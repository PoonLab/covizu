"""Call minimap2 to align sequence"""
import subprocess
import argparse
import re
import sys
import os
import json
import pkg_resources
import covizu
from covizu.utils import gisaid_utils


def apply_cigar(sequence, rpos, cigar):
    """
    Use CIGAR to pad sequence with gaps as required to
    align to reference.  Adapted from http://github.com/cfe-lab/MiCall
    """
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)

    if not is_valid:
        raise RuntimeError(f'Invalid CIGAR string: {cigar!r}.')
    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
    aligned = '-' * rpos
    left = 0
    for length, operation in tokens:
        length = int(length)
        if operation in 'M=X':
            aligned += sequence[left:(left + length)]
            left += length
        elif operation == 'D':
            aligned += '-' * length
        elif operation in 'SI':
            left += length  # soft clip

    return aligned


def minimap2(
        infile,
        ref,
        stream=False,
        path='minimap2',
        nthread=3,
        minlen=29000):
    """
    Wrapper function for minimap2.

    :param infile:  file object or StringIO
    :param ref:  str, path to FASTA with reference sequence(s)
    :param stream:  bool, if True then stream data from <infile> object via stdin
    :param path:  str, path to binary executable
    :param nthread:  int, number of threads for parallel execution of minimap2
    :param minlen:  int, filter genomes below minimum length; to accept all, set to 0.

    :yield:  query sequence name, reference index, CIGAR and original
             sequence
    """
    if stream:
        # input from StringIO in memory
        with subprocess.Popen(
            [path, '-t', str(nthread), '-a', '--eqx', ref, '-'], encoding='utf8',
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        ) as process:
            output = process.communicate(infile)[0]
            output = output.split('\n')
    else:
        # input read from file
        with subprocess.Popen(
            [path, '-t', str(nthread), '-a', '--eqx', ref, infile.name],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        ) as process:
            output = map(lambda x: x.decode('utf-8'), process.stdout)

    for out_line in output:
        if out_line == '' or out_line.startswith('@'):
            # split on \n leaves empty line; @ prefix header lines
            continue
        q_name, flag, rname, rpos, _, cigar, _, _, _, sequence = \
            out_line.strip('\n').split('\t')[:10]

        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            continue

        if len(sequence) < minlen:
            # reject sequence that is too short
            continue

        # validate CIGAR string
        if not re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar):
            raise RuntimeError(f'Invalid CIGAR string: {cigar!r}.')

        rpos = int(rpos) - 1  # convert to 0-index
        yield q_name, rpos, cigar, sequence


# return aligned sequence?
def output_fasta(iterable, outfile, reflen=0):
    """
    Stream output from minimap2 into FASTA file
    of aligned sequences.  CIGAR parsing code adapted from
    http://github.com/cfe-lab/MiCall

    :param iter:  generator from minimap2()
    :param outfile:  open file stream in write mode
    :param reflen:  int, length of reference genome to pad sequences;
                    defaults to no padding.
    """
    for q_name, rpos, cigar, sequence in iterable:
        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        aligned = '-' * rpos
        left = 0
        for length, operation in tokens:
            length = int(length)
            if operation in 'M=X':
                aligned += sequence[left:(left + length)]
                left += length
            elif operation == 'D':
                aligned += '-' * length
            elif operation in 'SI':
                left += length  # soft clip

        # pad on right
        aligned += '-' * (reflen - len(aligned))
        outfile.write(f'>{q_name}\n{aligned}\n')


def stream_fasta(iterable, reflen=0):
    """
    Stream output from minimap2 into list of tuples [(header, seq), ... , (header, seq)]
    of aligned sequences.  CIGAR parsing code adapted from
    http://github.com/cfe-lab/MiCall

    :param iter:  generator from minimap2()
    :param reflen:  int, length of reference genome to pad sequences;
                    defaults to no padding.
    :yield:  tuple, header and aligned sequence
    """
    for q_name, rpos, cigar, sequence in iterable:
        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        aligned = '-' * rpos
        left = 0
        for length, operation in tokens:
            length = int(length)
            if operation in 'M=X':
                aligned += sequence[left:(left + length)]
                left += length
            elif operation == 'D':
                aligned += '-' * length
            elif operation in 'SI':
                left += length  # soft clip

        # pad on right
        aligned += '-' * (reflen - len(aligned))
        yield q_name, aligned


def encode_diffs(iterable, reflen, alphabet='ACGT'):
    """
    Serialize differences of query sequences to reference
    genome, which comprise nucleotide substitutions, in-frame
    indels, and locations of missing data.
    NOTE: runs of 'N's are also represented by 'X' tokens in the CIGAR
    string.
    :param iter:  generator from minimap2()
    :param reflen:  length of reference genome
    """
    for q_name, rpos, cigar, sequence in iterable:
        differences = []
        missing_data = []
        pos = [0, rpos] # [index for query (left), rpos]
        if pos[1] > 0:
            # incomplete on left
            missing_data.append(tuple([0, pos[1]]))

        for length, operator in re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE):
            length = int(length)
            substr = sequence[pos[0]:(pos[0] + length)]
            if operator == 'X':
                # each nucleotide is a separate diff
                if 'N' in substr:
                    # for now, assume the whole substring is bs
                    missing_data.append(tuple([pos[1], pos[1] + length]))
                else:
                    # assume adjacent mismatches are independent substitutions
                    for i, nucleotide in enumerate(substr):
                        if nucleotide in alphabet:
                            differences.append(tuple(['~', pos[1] + i, nucleotide]))
                        else:
                            # skip ambiguous base calls, like "R"
                            missing_data.append(tuple([pos[1] + i, pos[1] + i + 1]))

                pos[0] += length
                pos[1] += length
            else:
                interpret_operator(operator, substr, pos, length, differences)

        # update missing if sequence incomplete on the right
        if pos[1] < reflen:
            missing_data.append(tuple([pos[1], reflen]))

        yield q_name, differences, missing_data

def interpret_operator(operator, substr, pos, length, differences):
    """interpret operator from encode_diff() to make pep8 compliance"""
    if operator == 'S':
        # discard soft clip
        pos[0] += length
    elif operator == 'I':
        # insertion relative to reference
        differences.append(tuple(['+', pos[1], substr]))
        pos[0] += length
    elif operator == 'D':
        # deletion relative to reference
        differences.append(tuple(['-', pos[1], length]))
        pos[1] += length
    elif operator == '=':
        # exact match
        pos[0] += length
        pos[1] += length
    elif operator == 'H':
        # hard clip, do nothing
        pass
    else:
        print(f"ERROR: unexpected operator {operator}")
        sys.exit()

def parse_args():
    """parse input args"""
    parser = argparse.ArgumentParser("Wrapper script for minimap2")
    parser.add_argument('infile', type=argparse.FileType(
        'r'), help="<input> path to query FASTA file or database (--stream)")
    parser.add_argument('-o', '--outfile',
                        type=argparse.FileType('w'),
                        required=False,
                        help="<output, optional> path to write output, "
                             "defaults to stdout")
    parser.add_argument('-a', '--align', action='store_true',
                        help="<option> output aligned sequences as FASTA")
    parser.add_argument('-t', '--thread', type=int, default=3,
                        help="<option> number of threads")
    parser.add_argument('-f', '--force-headers', action='store_true',
                        help="<option> use -f to force this script to accept "
                             "headers with spaces, which will be truncated "
                             "by minimap2")
    parser.add_argument('--filter', action='store_true',
                        help="<option> filter problematic sites")

    path = covizu.__path__[0]
    if '.egg' in path:
        ref_file = pkg_resources.resource_filename(
            'covizu', 'data/NC_045512.fa')
        vcf_file = pkg_resources.resource_filename(
            'covizu', "data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf")
    else:
        ref_file = os.path.join(path, "data", "NC_045512.fa")
        vcf_file = os.path.join(path,
                                "data/ProblematicSites_SARS-CoV2",
                                "problematic_sites_sarsCov2.vcf")

    parser.add_argument('--ref', type=str, default=ref_file,
                        help="<input> path to target FASTA (reference)")
    parser.add_argument(
        "--vcf",
        type=str,
        default=vcf_file,
        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
        "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument('--minlen', help="<option> minimum sequence length, "
                                         "defaults to 29000nt.",
                        type=int, default=29000)
    return parser.parse_args()


if __name__ == '__main__':
    from covizu.utils.seq_utils import convert_fasta, apply_features
    args = parse_args()
    if args.outfile is None:
        args.outfile = sys.stdout

    # check input headers for spaces
    if not args.force_headers:
        for line in args.infile:
            if line.startswith('>') and ' ' in line:
                print("WARNING: at least one FASTA header contains a space")
                print(line)
                print("Use '-f' to force this script to process the file")
                print(
                    "Otherwise use `sed -i 's/ /_/g' <file>` to replace all spaces in place.")
                sys.exit()

    mm2 = minimap2(args.infile, ref=args.ref, nthread=args.thread,
                   minlen=args.minlen)

    # get length of reference
    with open(args.ref, encoding='utf-8') as in_fasta:
        refseq = convert_fasta(in_fasta)[0][1]

    len_ref = len(refseq)
    if args.align:
        if args.filter:
            encoded = encode_diffs(mm2, reflen=len_ref)
            for qname, diffs, missing in gisaid_utils.filter_problematic(
                    encoded, vcf_file=args.vcf, encoded=True):
                seq = apply_features(diffs, missing, refseq=refseq)
                args.outfile.write(f">{qname}\n{seq}\n")
        else:
            output_fasta(mm2, reflen=len_ref, outfile=args.outfile)
    else:
        # serialize feature vectors as JSON
        res = []
        for qname, diffs, missing in encode_diffs(mm2, reflen=len_ref):
            res.append({'name': qname, 'diffs': diffs, 'missing': missing})
        serial = json.dumps(res).replace('},', '},\n')
        args.outfile.write(serial)
