import subprocess
import argparse
import re
import sys
from seq_utils import convert_fasta


def apply_cigar(seq, rpos, cigar):
    """
    Use CIGAR to pad sequence with gaps as required to
    align to reference.  Adapted from http://github.com/cfe-lab/MiCall
    """
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)

    if not is_valid:
        raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))
    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
    aligned = '-'*rpos
    left = 0
    for length, operation in tokens:
        length = int(length)
        if operation in 'M=X':
            aligned += seq[left:(left+length)]
            left += length
        elif operation == 'D':
            aligned += '-'*length
        elif operation in 'SI':
            left += length  # soft clip

    return aligned


def minimap2(fasta, ref, path='minimap2', nthread=3):
    """
    Wrapper function for minimap2
    :param fasta: str, path to FASTA with query sequences
    :param ref: str, path to FASTA with reference sequence(s)
    :param path: str, path to binary executable
    :yield:  query sequence name, reference index, CIGAR and original
             sequence
    """
    p = subprocess.Popen([path, '-t', str(nthread), '-a', '--eqx', ref, fasta],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in map(lambda x: x.decode('utf-8'), p.stdout):
        if line.startswith('@'):
            continue
        qname, flag, rname, rpos, _, cigar, _, _, _, seq = \
            line.strip().split()[:10]
        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            continue

        # validate CIGAR string
        is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
        if not is_valid:
            raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq


# return aligned sequence?
def output_fasta(iter, outfile, reflen=0):
    """
    Stream output from minimap2 into FASTA file
    of aligned sequences.  CIGAR parsing code adapted from
    http://github.com/cfe-lab/MiCall

    :param iter:  generator from minimap2()
    :param outfile:  open file stream in write mode
    :param reflen:  int, length of reference genome to pad sequences;
                    defaults to no padding.
    """
    for qname, rpos, cigar, seq in iter:
        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        aligned = '-' * rpos
        left = 0
        for length, operation in tokens:
            length = int(length)
            if operation in 'M=X':
                aligned += seq[left:(left + length)]
                left += length
            elif operation == 'D':
                aligned += '-' * length
            elif operation in 'SI':
                left += length  # soft clip

        # pad on right
        aligned += '-'*(reflen-len(aligned))
        outfile.write('>{}\n{}\n'.format(qname, aligned))


def compress_intervals(iseq):
    """
    Compress consecutive integers in a sequence into
    intervals.  Return as a serialized string.
    :param iseq:  list, integer sequence
    """
    pass


def encode_diffs(iter, reflen):
    """
    Serialize differences of query sequences to reference
    genome, which comprise nucleotide substitutions, in-frame
    indels, and locations of missing data.
    :param iter:  generator from minimap2()
    """
    for qname, rpos, cigar, seq in iter:
        diffs = []
        left = 0
        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        for length, operation in tokens:
            if operation == 'X':
                # each nucleotide is a separate diff
                pass                


def parse_args():
    parser = argparse.ArgumentParser("Wrapper script for minimap2")
    parser.add_argument('fasta', type=argparse.FileType('r'),
                        help="<input> path to query FASTA file")
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
    parser.add_argument('--ref', help="<input> path to target FASTA (reference)",
                        default='data/NC_045512.fa')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.outfile is None:
        args.outfile = sys.stdout

    # check input headers for spaces
    for line in args.fasta:
        if line.startswith('>') and ' ' in line:
            print("WARNING: at least one FASTA header contains a space")
            print(line)
            print("Use '-f' to force this script to process the file")
            print("Otherwise use `sed -i 's/ /_/g' <file>` to replace all spaces in place.")
            sys.exit()

    # get length of reference
    reflen = len(convert_fasta(open(args.ref))[0][1])
    mm2 = minimap2(args.fasta.name, ref=args.ref, nthread=args.thread)

    if args.align:
        output_fasta(mm2, reflen=reflen, outfile=args.outfile)
    else:
        encode_diffs(mm2, reflen=reflen)

