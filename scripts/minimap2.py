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


def minimap2(fasta, ref, path='minimap2', retseq=False):
    """
    Wrapper function for minimap2
    :param fasta: str, path to FASTA with query sequences
    :param ref: str, path to FASTA with reference sequence(s)
    :param path: str, path to binary executable
    """
    p = subprocess.Popen([path, '-a', '--eqx', ref, fasta],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in map(lambda x: x.decode('utf-8'), p.stdout):
        if line.startswith('@'):
            continue
        qname, _, rname, rpos, _, cigar, _, _, _, seq = \
            line.strip().split()[:10]
        if rname == '*':
            # did not map
            continue
        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq


# return aligned sequence?
def output_fasta(iter, reflen, outfile):
    for qname, rpos, cigar, seq in iter:
        aligned = apply_cigar(seq, rpos, cigar)
        # pad on right
        aligned += '-'*(reflen-len(aligned))
        outfile.write('>{}\n{}\n'.format(qname, aligned))


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
    parser.add_argument('--ref', help="<input> path to target FASTA (reference)",
                        default='../data/NC_045512.fa')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.outfile is None:
        args.outfile = sys.stdout

    mm2 = minimap2(args.fasta.name, ref=args.ref)
    if args.align:
        # get length of reference
        reflen = len(convert_fasta(open(args.ref))[0][1])
        output_fasta(mm2, reflen=reflen, outfile=args.outfile)
