"""
Rapid pairwise alignment of SARS-COV-2 sequences to reference,
eliminating insertions to enforce reference coordinate system.
"""

from gotoh2 import Aligner
import argparse
import re
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()


# default alignment settings
aligner = Aligner()

refseq = ''
with open("data/NC_045512.fa") as handle:
    _ = next(handle)
    for line in handle:
        refseq += line.strip()

# regex for terminal gaps
pat = re.compile("^(-*)[^-].+[^-](-*)$")


def len_terminal_gaps(seq):
    """
    Calculate lengths of terminal (prefix, suffix) gaps
    :param seq:  str, gapped nucleotide/amino acid sequence
    :return:  tuple, lengths of gap prefix and suffix
    """
    m = pat.findall(seq)
    return tuple(map(len, m[0]))


def iter_fasta (handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.
    """
    sequence = ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''
            h = line.strip('>\n')
        else:
            sequence += line.strip('\n').upper()
    yield h, sequence


def pair_align(ref, query, threshold=3.0):
    """
    Pairwise alignment of query against reference.  Any insertions in the query
    relative to the reference are excised and returned separately.  Filters
    sequences with low alignment scores (when scaled to sequence length, a good
    score should be around 5.0 - the minimum is 0).

    :param ref:  str, reference sequence
    :param query:  str, query sequence
    :return: str, list - aligned and trimmed query, list of insertions
    """
    aref, aquery, ascore = aligner.align(ref, query)
    if ascore / len(aref) < threshold:
        # bad alignment
        return None, None

    # if query is longer than reference, do not count terminal gaps as insertions
    left, right = len_terminal_gaps(aref)

    trim_seq = ''
    inserts = []
    for i in range(left, len(aref)-right):
        rn = aref[i]
        qn = aquery[i]
        if rn == '-':
            # insertion relative to reference
            inserts.append((i, qn))
            continue
        trim_seq += qn

    return trim_seq, inserts


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rapid pairwise alignment of SARS-COV-2 genomes against"
                    "reference sequence.  Requires about 8GB RAM."
    )
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="input, FASTA file of sequences to align")
    parser.add_argument("outfile",
                        help="output, FASTA file of aligned sequences")
    parser.add_argument("insfile",
                        help="output, text file of sequence insertions")
    parser.add_argument("--threshold", '-t', type=float, default=3.0,
                        help="optional, scaled alignment score threshold")
    return parser.parse_args()


def main():
    """
    Command-line interface
    """
    args = parse_args()

    outfile = open('{}.{}'.format(args.outfile, my_rank), 'w')
    insfile = open('{}.{}'.format(args.insfile, my_rank), 'w')

    for i, (h, s) in enumerate(iter_fasta(args.infile)):
        if i % nprocs != my_rank:
            continue

        print('Process {} of {} running {}'.format(my_rank, nprocs, h))

        trim_seq, inserts = pair_align(refseq, s, threshold=args.threshold)
        if trim_seq is None:
            print("{} failed alignment", h)
            continue

        outfile.write('>{}\n{}\n'.format(h, trim_seq))
        outfile.flush()
        for pos, ins in inserts:
            insfile.write('{},{},{}\n'.format(h, pos, ins))
            insfile.flush()

    outfile.close()
    insfile.close()


if __name__ == '__main__':
    main()
