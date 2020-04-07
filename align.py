"""
Rapid pairwise alignment of SARS-COV-2 sequences to reference,
eliminating insertions to enforce reference coordinate system.
"""

from gotoh2 import Aligner
import argparse
import re

# default alignment settings
aligner = Aligner()

refseq = ''
with open("data/NC_045512.fa") as handle:
    _ = next(handle)
    for line in handle:
        refseq += line.strip()

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


def pair_align(ref, query, threshold=0):
    """
    Pairwise alignment of query against reference.  Any insertions in the query
    relative to the reference are excised and returned separately.

    :param ref:  str, reference sequence
    :param query:  str, query sequence
    :return: str, list - aligned and trimmed query, list of insertions
    """
    aref, aquery, ascore = aligner.align(ref, query)
    print(ascore)

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


def main():
    """
    Command-line interface
    """
    parser = argparse.ArgumentParser(
        description="Rapid pairwise alignment of SARS-COV-2 genomes against"
                    "reference sequence."
    )
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="input, FASTA file of sequences to align")
    parser.add_argument("outfile", type=argparse.FileType('w'),
                        help="output, FASTA file of aligned sequences")
    parser.add_argument("insfile", type=argparse.FileType('w'),
                        help="output, text file of sequence insertions")
    parser.add_argument("--threshold", '-t', type=float, default=0.1,
                        help="optional, scaled alignment score threshold")
    args = parser.parse_args()

    for h, s in iter_fasta(args.infile):
        print(h)
        trim_seq, inserts = pair_align(refseq, s, threshold=args.threshold)
        if trim_seq is None:
            print("{} failed alignment", h)
            continue

        args.outfile.write('>{}\n{}\n'.format(h, trim_seq))

if __name__ == '__main__':
    main()
