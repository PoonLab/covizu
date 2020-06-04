import argparse
from gotoh2 import *


def update(srcfile, destfile, reffile, callback=None):
    """
    Incremental update of destination FASTA file by pairwise alignment
    of sequences from source FASTA file against a common reference.

    :param srcfile: open file stream to source FASTA in read mode
    :param destfile: open file stream to destination FASTA in append (r+) mode
    :param reffile: open file stream to FASTA with reference genome
    :param callback: optional, callback function for progress monitoring
    :return: int, number of genome sequences appended to destfile
    """
    # load SARS-COV-2 reference genome
    _, refseq = convert_fasta(reffile)[0]

    aligner = Aligner()
    count = update_alignment(
        ref=refseq, src=srcfile, dest=destfile, aligner=aligner,
        callback=callback
    )
    return count


def parse_args():
    """ Command-line interface """
    parser = argparse.ArgumentParser(
        description="Update previous FASTA alignment with new sequences by"
                    "pairwise alignment to a reference genome."
    )
    parser.add_argument('srcfile', type=argparse.FileType('r'),
                        help='input, FASTA file with new sequences')
    parser.add_argument('destfile', type=argparse.FileType('r+'),
                        help='input/output, FASTA file to append to')
    parser.add_argument('-ref', default=open('data/NC_045512.fa'),
                        type=argparse.FileType('r'),
                        help='Path to reference sequence.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    def progress(msg):
        print(msg)

    update(srcfile=args.srcfile, destfile=args.destfile, reffile=args.ref,
           callback=progress)
