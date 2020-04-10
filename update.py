import argparse
from gotoh2 import *

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

args = parser.parse_args()


def progress(msg):
    print(msg)


# load SARS-COV-2 reference genome
_, refseq = convert_fasta(args.ref)[0]

aligner = Aligner()
update = update_alignment(
    ref=refseq, src=args.srcfile, dest=args.destfile, aligner=aligner,
    callback=progress
)
