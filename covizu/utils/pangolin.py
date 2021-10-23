import sys
import lzma
import os
import csv
import subprocess
import tempfile
import argparse

tempdir = tempfile.gettempdir()
fieldnames = ['taxon', 'lineage', 'conflict', 'ambiguity_score', 'scorpio_call',
              'scorpio_support', 'scorpio_conflict', 'version', 'pangolin_version',
              'pangoLEARN_version', 'pango_version', 'status', 'note']


def pangolin(fasta):
    # write sequences to temporary file as input
    tmpfile = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    for h, s in fasta:
        tmpfile.write(">{}\n{}\n".format(h, s))
    tmpfile.close()

    # FIXME: can we capture the output stream?
    outfile = os.path.join(tempdir, 'mango{}.csv'.format(my_rank))
    subprocess.check_call(['pangolin', tmpfile.name, '--outfile', outfile],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # import results
    reader = csv.DictReader(open(outfile))
    results = [row for row in reader]
    return results


def iter_fasta(handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.
    :param handle:  open stream to FASTA file in read mode
    :yield tuples, (header, sequence)
    """
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    yield h, sequence


def chunk_fasta(fasta, chunk_size=1000):
    """
    Partition output stream by process number into chunks of a given size.
    """
    chunk = []
    for i, record in enumerate(fasta):
        if i % nprocs != my_rank:
            continue
        chunk.append(record)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    yield chunk


if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser("MPI wrapper for SARS-CoV-2 lineage classification with Pangolin")
    parser.add_argument('infile', type=str, help="<input> path to xz-compressed FASTA file")
    parser.add_argument('outfile', type=str, help="<output> path to write CSV output; will be modified "
                                                  "with MPI rank and .csv extension")
    parser.add_argument('--mpi', action='store_true', help="<option> run in MPI mode")
    parser.add_argument('-n', '--size', type=int, default=1000,
                        help="<option> chunk size for streaming from FASTA")
    parser.add_argument('--quiet', action='store_true', help="option, suppress console messages")
    args = parser.parse_args()

    if args.mpi:
        try:
            from mpi4py import MPI
        except ModuleNotFoundError:
            print("Script requires mpi4py - https://pypi.org/project/mpi4py/")
            sys.exit()

        comm = MPI.COMM_WORLD
        my_rank = comm.Get_rank()
        nprocs = comm.Get_size()
        outfile = open("{}.{}.csv".format(args.outfile, my_rank), 'w')
    else:
        my_rank = 0
        nprocs = 1
        outfile = sys.stdout

    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()

    # stream input from FASTA file
    handle = lzma.open(args.infile, 'rt')
    for fasta in chunk_fasta(iter_fasta(handle), chunk_size=args.size):
        if not args.quiet:
            sys.stderr.write("({}/{}) processing {} sequences\n".format(
                my_rank, nprocs, len(fasta)))
            sys.stderr.flush()
        results = pangolin(fasta)
        for row in results:
            writer.writerow(row)
    outfile.close()
