import sys
import lzma
import os
import csv
import subprocess
import tempfile
from mpi4py import MPI

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
nprocs = comm.Get_size()
tempdir = tempfile.gettempdir()
fieldnames = ['taxon', 'lineage', 'conflict', 'ambiguity_score', 'scorpio_call',
    'scorpio_support', 'scorpio_conflict', 'version', 'pangolin_version',
    'pangoLEARN_version', 'pango_version', 'status', 'note']


def pangolin(fasta):
    # write sequences to temporary file
    tmpfile = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    for h, s in fasta:
        tmpfile.write(">{}\n{}\n".format(h, s))
    tmpfile.close()

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
    if len(sys.argv) != 2:
        print("Usage: mpirun -np <number of processes> python3 mangolin.py <fasta.xz>")
        sys.exit()

    # create output file
    outfile = open("mangolin.{}.csv".format(my_rank), 'w')
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    handle = lzma.open(sys.argv[1], 'rt')
    for fasta in chunk_fasta(iter_fasta(handle)):
        sys.stdout.write("({}/{}) processing {} sequences\n".format(
            my_rank, nprocs, len(fasta)))
        results = pangolin(fasta)
        for row in results:
            writer.writerow(row)
    outfile.close()

