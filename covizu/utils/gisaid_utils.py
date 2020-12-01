import lzma
import json
from datetime import date
import argparse
import os

import covizu
from covizu.minimap2 import minimap2, encode_diffs
from covizu.utils.seq_utils import fromisoformat, convert_fasta


def load_gisaid(path, minlen=29000, mindate='2019-12-01'):
    """
    Read in GISAID feed as xz compressed JSON, applying some basic filters

    :param path:  str, path to xz-compressed JSON
    :param minlen:  int, minimum genome length
    :param mindate:  datetime.date, earliest reasonable sample collection date
    :yield:  dict, contents of each GISAID record
    """
    mindate = fromisoformat(mindate)
    with lzma.open(path, 'rb') as handle:
        for line in handle:
            record = json.loads(line)

            qname = record['covv_virus_name']
            if qname.split('/')[1][0].islower():
                # reject non-human isolates
                # FIXME: request host field
                continue

            seq = record['sequence'].replace('\n', '')
            if len(seq) < minlen:
                # reject sequences that are too short
                continue
            record['sequence'] = seq

            if record['covv_collection_date'].count('-') != 2:
                # reject sequences without complete collection date
                continue
            coldate = fromisoformat(record['covv_collection_date'])
            if coldate < mindate or coldate > date.today():
                # reject sequences with nonsense collection date
                continue

            yield record


def batch_fasta(gen, size=100):
    """
    Concatenate sequence records in stream into FASTA-formatted text
    """
    stdin = ''
    batch = []
    for i, record in enumerate(gen):
        sequence = record.pop('sequence')
        stdin += '>{}\n{}\n'.format(record['covv_virus_name'], sequence)
        batch.append(record)
        if i > 0 and i % size == 0:
            yield stdin, batch
            stdin = ''
            batch = []


def align_online(batcher, ref_file, binpath='minimap2', nthread=3, minlen=29000):
    """
    Stream output from JSON.xz file via load_gisaid() into minimap2
    via subprocess.
    :param batcher:  generator, returned by batch_fasta()
    :param ref_file:  str, path to reference genome (FASTA format)
    :param binpath:  str, path to minimap2 binary executable
    :param nthread:  int, number of threads to run minimap2
    :param minlen:  int, minimum genome length

    :return:  dict,
    """
    with open(ref_file) as handle:
        reflen = len(convert_fasta(handle)[0][1])

    for fasta, batch in batcher:
        mm2 = minimap2(fasta, ref_file, stream=True, path=binpath, nthread=nthread,
                       minlen=minlen)
        result = list(encode_diffs(mm2, reflen=reflen))
        yield result, batch


def parse_args():
    """ Command line help text"""
    parser = argparse.ArgumentParser("")
    parser.add_argument('infile', type=str, help="input, path to xz-compressed JSON")
    parser.add_argument('outfile', type=argparse.FileType('w'), help="output, path to write JSON")

    parser.add_argument('--minlen', type=int, default=29000, help='option, minimum genome length')
    parser.add_argument('--mindate', type=str, default='2019-12-01',
                        help='option, earliest possible sample collection date (ISO format, default '
                             '2019-12-01)')
    parser.add_argument('--ref', type=str, help="option, path to reference genome (FASTA)",
                        default=os.path.join(covizu.__path__[0], "data/NC_045512.fa"))
    parser.add_argument('--binpath', type=str, default='minimap2',
                        help="option, path to minimap2 binary executable file")
    parser.add_argument('--nthread', type=int, default=2, help='option, number of threads to run minimap2')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    loader = load_gisaid(args.infile, minlen=args.minlen, mindate=args.mindate)
    batcher = batch_fasta(loader, size=100)
    aligner = align_online(batcher, ref_file=args.ref, binpath=args.binpath, nthread=args.nthread, minlen=args.minlen)
    for result, batch in aligner:
        for aln, row in zip(result, batch):
            qname, diffs, missing = aln
            row.update({'diff': diffs, 'missing': missing})
            args.outfile.write(json.dumps(row)+'\n')

