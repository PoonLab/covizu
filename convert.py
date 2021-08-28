import argparse
import sys
from datetime import datetime, date
from csv import DictReader

from covizu.utils import seq_utils
from covizu.utils.progress_utils import Callback

import json
import gzip
import lzma
import codecs


def unescaped_str(arg_str):
    """ https://stackoverflow.com/questions/34145686/handling-argparse-escaped-character-as-option """
    return codecs.decode(str(arg_str), 'unicode_escape')


def parse_metadata(handle, set_region=None, delimiter=',', callback=None):
    """
    Parse CSV metadata input to dictionary

    :param handle:  open stream to input file in read mode
    :param delimiter:  set to '\t' for tab-separated value (TSV) files
    :param callback:  optional for progress monitoring
    :return:  dict, keyed by sequence name
    """
    reader = DictReader(handle, delimiter=delimiter)

    # check fieldnames
    fieldnames = [args.name, args.accession, args.coldate, args.country]
    if args.region:
        fieldnames.append(args.region)
    if args.division:
        fieldnames.append(args.division)

    for field in fieldnames:
        if field not in reader.fieldnames:
            if callback:
                callback("Missing fieldname {} in metadata CSV".format(field), level='ERROR')
                callback(reader.fieldnames, level='ERROR')
            sys.exit()

    metadata = {}
    for row in reader:
        metadata.update({row[args.name]: {
            'accession': row[args.accession],
            'coldate': row[args.coldate],
            'region': row.get(args.region, set_region),
            'country': row[args.country],
            'division': row.get(args.division, None),
            'lineage': row[args.lineage]
        }})

    return metadata


def combine(handle, metadata, minlen=29000, mindate='2019-12-01', callback=None):
    """
    Combine FASTA and metadata records into dictionary objects and do some basic filtering.

    :param handle:  file handle, open in read mode to FASTA
    :param metadata:  dict, returned from parse_metadata()
    :param minlen:  int, minimum genome length
    :param mindate:  str, earliest sample collection date in ISO-8601 format
    :param callback:  optional Callback() object

    :yield:  dict, combined genome and metadata for one record
    """
    mindate = seq_utils.fromisoformat(mindate)

    rejects = {'short': 0, 'baddate': 0}
    for label, seq in seq_utils.iter_fasta(handle):
        if label not in metadata:
            if callback:
                callback("Failed to retrieve metadata for genome {}".format(label), level='ERROR')
                sys.exit()
            sys.exit()

        if len(seq) < minlen:
            rejects['short'] += 1
            continue  # sequence is too short

        coldate = metadata[label]['coldate']
        if coldate.count('-') != 2:
            rejects['baddate'] += 1
            continue  # incomplete collection date
        dt = seq_utils.fromisoformat(coldate)
        if dt < mindate or dt > date.today():
            rejects['baddate'] += 1
            continue  # reject records with non-sensical collection date

        lineage = metadata[label]['lineage']
        if lineage is None:
            if callback:
                callback(
                    "Failed to retrieve lineage assignment for {}".format(label),
                    level='ERROR'
                )
            sys.exit()

        record = {
            'label': label,
            'accession': metadata[label]['accession'],
            'sequence': seq,
            'coldate': coldate,
            'region': metadata[label]['region'],
            'country': metadata[label]['country'],
            'division': metadata[label]['division'],
            'lineage': lineage
        }
        yield record

    if callback:
        callback("Rejected {short} short genomes, {baddate} records with bad dates".format(**rejects))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert and combine FASTA, Pangolin CSV and metadata files to "
                    "xz-compressed JSON, or append to an existing file."
    )

    parser.add_argument("infile", type=str,
                        help="input, path to FASTA file with genome sequences.")

    outputs = parser.add_mutually_exclusive_group()
    outputs.add_argument('-o', '--outfile', type=str,
                        help="path to output file to write xz-compressed JSON; otherwise "
                             "write to standard output to redirect to xz.")
    outputs.add_argument('-a', '--append', type=str,
                        help="path to optionally append output to existing xz-compressed JSON "
                             "file.")

    compression = parser.add_mutually_exclusive_group()
    compression.add_argument("--gz", action="store_true", help="FASTA is gzip-compressed.")
    compression.add_argument("--xz", action="store_true", help="FASTA is xz-compressed.")

    parser.add_argument("metadata", type=str,
                        help="input, CSV containing metadata, including Pango lineage")
    meta_compress = parser.add_mutually_exclusive_group()
    meta_compress.add_argument("--mgz", action="store_true", help="metadata file is gzip-compressed")
    meta_compress.add_argument("--mxz", action="store_true", help="metadata file is xz-compressed")

    parser.add_argument("--delimiter", type=unescaped_str, default='\t',
                        help="delimiter character for metadata TSV; "
                             "use ',' if comma-delimited")
    parser.add_argument("--name", type=str, default="strain",
                        help="column label for virus sample name in metadata; "
                             "required, must match FASTA")
    parser.add_argument("--accession", type=str, default="genbank_accession",
                        help="column label for accession number; required")
    parser.add_argument("--coldate", type=str, default="date",
                        help="column label for collection date in metadata CSV; required")
    parser.add_argument("--lineage", type=str, default="pango_lineage",
                        help="column label for Pango lineage classification; required")

    # geographical metadata
    parser.add_argument("--region", type=str, default=None,
                        help="column label for continent/region (e.g., Africa) in metadata CSV;"
                             "defaults to None (no region)")
    parser.add_argument("--set-region", type=str, default=None,
                        help="use if --region not specified; applies same region label to all records.")
    parser.add_argument("--country", type=str, default="country",
                        help="column label for country in metadata CSV")
    parser.add_argument("--division", type=str, default=None,
                        help="column label for country division (e.g., province) in metadata CSV; "
                             "defaults to None")

    parser.add_argument('--minlen', type=int, default=29000, help='minimum genome length (nt)')
    parser.add_argument('--mindate', type=str, default='2019-12-01',
                        help='earliest possible sample collection date (ISO format, default '
                             '2019-12-01')

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cb = Callback()
    callback = cb.callback

    # determine how to output results
    if args.outfile:
        outfile = lzma.open(args.outfile, 'wt')
    elif args.append:
        outfile = lzma.open(args.append, 'at')
    else:
        outfile = sys.stdout  # default
        sys.stderr.write("Streaming results to standard output, deactivating callback.")
        callback = None  # deactivate callback

    # process input metadata file
    if args.mgz:
        handle = gzip.open(args.metadata, 'rt')
    elif args.mxz:
        handle = lzma.open(args.metadata, 'rt')
    else:
        handle = open(args.metadata)
    metadata = parse_metadata(handle, set_region=args.set_region, delimiter=args.delimiter,
                              callback=callback)
    handle.close()

    # process input FASTA file
    if args.gz:
        handle = gzip.open(args.infile, 'rt')
    elif args.xz:
        handle = lzma.open(args.infile, 'rt')
    else:
        handle = open(args.infile)
    feed = combine(handle, metadata, minlen=args.minlen, mindate=args.mindate,
                   callback=callback)
    for record in feed:
        outfile.write(json.dumps(record)+'\n')
    handle.close()
    outfile.close()
