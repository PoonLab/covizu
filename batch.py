import argparse
import os
import covizu
from covizu import minimap2, clustering, treetime, beadplot
from covizu.utils.seq_utils import Callback, convert_fasta
from covizu.utils import db_utils
from tempfile import NamedTemporaryFile
import json


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--db", type=str, default='data/gsaid.db')
    parser.add_argument("--ref", type=argparse.FileType('r'),
                        default=open(os.path.join(covizu.__path__[0], "data/MT291829.fa"))),
    parser.add_argument("--thread", type=int, default=1,
                        help="option, number of threads for minimap2.")
    parser.add_argument('--misstol', type=int, default=300,
                        help="option, maximum tolerated number of missing bases per "
                             "genome (default 300).")
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='option, specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')
    parser.add_argument('--ft2bin', default='fasttree2',
                        help='option, path to fasttree2 binary executable')
    parser.add_argument('--ttbin', default='treetime',
                        help='option, path to treetime binary executable')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cb = Callback()

    # Generate time-scaled tree of Pangolin lineages
    cb.callback("Retrieving genomes")
    fasta = treetime.retrieve_genomes(args.db, ref_file=args.ref, misstol=args.misstol)

    cb.callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = treetime.fasttree(fasta, binpath=args.ft2bin)

    cb.callback("Reconstructing time-scaled tree with {}").format(args.ttbin)
    nexus_file = treetime.treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                          clock=args.clock)
    treetime.parse_nexus(nexus_file, fasta, date_tol=args.datetol)

    # Retrieve raw genomes from DB, align and extract features
    cb.callback("Retrieving raw genomes from database")
    tmpfile = NamedTemporaryFile()
    db_utils.dump_raw(tmpfile, db=args.db)
    mm2 = minimap2.minimap2(tmpfile.name, ref=args.ref, nthread=args.thread)

    cb.callback("Extracting features and serializing to JSON")
    reflen = len(convert_fasta(open(args.ref))[0][1])
    res = []
    for qname, diffs, missing in minimap2.encode_diffs(mm2, reflen=reflen):
        res.append({'name': qname, 'diffs': diffs, 'missing': missing})
    serial = json.dumps(res).replace('},', '},\n')
    args.outfile.write(serial)

    lineages = db_utils.dump_lineages(args.db)
