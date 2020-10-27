import argparse
import os
import covizu
from covizu import minimap2, clustering, treetime, beadplot
from covizu.utils import db_utils, seq_utils
from tempfile import NamedTemporaryFile
import json
import itertools


def parse_args():
    parser = argparse.ArgumentParser(description="CoVizu analysis pipeline automation")

    parser.add_argument("outfile", type=argparse.FileType('w'),
                        help="output, dest for JSON beadplot file")

    parser.add_argument("--db", type=str, default='data/gsaid.db',
                        help="input, path to sqlite3 database")
    parser.add_argument("--ref", type=str,
                        default=os.path.join(covizu.__path__[0], "data/MT291829.fa"),
                        help="input, path to FASTA file with reference genome"),

    parser.add_argument('--outdir', default='data/',
                        help='optional, directory to write TreeTime output files')
    parser.add_argument('--misstol', type=int, default=300,
                        help="option, maximum tolerated number of missing bases per "
                             "genome (default 300).")
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='option, specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')
    parser.add_argument('--datetol', type=float, default=0.1,
                        help='optional, exclude tips from time-scaled tree '
                             'with high discordance between estimated and '
                             'known sample collection dates (year units,'
                             'default: 0.1)')

    parser.add_argument('--ft2bin', default='fasttree2',
                        help='option, path to fasttree2 binary executable')
    parser.add_argument('--ttbin', default='treetime',
                        help='option, path to treetime binary executable')

    parser.add_argument('-mmt', "--mmthreads", type=int, default=1,
                        help="option, number of threads for minimap2.")
    parser.add_argument('-njt', "--njthreads", type=int, default=1,
                        help="option, number of threads for NJ reconstruction")
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")
    parser.add_argument("--cutoff", type=float, default=0.5,
                        help="Bootstrap cutoff for consensus tree (default 0.5). "
                             "Only used if --cons is specified.")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cb = seq_utils.Callback()

    # Generate time-scaled tree of Pangolin lineages
    cb.callback("Retrieving lineage genomes")
    fasta = treetime.retrieve_genomes(args.db, nthread=args.mmthreads, ref_file=args.ref,
                                      misstol=args.misstol)

    cb.callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = treetime.fasttree(fasta, binpath=args.ft2bin)

    cb.callback("Reconstructing time-scaled tree with {}".format(args.ttbin))
    nexus_file = treetime.treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                                   clock=args.clock)
    treetime.parse_nexus(nexus_file, fasta, date_tol=args.datetol)
    # output timetree.nwk

    # Retrieve raw genomes from DB, align and extract features
    cb.callback("Retrieving raw genomes from database")
    with NamedTemporaryFile(prefix="cvz_") as tmpfile:
        db_utils.dump_raw(outfile=tmpfile.name, db=args.db)
        mm2 = minimap2.minimap2(tmpfile, ref=args.ref, nthread=args.mmthreads)

        cb.callback("Extracting features and serializing to JSON")
        reflen = len(seq_utils.convert_fasta(open(args.ref))[0][1])
        features = []
        for qname, diffs, missing in minimap2.encode_diffs(mm2, reflen=reflen):
            features.append([qname, diffs, missing])

    # Neighbor-joining reconstruction
    lineages = db_utils.dump_lineages(args.db)
    result = []
    by_lineage = clustering.split_by_lineage(features, lineages)
    for lineage, feats in by_lineage.items():
        cb.callback('start {}, {} entries'.format(lineage, len(feats)))
        filtered = seq_utils.filter_outliers(feats)

        # FIXME: reduce data for debugging
        filtered = itertools.islice(filtered, 100)

        trees, labels = clustering.build_trees(
            filtered, nboot=args.nboot, threads=args.njthreads, callback=cb.callback
        )
        ctree = clustering.consensus(trees, cutoff=args.cutoff)
        label_dict = dict([(idx, lst) for idx, lst in enumerate(labels)])
        ctree = beadplot.annotate_tree(ctree, label_dict)
        result.append(beadplot.serialize_tree(ctree))

    args.outfile.write(json.dumps(result, indent=2))
