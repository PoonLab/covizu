import argparse
import os
import sys
import json
from datetime import datetime
import lzma
import covizu

from covizu.utils import seq_utils
from covizu.utils.batch_utils import *
from covizu.utils.progress_utils import Callback
from covizu.minimap2 import extract_features


def parse_args():
    parser = argparse.ArgumentParser(
        description="CoVizu analysis pipeline automation for execution on local files"
    )

    parser.add_argument("infile", type=str,
                        help="input, path to xz-compressed JSON file with genomes and metadata,"
                             " see convert.py")

    parser.add_argument("--outdir", type=str, default='data/',
                        help="option, path to write output files")
    parser.add_argument("--bylineage", type=str, default='data/by_lineage.json',
                        help="path to write JSON of features by lineage")

    parser.add_argument('--poisson-cutoff', type=float, default=0.001,
                        help='filtering outlying genomes whose distance exceeds the upper '
                             'quantile of Poisson distribution (molecular clock).  Default 0.001 '
                             'corresponds to 99.9%% cutoff.')
    parser.add_argument('--minlen', type=int, default=29000, help='minimum genome length (nt)')

    parser.add_argument('--batchsize', type=int, default=500,
                        help='number of records to batch process with minimap2')

    parser.add_argument("--ref", type=str,
                        default=os.path.join(covizu.__path__[0], "data/NC_045512.fa"),
                        help="path to FASTA file with reference genome")
    parser.add_argument('--mmbin', type=str, default='minimap2',
                        help="path to minimap2 binary executable")
    parser.add_argument('-mmt', "--mmthreads", type=int, default=8,
                        help="number of threads for minimap2.")

    parser.add_argument('--misstol', type=int, default=300,
                        help="maximum tolerated number of missing bases per "
                             "genome (default 300).")
    parser.add_argument("--vcf", type=str,
                        default=os.path.join(covizu.__path__[0], "data/problematic_sites_sarsCov2.vcf"),
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument('--ft2bin', default='fasttree2',
                        help='path to fasttree2 binary executable')

    parser.add_argument('--ttbin', default='treetime',
                        help='path to treetime binary executable')
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')
    parser.add_argument('--earliest', action='store_true',
                        help='option, use earliest sample per lineage for time-scaled '
                             'tree; otherwise defaults to most recent samples.')

    parser.add_argument('--datetol', type=float, default=0.1,
                        help='exclude tips from time-scaled tree '
                             'with high discordance between estimated and '
                             'known sample collection dates (year units,'
                             'default: 0.1)')

    parser.add_argument('--binpath', type=str, default='rapidnj',
                        help='path to RapidNJ binary executable')
    parser.add_argument('--mincount', type=int, default=500,
                        help='minimum number of variants in lineage '
                             'above which MPI processing will be used.')
    parser.add_argument('--machine_file', type=str, default='mfile',
                        help='path to machine file for MPI.')
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")

    parser.add_argument("--boot-cutoff", type=float, default=0.5,
                        help="Bootstrap cutoff for consensus tree (default 0.5). "
                             "Only used if --cons is specified.")

    return parser.parse_args()


def analyze_feed(handle, args, callback=None):
    """
    :param handle:  file stream in read mode, from lzma.open()
    :param args:  Namespace
    :param callback:  optional progress monitoring, see progress_utils.py
    """

    # check that user has loaded openmpi module
    try:
        subprocess.check_call(['mpirun', '-np', '2', 'ls'], stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        if callback:
            callback("mpirun not loaded - run `module load openmpi/gnu`", level='ERROR')
        sys.exit()

    # pre-processing feed
    feed = map(json.loads, handle)
    batcher = seq_utils.batch_fasta(feed, size=args.batchsize)
    aligned = extract_features(batcher, ref_file=args.ref, binpath=args.mmbin,
                               nthread=args.mmthreads, minlen=args.minlen)
    filtered = seq_utils.filter_problematic(aligned, vcf_file=args.vcf, cutoff=args.poisson_cutoff,
                                            callback=callback)
    by_lineage = sort_by_lineage(filtered, callback=callback)

    with open(args.bylineage, 'w') as handle:
        # export to file to process large lineages with MPI
        json.dump(by_lineage, handle)

    # reconstruct time-scaled tree
    timetree, residuals = build_timetree(by_lineage, args, callback)
    timestamp = datetime.now().isoformat().split('.')[0]
    nwk_file = os.path.join(args.outdir, 'timetree.{}.nwk'.format(timestamp))
    with open(nwk_file, 'w') as handle:
        Phylo.write(timetree, file=handle, format='newick')

    # generate beadplots and serialize to file
    result = make_beadplots(by_lineage, args, callback, t0=timestamp)
    outfile = open(os.path.join(args.outdir, 'clusters.{}.json'.format(timestamp)), 'w')
    outfile.write(json.dumps(result))  # serialize results to JSON
    outfile.close()

    # write data stats
    dbstat_file = os.path.join(args.outdir, 'dbstats.{}.json'.format(timestamp))
    with open(dbstat_file, 'w') as handle:
        nseqs = sum([len(rows) for rows in by_lineage.values()])
        val = {
            'lastupdate': timestamp.split('T')[0],
            'noseqs': nseqs,
            'lineages': {}
        }
        for lineage, samples in by_lineage.items():
            ndiffs = [len(x['diffs']) for x in samples]
            val['lineages'][lineage] = {
                'nsamples': len(samples),
                'lastcoldate': max(x['coldate'] for x in samples),
                'residual': residuals[lineage],
                'max_ndiffs': max(ndiffs),
                'mean_ndiffs': sum(ndiffs) / len(ndiffs)
            }
        json.dump(val, handle)

    if callback:
        callback("All done!")


if __name__ == '__main__':
    cb = Callback()
    args = parse_args()
    handle = lzma.open(args.infile, 'rb')
    analyze_feed(handle, args, callback=cb.callback)
