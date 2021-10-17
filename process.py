import argparse
import os
import sys
import json
from datetime import datetime
import lzma
import gzip
import covizu
from csv import DictReader
from itertools import chain

from covizu.utils import seq_utils
from covizu.utils.batch_utils import *
from covizu.utils.seq_utils import SC2Locator
from covizu.utils.progress_utils import Callback
from covizu.minimap2 import extract_features


vs_fields = {'name': "fasta header name",
             'accession': "specimen collector sample ID",
             'coldate': "sample collection date",
             'region': None,
             'country': "geo_loc_name (country)",
             'division': "geo_loc_name (state/province/territory)",
             'lineage': None
             }
op_fields = {'name': 'strain',
             'accession': "genbank_accession",
             'coldate': "date",
             'region': "region",
             'country': "country",
             'division': "division",
             'lineage': "pango_lineage"}


def merge_data(fasta_file, meta_file, fields, lineage_file=None, region=None, callback=None):
    """
    Append Pangolin lineage designations to metadata

    :param fasta_file:  str, path to xz-compressed FASTA file
    :param meta_file:  str, path to gz-compressed metadata file
    :param fields:  dict, map fieldnames to our keys
    :param lineage_file:  str, path to Pangolin CSV output; if None, assume lineages
                          are contained in metadata
    :param region:  str, default value if not present in metadata
    :param callback:  Callback function
    :yield:  dict, merged data for each record
    """
    # parse Pangolin output, if provided
    lineages = {}
    if lineage_file:
        reader = DictReader(open(lineage_file))
        if 'taxon' not in reader.fieldnames or 'lineage' not in reader.fieldnames:
            if callback:
                callback("Lineage CSV header does not match expected.", level='ERROR')
            sys.exit()
        for row in reader:
            lineages.update({row['taxon']: row['lineage']})

    # check metadata fieldnames
    handle = gzip.open(meta_file, 'rt')
    reader = DictReader(handle, delimiter='\t')
    for _, field in fields.items():
        if field is not None and field not in reader.fieldnames:
            if callback:
                callback("Missing fieldname {} in VirusSeq metadata".format(field), level='ERROR')
                callback(reader.fieldnames, level='ERROR')
            sys.exit()

    # parse metadata
    metadata = {}
    for row in reader:
        label = row[fields['name']]
        if lineage_file:
            lineage = lineages.get(label, None)
            if lineage is None:
                if callback:
                    callback("Failed to retrieve lineage assignment for {}".format(label), level='ERROR')
                sys.exit()
        else:
            lineage = row[fields['lineage']]

        metadata.update({label: {
            'accession': row[fields['accession']],
            'coldate': row[fields['coldate']],
            'region': row[fields['region']] if region is None else region,
            'country': row[fields['country']],
            'division': row[fields['division']],
            'lineage': lineage
        }})

    # parse FASTA
    handle = lzma.open(fasta_file, 'rt')
    for label, sequence in seq_utils.iter_fasta(handle):
        record = {'label': label, 'sequence': sequence}
        if label in metadata:
            record.update(metadata[label])
        else:
            if callback:
                callback("Failed to retrieve metadata for sequence {}".format(label), level='ERROR')
            sys.exit()
        yield record


def analyze_feed(feed, args, callback=None):
    """
    :param records:  generator returned from merge_data()
    :param args:  Namespace
    :param callback:  optional progress monitoring, see progress_utils.py
    """
    # check that user has loaded openmpi module
    if args.machine_file or args.np:
        try:
            subprocess.check_call(['mpirun', '-np', '2', 'ls'], stdout=subprocess.DEVNULL)
        except FileNotFoundError:
            if callback:
                callback("mpirun not loaded - run `module load openmpi/gnu`", level='ERROR')
            sys.exit()

    # pre-processing
    batcher = seq_utils.batch_fasta(feed, size=args.batchsize)
    aligned = extract_features(batcher, ref_file=args.ref, binpath=args.mmbin,
                               nthread=args.mmthreads, minlen=args.minlen)
    filtered = seq_utils.filter_problematic(aligned, vcf_file=args.vcf, cutoff=args.poisson_cutoff,
                                            callback=callback)
    by_lineage = sort_by_lineage(filtered, callback=callback)

    # reconstruct time-scaled tree
    timetree, residuals = build_timetree(by_lineage, args, callback)
    t0 = datetime.now()
    timestamp = t0.isoformat().split('.')[0]
    nwk_file = os.path.join(args.outdir, 'timetree.{}.nwk'.format(timestamp))
    with open(nwk_file, 'w') as handle:
        Phylo.write(timetree, file=handle, format='newick')

    # generate beadplots and serialize to file
    result = make_beadplots(by_lineage, args, callback, t0=t0.timestamp())
    clust_file = os.path.join(args.outdir, 'clusters.{}.json'.format(timestamp))
    with open(clust_file, 'w') as handle:
        json.dump(result, fp=handle)

    # get mutation info
    locator = SC2Locator()
    mutations = {}
    for lineage, features in get_mutations(by_lineage).items():
        annots = [locator.parse_mutation(f) for f in features]
        mutations.update({lineage: [a for a in annots if a is not None]})

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
                'mean_ndiffs': sum(ndiffs)/len(ndiffs),
                'mutations': mutations[lineage]
            }
        json.dump(val, handle)

    if callback:
        callback("All done!")



def parse_args():
    parser = argparse.ArgumentParser(
        description="CoVizu analysis pipeline automation for execution on local files"
    )

    parser.add_argument("vsfasta", type=str, help="input, path to VirusSeq xz-compressed FASTA")
    parser.add_argument("vspango", type=str, help="input, path to VirusSeq pangolin CSV")
    parser.add_argument("vsmeta", type=str, help="input, path to VirusSeq gzip metadata TSV")

    parser.add_argument("opfasta", type=str, help="input, path to Nextstrain xz-compresed FASTA")
    parser.add_argument("opmeta", type=str, help="input, path to Nextstrain gzip metadata TSV")

    parser.add_argument("--outdir", type=str, default='data/',
                        help="option, path to write output files")

    parser.add_argument('--poisson-cutoff', type=float, default=0.001,
                        help='filtering outlying genomes whose distance exceeds the upper '
                             'quantile of Poisson distribution (molecular clock).  Default 0.001 '
                             'corresponds to 99.9%% cutoff.')
    parser.add_argument('--minlen', type=int, default=29000, help='minimum genome length (nt)')

    parser.add_argument('--batchsize', type=int, default=2000,
                        help='number of records to batch process with minimap2')
    parser.add_argument("--max-variants", type=int, default=5000,
                        help="option, limit number of variants per lineage (default 5,000)")

    parser.add_argument("--ref", type=str,
                        default=os.path.join(covizu.__path__[0], "data/NC_045512.fa"),
                        help="path to FASTA file with reference genome")
    parser.add_argument('--mmbin', type=str, default='minimap2',
                        help="path to minimap2 binary executable")
    parser.add_argument('-mmt', "--mmthreads", type=int, default=16,
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

    parser.add_argument('--lineages', type=str,
                        default=os.path.join(covizu.__path__[0], "data/lineages.csv"),
                        help="optional, path to CSV file containing Pango lineage designations.")
    parser.add_argument('--ttbin', default='treetime',
                        help='path to treetime binary executable')
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')

    parser.add_argument('--datetol', type=float, default=0.1,
                        help='exclude tips from time-scaled tree '
                             'with high discordance between estimated and '
                             'known sample collection dates (year units,'
                             'default: 0.1)')

    parser.add_argument('--binpath', type=str, default='rapidnj',
                        help='path to RapidNJ binary executable')
    parser.add_argument('--mincount', type=int, default=5000,
                        help='minimum number of variants in lineage '
                             'above which MPI processing will be used.')

    # leave both to None for serial execution
    parser.add_argument('--machine_file', type=str, default=None,
                        help='optional, path to machine file for MPI.')
    parser.add_argument('-np', type=int, default=None,
                        help='optional, number of processes for MPI')

    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")

    parser.add_argument("--boot-cutoff", type=float, default=0.5,
                        help="Bootstrap cutoff for consensus tree (default 0.5). "
                             "Only used if --cons is specified.")

    return parser.parse_args()


if __name__ == '__main__':
    cb = Callback()
    args = parse_args()

    # parse VirusSeq
    virusseq = merge_data(fasta_file=args.vsfasta, meta_file=args.vsmeta, fields=vs_fields,
                          lineage_file=args.vspango, region='North America', callback=cb.callback)
    print(next(virusseq))

    # itertools.chain(virusseq, opendata)
    #analyze_feed(args, callback=cb.callback)
