import argparse
import os
from datetime import datetime, date
import lzma
import gzip
import covizu
from csv import DictReader
import itertools

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
             'lineage': None}

op_fields = {'name': 'strain',
             'accession': "genbank_accession",
             'coldate': "date",
             'region': "region",
             'country': "country",
             'division': "division",
             'lineage': "pango_lineage"}


def import_metadata(meta_file, fields, lineages=None, region=None, callback=None):
    """
    Import metadata from file and append Pangolin lineage classifications
    if necessary.

    :param meta_file:  str, path to gz-compressed metadata file
    :param fields:  dict, map fieldnames to our keys
    :param lineages:  dict, Pangolin classifications
    :param region:  str, default value if not present in metadata
    :param callback:  Callback function

    :yield:  dict, metadata keyed by sequence label
    """
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
        coldate = row[fields['coldate']]

        if lineages:
            lineage = lineages.get(label.upper(), None)
            if lineage is None:
                if callback:
                    callback("Failed to retrieve lineage assignment for {}".format(label), level='ERROR')
                sys.exit()
        else:
            lineage = row[fields['lineage']]

        metadata.update({label: {
            'accession': row[fields['accession']],
            'coldate': coldate,
            'region': row[fields['region']] if region is None else region,
            'country': row[fields['country']],
            'division': row[fields['division']],
            'lineage': lineage
        }})

    return metadata


def merge_data(fasta_file, metadata, minlen=29000, mindate=date(2019, 12, 1),
               callback=None, limit=None):
    """
    :param fasta_file:  str, path to xz-compressed FASTA file
    :param metadata:  dict, returned from import_metadata()
    :param minlen:  int, minimum sequence length
    :param mindate:  datetime.date object, earliest acceptable sample
                     collection date
    :param callback:  optional, progress_utils.Callback object
    :param limit:  int, stop iteration at count (DEBUGGING)
    :yield:  dict, record including label, sequence and metadata
    """
    # parse FASTA and do some basic QC
    handle = lzma.open(fasta_file, 'rt')
    count = 0
    for label, sequence in seq_utils.iter_fasta(handle):
        count += 1
        if limit and count > limit:
            break
        
        if len(sequence) < minlen:
            if callback:
                callback("Rejected short sequence: {}".format(label), level='WARN')
            continue

        record = {'label': label, 'sequence': sequence}
        if label in metadata:
            # validate collection date
            dt = metadata[label]['coldate']
            coldate = seq_utils.fromisoformat(dt)
            if coldate is None or coldate < mindate or coldate > date.today():
                if callback:
                    callback("Rejected record with bad date: {}".format(label), level="WARN")
                continue
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

    parser.add_argument("--vsfasta", required=True, type=str,
                        help="input, path to VirusSeq xz-compressed FASTA")
    parser.add_argument("--vsmeta", required=True, type=str,
                        help="input, path to VirusSeq gzip metadata TSV")
    parser.add_argument("--vspango", type=str, help="input, path to PANGO lineage classifications from the Viral AI database")

    parser.add_argument("--opfasta", required=True, type=str,
                        help="input, path to Nextstrain xz-compresed FASTA")
    parser.add_argument("--opmeta", required=True, type=str,
                        help="input, path to Nextstrain gzip metadata TSV")

    parser.add_argument("--outdir", type=str, default='data/',
                        help="option, path to write output files")
    parser.add_argument("--limit", type=int, help="option, stop at number for debugging")

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
                        default=os.path.join(covizu.__path__[0], "data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf"),
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument('--ft2bin', default='fasttree2',
                        help='path to fasttree2 binary executable')

    parser.add_argument('--lineages', type=str,
                        default=os.path.join(covizu.__path__[0], "data/pango-designation/lineages.csv"),
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

    # check that the user has included submodules
    if (not os.path.exists(os.path.join(covizu.__path__[0], "data/pango-designation/lineages.csv")) or 
            not os.path.exists(os.path.join(covizu.__path__[0], "data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf"))):
        try:
            subprocess.check_call("git submodule init; git submodule update", shell=True)
        except:
            cb.callback("Error adding the required submodules")
            sys.exit()

    # update submodules
    try:
        subprocess.check_call("git submodule foreach git pull origin master", shell=True)
    except:
        cb.callback("Could not update submodules", level='ERROR')

    lineages = {}
    handle = open(args.vspango)
    for row in DictReader(handle):
        lineages.update({row['isolate'].upper() : row['lineage']})

    cb.callback("importing VirusSeq metadata")
    metadata = import_metadata(meta_file=args.vsmeta, fields=vs_fields, lineages=lineages,
                               region='North America', callback=cb.callback)
    cb.callback("importing VirusSeq sequences")
    virusseq = merge_data(fasta_file=args.vsfasta, metadata=metadata,
                          callback=cb.callback, limit=args.limit)

    cb.callback("importing OpenData metadata")
    metadata = import_metadata(meta_file=args.opmeta, fields=op_fields,
                               callback=cb.callback)
    cb.callback("importing OpenData sequences")
    opendata = merge_data(fasta_file=args.opfasta, metadata=metadata,
                          callback=cb.callback, limit=args.limit)

    # run analysis
    cb.callback("starting main analysis")
    feed = itertools.chain(virusseq, opendata)
    analyze_feed(feed, args=args, callback=cb.callback)
