import argparse
import os
import sys
import json
from datetime import datetime, date
from csv import DictReader

import covizu
from covizu.utils import seq_utils, gisaid_utils
from covizu.utils.progress_utils import Callback
from covizu.utils.batch_utils import *
from covizu.utils.seq_utils import SC2Locator


def parse_args():
    parser = argparse.ArgumentParser(
        description="CoVizu analysis pipeline automation for execution on local files"
    )

    parser.add_argument("infile", type=str,
                        help="input, path to xz-compressed JSON; if not specified, "
                             "download xz file from GISAID provision feed.")
    parser.add_argument("pangolineages", type=argparse.FileType('r'),
                        help="input, CSV output generated by Pangolin")
    parser.add_argument("--outdir", type=str, default='data/',
                        help="option, path to write output files")

    parser.add_argument("--bylineage", type=str, default='data/by_lineage.json',
                        help="path to write JSON of features by lineage")

    parser.add_argument('--lineages', type=str,
                        default=os.path.join(covizu.__path__[0], "data/pango-designation/lineages.csv"),
                        help="optional, path to CSV file containing Pango lineage designations.")

    parser.add_argument('--minlen', type=int, default=29000, help='minimum genome length (nt)')
    parser.add_argument('--mindate', type=str, default='2019-12-01',
                        help='earliest possible sample collection date (ISO format, default '
                             '2019-12-01')
    parser.add_argument('--poisson-cutoff', type=float, default=0.001,
                        help='filtering outlying genomes whose distance exceeds the upper '
                             'quantile of Poisson distribution (molecular clock).  Default 0.001 '
                             'corresponds to 99.9%% cutoff.')

    parser.add_argument('--batchsize', type=int, default=500,
                        help='number of records to batch process with minimap2')
    parser.add_argument('--max-variants', type=int, default=5000,
                        help='option, limit number of variants per lineage (default 5000)')

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
                        default=os.path.join(covizu.__path__[0], "data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf"),
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument('--ft2bin', default='fasttree2',
                        help='path to fasttree2 binary executable')

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


def stream_local(path, lineage_file, minlen=29000, mindate='2019-12-01', callback=None):
    """ Convert local FASTA file to feed-like object - replaces load_gisaid() """
    mindate = seq_utils.fromisoformat(mindate)

    # parse CSV output from Pangolin
    reader = DictReader(lineage_file)
    if reader.fieldnames != ['taxon', 'lineage', 'conflict', 'ambiguity_score', 'scorpio_call', 'scorpio_support', 'scorpio_conflict', 'scorpio_notes', 'version', 'pangolin_version', 'scorpio_version', 'constellation_version', 'is_designated', 'qc_status', 'qc_notes', 'note']:
        if callback:
            callback("Lineage CSV header does not match expected.", level='ERROR')
        sys.exit()

    lineages = {}
    for row in reader:
        lineages.update({row['taxon']: row['lineage']})

    handle = open(path)
    rejects = {'short': 0, 'baddate': 0, 'nonhuman': 0}
    for header, seq in seq_utils.iter_fasta(handle):
        if len(seq) < minlen:
            rejects['short'] += 1
            continue  # sequence is too short

        # hCoV-19/Canada/Qc-L00240569/2020|EPI_ISL_465679|2020-03-27
        label, accn, coldate = header.split('|')
        country = label.split('/')[1]
        if country == '' or country[0].islower():
            rejects['nonhuman'] += 1
            continue

        if coldate.count('-') != 2:
            rejects['baddate'] += 1
            continue  # incomplete collection date
        dt = seq_utils.fromisoformat(coldate)
        if dt < mindate or dt > date.today():
            rejects['baddate'] += 1
            continue  # reject records with non-sensical collection date

        lineage = lineages.get(header, None)
        if lineage is None:
            if callback:
                callback(
                    "Failed to retrieve lineage assignment for {}".format(header),
                    level='ERROR'
                )
            sys.exit()

        record = {
            'covv_virus_name': label,
            'covv_accession_id': accn,
            'sequence': seq,
            'covv_collection_date': coldate,
            'covv_lineage': lineage,
            'covv_location': country
        }
        yield record

    if callback:
        callback("Rejected {short} short genomes\n         {baddate} records with bad "
                 "dates\n         {nonhuman} non-human genomes".format(**rejects))


def process_local(args, callback=None):
    """ Analyze genome sequences from local FASTA file """
    with open(args.ref) as handle:
        reflen = len(seq_utils.convert_fasta(handle)[0][1])

    loader = stream_local(args.infile, args.pangolineages, minlen=args.minlen,
                          mindate=args.mindate, callback=callback)
    batcher = gisaid_utils.batch_fasta(loader, size=args.batchsize)
    aligned = gisaid_utils.extract_features(batcher, ref_file=args.ref, binpath=args.mmbin,
                                            nthread=args.mmthreads, minlen=args.minlen)
    filtered = gisaid_utils.filter_problematic(aligned, vcf_file=args.vcf, cutoff=args.poisson_cutoff,
                                               callback=callback)
    return gisaid_utils.sort_by_lineage(filtered, callback=callback)


if __name__ == "__main__":
    args = parse_args()
    cb = Callback()

    # check that user has loaded openmpi module
    try:
        subprocess.check_call(['mpirun', '-np', '2', 'ls'], stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        cb.callback("mpirun not loaded - run `module load openmpi/gnu`", level='ERROR')
        sys.exit()

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
        cb.callback("Error updating submodules")
        sys.exit()

    by_lineage = process_local(args, cb.callback)
    with open(args.bylineage, 'w') as handle:
        # export to file to process large lineages with MPI
        json.dump(by_lineage, handle)

    # reconstruct time-scaled tree
    timetree, residuals = build_timetree(by_lineage, args, cb.callback)
    timestamp = datetime.now().isoformat().split('.')[0]
    nwk_file = os.path.join(args.outdir, 'timetree.{}.nwk'.format(timestamp))
    with open(nwk_file, 'w') as handle:
        Phylo.write(timetree, file=handle, format='newick')

    # generate beadplots and serialize to file
    result = make_beadplots(by_lineage, args, cb.callback, t0=cb.t0.timestamp())
    outfile = os.path.join(args.outdir, 'clusters.{}.json'.format(timestamp))
    with open(outfile, 'w') as handle:  # serialize results to JSON
        json.dump(result, fp=handle)

    # get mutation info
    locator = SC2Locator()
    mutations = {}
    for lineage, features in get_mutations(by_lineage).items():
        annots = {locator.parse_mutation(f) : freq for f, freq in features.items()}
        mutations.update({lineage: {a : freq for a, freq in annots.items() if a is not None}})

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
            samples = unpack_records(samples)
            ndiffs = [len(x['diffs']) for x in samples]
            val['lineages'][lineage] = {
                'nsamples': len(samples),
                'lastcoldate': max(x['covv_collection_date'] for x in samples),
                'residual': residuals.get(lineage, 0),
                'max_ndiffs': max(ndiffs),
                'mean_ndiffs': sum(ndiffs)/len(ndiffs),
                'mutations': mutations[lineage]
            }
        json.dump(val, handle)

    cb.callback("All done!")
