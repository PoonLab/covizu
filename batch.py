import argparse
import os
import sys
import json
import subprocess
from Bio import Phylo
from datetime import datetime

import covizu
from covizu import clustering, treetime, beadplot
from covizu.utils import gisaid_utils
from covizu.utils.progress_utils import Callback


def parse_args():
    parser = argparse.ArgumentParser(description="CoVizu analysis pipeline automation")

    parser.add_argument('--url', type=str, default=os.environ["GISAID_URL"],
                        help="URL to download provision file, defaults to environment variable.")
    parser.add_argument('--user', type=str, default=os.environ["GISAID_USER"],
                        help="GISAID username, defaults to environment variable.")
    parser.add_argument('--password', type=str, default=os.environ["GISAID_PSWD"],
                        help="GISAID password, defaults to environment variable.")

    parser.add_argument("--infile", type=str, default=None,
                        help="input, path to xz-compressed JSON; if not specified, "
                             "download xz file from GISAID provision feed.")
    parser.add_argument("--outfile", type=argparse.FileType('w'),
                        default='data/clusters.{}.json'.format(datetime.now().isoformat().split('.')[0]),
                        help="output, dest for JSON beadplot file")

    parser.add_argument("--bylineage", type=str, default='data/by_lineage.json',
                        help="option, path to write JSON of features by lineage")

    parser.add_argument('--minlen', type=int, default=29000, help='option, minimum genome length (nt)')
    parser.add_argument('--mindate', type=str, default='2019-12-01', 
                        help='option, earliest possible sample collection date (ISO format, default '
                              '2019-12-01')
   
    parser.add_argument('--batchsize', type=int, default=500,
                        help='option, number of records to batch process with minimap2')

    parser.add_argument("--ref", type=str,
                        default=os.path.join(covizu.__path__[0], "data/NC_045512.fa"),
                        help="option, path to FASTA file with reference genome")
    parser.add_argument('--mmbin', type=str, default='minimap2',
                        help="option, path to minimap2 binary executable")
    parser.add_argument('-mmt', "--mmthreads", type=int, default=8,
                        help="option, number of threads for minimap2.")

    parser.add_argument('--misstol', type=int, default=300,
                        help="option, maximum tolerated number of missing bases per "
                             "genome (default 300).")
    parser.add_argument("--vcf", type=str,
                        default=os.path.join(covizu.__path__[0], "data/problematic_sites_sarsCov2.vcf"),
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument('--ft2bin', default='fasttree2',
                        help='option, path to fasttree2 binary executable')

    parser.add_argument('--outdir', default='data/',
                        help='option, directory to write TreeTime output files')
    parser.add_argument('--ttbin', default='treetime',
                        help='option, path to treetime binary executable')
    parser.add_argument('--clock', type=float, default=8e-4,
                        help='option, specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')

    parser.add_argument('--datetol', type=float, default=0.1,
                        help='option, exclude tips from time-scaled tree '
                             'with high discordance between estimated and '
                             'known sample collection dates (year units,'
                             'default: 0.1)')

    parser.add_argument('--binpath', type=str, default='rapidnj',
                        help='option, path to RapidNJ binary executable')
    parser.add_argument('--mincount', type=int, default=500,
                        help='option, minimum number of variants in lineage '
                             'above which MPI processing will be used.')
    parser.add_argument('--machine_file', type=str, default='mfile',
                        help='option, path to machine file for MPI.')
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")

    parser.add_argument("--cutoff", type=float, default=0.5,
                        help="Bootstrap cutoff for consensus tree (default 0.5). "
                             "Only used if --cons is specified.")

    return parser.parse_args()


def process_feed(args, callback=None):
    """ Process feed data """
    if callback:
        callback("Processing GISAID feed data")
    loader = gisaid_utils.load_gisaid(args.infile, minlen=args.minlen, mindate=args.mindate)
    batcher = gisaid_utils.batch_fasta(loader, size=args.batchsize)
    aligned = gisaid_utils.extract_features(batcher, ref_file=args.ref, binpath=args.mmbin,
                                            nthread=args.mmthreads, minlen=args.minlen)
    filtered = gisaid_utils.filter_problematic(aligned, vcf_file=args.vcf, callback=callback)
    return gisaid_utils.sort_by_lineage(filtered, callback=callback)


def build_timetree(by_lineage, args, callback=None):
    """ Generate time-scaled tree of Pangolin lineages """
    fasta = treetime.retrieve_genomes(by_lineage, ref_file=args.ref)

    if callback:
        callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = treetime.fasttree(fasta, binpath=args.ft2bin)

    if callback:
        callback("Reconstructing time-scaled tree with {}".format(args.ttbin))
    nexus_file = treetime.treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                                   clock=args.clock, verbosity=0)

    # writes output to treetime.nwk at `nexus_file` path
    return treetime.parse_nexus(nexus_file, fasta, date_tol=args.datetol)


def beadplot_serial(lineage, features, args, callback=None):
    """ Compute distance matrices and reconstruct NJ trees """
    # bootstrap sampling and NJ tree reconstruction, serial mode
    trees, labels = clustering.build_trees(features, args, callback=callback)
    if trees is None:
        # lineage only has one variant, no meaningful tree
        beaddict = {'lineage': lineage, 'nodes': {}, 'edges': []}

        # use earliest sample as variant label
        intermed = [label.split('|')[::-1] for label in labels[0]]
        intermed.sort()
        variant = intermed[0][1]
        beaddict['nodes'].update({variant: []})

        for coldate, accn, label1 in intermed:
            beaddict['nodes'][variant].append({
                'accession': accn, 'label1': label1, 'country': label1.split('/')[1],
                'coldate': coldate
            })
        return beaddict

    # generate majority consensus tree
    ctree = clustering.consensus(iter(trees), cutoff=args.cutoff)

    # collapse polytomies and label internal nodes
    label_dict = dict([(str(idx), lst) for idx, lst in enumerate(labels)])
    ctree = beadplot.annotate_tree(ctree, label_dict, callback=callback)

    # convert to JSON format
    beaddict = beadplot.serialize_tree(ctree)
    beaddict.update({'lineage': lineage})
    return beaddict


def import_labels(handle, callback=None):
    """ Load map of genome labels to tip indices from CSV file """
    result = {}
    _ = next(handle)  # skip header line
    for line in handle:
        try:
            qname, idx = line.strip('\n').split(',')
        except ValueError:
            if callback:
                callback("import_labels() failed to parse line {}".format(line), level="ERROR")
            raise  # issue #206, sequence label contains delimiter

        if idx not in result:
            result.update({idx: []})
        result[idx].append(qname)
    return result


def make_beadplots(by_lineage, args, callback=None, t0=None):
    """
    Wrapper for beadplot_serial - divert to clustering.py in MPI mode if
    lineage has too many genomes.

    :param by_lineage:  dict, feature vectors stratified by lineage
    :param args:  Namespace, from argparse.ArgumentParser()
    :param t0:  float, datetime.timestamp.
    :return:  list, beadplot data by lineage
    """
    result = []
    for lineage, features in by_lineage.items():
        if callback:
            callback('start {}, {} entries'.format(lineage, len(features)))

        if len(features) < args.mincount:
            # serial processing
            if len(features) == 0:
                continue  # empty lineage, skip (should never happen)
            beaddict = beadplot_serial(lineage, features, args)
        else:
            # call out to MPI
            cmd = [
                "mpirun", "--machinefile", args.machine_file,
                "python3", "covizu/clustering.py",
                 args.bylineage, lineage,  # positional arguments <JSON file>, <str>
                 "--nboot", str(args.nboot), "--outdir", "data"
            ]
            if t0:
                cmd.extend(["--timestamp", str(t0)])
            subprocess.check_call(cmd)

            # import trees
            outfile = open('data/{}.nwk'.format(lineage))
            trees = Phylo.parse(outfile, 'newick')  # note this returns a generator

            # import label map
            with open('data/{}.labels.csv'.format(lineage)) as handle:
                label_dict = import_labels(handle)

            # generate beadplot data
            ctree = clustering.consensus(trees, cutoff=args.cutoff, callback=callback)
            outfile.close()  # done with Phylo.parse generator

            ctree = beadplot.annotate_tree(ctree, label_dict)
            beaddict = beadplot.serialize_tree(ctree)

        beaddict.update({'lineage': lineage})
        result.append(beaddict)

    return result


if __name__ == "__main__":
    args = parse_args()
    cb = Callback()

    # check that user has loaded openmpi module
    try:
        subprocess.check_call(['mpirun', '-np', '2', 'ls'], stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        cb.callback("mpirun not loaded - run `module load openmpi/gnu`", level='ERROR')
        sys.exit()

    # download xz file if not specified by user
    if args.infile is None:
        cb.callback("No input specified, downloading data from GISAID feed...")
        args.infile = gisaid_utils.download_feed(args.url, args.user, args.password)

    by_lineage = process_feed(args, cb.callback)
    with open(args.bylineage, 'w') as handle:
        # export to file to process large lineages with MPI
        json.dump(by_lineage, handle)

    timetree = build_timetree(by_lineage, args, cb.callback)

    # FIXME: this is fragile, fails if user specifies custom output file name
    head, tail = os.path.split(args.outfile.name)
    timestamp = tail.split('.')[1]
    nwk_file = os.path.join(head, 'timetree.{}.nwk'.format(timestamp))
    with open(nwk_file, 'w') as handle:
        Phylo.write(timetree, file=handle, format='newick')

    result = make_beadplots(by_lineage, args, cb.callback, t0=cb.t0.timestamp())
    args.outfile.write(json.dumps(result))  # serialize results to JSON

    # write data stats
    dbstat_file = os.path.join(head, 'dbstats.{}.json'.format(timestamp))
    with open(dbstat_file, 'w') as handle:
        nseqs = sum([len(rows) for rows in by_lineage.values()])
        val = {'lastupdate': timestamp.split('T')[0], 'noseqs': nseqs}
        json.dump(val, handle)

    # transfer output files to webserver
    server_root = 'filogeneti.ca:/var/www/html/covizu/data'
    subprocess.check_call(['scp', nwk_file, '{}/timetree.nwk'.format(server_root)])
    subprocess.check_call(['scp', args.outfile.name, '{}/clusters.json'.format(server_root)])
    subprocess.check_call(['scp', dbstat_file, '{}/dbstats.json'.format(server_root)])

    cb.callback("All done!")
