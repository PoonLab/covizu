import subprocess
from Bio import Phylo
from covizu import clustering, treetime, beadplot
import sys


def build_timetree(by_lineage, args, callback=None):
    """ Generate time-scaled tree of Pangolin lineages """
    fasta = treetime.retrieve_genomes(by_lineage, ref_file=args.ref, earliest=True)

    if callback:
        callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = treetime.fasttree(fasta, binpath=args.ft2bin)

    if callback:
        callback("Reconstructing time-scaled tree with {}".format(args.ttbin))
    nexus_file = treetime.treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                                   clock=args.clock, verbosity=0)

    # writes output to treetime.nwk at `nexus_file` path
    return treetime.parse_nexus(nexus_file, fasta)


def beadplot_serial(lineage, features, args, callback=None):
    """ Compute distance matrices and reconstruct NJ trees """
    # bootstrap sampling and NJ tree reconstruction, serial mode
    trees, labels = clustering.build_trees(features, args, callback=callback)
    if trees is None:
        # lineage only has one variant, no meaningful tree
        beaddict = {'lineage': lineage, 'nodes': {}, 'edges': []}

        # use earliest sample as variant label
        intermed = [(label['coldate'], label) for label in labels[0]]
        intermed.sort()
        variant = intermed[0]['accession']
        beaddict['nodes'].update({variant: []})

        for items in intermed:
            beaddict['nodes'][variant].append(items)
        return beaddict

    # generate majority consensus tree
    ctree = clustering.consensus(iter(trees), cutoff=args.boot_cutoff)

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


def make_beadplots(by_lineage, args, callback=None, t0=None, txtfile='minor_lineages.txt'):
    """
    Wrapper for beadplot_serial - divert to clustering.py in MPI mode if
    lineage has too many genomes.
    :param by_lineage:  dict, feature vectors stratified by lineage
    :param args:  Namespace, from argparse.ArgumentParser()
    :param t0:  float, datetime.timestamp.
    :param txtfile:  str, path to file to write minor lineage names
    :return:  list, beadplot data by lineage
    """
    # partition lineages into major and minor categories
    intermed = [(len(features), lineage) for lineage, features in by_lineage.items()
                if len(features) < args.mincount]
    intermed.sort(reverse=True)  # descending order
    minor = dict([(lineage, None) for _, lineage in intermed if lineage is not None])

    # export minor lineages to text file
    with open(txtfile, 'w') as handle:
        for lineage in minor:
            handle.write('{}\n'.format(lineage))

    # launch MPI job across minor lineages
    if callback:
        callback("start MPI on minor lineages")
    mpirun = ["mpirun"]
    if args.machine_file:
        mpirun.extend(["--machinefile", args.machine_file])
    elif args.np:
        mpirun.extend(["-np", str(args.np)])
    else:
        mpirun = [""]  # serial mode
        sys.exit()

    cmd = ["python3", "covizu/clustering.py",
           args.bylineage, txtfile,  # positional arguments <JSON file>, <str>
           "--mode", "flat",
           "--max-variants", str(args.max_variants),
           "--nboot", str(args.nboot), "--outdir", "data"]
    if t0:
        cmd.extend(["--timestamp", str(t0)])

    subprocess.check_call(mpirun+cmd)

    # process major lineages
    for lineage, features in by_lineage.items():
        if lineage in minor:
            continue
        if callback:
            callback('start {}, {} entries'.format(lineage, len(features)))
            # call out to MPI
            cmd = [
                "python3", "covizu/clustering.py",
                args.bylineage, lineage,  # positional arguments <JSON file>, <str>
                "--mode", "deep",
                "--max-variants", str(args.max_variants),
                "--nboot", str(args.nboot), "--outdir", "data", "--binpath", args.binpath
            ]
            if t0:
                cmd.extend(["--timestamp", str(t0)])

            subprocess.check_call(mpirun+cmd)

    # parse output files
    if callback:
        callback("Parsing output files")
    result = []
    for lineage in by_lineage:
        # import trees
        lineage_name = lineage.replace('/', '_')  # issue #297
        outfile = open('data/{}.nwk'.format(lineage_name))
        trees = Phylo.parse(outfile, 'newick')  # note this returns a generator

        # import label map
        with open('data/{}.labels.csv'.format(lineage_name)) as handle:
            label_dict = import_labels(handle)

        if len(label_dict) == 1:
            # handle case of only one variant
            # lineage only has one variant, no meaningful tree
            beaddict = {'lineage': lineage, 'nodes': {}, 'edges': []}

            # use earliest sample as variant label
            intermed = [label.split('|')[::-1] for label in label_dict['0']]
            intermed.sort()
            variant = intermed[0][3]  # coldate, country, region, _ACCESSION_, label
            beaddict['nodes'].update({variant: []})

            for values in intermed:
                beaddict['nodes'][variant].append([values])
        else:
            # generate beadplot data
            ctree = clustering.consensus(trees, cutoff=args.boot_cutoff, callback=callback)
            outfile.close()  # done with Phylo.parse generator

            ctree = beadplot.annotate_tree(ctree, label_dict)
            beaddict = beadplot.serialize_tree(ctree)
            beaddict.update({'sampled_variants': len(label_dict)})

            beaddict.update({'lineage': lineage})
        result.append(beaddict)

    return result


def get_mutations(by_lineage):
    """
    Extract common mutations from feature vectors for each lineage
    :param by_lineage:  dict, return value from process_feed()
    :return:  dict, common mutations by lineage
    """
    result = {}
    for lineage, samples in by_lineage.items():
        # enumerate features
        counts = {}
        for sample in samples:
            for diff in sample['diffs']:
                feat = tuple(diff)
                if feat not in counts:
                    counts.update({feat: 0})
                counts[feat] += 1
        # filter for mutations that occur in at least half of samples
        common = [feat for feat, count in counts.items() if count/len(samples) >= 0.5]
        result.update({lineage: common})
    return result


def sort_by_lineage(records, callback=None):
    """
    Resolve stream into a dictionary keyed by Pangolin lineage
    :param records:  generator, return value of extract_features()
    :return:  dict, lists of records keyed by lineage
    """
    result = {}
    for i, record in enumerate(records):
        if callback and i % 1000 == 0:
            callback('aligned {} records'.format(i))
        lineage = record['lineage']
        if lineage not in result:
            result.update({lineage: []})
        result[lineage].append(record)
    return result
