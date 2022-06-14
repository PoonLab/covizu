import subprocess
from Bio import Phylo
from covizu import clustering, beadplot
import sys
import json
import covizu.treetime


def unpack_records(records):
    """
    by_lineage is a nested dict with the inner dicts keyed by serialized
    mutation sets (diffs).  This function is used to reconstitute the mutations
    as a list of tuples, restoring the diffs entry for each record, and to return
    a list of dicts.
    Used in:
    - treetime:retrieve_genomes()
    - get_mutations()
    - batch.py, command line interface
    See issue: https://github.com/PoonLab/covizu/issues/387

    @param records:  dict, sets of genomes under a single lineage classification,
                     each set keyed by their shared set of mutations (diffs)
    @return:  a list of dicts, each dict representing a genome and its metadata
    """
    unpacked = []
    for key, variant in records.items():
        # reconstitute the mutations defining this variant
        diffs = []
        for mutation in key.split(','):
            typ, pos, alt = mutation.split('|')
            if typ == '-':
                alt = int(alt)  # number of nucleotides in indel
            diffs.append(tuple([typ, int(pos), alt]))

        for sample in variant:
            sample.update({'diffs': diffs})
            unpacked.append(sample)
    return unpacked


def build_timetree(by_lineage, args, callback=None):
    """ Generate time-scaled tree of Pangolin lineages """

    if callback:
        callback("Parsing Pango lineage designations")
    handle = open(args.lineages)
    header = next(handle)
    if header != 'taxon,lineage\n':
        if callback:
            callback("Error: {} does not contain expected header row 'taxon,lineage'".format(args.lineages))
        sys.exit()
    lineages = {}
    for line in handle:
        taxon, lineage = line.strip().split(',')
        lineages.update({taxon: lineage})

    # filter recombinants
    by_lineage = {lineage : records for lineage, records in by_lineage.items() if lineage[0] != 'X'}

    if callback:
        callback("Identifying lineage representative genomes")
    fasta = covizu.treetime.retrieve_genomes(by_lineage, known_seqs=lineages, ref_file=args.ref,
                                      earliest=True)

    if callback:
        callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = covizu.treetime.fasttree(fasta, binpath=args.ft2bin)

    if callback:
        callback("Reconstructing time-scaled tree with {}".format(args.ttbin))
    nexus_file = covizu.treetime.treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                                   clock=args.clock, verbosity=0)

    # writes output to treetime.nwk at `nexus_file` path
    return covizu.treetime.parse_nexus(nexus_file, fasta)


def beadplot_serial(lineage, features, args, callback=None):
    """ Compute distance matrices and reconstruct NJ trees """
    # bootstrap sampling and NJ tree reconstruction, serial mode
    trees, labels = clustering.build_trees(features, args, callback=callback)
    if trees is None:
        # lineage only has one variant, no meaningful tree
        beaddict = {'lineage': lineage, 'nodes': {}, 'edges': []}

        # use earliest sample as variant label
        intermed = [label.split('|')[::-1] for label in labels['0']]
        intermed.sort()
        variant = intermed[0][1]
        beaddict.update({'sampled_variants': len(labels)})
        beaddict['nodes'].update({variant: []})

        for coldate, accn, label1 in intermed:
            beaddict['nodes'][variant].append([coldate, accn, label1])
        return beaddict

    # generate majority consensus tree
    ctree = clustering.consensus(iter(trees), cutoff=args.boot_cutoff)

    # collapse polytomies and label internal nodes
    label_dict = dict([(str(idx), lst) for idx, lst in enumerate(labels)])
    ctree = beadplot.annotate_tree(ctree, label_dict, callback=callback)

    # convert to JSON format
    beaddict = beadplot.serialize_tree(ctree)
    beaddict.update({'lineage': lineage})
    beaddict.update({'sampled_variants': len(labels)})
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


def make_beadplots(by_lineage, args, callback=None, t0=None, txtfile='minor_lineages.txt',
                   recode_file="recoded.json"):
    """
    Wrapper for beadplot_serial - divert to clustering.py in MPI mode if
    lineage has too many genomes.

    :param by_lineage:  dict, feature vectors stratified by lineage
    :param args:  Namespace, from argparse.ArgumentParser()
    :param t0:  float, datetime.timestamp.
    :param txtfile:  str, path to file to write minor lineage names
    :param recode_file:  str, path to JSON file to write recoded lineage data

    :return:  list, beadplot data by lineage
    """

    # recode data into variants and serialize
    if callback:
        callback("Recoding features, compressing variants..")
    recoded = {}
    for lineage, records in by_lineage.items():
        union, labels, indexed = clustering.recode_features(records, limit=args.max_variants)

        # serialize tuple keys (features of union), #335
        union = dict([("{0}|{1}|{2}".format(*feat), idx) for feat, idx in union.items()])
        indexed = [list(s) for s in indexed]  # sets cannot be serialized to JSON, #335
        recoded.update({lineage: {'union': union, 'labels': labels,
                                  'indexed': indexed}})

    with open(recode_file, 'w') as handle:
        json.dump(recoded, handle)

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
    cmd = ["mpirun", "--machinefile", args.machine_file, "python3", "covizu/clustering.py",
           recode_file, txtfile,  # positional arguments <JSON file>, <str>
           "--mode", "flat",
           "--max-variants", str(args.max_variants),
           "--nboot", str(args.nboot),
           "--outdir", args.outdir,
           "--binpath", args.binpath  # RapidNJ
           ]
    if t0:
        cmd.extend(["--timestamp", str(t0)])
    subprocess.check_call(cmd)

    # process major lineages
    for lineage, features in by_lineage.items():
        if lineage in minor:
            continue

        if callback:
            callback('start {}, {} entries'.format(lineage, len(features)))

        cmd = [
            "mpirun", "--machinefile", args.machine_file, "python3", "covizu/clustering.py",
            recode_file, lineage,  # positional arguments <JSON file>, <str>
            "--mode", "deep",
            "--max-variants", str(args.max_variants),
            "--nboot", str(args.nboot),
            "--outdir", args.outdir,
            "--binpath", args.binpath
        ]
        if t0:
            cmd.extend(["--timestamp", str(t0)])
        subprocess.check_call(cmd)

    # parse output files
    if callback:
        callback("Parsing output files")
    result = []
    for lineage in recoded:
        # import trees
        lineage_name = lineage.replace('/', '_')  # issue #297
        outfile = open('data/{}.nwk'.format(lineage_name))
        trees = Phylo.parse(outfile, 'newick')  # note this returns a generator

        label_dict = recoded[lineage]['labels']

        if len(label_dict) == 1:
            # handle case of only one variant
            # lineage only has one variant, no meaningful tree
            beaddict = {'nodes': {}, 'edges': []}

            # use earliest sample as variant label
            intermed = [label.split('|')[::-1] for label in label_dict['0']]
            intermed.sort()
            variant = intermed[0][1]
            beaddict['nodes'].update({variant: []})

            for coldate, accn, label1 in intermed:
                beaddict['nodes'][variant].append([coldate, accn, label1])
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
    for lineage, records in by_lineage.items():
        samples = unpack_records(records)

        # enumerate features
        counts = {}
        for sample in samples:
            for diff in sample['diffs']:
                feat = tuple(diff)
                if feat not in counts:
                    counts.update({feat: 0})
                counts[feat] += 1

        # filter for mutations that occur in at least half of samples
        common = dict([(feat, count/len(samples)) for feat, count in counts.items()
                       if count/len(samples) >= 0.5])
        result.update({lineage: common})

    return result

