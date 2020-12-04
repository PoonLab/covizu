import random
import json

from covizu.utils import db_utils, seq_utils
from covizu.utils.progress_utils import Callback

import argparse
import tempfile
import subprocess

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from io import StringIO
from multiprocessing import Pool

import sys
import os

from covizu.utils.seq_utils import load_vcf, filter_problematic


sys.setrecursionlimit(20000)  # fix for issue #127, default limit 1000


def import_json(path, vcf_file, max_missing=600, callback=None):
    """
    Read genome features (genetic differences from reference) from JSON file.
    For command line execution with JSON file input.

    :param path:  str, relative or absolute path to JSON input file
    :param max_missing:  int, maximum tolerated number of uncalled bases
    :return:  list, [qname, diffs, missing]
    """
    with open(path) as fp:
        obj = json.load(fp)

    # remove features known to be problematic - note entries are dicts
    mask = load_vcf(vcf_file)
    filtered = filter_problematic(obj, mask=mask, callback=callback)

    # remove genomes with too many uncalled bases
    count = len(filtered)
    features = [row for row in filtered if seq_utils.total_missing(row) < max_missing]

    if callback:
        callback("dropped {} records with uncalled bases in excess "
                 "of {}".format(count - len(features), max_missing))
    return features


def split_by_lineage(features, lineages):
    """
    Partition feature list by Pangolin lineage assignments
    :param features:  list, feature vectors by GISAID record
    :param lineages:  dict, lineage assignment keyed by accession, from dump_lineages()
    :return:  dict, feature lists keyed by lineage
    """
    result = {}
    for row in features:
        qname = row[0]
        accn = qname.split('|')[1]
        val = lineages.get(accn, None)
        if val is None:
            print("Error in clustering::split_by_lineage(), no lineage assignment"
                  " for accession {}".format(accn))
            continue
        lineage = val['lineage']

        if lineage not in result:
            result.update({lineage: []})
        result[lineage].append(row)
    return result


def sample_with_replacement(population, k):
    """ Sample <k> items from iterable <population> with replacement """
    n = len(population)
    pop = list(population)
    return [pop[int(n * random.random())] for _ in range(k)]


def build_trees(records, nboot=100, binpath='rapidnj', threads=1, callback=None):
    """
    Recode feature vectors with integer indices based on set union.
    Pass results to bootstrap() to reconstruct trees by neighbor-joining method.

    :param records:  list, dict for each record
    :param nboot:  int, number of bootstrap samples
    :param binpath:  str, path to RapidNJ binary executable file
    :param threads:  int, number of threads for bootstrap
    :param callback:  optional, function for progress monitoring
    :return:  list, list; Phylo.BaseTree objects and labels associated with tips
    """
    # compress genomes with identical feature vectors
    fvecs = {}
    for record in records:
        label = '|'.join([record['covv_virus_name'], record['covv_accession_id'],
                          record['covv_collection_date']])
        key = tuple(record['diffs'])
        if key not in fvecs:
            fvecs.update({key: []})
        fvecs[key].append(label)

    # generate union of all features
    if callback:
        callback("Reduced to {} variants; generating feature set union".format(len(fvecs)))
    union = {}
    labels = []
    indexed = []
    for fvec in fvecs:
        labels.append(fvecs[fvec])
        for feat in fvec:
            if feat not in union:
                union.update({feat: len(union)})
        # recode feature vectors as integers (0-index)
        indexed.append(set(union[feat] for feat in fvec))

    # generate distance matrices from bootstrap samples
    if callback:
        callback("Reconstructing trees ({} threads)".format(threads))
    if threads == 1:
        trees = [bootstrap(union, indexed, binpath) for _ in range(nboot)]
    else:
        with Pool(threads) as pool:
            results = [pool.apply_async(bootstrap, [union, indexed, binpath]) for _ in range(nboot)]
            trees = [r.get() for r in results]

    return trees, labels


def bootstrap(union, indexed, binpath='rapidnj', callback=None):
    """
    Sample features from set union at random with replacement.  We use the
    result to weight the symmetric differences when calculating pairwise
    distances.  Pass the resulting distance matrix to RapidNJ to reconstruct
    a tree.

    :param union:  set, all observed genetic differences from reference (features)
    :param indexed:  list, feature vectors encoded as integers
    :param binpath:  str, path to RapidNJ binary executable
    :param callback:  function, optional for progress monitoring
    """
    sample = [int(len(union) * random.random()) for _ in range(len(union))]
    weights = dict([(y, sample.count(y)) for y in sample])
    n = len(indexed)

    outfile = tempfile.NamedTemporaryFile('w', prefix="cvz_boot_")
    outfile.write('{0:>5}\n'.format(n))
    for i in range(n):
        outfile.write('{}'.format(i))
        for j in range(n):
            if i == j:
                outfile.write(' {0:>2}'.format(0))
            else:
                # symmetric difference
                sd = tuple(indexed[i] ^ indexed[j])
                d = sum(weights.get(y, 0) for y in sd)
                outfile.write(' {0:>2}'.format(d))
        outfile.write('\n')
    outfile.flush()

    # call RapidNJ on temp file
    cmd = [binpath, outfile.name, '-i', 'pd', '--no-negative-length']
    stdout = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
    handle = StringIO(stdout.decode('utf-8'))
    outfile.close()

    return Phylo.read(handle, 'newick')


def label_nodes(tree, tip_index):
    """
    Attempt a faster annotation of internal nodes with tip labels by
    preorder traversal.  See issue #158.
    :param tree:  Phylo.BaseTree object
    :param tip_index:  dict, integer indices keyed by tip name
    :return:  Phylo.BaseTree annotated with `tip_index` attribute
    """
    for node in tree.get_nonterminals(order='postorder'):
        tips = []
        for child in node.clades:
            if child.is_terminal():
                tips.append(tip_index[child.name])
            else:
                tips.extend(child.tip_index)
        tips.sort()
        node.tip_index = tips
    return tree


def consensus(trees, cutoff=0.5, callback=None):
    """
    Generate a consensus tree by counting splits and using the splits with
    frequencies above the cutoff to resolve a star tree.
    :param trees:  iterable containing Phylo.BaseTree objects
    :param cutoff:  float, bootstrap threshold (default 0.5)
    :return:  Phylo.BaseTree
    """
    if type(trees) is not list:
        # resolve generator object
        trees = list(trees)

    count = len(trees)

    # store terminal labels and branch lengths
    tip_index = {}
    for i, tip in enumerate(trees[0].get_terminals()):
        tip_index.update({tip.name: i})

    if callback:
        callback("Recording splits and branch lengths")
    splits = {}
    terminals = dict([(tn, []) for tn in tip_index.keys()])
    for phy in trees:
        # record terminal branch lengths
        for tip in phy.get_terminals():
            terminals[tip.name].append(tip.branch_length)

        # record splits in tree
        phy = label_nodes(phy, tip_index)
        for node in phy.get_nonterminals():
            key = tuple(node.tip_index)
            if key not in splits:
                splits.update({key: []})
            splits[key].append(node.branch_length)

    # filter splits by frequency threshold
    intermed = [(len(k), k, v) for k, v in splits.items() if len(v)/count >= cutoff]
    intermed.sort()

    # construct consensus tree
    if callback:
        callback("Building consensus tree")
    orphans = dict(
        [(tip_index[tname], Clade(name=tname, branch_length=sum(tdata)/len(tdata)))
        for tname, tdata in terminals.items()]
    )

    for _, key, val in intermed:
        # average branch lengths across relevant trees
        if all([v is None for v in splits[key]]):
            bl = None
        else:
            bl = sum(splits[key]) / len(splits[key])
        support = len(val) / count
        node = Clade(branch_length=bl, confidence=support)

        for child in key:
            branch = orphans.pop(child, None)
            if branch:
                node.clades.append(branch)

        # use a single tip name to label ancestral node
        newkey = tip_index[node.get_terminals()[0].name]
        orphans.update({newkey: node})

    return orphans.popitem()[1]


def parse_args():
    """ Command-line interface """
    parser = argparse.ArgumentParser(
        description="Partition genomes by Pangolin lineage, generate distance"
                    "matrices and use neighbor-joining."
    )
    parser.add_argument("json", type=str,
                        help="input, JSON file generated by minimap2.py")
    parser.add_argument("-o", "--outdir", type=str, default=os.getcwd(),
                        help="output, directory to export trees (one file per lineage). "
                             "Defaults to current working directory.")
    parser.add_argument("--cons", action='store_true',
                        help="Generate consensus trees instead of exporting "
                             "bootstrap trees to files.")
    parser.add_argument("--db", type=str, default="data/gsaid.db",
                        help="Path to sqlite3 database")
    parser.add_argument("--vcf", type=str, default="data/problematic_sites_sarsCov2.vcf",
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads to run, default 1.")
    parser.add_argument("--maxsize", type=int, default=1000,
                        help="Write data to disk for lineages above this "
                             "threshold; otherwise work in RAM.  Override "
                             "with `--threads 1`.")
    parser.add_argument("--cutoff", type=float, default=0.5,
                        help="Bootstrap cutoff for consensus tree (default 0.5). "
                             "Only used if --cons is specified.")
    return parser.parse_args()


if __name__ == "__main__":
    # command-line execution
    args = parse_args()
    cb = Callback()

    cb.callback('loading lineage classifications from database')
    lineages = db_utils.dump_lineages(args.db)

    cb.callback('loading JSON')
    features = import_json(args.json, vcf_file=args.vcf, callback=cb.callback)

    by_lineage = split_by_lineage(features, lineages)
    for lineage, lfeatures in by_lineage.items():
        cb.callback('start {}, {} entries'.format(lineage, len(lfeatures)))

        # calculate symmetric difference matrix and run NJ on bootstrap samples
        filtered = seq_utils.filter_outliers(lfeatures)
        trees, labels = build_trees(filtered, nboot=args.nboot, callback=cb.callback)

        # export labels as CSV
        outfile = os.path.join(args.outdir, "{}.labels.csv".format(lineage))
        with open(outfile, 'w') as handle:
            handle.write("name,index\n")
            for i, names in enumerate(labels):
                for nm in names:
                    handle.write('{},{}\n'.format(nm, i))

        # export trees
        outfile = os.path.join(
            args.outdir, '{}.{}.nwk'.format(lineage, 'cons' if args.cons else 'boot')
        )
        if trees is None:
            # lineage only has one variant, no meaningful tree
            with open(outfile, 'w') as handle:
                handle.write('({}:0);\n'.format(labels[0][0]))
            continue

        if args.cons:
            trees = consensus(trees, cutoff=args.cutoff)
            cb.callback('built consensus tree')

        Phylo.write(trees, file=outfile, format='newick')
