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


def get_sym_diffs(features, use_file=True):
    """
    Calculate symmetric differences of feature vectors.

    :param features:  list, feature vectors from load_json()/split_by_lineage()
    :param use_file:  bool, if True, write matrix to file
    :return:  dict, symmetric differences indexed by (i,j) tuples OR
              str, path to output file;
              int, number of features in set union
              list, nested list of genome labels by index
    """
    # compress genomes with identical feature vectors
    fvecs = {}
    for qname, diffs, missing in features:
        key = tuple(diffs)
        if key not in fvecs:
            fvecs.update({key: []})
        fvecs[key].append(qname)

    # generate union of all features
    union = {}
    for fvec in fvecs:
        for feat in fvec:
            if feat not in union:
                union.update({feat: len(union)})

    # recode feature vectors as integers
    indexed = []
    labels = []
    for fvec in fvecs:
        indexed.append(set(union[feat] for feat in fvec))
        labels.append(fvecs[fvec])

    # calculate symmetric differences
    n = len(fvecs)
    if use_file:
        # write integer tuples to temporary CSV file
        handle = tempfile.NamedTemporaryFile('w', prefix="cvz_sd_", delete=False)
        for i in range(n):
            for j in range(n):
                sdiff = tuple(indexed[i] ^ indexed[j])
                handle.write(','.join(map(str, sdiff)))
                handle.write('\n')

        handle.close()
        return handle.name, len(union), labels
    else:
        # store tuples in memory
        sym_diffs = {}
        for i in range(n-1):
            for j in range(i+1, n):
                sdiff = tuple(indexed[i] ^ indexed[j])
                sym_diffs[(i, j)] = sdiff
        return sym_diffs, len(union), labels


def bootstrap(sym_diffs, n, m, binpath='rapidnj', callback=None):
    """
    Use nonparametric bootstrap sampling of the feature set union to convert symmetric
    differences to distances by re-weighting features.

    :param sym_diffs:  dict or TemporaryNamedFile; if dict, symmetric difference subsets
                       keyed by (i, j) tuples; else file handle to CSV
    :param n:  int, number of unique feature vectors
    :param m:  int, size of feature set union
    :param binpath:  str, path to rapidNJ binary executable

    :return:  Biopython.Phylo object
    """
    sample = [int(m*random.random()) for _ in range(m)]
    weights = dict([(y, sample.count(y)) for y in set(sample)])

    # write directly to file to save memory
    outfile = tempfile.NamedTemporaryFile('w', prefix="cvz_boot_")
    outfile.write('{0:>5}\n'.format(n))

    if type(sym_diffs) is dict:
        for i in range(n):
            outfile.write('{}'.format(i))
            for j in range(n):
                d = 0
                if i != j:
                    key = (i, j) if i < j else (j, i)
                    d = sum(weights.get(y, 0) for y in sym_diffs[key])

                outfile.write(' {0:>2}'.format(d))
            outfile.write('\n')
    else:
        # retrieve symmetric differences from file
        infile = open(sym_diffs)
        for ij, line in enumerate(infile):
            i = ij // n  # row number
            j = ij % n  # column number
            if j == 0:
                outfile.write('{}'.format(i))

            sym_diff = line.strip().split(',')
            if len(sym_diff) == 1 and sym_diff[0] == '':
                d = 0  # split on blank line returns ['']
            else:
                d = sum(weights.get(int(y), 0) for y in sym_diff)

            outfile.write(' {0:>2}'.format(d))

            if j == n-1:
                outfile.write('\n')

        infile.close()
    outfile.flush()

    if callback:
        callback('generated dist matrix')

    # call rapidNJ
    cmd = [binpath, outfile.name, '-i', 'pd', '--no-negative-length']

    stdout = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
    handle = StringIO(stdout.decode('utf-8'))
    phy = Phylo.read(handle, 'newick')

    if callback:
        callback('rapidNJ complete')
    outfile.close()  # delete

    return phy


def build_trees(features, nboot=100, threads=1, use_file=True, callback=None):
    """
    Main function for clustering script.  Calculate matrices of symmetric differences
    between feature vectors as a distance measure.  Use bootstrap sampling to re-weight
    distances and then use the resulting matrices to reconstruct neighbor-joining trees.

    :param features:  list, tuples corresponding to features (genetic differences)
                      of genomes relative to a reference - generated by minimap2.py
    :param nboot:  int, number of bootstrap samples
    :param threads:  int, number of threads for parallel processing
    :param use_file:  bool, if True, write symmetric differences matrix to a temporary
                      file to minimize RAM consumption.
    :param callback:  optional progress monitoring callback function
    """

    if callback is None:
        def callback(msg):
            pass  # ignore calls to callback()

    # compute symmetric differences
    sym_diffs, m, labels = get_sym_diffs(features, use_file=use_file)
    callback('computed symmetric differences ({} variants)'.format(len(labels)))

    if len(labels) == 1:
        # data collapsed to a single variant
        return None, labels

    # bootstrap sampling and tree reconstruction
    if threads == 1:
        callback('launching single-threaded mode')
        trees = []
        for _ in range(nboot):
            # automatically detects if sym_diffs is a dict in memory or file
            phy = bootstrap(sym_diffs, n=len(labels), m=m)
            trees.append(phy)
    else:
        callback('launching pool with {} threads'.format(threads))
        with Pool(threads) as pool:
            results = [pool.apply_async(bootstrap, [sym_diffs, len(labels), m])
                       for _ in range(nboot)]
            trees = [r.get() for r in results]

    # clean up temporary files
    if use_file:
        os.remove(sym_diffs)

    callback('built NJ trees from {} bootstraps'.format(nboot))
    return trees, labels


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
        trees, labels = build_trees(
            filtered, nboot=args.nboot, threads=args.threads, callback=cb.callback
        )

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
