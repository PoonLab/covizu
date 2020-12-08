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
from multiprocessing import Pool, Manager

import sys
import os

from covizu.utils.seq_utils import load_vcf, filter_problematic


sys.setrecursionlimit(20000)  # fix for issue #127, default limit 1000


def recode_features(records, callback=None):
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
        key = tuple([tuple(x) for x in record['diffs']])
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

    return union, labels, indexed


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

    if callback:
        callback("Writing distance matrix to temporary file")

    outfile = tempfile.NamedTemporaryFile('w', prefix="cvz_boot_")
    outfile.write('{0:>5}\n'.format(n))
    for i in range(n):
        if callback and i % 100 == 0:
            callback("  row {} of {}".format(i, n))
        outfile.write('{}'.format(i))
        for j in range(n):
            if i == j:
                outfile.write(' {0:>2}'.format(0))
            else:
                sd = indexed[i] ^ indexed[j]  # symmetric difference
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


def build_trees(records, callback=None):
    """ Serial mode """
    union, labels, indexed = recode_features(records, callback=callback)
    trees = [bootstrap(union, indexed, args.binpath, callback=callback)
             for _ in range(args.nboot)]
    return trees, labels


def parse_args():
    """ Command-line interface """
    parser = argparse.ArgumentParser(
        description="Partition genomes by Pangolin lineage, generate distance"
                    "matrices and use neighbor-joining."
    )
    parser.add_argument("json", type=str,
                        help="input, JSON file")
    parser.add_argument("lineage", type=str, help="input, lineage to process")

    parser.add_argument("-o", "--outdir", type=str, default=os.getcwd(),
                        help="output, directory to export trees (one file per lineage). "
                             "Defaults to current working directory.")
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")
    parser.add_argument("-t", "--threads", type=int, default=2,
                        help="Number of threads to run, default 2.")

    return parser.parse_args()


if __name__ == "__main__":
    try:
        from mpi4py import MPI
    except ModuleNotFoundError:
        print("Script requires mpi4py - https://pypi.org/project/mpi4py/")
        sys.exit()

    comm = MPI.COMM_WORLD
    my_rank = comm.Get_rank()
    nprocs = comm.Get_size()

    # command-line execution
    args = parse_args()
    cb = Callback()

    cb.callback('loading JSON')
    with open(args.json) as handle:
        by_lineage = json.load(handle)

    records = by_lineage.get(args.lineage, None)
    if records is None:
        cb.callback("ERROR: JSON did not contain lineage {}".format(args.lineage))
        sys.exit()

    # generate distance matrices from bootstrap samples
    union, labels, indexed = recode_features(records, callback=cb.callback)
    trees = []
    for bn in range(args.nboot):
        if bn % nprocs != my_rank:
            continue
        phy = bootstrap(union, indexed, args.binpath, callback=cb.callback)
        trees.append(phy)
    comm.Barrier()
    result = comm.gather(trees, root=0)

    if my_rank == 0:
        trees = [phy for batch in result for phy in batch]  # flatten nested lists

        # export map of sequence labels to tip indices
        outfile = os.path.join(args.outdir, "{}.labels.csv".format(args.lineage))
        with open(outfile, 'w') as handle:
            handle.write("name,index\n")
            for i, names in enumerate(labels):
                for nm in names:
                    handle.write('{},{}\n'.format(nm, i))

        # export trees
        outfile = os.path.join(
            args.outdir, '{}.{}.nwk'.format(
                args.lineage, 'cons' if args.cons else 'boot'
            )
        )
        if trees is None:
            # lineage only has one variant, no meaningful tree
            with open(outfile, 'w') as handle:
                handle.write('({}:0);\n'.format(labels[0][0]))
        else:
            Phylo.write(trees, file=outfile, format='newick')
