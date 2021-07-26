import random
import json
import argparse
import tempfile
import subprocess
from io import StringIO
import sys
import os

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from covizu.utils.progress_utils import Callback


sys.setrecursionlimit(20000)  # fix for issue #127, default limit 1000


def recode_features(records, callback=None, limit=10000):
    """
    Recode feature vectors with integer indices based on set union.
    Pass results to bootstrap() to reconstruct trees by neighbor-joining method.

    :param records:  list, dict for each record
    :param callback:  optional, function for progress monitoring
    :param limit:  int, maximum number of variants to prevent memory allocation crashes
    :return:  dict, key-value pairs of all features indexed by integers
              list, nested list of labels by variant (identical feature vectors)
              list, sets of feature vectors encoded by integers, by variant
    """
    # compress genomes with identical feature vectors
    fvecs = {}
    for record in records:
        label = '{label}|{accession}|{region}|{country}|{division}|{coldate}'.format(**record)
        key = tuple([tuple(x) for x in record['diffs']])
        if key not in fvecs:
            fvecs.update({key: []})
        fvecs[key].append(label)

    # limit to N most recently-sampled feature vectors
    intermed = [(max([l.split('|')[-1] for l in label]), key) for key, label in fvecs.items()]
    intermed.sort(reverse=True)

    # generate union of all features
    if callback:
        callback("Reduced to {} variants; generating feature set union".format(len(fvecs)))
    union = {}
    labels = []
    indexed = []
    for count, item in enumerate(intermed):
        _, key = item  # discard max coldate
        labels.append(fvecs[key])
        if count < limit:
            for feat in key:
                if feat not in union:
                    union.update({feat: len(union)})
            # recode feature vectors as integers (0-index)
            indexed.append(set(union[feat] for feat in key))

    return union, labels, indexed


def bootstrap(union, indexed, binpath='rapidnj', callback=None, callfreq=1000):
    """
    Sample features from set union at random with replacement.  We use the
    result to weight the symmetric differences when calculating pairwise
    distances.  Pass the resulting distance matrix to RapidNJ to reconstruct
    a tree.

    :param union:  set, all observed genetic differences from reference (features)
    :param indexed:  list, feature vectors encoded as integers
    :param binpath:  str, path to RapidNJ binary executable
    :param callback:  function, optional for progress monitoring
    :param callfreq:  int, sampling interval for callback

    :return:  Bio.Phylo.BaseTree object
    """
    sample = [int(len(union) * random.random()) for _ in range(len(union))]
    weights = dict([(y, sample.count(y)) for y in sample])
    n = len(indexed)

    if callback:
        callback("Writing distance matrix to temporary file")

    # TODO: this is the slowest step - port to C? cache results to traverse half matrix?
    outfile = tempfile.NamedTemporaryFile('w', prefix="cvz_boot_")
    outfile.write('{0:>5}\n'.format(n))
    for i in range(n):
        if callback and i % callfreq == 0:
            callback("  row {} of {}".format(i, n), level='DEBUG')
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
    :param callback:  function, optional callback
    :return:  Phylo.BaseTree
    """
    ntrees = 1
    tree = next(trees)

    # store terminal labels and branch lengths
    tip_index = {}
    for i, tip in enumerate(tree.get_terminals()):
        tip_index.update({tip.name: i})
    ntips = len(tip_index)

    if callback:
        callback("Recording splits and branch lengths")
    splits = {}
    terminals = dict([(tn, 0) for tn in tip_index.keys()])

    while True:
        # record terminal branch lengths
        for tip in tree.get_terminals():
            terminals[tip.name] += tip.branch_length

        # record splits in tree
        tree = label_nodes(tree, tip_index)  # aggregates tip indices down tree
        for node in tree.get_nonterminals():
            key = ','.join(map(str, node.tip_index))
            if key not in splits:
                splits.update({key: {'sum': 0., 'count': 0}})

            if node.branch_length is not None:
                # None interpreted as zero length (e.g., root branch)
                splits[key]['sum'] += node.branch_length
            splits[key]['count'] += 1
        try:
            tree = next(trees)
            if callback:
                callback(".. {} completed ".format(ntrees), level="DEBUG")
            ntrees += 1
        except StopIteration:
            if callback:
                callback("... done", level='DEBUG')
            break

    # filter splits by frequency (support) threshold
    intermed = [(k.count(',')+1, k, v) for k, v in splits.items() if v['count']/ntrees >= cutoff]
    intermed.sort()  # sort by level (tips to root)
    del splits  # free some memory

    # construct consensus tree
    if callback:
        callback("Building consensus tree")
    orphans = dict([
        (tip_index[tname], Clade(name=tname, branch_length=totlen/ntrees))
        for tname, totlen in terminals.items()
    ])

    for _, key, val in intermed:
        # average branch lengths across relevant trees
        bl = val['sum'] / val['count']
        support = val['count'] / ntrees
        node = Clade(branch_length=bl, confidence=support)

        for child in map(int, key.split(',')):
            branch = orphans.pop(child, None)
            if branch:
                node.clades.append(branch)

        # use a single tip name to label ancestral node
        newkey = tip_index[node.get_terminals()[0].name]
        orphans.update({newkey: node})

    return orphans.popitem()[1]


def build_trees(records, args, callback=None):
    """
    Serial mode, called from batch.py

    :param records:  list, feature vectors as dicts
    :param args:  Namespace, from argparse.ArgumentParser
    :param callback:  function, optional for progress monitoring
    """
    union, labels, indexed = recode_features(records, callback=callback)
    if len(indexed) == 1:
        # only one variant, no meaningful tree
        trees = None
    else:
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
    parser.add_argument("lineage", type=str,
                        help="input, name of lineage to process ('deep' mode) or path to "
                             "text file of minor lineage names ('flat' mode)")

    parser.add_argument("--mode", type=str, default='deep',
                        help="'flat' mode distributes many lineages across MPI nodes, whereas"
                             "'deep' mode distributes bootstrap replicates for a single lineage "
                             "across MPI nodes.  Defaults to 'deep'.")

    parser.add_argument("-o", "--outdir", type=str, default=os.getcwd(),
                        help="output, directory to export trees (one file per lineage). "
                             "Defaults to current working directory.")
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")
    parser.add_argument("--binpath", type=str, default='rapidnj',
                        help="Path to RapidNJ binary executable.")
    parser.add_argument("--timestamp", type=float, default=None,
                        help="option, timestamp to set callback function")
    parser.add_argument("--max-variants", type=int, default=10000,
                        help="option, limit number of variants per lineage, prioritizing the "
                             "most recently sampled variants")

    return parser.parse_args()


if __name__ == "__main__":
    """
    Called by batch.py via subprocess to handle lineages with excessive
    numbers of genomes, to process via MPI
    """
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
    cb = Callback(t0=args.timestamp, my_rank=my_rank, nprocs=nprocs)

    # import lineage data from file
    # FIXME: this consumes a few minutes to load the entire data set
    cb.callback('loading JSON')
    with open(args.json) as handle:
        by_lineage = json.load(handle)

    if args.mode == 'deep':
        records = by_lineage.get(args.lineage, None)
        if records is None:
            cb.callback("ERROR: JSON did not contain lineage {}".format(args.lineage))
            sys.exit()

        # generate distance matrices from bootstrap samples [[ MPI ]]
        union, labels, indexed = recode_features(records, callback=cb.callback, limit=args.max_variants)

        # export map of sequence labels to tip indices
        lineage_name = args.lineage.replace('/', '_')  # issue #297
        if my_rank == 0:
            csvfile = os.path.join(args.outdir, "{}.labels.csv".format(lineage_name))
            with open(csvfile, 'w') as handle:
                handle.write("name,index\n")
                for i, names in enumerate(labels):
                    for nm in names:
                        handle.write('{},{}\n'.format(nm, i))

        outfile = os.path.join(args.outdir, '{}.nwk'.format(lineage_name))
        if len(indexed) == 1:
            # lineage only has one variant, no meaningful tree
            if my_rank == 0:
                with open(outfile, 'w') as handle:
                    handle.write('({}:0);\n'.format(labels[0][0]))
        else:
            # MPI processing
            trees = []
            for bn in range(args.nboot):
                if bn % nprocs == my_rank:
                    phy = bootstrap(union, indexed, args.binpath)
                    trees.append(phy)
            comm.Barrier()  # wait for other processes to finish
            result = comm.gather(trees, root=0)

            # only head node exports trees
            if my_rank == 0:
                trees = [phy for batch in result for phy in batch]  # flatten nested lists
                Phylo.write(trees, file=outfile, format='newick')

    elif args.mode == 'flat':
        # load list of lineages from text file
        minor_lineages = []
        with open(args.lineage) as handle:
            for line in handle:
                minor_lineages.append(line.strip())

        for li, lineage in enumerate(minor_lineages):
            if li % nprocs != my_rank:
                continue

            records = by_lineage.get(lineage, None)
            if records is None:
                cb.callback("ERROR: JSON did not contain lineage {}".format(args.lineage))
                sys.exit()

            union, labels, indexed = recode_features(records, callback=cb.callback, limit=args.max_variants)

            # export map of sequence labels to tip indices
            lineage_name = lineage.replace('/', '_')  # issue #297
            csvfile = os.path.join(args.outdir, "{}.labels.csv".format(lineage_name))
            with open(csvfile, 'w') as handle:
                handle.write("name,index\n")
                for i, names in enumerate(labels):
                    for nm in names:
                        handle.write('{},{}\n'.format(nm, i))

            outfile = os.path.join(args.outdir, '{}.nwk'.format(lineage_name))
            if len(indexed) == 1:
                # lineage only has one variant, no meaningful tree
                with open(outfile, 'w') as handle:
                    handle.write('({}:0);\n'.format(labels[0][0]))
            else:
                trees = [bootstrap(union, indexed, args.binpath)
                         for _ in range(args.nboot)]
                Phylo.write(trees, file=outfile, format='newick')
    else:
        cb.callback("Unexpected mode argument {} in clustering.py".format(args.mode))
        sys.exit()
