"""what clustering does"""
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

    :param records:  dict, samples keyed by unique mutation set
    :param callback:  optional, function for progress monitoring
    :param limit:  int, maximum number of variants to prevent memory allocation crashes
    :return:  dict, key-value pairs of all features indexed by integers
              dict, lists of labels by variant (identical feature vectors), keyed by index
              list, sets of feature vectors encoded by integers, by variant
    """
    # compress genomes with identical feature vectors
    fvecs = {}
    for muts, variant in records.items():
        key = tuple(tuple(x.split('|')) for x in muts.split(','))
        if key not in fvecs:
            fvecs.update({key: []})
        for sample in variant:
            label = (f"{sample['covv_virus_name']}|{sample['covv_location']}|"
                    f"{sample['covv_accession_id']}|{sample['covv_collection_date']}")

            fvecs[key].append(label)

    # limit to N most recently-sampled feature vectors
    intermed = [(max(label.split('|')[-1] for label in labels), key)
                for key, labels in fvecs.items()]
    intermed.sort(reverse=True)

    # generate union of all features
    if callback:
        callback(
            f"Reduced to {len(fvecs)} variants; generating feature set union")
    intermed_out = [{}, {}, []] # [{union}, {labels}, [idx]]
    for count, item in enumerate(intermed):
        fvec = item[1]
        intermed_out[1].update({str(count): fvecs[fvec]})
        if count < limit:
            for feat in fvec:
                if feat not in intermed_out[0]:
                    intermed_out[0].update({feat: len(intermed_out[0])})
            # recode feature vectors as integers (0-index)
            intermed_out[2].append(set(intermed_out[0][feat] for feat in fvec))

    return intermed_out


def bootstrap(input_union, idxed, binpath='rapidnj', callback=None, callfreq=1000):
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
    sample = [int(len(input_union) * random.random()) for _ in range(len(input_union))]
    weights = {y: sample.count(y) for y in sample}
    num_vec = len(idxed)

    # TODO: this is the slowest step - port to C? cache results to traverse
    # half matrix?
    with tempfile.NamedTemporaryFile('w', prefix="cvz_boot_") as temp_out:
        temp_out.write(f'{num_vec:>5}\n')
        for i in range(num_vec):
            if callback and i % callfreq == 0:
                callback(f"  row {i} of {num_vec}", level='DEBUG')
            temp_out.write(f'{i}')
            for j in range(num_vec):
                if i == j:
                    temp_out.write(f' {0:>2}')
                else:
                    sym_diff = idxed[i] ^ idxed[j]  # symmetric difference
                    difference = sum(weights.get(y, 0) for y in sym_diff)
                    temp_out.write(f' {difference:>2}')
            temp_out.write('\n')
        temp_out.flush()

        # call RapidNJ on temp file

        stdout = subprocess.check_output([
            binpath, temp_out.name, '-i', 'pd', '--no-negative-length'],
            stderr=subprocess.DEVNULL)
        out_handle = StringIO(stdout.decode('utf-8'))


    return Phylo.read(out_handle, 'newick')


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


def consensus(input_trees, cutoff=0.5, callback=None):
    """
    Generate a consensus tree by counting splits and using the splits with
    frequencies above the cutoff to resolve a star tree.
    :param trees:  iterable containing Phylo.BaseTree objects
    :param cutoff:  float, bootstrap threshold (default 0.5)
    :param callback:  function, optional callback
    :return:  Phylo.BaseTree
    """
    ntrees = 1
    tree = next(input_trees)

    # store terminal labels and branch lengths
    tip_index = {}
    for i, tip in enumerate(tree.get_terminals()):
        tip_index.update({tip.name: i})


    if callback:
        callback("Recording splits and branch lengths", level='DEBUG')
    splits = {}
    terminals = {tn: 0 for tn in tip_index}

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
            tree = next(input_trees)
            if callback:
                callback(f".. {ntrees} completed ", level="DEBUG")
            ntrees += 1
        except StopIteration:
            if callback:
                callback("... done", level='DEBUG')
            break

    # filter splits by frequency (support) threshold
    intermed = sorted([(k.count(',') + 1, k, v)
                      for k, v in splits.items() if v['count'] / ntrees >= cutoff])
    del splits  # free some memory

    # construct consensus tree
    if callback:
        callback("Building consensus tree", level='DEBUG')
    orphans = {
        tip_index[tname]: Clade(name=tname, branch_length=totlen / ntrees)
        for tname, totlen in terminals.items()
    }

    parse_intermed(intermed, ntrees, tip_index, orphans)

    return orphans.popitem()[1]

def parse_intermed(input_intermed, num_trees, tip_ind, in_orphans):
    """make pep8 compliant for consensus"""
    for _, key, val in input_intermed:
        # average branch lengths across relevant trees
        len_branch = val['sum'] / val['count']
        support = val['count'] / num_trees
        node = Clade(branch_length=len_branch, confidence=support)

        for child in map(int, key.split(',')):
            branch = in_orphans.pop(child, None)
            if branch:
                node.clades.append(branch)

        # use a single tip name to label ancestral node
        newkey = tip_ind[node.get_terminals()[0].name]
        in_orphans.update({newkey: node})


def build_trees(records, input_args, callback=None):
    """
    Serial mode, called from batch.py

    :param records:  list, feature vectors as dicts
    :param args:  Namespace, from argparse.ArgumentParser
    :param callback:  function, optional for progress monitoring
    """
    recode_union, recode_labels, idxed = recode_features(records, callback=callback)
    if len(idxed) == 1:
        # only one variant, no meaningful tree
        return None, recode_labels
    out_trees = [bootstrap(recode_union, idxed, input_args.binpath, callback=callback)
            for _ in range(input_args.nboot)]
    return out_trees, recode_labels


def parse_args():
    """ Command-line interface """
    parser = argparse.ArgumentParser(
        description="Partition genomes by Pangolin lineage, generate distance"
                    "matrices and use neighbor-joining."
    )
    parser.add_argument("json", type=str,
                        help="input, JSON file")
    parser.add_argument(
        "lineage",
        type=str,
        help="input, name of lineage to process ('deep' mode) or path to "
        "text file of minor lineage names ('flat' mode)")

    parser.add_argument(
        "--mode",
        type=str,
        default='deep',
        help="'flat' mode distributes many lineages across MPI nodes, whereas"
        "'deep' mode distributes bootstrap replicates for a single lineage "
        "across MPI nodes.  Defaults to 'deep'.")

    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default=os.getcwd(),
        help="output, directory to export trees (one file per lineage). "
        "Defaults to current working directory.")
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")
    parser.add_argument("--binpath", type=str, default='rapidnj',
                        help="Path to RapidNJ binary executable.")
    parser.add_argument("--timestamp", type=float, default=None,
                        help="option, timestamp to set callback function")
    parser.add_argument(
        "--max-variants",
        type=int,
        default=10000,
        help="option, limit number of variants per lineage, prioritizing the "
        "most recently sampled variants")

    return parser.parse_args()


def unpack_recoded(input_recoded, input_lineage, callback=None):
    """
    Recover dictionary from JSON.
    :param recoded:  dict, directly returned from json.load
    :param lineage:  str, PANGO lineage specifier
    :param callback:  optional callback function
    """
    rdata = input_recoded.get(input_lineage, None)
    if rdata is None:
        if callback:
            callback(f"ERROR: JSON did not contain lineage {args.lineage}")
        sys.exit()

    union1 = rdata['union']  # unpack JSON data
    union2 = {}
    for feat, idx in union1.items():
        typ, pos, length = feat.split('|')
        if typ == '-':
            # length of deletion is an integer
            key = tuple([typ, int(pos), int(length)])
        else:
            key = tuple([typ, int(pos), length])
        union2.update({key: idx})

    recoded_indexed = [set(l) for l in rdata['indexed']]  # see #335
    return union2, rdata['labels'], recoded_indexed

#   Called by batch.py via subprocess to handle lineages with excessive
#   numbers of genomes, to process via MPI
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
    cb = Callback(t0=args.timestamp, my_rank=my_rank, nprocs=nprocs)

    # import lineage data from file
    with open(args.json, encoding='utf-8') as handle:
        recoded = json.load(handle)

    if args.mode == 'deep':
        union, labels, indexed = unpack_recoded(
            recoded, args.lineage, callback=cb.callback)

        # export map of sequence labels to tip indices
        lineage_name = args.lineage.replace('/', '_')  # issue #297

        outfile = os.path.join(args.outdir, f'{lineage_name}.nwk')
        if len(indexed) == 1:
            # lineage only has one variant, no meaningful tree
            if my_rank == 0:
                with open(outfile, 'w', encoding='utf-8') as handle:
                    handle.write(f"({labels['0'][0]}:0);\n")
        else:
            # MPI processing
            trees = []
            for bn in range(args.nboot):
                if bn % nprocs == my_rank:
                    phy = bootstrap(
                        union, indexed, args.binpath, callback=cb.callback)
                    trees.append(phy)
            comm.Barrier()  # wait for other processes to finish
            result = comm.gather(trees, root=0)

            # only head node exports trees
            if my_rank == 0:
                # flatten nested lists
                trees = [phy for batch in result for phy in batch]
                Phylo.write(trees, file=outfile, format='newick')

    elif args.mode == 'flat':
        # load list of lineages from text file
        minor_lineages = []
        with open(args.lineage, encoding='utf-8') as handle:
            for line in handle:
                minor_lineages.append(line.strip())

        for li, lineage in enumerate(minor_lineages):
            if li % nprocs != my_rank:
                continue
            cb.callback(f"starting {lineage}")
            union, labels, indexed = unpack_recoded(
                recoded, lineage, callback=cb.callback)

            lineage_name = lineage.replace('/', '_')  # issue #297
            outfile = os.path.join(args.outdir, f'{lineage_name}.nwk')
            if len(indexed) == 1:
                # lineage only has one variant, no meaningful tree
                with open(outfile, 'w', encoding='utf-8') as handle:
                    handle.write(f"({labels['0'][0]}:0);\n")
            else:
                trees = [
                    bootstrap(
                        union,
                        indexed,
                        args.binpath,
                        callback=cb.callback) for _ in range(
                        args.nboot)]
                Phylo.write(trees, file=outfile, format='newick')
    else:
        cb.callback(f"Unexpected mode argument {args.mode} in clustering.py")
        sys.exit()
