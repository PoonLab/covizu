import random
from Bio import Phylo
import glob
import os
from covizu import clustering 
from covizu.utils.progress_utils import Callback
from covizu.utils.batch_utils import manage_collapsed_nodes
import json
import csv
import argparse
from Bioplus import prunetree

def parse_args():
    parser = argparse.ArgumentParser(description='Prune tree')
    parser.add_argument('recoded', type=str, help='Path to recoded JSON')
    parser.add_argument('--prune', type=str, help='Path to write pruned trees', default=None)
    parser.add_argument('--ctree', type=str, help='Path to consensus tree', default=None)
    parser.add_argument('--bootstrap', type=str, help='Path to bootstrap trees', default=None)
    parser.add_argument('--write-ctree', type=str, help='Path to write consensus tree', default=None)
    parser.add_argument('--write-labels',  type=str, help='Path to write labels', default=None)
    return parser.parse_args()


def write_labels(outpath, clabel_dict, lineage):
    # Write labels CSV
    with open(f'{outpath}/labels.{lineage}.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        for key, value in clabel_dict.items():
            writer.writerow([key, value])


if __name__ == '__main__':
    args = parse_args()
    cb = Callback()

    with open(args.recoded, 'r') as f:
        recoded = json.load(f)

    if (args.bootstrap):
        nwk_trees = glob.glob(args.bootstrap + '/*.nwk')
    elif (args.ctree):
        nwk_trees = glob.glob(args.ctree + '/*.nwk')
    else:
        cb.callback("No trees provided")
        exit(-1)

    for f in nwk_trees:
        lineage = os.path.basename(f).split('.nwk')[0]
        label_dict = recoded[lineage]['labels']
        if (len(label_dict) == 1):
            cb.callback(f"Skipping {lineage}")
            continue

        trees = None
        with open(f, encoding='utf-8') as nwk_file:
            if args.ctree:
                tree = Phylo.read(nwk_file, 'newick')
            else:
                trees = Phylo.parse(nwk_file, "newick")
                tree = clustering.consensus(trees, cutoff=0.5)

            clabel_dict = manage_collapsed_nodes(label_dict, tree)

            if args.prune:
                if os.path.exists(f"{args.prune}/{lineage}.n500.nwk"):
                    continue
             
                tree = prunetree.prune_tree(tree, 500)

            for tip in tree.get_terminals():
                if tip.name not in clabel_dict:
                    continue
                count = len(clabel_dict[tip.name])
                if count < 2:
                    continue
                tip.name += '_' 
                tip.split(n=count, branch_length=0.)
                tip.name = None  # remove internal label


            if args.write_ctree:
                # Write consensus tree
                Phylo.write(tree, f'{args.write_ctree}/{lineage}.nwk', 'newick')
            
            if args.write_labels:
                write_labels(args.write_labels, clabel_dict, lineage)
            
            if args.prune:
                Phylo.write(tree, f'{args.prune}/{lineage}.n500.nwk', 'newick')
