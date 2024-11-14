import random
from Bio import Phylo
import glob
import os
from covizu import clustering 
from covizu.utils.progress_utils import Callback
from covizu.utils.batch_utils import manage_collapsed_nodes
import json
import csv


cb = Callback()

with open('2024-08-13/recode.json', 'r') as f:
    recoded = json.load(f)

# Bootstrap trees
nwk_trees = glob.glob('2024-08-13/*.nwk')

# Consensus trees
# nwk_trees = glob.glob('ctree/*.nwk')

for f in nwk_trees:
    lineage = os.path.basename(f).split('.nwk')[0]
    label_dict = recoded[lineage]['labels']
    if (len(label_dict) == 1):
        cb.callback(f"Skipping {lineage}")
        continue

    trees = None
    with open(f, encoding='utf-8') as outfile:
        trees = Phylo.parse(outfile, "newick")

        tree = clustering.consensus(trees, cutoff=0.5)

        clabel_dict = manage_collapsed_nodes(label_dict, tree)

        # Write labels CSV
        # with open(f'labels.{lineage}.csv', 'w') as csvfile:
        #     writer = csv.writer(csvfile)
        #     for key, value in clabel_dict.items():
        #         writer.writerow([key, value])

        all_tips = tree.get_terminals()

        if os.path.exists("pruned/{}.n500.nwk".format(lineage)):
            continue

        cb.callback(f"Lineage {lineage} has {len(all_tips)} tips")
        if len(all_tips) > 500:
            selected_tips = set(random.sample(all_tips, 500))
            
            tips_to_remove = [tip for tip in all_tips if tip not in selected_tips]
            
            for tip in tips_to_remove:
                tree.prune(tip)

        for tip in tree.get_terminals():
            if tip.name not in clabel_dict:
                continue
            count = len(clabel_dict[tip.name])
            if count < 2:
                continue
            tip.name += '_' 
            tip.split(n=count, branch_length=0.)
            tip.name = None  # remove internal label

        Phylo.write(tree, "pruned/{}.n500.nwk".format(lineage), "newick")
