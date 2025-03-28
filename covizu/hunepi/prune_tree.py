import random
from Bio import Phylo
from glob import glob
import os
from covizu import clustering 
from covizu.utils.progress_utils import Callback
from covizu.utils.batch_utils import manage_collapsed_nodes
import json
import csv


cb = Callback()

with open('iss539/recode.json', 'r') as f:
    recoded = json.load(f)

# Consensus trees
for f in glob('iss539/ctree/*.nwk'):
    lineage = os.path.basename(f).split('.nwk')[0]
    label_dict = recoded[lineage]['labels']
    if (len(label_dict) == 1):
        cb.callback(f"Skipping {lineage}")
        continue

    nwkfile = open(f, encoding='utf-8')
    tree = Phylo.read(nwkfile, "newick")
    #tree = clustering.consensus(trees, cutoff=0.5)

    #clabel_dict = manage_collapsed_nodes(label_dict, tree)

    # Write labels CSV
    # with open(f'labels.{lineage}.csv', 'w') as csvfile:
    #     writer = csv.writer(csvfile)
    #     for key, value in clabel_dict.items():
    #         writer.writerow([key, value])

    all_tips = tree.get_terminals()

    if os.path.exists("iss539/pruned/{}.n500.nwk".format(lineage)):
        continue

    cb.callback(f"Lineage {lineage} has {len(all_tips)} tips")
    
    # expand tips with multiple samples
    for tip in all_tips:
        if tip.name not in label_dict:
            continue
        count = len(label_dict[tip.name])
        if count < 2:
            continue
        tip.name += '_' 
        tip.split(n=count, branch_length=0.)
        tip.name = None  # remove internal label
        cb.callback(count)
    
    all_tips = tree.get_terminals()
    cb.callback(f"... expanded to {len(all_tips)} tips.")
    
    if len(all_tips) > 500:
        # random sample without replacement
        selected_tips = set(random.sample(all_tips, 500))
        for tip in all_tips:
            if tip in selected_tips:
                continue
            tree.prune(tip)

    Phylo.write(tree, "iss539/pruned/{}.n500.nwk".format(lineage), "newick")
    break
