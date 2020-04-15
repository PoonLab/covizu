from Bio import Phylo
cutoff = 0.0006

phy = Phylo.read('data/clusters.ft2.nwk', format='newick')

targets = []
for clade in phy.find_clades(terminal=True):
    # print('{},{}'.format(clade.name, clade.branch_length))
    # continue
    if clade.branch_length > cutoff:
        targets.append(clade)

for clade in targets:
    phy.prune(clade)

# remove confidence - Biopython doesn't write them out correctly
for node in phy.get_nonterminals():
    node.confidence = None

Phylo.write(phy, file='data/clusters.pruned.nwk', format='newick')
