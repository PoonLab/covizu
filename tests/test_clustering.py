import os
import unittest
from covizu import clustering
from Bio import Phylo as phy
from Bio.Phylo.BaseTree import Clade

# set up
basepath = os.path.dirname(os.path.abspath(__file__))
tree1 = phy.read(os.path.join(basepath, 'test_tree1.nwk'), format='newick', rooted=True)
tip_index1 = {}
for i, tip in enumerate(tree1.get_terminals()):
    tip_index1.update({tip.name: i})

tree2 = phy.read(os.path.join(basepath, 'test_tree2.nwk'), format='newick', rooted=True)
tip_index2 = {}
for i, tip in enumerate(tree2.get_terminals()):
    tip_index2.update({tip.name: i})

class Test_Label_Nodes(unittest.TestCase):
    def test_label_nodes(self):
        output_tree1 = clustering.label_nodes(tree1, tip_index1)
        self.assertEqual(tree1, output_tree1)

    def test_label_nodes2(self):
        output_tree2 = clustering.label_nodes(tree2, tip_index2)
        self.assertEqual(tree2, output_tree2)

class Test_Consensus(unittest.TestCase):
    def setUp(self):
        self.expected = [Clade(branch_length=0.55, confidence=1.0),
                         Clade(branch_length=0.95, confidence=1.0),
                         Clade(branch_length=1.0, confidence=1.0)]

    def test_consensus(self):
        out_tree = clustering.consensus(iter([tree1, tree2])).clades
        for ind, clade in enumerate(self.expected):
            self.assertEqual(clade.branch_length, out_tree[ind].branch_length)
            self.assertEqual(clade.confidence, out_tree[ind].confidence)
