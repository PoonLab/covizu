import unittest
from covizu import beadplot
from io import StringIO
from Bio import Phylo


def phylo_from_str(nwk):
    """ Returns a Phylo.BaseTree object given Newick string """
    handle = StringIO()
    handle.write(nwk)
    handle.seek(0)
    return Phylo.read(handle, "newick")


def tip_names(tree):
    names = [tip.name for tip in tree.get_terminals()]
    names.sort()
    return tuple(names)

def clades(tree):
    return dict([
        (tip_names(node), node.branch_length)
        for node in tree.get_nonterminals()
    ])

def compare_phylo(tree1, tree2):
    """ Compare Phylo.BaseTree objects """
    return clades(tree1) == clades(tree2)


class TestParseLabels(unittest.TestCase):
    def setUp(self):
        self.expected = {
            '0': [
                'hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25',
                'hCoV-19/HongKong/HKU-200723-093/2020|EPI_ISL_497860|2020-01-25'
            ],
            '1': [
                'hCoV-19/Finland/14M82/2020|EPI_ISL_418411|2020-03-14'
            ],
            '2': [
                'hCoV-19/Australia/VIC129/2020|EPI_ISL_419732|2020-03-20'
            ]
        }

    def testParse(self):
        handle = StringIO()
        handle.write(
            "name,index\n"
            "hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25,0\n"
            "hCoV-19/HongKong/HKU-200723-093/2020|EPI_ISL_497860|2020-01-25,0\n"
            "hCoV-19/Finland/14M82/2020|EPI_ISL_418411|2020-03-14,1\n"
            "hCoV-19/Australia/VIC129/2020|EPI_ISL_419732|2020-03-20,2\n"
        )
        handle.seek(0)
        result = beadplot.parse_labels(handle)
        self.assertEqual(self.expected, result)

    def testBadParse(self):
        handle = StringIO()
        handle.write(
            "label,index\n"
            "hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25,0\n"
        )
        handle.seek(0)

        with self.assertRaises(KeyError):
            beadplot.parse_labels(handle)

    def testReverse(self):
        handle = StringIO()
        handle.write(
            "index,name\n"
            "0,hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25\n"
            "0,hCoV-19/HongKong/HKU-200723-093/2020|EPI_ISL_497860|2020-01-25\n"
            "1,hCoV-19/Finland/14M82/2020|EPI_ISL_418411|2020-03-14\n"
            "2,hCoV-19/Australia/VIC129/2020|EPI_ISL_419732|2020-03-20\n"
        )
        handle.seek(0)
        result = beadplot.parse_labels(handle)
        self.assertEqual(self.expected, result)


class TestGetParents(unittest.TestCase):
    def testGetParents(self):
        tree = phylo_from_str("((A,B)D,C)E;\n")
        result = beadplot.get_parents(tree)

        # it's difficult to compare Clade objects, so extract names
        result2 = {}
        for child, parent in result.items():
            result2.update({child.name: parent.name})

        expected = {'A': 'D', 'B': 'D', 'D': 'E', 'C': 'E'}
        self.assertEqual(expected, result2)


class TestCollapsePolytomies(unittest.TestCase):
    def testSimple(self):
        tree = phylo_from_str("((A:1,B:1):0,C:1):0;")
        result = beadplot.collapse_polytomies(tree)
        expected = phylo_from_str("(A:1,B:1,C:1):0;")
        self.assertTrue(compare_phylo(expected, result))

        tree2 = phylo_from_str("((A:1,B:1):0,(C:1,D:1):1):0;")
        result = beadplot.collapse_polytomies(tree2)
        expected = phylo_from_str("(A:1,B:1,(C:1,D:1):1):0;")
        self.assertTrue(compare_phylo(expected, result))

    def testLabeling(self):
        tree3 = phylo_from_str("((A:1,B:1):1,C:0):0;")
        result = beadplot.collapse_polytomies(tree3)
        expected = phylo_from_str("(((A:1,B:1):1)C:0);")
        self.assertTrue(compare_phylo(expected, result))

        tree4 = phylo_from_str("((A:1,B:1):1,(C:0,D:0):1):0;")
        result = beadplot.collapse_polytomies(tree4)
        expected = phylo_from_str("((A:1,B:1):1,C|D:1):0;")
        self.assertTrue(compare_phylo(expected, result))