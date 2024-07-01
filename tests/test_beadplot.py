import unittest
from io import StringIO
from covizu import beadplot
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

    def test_parse(self):
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

    def test_bad_parse(self):
        handle = StringIO()
        handle.write(
            "label,index\n"
            "hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25,0\n"
        )
        handle.seek(0)

        with self.assertRaises(KeyError):
            beadplot.parse_labels(handle)

    def test_reverse(self):
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
    def test_get_parents(self):
        tree = phylo_from_str("((A,B)D,C)E;\n")
        result = beadplot.get_parents(tree)

        # it's difficult to compare Clade objects, so extract names
        result2 = {}
        for child, parent in result.items():
            result2.update({child.name: parent.name})

        expected = {'A': 'D', 'B': 'D', 'D': 'E', 'C': 'E'}
        self.assertEqual(expected, result2)


class TestCollapsePolytomies(unittest.TestCase):
    def test_simple(self):
        tree = phylo_from_str("((A:1,B:1):0,C:1):0;")
        result = beadplot.collapse_polytomies(tree)
        expected = phylo_from_str("(A:1,B:1,C:1):0;")
        self.assertTrue(compare_phylo(expected, result))

        tree2 = phylo_from_str("((A:1,B:1):0,(C:1,D:1):1):0;")
        result = beadplot.collapse_polytomies(tree2)
        expected = phylo_from_str("(A:1,B:1,(C:1,D:1):1):0;")
        self.assertTrue(compare_phylo(expected, result))

    def test_labeling(self):
        tree3 = phylo_from_str("((A:1,B:1):1,C:0):0;")
        result = beadplot.collapse_polytomies(tree3)
        expected = phylo_from_str("(((A:1,B:1):1)C:0);")
        self.assertTrue(compare_phylo(expected, result))

        tree4 = phylo_from_str("((A:1,B:1):1,(C:0,D:0):1):0;")
        result = beadplot.collapse_polytomies(tree4)
        expected = phylo_from_str("((A:1,B:1):1,C|D:1):0;")
        self.assertTrue(compare_phylo(expected, result))


class TestIssues(unittest.TestCase):
    def test_issue150(self):
        # consensus of B.13 bootstrap trees
        tree = phylo_from_str(
            "(0:0.00000,(((1:0.00000,8:1.86735)0.71:1.54286,3:0.00000)"
            "0.72:1.43662,7:0.93878)0.88:2.02326,2:1.14286,4:0.95918,5:"
            "1.01020,6:2.19388)1.00:0.00000;"
        )
        labels = {
            '0': ['hCoV-19/USA/WI-UW-306/2020|EPI_ISL_436600|2020-04-14',
                  'hCoV-19/USA/WI-UW-278/2020|EPI_ISL_436572|2020-04-03',
                  'hCoV-19/USA/WI-UW-279/2020|EPI_ISL_436573|2020-04-03',
                  'hCoV-19/USA/WI-UW-287/2020|EPI_ISL_436581|2020-04-06',
                  'hCoV-19/USA/WI-UW-289/2020|EPI_ISL_436583|2020-04-08',
                  'hCoV-19/USA/WI-UW-291/2020|EPI_ISL_436585|2020-04-08',
                  'hCoV-19/USA/WI-UW-299/2020|EPI_ISL_436593|2020-04-13'],
            '1': ['hCoV-19/USA/WI-UW-310/2020|EPI_ISL_436604|2020-04-15'],
            '2': ['hCoV-19/USA/WI-UW-314/2020|EPI_ISL_436608|2020-04-16'],
            '3': ['hCoV-19/USA/WI-UW-315/2020|EPI_ISL_436609|2020-04-17',
                  'hCoV-19/USA/WI-UW-282/2020|EPI_ISL_436576|2020-04-06',
                  'hCoV-19/USA/WI-UW-286/2020|EPI_ISL_436580|2020-04-06'],
            '4': ['hCoV-19/USA/WI-UW-281/2020|EPI_ISL_436575|2020-04-06',
                  'hCoV-19/USA/WI-UW-292/2020|EPI_ISL_436586|2020-04-08',
                  'hCoV-19/USA/WI-UW-305/2020|EPI_ISL_436599|2020-04-14'],
            '5': ['hCoV-19/USA/WI-UW-294/2020|EPI_ISL_436588|2020-04-09'],
            '6': ['hCoV-19/USA/WI-UW-298/2020|EPI_ISL_436592|2020-04-13'],
            '7': ['hCoV-19/USA/WI-UW-301/2020|EPI_ISL_436595|2020-04-13'],
            '8': ['hCoV-19/USA/WI-UW-302/2020|EPI_ISL_436596|2020-04-13']
        }
        # this raised an exception
        tree = beadplot.annotate_tree(tree, labels)

if __name__ == '__main__':
    unittest.main()

