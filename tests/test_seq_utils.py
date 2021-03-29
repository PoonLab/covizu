import unittest
from covizu.utils.seq_utils import *
from io import StringIO


class TestIterFasta(unittest.TestCase):
    def setUp(self):
        self.expected = \
            [
                ('hCoV-19/Canada/Qc-L00240569/2020|EPI_ISL_465679|2020-03-27',
                    'GGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGA'
                 ),
                ('hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25',
                    "CATCTACAGATACTTGTTTTGCTAACAAACATGCTGATTTTGACACATGGTTTAGCCAGCGTGGTGGTAGTTATACTAATTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAA"
                 ),
                ('hCoV-19/HongKong/HKU-200723-093/2020|EPI_ISL_497860|2020-01-25',
                    "GTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTTGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTT"
                 )
            ]

    def testIterFasta(self):
        handle = StringIO()
        handle.write(
            ">hCoV-19/Canada/Qc-L00240569/2020|EPI_ISL_465679|2020-03-27\n"
            "GGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTG\n"
            "TGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGA\n"
            ">hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25\n"
            "CATCTACAGATACTTGTTTTGCTAACAAACATGCTGATTTTGACACATGGTTTAGCCAGCGTGGTGGTAGTTATACTAAT\n"
            "TACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAA\n"
            ">hCoV-19/HongKong/HKU-200723-093/2020|EPI_ISL_497860|2020-01-25\n"
            "GTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTTGTCCG\n"
            "GGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTT\n"
        )
        handle.seek(0)
        result = list(iter_fasta(handle))
        self.assertEqual(self.expected, result)


class TestConvertFasta(unittest.TestCase):
    def setUp(self):
        self.expected = \
            [
                ['NC_045512.2, Complete genome',
                 'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC'
                ],
                ['NC_045512.3, Test genome',
                 "CTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAA"
                ]
            ]

    def testConvertFasta(self):
        handle = StringIO()
        handle.write(
            ">NC_045512.2, Complete genome\n"
            "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA\n"
            "CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC\n"
            ">NC_045512.3, Test genome\n"
            "CTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAA\n"
            "ATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAA\n"
        )
        handle.seek(0)
        result = convert_fasta(handle)
        # Test case fails if line starts with a #
        self.assertEqual(self.expected, result)


class TestTotalMissing(unittest.TestCase):

    def testMissing(self):
        row = {
            'missing': [(0,6), (5203,5222), (29844,29903)]
        }
        res = total_missing(row)
        self.assertEqual(84, res)

    def testListRowMissing(self):
        row = [["test"], ["test"], [(0,6), (5203,5222), (29844,29903)]]
        res = total_missing(row)
        self.assertEqual(84, res)


class TestQPois(unittest.TestCase):
    def setUp(self):
        self.qp = QPois(quantile=1-0.005, rate=0.0655, maxtime=1e3, origin='2019-12-01')

    def testIsOutlier(self):
        coldate = '2021-01-18'
        ndiffs = 30
        result = self.qp.is_outlier(coldate, ndiffs)
        self.assertEqual(result, False)

    def testIsOutlierManyDiff(self):
        coldate = '2021-01-18'
        ndiffs = 42
        result = self.qp.is_outlier(coldate, ndiffs)
        self.assertEqual(result, True)

    def testIsOutlierFewDiff(self):
        coldate = '2021-01-18'
        ndiffs = 14
        result = self.qp.is_outlier(coldate, ndiffs)
        self.assertEqual(result, True)


class TestApplyFeatures(unittest.TestCase):
    def setUp(self):
        self.expected = 'NNNNNTGTTTNTNNCCTTCCCAGGTTNTNTAGCAAACAAC'

    def testApplyFeatures(self):
        diffs = [('~', 5, 'T'), ('~', 7, 'T'), ('~', 9, 'T'), ('~', 11, 'T'), ('~', 25, 'T'), ('~', 27, 'T'), ('~', 29, 'T'), ('~', 31, 'G'), ('~', 35, 'A')]
        missing = [(0, 6), (10, 14), (25, 29)]
        refseq = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAAC"
        result = apply_features(diffs, missing, refseq)
        self.assertEqual(self.expected, result)


if __name__ == '__main__':
    unittest.main()
