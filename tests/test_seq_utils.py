import unittest
from io import StringIO
from covizu.utils.seq_utils import iter_fasta, convert_fasta, total_missing, QPois, apply_features


class TestIterFasta(unittest.TestCase):
    def setUp(self):
        self.expected = \
            [
                ('hCoV-19/Canada/Qc-L00240569/2020|EPI_ISL_465679|2020-03-27',
                    ('GGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTG'
                    'TTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTG'
                    'CACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGA')
                 ),
                ('hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25',
                    ('CATCTACAGATACTTGTTTTGCTAACAAACATGCTGATTTTGACACATGGTTTAG'
                    'CCAGCGTGGTGGTAGTTATACTAATTACAGGTTCGCGACGTGCTCGTACGTGGCT'
                    'TTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAA')
                 ),
                ('hCoV-19/HongKong/HKU-200723-093/2020|EPI_ISL_497860|2020-01-25',
                    ('GTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATC'
                    'AGCACATCTAGGTTTTGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTG'
                    'GTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTT')
                 )
            ]

    def test_iter_fasta(self):
        handle = StringIO()
        handle.write(
            ">hCoV-19/Canada/Qc-L00240569/2020|EPI_ISL_465679|2020-03-27\n"
            "GGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTA"
            "AACGAACTTTAAAATCTG\n"
            "TGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTG"
            "TCGTTGACAGGACACGA\n"
            ">hCoV-19/HongKong/HKPU6_2101/2020|EPI_ISL_417178|2020-01-25\n"
            "CATCTACAGATACTTGTTTTGCTAACAAACATGCTGATTTTGACACATGGTTTAGCCAGCGTG"
            "GTGGTAGTTATACTAAT\n"
            "TACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGG"
            "CACGTCAACATCTTAAA\n"
            ">hCoV-19/HongKong/HKU-200723-093/2020|EPI_ISL_497860|2020-01-25\n"
            "GTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCA"
            "CATCTAGGTTTTGTCCG\n"
            "GGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCA"
            "ACTCAGTTTGCCTGTTT\n"
        )
        handle.seek(0)
        result = list(iter_fasta(handle))
        self.assertEqual(self.expected, result)


class TestConvertFasta(unittest.TestCase):
    def setUp(self):
        self.expected = \
            [
                ['NC_045512.2, Complete genome',
                 'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTC'
                 'TCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC'
                ],
                ['NC_045512.3, Test genome',
                 "CTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGC"
                 "AAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAA"
                ]
            ]

    def test_convert_fasta(self):
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

    def test_missing(self):
        row = {
            'missing': [(0,6), (5203,5222), (29844,29903)]
        }
        res = total_missing(row)
        self.assertEqual(84, res)

    def test_list_row_missing(self):
        row = [["test"], ["test"], [(0,6), (5203,5222), (29844,29903)]]
        res = total_missing(row)
        self.assertEqual(84, res)


class TestApplyFeatures(unittest.TestCase):
    def setUp(self):
        self.expected = 'NNNNNTGTTTNTNNCCTTCCCAGGTTNTNTAGCAAACAAC'

    def test_apply_features(self):
        diffs = [('~', 5, 'T'), ('~', 7, 'T'), ('~', 9, 'T'), ('~', 11, 'T'), ('~', 25, 'T'),
                 ('~', 27, 'T'), ('~', 29, 'T'), ('~', 31, 'G'), ('~', 35, 'A')]
        missing = [(0, 6), (10, 14), (25, 29)]
        refseq = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAAC"
        result = apply_features(diffs, missing, refseq)
        self.assertEqual(self.expected, result)


class TestQPois(unittest.TestCase):
    def setUp(self):
        self.poisson_expected = QPois(quantile=1-0.005, rate=0.0655, maxtime=1e3, origin='2019-12-01')

    def test_is_outlier(self):
        coldate = '2021-01-18'
        ndiffs = 30
        result = self.poisson_expected.is_outlier(coldate, ndiffs)
        self.assertEqual(result, False)

    def test_is_outlier_many_diffs(self):
        coldate = '2021-01-18'
        ndiffs = 42
        result = self.poisson_expected.is_outlier(coldate, ndiffs)
        self.assertEqual(result, True)

    def test_is_outlier_few_diffs(self):
        coldate = '2021-01-18'
        ndiffs = 14
        result = self.poisson_expected.is_outlier(coldate, ndiffs)
        self.assertEqual(result, True)

# Load VCF

# Filter Problematic

# SC2Locator

if __name__ == '__main__':
    unittest.main()
