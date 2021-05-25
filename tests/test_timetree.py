import unittest
from covizu.treetime import *

# Parse Nexus File - Nexus to Newick

class TestRetrieveGenome(unittest.TestCase):
    def setUp(self):
        self.expected = \
        {
            '|B.1.1.171|2020-03-27': 'NNNNNNGGTTTATACCTTCCCAGGTAAC-AACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTC',
            '|B.1.265|2020-03-27': 'NNNNNNGGTTTATACCTTCCCAGGTAACAACCAAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTC',
            '|B.1.2|2020-03-27': 'NNNNNNGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTC'
         }
    def testRetrieveGenome(self):
        by_lineage = \
            {'B.1.1.171':
                [
                    {
                        'covv_virus_name': 'hCoV-19/Canada/Qc-L00240569/2020',
                        'covv_accession_id': 'EPI_ISL_465679', 'covv_collection_date': '2020-03-27',
                        'covv_lineage': 'B.1.1.171',
                        'diffs': [('-', 28, 1)],
                        'missing': [(0, 6)]
                    }
                ],
                'B.1.265':
                    [
                        {
                            'covv_virus_name': 'hCoV-19/Canada/Qc-L00240594/2020',
                            'covv_accession_id': 'EPI_ISL_465680',
                            'covv_collection_date': '2020-03-27',
                            'covv_lineage': 'B.1.265',
                            'diffs': [('~', 30, 'C'), ('~', 32, 'A')],
                            'missing': [(0, 6)]
                        }
                    ],
                'B.1.2':
                    [
                        {
                            'covv_virus_name': 'hCoV-19/Canada/Qc-L00240624/2020',
                            'covv_accession_id': 'EPI_ISL_465681',
                            'covv_collection_date': '2020-03-27',
                            'covv_lineage': 'B.1.2',
                            'diffs': [('+', 13, 'A')],
                            'missing': [(0, 6)]
                        }
                    ]
            }
        result = retrieve_genomes(by_lineage, 'NC_Test.fa')
        self.assertEqual(self.expected, result)


if __name__ == '__main__':
    unittest.main()
