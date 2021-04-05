import unittest
from covizu.utils.gisaid_utils import *


class TestBatchFasta(unittest.TestCase):
    def setUp(self):
        self.expected = \
            [('>hCoV-19/Canada/Qc-L00240569/2020\nGGTTTATACCTTCCCAGGTAACAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC\n>hCoV-19/Canada/Qc-L00240594/2020\nATTAAAGGTTTATACCTTCCCAGGTAACAACCAAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC\n',
                [{'covv_virus_name': 'hCoV-19/Canada/Qc-L00240569/2020', 'covv_accession_id': 'EPI_ISL_465679', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.1.171'},
                {'covv_virus_name': 'hCoV-19/Canada/Qc-L00240594/2020', 'covv_accession_id': 'EPI_ISL_465680', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.265'}]
              ),
             ('>hCoV-19/Canada/Qc-L00240624/2020\nATTAAAGGTTTATAACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC\n',
                 [{'covv_virus_name': 'hCoV-19/Canada/Qc-L00240624/2020', 'covv_accession_id': 'EPI_ISL_465681', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.2'}]
             )
            ]

    def testBatchFasta(self):
        gen = \
            [
                {'covv_virus_name':'hCoV-19/Canada/Qc-L00240569/2020',
                 'covv_accession_id': 'EPI_ISL_465679',
                 'covv_collection_date': '2020-03-27',
                 'covv_lineage': 'B.1.1.171',
                 'sequence': 'GGTTTATACCTTCCCAGGTAACAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC'},
                {'covv_virus_name':'hCoV-19/Canada/Qc-L00240594/2020',
                 'covv_accession_id': 'EPI_ISL_465680',
                 'covv_collection_date': '2020-03-27',
                 'covv_lineage': 'B.1.265',
                 'sequence': 'ATTAAAGGTTTATACCTTCCCAGGTAACAACCAAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC'},
                {'covv_virus_name':'hCoV-19/Canada/Qc-L00240624/2020',
                 'covv_accession_id': 'EPI_ISL_465681',
                 'covv_collection_date': '2020-03-27',
                 'covv_lineage': 'B.1.2',
                 'sequence': 'ATTAAAGGTTTATAACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC'}
            ]
        result = list(batch_fasta(gen, 2))
        self.assertEqual(self.expected, result)


class TestExtractFeatures(unittest.TestCase):
    def setUp(self):
        self.batcher = \
            [('>hCoV-19/Canada/Qc-L00240569/2020\nGGTTTATACCTTCCCAGGTAACAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC\n>hCoV-19/Canada/Qc-L00240594/2020\nATTAAAGGTTTATACCTTCCCAGGTAACAACCAAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC\n',
              [{'covv_virus_name': 'hCoV-19/Canada/Qc-L00240569/2020', 'covv_accession_id': 'EPI_ISL_465679', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.1.171'},
               {'covv_virus_name': 'hCoV-19/Canada/Qc-L00240594/2020', 'covv_accession_id': 'EPI_ISL_465680', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.265'}]
              ),
             ('>hCoV-19/Canada/Qc-L00240624/2020\nATTAAAGGTTTATAACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC\n',
              [{'covv_virus_name': 'hCoV-19/Canada/Qc-L00240624/2020', 'covv_accession_id': 'EPI_ISL_465681', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.2'}]
              )
             ]
        self.expected = \
            [
                {'covv_virus_name': 'hCoV-19/Canada/Qc-L00240569/2020', 'covv_accession_id': 'EPI_ISL_465679', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.1.171', 'diffs': [('-', 28, 1)], 'missing': [(0, 6), (75, 29903)]},
                {'covv_virus_name': 'hCoV-19/Canada/Qc-L00240594/2020', 'covv_accession_id': 'EPI_ISL_465680', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.265', 'diffs': [('~', 30, 'C'), ('~', 32, 'A')], 'missing': [(75, 29903)]},
                {'covv_virus_name': 'hCoV-19/Canada/Qc-L00240624/2020', 'covv_accession_id': 'EPI_ISL_465681', 'covv_collection_date': '2020-03-27', 'covv_lineage': 'B.1.2', 'diffs': [('+', 13, 'A')], 'missing': [(75, 29903)]}
            ]

    def testExtractFeatures(self):
        # Delete A in the 28th position, missing ATTAAA, missing 75-29903
        # Substitution in the 30th position C, and 32nd position A
        # Add A in the 13th position
        result = list(extract_features(self.batcher, '/vagrant/covizu/covizu/data/NC_045512.fa', binpath='minimap2', nthread=3, minlen=40))
        self.assertEqual(self.expected, result)

# Test Filter Problematic

if __name__ == '__main__':
    unittest.main()
