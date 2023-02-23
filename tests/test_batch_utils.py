import unittest
from argparse import Namespace
from covizu.utils.batch_utils import beadplot_serial, import_labels, make_beadplots, get_mutations

# Build TimeTree - Tests in test_timetree.py


class TestBeadplotSerial(unittest.TestCase):
    def setUp(self):
        self.expected = {
            'lineage': 'B.1.1.171',
            'nodes': {
                'EPI_ISL_465679': [['2020-03-27', 'EPI_ISL_465679', 'hCoV-19/Canada/Qc-L00240569/2020']]
            },
            'edges': [],
            'sampled_variants': 1
        }

    def testBeadplotSerial(self):
        # Lineage with only one variant

        lineage = 'B.1.1.171'
        features = {
            '~|240|T,~|1436|T,~|3036|T,~|3713|T,~|5883|T,~|14407|T,~|20543|T,~|23402|G,~|28880|A,~|28881|A,~|28882|C': [
                {
                    'covv_virus_name': 'hCoV-19/Canada/Qc-L00240569/2020',
                    'covv_accession_id': 'EPI_ISL_465679',
                    'covv_collection_date': '2020-03-27',
                    'covv_lineage': 'B.1.1.171',
                    'covv_location': 'North America / Canada / London',
                    'diffs' : [
                        tuple(['~', 240, 'T']),
                        tuple(['~', 1436, 'T']),
                        tuple(['~', 3036, 'T']),
                        tuple(['~', 3713, 'T']),
                        tuple(['~', 5883, 'T']),
                        tuple(['~', 14407, 'T']),
                        tuple(['~', 20543, 'T']),
                        tuple(['~', 23402, 'G']),
                        tuple(['~', 28880, 'A']),
                        tuple(['~', 28881, 'A']),
                        tuple(['~', 28882, 'C'])
                    ],
                    'missing' : [
                        tuple([0, 6]),
                        tuple([5203, 5222]),
                        tuple([29844, 29903])
                    ]
                }
            ]
        }

        args = Namespace(boot_cuttoff=0.5)
        result = beadplot_serial(lineage, features, args)
        self.assertEqual(self.expected, result)


class TestGetMutations(unittest.TestCase):
    def setUp(self):
        self.expected = {
            'B.1.1.171': [('~', 240, 'T'), ('~', 1436, 'T')],
            'B.1.265': [('~', 240, 'T'), ('~', 1436, 'T')]
        }

    def testGetMutations(self):
        by_lineage = {
            'B.1.1.171': [{
                'covv_virus_name' : '',
                'covv_accession_id' : '',
                'covv_collection_date': '',
                'covv_lineage' : 'B.1.1.171',
                'diffs' : [
                    tuple(['~', 240, 'T']),
                    tuple(['~', 1436, 'T']),
                ],
                'missing' : []
            }],
            'B.1.265': [
                {
                    'covv_virus_name' : '',
                    'covv_accession_id' : '',
                    'covv_collection_date': '',
                    'covv_lineage' : 'B.1.265',
                    'diffs' : [
                        tuple(['~', 240, 'T']),
                        tuple(['~', 240, 'T']),
                        tuple(['~', 240, 'T']),
                        tuple(['~', 1436, 'T']),
                    ],
                    'missing' : []
                },
                {
                    'covv_virus_name' : '',
                    'covv_accession_id' : '',
                    'covv_collection_date': '',
                    'covv_lineage' : 'B.1.1.171',
                    'diffs' : [
                        tuple(['~', 240, 'T']),
                        tuple(['~', 1436, 'T']),
                    ],
                    'missing' : []
                }
            ]
        }

        results = get_mutations(by_lineage)
        self.assertEqual(self.expected, results)


if __name__ == '__main__':
    unittest.main()
