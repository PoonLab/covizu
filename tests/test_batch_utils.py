import unittest
from argparse import Namespace
from covizu.utils.batch_utils import beadplot_serial, import_labels, make_beadplots


class TestBeadplotSerial(unittest.TestCase):
    def setUp(self):
        self.expected = {
            'lineage' : 'B.1.1.171',
            'nodes' : {
                'EPI_ISL_465679' :
                [{
                    'accession' : 'EPI_ISL_465679',
                    'label1' : 'hCoV-19/Canada/Qc-L00240569/2020',
                    'country' : 'Canada',
                    'coldate' : '2020-03-27'
                }]
            },
            'edges' : []
        }

    def testBeadplotSerial(self):
        # Lineage with only one variant

        lineage = 'B.1.1.171'
        features = [
            {
                'covv_virus_name' : 'hCoV-19/Canada/Qc-L00240569/2020',
                'covv_accession_id' : 'EPI_ISL_465679',
                'covv_collection_date': '2020-03-27',
                'covv_lineage' : 'B.1.1.171',
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
        args = Namespace(boot_cuttoff=0.5)
        result = beadplot_serial(lineage, features, args)
        self.assertEqual(self.expected, result)


if __name__ == '__main__':
    unittest.main()
