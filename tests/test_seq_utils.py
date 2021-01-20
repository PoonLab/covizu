import unittest
from covizu.utils.seq_utils import QPois


class TestQPois(unittest.TestCase):
    def setUp(self):
        self.qp = QPois(quantile=1-0.005, rate=0.0655, maxtime=1e3, origin='2019-12-01')

    def testIsOutlier(self):
        coldate = '2020-03-10'
        ndiffs = 10
        result = self.qp.is_outlier(coldate, ndiffs)
        self.assertEqual(result, False)

    def testIsOutlierManyDiff(self):
        coldate = '2020-03-10'
        ndiffs = 500
        result = self.qp.is_outlier(coldate, ndiffs)
        self.assertEqual(result, True)

    def testIsOutlierFewDiff(self):
        coldate = '2020-03-10'
        ndiffs = 1
        result = self.qp.is_outlier(coldate, ndiffs)
        self.assertEqual(result, True)

if __name__ == '__main__':
    unittest.main()
