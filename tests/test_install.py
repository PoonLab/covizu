import unittest
import os
import covizu


class TestInstalledFiles(unittest.TestCase):
    def test_has_vcf_file(self):
        path = os.path.join(covizu.__path__[0], 'data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf')
        self.assertTrue(os.path.exists(path))
