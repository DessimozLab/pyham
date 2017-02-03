import collections
import unittest
from ham import ham
from ham import utils
from ham import Gene, ExtantGenome

class QueryTest(unittest.TestCase):

    def setUp(self):
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        self.ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type='orthoxml')

    def test_get_all_top_level_hogs(self):
        toplevel_hogs = self.ham_analysis.get_all_top_level_hogs()
        self.assertEqual(len(toplevel_hogs), 3)

    def test_get_all_extant_genes_dict(self):
        self.assertEqual(len(self.ham_analysis.get_all_extant_genes_dict()), 19)

    def test_get_ancestral_genomes(self):
        self.assertEqual(len(self.ham_analysis.get_ancestral_genomes()), 5)

    def test_get_extant_genome_by_name(self):
        query_name = "HUMAN"
        a = self.ham_analysis._get_extant_genome_by_name(name=query_name)
        self.assertIsInstance(a, ExtantGenome)
        self.assertEqual(query_name, a.name)

    def test_get_all_genes_of_hog(self):
        # TODO
        pass

if __name__ == "__main__":
    unittest.main()