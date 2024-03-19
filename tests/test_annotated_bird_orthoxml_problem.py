import unittest
from pyham import utils
from pyham import ham
from unittest import skip
import os

class BirdExampleTests(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/birds.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/birds.orthoxml')
        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)
        self.hog = self.ham_analysis.get_list_top_level_hogs()[1]

    def test_parent_of_gene1703732_is_specific_subhog(self):
        self.assertEqual(self.ham_analysis.get_gene_by_id(1703732).parent.hog_id, "HOG:0011319.2b.2b.2b")







if __name__ == "__main__":
    unittest.main()
