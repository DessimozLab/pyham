import unittest
from pyham import ham, utils, iham
import os

class IHAMTest(unittest.TestCase):

    def setUp(self):

        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                    use_internal_name=True)

    def test_ok(self):

        hog_3 = self.ham_analysis.get_hog_by_id(1)
        self.ham_analysis.create_iHam(hog=hog_3)


