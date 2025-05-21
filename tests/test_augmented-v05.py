import unittest
from pyham import utils
from pyham import ham
import os

class AugmentedExampleTests(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/augmented-test-v05.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/augmented-test-v05.orthoxml')
        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)

    def test_all_hogs_at_Eukaryota_have_taxid_2759_in_hogid(self):
        hogs_with_wrong_taxid = [hog for hog in self.ham_analysis.get_ancestral_genome_by_name("Eukaryota").genes if hog.hog_id.split('_')[1] != "2759"]
        self.assertEqual(hogs_with_wrong_taxid, [], "All HOGs at Eukaryota should have taxid 2759 in hog_id")


if __name__ == "__main__":
    unittest.main()
