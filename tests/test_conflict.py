import unittest
from pyham import utils
from pyham import ham
from unittest import skip
import os

@skip
class HAMTestConflict1(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/Conflict1/tree.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/Conflict1/hog.orthoxml')
        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)
        self.hog = self.ham_analysis.get_list_top_level_hogs()[0]

    def test_only_one_hog(self):
        self.assertEqual(len(self.ham_analysis.get_list_top_level_hogs()), 1)

    def test_root(self):
        self.assertEqual(self.hog.genome.name, "Chordata")

    def test_hog_structure(self):
        self.assertEqual(len(self.hog.children), 2)

        hogs_Euteleostomi = []
        for j in self.hog.get_all_descendant_hogs():
            if j.genome.name == "Euteleostomi":
                hogs_Euteleostomi.append(j)

        self.assertEqual(len(hogs_Euteleostomi), 1)

        self.assertEqual(len(hogs_Euteleostomi[0].duplications), 1)

        self.assertEqual(len(hogs_Euteleostomi[0].duplications[0].children), 2)








if __name__ == "__main__":
    unittest.main()
