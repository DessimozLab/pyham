import shutil
import unittest
from pyham import ham, utils, TreeProfile
import tempfile
import os
from unittest import skip

class TreeProfileTest(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                    use_internal_name=True)
        self.ham_analysis_no_name = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                            use_internal_name=False)

        # gene, retained, duplicated, gain, lost
        self.exp_full = {"Mammalia":[3,2,0,1,0], "Euarchontoglires":[4,2,2,0,0], "Primates":[4,4,0,0,0], "Rodents":[4,4,0,0,0], "Vertebrata":[2,None,None,None,None] }

        self.exp_full_nn = {"HUMAN/PANTR/MOUSE/RATNO/CANFA": [3, 2, 0, 1, 0], "HUMAN/PANTR/MOUSE/RATNO": [4, 2, 2, 0, 0], "HUMAN/PANTR": [4, 4, 0, 0, 0],
                          "MOUSE/RATNO": [4, 4, 0, 0, 0], "XENTR/HUMAN/PANTR/MOUSE/RATNO/CANFA": [2, None, None, None, None],
                            "HUMAN": [4, 3, 0, 1, 1], "PANTR": [4, 4, 0, 0, 0], "MOUSE": [4, 4, 0, 0, 0], "RATNO": [2, 1, 0, 1, 3],
                            "CANFA": [3, 3, 0, 0,0], "XENTR": [2, 2, 0, 0, 0]}

        # gene, retained, duplicated, lost
        self.expected_level_1 = {"Mammalia": [1, 1, 0, 0], "Euarchontoglires": [1, 1, 0, 0], "Primates": [1, 1, 0, 0],
                            "Rodents": [1, 1, 0, 0], "Vertebrata": [1, None, None, None], "HUMAN":[1,1,0,0],
                            "PANTR":[1,1,0,0], "MOUSE":[1,1,0,0], "RATNO":[1,1,0,0], "CANFA":[1,1,0,0], "XENTR":[1,1,0,0]}

        self.expected_level_2 = {"Mammalia": [1, None, None, None], "Euarchontoglires": [1, 1, 0, 0],
                            "Primates": [1, 1, 0, 0],
                            "Rodents": [1, 1, 0, 0], "HUMAN":[1,1,0,0], "PANTR":[1,1,0,0], "MOUSE":[1,1,0,0], "RATNO":[0,0,0,1], "CANFA":[1,1,0,0]}

        self.expected_level_3 = {"Mammalia": [1, 1, 0, 0], "Euarchontoglires": [2, 0, 2, 0],
                            "Primates": [2, 2, 0, 0],
                            "Rodents": [2, 2, 0, 0], "Vertebrata": [1, None, None, None], "HUMAN":[1,1,0,1],
                            "PANTR":[2,2,0,0], "MOUSE":[2,2,0,0], "RATNO":[0,0,0,2], "CANFA":[1,1,0,0], "XENTR":[1,1,0,0]}

    def test_incorrect_hog(self):
        with self.assertRaises(TypeError):
            TreeProfile(self.ham_analysis, hog="")

    def test_tp_hog_at_root(self):

        for hog in self.ham_analysis.get_list_top_level_hogs():
            tp_h = TreeProfile(self.ham_analysis, hog=hog)

            if hog.hog_id == "1":
                exp = self.expected_level_1
            elif hog.hog_id == "2":
                exp = self.expected_level_2
            elif hog.hog_id == "3":
                exp = self.expected_level_3

            for lvl in tp_h.treemap.traverse():
                if lvl.name in exp.keys():
                    observed_level = [lvl.nbr_genes, lvl.retained, lvl.dupl, lvl.lost]
                    self.assertListEqual(exp[lvl.name], observed_level)

    def test_tp_hog_at_level(self):
        mammalia = self.ham_analysis.get_ancestral_genome_by_name("Mammalia")

        for hog in mammalia.genes:

            tp_h = TreeProfile(self.ham_analysis, hog=hog)
            if hog.get_top_level_hog().hog_id == "1":
                exp = self.expected_level_1
            elif hog.get_top_level_hog().hog_id == "2":
                exp = self.expected_level_2
            elif hog.get_top_level_hog().hog_id == "3":
                exp = self.expected_level_3

            for lvl in tp_h.treemap.traverse():
                if lvl.name in exp.keys() and not lvl.is_root():
                    observed_level = [lvl.nbr_genes, lvl.retained, lvl.dupl, lvl.lost]
                    self.assertListEqual(exp[lvl.name], observed_level)

    def test_tp_full(self):
        tp = TreeProfile(self.ham_analysis)
        for lvl in tp.treemap.traverse():
            if lvl.name in self.exp_full.keys():
                observed_level = [lvl.nbr_genes, lvl.retained, lvl.dupl, lvl.gain, lvl.lost]
                self.assertListEqual(self.exp_full[lvl.name], observed_level)

        tp_nn = TreeProfile(self.ham_analysis_no_name)
        for lvl in tp_nn.treemap.traverse():
            if lvl.name in self.exp_full_nn.keys():
                observed_level = [lvl.nbr_genes, lvl.retained, lvl.dupl, lvl.gain, lvl.lost]
                self.assertListEqual(self.exp_full_nn[lvl.name], observed_level)


class ServerBasedTreeProfileTest(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_non_luca_root_hog_works_from_omabrowser(self):
        analysis = ham.Ham(query_database='P53_RAT', use_data_from='oma')

        fn = os.path.join(self.tmpdir, 'tree_profile.html')
        analysis.create_tree_profile(outfile=fn)
        with open(fn, 'rt') as fh:
            html = fh.read()
        self.assertIn('treeData', html)
