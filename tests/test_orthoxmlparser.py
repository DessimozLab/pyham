import collections
import unittest
from ham import ham
from ham import utils

class OrthoXMLParserTest(unittest.TestCase):

    def setUp(self):
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type='orthoxml')
        self.hogs = ham_analysis.get_all_top_level_hogs()
        self.genes = ham_analysis.get_all_extant_genes_dict()

    def test_numberOfGenesPerSpecies(self):
        expected_cnts = dict(HUMAN=4, PANTR=4, MOUSE=4, RATNO=2,
                            CANFA=3, XENTR=2)
        observed_cnts = collections.defaultdict(int)
        for g in self.genes.values():
            observed_cnts[g.genome.name] += 1
        self.assertDictEqual(observed_cnts, expected_cnts)

    def test_scores_on_toplevel(self):
        self.assertEqual(self.hogs["1"].score('consistency'), 1)
        with self.assertRaises(KeyError):
            self.hogs["2"].score('consistency')
        with self.assertRaises(KeyError):
            self.hogs["1"].score('coverage')

    def test_hog_hierarchy(self):
        # TODO
        # test that the hogs reconstruction fit to the schema with all missing taxon, etc..
        pass

if __name__ == "__main__":
    unittest.main()