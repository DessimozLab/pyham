import unittest
from pyham import utils
from pyham import ham
from pathlib import Path
import os

class AncestralTomatoGenomeContentTest(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/tomato.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/tom.orthoxml')
        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)

    def test_ancestral_genome(self):
        anc_genome = self.ham_analysis.get_ancestral_genome_by_name("NODE_29")
        self.assertIsNotNone(anc_genome, "Ancestral genome should not be None")
        self.assertEqual(len(anc_genome.genes), 1, "Ancestral genome should have genes")


class AncestralMetazoaGenomeContentTest(unittest.TestCase):
    def setUp(self):
        data = Path(__file__).parent / 'data'
        nwk_path = data / 'euk-HOG0001687.nwk'
        orthoxml_path = data / 'euk-HOG0001687.orthoxml'
        self.ham_analysis = ham.Ham(tree_file=nwk_path, hog_file=orthoxml_path,
                                    tree_format="newick", use_internal_name=True)

    def test_ancestral_genome_empty_at_treeroot(self):
        anc_genome = self.ham_analysis.get_ancestral_genome_by_name("internal_0")
        self.assertEqual(len(anc_genome.genes), 0, "Ancestral genome at root should be empty")

    def test_ancestral_genome_at_hog_root(self):
        anc_genome = self.ham_analysis.get_ancestral_genome_by_name("42289")
        self.assertEqual(len(anc_genome.genes), 1, "Ancestral genome at hog have one HOG")
        self.assertIsNone(anc_genome.genes[0].parent)


if __name__ == "__main__":
    unittest.main()