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
        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                    use_internal_name=True)

    def test_ancestral_genome(self):
        anc_genome = self.ham_analysis.get_ancestral_genome_by_name("NODE_29")
        self.assertIsNotNone(anc_genome, "Ancestral genome should not be None")
        self.assertEqual(len(anc_genome.genes), 1, "Ancestral genome should have genes")

    def test_fetch_subhog_by_id(self):
        query_hog_id = "HOG:0027790.3c.9ba.5i.8ai_69"
        subhog = self.ham_analysis.get_hog_by_id(query_hog_id)
        self.assertEqual(subhog.hog_id, query_hog_id, "SubHOG should match the queried ID")

    def test_fetch_roothog_by_id(self):
        query_hog_id = "HOG:0027790_45"
        root_hog = self.ham_analysis.get_hog_by_id(query_hog_id)
        self.assertEqual(root_hog.hog_id, query_hog_id, "Root HOG should match the queried ID")
        self.assertIsNone(root_hog.parent, "Root HOG should not have a parent")


class PyHAMInitWithDifferentTaxonomyTests(unittest.TestCase):
    def setUp(self):
        self.nwk_path = os.path.join(os.path.dirname(__file__), './data/tomato.nwk')
        self.tree_str = utils.get_newick_string(self.nwk_path, type="nwk")
        self.orthoxml_path = os.path.join(os.path.dirname(__file__), './data/tom.orthoxml')

    def test_init_with_newick_taxonomy(self):
        analysis = ham.Ham(tree_file=self.tree_str, hog_file=self.orthoxml_path)
        self.assertIsNotNone(analysis.taxonomy, "Taxonomy should be initialized")

    def test_init_with_orthoxml_taxonomy(self):
        analysis = ham.Ham(hog_file=self.orthoxml_path)
        self.assertIsNotNone(analysis.taxonomy, "Taxonomy should be initialized with internal names")

    def test_taxonomy_from_newick_identical(self):
        analysis_newick = ham.Ham(tree_file=self.tree_str, hog_file=self.orthoxml_path)
        analysis_orthoxml = ham.Ham(hog_file=self.orthoxml_path)
        self.assertEqual(analysis_newick.taxonomy, analysis_orthoxml.taxonomy,
                         "Taxonomy from Newick and Orthoxml should be identical")


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
