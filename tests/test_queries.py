import collections
import unittest
from ham import ham
from ham import utils
from ham import Gene, ExtantGenome

class QueryTest(unittest.TestCase):

    def setUp(self):
        nwk_path = './tests/data/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/data/simpleEx.orthoxml'
        self.ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')

    def test_get_all_top_level_hogs(self):
        toplevel_hogs = self.ham_analysis.get_dict_top_level_hogs()
        self.assertEqual(len(toplevel_hogs), 3)

    def test_get_all_extant_genes_dict(self):
        self.assertEqual(len(self.ham_analysis.get_dict_extant_genes()), 19)

    def test_get_ancestral_genomes(self):
        self.assertEqual(len(self.ham_analysis.get_list_ancestral_genomes()), 5)

    def test_get_extant_genome_by_name(self):
        query_name = "HUMAN"
        a = self.ham_analysis._get_extant_genome_by_name(name=query_name)
        self.assertIsInstance(a, ExtantGenome)
        self.assertEqual(query_name, a.name)

    def test_get_all_genes_of_hog(self):
        # TODO
        pass

    def test_get_mrca_ancestral_genome_from_genome_set(self):
        human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        mouse = self.ham_analysis._get_extant_genome_by_name(name="MOUSE")

        with self.assertRaises(ValueError):
            self.ham_analysis._get_ancestral_genome_by_mrca_of_genome_set(set([mouse]))

        with self.assertRaises(TypeError):
            self.ham_analysis._get_ancestral_genome_by_mrca_of_genome_set(set([mouse, "human"]))

        mrca_euch = self.ham_analysis._get_ancestral_genome_by_mrca_of_genome_set(set([human, mouse]))
        self.assertEqual("Euarchontoglires", mrca_euch.taxon.name)

        mrca_euch2 = self.ham_analysis._get_ancestral_genome_by_mrca_of_genome_set(set([mrca_euch, mouse]))
        self.assertEqual("Euarchontoglires", mrca_euch2.taxon.name)


if __name__ == "__main__":
    unittest.main()