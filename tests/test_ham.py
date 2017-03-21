import unittest
from ham import utils
from ham import ham

class HAMTest(unittest.TestCase):

    def setUp(self):
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        self.ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type='orthoxml')
        self.hogs = self.ham_analysis.get_all_top_level_hogs()
        self.genes = self.ham_analysis.get_all_extant_genes_dict()

        self.human = self.ham_analysis.get_extant_genome_by_name(name="HUMAN")
        self.frog = self.ham_analysis.get_extant_genome_by_name(name="XENTR")
        self.mouse = self.ham_analysis.get_extant_genome_by_name(name="MOUSE")
        self.rat = self.ham_analysis.get_extant_genome_by_name(name="RATNO")
        self.chimp = self.ham_analysis.get_extant_genome_by_name(name="PANTR")
        self.vertebrates = self.ham_analysis.get_mrca_ancestral_genome_from_genome_set({self.human, self.frog})

    def test_compare_level_correct_input(self):

        # launch wrong analysis
        with self.assertRaises(TypeError):
            self.ham_analysis.compare_genomes({self.human,self.chimp}, "xxx")

        # launch vertical analysis with more than 2 genomes
        with self.assertRaises(TypeError):
            self.ham_analysis.compare_genomes({self.human,self.chimp, self.rat}, "vertical")

        # launch lateral analysis with 1 genomes
        with self.assertRaises(TypeError):
            self.ham_analysis.compare_genomes({self.human}, "lateral")

        # launch lateral analysis with 0 genomes
        with self.assertRaises(TypeError):
            self.ham_analysis.compare_genomes(set(), "lateral")


if __name__ == "__main__":
    unittest.main()
