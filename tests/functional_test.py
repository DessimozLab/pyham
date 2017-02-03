import unittest

from ham import ham
from ham import utils
import logging


class SetUpHamAnalysis(unittest.TestCase):

    def test_load_taxonomy_from_nwk_file_and_all_hogs_from_orthoxml_file_no_filter(self):

        # load the logger
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

        # Clement select a nwk file as a taxonomy reference
        nwk_path = './tests/simpleEx.nwk'
        # And extract the newick tree as a string
        tree_str = utils.get_newick_string(nwk_path, type="nwk")

        # then clement select his favorite orthoXML file
        orthoxml_path = './tests/simpleEx.orthoxml'

        # Clement create the HAM object that will be the kernel of all analysis
        ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type='orthoxml')

        # And verifying if all tree elements are created
        self.assertEqual(ham_analysis.taxonomy.newick_str,
                         "(XENTR, (((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, "
                         "CANFA)Mammalia)Vertebrata;")

        # After he get all the top level hogs
        self.assertEqual(len(ham_analysis.toplevel_hogs), 3)
        self.assertEqual(len(ham_analysis.extant_gene_map), 19)

        # Clement is curious to look at the species present within this taxonomy
        extant_genomes = ham_analysis.get_extant_genomes()

        # then look if its 6 species are present
        extant_genomes_name = set(ext_genome.name for ext_genome in extant_genomes)
        self.assertEqual(extant_genomes_name, {'RATNO', 'HUMAN', 'CANFA', 'PANTR', 'XENTR', 'MOUSE'})
        self.assertEqual(len(ham_analysis.taxonomy.leaves), 6)
        self.assertEqual(len(ham_analysis.taxonomy.internal_nodes), 5)

        # Clement is happy the parsing went well he celebrates...

    def test_single_hog_analysis(self):
        pass

    def test_ancestral_genome_analysis(self):
        pass

    def test_comparative_analysis(self):
        pass


if __name__ == "__main__":
    unittest.main()

