import unittest

from ham import ham
from ham import queries
from ham import utils

class SetUpHamAnalysis(unittest.TestCase):

    def test_load_taxonomy_from_nwk_file_and_all_hogs_from_orthoxml_file_no_filter(self):


        # Clement select a nwk file as a taxonomy reference
        nwk_path = './tests/simpleEx.nwk'
        # And extract the newick tree as a string
        tree_str = utils.get_newick_string(nwk_path, type="nwk")

        # Then he build the taxonomy objects induced by the newick tree.
        taxonomy = ham.build_taxonomy_and_ancestral_genomes(tree_str)
        # And verifying if all tree elements are created
        self.assertEqual(taxonomy.newick_str, "(XENTR, (((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, CANFA)Mammalia)Vertebrata;")
        self.assertEqual(len(taxonomy.internal_nodes), 5 )
        self.assertEqual(len(taxonomy.leaves), 6 )

        # Clement is curious to look at the species present within this taxonomy
        extant_genomes = queries.get_extant_genomes(taxonomy)
        # then look if its 6 species are present
        extant_genomes_name = set(ext_genome.taxon.name for ext_genome in extant_genomes)
        self.assertEqual(extant_genomes_name, {'RATNO', 'HUMAN', 'CANFA', 'PANTR', 'XENTR', 'MOUSE'})

        # After he get its freshly built ancestral genomes
        ancestral_genomes = queries.get_ancestral_genomes(taxonomy)
        # And inspect if they are all present
        ancestral_genomes_name = set(anc_genome.taxon.name for anc_genome in ancestral_genomes)
        self.assertEqual(ancestral_genomes_name, {'Vertebrata', 'Primates', 'Mammalia', 'Rodents', 'Euarchontoglires'})

        # then clement select his favorite orthoXML file
        orthoxml_path = './tests/simpleEx.orthoxml'
        orthoxml_file = open(orthoxml_path, 'r')

        # And build the Abstract_gene objects from file object
        hogs, genes = ham.build_hogs_and_genes(orthoxml_file, taxonomy)
        orthoxml_file.close()
        # After he get all the top level hogs
        top_hogs = queries.get_top_level_hogs(hogs)
        self.assertEqual(len(top_hogs),3 )
        self.assertEqual(len(genes),17 )

        # To finish, he resolve taxonomy problem that can occured
        # due to incomplete lineage name within hogs hierarchy
        ham.resolve_taxonomy_and_hogs() # TODO


if __name__ == "__main__":
    unittest.main()

