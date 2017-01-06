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
        taxonomy = ham.build_taxonomy(tree_str)
        # And verifying if all tree elements are created
        self.assertEqual(taxonomy.newick_str,
                         "(XENTR, (((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, "
                         "CANFA)Mammalia)Vertebrata;")



        # then clement select his favorite orthoXML file
        orthoxml_path = './tests/simpleEx.orthoxml'
        with open(orthoxml_path, 'r') as orthoxml_file:
            # And build the Abstract_gene objects from file object
            top_hogs, genes = ham.build_hogs_and_genes(orthoxml_file, taxonomy)

        # After he get all the top level hogs
        self.assertEqual(len(top_hogs), 3)
        self.assertEqual(len(genes), 19)

        # Clement is curious to look at the species present within this taxonomy
        extant_genomes = queries.get_extant_genomes(taxonomy)
        # then look if its 6 species are present
        extant_genomes_name = set(ext_genome.name for ext_genome in extant_genomes)
        self.assertEqual(extant_genomes_name, {'RATNO', 'HUMAN', 'CANFA', 'PANTR', 'XENTR', 'MOUSE'})
        self.assertEqual(len(taxonomy.leaves), 6)


        # After he get its freshly built ancestral genomes
        #ancestral_genomes = queries.get_ancestral_genomes(taxonomy)
        # And inspect if they are all present
        #ancestral_genomes_name = set(anc_genome.name for anc_genome in ancestral_genomes)
        #self.assertEqual(ancestral_genomes_name,
                        # {'Vertebrata', 'Primates', 'Mammalia', 'Rodents', 'Euarchontoglires'})
        self.assertEqual(len(taxonomy.internal_nodes), 5)

        '''
        # To finish, he resolve taxonomy problem that can occured
        # due to incomplete lineage name within hogs hierarchy
        ham.resolve_taxonomy_and_hogs() # TODO
        '''

        for node in taxonomy.tree.traverse("postorder"):

            print(node.genome)
            print(node.genome.genes)




if __name__ == "__main__":
    unittest.main()

