import unittest

from ham import ham

class SetUpHamAnalysis(unittest.TestCase):

    def test_load_taxonomy_from_nwk_file_and_all_hogs_from_orthoxml_file_no_filter(self):

        # Clement select a nwk file as a taxonomy reference
        nwk_path = './tests/simpleEx.nwk'

        # And extract the newick tree as a string
        nwk_file = open(nwk_path, 'r')
        tree_str = nwk_file.read()
        nwk_file.close()

        # Then he build the taxonomy/ancestral genomes objects induced by the newick tree
        taxonomy, ancestral_genomes = ham.build_taxonomy_and_ancestral_genomes(tree_str)

        # then clement select his favorite orthoXML file
        orthoxml_path = './tests/simpleEx.orthoxml'
        orthoxml_file = open(orthoxml_path, 'r')

        # And build the Abstract_gene objects from file object
        hogs, genes = ham.build_hogs_and_genes(orthoxml_file)
        orthoxml_file.close()

        # To finish, he resolve taxonomy problem that can occured
        # due to incomplete lineage name within hogs hierarchy
        ham.resolve_taxonomy_and_hogs(taxonomy, hogs, genes)


if __name__ == "__main__":
    unittest.main()

