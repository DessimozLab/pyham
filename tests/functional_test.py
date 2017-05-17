import unittest
from pyham import ham
from pyham import utils
import logging
import os

def _str_array(array):
    array_converted = []
    for e in array:
        array_converted.append(str(e))
    return set(array_converted)


def _str_dict_one_value(dict):
    for kk in dict.keys():
        dict[str(kk)] = dict.pop(kk)
    for k, v in dict.items():
        dict[k] = str(v)
    return dict


def _str_dict_array_value(dict):
    for kk in dict.keys():
        dict[str(kk)] = dict.pop(kk)
    for k, vs in dict.items():
        array = []
        for v in vs:
            array.append(str(v))
        dict[k] = set(array)
    return dict


class HamAnalysis(unittest.TestCase):

    def test_load_taxonomy_from_nwk_file_and_all_hogs_from_orthoxml_file_no_filter(self):

        # load the logger
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

        # Clement select a nwk file as a taxonomy reference
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        # And extract the newick tree as a string
        tree_str = utils.get_newick_string(nwk_path, type="nwk")

        # then clement select his favorite orthoXML file
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        # Clement create the Ham object that will be the kernel of all analysis
        ham_analysis = ham.Ham(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')

        # And verifying if all tree elements are created
        self.assertEqual(ham_analysis.taxonomy.newick_str,
                         "(XENTR, (((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, "
                         "CANFA)Mammalia)Vertebrata;")

        # After he get all the top level hogs
        self.assertEqual(len(ham_analysis.top_level_hogs), 3)
        self.assertEqual(len(ham_analysis.extant_gene_map), 19)

        # Clement is curious to look at the species present within this taxonomy
        extant_genomes = ham_analysis.get_list_extant_genomes()


        # then look if its 6 species are present
        extant_genomes_name = set(ext_genome.name for ext_genome in extant_genomes)
        self.assertEqual(extant_genomes_name, {'RATNO', 'HUMAN', 'CANFA', 'PANTR', 'XENTR', 'MOUSE'})
        self.assertEqual(len(ham_analysis.taxonomy.leaves), 6)
        self.assertEqual(len(ham_analysis.taxonomy.internal_nodes), 5)

        # Clement is happy the parsing went well, now he celebrates...

    def test_single_hog_analysis(self):
        pass

    def test_ancestral_genome_analysis(self):
        pass

    def test_lineage_comparative_analysis(self):

        # Clement initialise the pyham analyzer objet as it's explained in the documentation
        logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        ham_analysis = ham.Ham(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)

        # Then clement is interest to investigate on what happened between the ancestral genomes of vertebrates
        # and the extent genomes of the mouse.

        # First he select the related genomes objectt via their name or mrca.
        mouse = ham_analysis.get_extant_genome_by_name(name="MOUSE")
        frog = ham_analysis.get_extant_genome_by_name(name="XENTR")
        vertebrates = ham_analysis.get_ancestral_genome_by_mrca_of_genome_set({mouse, frog})

        # Then, he compare the two genomes of interest
        vertical_map_mouse_vs_vert = ham_analysis.compare_genomes_vertically(mouse, vertebrates)

        # Now he is interest by the HOG that have stay single copy between these two levels
        self.assertDictEqual({'<HOG(1)>': 'Gene(31)'}, _str_dict_one_value(vertical_map_mouse_vs_vert.get_identical()))

        # ... and at the duplicated genes
        self.assertDictEqual({'<HOG(3)>': {'Gene(34)', 'Gene(33)'}}, _str_dict_array_value(vertical_map_mouse_vs_vert.get_duplicated()))

        # Clement is curious and want to look if there is gene that have been lost
        self.assertSetEqual(set(), _str_array(vertical_map_mouse_vs_vert.get_lost()))

    def test_load_taxonomy_from_nwk_file_and_from_orthoxml_file_with_filter_hog2(self):

        # load the logger
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

        # Clement select a nwk file as a taxonomy reference
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        # And extract the newick tree as a string
        tree_str = utils.get_newick_string(nwk_path, type="nwk")

        # Clement will now setup the filter object
        f = ham.ParserFilter()
        f.add_hogs_via_hogId([2])

        # Clement check that the filter contained all information
        self.assertEqual(set(f.HOGId_filter), {'2'})
        self.assertEqual(set(f.GeneExtId_filter), set())
        self.assertEqual(set(f.GeneIntId_filter), set())

        # then clement select his favorite orthoXML file
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        # Clement create the Ham object that will be the kernel of all analysis
        ham_analysis = ham.Ham(tree_str, orthoxml_path, filter_object=f, use_internal_name=True)
        self.assertEqual(f, ham_analysis.filter_obj)

        # Clement check that what the filter understood was good
        self.assertSetEqual(set(ham_analysis.filter_obj.geneUniqueId), {'2', '32', '22', '12'})
        self.assertSetEqual(set(ham_analysis.filter_obj.hogsId), {'2'})

        # Clement check that the parsed informatio is correct
        self.assertEqual(len(ham_analysis.top_level_hogs), 1)
        self.assertEqual(len(ham_analysis.extant_gene_map), 4)

    def test_treeProfile_on_hog3(self):

        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')
        ham_analysis = ham.Ham(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)

        # Clement get the toplevel hog and build the tree profile on them
        #for hog in ham_analysis.get_list_top_level_hogs():
            #ham_analysis.create_tree_profile(hog, outfile="./hog{}.png".format(hog.hog_id), export_with_histogram=True)

    def test_treeProfile_on_full_setup(self):

        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')
        ham_analysis = ham.Ham(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')

        # Clement run the create_tree_profile for the whole genomic setup
        #ham_analysis.create_tree_profile(outfile="./tp.png", export_with_histogram=True)
        #ham_analysis.create_tree_profile(outfile="./tp2.png", export_with_histogram=True)



if __name__ == "__main__":
    unittest.main()

