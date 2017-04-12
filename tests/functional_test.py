import unittest
from unittest import skip
from ham import ham
from ham import utils
import logging

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
        nwk_path = './tests/simpleEx.nwk'
        # And extract the newick tree as a string
        tree_str = utils.get_newick_string(nwk_path, type="nwk")

        # then clement select his favorite orthoXML file
        orthoxml_path = './tests/simpleEx.orthoxml'

        # Clement create the HAM object that will be the kernel of all analysis
        ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')

        # And verifying if all tree elements are created
        self.assertEqual(ham_analysis.taxonomy.newick_str,
                         "(XENTR, (((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, "
                         "CANFA)Mammalia)Vertebrata;")

        # After he get all the top level hogs
        self.assertEqual(len(ham_analysis.toplevel_hogs), 3)
        self.assertEqual(len(ham_analysis.extant_gene_map), 19)

        # Clement is curious to look at the species present within this taxonomy
        extant_genomes = ham_analysis.get_all_extant_genomes()


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

        # Clement initialise the ham analyzer objet as it's explained in the documentation
        logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')

        # Then clement is interest to investigate on what happened between the ancestral genomes of vertebrates
        # and the extent genomes of the mouse.

        # First he select the related genomes objectt via their name or mrca.
        mouse = ham_analysis.get_extant_genome_by_name(name="MOUSE")
        frog = ham_analysis.get_extant_genome_by_name(name="XENTR")
        vertebrates = ham_analysis.get_mrca_ancestral_genome_from_genome_set({mouse, frog})

        # Then, he compare the two genomes of interest
        vertical_map_mouse_vs_vert = ham_analysis.compare_genomes({mouse, vertebrates}, analysis='vertical')

        # Now he is interest by the HOG that have stay single copy between these two levels
        self.assertDictEqual({'<HOG(1)>': 'Gene(31)'}, _str_dict_one_value(vertical_map_mouse_vs_vert.get_single()))

        # ... and at the duplicated genes
        self.assertDictEqual({'<HOG(3)>': {'Gene(34)', 'Gene(33)'}}, _str_dict_array_value(vertical_map_mouse_vs_vert.get_duplicated()))

        # Clement is curious and want to look if there is gene that have been lost
        self.assertSetEqual(set(), _str_array(vertical_map_mouse_vs_vert.get_lost()))

    def test_multi_genome_comparative_analysis(self):

        # Clement initialise the ham analyzer objet as it's explained in the documentation
        logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')

        '''
        # ... UNDER CONSTRUCTION ... #

        # Then clement is interest to investigate on the difference/similarities that exists between the primates and the rodents.

        # First he select the related genomes object via their name and mrca.
        mouse = ham_analysis.get_extant_genome_by_name(name="MOUSE")
        rat = ham_analysis.get_extant_genome_by_name(name="RATNO")
        human = ham_analysis.get_extant_genome_by_name(name="HUMAN")
        chimp = ham_analysis.get_extant_genome_by_name(name="PANTR")

        rodents = ham_analysis.get_mrca_ancestral_genome_from_genome_set({mouse, rat})
        primates = ham_analysis.get_mrca_ancestral_genome_from_genome_set({human, chimp})

        # Then, he compare the two genomes of interest through their mrca (Euarchontoglires)
        lateral_map_rodents_vs_primates = ham_analysis.compare_genomes({rodents, primates}, analysis='lateral')

        # Now he is interest by the HOG that have stay single in the both taxon
        s = lateral_map_rodents_vs_primates.get_single()
        print(lateral_map_rodents_vs_primates.ancestor.taxon)
        for h in s:
            print(h)
            print("\t {} ".format(s[h]))

        # ... and at the duplicated genes

        # Clement is curious and want to look if there is gene that have been lost

        # ... UNDER CONSTRUCTION ... #

        '''

    def test_load_taxonomy_from_nwk_file_and_from_orthoxml_file_with_filter_hog2(self):

        # load the logger
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

        # Clement select a nwk file as a taxonomy reference
        nwk_path = './tests/simpleEx.nwk'
        # And extract the newick tree as a string
        tree_str = utils.get_newick_string(nwk_path, type="nwk")

        # Clement will now setup the filter object
        f = ham.ParserFilter()
        f.add_hogs_via_hogId(["2"])

        # Clement check that the filter contained all information
        self.assertEqual(set(f.HOGId_filter), {"2"})
        self.assertEqual(set(f.GeneExtId_filter), set())
        self.assertEqual(set(f.GeneIntId_filter), set())

        # then clement select his favorite orthoXML file
        orthoxml_path = './tests/simpleEx.orthoxml'

        # Clement create the HAM object that will be the kernel of all analysis
        ham_analysis = ham.HAM(tree_str, orthoxml_path, filterObject=f)
        self.assertEqual(f, ham_analysis.filterObj)

        # Clement check that what the filter understood was good
        self.assertSetEqual(set(ham_analysis.filterObj.geneUniqueId), {'2', '32', '22', '12'})
        self.assertSetEqual(set(ham_analysis.filterObj.hogsId), {'2'})

        # Clement check that the parsed informatio is correct
        self.assertEqual(len(ham_analysis.toplevel_hogs), 1)
        self.assertEqual(len(ham_analysis.extant_gene_map), 4)

    def test_treeProfile_on_hog3(self):

        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')


        # Clement get the hog2 of interest
        hog3 = ham_analysis.get_hog_by_id("3")
        self.assertEqual(str(hog3), "<HOG(3)>")

        # and run the treeProfile on it
        tp_hog3 = ham_analysis.treeProfile(hog3)
        #print(tp_hog3.dirty_display())


        '''
        for hog in ham_analysis.get_all_top_level_hogs().values():
            tp_hog = ham_analysis.treeProfile(hog)
            print(hog)
            print(tp_hog.dirty_display())

        '''

    def test_treeProfile_on_full_setup(self):

        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')

        # Clement run the treeProfile for the whole genomic setup
        tp = ham_analysis.treeProfile(outfile="./tp.png", export_with_histogram=True)
        #tp.dirty_display()



if __name__ == "__main__":
    unittest.main()

