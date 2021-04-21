import unittest
from pyham import utils
from pyham import ham
import ete3
import os


# This helps to convert elements of list/dictionary to string in order to make easier assertEqual test.
def _str_dict_one_value(dict):
    return {str(key): str(val) for key, val in dict.items()}


def _str_array(array):
    return set(str(e) for e in array)


class HAMTestSetUp(unittest.TestCase):

    def setUp(self):
        self.nwk_str_empty = ""
        self.nwk_str_wrong = "(A,B)x0.4"
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        self.nwk_str = utils.get_newick_string(nwk_path, type="nwk")

        self.phyloxml_file = os.path.join(os.path.dirname(__file__), './data/simpleEx.phyloxml')

        self.orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        with open(self.orthoxml_path, 'r') as orthoxml_file:
            data = orthoxml_file.read()
            self.orthoxml_string = data

    def test_load_from_server(self):

        # Not valid use_data_from
        with self.assertRaises(TypeError):
            ham.Ham(use_data_from='xxx')

        #  Valid use_data_from but missing query_database
        with self.assertRaises(TypeError):
            ham.Ham(use_data_from='oma')

    '''
    def test_load_from_oma(self):

        ham.Ham(query_database="HUMAN2", use_data_from='oma')

        with self.assertRaises(TypeError):
            ham.Ham(query_database = "1bshwls", use_data_from = 'oma')

        with self.assertRaises(TypeError):
            ham.Ham(use_data_from = 'oma')
    '''

    def test_wrong_newick_str(self):

        with self.assertRaises(ete3.parser.newick.NewickError):
            ham.Ham(self.nwk_str_empty, self.orthoxml_path, type_hog_file='orthoxml')

        with self.assertRaises(ete3.parser.newick.NewickError):
            ham.Ham(self.nwk_str_wrong, self.orthoxml_path, type_hog_file='orthoxml')

    def test_phyloxml_tag(self):

        with self.assertRaises(TypeError):
            ham.Ham(self.phyloxml_file, self.orthoxml_path, type_hog_file='orthoxml', tree_format='phyloxml', phyloxml_leaf_name_tag='None')

        ham.Ham(self.phyloxml_file, self.orthoxml_path, type_hog_file='orthoxml', tree_format='phyloxml',  phyloxml_leaf_name_tag='clade_name')

    def test_wrong_type_hog_file(self):

        with self.assertRaises(TypeError):
            ham.Ham(self.nwk_str, self.orthoxml_path, type_hog_file='xs')

        with self.assertRaises(TypeError):
            ham.Ham(self.nwk_str, self.orthoxml_path, type_hog_file='')

        with self.assertRaises(TypeError):
            ham.Ham(self.nwk_str, self.orthoxml_path, type_hog_file=None)

    def test_wrong_filter(self):
        with self.assertRaises(TypeError):
            ham.Ham(self.nwk_str, self.orthoxml_path, filter_object="x")

    def test_orthoxml_as_string(self):

        self.ham_analysis = ham.Ham(tree_file=self.nwk_str, hog_file=self.orthoxml_string, type_hog_file='orthoxml', orthoXML_as_string = True)

        with self.assertRaises(IOError):
            self.ham_analysis = ham.Ham(tree_file=self.nwk_str, hog_file=self.orthoxml_string, type_hog_file='orthoxml')


class HAMTest(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")

        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')
        self.hogs = self.ham_analysis.get_dict_top_level_hogs()
        self.genes = self.ham_analysis.get_dict_extant_genes()

        self.human = self.ham_analysis.get_extant_genome_by_name(name="HUMAN")
        self.frog = self.ham_analysis.get_extant_genome_by_name(name="XENTR")
        self.mouse = self.ham_analysis.get_extant_genome_by_name(name="MOUSE")
        self.rat = self.ham_analysis.get_extant_genome_by_name(name="RATNO")
        self.chimp = self.ham_analysis.get_extant_genome_by_name(name="PANTR")
        self.vertebrates = self.ham_analysis.get_ancestral_genome_by_mrca_of_genome_set({self.human, self.frog})

    def test_compare_level_correct_input_vertical(self):

        # launch vertical analysis with incorrect genomes
        with self.assertRaises(TypeError):
            self.ham_analysis.compare_genomes_vertically(None, self.chimp)

    def test_compare_level_correct_input_lateral(self):

        # launch lateral analysis analysis with incorrect genomes
        with self.assertRaises(TypeError):
            self.ham_analysis.compare_genomes_lateral(None, self.chimp)


class HAMTestQuery(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        nwk_str = utils.get_newick_string(nwk_path, type="nwk")

        nwk_path_no_name = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        nwk_str_no_name = utils.get_newick_string(nwk_path_no_name, type="nwk")

        phyloxml_file = os.path.join(os.path.dirname(__file__), './data/simpleEx.phyloxml')

        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        with open(orthoxml_path, 'r') as orthoxml_file:
            data = orthoxml_file.read()
            self.orthoxml_string = data

        # using newick with name on both internal nodes and leaves
        self.h = ham.Ham(nwk_str, orthoxml_path, use_internal_name=True)

        # using newick with name only at leaves
        self.hn = ham.Ham(nwk_str_no_name, orthoxml_path)

        # using phyloxml file with name
        self.hpx = ham.Ham(phyloxml_file, orthoxml_path,  use_internal_name=True, tree_format='phyloxml', phyloxml_internal_name_tag='taxonomy_scientific_name', phyloxml_leaf_name_tag='taxonomy_code' )

        # using newick with name on both internal nodes and leaves and filter for HOG2
        self.filter_genome = {"HUMAN", "MOUSE", "CANFA", "PANTR"}
        self.filter_genes = {'2', '32', '22', '12'}
        self.filter_genes_ext = {'HUMAN2', 'MOUSE2', 'CANFA2', 'PANTR2'}
        self.filter_hogs = {'2'}
        self.no_filter_genome = {"XENTR", "RATNO"}
        self.no_filter_genes = {'1', '11', '21', '31', '41', '51', '3', '13', '23', '33', '53', '34', '14'}
        self.no_filter_genes_ext = {'HUMAN1', 'PANTR1', 'CANFA1', 'MOUSE1', 'RATNO1', 'XENTR1', 'HUMAN3',
                                    'PANTR3', 'CANFA3', 'MOUSE3', 'XENTR3', 'MOUSE4', 'PANTR4'}
        self.no_filter_hogs = {'1','3'}
        f = ham.ParserFilter()
        f.add_hogs_via_hogId([2])
        self.hf = ham.Ham(nwk_str, orthoxml_path, filter_object=f, use_internal_name=True)

        # test that filter work with string
        self.hstring = ham.Ham(nwk_str, self.orthoxml_string, use_internal_name=True,
                                orthoXML_as_string=True)

        # test that filter work with string
        self.hfstring = ham.Ham(nwk_str, self.orthoxml_string, filter_object=f, use_internal_name=True, orthoXML_as_string = True)

    # HOG

    def test_get_hog_by_id(self):

        # If Id not exist
        with self.assertRaises(KeyError):
            self.h.get_hog_by_id("d")

        # Get Hog with str(id)
        hog3 = self.h.get_hog_by_id("3")
        self.assertEqual(str(hog3), "<HOG(3)>")

        # Get Hog with int(id)
        hog3 = self.h.get_hog_by_id(3)
        self.assertEqual(str(hog3), "<HOG(3)>")

        hog3 = self.hstring.get_hog_by_id(3)
        self.assertEqual(str(hog3), "<HOG(3)>")

        hog3 = self.hpx.get_hog_by_id(3)
        self.assertEqual(str(hog3), "<HOG(3)>")


        ###############
        # With filter #
        ###############

        for hog_id in self.filter_hogs:
            hog = self.hf.get_hog_by_id(hog_id)
            self.assertEqual(str(hog), "<HOG({})>".format(hog_id))

            hog = self.hfstring.get_hog_by_id(hog_id)
            self.assertEqual(str(hog), "<HOG({})>".format(hog_id))

        for hog_id in self.no_filter_hogs:
            with self.assertRaises(KeyError):
                self.hf.get_hog_by_id(hog_id)

            with self.assertRaises(KeyError):
                self.hfstring.get_hog_by_id(hog_id)

    def test_get_hog_by_gene(self):

        gene3 = self.h.get_gene_by_id("3")

        # Get Hog with Gene object
        hog3 = self.h.get_hog_by_gene(gene3)
        self.assertEqual(str(hog3), "<HOG(3)>")

        hog3 = self.hpx.get_hog_by_gene(gene3)
        self.assertEqual(str(hog3), "<HOG(3)>")

        # If gene argument is not a Gene object
        with self.assertRaises(KeyError):
            self.h.get_hog_by_gene("gene")

        ###############
        # With filter #
        ###############

        for gene_id in self.filter_genes:
            gene = self.hf.get_gene_by_id(gene_id)
            hog = self.hf.get_hog_by_gene(gene)
            self.assertEqual(str(hog), "<HOG({})>".format(hog.hog_id))

        for gene_id in self.no_filter_genes:
            with self.assertRaises(KeyError):
                gene = self.hf.get_gene_by_id(gene_id)
                self.hf.get_hog_by_gene(gene)

    def test_get_list_top_level_hogs(self):
        hogs = self.h.get_list_top_level_hogs()
        self.assertSetEqual(_str_array(hogs), {"<HOG(2)>", "<HOG(1)>", "<HOG(3)>"})

        hogs = self.hpx.get_list_top_level_hogs()
        self.assertSetEqual(_str_array(hogs), {"<HOG(2)>", "<HOG(1)>", "<HOG(3)>"})

        ###############
        # With filter #
        ###############

        hogs = self.hf.get_list_top_level_hogs()
        self.assertSetEqual(_str_array(hogs), {"<HOG(2)>"})

    def test_get_dict_top_level_hogs(self):
        hogs = self.h.get_dict_top_level_hogs()
        self.assertDictEqual(_str_dict_one_value(hogs), {"2": "<HOG(2)>", "1": "<HOG(1)>", "3": "<HOG(3)>"})

        hogs = self.hpx.get_dict_top_level_hogs()
        self.assertDictEqual(_str_dict_one_value(hogs), {"2": "<HOG(2)>", "1": "<HOG(1)>", "3": "<HOG(3)>"})

        ###############
        # With filter #
        ###############

        hogs = self.hf.get_dict_top_level_hogs()
        self.assertDictEqual(_str_dict_one_value(hogs), {"2": "<HOG(2)>"})

        hogs = self.hfstring.get_dict_top_level_hogs()
        self.assertDictEqual(_str_dict_one_value(hogs), {"2": "<HOG(2)>"})

    # Gene

    def test_get_gene_by_id(self):

        # If Id not exist
        with self.assertRaises(KeyError):
            self.h.get_gene_by_id("d")

        # Get Gene with str(id)
        gene3 = self.h.get_gene_by_id("3")
        self.assertEqual(str(gene3), "Gene(3)")

        gene3 = self.hpx.get_gene_by_id("3")
        self.assertEqual(str(gene3), "Gene(3)")

        # Get Gene with int(id)
        gene3 = self.h.get_gene_by_id(3)
        self.assertEqual(str(gene3), "Gene(3)")

        ###############
        # With filter #
        ###############

        for gene_id in self.filter_hogs:
            gene = self.hf.get_gene_by_id(gene_id)
            self.assertEqual(str(gene), "Gene({})".format(gene_id))

        for gene_id in self.no_filter_genes:
            with self.assertRaises(KeyError):
                self.hf.get_gene_by_id(gene_id)

    def test_get_genes_by_external_id(self):

        # If Id not exist
        with self.assertRaises(KeyError):
            self.h.get_genes_by_external_id("d")

        # Get Gene with protId
        gene12 = self.h.get_genes_by_external_id("PANTR2")[0]
        self.assertEqual(str(gene12), "Gene(12)")

        # Get Gene with geneId
        gene12 = self.h.get_genes_by_external_id("PANTRg2")[0]
        self.assertEqual(str(gene12), "Gene(12)")

        gene12 = self.hstring.get_genes_by_external_id("PANTRg2")[0]
        self.assertEqual(str(gene12), "Gene(12)")

        # Get Gene with protId
        gene12 = self.hpx.get_genes_by_external_id("PANTR2")[0]
        self.assertEqual(str(gene12), "Gene(12)")

        # Get Gene with geneId
        gene12 = self.hpx.get_genes_by_external_id("PANTRg2")[0]
        self.assertEqual(str(gene12), "Gene(12)")


        ###############
        # With filter #
        ###############

        for gene_id in self.filter_genes_ext:
            genes = self.hf.get_genes_by_external_id(gene_id)
            self.assertEqual(len(genes), 1)
            self.assertEqual(genes[0].prot_id, gene_id)

        for gene_id in self.no_filter_genes_ext:
            with self.assertRaises(KeyError):
                self.hf.get_genes_by_external_id(gene_id)

    def test_get_list_extant_genes(self):

        genes = self.h.get_list_extant_genes()
        expected = {'Gene(33)', 'Gene(14)', 'Gene(31)', 'Gene(51)', 'Gene(13)', 'Gene(11)', 'Gene(12)', 'Gene(23)',
                    'Gene(21)', 'Gene(2)', 'Gene(34)', 'Gene(1)', 'Gene(32)', 'Gene(5)', 'Gene(22)', 'Gene(3)',
                    'Gene(41)', 'Gene(53)', 'Gene(43)'}
        self.assertSetEqual(_str_array(genes), expected)

        genes = self.hpx.get_list_extant_genes()
        expected = {'Gene(33)', 'Gene(14)', 'Gene(31)', 'Gene(51)', 'Gene(13)', 'Gene(11)', 'Gene(12)', 'Gene(23)',
                    'Gene(21)', 'Gene(2)', 'Gene(34)', 'Gene(1)', 'Gene(32)', 'Gene(5)', 'Gene(22)', 'Gene(3)',
                    'Gene(41)', 'Gene(53)', 'Gene(43)'}
        self.assertSetEqual(_str_array(genes), expected)

        ###############
        # With filter #
        ###############

        genes = self.hf.get_list_extant_genes()
        self.assertSetEqual(set([g.unique_id for g in genes]),self.filter_genes)

        genes = self.hfstring.get_list_extant_genes()
        self.assertSetEqual(set([g.unique_id for g in genes]), self.filter_genes)

    def test_get_dict_extant_genes(self):

        genes = self.h.get_dict_extant_genes()
        expected = {'34': 'Gene(34)', '33': 'Gene(33)', '2': 'Gene(2)', '12': 'Gene(12)', '5': 'Gene(5)', '11':
            'Gene(11)', '1': 'Gene(1)', '22': 'Gene(22)', '3': 'Gene(3)', '41': 'Gene(41)', '51': 'Gene(51)', '53':
            'Gene(53)', '14': 'Gene(14)', '13': 'Gene(13)', '43': 'Gene(43)', '32': 'Gene(32)', '23': 'Gene(23)', '21':
            'Gene(21)', '31': 'Gene(31)'}
        self.assertDictEqual(_str_dict_one_value(genes), expected)

        genes = self.hpx.get_dict_extant_genes()
        expected = {'34': 'Gene(34)', '33': 'Gene(33)', '2': 'Gene(2)', '12': 'Gene(12)', '5': 'Gene(5)', '11':
            'Gene(11)', '1': 'Gene(1)', '22': 'Gene(22)', '3': 'Gene(3)', '41': 'Gene(41)', '51': 'Gene(51)', '53':
                        'Gene(53)', '14': 'Gene(14)', '13': 'Gene(13)', '43': 'Gene(43)', '32': 'Gene(32)',
                    '23': 'Gene(23)', '21':
                        'Gene(21)', '31': 'Gene(31)'}
        self.assertDictEqual(_str_dict_one_value(genes), expected)

        ###############
        # With filter #
        ###############

        genes = self.hf.get_dict_extant_genes()
        self.assertDictEqual(_str_dict_one_value(genes),{'32': 'Gene(32)', '22': 'Gene(22)', '2': 'Gene(2)', '12': 'Gene(12)'})

    # ExtantGenome

    def test_get_extant_genome_by_name(self):

        with self.assertRaises(KeyError):
            self.h.get_extant_genome_by_name("")

        with self.assertRaises(KeyError):
            self.h.get_extant_genome_by_name("abc")

        with self.assertRaises(KeyError):
            self.h.get_extant_genome_by_name(1)

        h = self.h.get_extant_genome_by_name("HUMAN")
        self.assertEqual(h.name, "HUMAN")

        h = self.hpx.get_extant_genome_by_name("HUMAN")
        self.assertEqual(h.name, "HUMAN")
        ###############
        # With filter #
        ###############

        for gname in self.filter_genome:
            genome = self.hf.get_extant_genome_by_name(gname)
            self.assertEqual(str(genome), gname)

        for gname in self.no_filter_genome:
            genome = self.hf.get_extant_genome_by_name(gname)
            self.assertEqual(len(genome.genes), 0)

    def test_get_list_extant_genomes(self):

        leg = self.h.get_list_extant_genomes()
        expected = {'RATNO', 'HUMAN', 'MOUSE', 'XENTR', 'PANTR', 'CANFA'}
        self.assertSetEqual(_str_array(leg), expected)

        leg = self.hpx.get_list_extant_genomes()
        expected = {'RATNO', 'HUMAN', 'MOUSE', 'XENTR', 'PANTR', 'CANFA'}
        self.assertSetEqual(_str_array(leg), expected)

    # TaxNode

    def test_get_taxon_by_name(self):

        with self.assertRaises(KeyError):
            self.h.get_taxon_by_name("")

        # get Taxnode at internal node
        mammals = self.h.get_taxon_by_name("Mammalia")
        self.assertEqual(mammals.name, "Mammalia")

        rodents = self.hn.get_taxon_by_name("MOUSE/RATNO")
        self.assertEqual(rodents.name, "MOUSE/RATNO")

        # get Taxnode at leaf
        mouse = self.hn.get_taxon_by_name("MOUSE")
        self.assertEqual(mouse.name, "MOUSE")


        # get Taxnode at leaf
        mouse = self.hpx.get_taxon_by_name("MOUSE")
        self.assertEqual(mouse.name, "MOUSE")

    # AncestralGenome

    def test_get_ancestral_genome_by_taxon(self):

        mammals_tax = self.h.get_taxon_by_name("Mammalia")
        mammals = self.h.get_ancestral_genome_by_taxon(mammals_tax)
        self.assertEqual(mammals.name, "Mammalia")

        mammals_tax = self.hpx.get_taxon_by_name("Mammalia")
        mammals = self.hpx.get_ancestral_genome_by_taxon(mammals_tax)
        self.assertEqual(mammals.name, "Mammalia")

        ###############
        # With filter #
        ###############

        with self.assertRaises(KeyError):
            vertebrate_tax = self.hf.get_taxon_by_name("Vertebrata")
            self.h.get_ancestral_genome_by_taxon(vertebrate_tax)

    def test_get_ancestral_genome_by_name(self):

        mammals = self.h.get_ancestral_genome_by_name("Mammalia")
        self.assertEqual(mammals.name, "Mammalia")

        mammals = self.hpx.get_ancestral_genome_by_name("Mammalia")
        self.assertEqual(mammals.name, "Mammalia")

        rodents = self.hn.get_ancestral_genome_by_name("MOUSE/RATNO")
        self.assertEqual(rodents.name, "MOUSE/RATNO")

        ###############
        # With filter #
        ###############

        with self.assertRaises(KeyError):
            self.hf.get_ancestral_genome_by_name("Vertebrata")

    def test_get_list_ancestral_genomes(self):

        lag = self.h.get_list_ancestral_genomes()
        self.assertSetEqual(_str_array(lag), {"Euarchontoglires", "Vertebrata", "Mammalia", "Primates", "Rodents"})

        lag_no_name = self.hn.get_list_ancestral_genomes()
        self.assertSetEqual(_str_array(lag_no_name), {"HUMAN/PANTR/MOUSE/RATNO", "XENTR/HUMAN/PANTR/MOUSE/RATNO/CANFA", "HUMAN/PANTR/MOUSE/RATNO/CANFA", "HUMAN/PANTR", "MOUSE/RATNO"})

    def test_get_ancestral_genome_by_mrca_of_genome_set(self):

        human = self.h.get_extant_genome_by_name(name="HUMAN")
        mouse = self.h.get_extant_genome_by_name(name="MOUSE")

        with self.assertRaises(ValueError):
            self.h.get_ancestral_genome_by_mrca_of_genome_set({mouse})

        with self.assertRaises(TypeError):
            self.h.get_ancestral_genome_by_mrca_of_genome_set({mouse, "human"})

        euarchontoglires = self.h.get_ancestral_genome_by_mrca_of_genome_set({human, mouse})
        self.assertEqual("Euarchontoglires", euarchontoglires.name)

        # If one of the genomes is the mrca
        euarchontoglires2 = self.h.get_ancestral_genome_by_mrca_of_genome_set({euarchontoglires, mouse})
        self.assertEqual("Euarchontoglires", euarchontoglires2.name)

        ###############
        # With filter #
        ###############

        human = self.hf.get_extant_genome_by_name(name="HUMAN")
        mouse = self.hf.get_extant_genome_by_name(name="MOUSE")

        euarchontoglires = self.hf.get_ancestral_genome_by_mrca_of_genome_set({human, mouse})
        self.assertEqual("Euarchontoglires", euarchontoglires.name)

        with self.assertRaises(KeyError):
            frog = self.hf.get_extant_genome_by_name(name="XENTR")
            self.hf.get_ancestral_genome_by_mrca_of_genome_set({human, frog})


class HAMTestPrivate(unittest.TestCase):

    def setUp(self):

        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        nwk_str = utils.get_newick_string(nwk_path, type="nwk")

        nwk_path_no_name = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        nwk_str_no_name = utils.get_newick_string(nwk_path_no_name, type="nwk")

        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        # using newick with name on both internal nodes and leaves
        self.h = ham.Ham(nwk_str, orthoxml_path)

        # using newick with name only at leaves
        self.hn = ham.Ham(nwk_str_no_name, orthoxml_path)

        # using newick with name on both internal nodes and leaves and filter for HOG2

        self.filter_genome = {"HUMAN", "MOUSE", "CANFA", "PANTR"}
        self.filter_genes = {'2', '32', '22', '12'}
        self.filter_genes_ext = {'HUMAN2', 'MOUSE2', 'CANFA2', 'PANTR2'}
        self.filter_hogs = {'2'}
        self.no_filter_genome = {"XENTR", "RATNO"}
        self.no_filter_genes = {'1', '11', '21', '31', '41', '51', '3', '13', '23', '33', '53', '34', '14'}
        self.no_filter_genes_ext = {'HUMAN1', 'PANTR1', 'CANFA1', 'MOUSE1', 'RATNO1', 'XENTR1', 'HUMAN3',
                                    'PANTR3', 'CANFA3', 'MOUSE3', 'XENTR3', 'MOUSE4', 'PANTR4'}
        self.no_filter_hogs = {'1','3'}
        f = ham.ParserFilter()
        f.add_hogs_via_hogId([2])
        self.hf = ham.Ham(nwk_str, orthoxml_path, filter_object=f)

if __name__ == "__main__":
    unittest.main()
