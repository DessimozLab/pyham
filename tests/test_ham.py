import unittest
from ham import utils
from ham import ham
import ete3

# This helps to convert elements of list/dictionary to string in order to make easier assertEqual test.
def _str_dict_one_value(dict):
    for kk in dict.keys():
        dict[str(kk)] = dict.pop(kk)
    for k, v in dict.items():
        dict[k] = str(v)
    return dict
def _str_array(array):
    array_converted = []
    for e in array:
        array_converted.append(str(e))
    return set(array_converted)


class HAMTestSetUp(unittest.TestCase):

    def setUp(self):
        self.nwk_str_empty = ""
        self.nwk_str_wrong = "(A,B)x0.4"
        self.nwk_str = utils.get_newick_string('./tests/simpleEx.nwk', type="nwk")
        self.orthoxml_path = './tests/simpleEx.orthoxml'

    def test_wrong_newick_str(self):

        with self.assertRaises(ete3.parser.newick.NewickError):
            ham.HAM(self.nwk_str_empty, self.orthoxml_path, type_hog_file='orthoxml')

        with self.assertRaises(ete3.parser.newick.NewickError):
            ham.HAM(self.nwk_str_wrong, self.orthoxml_path, type_hog_file='orthoxml')

    def test_wrong_type_hog_file(self):

        with self.assertRaises(TypeError):
            ham.HAM(self.nwk_str, self.orthoxml_path, type_hog_file='xs')

        with self.assertRaises(TypeError):
            ham.HAM(self.nwk_str, self.orthoxml_path, type_hog_file='')

        with self.assertRaises(TypeError):
            ham.HAM(self.nwk_str, self.orthoxml_path, type_hog_file=None)


class HAMTest(unittest.TestCase):

    def setUp(self):
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        self.ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')
        self.hogs = self.ham_analysis.get_dict_top_level_hogs()
        self.genes = self.ham_analysis.get_dict_extant_genes()

        self.human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        self.frog = self.ham_analysis._get_extant_genome_by_name(name="XENTR")
        self.mouse = self.ham_analysis._get_extant_genome_by_name(name="MOUSE")
        self.rat = self.ham_analysis._get_extant_genome_by_name(name="RATNO")
        self.chimp = self.ham_analysis._get_extant_genome_by_name(name="PANTR")
        self.vertebrates = self.ham_analysis.get_ancestral_genome_by_mrca_of_genome_set({self.human, self.frog})

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


class HAMTestQuery(unittest.TestCase): # todo tets with filter

    def setUp(self):

        nwk_path = './tests/simpleEx.nwk'
        nwk_str = utils.get_newick_string(nwk_path, type="nwk")

        nwk_path_no_name = './tests/simpleExNoName.nwk'
        nwk_str_no_name = utils.get_newick_string(nwk_path_no_name, type="nwk")

        orthoxml_path = './tests/simpleEx.orthoxml'

        self.h = ham.HAM(nwk_str, orthoxml_path)
        self.hn = ham.HAM(nwk_str_no_name, orthoxml_path)

    # HOG

    def test_get_hog_by_id(self):

        # If Id not exist
        with self.assertRaises(TypeError):
            self.h.get_hog_by_id("d")

        # Get Hog with str(id)
        hog3 = self.h.get_hog_by_id("3")
        self.assertEqual(str(hog3), "<HOG(3)>")

        # Get Hog with int(id)
        hog3 = self.h.get_hog_by_id(3)
        self.assertEqual(str(hog3), "<HOG(3)>")

    def test_get_hog_by_gene(self):

        gene3 = self.h.get_gene_by_id("3")

        # Get Hog with Gene object
        hog3 = self.h.get_hog_by_gene(gene3)
        self.assertEqual(str(hog3), "<HOG(3)>")

        # If gene argument is not a Gene object
        with self.assertRaises(TypeError):
            self.h.get_hog_by_gene("gene")

    def test_get_list_top_level_hogs(self):
        hogs = self.h.get_list_top_level_hogs()
        self.assertSetEqual(_str_array(hogs), {"<HOG(2)>", "<HOG(1)>", "<HOG(3)>"})

    def test_get_dict_top_level_hogs(self):
        hogs = self.h.get_dict_top_level_hogs()
        self.assertDictEqual(_str_dict_one_value(hogs), {"2": "<HOG(2)>", "1": "<HOG(1)>", "3": "<HOG(3)>"})

    # Gene

    def test_get_gene_by_id(self):

        # If Id not exist
        with self.assertRaises(TypeError):
            self.h.get_gene_by_id("d")

        # Get Gene with str(id)
        gene3 = self.h.get_gene_by_id("3")
        self.assertEqual(str(gene3), "Gene(3)")

        # Get Gene with int(id)
        gene3 = self.h.get_gene_by_id(3)
        self.assertEqual(str(gene3), "Gene(3)")

    def test_get_genes_by_external_id(self):

        # If Id not exist
        with self.assertRaises(TypeError):
            self.h.get_genes_by_external_id("d")

        # Get Gene with protId
        gene12 = self.h.get_genes_by_external_id("PANTR2")[0]
        self.assertEqual(str(gene12), "Gene(12)")

        # Get Gene with geneId
        gene12 = self.h.get_genes_by_external_id("PANTRg2")[0]
        self.assertEqual(str(gene12), "Gene(12)")

    def test_get_list_extant_genes(self):

        genes = self.h.get_list_extant_genes()
        expected = {'Gene(33)', 'Gene(14)', 'Gene(31)', 'Gene(51)', 'Gene(13)', 'Gene(11)', 'Gene(12)', 'Gene(23)',
                    'Gene(21)', 'Gene(2)', 'Gene(34)', 'Gene(1)', 'Gene(32)', 'Gene(5)', 'Gene(22)', 'Gene(3)',
                    'Gene(41)', 'Gene(53)', 'Gene(43)'}
        self.assertSetEqual(_str_array(genes), expected)

    def test_get_dict_extant_genes(self):

        genes = self.h.get_dict_extant_genes()
        expected = {'34': 'Gene(34)', '33': 'Gene(33)', '2': 'Gene(2)', '12': 'Gene(12)', '5': 'Gene(5)', '11':
            'Gene(11)', '1': 'Gene(1)', '22': 'Gene(22)', '3': 'Gene(3)', '41': 'Gene(41)', '51': 'Gene(51)', '53':
            'Gene(53)', '14': 'Gene(14)', '13': 'Gene(13)', '43': 'Gene(43)', '32': 'Gene(32)', '23': 'Gene(23)', '21':
            'Gene(21)', '31': 'Gene(31)'}
        self.assertDictEqual(_str_dict_one_value(genes), expected)

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

    def test_get_list_extant_genomes(self):

        leg = self.h.get_list_extant_genomes()
        expected = {'RATNO', 'HUMAN', 'MOUSE', 'XENTR', 'PANTR', 'CANFA'}
        self.assertSetEqual(_str_array(leg), expected)

    # AncestralGenome

    def test_get_ancestral_genome_by_taxon(self):
        pass

    def test_get_list_ancestral_genomes(self):
        pass


if __name__ == "__main__":
    unittest.main()
