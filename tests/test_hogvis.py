import unittest
from pyham import ham, utils, hogvis
import os

class HOGVisTest(unittest.TestCase):
    def setUp(self):

        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/hogvisEx.orthoxml')

        self.ham_analysis = ham.Ham(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)
        hogs = self.ham_analysis.get_dict_top_level_hogs()
        hog = hogs["3"]

        self.hogvis = hogvis.Hogvis(self.ham_analysis, hog )
        self.map = self.hogvis.ogs_mapper

    def test_levelMap(self):
        self.assertEqual(self.map.nr_subhogs_on_level('Euarchontoglires'), 2)
        self.assertEqual(self.map.nr_subhogs_on_level('Primates'), 2)
        self.assertEqual(self.map.nr_subhogs_on_level('Mammalia'), 1)
        self.assertEqual(self.map.nr_subhogs_on_level('LUCA'), 0)

    def test_og_id_equality(self):
        primates_hogs = self.ham_analysis.taxonomy.tree.search_nodes(name="Primates")[0].genome.genes
        primate_hogs_from_map = self.map.levels['Primates']
        glires = self.map.levels['Euarchontoglires']
        self.assertEqual(primates_hogs, primate_hogs_from_map)
        self.assertNotEqual(glires[0], glires[1])

    def test_per_species_mapping(self):
        mapping = self.hogvis.hog.get_all_descendant_genes_clustered_by_species()
        mouse_geneids = [int(z.unique_id) for z in mapping[self.ham_analysis.get_extant_genome_by_name(name="MOUSE")]]
        self.assertEqual([33, 34], mouse_geneids)

    def test_mapping(self):
        expected_HUMAN = {'HUMAN': [[3], [5]], 'Primates': [[3, 5],[]],
                          'Euarchontoglires': [[3, 5], []],
                          'Mammalia': [[3, 5]], 'Vertebrata': [[3, 5]]}
        expected_MOUSE = {'MOUSE': [[33], [34]], 'Euarchontoglires': [[33], [34]],
                          'Mammalia': [[33, 34]], 'Vertebrata': [[33, 34]], 'Rodents': [[33], [34]]}
        per_species = self.hogvis._get_per_species_structure()
        self.assertDictEqual(expected_HUMAN, per_species['HUMAN'])
        self.assertDictEqual(expected_MOUSE, per_species['MOUSE'])

    def test_xrefs(self):
        expected_of_human5 = {'id': "5", 'protId': "HUMAN5", 'geneId': "HUMANg5"}
        xrefs = self.hogvis.xrefs
        self.assertEqual(expected_of_human5, xrefs["5"])

class HOGVisTestNoName(unittest.TestCase):
    def setUp(self):

        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleExNoName.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/hogvisEx.orthoxml')

        self.ham_analysis = ham.Ham(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml')
        hogs = self.ham_analysis.get_dict_top_level_hogs()
        hog = hogs["3"]

        self.hogvis = hogvis.Hogvis(self.ham_analysis, hog )
        self.map = self.hogvis.ogs_mapper

    def test_levelMap(self):
        self.assertEqual(self.map.nr_subhogs_on_level('HUMAN/PANTR/MOUSE/RATNO'), 2)
        self.assertEqual(self.map.nr_subhogs_on_level('HUMAN/PANTR'), 2)
        self.assertEqual(self.map.nr_subhogs_on_level('HUMAN/PANTR/MOUSE/RATNO/CANFA'), 1)
        self.assertEqual(self.map.nr_subhogs_on_level('LUCA'), 0)

    def test_og_id_equality(self):
        HUMAN = self.ham_analysis.get_extant_genome_by_name(name="HUMAN")
        RATNO = self.ham_analysis.get_extant_genome_by_name(name="PANTR")
        primates = self.ham_analysis.get_ancestral_genome_by_mrca_of_genome_set({HUMAN, RATNO})
        primates_hogs = primates.genes
        primate_hogs_from_map = self.map.levels['HUMAN/PANTR']
        glires = self.map.levels['HUMAN/PANTR/MOUSE/RATNO']
        self.assertEqual(primates_hogs, primate_hogs_from_map)
        self.assertNotEqual(glires[0], glires[1])

    def test_per_species_mapping(self):
        mapping = self.hogvis.hog.get_all_descendant_genes_clustered_by_species()
        mouse_geneids = [int(z.unique_id) for z in mapping[self.ham_analysis.get_extant_genome_by_name(name="MOUSE")]]
        self.assertEqual([33, 34], mouse_geneids)

    def test_mapping(self):
        expected_HUMAN = {'HUMAN': [[3], [5]], 'HUMAN/PANTR': [[3, 5],[]],
                          'HUMAN/PANTR/MOUSE/RATNO': [[3, 5], []],
                          'HUMAN/PANTR/MOUSE/RATNO/CANFA': [[3, 5]], 'XENTR/HUMAN/PANTR/MOUSE/RATNO/CANFA': [[3, 5]]}
        expected_MOUSE = {'MOUSE': [[33], [34]], 'HUMAN/PANTR/MOUSE/RATNO': [[33], [34]],
                          'HUMAN/PANTR/MOUSE/RATNO/CANFA': [[33, 34]], 'XENTR/HUMAN/PANTR/MOUSE/RATNO/CANFA': [[33, 34]], 'MOUSE/RATNO': [[33], [34]]}
        per_species = self.hogvis._get_per_species_structure()
        self.assertDictEqual(expected_HUMAN, per_species['HUMAN'])
        self.assertDictEqual(expected_MOUSE, per_species['MOUSE'])

    def test_xrefs(self):
        expected_of_human5 = {'id': "5", 'protId': "HUMAN5", 'geneId': "HUMANg5"}
        xrefs = self.hogvis.xrefs
        self.assertEqual(expected_of_human5, xrefs["5"])

if __name__ == "__main__":
    unittest.main()