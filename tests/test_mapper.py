import unittest
from ham import ham
from ham import utils
from ham import HOGsMap

class MapperTest(unittest.TestCase):

    def setUp(self):
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        self.ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type='orthoxml')

    def _get_identifier(self, item):
        if isinstance(item, ham.abstractgene.Gene):
            return item.unique_id
        elif isinstance(item, ham.abstractgene.HOG):
            if item.hog_id != None:
                return item.hog_id
            else:
                return item.genome.taxon.name
        elif item == None:
            return item
        else:
            raise TypeError("expect subclass obj of '{}', got {}"
                        .format(ham.abstractgene.AbstractGene.__name__,
                                type(item).__name__))

    def test_set_ancestor_and_descendants(self):

        # genomes not on the same lineage
        human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        mouse = self.ham_analysis._get_extant_genome_by_name(name="MOUSE")

        map = HOGsMap(self.ham_analysis, {human, mouse})
        self.assertEqual("Euarchontoglires", map.ancestor.taxon.name)
        self.assertEqual({human:human, mouse:mouse}, map.descendants)

        # genomes on the same lineage
        human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        mouse = self.ham_analysis._get_extant_genome_by_name(name="MOUSE")
        euarchontoglires = self.ham_analysis._get_mrca_ancestral_genome_from_genome_set({human, mouse})

        map = HOGsMap(self.ham_analysis, {human, euarchontoglires})
        self.assertEqual("Euarchontoglires", map.ancestor.taxon.name)
        self.assertEqual({human:human}, map.descendants)

    def test_UpMaps(self):

        def _convert_map(single_mapUp):
            observed_map = {}
            for source, target in single_mapUp.items():
                observed_map[self._get_identifier(source)] = self._get_identifier(target)
            return observed_map

        # an extant genome and its ancestor
        human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        frog = self.ham_analysis._get_extant_genome_by_name(name="XENTR")
        vertebrates = self.ham_analysis._get_mrca_ancestral_genome_from_genome_set({human, frog})

        map = HOGsMap(self.ham_analysis, {human, vertebrates})
        expected_map = {'1': '1', '2': None, '3': '3'}
        observed_map = _convert_map(map.upMaps[human])
        self.assertDictEqual(expected_map, observed_map)

        # two extant genomes and their MRCA
        human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        mouse = self.ham_analysis._get_extant_genome_by_name(name="MOUSE")
        map = HOGsMap(self.ham_analysis, {human, mouse})

        expected_map_human = {'1': 'Euarchontoglires', '2': 'Euarchontoglires', '3': 'Euarchontoglires'}
        observed_map_human = _convert_map(map.upMaps[human])
        self.assertDictEqual(expected_map_human, observed_map_human)

        expected_map_mouse = {'31': 'Euarchontoglires', '32': 'Euarchontoglires', '33': 'Euarchontoglires', '34': 'Euarchontoglires'}
        observed_map_mouse = _convert_map(map.upMaps[mouse])
        self.assertDictEqual(expected_map_mouse, observed_map_mouse)

        # an extant genomes, an ancestral genome and their MRCA
        human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        mouse = self.ham_analysis._get_extant_genome_by_name(name="MOUSE")
        rat = self.ham_analysis._get_extant_genome_by_name(name="RATNO")
        rodents = self.ham_analysis._get_mrca_ancestral_genome_from_genome_set({mouse, rat})
        map = HOGsMap(self.ham_analysis, {human, rodents})

        expected_map_human = {'1': 'Euarchontoglires', '2': 'Euarchontoglires', '3': 'Euarchontoglires'}
        observed_map_human = _convert_map(map.upMaps[human])
        self.assertDictEqual(expected_map_human, observed_map_human)

        expected_map_rodents = {'Rodents': 'Euarchontoglires', 'Rodents': 'Euarchontoglires', 'Rodents': 'Euarchontoglires', 'Rodents': 'Euarchontoglires'}
        observed_map_rodents = _convert_map(map.upMaps[rodents])
        self.assertDictEqual(expected_map_rodents, observed_map_rodents)


if __name__ == "__main__":
    unittest.main()