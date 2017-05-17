import unittest
from pyham import ham
from pyham import utils
from pyham import HOGsMap, MapLateral, MapVertical
import os

##########################################################################################
####  ATTENTION THIS IS QUICK AND DIRTY UNIT TEST TO DEBUG NEED TO BE REDONE     #########
##########################################################################################

def _str_array(array):
    array_converted = []
    for e in array:
        array_converted.append(str(e))
    return array_converted


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


class MapperTestCases:

    class MapperTest(unittest.TestCase):
        def setUp(self):

            nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
            tree_str = utils.get_newick_string(nwk_path, type="nwk")
            orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

            self.ham_analysis = ham.Ham(newick_str=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)

            self.human = self.ham_analysis.get_extant_genome_by_name(name="HUMAN")
            self.frog = self.ham_analysis.get_extant_genome_by_name(name="XENTR")
            self.mouse = self.ham_analysis.get_extant_genome_by_name(name="MOUSE")
            self.rat = self.ham_analysis.get_extant_genome_by_name(name="RATNO")
            self.chimp = self.ham_analysis.get_extant_genome_by_name(name="PANTR")
            self.vertebrates = self.ham_analysis.get_ancestral_genome_by_mrca_of_genome_set({self.human, self.frog})
            self.rodents = self.ham_analysis.get_ancestral_genome_by_mrca_of_genome_set({self.mouse, self.rat})
            self.primates = self.ham_analysis.get_ancestral_genome_by_mrca_of_genome_set({self.human, self.chimp})
            self.euarchontoglires = self.ham_analysis.get_ancestral_genome_by_mrca_of_genome_set(
                {self.human, self.mouse})


class HOGsMapTest(MapperTestCases.MapperTest):

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

    def _get_topLevel_id(self, hog):
        current_hog = hog
        while current_hog.parent is not None:
            current_hog = current_hog.parent
        return current_hog.hog_id

    def test_set_ancestor_and_descendant(self):

        # two genomes on the same lineage
        map = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        self.assertEqual("Euarchontoglires", map.ancestor.taxon.name)
        self.assertEqual(self.human, map.descendant)

    def test_add_genome_not_on_lineage(self):

        with self.assertRaises(TypeError):
            map2 = HOGsMap(self.ham_analysis, self.frog, self.primates)

    def test_UpMap(self):

        def _convert_map(single_mapUp):
            observed_map = {}
            for source, target in single_mapUp.items():
                observed_map[self._get_identifier(source)] = self._get_identifier(target[0])
            return observed_map

        # an extant genome (human) and its ancestor (Vertebrates)
        map = HOGsMap(self.ham_analysis, self.human, self.vertebrates)
        expected_map = {'1': '1', '2': None, '3': '3', "5":None}
        observed_map = _convert_map(map.upMap)
        self.assertDictEqual(expected_map, observed_map)

    def test_buildEventClusters(self):

        def convert_LOSS(LOSS):
            return set([str(hog_ancestor) for hog_ancestor in LOSS])

        def convert_GAIN(GAIN):
            return set([str(new_hog) for new_hog in GAIN])

        def convert_IDENTICAL(IDENTICAL):
            cIDENTICAL = set()
            for hog_ancestor, hog_descendant in IDENTICAL.items():
                x = [str(hog_ancestor), str(hog_descendant)]
                cIDENTICAL.add(frozenset(x))
            return cIDENTICAL

        def convert_DUPLICATE(DUPLICATE):
            cDUPLICATE = set()
            for hog_ancestor, hog_descendants in DUPLICATE.items():
                x = [str(hog_ancestor)]
                for g in hog_descendants:
                    x.append(str(g))
                cDUPLICATE.add(frozenset(x))
            return cDUPLICATE

        # an extant genome (human) and its ancestor (Vertebrates)
        map = HOGsMap(self.ham_analysis, self.human, self.vertebrates)

        expected_LOSS = set()
        self.assertSetEqual(expected_LOSS, convert_LOSS(map.LOSS))

        expected_GAIN = set()
        expected_GAIN.add("Gene(2)")
        expected_GAIN.add("Gene(5)")
        self.assertSetEqual(expected_GAIN, convert_GAIN(map.GAIN))

        expected_IDENTICAL = set()
        expected_IDENTICAL.add(frozenset(["<HOG(1)>", "Gene(1)"]))
        self.assertSetEqual(expected_IDENTICAL, convert_IDENTICAL(map.IDENTICAL))

        expected_DUPLICATE = set()
        expected_DUPLICATE.add(frozenset(["<HOG(3)>", "Gene(3)"]))
        self.assertSetEqual(expected_DUPLICATE, convert_DUPLICATE(map.DUPLICATE))


class VerticalMapperTest(MapperTestCases.MapperTest):
    def test_create_correctly_vertical_map(self):
        vertical_map = MapVertical(self.ham_analysis)

        map = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        vertical_map.add_map(map)

        self.assertEqual(self.euarchontoglires, vertical_map.ancestor)
        self.assertEqual(self.human, vertical_map.descendant)
        self.assertEqual(self.ham_analysis, vertical_map.HAM)

    def test_cannot_add_more_than_one_map(self):
        vertical_map = MapVertical(self.ham_analysis)

        map = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        vertical_map.add_map(map)

        map2 = HOGsMap(self.ham_analysis, self.human, self.vertebrates)
        with self.assertRaises(TypeError):
            vertical_map.add_map(map2)

    def test_can_only_add_HOGMAP_as_map(self):
        vertical_map = MapVertical(self.ham_analysis)

        with self.assertRaises(TypeError):
            vertical_map.add_map("111")

        with self.assertRaises(TypeError):
            vertical_map.add_map("")

    def test_can_only_create_vertical_map_with_ham_object(self):
        with self.assertRaises(TypeError):
            MapVertical()

    def test_get_lost(self):
        vertical_map = MapVertical(self.ham_analysis)
        map = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        vertical_map.add_map(map)

        loss = vertical_map.get_lost()
        self.assertEqual(["<HOG(3.E.2)>"], _str_array(loss))

        vertical_map2 = MapVertical(self.ham_analysis)
        map2 = HOGsMap(self.ham_analysis, self.rat, self.euarchontoglires)
        vertical_map2.add_map(map2)

        loss2 = vertical_map2.get_lost()
        self.assertEqual({"<HOG(2.E)>", "<HOG(3.E.2)>", "<HOG(3.E.1)>"}, set(_str_array(loss2)))

    def test_get_gained(self):
        vertical_map = MapVertical(self.ham_analysis)
        map = HOGsMap(self.ham_analysis, self.human, self.vertebrates)
        vertical_map.add_map(map)

        gain = vertical_map.get_gained()
        self.assertEqual({"Gene(2)", "Gene(5)"}, set(_str_array(gain)))

    def test_get_single(self):
        vertical_map = MapVertical(self.ham_analysis)
        map = HOGsMap(self.ham_analysis, self.human, self.vertebrates)
        vertical_map.add_map(map)

        single = vertical_map.get_identical()
        self.assertDictEqual({'<HOG(1)>': 'Gene(1)'}, _str_dict_one_value(single))

    def test_get_duplicated(self):
        vertical_map = MapVertical(self.ham_analysis)
        map = HOGsMap(self.ham_analysis, self.mouse, self.vertebrates)
        vertical_map.add_map(map)

        duplicate = vertical_map.get_duplicated()
        self.assertDictEqual({'<HOG(3)>': {'Gene(33)', 'Gene(34)'}}, _str_dict_array_value(duplicate))


class LateralMapperTest(MapperTestCases.MapperTest):

    def test_create_correctly_lateral_map(self):
        lateral_map = MapLateral(self.ham_analysis)

        map = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        lateral_map.add_map(map)

        map2 = HOGsMap(self.ham_analysis, self.rodents, self.euarchontoglires)
        lateral_map.add_map(map2)

        self.assertEqual(self.euarchontoglires, lateral_map.ancestor)
        self.assertEqual([self.human, self.rodents], lateral_map.descendants)
        self.assertEqual(self.ham_analysis, lateral_map.HAM)

    def test_cannot_add_map_with_different_ancestor(self):
        lateral_map = MapLateral(self.ham_analysis)

        map = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        lateral_map.add_map(map)

        map2 = HOGsMap(self.ham_analysis, self.human, self.vertebrates)
        with self.assertRaises(TypeError):
            lateral_map.add_map(map2)

    def test_can_only_add_HOGMAP_as_map(self):
        lateral_map = MapLateral(self.ham_analysis)

        with self.assertRaises(TypeError):
            lateral_map.add_map("111")

        with self.assertRaises(TypeError):
            lateral_map.add_map("")

    def test_can_only_create_vertical_map_with_ham_object(self):
        with self.assertRaises(TypeError):
            lateral_map = MapLateral()

    def test_get_lost(self):
        lateral_map = MapLateral(self.ham_analysis)

        map_human_euarc = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        lateral_map.add_map(map_human_euarc)

        map_rat_euarc = HOGsMap(self.ham_analysis, self.rat, self.euarchontoglires)
        lateral_map.add_map(map_rat_euarc)

        loss = lateral_map.get_lost()

        for hog in loss.keys():
            if str(hog) == "<HOG(3.E.2)>":
                H3_E_2 = hog
            elif str(hog) == "<HOG(3.E.1)>":
                H3_E_1 = hog
            elif str(hog) == "<HOG(2.E)>":
                H2_E = hog

        self.assertEqual({self.rat}, set(loss[H2_E]))
        self.assertEqual({self.rat}, set(loss[H3_E_1]))
        self.assertEqual({self.rat, self.human}, set(loss[H3_E_2]))

        # trying if lazy property work
        loss2 = lateral_map.get_lost()
        self.assertEqual({self.rat}, set(loss2[H2_E]))
        self.assertEqual({self.rat}, set(loss2[H3_E_1]))
        self.assertEqual({self.rat, self.human}, set(loss2[H3_E_2]))

    def test_get_gained(self):

        lateral_map = MapLateral(self.ham_analysis)

        map_frog_vertebrate = HOGsMap(self.ham_analysis, self.frog, self.vertebrates)
        lateral_map.add_map(map_frog_vertebrate)

        map_human_euarc = HOGsMap(self.ham_analysis, self.human, self.vertebrates)
        lateral_map.add_map(map_human_euarc)

        gain = lateral_map.get_gained()
        self.assertEqual(set(), set(_str_array(gain[self.frog])))
        self.assertEqual({"Gene(2)", "Gene(5)" }, set(_str_array(gain[self.human])))

    def test_get_single(self):
        lateral_map = MapLateral(self.ham_analysis)

        map_human_euarc = HOGsMap(self.ham_analysis, self.human, self.euarchontoglires)
        lateral_map.add_map(map_human_euarc)

        map_rat_euarc = HOGsMap(self.ham_analysis, self.rat, self.euarchontoglires)
        lateral_map.add_map(map_rat_euarc)

        single = lateral_map.get_identical()

        for hog in single.keys():
            if str(hog) == "<HOG(3.E.1)>":
                H3_E_1 = hog
            elif str(hog) == "<HOG(1.M.E)>":
                H1_M_E = hog
            elif str(hog) == "<HOG(2.E)>":
                H2_E = hog

        self.assertDictEqual({str(self.human): "Gene(3)"},_str_dict_one_value(single[H3_E_1]))
        self.assertDictEqual({str(self.human): "Gene(1)", str(self.rat): "Gene(41)"},_str_dict_one_value(single[H1_M_E]))
        self.assertDictEqual({str(self.human): "Gene(2)"},_str_dict_one_value(single[H2_E]))

    def test_get_duplicated(self):
        lateral_map = MapLateral(self.ham_analysis)

        map_human_vert = HOGsMap(self.ham_analysis, self.human, self.vertebrates)
        lateral_map.add_map(map_human_vert)

        map_mouse_vert = HOGsMap(self.ham_analysis, self.mouse, self.vertebrates)
        lateral_map.add_map(map_mouse_vert)

        duplicate = lateral_map.get_duplicated()

        for hog in duplicate.keys():
            if str(hog) == "<HOG(3)>":
                H = hog

        self.assertDictEqual({str(self.human): {"Gene(3)"}, str(self.mouse):{ 'Gene(33)', 'Gene(34)'}}, _str_dict_array_value(duplicate[H]))

if __name__ == "__main__":
    unittest.main()
