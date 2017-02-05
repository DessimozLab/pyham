import unittest
from unittest import skip
from ham import ham
from ham import utils
from ham import HOGsMap



##########################################################################################
####  ATTENTION THIS I QUICK AND DIRTY UNIT TEST TO DEBUG NEED TO BE REDO        #########
##########################################################################################


class MapperTest(unittest.TestCase):

    def setUp(self):
        nwk_path = './tests/simpleEx.nwk'
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = './tests/simpleEx.orthoxml'
        self.ham_analysis = ham.HAM(newick_str=tree_str, hog_file=orthoxml_path, type='orthoxml')

        self.human = self.ham_analysis._get_extant_genome_by_name(name="HUMAN")
        self.frog = self.ham_analysis._get_extant_genome_by_name(name="XENTR")
        self.mouse = self.ham_analysis._get_extant_genome_by_name(name="MOUSE")
        self.rat = self.ham_analysis._get_extant_genome_by_name(name="RATNO")
        self.chimp = self.ham_analysis._get_extant_genome_by_name(name="PANTR")
        self.vertebrates = self.ham_analysis._get_mrca_ancestral_genome_from_genome_set({self.human, self.frog})
        self.rodents = self.ham_analysis._get_mrca_ancestral_genome_from_genome_set({self.mouse, self.rat})
        self.primates = self.ham_analysis._get_mrca_ancestral_genome_from_genome_set({self.human, self.chimp})
        self.euarchontoglires = self.ham_analysis._get_mrca_ancestral_genome_from_genome_set({self.human, self.mouse})

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

    def test_set_ancestor_and_descendants(self):

        # genomes (two extant genomes) not on the same lineage
        map = HOGsMap(self.ham_analysis, {self.human, self.mouse})
        self.assertEqual("Euarchontoglires", map.ancestor.taxon.name)
        self.assertEqual({self.human:self.human, self.mouse:self.mouse}, map.descendants)

        # two genomes on the same lineage
        map = HOGsMap(self.ham_analysis, {self.human, self.euarchontoglires})
        self.assertEqual("Euarchontoglires", map.ancestor.taxon.name)
        self.assertEqual({self.human:self.human}, map.descendants)

        # genomes (one extant genomes and one ancestral genome) not on the same lineage
        map = HOGsMap(self.ham_analysis, {self.human, self.rodents})
        self.assertEqual("Euarchontoglires", map.ancestor.taxon.name)
        self.assertEqual({self.human:self.human, self.rodents:self.rodents}, map.descendants)

        # genomes (two ancestral genome) not on the same lineage
        map = HOGsMap(self.ham_analysis, {self.rodents, self.primates})
        self.assertEqual("Euarchontoglires", map.ancestor.taxon.name)
        self.assertEqual({self.primates:self.primates, self.rodents:self.rodents}, map.descendants)

    def test_UpMaps(self):

        def _convert_map(single_mapUp):
            observed_map = {}
            for source, target in single_mapUp.items():
                observed_map[self._get_identifier(source)] = self._get_identifier(target[0])
            return observed_map

        # an extant genome (human) and its ancestor (Vertebrates)
        map = HOGsMap(self.ham_analysis, {self.human, self.vertebrates})
        expected_map = {'1': '1', '2': None, '3': '3'}
        observed_map = _convert_map(map.upMaps[self.human])
        self.assertDictEqual(expected_map, observed_map)

         # two extant genomes(human,mouse) and their MRCA(Euarchontoglires)
        map = HOGsMap(self.ham_analysis, {self.human, self.mouse})

        expected_map_human = {'1': 'Euarchontoglires', '2': 'Euarchontoglires', '3': 'Euarchontoglires'}
        observed_map_human = _convert_map(map.upMaps[self.human])
        self.assertDictEqual(expected_map_human, observed_map_human)

        expected_map_mouse = {'31': 'Euarchontoglires', '32': 'Euarchontoglires', '33': 'Euarchontoglires', '34': 'Euarchontoglires'}
        observed_map_mouse = _convert_map(map.upMaps[self.mouse])
        self.assertDictEqual(expected_map_mouse, observed_map_mouse)

        # an extant genomes, an ancestral genome and their MRCA
        map = HOGsMap(self.ham_analysis, {self.human, self.rodents})

        expected_map_human = {'1': 'Euarchontoglires', '2': 'Euarchontoglires', '3': 'Euarchontoglires'}
        observed_map_human = _convert_map(map.upMaps[self.human])
        self.assertDictEqual(expected_map_human, observed_map_human)

        expected_map_rodents = {72: '1', 32: '2', 34: "3", 33: '3'}
        observed_map_rodents = {}
        for hog_rodents, hog_euarch in map.upMaps[self.rodents].items(): #dirty but it's work (the trick if to sum up the children id..)
            sum_child = 0
            for child in hog_rodents.children:
                sum_child += int(child.unique_id)
            observed_map_rodents[sum_child]= self._get_topLevel_id(hog_euarch[0])
        self.assertDictEqual(expected_map_rodents, observed_map_rodents)

    @skip
    def test_DownMap(self):

        def _convert_map(downMap):
            observed_map = {}
            for hog_top, value in downMap.items():
                item = {}
                for genome, list_hogs in value.items():
                    list_cleaned = []
                    for e in list_hogs:
                        list_cleaned.append(self._get_identifier(e))
                    item[genome]=list_cleaned
                observed_map[self._get_identifier(hog_top)] = item
            return observed_map

        # an extant genome and its ancestor
        map = HOGsMap(self.ham_analysis, {self.human, self.vertebrates})
        expected_map = {'1': {self.human:['1']}, '3': {self.human:['3']}}
        observed_map = _convert_map(map.downMap)
        self.assertDictEqual(expected_map, observed_map)


        # two extant genomes and their MRCA
        map = HOGsMap(self.ham_analysis, {self.human, self.mouse})
        expected_map = {'<HOG(1)>': {self.human:['1'], self.mouse:['31']}, '<HOG(2)>': {self.human:['2'], self.mouse:['32']}, '<HOG(3)>': {self.human:['3'], self.mouse:['33']}, '<HOG()>': {self.human:[], self.mouse:['34']}}
        observed_map = {}
        for hog_top, value in map.downMap.items(): #dirty but it's work
                item = {}
                for genome, list_hogs in value.items():
                    list_cleaned = []
                    for e in list_hogs:
                        list_cleaned.append(self._get_identifier(e))
                    item[genome]=list_cleaned
                if len(item[self.human]) > 0:
                    observed_map["<HOG({})>".format(item[self.human][0])] = item
                else:
                    observed_map["<HOG()>"] = item
        self.assertDictEqual(expected_map, observed_map)

        # an extant genomes, an ancestral genome and their MRCA
        map = HOGsMap(self.ham_analysis, {self.human, self.rodents})
        expected_map = {'<HOG(1)>': {self.human:['1'], self.rodents:72}, '<HOG(2)>': {self.human:['2'], self.rodents:32}, '<HOG(3)>': {self.human:['3'], self.rodents:33}, '<HOG()>': {self.human:[], self.rodents:34}}
        observed_map = {}
        for hog_top, value in map.downMap.items(): #dirty but it's work
                item = {}
                for genome, list_hogs in value.items():

                    if genome == self.human:
                        list_cleaned = []
                        for e in list_hogs:
                            list_cleaned.append(self._get_identifier(e))
                        item[genome]=list_cleaned
                    else:
                        child_sum = 0
                        for h in list_hogs:
                            for child in h.children:
                                child_sum += int(child.unique_id)
                        item[genome]=child_sum
                if len(item[self.human]) > 0:
                    observed_map["<HOG({})>".format(item[self.human][0])] = item
                else:
                    observed_map["<HOG()>"] = item
        self.assertDictEqual(expected_map, observed_map)

    def test_buildEventClusters(self):

        def convert_LOSS(LOSS):
            cLOSS = {}
            for hog_ancestor, descendants in LOSS.items():
                cLOSS[str(hog_ancestor)] = descendants
            return cLOSS

        def convert_GAIN(GAIN):
            cGAIN = {}
            for genome, list_genes in GAIN.items():
                clist = []
                for g in list_genes:
                    clist.append(str(g))
                cGAIN[genome]=clist
            return cGAIN

        def convert_SINGLE(SINGLE):
            cSINGLE = set()
            for hog_ancestor, descendant in SINGLE.items():
                x = []
                for genome, gene in descendant.items():
                    x = x + [str(hog_ancestor), genome, str(gene)]

                cSINGLE.add(frozenset(x))
            return cSINGLE

        def convert_DUPLICATE(DUPLICATE):
            cDUPLICATE = set()
            for hog_ancestor, descendant in DUPLICATE.items():
                x = []
                for genome, genes in descendant.items():
                    x = x + [str(hog_ancestor), genome]
                    for g in genes:
                        x.append(str(g))
                cDUPLICATE.add(frozenset(x))
            return cDUPLICATE

        # an extant genome (human) and its ancestor (Vertebrates)
        map = HOGsMap(self.ham_analysis, {self.human, self.vertebrates})

        expected_LOSS = {}
        self.assertDictEqual(expected_LOSS, convert_LOSS(map.LOSS))

        expected_GAIN = {self.human: ["Gene(2)"]}
        self.assertDictEqual(expected_GAIN, convert_GAIN(map.GAIN))

        expected_SINGLE = set()
        expected_SINGLE.add(frozenset(["<HOG(1)>",self.human,"Gene(1)"]))
        self.assertSetEqual(expected_SINGLE, convert_SINGLE(map.SINGLE))

        expected_DUPLICATE = set()
        expected_DUPLICATE.add(frozenset(["<HOG(3)>",self.human,"Gene(3)"]))
        self.assertSetEqual(expected_DUPLICATE, convert_DUPLICATE(map.DUPLICATE))

        # two extant genomes(human,mouse) and their MRCA(Euarchontoglires)
        map = HOGsMap(self.ham_analysis, {self.human, self.mouse})

        expected_LOSS = {"<HOG()>": [self.human]}
        self.assertDictEqual(expected_LOSS, convert_LOSS(map.LOSS))

        expected_GAIN = {self.human:[], self.mouse:[]}
        self.assertDictEqual(expected_GAIN, convert_GAIN(map.GAIN))

        expected_SINGLE = set()
        expected_SINGLE.add(frozenset(["<HOG()>", self.human,"Gene(1)",self.mouse,"Gene(31)"]))
        expected_SINGLE.add(frozenset(["<HOG()>",self.human,"Gene(2)", self.mouse,"Gene(32)"]))
        expected_SINGLE.add(frozenset(["<HOG()>",self.human,"Gene(3)", self.mouse,"Gene(33)"]))
        expected_SINGLE.add(frozenset(["<HOG()>", self.mouse,"Gene(34)"]))
        self.assertSetEqual(expected_SINGLE, convert_SINGLE(map.SINGLE))

        expected_DUPLICATE = set()
        self.assertSetEqual(expected_DUPLICATE, convert_DUPLICATE(map.DUPLICATE))


if __name__ == "__main__":
    unittest.main()