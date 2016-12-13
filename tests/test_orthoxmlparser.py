__author__ = 'admin'

import unittest
from ham import parsers
from xml.etree.ElementTree import XMLParser
from ham import queries
from ham import ham
from ham import utils


class orthoxmlParser(unittest.TestCase):

    def setUp(self):

        taxonomy = ham.build_taxonomy_and_ancestral_genomes(utils.get_newick_string('./tests/simpleEx.nwk', type="nwk"))

        factory = parsers.orthoxmlParser(taxonomy)
        parser = XMLParser(target=factory)

        with open('./tests/simpleEx.orthoxml', 'r') as orthoxml_file:
            for line in orthoxml_file:
                parser.feed(line)

        self.hogs = factory.hogs
        self.genes = factory.genes

    def test_nrOfToplevelFamilies(self):
        self.assertEqual(len(queries.get_top_level_hogs(self.hogs)), 3)

    def test_numberOfGenesPerSpecies(self): ## Be carefull here expected counts are based on how many genes are present in the mappin AND within an HOGs

        expectedCnts = dict(HUMAN=3, PANTR=4, MOUSE=4, RATNO=1,
                            CANFA=3, XENTR=2)

        observedCnts = dict(HUMAN=0, PANTR=0, MOUSE=0, RATNO=0,
                            CANFA=0, XENTR=0)

        for g in self.genes:
            observedCnts[g.species] += 1

        for species in expectedCnts.keys():
            self.assertEqual(observedCnts[species], expectedCnts[species], "number of genes not correct for "+species)


if __name__ == "__main__":
    unittest.main()