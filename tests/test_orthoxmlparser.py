import collections
import unittest
from ham import parsers
from xml.etree.ElementTree import XMLParser
from ham import ham
from ham import utils


class OrthoXMLParserTest(unittest.TestCase):

    def setUp(self):
        taxonomy = ham.build_taxonomy_and_ancestral_genomes(utils.get_newick_string('./tests/simpleEx.nwk', type="nwk"))
        factory = parsers.OrthoXMLParser(taxonomy)
        parser = XMLParser(target=factory)
        with open('./tests/simpleEx.orthoxml', 'r') as orthoxml_file:
            for line in orthoxml_file:
                parser.feed(line)

        self.hogs = factory.toplevel_hogs
        self.genes = factory.extant_gene_map

    def test_numberOfGenesPerSpecies(self):
        expected_cnts = dict(HUMAN=4, PANTR=4, MOUSE=4, RATNO=2,
                            CANFA=3, XENTR=2)
        observed_cnts = collections.defaultdict(int)
        for g in self.genes.values():
            observed_cnts[g.taxon.name] += 1
        self.assertDictEqual(observed_cnts, expected_cnts)


if __name__ == "__main__":
    unittest.main()