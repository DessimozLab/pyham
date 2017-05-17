import collections
import unittest
from pyham import ham
from pyham import utils
import os

class OrthoXMLParserTest(unittest.TestCase):
    def _get_identifier(self, item):
        if isinstance(item, ham.abstractgene.Gene):
            return item.unique_id
        elif isinstance(item, ham.abstractgene.HOG):
            return item.genome.taxon.name
        else:
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(ham.abstractgene.AbstractGene.__name__,
                                    type(item).__name__))

    def _get_child_by_identifier(self, item, query):
        founded_children = []
        for child in item.children:
            if self._get_identifier(child) == query:
                founded_children.append(child)
        if len(founded_children) == 1:
            return founded_children[0]
        return founded_children

    def _check_children_consistency(self, hog, expected_members):
        expected_members = sorted(list(expected_members))
        observed_members = sorted(list(self._get_identifier(item) for item in hog.children))
        self.assertEqual(expected_members, observed_members)

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        nwk_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        self.ham_analysis = ham.Ham(newick_str=nwk_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)
        self.hogs = self.ham_analysis.get_dict_top_level_hogs()
        self.genes = self.ham_analysis.get_dict_extant_genes()

    def test_numberOfGenesPerSpecies(self):
        expected_cnts = dict(HUMAN=4, PANTR=4, MOUSE=4, RATNO=2,
                             CANFA=3, XENTR=2)
        observed_cnts = collections.defaultdict(int)
        for g in self.genes.values():
            observed_cnts[g.genome.name] += 1
        self.assertDictEqual(observed_cnts, expected_cnts)

    def test_number_hog_per_ancestral_genome(self):
        ags = self.ham_analysis.get_list_ancestral_genomes()
        expected_numbers = {'Vertebrata': 2, 'Mammalia': 3, 'Euarchontoglires': 4, 'Rodents': 4, 'Primates': 4}
        observed_numbers = {'Vertebrata': 0, 'Mammalia': 0, 'Euarchontoglires': 0, 'Rodents': 0, 'Primates': 0}
        for ag in ags:
            observed_numbers[ag.taxon.name] += len(ag.genes)
        self.assertDictEqual(expected_numbers, observed_numbers)

    def test_scores_on_toplevel(self):
        self.assertEqual(self.hogs["1"].score('consistency'), 1)
        with self.assertRaises(KeyError):
            self.hogs["2"].score('consistency')
        with self.assertRaises(KeyError):
            self.hogs["1"].score('coverage')

    def test_hog_membership(self):

        hog1 = self.hogs["1"]
        hog2 = self.hogs["2"]
        hog3 = self.hogs["3"]

        expectedMembers_1 = {'51', '21', '1', '11', '31', '41'}
        expectedMembers_2 = {'22', '32', '2', '12'}
        expectedMembers_3 = {'3', '13', '23', '33', '53', '14', '34'}

        hog1_genes = hog1.get_all_descendant_genes()
        hog2_genes = hog2.get_all_descendant_genes()
        hog3_genes = hog3.get_all_descendant_genes()

        members_1 = set(x.unique_id for x in hog1_genes)
        members_2 = set(x.unique_id for x in hog2_genes)
        members_3 = set(x.unique_id for x in hog3_genes)

        self.assertSetEqual(expectedMembers_1, members_1)
        self.assertSetEqual(expectedMembers_2, members_2)
        self.assertSetEqual(expectedMembers_3, members_3)

    def test_simple_hog_structure(self):

        hog1 = self.hogs["1"]

        # Vertebrata
        self.assertEqual("Vertebrata", hog1.genome.taxon.name)
        self._check_children_consistency(hog1, ["51", "Mammalia"])

        # Mammalia
        mammalia = self._get_child_by_identifier(hog1, "Mammalia")
        self._check_children_consistency(mammalia, ["21", "Euarchontoglires"])

        # Euarchontoglires
        euarchontoglires = self._get_child_by_identifier(mammalia, "Euarchontoglires")
        self._check_children_consistency(euarchontoglires, ["Rodents", "Primates"])

        # Rodents
        rodents = self._get_child_by_identifier(euarchontoglires, "Rodents")
        self._check_children_consistency(rodents, ["41", "31"])

        # Primates
        primates = self._get_child_by_identifier(euarchontoglires, "Primates")
        self._check_children_consistency(primates, ["1", "11"])

    def test_incomplete_lineage_taxon_expansion(self):

        hog2 = self.hogs["2"]

        # Mammalia
        self.assertEqual("Mammalia", hog2.genome.taxon.name)
        self._check_children_consistency(hog2, ["22", "Euarchontoglires"])

        # Euarchontoglires
        euarchontoglires = self._get_child_by_identifier(hog2, "Euarchontoglires")
        self._check_children_consistency(euarchontoglires, ["Rodents", "Primates"])

        # Primates
        primates = self._get_child_by_identifier(euarchontoglires, "Primates")
        self._check_children_consistency(primates, ["2", "12"])

        # Rodents
        rodents = self._get_child_by_identifier(euarchontoglires, "Rodents")
        self._check_children_consistency(rodents, {"32"})

    def test_hog_with_duplication(self):

        hog3 = self.hogs["3"]

        # Vertebrata
        self.assertEqual("Vertebrata", hog3.genome.taxon.name)
        self._check_children_consistency(hog3, ["53", "Mammalia"])
        self.assertFalse(hog3.arose_by_duplication)


        # Mammalia
        mammalia = self._get_child_by_identifier(hog3, "Mammalia")
        self._check_children_consistency(mammalia, ["23", "Euarchontoglires", "Euarchontoglires"])
        self.assertFalse(mammalia.arose_by_duplication)

        # Euarchontoglires
        euarchontoglires = self._get_child_by_identifier(mammalia, "Euarchontoglires")
        for euarchontoglire in euarchontoglires:
            self.assertTrue(euarchontoglire.arose_by_duplication)
            self._check_children_consistency(euarchontoglire, ["Primates", "Rodents"])

            # Primates and Rodents
            primates = self._get_child_by_identifier(euarchontoglire, "Primates")
            self.assertFalse(primates.arose_by_duplication)
            rodents = self._get_child_by_identifier(euarchontoglire, "Rodents")
            self.assertFalse(rodents.arose_by_duplication)

            if len(primates.children) == 2:
                self._check_children_consistency(rodents, ["33"])
                self._check_children_consistency(primates, ["3", "13"])
            else:
                self._check_children_consistency(rodents, ["34"])
                self._check_children_consistency(primates, ["14"])

"""
Tests for FilterOrthoXMLParser are made indirectly inside Ham unit tests.
"""

if __name__ == "__main__":
    unittest.main()
