import collections
import unittest
from pyham import ham
from pyham import abstractgene
from pyham import utils
import os
import numpy as np
from unittest import skip


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

        self.ham_analysis = ham.Ham(tree_file=nwk_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                    use_internal_name=True)
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


class OrthoXMLParserTest_complexParalogs(unittest.TestCase):
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
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx_complexParalogs.orthoxml')
        orthoxml_path_v3 = os.path.join(os.path.dirname(__file__), './data/SimpleEx_complexParalogGetHOGS3format.xml')

        self.ham_analysis = ham.Ham(tree_file=nwk_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                    use_internal_name=True)
        self.hogs = self.ham_analysis.get_dict_top_level_hogs()
        self.genes = self.ham_analysis.get_dict_extant_genes()

        self.ham_analysis_v3 = ham.Ham(tree_file=nwk_str, hog_file=orthoxml_path_v3, type_hog_file='orthoxml',
                                       use_internal_name=True)
        self.hogs_v3 = self.ham_analysis_v3.get_dict_top_level_hogs()
        self.genes_v3 = self.ham_analysis_v3.get_dict_extant_genes()

    def test_consistency_number(self):

        from pyham import TreeProfile

        # to compute all HOGMAP
        tp = TreeProfile(self.ham_analysis)

        # Verify for each HOGMAP the consistency -- that the parser did his job
        for hg in self.ham_analysis.HOGMaps.values():
            self.assertTrue(hg.consistent)

    def test_numberOfGenesPerSpecies(self):

        def _test_one(genes):
            expected_cnts = dict(HUMAN=6, PANTR=9, MOUSE=8, RATNO=4,
                                 CANFA=7, XENTR=4)
            observed_cnts = collections.defaultdict(int)
            for g in genes.values():
                observed_cnts[g.genome.name] += 1
            self.assertDictEqual(observed_cnts, expected_cnts)

        _test_one(self.genes)
        _test_one(self.genes_v3)

    def test_number_hog_per_ancestral_genome(self):

        def _test_one(ham):
            ags = ham.get_list_ancestral_genomes()
            expected_numbers = {'Vertebrata': 4, 'Mammalia': 6, 'Euarchontoglires': 8, 'Rodents': 5, 'Primates': 7}
            observed_numbers = {'Vertebrata': 0, 'Mammalia': 0, 'Euarchontoglires': 0, 'Rodents': 0, 'Primates': 0}
            for ag in ags:
                observed_numbers[ag.taxon.name] += len(ag.genes)

        _test_one(self.ham_analysis)
        _test_one(self.ham_analysis_v3)

    def test_scores_on_toplevel(self):

        def _test_one(hogs):
            self.assertEqual(hogs["1"].score('consistency'), 1)
            with self.assertRaises(KeyError):
                hogs["2"].score('consistency')
            with self.assertRaises(KeyError):
                hogs["1"].score('coverage')

        _test_one(self.hogs)
        _test_one(self.hogs_v3)

    def test_hog_membership(self):

        def _test_one(hogs):
            hog1 = hogs["1"]
            hog2 = hogs["2"]
            hog3 = hogs["3"]
            hog4 = hogs["4"]
            hog5 = hogs["5"]
            hog6 = hogs["6"]

            expectedMembers_1 = {'51', '21', '1', '11', '31', '41'}
            expectedMembers_2 = {'22', '32', '2', '12'}
            expectedMembers_3 = {'13', '23', '53', '14'}
            expectedMembers_4 = {'15', '16', '24', '25', '54'}
            expectedMembers_5 = {'26', '6', '17', '35', '44', '7', '18', '36', '37', '45'}
            expectedMembers_6 = {'19', '38', '27', '56'}

            hog1_genes = hog1.get_all_descendant_genes()
            hog2_genes = hog2.get_all_descendant_genes()
            hog3_genes = hog3.get_all_descendant_genes()
            hog4_genes = hog4.get_all_descendant_genes()
            hog5_genes = hog5.get_all_descendant_genes()
            hog6_genes = hog6.get_all_descendant_genes()

            members_1 = set(x.unique_id for x in hog1_genes)
            members_2 = set(x.unique_id for x in hog2_genes)
            members_3 = set(x.unique_id for x in hog3_genes)
            members_4 = set(x.unique_id for x in hog4_genes)
            members_5 = set(x.unique_id for x in hog5_genes)
            members_6 = set(x.unique_id for x in hog6_genes)

            self.assertSetEqual(expectedMembers_1, members_1)
            self.assertSetEqual(expectedMembers_2, members_2)
            self.assertSetEqual(expectedMembers_3, members_3)
            self.assertSetEqual(expectedMembers_4, members_4)
            self.assertSetEqual(expectedMembers_5, members_5)
            self.assertSetEqual(expectedMembers_6, members_6)

        _test_one(self.hogs)
        _test_one(self.hogs_v3)

    def test_simple_hog_structure(self):

        def _test_one(hog):
            # Vertebrata
            self.assertEqual("Vertebrata", hog.genome.taxon.name)
            self._check_children_consistency(hog, ["51", "Mammalia"])

            # Mammalia
            mammalia = self._get_child_by_identifier(hog, "Mammalia")
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

        _test_one(self.hogs["1"])
        _test_one(self.hogs_v3["1"])

    def test_incomplete_lineage_taxon_expansion(self):

        def _test_one(hog2):
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

        _test_one(self.hogs["2"])
        _test_one(self.hogs_v3["2"])

    def test_hog_with_duplication(self):

        def _test_one(hog3):
            # Vertebrata
            self.assertEqual("Vertebrata", hog3.genome.taxon.name)
            self._check_children_consistency(hog3, ["53", "Mammalia"])
            self.assertFalse(hog3.arose_by_duplication)

            # Mammalia
            mammalia = self._get_child_by_identifier(hog3, "Mammalia")
            self._check_children_consistency(mammalia, ["23", "Euarchontoglires"])
            self.assertFalse(mammalia.arose_by_duplication)

            # Euarchontoglires
            euarchontoglire = self._get_child_by_identifier(mammalia, "Euarchontoglires")
            self.assertFalse(euarchontoglire.arose_by_duplication)
            self._check_children_consistency(euarchontoglire, ["Primates"])

            # Primates
            primates = self._get_child_by_identifier(euarchontoglire, "Primates")
            self.assertFalse(primates.arose_by_duplication)
            self._check_children_consistency(primates, ["13", "14"])
            g13 = self.ham_analysis.get_gene_by_id("13")
            self.assertTrue(g13.arose_by_duplication)
            g14 = self.ham_analysis.get_gene_by_id("14")
            self.assertTrue(g14.arose_by_duplication)

        _test_one(self.hogs["3"])
        _test_one(self.hogs_v3["3"])

    def test_hog_with_sister_duplication(self):

        def _test_one(hog4):
            # Vertebrata
            self.assertEqual("Vertebrata", hog4.genome.taxon.name)
            self._check_children_consistency(hog4, ["54", "Mammalia"])
            self.assertFalse(hog4.arose_by_duplication)

            # Mammalia
            mammalia = self._get_child_by_identifier(hog4, "Mammalia")
            self._check_children_consistency(mammalia, ["24", "25", "Euarchontoglires"])
            self.assertFalse(mammalia.arose_by_duplication)
            g24 = self.ham_analysis.get_gene_by_id("24")
            self.assertTrue(g24.arose_by_duplication)
            g25 = self.ham_analysis.get_gene_by_id("25")
            self.assertTrue(g25.arose_by_duplication)

            # Euarchontoglires
            euarchontoglire = self._get_child_by_identifier(mammalia, "Euarchontoglires")
            self.assertFalse(euarchontoglire.arose_by_duplication)
            self._check_children_consistency(euarchontoglire, ["Primates"])

            # Primates
            primates = self._get_child_by_identifier(euarchontoglire, "Primates")
            self.assertFalse(primates.arose_by_duplication)
            self._check_children_consistency(primates, ["15", "16"])
            g15 = self.ham_analysis.get_gene_by_id("15")
            self.assertTrue(g15.arose_by_duplication)
            g16 = self.ham_analysis.get_gene_by_id("16")
            self.assertTrue(g16.arose_by_duplication)

        _test_one(self.hogs["4"])
        _test_one(self.hogs_v3["4"])

    def test_hog_with_nested_duplication(self):

        def _test_one(hog5):
            # Mammalia
            self.assertEqual("Mammalia", hog5.genome.taxon.name)
            self._check_children_consistency(hog5, ["26", "Euarchontoglires", "Euarchontoglires"])
            self.assertFalse(hog5.arose_by_duplication)

            # Euarchontoglires
            euarchontoglires = self._get_child_by_identifier(hog5, "Euarchontoglires")

            if self.ham_analysis.get_gene_by_id("35") in euarchontoglires[0].get_all_descendant_genes():
                euarchontoglires_A = euarchontoglires[1]
                euarchontoglires_B = euarchontoglires[0]
            else:
                euarchontoglires_A = euarchontoglires[0]
                euarchontoglires_B = euarchontoglires[1]

            # Euarchontoglires paralogs A
            self.assertNotEqual(False, euarchontoglires_A.arose_by_duplication)
            self._check_children_consistency(euarchontoglires_A, ["Primates", "Rodents"])

            primates = self._get_child_by_identifier(euarchontoglires_A, "Primates")
            self.assertFalse(primates.arose_by_duplication)
            self._check_children_consistency(primates, ["6", "17"])

            rodents = self._get_child_by_identifier(euarchontoglires_A, "Rodents")
            self.assertFalse(rodents.arose_by_duplication)
            self._check_children_consistency(rodents, ["36", "37", "44"])
            g36 = self.ham_analysis.get_gene_by_id("36")
            self.assertNotEqual(False, g36.arose_by_duplication)
            g37 = self.ham_analysis.get_gene_by_id("37")
            self.assertNotEqual(False, g37.arose_by_duplication)

            # Euarchontoglires paralogs B
            self.assertNotEqual(False, euarchontoglires_B.arose_by_duplication)
            self._check_children_consistency(euarchontoglires_B, ["Primates", "Rodents"])

            primates = self._get_child_by_identifier(euarchontoglires_B, "Primates")
            self.assertFalse(primates.arose_by_duplication)
            self._check_children_consistency(primates, ["7", "18"])

            rodents = self._get_child_by_identifier(euarchontoglires_B, "Rodents")
            self.assertFalse(rodents.arose_by_duplication)
            self._check_children_consistency(rodents, ["35", "45"])

        _test_one(self.hogs["5"])
        _test_one(self.hogs_v3["5"])

    def test_hog_with_losses_duplication(self):

        def _test_one(hog6):
            # Vertebrata
            self.assertEqual("Vertebrata", hog6.genome.taxon.name)
            self._check_children_consistency(hog6, ["56", "Mammalia"])
            self.assertFalse(hog6.arose_by_duplication)

            # Mammalia
            mammalia = self._get_child_by_identifier(hog6, "Mammalia")
            self._check_children_consistency(mammalia, ["27", "Euarchontoglires", "Euarchontoglires"])
            self.assertFalse(mammalia.arose_by_duplication)

            # Euarchontoglires
            euarchontoglires = self._get_child_by_identifier(mammalia, "Euarchontoglires")

            if self.ham_analysis.get_gene_by_id("38") in euarchontoglires[0].get_all_descendant_genes():
                euarchontoglires_A = euarchontoglires[1]
                euarchontoglires_B = euarchontoglires[0]
            else:
                euarchontoglires_A = euarchontoglires[0]
                euarchontoglires_B = euarchontoglires[1]


            # Euarchontoglires paralogs A
            self.assertNotEqual(False, euarchontoglires_A.arose_by_duplication)
            self._check_children_consistency(euarchontoglires_A, ["Primates"])

            primates = self._get_child_by_identifier(euarchontoglires_A, "Primates")
            self.assertFalse(primates.arose_by_duplication)
            self._check_children_consistency(primates, ["19"])
            g19 = self.ham_analysis.get_gene_by_id("19")
            self.assertEqual(False, g19.arose_by_duplication)

            # Euarchontoglires paralogs B
            self.assertNotEqual(False, euarchontoglires_B.arose_by_duplication)
            self._check_children_consistency(euarchontoglires_B, ["Rodents"])

            rodents = self._get_child_by_identifier(euarchontoglires_B, "Rodents")
            self.assertFalse(rodents.arose_by_duplication)
            self._check_children_consistency(rodents, ["38"])
            g38 = self.ham_analysis.get_gene_by_id("38")
            self.assertEqual(False, g38.arose_by_duplication)

        _test_one(self.hogs["6"])
        _test_one(self.hogs_v3["6"])

    def test_duplication_children_all_in_hog_children(self):

        for top_hog in self.ham_analysis.get_list_top_level_hogs():

            for sub_hog in top_hog.get_all_descendant_hogs():

                for dup in sub_hog.duplications:

                    for e in dup.children:
                        self.assertTrue(e in sub_hog.children)


class OrthoXMLParser_Weird_Element_In_duplicationChildren(unittest.TestCase):

    def setUp(self):

        nwk_path = os.path.join(os.path.dirname(__file__), './data/parser/tree.newick')
        nwk_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/parser/conflict_duplication_children.orthoxml')

        self.ham_analysis = ham.Ham(tree_file=nwk_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                    use_internal_name=True)

    def test_duplication_children_all_in_hog_children(self):

        for top_hog in self.ham_analysis.get_list_top_level_hogs():

            for sub_hog in top_hog.get_all_descendant_hogs():

                for dup in sub_hog.duplications:

                    for e in dup.children:
                        self.assertTrue(e in sub_hog.children)


"""
Tests for FilterOrthoXMLParser are made indirectly inside Ham unit tests.
"""

if __name__ == "__main__":
    unittest.main()
