import unittest
from pyham import taxonomy, EvolutionaryConceptError


class HAMTaxonomy(unittest.TestCase):

    def setUp(self):
        self.newick_str = "((HUMAN, PANTR)Primates,(MOUSE, RATNO)Rodents)Euarchontoglires;"
        self.newick_str_support = "((HUMAN, PANTR)1:0.1,(MOUSE, RATNO)1:0.1)1:0.1;"
        self.newick_str_non_unique = "((HUMAN, HUMAN)Primates,(MOUSE, RATNO)Rodents)Euarchontoglires;"

        self.expected_name_native = {"Primates", "Rodents", "Euarchontoglires"}
        self.expected_name_concat = {"HUMAN/PANTR", "MOUSE/RATNO", "HUMAN/PANTR/MOUSE/RATNO"}

    def test_non_unique_leaf_names(self):
        with self.assertRaises(KeyError):
            taxonomy.Taxonomy(self.newick_str_non_unique)

        with self.assertRaises(KeyError):
            taxonomy.Taxonomy(self.newick_str_non_unique, use_internal_name=False)

    def test_use_internal_name(self):

        # using the normal newick
        t = taxonomy.Taxonomy(self.newick_str, use_internal_name=True)
        observed_name = {node.name for node in t.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_native, observed_name)

    def test_dont_use_internal_name(self):

        # using the normal newick
        t = taxonomy.Taxonomy(self.newick_str, use_internal_name=False)
        observed_name = {node.name for node in t.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_concat, observed_name)

        # using newick  with support values as internal names
        t_support = taxonomy.Taxonomy(self.newick_str_support, use_internal_name=False)
        observed_name = {node.name for node in t_support.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_concat, observed_name)