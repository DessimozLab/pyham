import unittest
from pyham import taxonomy, EvolutionaryConceptError
import os


class HAMTaxonomy(unittest.TestCase):

    def setUp(self):
        self.newick_str = "((HUMAN, PANTR)Primates,(MOUSE, RATNO)Rodents)Euarchontoglires;"
        self.newick_str_support = "((HUMAN, PANTR)1:0.1,(MOUSE, RATNO)1:0.1)1:0.1;"
        self.newick_str_non_unique = "((HUMAN, HUMAN)Primates,(MOUSE, RATNO)Rodents)Euarchontoglires;"
        self.expected_name_native_ns = {"Primates", "Rodents", "Euarchontoglires"}
        self.expected_name_concat_ns = {"HUMAN/PANTR", "MOUSE/RATNO", "HUMAN/PANTR/MOUSE/RATNO"}

        self.newick_file =os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        self.expected_name_native_nf = {"Primates", "Rodents", "Euarchontoglires","Mammalia","Vertebrata"}
        self.expected_name_concat_nf = {"HUMAN/PANTR", "MOUSE/RATNO", "HUMAN/PANTR/MOUSE/RATNO", "HUMAN/PANTR/MOUSE/RATNO/CANFA", "XENTR/HUMAN/PANTR/MOUSE/RATNO/CANFA"}


        self.phyloxml_file = os.path.join(os.path.dirname(__file__), './data/simpleEx.phyloxml')
        self.phyloxml_file_no_int_name = os.path.join(os.path.dirname(__file__), './data/simpleExNoName.phyloxml')
        self.phyloxml_file_no_clade_name = os.path.join(os.path.dirname(__file__), './data/simpleExNoCladeName.phyloxml')
        self.expected_name_concat_nf2 = {"PANTR/HUMAN", "RATNO/MOUSE", "RATNO/MOUSE/PANTR/HUMAN",
                                        "RATNO/MOUSE/PANTR/HUMAN/CANFA", "XENTR/RATNO/MOUSE/PANTR/HUMAN/CANFA"}

        self.set_species_name = {"HUMAN", "PANTR", "MOUSE", "RATNO", "CANFA", "XENTR"}
        self.set_species_sciname = {"Homo Sapiens", "Chimp", "Mus Musculus", "Ratus Norvegicus", "Canis Familiaris", "Xenopus Tropicallis"}

    def test_non_unique_leaf_names(self):
        with self.assertRaises(KeyError):
            taxonomy.Taxonomy(self.newick_str_non_unique)

        with self.assertRaises(KeyError):
            taxonomy.Taxonomy(self.newick_str_non_unique, use_internal_name=False)


    def test_use_internal_name(self):

        # using the normal newick
        t = taxonomy.Taxonomy(self.newick_str, use_internal_name=True)
        observed_name = {node.name for node in t.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_native_ns, observed_name)

        # using the file newick
        t2 = taxonomy.Taxonomy(self.newick_file, use_internal_name=True, tree_format='newick')
        observed_name = {node.name for node in t2.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_native_nf, observed_name)

        # using the file phyloxml with sciname
        t3 = taxonomy.Taxonomy(self.phyloxml_file, use_internal_name=True, tree_format='phyloxml', phyloxml_internal_name_tag='taxonomy_scientific_name', phyloxml_leaf_name_tag='taxonomy_scientific_name')
        observed_name = {node.name for node in t3.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_native_nf, observed_name)

        # using the file phyloxml with clade name
        t4 = taxonomy.Taxonomy(self.phyloxml_file, use_internal_name=True, tree_format='phyloxml',
                               phyloxml_internal_name_tag='clade_name',
                               phyloxml_leaf_name_tag='taxonomy_scientific_name')
        observed_name = {node.name for node in t4.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_native_nf, observed_name)


    def test_dont_use_internal_name(self):

        # using the normal newick
        t = taxonomy.Taxonomy(self.newick_str, use_internal_name=False)
        observed_name = {node.name for node in t.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_concat_ns, observed_name)

        # using the file newick
        t2 = taxonomy.Taxonomy(self.newick_file, use_internal_name=False, tree_format='newick')
        observed_name = {node.name for node in t2.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_concat_nf, observed_name)

        # using the file phyloxml
        t3 = taxonomy.Taxonomy(self.phyloxml_file, use_internal_name=False, tree_format='phyloxml', phyloxml_leaf_name_tag='taxonomy_code')
        observed_name = {node.name for node in t3.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_concat_nf2, observed_name)


        # using the file phyloxml with phylogeny code
        t6 = taxonomy.Taxonomy(self.phyloxml_file_no_int_name, use_internal_name=False, tree_format='phyloxml', phyloxml_leaf_name_tag='taxonomy_code')
        observed_name = {node.name for node in t6.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_concat_nf2, observed_name)

        # using newick  with support values as internal names
        t_support = taxonomy.Taxonomy(self.newick_str_support, use_internal_name=False)
        observed_name = {node.name for node in t_support.tree.traverse() if node.is_leaf() is False}
        self.assertSetEqual(self.expected_name_concat_ns, observed_name)


    def test_all_correct_name(self):
        # using the file newick
        t2 = taxonomy.Taxonomy(self.newick_file, use_internal_name=True, tree_format='newick')
        observed_name = {node.name for node in t2.tree.traverse() if node.is_leaf() is True}
        self.assertSetEqual(self.set_species_name, observed_name)

        # using the file phyloxml with clade name
        t4 = taxonomy.Taxonomy(self.phyloxml_file, use_internal_name=True, tree_format='phyloxml', phyloxml_internal_name_tag='taxonomy_scientific_name', phyloxml_leaf_name_tag='clade_name')
        observed_name = {node.name for node in t4.tree.traverse() if node.is_leaf() is True}
        self.assertSetEqual(self.set_species_name, observed_name)

        # using the file phyloxml with phylogeny sciname
        t5 = taxonomy.Taxonomy(self.phyloxml_file, use_internal_name=True, tree_format='phyloxml', phyloxml_internal_name_tag='taxonomy_scientific_name' ,
                               phyloxml_leaf_name_tag='taxonomy_scientific_name')
        observed_name = {node.name for node in t5.tree.traverse() if node.is_leaf() is True}
        self.assertSetEqual(self.set_species_sciname, observed_name)

        # using the file phyloxml with phylogeny code
        t6 = taxonomy.Taxonomy(self.phyloxml_file, use_internal_name=True, tree_format='phyloxml', phyloxml_internal_name_tag='taxonomy_scientific_name',
                               phyloxml_leaf_name_tag='taxonomy_code')
        observed_name = {node.name for node in t6.tree.traverse() if node.is_leaf() is True}
        self.assertSetEqual(self.set_species_name, observed_name)


