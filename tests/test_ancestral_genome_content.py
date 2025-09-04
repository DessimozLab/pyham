import unittest
import pytest
from pyham import utils
from pyham import ham
from pathlib import Path
import os

@pytest.fixture(params=[
    ("tom.orthoxml", None),
    ("tom.orthoxml", "tomato.nwk"),
])
def tomato_ham(request):
    oxml, nwk = request.param
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data', oxml)
    tree_path = os.path.join(os.path.dirname(__file__), './data', nwk) if nwk else None
    return ham.Ham(tree_file=tree_path, hog_file=orthoxml_path, tree_format="newick", use_internal_name=True)


def test_ancestral_genome(tomato_ham):
    anc_genome = tomato_ham.get_ancestral_genome_by_name("NODE_29")
    assert anc_genome is not None, "Ancestral genome should not be None"
    assert len(anc_genome.genes) == 262, "Ancestral genome should have genes"


def test_fetch_subhog_by_id(tomato_ham):
    query_hog_id = "HOG:0027790.3c.9ba.5i.8ai_69"
    subhog = tomato_ham.get_hog_by_id(query_hog_id)
    assert subhog.hog_id == query_hog_id, "SubHOG should match the queried ID"


def test_fetch_roothog_by_id(tomato_ham):
    query_hog_id = "HOG:0027790_45"
    root_hog = tomato_ham.get_hog_by_id(query_hog_id)
    assert root_hog.hog_id == query_hog_id, "Root HOG should match the queried ID"
    assert root_hog.parent is None, "Root HOG should not have a parent"


def test_fetch_nonexistent_subhog(tomato_ham):
    query_hog_id = "HOG:0027790.3c.9b.88a_69"
    with pytest.raises(KeyError):
        tomato_ham.get_hog_by_id(query_hog_id)


def test_tomato_hog_follow_species_tree(tomato_ham):
    for gene in tomato_ham.get_list_extant_genes():
        hog = gene
        while hog.parent:
            assert hog.genome.taxon.parent == hog.parent.genome.taxon, "HOG should follow species tree"
            hog = hog.parent


####
# tests on hog_1074943 (that failed while building edgehog)
@pytest.fixture(params=[
    ("hog_1074943.augmented.orthoxml", None),
    ("hog_1074943.augmented.orthoxml", "hog_1074943.nwk")
])
def hog_1074943_ham(request):
    oxml, nwk = request.param
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data', oxml)
    tree_path = os.path.join(os.path.dirname(__file__), './data', nwk) if nwk else None
    return ham.Ham(tree_file=tree_path, hog_file=orthoxml_path, tree_format="newick", use_internal_name=True)


def test_root_level(hog_1074943_ham):
    assert len(hog_1074943_ham.get_list_top_level_hogs()) == 1, "HOG 1074943 should have one top level hog"
    hog = hog_1074943_ham.get_list_top_level_hogs()[0]
    assert hog.hog_id == "HOG:1074943_0", "Top level HOG ID should match"
    assert hog.genome.name == "LUCA", "Top level HOG level should be LUCA"


def test_hog_follow_species_tree(hog_1074943_ham):
    for gene in hog_1074943_ham.get_list_extant_genes():
        hog = gene
        while hog.parent:
            assert hog.genome.taxon.parent == hog.parent.genome.taxon, "HOG should follow species tree"
            hog = hog.parent


class PyHAMInitWithDifferentTaxonomyTests(unittest.TestCase):
    def setUp(self):
        self.nwk_path = os.path.join(os.path.dirname(__file__), './data/tomato.nwk')
        self.tree_str = utils.get_newick_string(self.nwk_path, type="nwk")
        self.orthoxml_path = os.path.join(os.path.dirname(__file__), './data/tom.orthoxml')

    def test_init_with_newick_taxonomy(self):
        analysis = ham.Ham(tree_file=self.tree_str, hog_file=self.orthoxml_path)
        self.assertIsNotNone(analysis.taxonomy, "Taxonomy should be initialized")

    def test_init_with_orthoxml_taxonomy(self):
        analysis = ham.Ham(hog_file=self.orthoxml_path)
        self.assertIsNotNone(analysis.taxonomy, "Taxonomy should be initialized with internal names")

    def test_taxonomy_from_newick_identical(self):
        analysis_newick = ham.Ham(tree_file=self.tree_str, hog_file=self.orthoxml_path, use_internal_name=True)
        analysis_orthoxml = ham.Ham(hog_file=self.orthoxml_path)
        rf = analysis_orthoxml.taxonomy.tree.robinson_foulds(analysis_newick.taxonomy.tree)
        self.assertEqual(rf[0], 0, "Taxonomy from Newick and Orthoxml should be identical")


class AncestralMetazoaGenomeContentTest(unittest.TestCase):
    def setUp(self):
        data = Path(__file__).parent / 'data'
        nwk_path = data / 'euk-HOG0001687.nwk'
        orthoxml_path = data / 'euk-HOG0001687.orthoxml'
        self.ham_analysis = ham.Ham(tree_file=nwk_path, hog_file=orthoxml_path,
                                    tree_format="newick", use_internal_name=True)

    def test_ancestral_genome_empty_at_treeroot(self):
        anc_genome = self.ham_analysis.get_ancestral_genome_by_name("internal_0")
        self.assertEqual(len(anc_genome.genes), 0, "Ancestral genome at root should be empty")

    def test_ancestral_genome_at_hog_root(self):
        anc_genome = self.ham_analysis.get_ancestral_genome_by_name("42289")
        self.assertEqual(len(anc_genome.genes), 1, "Ancestral genome at level 42289 should have one HOG")
        self.assertIsNone(anc_genome.genes[0].parent)

    def test_level_7145_has_3_ancestral_genes(self):
        anc_genome = self.ham_analysis.get_ancestral_genome_by_name("7145")
        self.assertEqual(len(anc_genome.genes), 3, "Ancestral genome at level 7145 should have 3 HOGs")


if __name__ == "__main__":
    unittest.main()
