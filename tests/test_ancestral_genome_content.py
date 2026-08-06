import unittest
import pytest
from pyham import utils
from pyham import ham
from pyham.abstractgene import TaxonomicConflictError
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


@pytest.fixture(params=[
     ("fam_402997.augmented.orthoxml", None, "HOG:0402997_5302", "Agaricomycotina"),
     ("fam_402997.augmented.orthoxml", "fam_402997.nwk", "HOG:0402997_5302", "Agaricomycotina"),
     ("fam_800112.augmented.orthoxml", None, "HOG:0800112_7711", "Chordata"),
     ("fam_800112.augmented.orthoxml", "fam_402997.nwk", "HOG:0800112_7711", "Chordata"),
 ])
def hog_fam_ham(request):
     oxml, nwk, toplevel_hogid, root_level = request.param
     orthoxml_path = os.path.join(os.path.dirname(__file__), './data', oxml)
     tree_path = os.path.join(os.path.dirname(__file__), './data', nwk) if nwk else None
     ham_obj = ham.Ham(tree_file=tree_path, hog_file=orthoxml_path, tree_format="newick", use_internal_name=True)
     return ham_obj, toplevel_hogid, root_level


def test_root_level_as_expected(hog_fam_ham):
     ham_obj, toplevel_hogid, root_level = hog_fam_ham
     assert len(ham_obj.get_list_top_level_hogs()) == 1, "HOG 1074943 should have one top level hog"
     hog = ham_obj.get_list_top_level_hogs()[0]
     assert hog.hog_id == toplevel_hogid, "Top level HOG ID should match"
     assert hog.genome.name == root_level, f"Top level HOG level should be {root_level}"


####
# tests on hog_1074943 (that failed while building edgehog)
@pytest.fixture
def hog_1074943_ham():
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/hog_1074943.augmented.orthoxml')
    return ham.Ham(hog_file=orthoxml_path, use_internal_name=True)


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


####
# hog_1074943 against the external, more finely-resolved species tree: this exposes two
# real taxonomic conflicts between the orthoxml's own TaxRange claims and the given tree
# (see pyham/parsers.py's TaxonomicConflict machinery).
def _hog_1074943_with_external_tree(**kwargs):
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/hog_1074943.augmented.orthoxml')
    tree_path = os.path.join(os.path.dirname(__file__), './data/hog_1074943.nwk')
    return ham.Ham(tree_file=tree_path, hog_file=orthoxml_path, tree_format="newick", use_internal_name=True, **kwargs)


def test_hog_1074943_conflicts_collected_by_default():
    with pytest.raises(TaxonomicConflictError) as exc:
        _hog_1074943_with_external_tree()

    conflicts = exc.value.conflicts
    assert len(conflicts) == 2
    assert {c.kind for c in conflicts} == {"same_rank", "disjoint"}

    same_rank = next(c for c in conflicts if c.kind == "same_rank")
    assert same_rank.family_id == "HOG:1074943"
    assert same_rank.hog_id == "HOG:1074943.2a.7b.4a_35718"
    assert same_rank.resolved_level == "Chaetomiaceae"
    assert [label for label, species, detail in same_rank.offending_members] == ["HOG:1074943.2a.7b.4a_5149"]
    assert [label for label, species in same_rank.sibling_members] == ["HOG:1074943.2a.7b.4a_78579"]

    disjoint = next(c for c in conflicts if c.kind == "disjoint")
    assert disjoint.family_id == "HOG:1074943"
    assert disjoint.resolved_level == "Aspergillus subgen. Circumdati"
    assert [label for label, species, detail in disjoint.offending_members] == ["ASPPS06251"]


def test_hog_1074943_fail_fast_stops_at_first_conflict():
    with pytest.raises(TaxonomicConflictError) as exc:
        _hog_1074943_with_external_tree(fail_fast=True)

    conflicts = exc.value.conflicts
    assert len(conflicts) == 1
    assert conflicts[0].kind == "same_rank"
    assert conflicts[0].hog_id == "HOG:1074943.2a.7b.4a_35718"


####
# a small synthetic fixture for the "inverted" conflict kind (a HOG's claimed level is
# itself a descendant of one of its children's claimed level) -- no real example of this
# is known, unlike "same_rank" and "disjoint" above.
def test_inverted_level_conflict():
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/inverted_level.orthoxml')
    tree_path = os.path.join(os.path.dirname(__file__), './data/inverted_level.nwk')

    with pytest.raises(TaxonomicConflictError) as exc:
        ham.Ham(tree_file=tree_path, hog_file=orthoxml_path, tree_format="newick", use_internal_name=True)

    conflicts = exc.value.conflicts
    assert len(conflicts) == 1
    conflict = conflicts[0]
    assert conflict.kind == "inverted"
    assert conflict.family_id == "HOG:0000001"
    assert conflict.hog_id == "HOG:0000001_sub"
    assert conflict.resolved_level == "Sub"
    assert conflict.offending_members == [("HOG:0000001_family_inner", ["SP3"], None)]
    assert {label for label, species in conflict.sibling_members} == {"SP1_1", "SP2_1"}


####
# tests for HOG children that share a missing intermediate taxonomic level, i.e.
# the reference species tree resolves part of the taxonomy more finely than the
# orthoxml's own HOG structure (regression test for the shared_missing_level fixture)
@pytest.fixture
def shared_missing_level_ham():
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/shared_missing_level.orthoxml')
    tree_path = os.path.join(os.path.dirname(__file__), './data/shared_missing_level.nwk')
    return ham.Ham(tree_file=tree_path, hog_file=orthoxml_path, tree_format="newick", use_internal_name=True)


def test_shared_missing_level_creates_single_shared_hog(shared_missing_level_ham):
    h = shared_missing_level_ham

    family_hog = h.get_list_top_level_hogs()[0]
    assert {c.genome.name for c in family_hog.children} == {"SP6", "DeepA"}

    deep_a = h.get_ancestral_genome_by_name("DeepA")
    assert len(deep_a.genes) == 1, "DeepA should hold exactly one HOG for this family, not one per child"
    deep_a_hog = deep_a.genes[0]
    assert deep_a_hog.parent is family_hog
    assert {c.genome.name for c in deep_a_hog.children} == {"SP4", "Clade", "DeepB"}
    assert deep_a_hog._missing_in_xml.name == "Family"

    clade = h.get_ancestral_genome_by_name("Clade")
    assert len(clade.genes) == 1, "Clade should hold exactly one HOG for this family, not one per child"
    clade_hog = clade.genes[0]
    assert clade_hog.parent is deep_a_hog
    assert {c.genome.name for c in clade_hog.children} == {"SP1", "SP2", "SP3"}
    assert clade_hog._missing_in_xml.name == "DeepA"

    deep_b = h.get_ancestral_genome_by_name("DeepB")
    assert len(deep_b.genes) == 1
    deep_b_hog = deep_b.genes[0]
    assert deep_b_hog.parent is deep_a_hog
    assert len(deep_b_hog.children) == 1 and deep_b_hog.children[0].genome.name == "SP5"
    assert deep_b_hog._missing_in_xml.name == "DeepA"


def test_shared_missing_level_follows_species_tree(shared_missing_level_ham):
    for gene in shared_missing_level_ham.get_list_extant_genes():
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
