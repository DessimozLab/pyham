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


####
# tests for HOG id scheme auto-detection and missing-level id synthesis under LOFT_TAXID
# (regression test for pyham/id_formats.py). `loft_taxid_gap.orthoxml` embeds its own
# taxonomy (GrandRoot > Root > Mid > SubMid > SubSubMid > SP1/SP2, with SP3 a sibling of
# Mid under Root) and a single duplicated family whose HOGs only explicitly claim
# GrandRoot and Mid -- so both the "several stacked missing levels" gap (Mid -> SubMid ->
# SubSubMid) and the "two duplication branches missing the same level" gap (both branches
# missing Root) get synthesized.
def _all_hogs(ham_obj):
    from pyham import abstractgene
    hogs = []
    for taxon in ham_obj.taxonomy.tree.traverse():
        genome = taxon.props.get('genome')
        if genome is None:
            continue
        hogs.extend(g for g in genome.genes if isinstance(g, abstractgene.HOG))
    return hogs


def test_loft_taxid_scheme_is_auto_detected():
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/loft_taxid_gap.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)
    assert h.id_schema == "LOFT_TAXID"


def test_loft_taxid_missing_levels_get_distinct_resolvable_ids():
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/loft_taxid_gap.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)

    hogs = _all_hogs(h)
    ids = [hog.hog_id for hog in hogs]
    assert len(ids) == len(set(ids)), "every HOG (real or synthesized) should have a unique id"

    assert "HOG:0000005_50" in ids  # GrandRoot, real (top-level)
    assert "HOG:0000005.1a_200" in ids  # Mid, real
    # stacked missing levels below Mid each get their own distinct, taxid-suffixed id
    assert "HOG:0000005.1a_300" in ids  # SubMid
    assert "HOG:0000005.1a_350" in ids  # SubSubMid
    # both duplication branches are missing "Root"; each gets its own synthesized hog,
    # correctly labeled as sibling branches (.1a/.1b) of the same duplication rather
    # than colliding on a shared, dot-chain-less id
    assert "HOG:0000005.1a_100" in ids  # Root, ancestor of the real .1a/Mid lineage
    assert "HOG:0000005.1b_100" in ids  # Root, wrapper for the bare-geneRef duplication branch

    # every real or synthesized id resolves through the public lookup API.
    for hog_id in ids:
        assert h.get_hog_by_id(hog_id).hog_id == hog_id


def test_get_hog_by_id_resolves_bare_fam_taxid_id():
    # regression test: get_hog_by_id used to only delegate to the taxid-aware
    # HOG.find_by_id when the queried id had a dot-chain (`subhog`); a bare
    # "<fam>_<taxid>" id (no dot-chain, e.g. a synthesized level with no
    # duplication above it) silently fell through to `return roothog`, ignoring
    # the taxid entirely and returning the wrong HOG. Uses a small inline fixture
    # (no duplication at all) so this stays decoupled from how duplication branches
    # get labeled elsewhere.
    orthoxml = """<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="pyham test fixture" originVersion="0.1">
 <species name="SP1" NCBITaxId="1"><database name="SP1fake" version="0.1"><genes><gene id="1" protId="SP1_1" geneId="SP1g1" /></genes></database></species>
 <species name="SP2" NCBITaxId="2"><database name="SP2fake" version="0.1"><genes><gene id="2" protId="SP2_1" geneId="SP2g1" /></genes></database></species>
 <species name="SP3" NCBITaxId="3"><database name="SP3fake" version="0.1"><genes><gene id="3" protId="SP3_1" geneId="SP3g1" /></genes></database></species>
 <taxonomy>
  <taxon id="10" name="Family">
   <taxon id="20" name="Mid">
    <taxon id="1" name="SP1"/>
    <taxon id="2" name="SP2"/>
   </taxon>
   <taxon id="3" name="SP3"/>
  </taxon>
 </taxonomy>
 <groups>
  <orthologGroup id="HOG:0000009_10" taxonId="10">
   <geneRef id="1"/>
   <geneRef id="2"/>
   <geneRef id="3"/>
  </orthologGroup>
 </groups>
</orthoXML>
"""
    h = ham.Ham(hog_file=orthoxml, orthoXML_as_string=True, use_internal_name=True, id_schema="LOFT_TAXID")

    hog = h.get_hog_by_id("HOG:0000009_20")
    assert hog.hog_id == "HOG:0000009_20"
    assert hog.genome.name == "Mid"

    root_hog = h.get_hog_by_id("HOG:0000009_10")
    assert root_hog.genome.name == "Family"
    assert root_hog is not hog


def test_get_hog_by_id_with_dot_chain_but_no_taxid_does_not_crash():
    # regression test: querying a dot-chain id with no taxid suffix used to crash
    # with TypeError (int(None)) instead of doing a best-effort, taxid-less lookup.
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/loft_taxid_gap.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)

    hog = h.get_hog_by_id("HOG:0000005.1a")
    assert hog.hog_id.startswith("HOG:0000005.1a")


def test_generic_scheme_reproduces_old_colliding_behavior():
    # id_schema='GENERIC' (today's pre-existing "no adjustment" behavior) collapses the
    # same gaps down to just 2 distinct ids across 6 HOGs -- the exact bug this feature
    # fixes for files that actually use LOFT_TAXID.
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/loft_taxid_gap.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True, id_schema="GENERIC")

    ids = [hog.hog_id for hog in _all_hogs(h)]
    assert len(ids) == 6


####
# tests for backfilling LOFT dot-chain ids on duplication branches that have no id of
# their own (a bare geneRef, or an orthologGroup element missing its id attribute) --
# regression tests for OrthoXMLParser._assign_duplication_branch_ids, using two real
# FastOMA-exported families that exhibit the bug at different scales.
def _assert_all_hog_ids_unique(h):
    ids = [hog.hog_id for hog in _all_hogs(h)]
    duplicates = {i for i in ids if ids.count(i) > 1}
    assert not duplicates, f"colliding HOG ids: {duplicates}"
    return ids


def test_paralog_group_some_child_no_id_has_no_collisions():
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/paralog_group_some_child_no_id.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)
    assert h.id_schema == "LOFT_TAXID"

    ids = _assert_all_hog_ids_unique(h)
    for hog_id in ids:
        assert h.get_hog_by_id(hog_id).hog_id == hog_id


def test_paralog_group_second_maize_copy_becomes_sibling_branch():
    # the motivating example: a bare geneRef for a second maize paralog, sibling to the
    # real HOG:0000001.3a.1a_71, must become .3a.1b_71 -- same duplication (idx 1), next
    # letter, same taxid (both branches share the exact same duplication origin) -- not
    # a confusing near-duplicate of its sibling's id.
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/paralog_group_some_child_no_id.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)

    gene = h.extant_gene_map['1051020597']
    assert gene.parent.hog_id == "HOG:0000001.3a.1b_71"
    assert gene.parent.genome.name == "NODE_51"


def test_paralog_group_zero_depth_branches_get_loft_only_labels():
    # bare genes that are immediate children of the duplication (no vertical gap at
    # all) get labeled directly via set_LOFT (dot-chain only, no taxid -- matching how
    # real LOFT="..." geneRef attributes look), instead of staying unlabeled.
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/paralog_group_some_child_no_id.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)

    labels = {gid: getattr(h.extant_gene_map[gid], 'hog_id', None)
              for gid in ('1009021332', '1009015113', '1009024035')}
    assert all(label is not None for label in labels.values())
    assert len(set(labels.values())) == 3, "each sibling branch should get its own label"
    for label in labels.values():
        assert '_' not in label, "a zero-depth branch label should not carry a taxid suffix"


def test_nested_duplication_without_ids_has_no_collisions():
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/nested_duplication_without_ids.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)
    assert h.id_schema == "LOFT_TAXID"

    ids = _assert_all_hog_ids_unique(h)
    for hog_id in ids:
        assert h.get_hog_by_id(hog_id).hog_id == hog_id


def test_nested_duplication_third_sibling_branch_gets_next_letter():
    # HOG:0078342.6a_124 and HOG:0078342.6b_128 are real siblings of duplication "6";
    # a bare geneRef third sibling must become .6c, not collide with either.
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/nested_duplication_without_ids.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)

    gene = h.extant_gene_map['1028024909']
    assert gene.parent.hog_id == "HOG:0078342.6c_126"


def test_nested_duplication_fresh_index_avoids_real_sibling_index():
    # HOG:0078342.6a_124 directly contains two independent paralogGroups: one entirely
    # unlabeled, one with a real sibling id HOG:0078342.6a.5b_126 (revealing index 5).
    # The fully-unlabeled one must mint its own fresh index rather than colliding with 5.
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/nested_duplication_without_ids.orthoxml')
    h = ham.Ham(hog_file=orthoxml_path, use_internal_name=True)

    labels = [getattr(h.extant_gene_map[gid], 'hog_id', None)
              for gid in ('1029010392', '1029023514', '1029006173', '1029006174')]
    assert all(label is not None for label in labels)
    assert len(set(labels)) == 4
    for label in labels:
        idx = int(label.split('.')[1].rstrip('abcdefghijklmnopqrstuvwxyz'))
        assert idx != 5, "a freshly-minted index must not collide with the real sibling index 5"


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
