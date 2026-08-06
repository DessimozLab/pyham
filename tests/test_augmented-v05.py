import unittest
import pytest
from pyham import utils
from pyham import ham
from pyham.abstractgene import TaxonomicConflictError
import os

class AugmentedExampleTests(unittest.TestCase):

    def setUp(self):
        # embedded taxonomy is self-consistent by construction (no external tree given),
        # so this builds without any taxonomic conflict.
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/augmented-test-v05.orthoxml')
        self.ham_analysis = ham.Ham(hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True)

    def test_all_hogs_at_Eukaryota_have_taxid_2759_in_hogid(self):
        hogs_with_wrong_taxid = [hog for hog in self.ham_analysis.get_ancestral_genome_by_name("Eukaryota").genes if hog.hog_id.split('_')[1] != "2759"]
        self.assertEqual(hogs_with_wrong_taxid, [], "All HOGs at Eukaryota should have taxid 2759 in hog_id")


####
# The external species tree resolves Ogataea parapolymorpha in a lineage disjoint from
# "Saccharomycetales", conflicting with (almost) every family's own TaxRange claim that
# groups it there -- a systemic, whole-file taxonomic conflict.
def _augmented_v05_with_external_tree(**kwargs):
    nwk_path = os.path.join(os.path.dirname(__file__), './data/augmented-test-v05.nwk')
    tree_str = utils.get_newick_string(nwk_path, type="nwk")
    orthoxml_path = os.path.join(os.path.dirname(__file__), './data/augmented-test-v05.orthoxml')
    return ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml', use_internal_name=True, **kwargs)


def test_augmented_v05_conflicts_collected_by_default():
    with pytest.raises(TaxonomicConflictError) as exc:
        _augmented_v05_with_external_tree()

    conflicts = exc.value.conflicts
    assert len(conflicts) == 3359
    assert {c.kind for c in conflicts} == {"disjoint"}
    assert len({c.family_id for c in conflicts}) == 3237
    assert conflicts[0].family_id == "HOG:0002640"
    assert conflicts[0].resolved_level == "Saccharomycetales"


def test_augmented_v05_fail_fast_stops_at_first_conflict():
    with pytest.raises(TaxonomicConflictError) as exc:
        _augmented_v05_with_external_tree(fail_fast=True)

    conflicts = exc.value.conflicts
    assert len(conflicts) == 1
    assert conflicts[0].kind == "disjoint"
    assert conflicts[0].family_id == "HOG:0002640"


if __name__ == "__main__":
    unittest.main()
