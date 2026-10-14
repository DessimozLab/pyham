import contextlib
import io
import os
import tempfile
import unittest

from pyham import Ham, utils
from pyham.abstractgene import EvolutionaryConceptError
from pyham.validation import validate_input, InputValidationError, main


DATA = os.path.join(os.path.dirname(__file__), 'data')


def data(*parts):
    return os.path.join(DATA, *parts)


def codes(report):
    return {issue.code for issue in report.issues}


def issue(report, code):
    found = [i for i in report.issues if i.code == code]
    assert len(found) == 1, f"expected one {code!r} issue, got {[i.code for i in report.issues]}"
    return found[0]


SIMPLE_NWK = utils.get_newick_string(data('simpleEx.nwk'), type="nwk")
with open(data('simpleEx.orthoxml')) as fh:
    SIMPLE_ORTHOXML = fh.read()


def orthoxml(species, groups):
    """species: {name: [gene ids]}, groups: the inner XML of <groups>."""
    blocks = []
    for name, genes in species.items():
        gene_xml = "".join(f'<gene id="{g}" protId="{name}{g}"/>' for g in genes)
        blocks.append(f'<species name="{name}" NCBITaxId="1"><database name="d" version="1">'
                      f'<genes>{gene_xml}</genes></database></species>')
    return ('<?xml version="1.0" encoding="UTF-8"?>'
            '<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="test" originVersion="1">'
            + "".join(blocks) + f"<groups>{groups}</groups></orthoXML>")


class ValidInputTest(unittest.TestCase):

    def test_clean_input_has_no_errors(self):
        report = validate_input(SIMPLE_NWK, data('simpleEx.orthoxml'), use_internal_name=True)
        self.assertTrue(report.ok)
        self.assertEqual(report.errors, [])
        self.assertEqual(report.warnings, [])
        self.assertEqual(report.stats['families'], 3)
        self.assertEqual(report.stats['families_built'], 3)
        self.assertIn("all 3 families built", report.trial_build)
        report.raise_for_errors()

    def test_report_text(self):
        text = str(validate_input(SIMPLE_NWK, data('simpleEx.orthoxml'), use_internal_name=True))
        self.assertIn("pyHam input check: 0 errors, 0 warnings, 1 note", text)
        self.assertIn("6 species, 19 genes, 3 families", text)
        self.assertIn("singleton-genes", text)


class SpeciesTreeTest(unittest.TestCase):

    def test_unreadable_newick(self):
        report = validate_input("((Synechococcus sp. (strain JA-2-3B'a),B),C);", data('simpleEx.orthoxml'))
        self.assertIn('tree-unreadable', codes(report))
        self.assertIn("must be quoted", issue(report, 'tree-unreadable').hint)
        self.assertIn("skipped", report.trial_build)

    def test_support_values_as_internal_names(self):
        # issue #26: numeric internal labels are support values, not taxon names
        report = validate_input("((A,B)1,(C,D)1)1;", data('simpleEx.orthoxml'), use_internal_name=True,
                                trial_build=False)
        dup = issue(report, 'duplicate-internal-names')
        self.assertEqual(dup.items, ['"1" (3 times)'])
        self.assertIn("support values", dup.hint)
        self.assertIn("use_internal_name=False", dup.hint)

    def test_phyloxml_without_taxonomy_elements(self):
        # issue #25: names only in <clade><name>, but the default tag reads <taxonomy><scientific_name>
        phyloxml = ('<phyloxml xmlns="http://www.phyloxml.org"><phylogeny rooted="true"><clade><name>Mammalia</name>'
                    '<clade><name>HUMAN</name></clade><clade><name>MOUSE</name></clade></clade></phylogeny></phyloxml>')
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'tree.phyloxml')
            with open(path, 'w') as fh:
                fh.write(phyloxml)
            report = validate_input(path, data('simpleEx.orthoxml'), use_internal_name=True, trial_build=False)
        unnamed = issue(report, 'unnamed-leaves')
        self.assertEqual(unnamed.items, ["leaf with clade_name='HUMAN'", "leaf with clade_name='MOUSE'"])
        self.assertIn("phyloxml_leaf_name_tag='clade_name'", unnamed.hint)
        self.assertIn("phyloxml_internal_name_tag='clade_name'", issue(report, 'unnamed-internal-nodes').hint)


class SpeciesMappingTest(unittest.TestCase):

    def test_names_differing_only_in_case_or_whitespace(self):
        xml = SIMPLE_ORTHOXML.replace('name="HUMAN"', 'name="Human"').replace('name="PANTR"', 'name="PANTR "')
        report = validate_input(SIMPLE_NWK, xml, orthoXML_as_string=True, use_internal_name=True)
        missing = issue(report, 'species-not-in-tree')
        self.assertEqual(missing.data['species'], ["Human", "PANTR "])
        self.assertEqual(missing.data['suggestions'], {"Human": "HUMAN", "PANTR ": "PANTR"})
        self.assertIn("letter case", missing.hint)
        self.assertIn("skipped", report.trial_build)

    def test_ham_reports_all_missing_species_at_once(self):
        xml = SIMPLE_ORTHOXML.replace('name="HUMAN"', 'name="Human"').replace('name="PANTR"', 'name="PANTR "')
        with self.assertRaises(InputValidationError) as ctx:
            Ham(SIMPLE_NWK, xml, orthoXML_as_string=True, use_internal_name=True)
        self.assertIsInstance(ctx.exception, KeyError)  # what Ham raised before, kept for compatibility
        message = str(ctx.exception)
        self.assertIn('"Human": the tree has "HUMAN"', message)
        self.assertIn('"PANTR ": the tree has "PANTR"', message)
        self.assertIn("pyham.validate_input()", message)
        self.assertNotIn("Trial build", message)

    def test_species_missing_from_pruned_tree(self):
        report = validate_input(data('issue25', 'speciestree.nwk'), data('issue25', 'hog2.orthoxml'),
                                use_internal_name=True)
        missing = issue(report, 'species-not-in-tree')
        self.assertEqual(missing.data['species'], ["Syncephalastrum_racemosum"])
        self.assertIn("skipped", report.trial_build)


class LevelTest(unittest.TestCase):

    def test_levels_missing_from_pruned_tree(self):
        # issue #25's data: the tree lacks Mucorales/Mucorineae, so a duplication ends up above the root
        report = validate_input(data('issue25', 'speciestree.nwk'), data('issue25', 'hog.orthoxml'),
                                use_internal_name=True)
        levels = issue(report, 'level-not-in-tree')
        self.assertEqual(levels.data['levels'], {"Mucorineae": 1, "Mucorales": 1})
        failed = issue(report, 'family-build-failed')
        self.assertIn("would predate the tree", failed.items[0])

    def test_suggests_use_internal_name(self):
        tree = data('issue25', 'species_tree.nwk_scientific.nwk')
        report = validate_input(tree, data('issue25', 'hog.orthoxml'))
        levels = issue(report, 'level-not-in-tree')
        self.assertIn("use_internal_name=True", levels.hint)

        report = validate_input(tree, data('issue25', 'hog.orthoxml'), use_internal_name=True)
        self.assertTrue(report.ok, str(report))
        self.assertNotIn('level-not-in-tree', codes(report))

    def test_no_level_annotations_note(self):
        report = validate_input(data('issue19', 'insecta.nwk'), data('issue19', 'eglps.orthoxml'))
        self.assertEqual(issue(report, 'no-level-annotations').severity, 'note')


class TrialBuildTest(unittest.TestCase):

    def test_unreconciled_gene_tree(self):
        # issue #19: an orthologGroup puts PHUM outside (LIPO, MGAL), but PHUM and MGAL are sisters in the tree
        report = validate_input(data('issue19', 'insecta.nwk'), data('issue19', 'eglps.orthoxml'))
        failed = issue(report, 'family-build-failed')
        self.assertEqual(len(failed.items), 1)
        self.assertIn("contradicts the species tree", failed.items[0])
        self.assertIn("reconcile", failed.hint)

    def test_ham_names_family_and_explains(self):
        with self.assertRaises(EvolutionaryConceptError) as ctx:
            Ham(data('issue19', 'insecta.nwk'), data('issue19', 'eglps.orthoxml'), tree_format='newick')
        self.assertIn("contradicts the species tree", str(ctx.exception))
        if hasattr(ctx.exception, '__notes__'):
            self.assertIn("family #1", ctx.exception.__notes__[0])

    def test_taxonomic_conflicts_across_whole_file(self):
        with open(data('augmented-test-v05.nwk')) as fh:
            tree = fh.read()
        report = validate_input(tree, data('augmented-test-v05.orthoxml'), use_internal_name=True)
        conflicts = issue(report, 'taxonomic-conflicts')
        self.assertEqual(len({c.family_id for c in conflicts.data['conflicts']}), 3237)
        self.assertEqual(report.stats['families_built'], report.stats['families'] - 3237)
        self.assertEqual(issue(report, 'same-level-hogs-collapsed').severity, 'note')


class OrthoXMLContentTest(unittest.TestCase):

    def test_gene_and_species_integrity(self):
        xml = orthoxml({"HUMAN": ["1", "2"], "MOUSE": ["3", "2"], "PANTR": ["5"], "CANFA": ["6"]},
                       '<orthologGroup id="1"><geneRef id="1"/><geneRef id="3"/></orthologGroup>'
                       '<orthologGroup id="2"><geneRef id="1"/><geneRef id="99"/></orthologGroup>')
        xml = xml.replace('<species name="CANFA"', '<species name="PANTR"')
        report = validate_input(SIMPLE_NWK, xml, orthoXML_as_string=True, use_internal_name=True)
        self.assertEqual(issue(report, 'duplicate-gene-ids').items, ["gene id 2 (in HUMAN, and again in MOUSE)"])
        self.assertEqual(issue(report, 'unknown-gene-reference').items, ["geneRef 99 in family 2"])
        self.assertEqual(issue(report, 'gene-referenced-twice').items, ["gene 1 (HUMAN) in families 1, 2"])
        self.assertEqual(issue(report, 'duplicate-species').items, ['"PANTR" (2 times)'])
        self.assertIn("skipped", report.trial_build)

    def test_unreadable_orthoxml(self):
        report = validate_input(SIMPLE_NWK, "<orthoXML><groups>", orthoXML_as_string=True)
        self.assertIn('orthoxml-unreadable', codes(report))


class CommandLineTest(unittest.TestCase):

    def run_main(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            status = main(list(argv))
        return status, out.getvalue()

    def test_exit_status(self):
        status, out = self.run_main(data('issue25', 'species_tree.nwk_scientific.nwk'),
                                    data('issue25', 'hog.orthoxml'), '--use-internal-name')
        self.assertEqual(status, 0, out)
        status, out = self.run_main(data('issue19', 'insecta.nwk'), data('issue19', 'eglps.orthoxml'))
        self.assertEqual(status, 1)
        self.assertIn("family-build-failed", out)


if __name__ == "__main__":
    unittest.main()
