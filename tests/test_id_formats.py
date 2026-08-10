import unittest
from xml.etree.ElementTree import XMLParser
from pyham import id_formats
from pyham import parsers


class FakeTaxon:
    def __init__(self, props=None, name=None):
        self.props = props or {}
        self.name = name


class DetectIdSchemeTest(unittest.TestCase):

    def test_loft_taxid_ids(self):
        samples = ["HOG:0000001.1b_2", "HOG:0000001.2a_10", "HOG:0402997.2b.1a_47424", "HOG:0402997", "HOG:F002251.53b.5az_9606"]
        self.assertEqual(id_formats.detect_id_scheme(samples), id_formats.LOFT_TAXID)

    def test_loft_ids_without_taxid(self):
        samples = ["HOG:00000001.1a.1a", "HOG:00000001.1b", "HOG:00000001"]
        self.assertEqual(id_formats.detect_id_scheme(samples), id_formats.LOFT)

    def test_generic_dotted_but_not_loft_shaped(self):
        # dots present, but not "<digits><letters>" duplication-branch pairs.
        samples = ["1.M.E.P", "2.E.P", "HOG_0010360_sub100"]
        self.assertEqual(id_formats.detect_id_scheme(samples), id_formats.GENERIC)

    def test_generic_plain_incrementing_ids(self):
        samples = ["1", "2", "3", "4", "5", "6"]
        self.assertEqual(id_formats.detect_id_scheme(samples), id_formats.GENERIC)

    def test_generic_when_no_nested_ids_seen(self):
        # only bare top-level ids (no subhog part) says nothing about the nested convention.
        samples = ["HOG:0000001", "HOG:0000002", "HOG:0000003"]
        self.assertEqual(id_formats.detect_id_scheme(samples), id_formats.GENERIC)

    def test_empty_input(self):
        self.assertEqual(id_formats.detect_id_scheme([]), id_formats.GENERIC)


class MakeMissingHogIdTest(unittest.TestCase):

    def test_loft_taxid_appends_taxid_from_id_prop(self):
        taxon = FakeTaxon(props={'id': 47424}, name='Marasmiineae')
        result = id_formats.make_missing_hog_id(id_formats.LOFT_TAXID, "HOG:0402997.2b.1a_10", taxon)
        self.assertEqual(result, "HOG:0402997.2b.1a_47424")

    def test_loft_taxid_falls_back_to_taxon_id_prop(self):
        taxon = FakeTaxon(props={'taxon_id': '99'}, name='Clade')
        result = id_formats.make_missing_hog_id(id_formats.LOFT_TAXID, "HOG:1.1a_1", taxon)
        self.assertEqual(result, "HOG:1.1a_99")

    def test_loft_taxid_without_numeric_taxid_keeps_anchor(self):
        # no 'id'/'taxon_id' prop at all (e.g. a plain newick-derived tree) -- returning a
        # non-unique id is preferable to fabricating one that can't be parsed back out.
        taxon = FakeTaxon(props={}, name='SomeClade')
        result = id_formats.make_missing_hog_id(id_formats.LOFT_TAXID, "HOG:1.1a_1", taxon)
        self.assertEqual(result, "HOG:1.1a_1")

    def test_loft_scheme_keeps_anchor_unchanged(self):
        taxon = FakeTaxon(props={'id': 47424}, name='Marasmiineae')
        result = id_formats.make_missing_hog_id(id_formats.LOFT, "HOG:0402997.2b.1a_10", taxon)
        self.assertEqual(result, "HOG:0402997.2b.1a_10")

    def test_generic_scheme_keeps_anchor_unchanged(self):
        taxon = FakeTaxon(props={'id': 47424}, name='Marasmiineae')
        result = id_formats.make_missing_hog_id(id_formats.GENERIC, "1", taxon)
        self.assertEqual(result, "1")

    def test_none_anchor_stays_none(self):
        taxon = FakeTaxon(props={'id': 1}, name='X')
        self.assertIsNone(id_formats.make_missing_hog_id(id_formats.LOFT_TAXID, None, taxon))

    def test_unparseable_anchor_returned_unchanged(self):
        taxon = FakeTaxon(props={'id': 1}, name='X')
        result = id_formats.make_missing_hog_id(id_formats.LOFT_TAXID, "not-a-hog-id!", taxon)
        self.assertEqual(result, "not-a-hog-id!")


class IDSchemeSnifferTest(unittest.TestCase):
    """A single family can hold arbitrarily many HOGs (e.g. a giant root family spanning
    most of a file) -- the sniffer must stop once it has enough samples without waiting
    for that family to close, or it defeats the point of bounding the sniff at all."""

    NS = "{http://orthoXML.org/2011/}orthologGroup"

    def _feed(self, sniffer, xml_fragments):
        parser = XMLParser(target=sniffer)
        for fragment in xml_fragments:
            parser.feed(fragment)
            if sniffer.done:
                break
        return sniffer

    def test_stops_mid_family_once_max_ids_reached(self):
        sniffer = parsers.IDSchemeSniffer()
        sniffer.MAX_IDS = 5
        fragments = ['<orthoXML xmlns="http://orthoXML.org/2011/"><groups>',
                     '<orthologGroup id="HOG:1_1">']
        fragments += ['<orthologGroup id="HOG:1.{0}a_{0}">'.format(i) for i in range(10)]
        self._feed(sniffer, fragments)
        self.assertTrue(sniffer.done)
        self.assertEqual(len(sniffer.samples), 5)

    def test_stops_after_max_families_of_small_families(self):
        sniffer = parsers.IDSchemeSniffer()
        sniffer.MAX_FAMILIES = 2
        fragments = ['<orthoXML xmlns="http://orthoXML.org/2011/"><groups>']
        for fam in range(5):
            fragments.append('<orthologGroup id="HOG:{0}_1"><geneRef id="1"/></orthologGroup>'.format(fam))
        self._feed(sniffer, fragments)
        self.assertTrue(sniffer.done)
        self.assertEqual(sniffer.families_seen, 2)
        self.assertEqual(sniffer.samples, ["HOG:0_1", "HOG:1_1"])


if __name__ == '__main__':
    unittest.main()
