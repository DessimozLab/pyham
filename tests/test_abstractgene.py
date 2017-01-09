import unittest
from ham import Gene, HOG, AbstractGene, EvolutionaryConceptError


class GeneTest(unittest.TestCase):
    def test_id_required(self):
        with self.assertRaises(TypeError):
            Gene()
        self.assertIsInstance(Gene('324'), AbstractGene)

    def test_cannot_add_multiple_taxon_ranges(self):
        g = Gene(id="423")
        g.set_genome("test")
        with self.assertRaises(EvolutionaryConceptError):
            g.set_genome("bla")
        # but same twice should be ok
        g.set_genome("test")


class HogTest(unittest.TestCase):
    def test_hog_can_be_nested_and_traversed_up_and_down(self):
        a = HOG(id="a")
        b = HOG(id="b")
        a.add_child(b)
        self.assertIn(b, b.parent.children)

    def test_cannot_add_any_type_as_child(self):
        a = HOG()
        with self.assertRaises(TypeError):
            a.add_child("test")

    def test_cannot_be_child_of_self(self):
        a = HOG()
        with self.assertRaises(EvolutionaryConceptError):
            a.add_child(a)

    def test_several_taxon_ranges_possible(self):
        a = HOG()
        taxranges = ("Primates", "Mammalia")
        for tax in taxranges:
            a.set_genome(tax)
        self.assertEqual(set(taxranges), a.genome)

    def test_represenation_of_Hog(self):
        a = HOG(id="441.2a", bla="don't know")
        self.assertEqual("{!r}".format(a), "<HOG(441.2a)>")

    def test_add_score(self):
        a = HOG()
        a.score('testscore', 0.932)
        self.assertAlmostEqual(a.score('testscore'), 0.932)

    def test_unknown_score_raises_keyerror(self):
        a = HOG()
        with self.assertRaises(KeyError):
            a.score('testscore')