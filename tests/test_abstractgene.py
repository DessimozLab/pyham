import unittest
from ham import Gene, HOG, AbstractGene, AncestralGenome, ExtantGenome

from ham.abstractgene import EvolutionaryConceptError as ECE

class GeneTest(unittest.TestCase):
    def test_id_required(self):
        with self.assertRaises(TypeError):
            Gene()
        self.assertIsInstance(Gene('324'), AbstractGene)

    def test_cannot_add_multiple_taxon_ranges(self):
        g = Gene(id="423")
        # invalid genome
        with self.assertRaises(TypeError):
            g.set_genome("wrong")
        # valid
        b = ExtantGenome(name="HUMAN", NCBITaxId="9601")
        g.set_genome(b)
        # cannot reassign a genome is already set
        c = ExtantGenome(name="MONDO", NCBITaxId="96")
        with self.assertRaises(ECE):
            g.set_genome(c)
        # but same twice should be ok
        g.set_genome(b)


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
        with self.assertRaises(ECE):
            a.add_child(a)

    def test_only_one_genome_possible(self):
        h = HOG()
        # invalid genome
        with self.assertRaises(TypeError):
            h.set_genome("wrong")
        # valid
        a = AncestralGenome()
        h.set_genome(a)
        # cannot reassign a genome is already set
        b = AncestralGenome()
        with self.assertRaises(ECE):
            h.set_genome(b)
        # but same twice should be ok
        h.set_genome(a)


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


if __name__ == "__main__":
    unittest.main()