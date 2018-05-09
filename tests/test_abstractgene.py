import unittest
from pyham import Gene, HOG, AbstractGene, AncestralGenome, ExtantGenome, utils, ham
from pyham.abstractgene import EvolutionaryConceptError as ECE
import os

class GeneTest(unittest.TestCase):

    def test_id_required(self):
        with self.assertRaises(TypeError):
            Gene()
        self.assertIsInstance(Gene('324'), AbstractGene)

    def test_cannot_add_multiple_more_than_one_extant_genome(self):
        g = Gene(id="423")

        # invalid genome
        with self.assertRaises(TypeError):
            g.set_genome("wrong")

        # invalid ancestral genome
        with self.assertRaises(TypeError):
            a = AncestralGenome()
            g.set_genome(a)

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

    def test_remove_child(self):
        a = HOG(id="a")
        b = HOG(id="b")
        c = HOG(id="b")

        self.assertListEqual(a.children, [])

        a.add_child(b)
        a.add_child(c)
        self.assertListEqual(a.children, [b,c])

        a.remove_child(b)
        self.assertListEqual(a.children, [c])

    def test_only_one_genome_possible(self):
        h = HOG()

        # invalid genome
        with self.assertRaises(TypeError):
            h.set_genome("wrong")

        # valid ancestral genome
        a = AncestralGenome()
        h.set_genome(a)

        # invalid extant genome
        hum = ExtantGenome("HUMAN", "1")
        with self.assertRaises(TypeError):
            h.set_genome(hum)

        # cannot reassign a genome is already set
        b = AncestralGenome()
        with self.assertRaises(ECE):
            h.set_genome(b)

        # but same twice should be ok
        h.set_genome(a)

    def test_representation_of_Hog(self):
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


class AbstractGeneTest(unittest.TestCase):

    def setUp(self):
        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        nwk_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        self.h = ham.Ham(nwk_str, orthoxml_path, use_internal_name=True)

    def test_search_ancestor_hog_in_ancestral_genome(self):

        # there is an ancestor founded
        gene2 = self.h.get_gene_by_id("2")
        hog2 = self.h.get_hog_by_id("2")

        ancestor, is_dup = gene2.search_ancestor_hog_in_ancestral_genome(hog2.genome)
        self.assertEqual(ancestor, hog2)
        self.assertEqual(is_dup, False)

        # there is no ancestor founded
        hog1 = self.h.get_hog_by_id(1)

        ancestor, is_dup = gene2.search_ancestor_hog_in_ancestral_genome(hog1.genome)
        self.assertEqual(ancestor, None)
        self.assertEqual(is_dup, False)

        # there is a duplication
        gene3 = self.h.get_gene_by_id("3")
        hog3 = self.h.get_hog_by_id(3)

        ancestor, is_dup = gene3.search_ancestor_hog_in_ancestral_genome(hog3.genome)
        self.assertEqual(ancestor, hog3)
        self.assertEqual(is_dup, True)

        # level given is the same as hog
        hog3 = self.h.get_hog_by_id(3)

        ancestor, is_dup = hog3.search_ancestor_hog_in_ancestral_genome(hog3.genome)
        self.assertEqual(ancestor, None)
        self.assertEqual(is_dup, False)

    def test_get_top_level_hog(self):

        # there is an top level hog founded
        gene2 = self.h.get_gene_by_id("2")
        hog2 = self.h.get_hog_by_id("2")

        ancestor = gene2.get_top_level_hog()
        self.assertEqual(ancestor, hog2)

        # hog is already top level
        hog2 = self.h.get_hog_by_id("2")

        ancestor = hog2.get_top_level_hog()
        self.assertEqual(ancestor, hog2)

        # gene is singleton
        singleton = self.h.get_gene_by_id("5")

        ancestor = singleton.get_top_level_hog()
        self.assertEqual(ancestor, singleton)

    def test_get_at_level(self):
        gene1 = self.h.get_gene_by_id("1")

        hog3 = self.h.get_hog_by_id("3")
        hog2 = self.h.get_hog_by_id("2")

        vert = self.h.get_ancestral_genome_by_name("Vertebrata")
        rodents = self.h.get_ancestral_genome_by_name("Rodents")
        euarch = self.h.get_ancestral_genome_by_name("Euarchontoglires")

        # Level is itself genomes
        with self.assertRaises(KeyError):
            gene1.get_at_level(gene1.genome)
        with self.assertRaises(KeyError):
            hog3.get_at_level(hog3.genome)

        # Wrong input Type or not in this family
        with self.assertRaises(TypeError):
            gene1.get_at_level("")
        with self.assertRaises(KeyError):
            hog2.get_at_level(vert)

        gal_1 = gene1.get_at_level(rodents)[0]
        self.assertEqual(gal_1.genome, rodents)

        gal_2 = hog3.get_at_level(euarch)
        self.assertEqual(len(gal_2), 2)


if __name__ == "__main__":
    unittest.main()
