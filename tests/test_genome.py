__author__ = 'admin'

import unittest
from pyham import Gene, HOG, AbstractGene, EvolutionaryConceptError
from unittest import skip
import pyham.genome as g
import ete3


class ExtantGenomeTest(unittest.TestCase):

    def test_name_and_NCBITaxId_required(self):

        # missing input should raises error
        with self.assertRaises(TypeError):
            a = g.ExtantGenome()

        # missing NCBITaxId will put -1 as default
        b = g.ExtantGenome(name="HUMAN")
        self.assertEqual(b.taxid, "-1")

        # missing name should raises error
        with self.assertRaises(TypeError):
            c = g.ExtantGenome(NCBITaxId="9601")

        # valid
        self.assertIsInstance(g.ExtantGenome(name="HUMAN", NCBITaxId="9601"), g.ExtantGenome)


class GenomeTest(unittest.TestCase):

    def test_cannot_add_anything_as_taxon_range(self):
        a1 = g.ExtantGenome(name="HUMAN", NCBITaxId="9601")
        a2 = g.AncestralGenome()
        b = Gene(id="423")

        # wrong type of input
        with self.assertRaises(TypeError):
            a1.set_taxon(b)
        with self.assertRaises(TypeError):
            a2.set_taxon(b)

        # valid
        c = ete3.Tree("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", format=1)
        a1.set_taxon(c.get_tree_root())
        a2.set_taxon(c.get_tree_root())

        # try to wrongly reassign a taxonomic range to a extant genomes
        d = ete3.Tree("(E:1,D:1)Internal;", format=1)
        with self.assertRaises(EvolutionaryConceptError):
            a1.set_taxon(d.get_tree_root())
        with self.assertRaises(EvolutionaryConceptError):
            a2.set_taxon(d.get_tree_root())

        # reassign the same taxonomy
        a1.set_taxon(c.get_tree_root())
        a2.set_taxon(c.get_tree_root())


if __name__ == "__main__":
    unittest.main()