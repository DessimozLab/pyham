__author__ = 'admin'

from abc import ABCMeta


class Genome(metaclass=ABCMeta):
    def __init__(self):
        self.genes = []
        self.taxon = None

    def add_gene(self, gene):
        self.genes.append(gene)
        gene.set_taxon_range(self)


class AncestralGenome(Genome):

    def __init__(self):
        self.parent = None
        self.children = []


class ExtantGenome(Genome):

    def __init__(self, name, NCBITaxId, **kwargs):
        super(ExtantGenome, self).__init__()
        self.parent = None
        self.name = name
        self.taxid = NCBITaxId
