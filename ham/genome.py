__author__ = 'admin'

from abc import ABCMeta, abstractmethod


class Genome(metaclass=ABCMeta):
    def __init__(self):
        self.genes = []
        self.taxon = None

    def add_gene(self, gene):
        self.genes.append(gene)
        gene.set_genome(self)



class AncestralGenome(Genome):

    def __init__(self):
        super(AncestralGenome, self).__init__()
        self.parent = None
        self.children = []


class ExtantGenome(Genome):

    def __init__(self, name, NCBITaxId, **kwargs):
        super(ExtantGenome, self).__init__()
        self.parent = None
        self.name = name
        self.taxid = NCBITaxId
