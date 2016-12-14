__author__ = 'admin'

from abc import ABCMeta, abstractmethod


class AbstractGene(metaclass=ABCMeta):
    def __init__(self):
        self.parent = None
        self.taxon = None

    @abstractmethod
    def set_taxon_range(self, tax):
        """Set the taxonomic range"""
        pass


class HOG(AbstractGene):

    def __init__(self, id=None, **kwargs):
        super(HOG, self).__init__()
        self.hog_id = id
        self.depth = 0
        self.children = []

    def add_child(self, elem):
        if not isinstance(elem, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(elem).__name__))
        self.children.append(elem)
        elem.parent = self

    def append(self, elem):
        self.add_child(elem)

    def set_taxon_range(self, tax):
        """A HOG can potentially belong to several taxonomic ranges"""
        if self.taxon is None:
            self.taxon = set([])
        self.taxon.add(tax)


class Gene(AbstractGene):
    def __init__(self, id, geneId=None, protId=None, transcriptId=None, **kwargs):
        super(Gene, self).__init__()
        self.unique_id = id
        self.gene_id = geneId
        self.prot_id = protId
        self.transcript_id = transcriptId

    def set_taxon_range(self, tax):
        if self.taxon is not None and tax != self.taxon:
            raise KeyError("Gene can only belong to one genome")
        self.taxon = tax


