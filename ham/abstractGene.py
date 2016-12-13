__author__ = 'admin'

from abc import ABCMeta
from xml.etree.ElementTree import XMLParser
from . import parsers


class AbstractGene(metaclass=ABCMeta):
    @classmethod
    def create_abstract_gene_from_orthoxml(cls, orthoxml):
        target = parsers.orthoxmlParser()
        parser = XMLParser(target=target)
        with open(orthoxml, 'rt') as f:
            for line in f:
                parser.feed(line)


class HOG(AbstractGene):
    instances = []

    def __init__(self):
        self.hog_id = None
        self.parent = None
        self.taxon = None
        self.depth = 0
        self.children = []
        HOG.instances.append(self)


class Gene(AbstractGene):
    instances = []

    def __init__(self):
        self.unique_id = None
        self.parent = None
        self.species = None
        self.mapping = {}
        Gene.instances.append(self)
