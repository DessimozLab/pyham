__author__ = 'admin'

from abc import ABCMeta


class AbstractGene(metaclass=ABCMeta):
    pass


class HOG(AbstractGene):

    def __init__(self):
        self.hog_id = None
        self.parent = None
        self.taxon = None
        self.depth = 0
        self.children = []


class Gene(AbstractGene):

    def __init__(self):
        self.unique_id = None
        self.parent = None
        self.species = None
        self.mapping = {}
