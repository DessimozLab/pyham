__author__ = 'admin'

from abc import ABCMeta, abstractmethod


class Genome(metaclass=ABCMeta):
    pass


class AncestralGenome(Genome):

    def __init__(self):
        self.taxon = None
        self.hogs = []
        self.parent = None
        self.children = []
        pass


class ExtantGenome(Genome):

    def __init__(self):
        self.taxon = None
        self.genes = []
        self.parent = None
        pass
