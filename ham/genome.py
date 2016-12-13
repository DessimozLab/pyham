__author__ = 'admin'

from abc import ABCMeta


class Genome(metaclass=ABCMeta):
    pass


class AncestralGenome(Genome):

    def __init__(self):
        pass


class ExtantGenome(Genome):

    def __init__(self):
        pass
