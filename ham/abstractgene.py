__author__ = 'admin'

from abc import ABCMeta, abstractmethod
import numbers


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
        self.children = []

    def add_child(self, elem):
        if not isinstance(elem, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(elem).__name__))
        if self == elem:
            raise EvolutionaryConceptError('Cannot be a child of itself')
        self.children.append(elem)
        elem.parent = self

    def append(self, elem):
        self.add_child(elem)

    def score(self, score_id, value=None):
        """get or set a score attribute for this HOG.

        If no value is passed, the stored value is returned or a KeyError is raised
        if this score has not been stored.

        :param score_id: score to be set or returned.
        :param value: (optional) A numeric value for the score of this HOG.
        :raises KeyError: if score_id is unknown and accessed."""
        if value is None:
            try:
                return self.scores[score_id]
            except (AttributeError, KeyError):
                raise KeyError("Score '{}' not set".format(score_id))

        if not isinstance(value, numbers.Number):
            raise ValueError("value must be numeric")
        scores = getattr(self, 'scores', {})
        scores[score_id] = value
        self.scores = scores

    def set_taxon_range(self, tax):
        """A HOG can potentially belong to several taxonomic ranges"""
        if self.taxon is None:
            self.taxon = set([])
        self.taxon.add(tax)

    def __repr__(self):
        return "<{}({})>".format(self.__class__.__name__, self.hog_id if self.hog_id else "")


class Gene(AbstractGene):
    def __init__(self, id, geneId=None, protId=None, transcriptId=None, **kwargs):
        super(Gene, self).__init__()
        self.unique_id = id
        self.gene_id = geneId
        self.prot_id = protId
        self.transcript_id = transcriptId

    def set_taxon_range(self, tax):
        if self.taxon is not None and tax != self.taxon:
            raise EvolutionaryConceptError("Gene can only belong to one genome")
        self.taxon = tax

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, self.unique_id)


class EvolutionaryConceptError(Exception):
    pass