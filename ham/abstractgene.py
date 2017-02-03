__author__ = 'admin'

from abc import ABCMeta, abstractmethod
import numbers
from ham.genome import ExtantGenome, AncestralGenome

class AbstractGene(metaclass=ABCMeta):
    """AbstractGene class for gene entities.

    Attributes:
        parent      direct parent AbstractGene (can only be a HOG) of this AbstractGene
        genome      Genome containing this AbstractGene.

    """
    def __init__(self):
        self.parent = None
        self.genome = None

    @abstractmethod
    def set_genome(self, genome):
        """Set the genome attribute using the """
        pass


class HOG(AbstractGene):
    """HOG class for HOGs (inherit from AbstractGene).

    Attributes:
        hog_id      internal unique id
        children    AbstractGene objects that have descended from this HOG.

    """

    def __init__(self, id=None, is_paralog=False, **kwargs):
        super(HOG, self).__init__()
        self.hog_id = id
        self.children = []
        self.is_paralog = is_paralog


    # add a AbstractGene to the hog children
    def add_child(self, elem):
        if not isinstance(elem, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(elem).__name__))
        if self == elem:
            raise EvolutionaryConceptError('Cannot be a child of itself')
        self.children.append(elem)
        elem.parent = self

    # safely remove a child from the HOG children and its parent reference to the HOG
    def remove_child(self, elem):
        if not isinstance(elem, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(elem).__name__))
        if elem in self.children:
            self.children.remove(elem)
            elem.parent = None
        else:
            raise ValueError("element not found in the hog children")

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

    def set_genome(self, genome):
        if not isinstance(genome, AncestralGenome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AncestralGenome.__name__,
                                    type(genome).__name__))
        else:
            if self.genome is not None and genome != self.genome:
                raise EvolutionaryConceptError("HOG can only be mapped to one ancestral genome")
            self.genome = genome


    def visit(self, elem, function_leaf=None, function_postfix=None, function_prefix=None):

        if function_prefix != None:
            elem = function_prefix(self, elem)

        for child in self.children:
            if isinstance(child, Gene):
                if function_leaf != None:
                    elem = function_leaf(self, child, elem)
            else:
                child.visit(elem, function_leaf=function_leaf, function_postfix=function_postfix, function_prefix=function_prefix)
                if function_postfix != None:
                    elem = function_postfix(self, child, elem)
        return elem

    def __repr__(self):
        return "<{}({})>".format(self.__class__.__name__, self.hog_id if self.hog_id else "")


class Gene(AbstractGene):
    """Gene class for extant genes (inherit from AbstractGene).

    Attributes:
        unique_id  internal unique id
        gene_id    id used to mapped to external ids.
        prot_id    id used to mapped to external ids.
        transcript_id    id used to mapped to external ids.

    """
    def __init__(self, id, geneId=None, protId=None, transcriptId=None, **kwargs):
        super(Gene, self).__init__()
        self.unique_id = id
        self.gene_id = geneId
        self.prot_id = protId
        self.transcript_id = transcriptId

    def set_genome(self, genome):
        if not isinstance(genome, ExtantGenome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(ExtantGenome.__name__,
                                    type(genome).__name__))
        else:
            if self.genome is not None and genome != self.genome:
                raise EvolutionaryConceptError("Gene can only belong to one genome")
            self.genome = genome

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, self.unique_id)


class EvolutionaryConceptError(Exception):
    pass
