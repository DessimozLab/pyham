from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import super
from builtins import object
from future import standard_library
standard_library.install_aliases()
import numbers
from .genome import ExtantGenome, AncestralGenome, Genome
from .iham import IHAM
import abc, six

@six.add_metaclass(abc.ABCMeta)
class AbstractGene(object):
    """  
    AbstractGene is an abstract class representing extant or ancestral genes. An AbstractGene is defined by an unique
    id, the genome it belongs to and  its parent AbstractGene.
    
    Attributes:
        | parent (:obj:`pyham.abstractgene.HOG`): Direct parent HOG this gene is descendant to. For and only for top level HOG, parent
        attribute is set to  None.
        | genome (:obj:`pyham.genome.Genome`): Related Genome object.
        arose_by_duplication: pyham.abstractgene.DuplicationNode if this AbstractGene arose by a duplication from its parent otherwise False.
        
    """

    def __init__(self, arose_by_duplication=False, **kwargs):

        """
        Attributes:
            | arose_by_duplication (:obj:`Bool`): True if this AbstractGene arose by a duplication from its parent.
            | ** kwargs: dictionary of attribute and value required to create the: obj:`pyham.abstractgene.AbstractGene`.
        """

        self.parent = None
        self.genome = None
        self._properties = {}
        self.arose_by_duplication = arose_by_duplication

    @abc.abstractmethod
    def set_genome(self, genome):
        """Set the genome attribute using the given :obj:`Genome`."""
        pass

    @abc.abstractmethod
    def is_singleton(self):
        """ Return true if singletons otherwise false"""
        pass

    def search_ancestor_hog_in_ancestral_genome(self, ancestral_genome):

        """ 
        Get the :obj:`HOG` related to the query :obj:`AncestralGenome` from which this :obj:`AbstractGene` is
        descendant otherwise return None. In addition, it returns a boolean whether a duplication occurs in between the two 
        :obj:`AbstractGene`.
        
            Args:
                ancestral_genome (:obj:`pyham.genome.AncestralGenome`): Ancestral genome of interest.

            Returns:
                return :obj:`pyham.abstractgene.HOG` or None and :obj:`bool`

        """

        found = None
        current_hog = self
        paralog = current_hog.arose_by_duplication

        # Until we reach the top level HOG or the ancestor HOG is founded
        while current_hog.parent is not None and found is None:
            current_hog = current_hog.parent
            if current_hog.genome == ancestral_genome:
                return current_hog, paralog
            if current_hog.arose_by_duplication != False:
                paralog = True
        return found, paralog

    def get_top_level_hog(self):

        """ 
        Get the related top level :obj:`pyham.abstractgene.HOG`. If :obj:`pyham.abstractgene.HOG` return self, the :obj:`pyham.abstractgene.HOG` is already top level. If :obj:`pyham.abstractgene.Gene` return self, the :obj:`pyham.abstractgene.Gene` is a singleton.

            Returns:
                return :obj:`pyham.abstractgene.HOG`.
        """

        current_hog = self

        while current_hog.parent is not None:
            current_hog = current_hog.parent

        return current_hog

    def get_at_level(self, genome):

        """ 
            Get for the query HOG, the ancestor or children AbstractGene in the genome of interest
            
            Args:
                genome (:obj:`pyham.genome.Genome`): Genome of interest.
                
            Returns:
                list of  :obj:`pyham.abstractgene.AbstractGene`.
                
            Raises:
                TypeError: if genome is not a Genome.
                KeyError: If no HOG matches the query genome or it return itself. 
        """


        if not isinstance(genome, Genome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Genome.__name__,
                                    type(genome).__name__))

        tl = self.get_top_level_hog()

        def pre(node, elem):
            if node.genome == genome:
                elem.append(node)
            return elem

        def leaf(node, leaf, elem):
            if leaf.genome == genome:
                elem.append(leaf)
            return elem

        rlist = tl.visit([], function_extant_gene=leaf, function_prefix=pre)

        if not rlist:
            raise KeyError("Level {} not found within this HOG".format(genome))

        for e in rlist:
            if e == self:
                raise KeyError("get level cannot return itself".format())

        return rlist

    def add_property(self, name, value):
        self._properties[name] = value

    def __getitem__(self, item):
        return self._properties[item]


class HOG(AbstractGene):

    """  
    HOG class represents ancestral genes. An AbstractGene is defined the genome it belongs to and its child/parent
    AbstractGene. Only top level HOG doesn't have parent attribute (set to None).

    Attributes:
        | hog_id (:obj:`str`): hog id. Defaults is None.
        | children (:obj:`list` of :obj:`pyham.abstractgene.AbstractGene`): A list of direct descendants AbstractGene.
        | hogvis (:obj:`pyham.IHAM`): :obj:`pyham.IHAM` object of this HOG.
        | duplications (:obj:`list` of :obj:`pyham.abstractgene.DuplicationNode`): list of all duplication node child of this HOG.

    """

    def __init__(self, id=None, **kwargs):
        super(HOG, self).__init__(**kwargs)
        self.hog_id = id
        self.children = []
        self.hogvis = None
        self.duplications = []

    def add_child(self, child_to_add):

        """  
            Add in the children attribute a given :obj:`pyham.abstractgene.AbstractGene` and update parent attribute in child AbstractGene.

            Attributes:
                child_to_add (:obj:`pyham.abstractgene.AbstractGene`): add_child to add.

        """

        if not isinstance(child_to_add, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(child_to_add).__name__))

        if self == child_to_add:
            raise EvolutionaryConceptError('Cannot be a child of itself')

        self.children.append(child_to_add)
        child_to_add.parent = self

    def remove_child(self, child_to_remove):

        """  
            Remove from the children attribute a given :obj:`AbstractGene` and update parent attribute in child
            AbstractGene.

            Attributes:
                child_to_add (:obj:`AbstractGene`): add_child to add.

        """

        if not isinstance(child_to_remove, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(child_to_remove).__name__))

        if child_to_remove in self.children:
            self.children.remove(child_to_remove)
            child_to_remove.parent = None

        else:
            raise ValueError("element not found in the hog children")

    def score(self, score_id, value=None):
        """  
        Get or set a score attribute for this HOG.

        If no value is passed, the stored value is returned or a KeyError is raised if this score has not been stored.

        Attributes:
            score_id (:obj:`ExtantGenome`): score to be set or returned.
            value (optional): A numeric value for the score of this HOG.
            score_id (:obj:`ExtantGenome`): score to be set or returned.

        Raises:
            KeyError: if score_id is unknown and accessed.
        """

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

        """  
        This method set :obj:`pyham.genome.AncestralGenome` given as parameter as genome attribute.

        Attributes:
            | genome (:obj:`pyham.genome.AncestralGenome`): Ancestral genome to be set as genome.
        """

        # Check that genome is an :obj:`AncestralGenome`
        if not isinstance(genome, AncestralGenome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AncestralGenome.__name__,
                                    type(genome).__name__))

        if self.genome is not None and genome != self.genome:
            raise EvolutionaryConceptError("HOG can only be mapped to one ancestral genome")

        self.genome = genome

    def visit(self, elem, function_extant_gene=None, function_postfix=None, function_prefix=None):

        """  
        Recursive function to traverse the HOG nested children hierarchy. 
        
        An arbitrary object "elem" is return at each visit() and is pass/returned from each recursive state function 
        calls. 
        
        There is tree recursive state functions:
            - function_extant_gene: Function called at each ExtantGene child. Its take as argument the child ExtantGene
             and the "elem" and return processed "elem".
            - function_postfix: Function called after each HOG child visit(). Its take as argument the child HOG and 
            the "elem" and return processed "elem".
            - function_prefix: First function called that take as argument the current HOG and the "elem" and return
             processed "elem".

        Attributes:
            | elem (arbitrary object): object that is pass through visit recursive and recursive state functions calls.
            | function_extant_gene (callback function): Callback function for ExtantGenes (leaves).
            | function_postfix (callback function): Callback function for after child visit current HOG.
            | function_prefix (callback function): Callback function for current HOG.
        
        Returns:
            return arbitrary
        """

        if function_prefix is not None:
            elem = function_prefix(self, elem)

        for child in self.children:
            if isinstance(child, Gene):
                if function_extant_gene is not None:
                    elem = function_extant_gene(self, child, elem)
            else:
                child.visit(elem, function_extant_gene=function_extant_gene, function_postfix=function_postfix,
                            function_prefix=function_prefix)
                if function_postfix is not None:
                    elem = function_postfix(self, child, elem)
        return elem

    #  TODO: UT

    def get_all_descendant_genes(self):

        """ 
        Get all :obj:`Gene` present in this HOG.

            Returns:
                list of :obj:`pyham.abstractgene.Gene`

        """

        def append_leaf(current, child, list):
            list.append(child)
            return list

        return self.visit([], function_extant_gene=append_leaf)

    def get_all_descendant_genes_clustered_by_species(self):
        """ 
        Get all :obj:`pyham.abstractgene.Gene` present in this HOG clustered by species :obj:`pyham.genomeExtantGenome`.

            Returns:
                Dictionary of :obj:`pyham.genome.ExtantGenome` map their list of :obj:`pyham.abstractgene.Gene`.

        """
  
        def append_add_lead_with_genome(current, child, dict):
            dict.setdefault(child.genome, []).append(child)
            return dict

        return self.visit({}, function_extant_gene=append_add_lead_with_genome)

    def get_all_descendant_hogs(self): # TODO: self is also returned with it
        
        """ 
        Get all :obj:`pyham.abstractgene.HOG` present in this HOG hierarchy.

            Returns:
                Dictionary of :obj:`pyham.genomeExtantGenome` map their list of :obj:`pyham.abstractgene.Gene`.

        """

        def append_child(current, list):
            list.append(current)
            return list
        
        return self.visit([], function_prefix=append_child)

    def get_all_descendant_hog_levels(self): # TODO: self is also returned with it

        """ 
        Get all :obj:`pyham.genome.Genome` present in this HOG hierarchy.

            Returns:
                List of :obj:`pyham.genome.Genome`.

        """

        def append_current_genome(current,list):
            list.append(current.genome)
            return list

        return self.visit([], function_prefix=append_current_genome)

    def get_hog_vis(self, newick_str):

        """ Lazy getter of the :obj:`pyham.abstractgene.HOG` IHAM.

            Args:
                newick_str (:obj:`str`): newick species tree used by the create_iHam.

            Returns:
                :obj:`pyham.IHAM` of this HOG.

        """

        if self.hogvis is None:
            self.hogvis = IHAM(newick_str, self)

        return self.hogvis

    #  TODO: end UT

    def is_singleton(self):
        return False

    def __repr__(self):
        return "<{}({})>".format(self.__class__.__name__, self.hog_id if self.hog_id else "")


class Gene(AbstractGene):
    """Gene class for extant genes (inherit from AbstractGene).

    Attributes:
        | unique_id  internal unique id
        | gene_id    id used to mapped to external ids.
        | prot_id    id used to mapped to external ids.
        | transcript_id    id used to mapped to external ids.

    """

    def __init__(self, id, geneId=None, protId=None, transcriptId=None, **kwargs):
        super(Gene, self).__init__(**kwargs)
        self.unique_id = id
        self.gene_id = geneId
        self.prot_id = protId
        self.transcript_id = transcriptId

    def set_genome(self, genome):
        """  
        This method set :obj:`pyham.genome.ExtantGenome` given as parameter as genome attribute.

        Attributes:
            genome (:obj:`pyham.genome.ExtantGenome`): Extant genome to be set as genome.
             
        Raises:
            TypeError: If genome is not :obj:`pyham.genome.ExtantGenome`.
            EvolutionaryConceptError: If we try to add more than one genome.
         """

        # Check that genome is an :obj:`ExtantGenome`
        if not isinstance(genome, ExtantGenome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(ExtantGenome.__name__,
                                    type(genome).__name__))

        if self.genome is not None and genome != self.genome:
            raise EvolutionaryConceptError("Gene can only belong to one genome")

        self.genome = genome

    def get_dict_xref(self): #              <-- TODO: UT

        """ Get the dictionary of cross references.

            Returns:
                Dictionary of external id tag map with their value.
        """

        xref = {}

        if self.unique_id is not None:
            xref["id"] = self.unique_id
        if self.gene_id is not None:
            xref["geneId"] = self.gene_id
        if self.prot_id is not None:
            xref["protId"] = self.prot_id
        if self.transcript_id is not None:
            xref["transcriptId"] = self.transcript_id

        return xref

    def is_singleton(self):
        if self.parent is None:
            return True
        else:
            return False

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, self.unique_id)


class EvolutionaryConceptError(Exception):
    pass


class DuplicationNode(object):
    """
        This object link together all abstract genes that emerges from the same duplication event. Its composed of a set of genes
        all descendant from this duplication event and their direct related parent hog. In addition it contained the parent of the MRCA of all 
        children genes (otherwise MRCA return genome father duplication).
        
        Attributes:
        | parent (:obj:`pyham.abstractgene.HOG`): Direct parent HOG this DuplicationNode.
        | children: List of all duplicated genes
        | MRCA: Parent of the MRCA of all children gene that contained the parent HOG.
            
    """

    def __init__(self,ham_object):
        self.ham = ham_object
        self.MRCA = None
        self.children = []
        self.parent = None

    def set_parent(self, hog):
        """
                Args:
                    | hog (:obj:`pyham.abstractegene.HOG`): parent :obj:`pyham.abstractegene.HOG` of this duplication.

        """

        if not isinstance(hog, HOG):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(HOG.__name__,
                                    type(hog).__name__))

        hog.duplications.append(self)
        self.parent = hog

    def set_MRCA(self):
        """
            Compute the MRCA of all genes genomes.
        """

        children_genomes = set([child.genome for child in self.children])

        if len(children_genomes) < 2:
            self.MRCA = self.ham._get_ancestral_genome_by_taxon(list(children_genomes)[0].taxon.up)

        else:
            self.MRCA = self.ham._get_ancestral_genome_by_mrca_of_genome_set(children_genomes)
            self.MRCA = self.ham._get_ancestral_genome_by_taxon(self.MRCA.taxon.up)

    def add_child(self, child):
        """
            Add child to list of children
        """

        if not isinstance(child, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(child).__name__))

        self.children.append(child)
        child.arose_by_duplication = self

    def remove_child(self, child_to_remove):
        """
            remove child to list of children
        """

        if not isinstance(child_to_remove, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(child_to_remove).__name__))

        if child_to_remove in self.children:
            self.children.remove(child_to_remove)
            child_to_remove.arose_by_duplication = False

        else:
            raise ValueError("element not found in the duplication children")


