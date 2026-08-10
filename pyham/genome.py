import abc
from ete4 import Tree


class Genome(metaclass=abc.ABCMeta):
    """  
    Genome is an abstract class representing extant or ancestral genomes. An Genome is defined by a unique taxon (node)
    in the Taxonomy.tree.

    Attributes:
        | taxon (:obj:`ete3.TreeNode` of the :obj:`pyham.taxonomy`.tree): corresponding taxon.
        | name (:obj:`str`): Name of the Genome. Get from the newick tree if specified otherwise build it by concatenating all children genome names.
        | genes (:obj:`list`): list of :obj:`pyham.abstractgene.AbstractGene` related to this Genome.

    """

    def __init__(self):
        self.genes = []
        self.taxon = None
        self.name = None

    def add_gene(self, gene):

        """  
        This method adds an AbstractGene to the genes attributes and updates the AbstractGene.genome attribute.

        Attributes:
            | taxon (Taxonomy.tree node): corresponding taxon.
            | name (:obj:`str`): Name of the Genome. Get from the newick tree if specified otherwise build it by concatenating all children genome names.
            | genes (:obj:`list`): list of :obj:`AbstractGene` related to this Genome.
        
        Raises:
            TypeError: if gene is not an AbstractGene.
        """
        from .abstractgene import AbstractGene
        if not isinstance(gene, AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(AbstractGene.__name__,
                                    type(gene).__name__))

        self.genes.append(gene)
        gene.set_genome(self)

    def set_taxon(self, taxon):

        """  
        This method sets the taxon of the Genome.

        Attributes:
            taxon (Taxonomy.tree node): corresponding taxon.

        Raises:
            TypeError: if taxon is not an ete4.Tree.
            EvolutionaryConceptError: if the genome already has a taxon set.
        """

        if not isinstance(taxon, Tree):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Tree.__name__,
                                    type(taxon).__name__))
        if self.taxon is not None and taxon != self.taxon:
            raise EvolutionaryConceptError("only one taxon can refer to one genome")
        self.taxon = taxon

    @abc.abstractmethod
    def get_number_genes(self):
        """ Return the number of genes"""
        pass

    def __str__(self):
        return self.name


class AncestralGenome(Genome):

    """  
    AncestralGenome class representing ancestral genomes.

    Attributes:
        ancestral_clustering (dict): a dictionary that map each of this ancestral genome HOGs to its list of descendant
        extant genes.

    """

    def __init__(self):
        super(AncestralGenome, self).__init__()
        self.ancestral_clustering = None

    def get_ancestral_clustering(self):

        """ Lazy getter of the ancestral_clustering attribute.

            Returns:
                :obj:`dict`.

        """

        if self.ancestral_clustering is None:
            self.ancestral_clustering = {}
            for hog in self.genes:
                self.ancestral_clustering[hog] = hog.get_all_descendant_genes()
        return self.ancestral_clustering

    def get_number_genes(self):
        return len(self.genes)


class ExtantGenome(Genome):

    """  
    ExtantGenome class representing extant genomes. An ExtantGenome is defined by a unique name and a NCBI taxId.

    Attributes:
        | name (:obj:`str`): unique genome name. This is the unique identifier of ExtantGenome.
        | taxid (dict): NCBI taxId.

    """

    def __init__(self, name, NCBITaxId='-1', **kwargs):
        super(ExtantGenome, self).__init__()
        self.name = name
        self.taxid = NCBITaxId
        self.taxon_id = kwargs.get('taxon_id', None)

    def get_number_genes(self, singleton=True): #  TODO: UT

        """ Get the number of genes of this ExtantGenome.
        
            Args:
                | singleton (:obj:`Bool`): boolean to take into account singleton or not.

            Returns:
                :obj:`int`.

        """

        if singleton:
            return len(self.genes)
        else:
            nbr = 0
            for g in self.genes:
                if not g.is_singleton():
                    nbr += 1
            return nbr


class EvolutionaryConceptError(Exception):
    pass
