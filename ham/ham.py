from xml.etree.ElementTree import XMLParser
from . import taxonomy as tax
from . import genome
from ham import parsers
from ham import mapper
import logging
import sys
from . import abstractgene
import copy
logger = logging.getLogger(__name__)


def _build_hogs_and_genes(file_object, taxonomy=None):
    """
    build AbstractGene.HOG and AbstractGene.Gene using a given orthoXML file object
    :param file_object: orthoXML file object
    :param taxonomy: Taxonomy object used to map taxon with ete3 node
    :return: a set of AbstractGene.HOG and a set of  AbstractGene.Gene
    """

    factory = parsers.OrthoXMLParser(taxonomy)
    parser = XMLParser(target=factory)

    for line in file_object:
        parser.feed(line)

    return factory.toplevel_hogs, factory.extant_gene_map


class HAM(object):
    def __init__(self, newick_str=None, hog_file=None, type="orthoxml"):

        if newick_str is None or hog_file is None:
            logger.debug('newick_str: {}, hogs_file: {}'.format(newick_str, hog_file))
            sys.exit('Both newick string or hogs file are required to create HAM object')

        self.newick_str = newick_str
        self.hog_file = hog_file
        self.hog_file_type = type
        self.toplevel_hogs = None
        self.extant_gene_map = None
        self.HOGMaps = {} # In order to not recompute two time a hog map between two level we are storing them globally

        self.taxonomy = tax.Taxonomy(self.newick_str)
        logger.info('Build taxonomy: completed.'.format(self.newick_str))

        if type == "orthoxml":
            with open(self.hog_file, 'r') as orthoxml_file:
                self.toplevel_hogs, self.extant_gene_map = _build_hogs_and_genes(orthoxml_file, self)
            logger.info('Parse Orthoxml: {} top level hogs and {} extant genes extract.'.format(len(self.toplevel_hogs),
                                                                                                len(
                                                                                                    self.extant_gene_map)))

        elif type == "hdf5":
            # Looping through all orthoXML within the hdf5
            #   for each run self.build_...
            #       update self.toplevel_hogs and self.extant_gene_map for each
            pass

        logger.info(
            'Set up HAM analysis: ready to go with {} hogs founded within {} species.'.format(len(self.toplevel_hogs),
                                                                                              len(
                                                                                                  self.taxonomy.leaves)))

    def get_HOGMap(self, genome_pair_set):
        """
        if not already computed, do the hog mapping between a pair of levels.
        :param genome_pair_set:
        :return:
        """
        f = frozenset(genome_pair_set)
        if f in self.HOGMaps.keys():
            return self.HOGMaps[f]
        else:
            m = mapper.HOGsMap(self, genome_pair_set)
            self.HOGMaps[f] = m
            return m

    def get_all_top_level_hogs(self):
        return self.toplevel_hogs

    def get_all_genes_of_hog(self, hog):

        def append_child(current, child, list):
            list.append(child)
            return list

        return hog.visit([], function_leaf=append_child)

    def get_all_extant_genes_dict(self):
        return self.extant_gene_map

    def get_extant_genomes(self):
        """
        return all Genome.ExtantGenome
        :return: set of Genome.ExtantGenome
        """
        return set(leaf.genome for leaf in self.taxonomy.leaves)

    def get_ancestral_genomes(self):
        """
        return all Genome.AncestralGenome
        :return: set of Genome.AncestralGenome
        """
        return set(internal_node.genome for internal_node in self.taxonomy.internal_nodes)

    def get_ancestral_genome(self, tax_node):

        if "genome" in tax_node.features:
                return tax_node.genome

        else:
            ancestral_genome = genome.AncestralGenome()
            self.taxonomy.add_ancestral_genome_to_node(tax_node, ancestral_genome)

            return ancestral_genome

    def _get_extant_genome_by_name(self,**kwargs):

        nodes_founded = self.taxonomy.tree.search_nodes(name=kwargs['name'])

        if len(nodes_founded) == 1:
            if "genome" in nodes_founded[0].features:
                return nodes_founded[0].genome

            else:
                extant_genome = genome.ExtantGenome(**kwargs)
                self.taxonomy.add_extant_genome_to_node(nodes_founded[0], extant_genome)
                return extant_genome
        else:
            raise ValueError('{} node(s) founded for the species name: {}'.format(len(nodes_founded), kwargs['name']))

    def _get_mrca_ancestral_genome_from_genome_set(self, genome_set):
        """
        return the MRCA of all the genomes present in genome_set.
        :param genome_set: set of ancestral genomes.
        :return: ancestral genome
        """
        if len(genome_set) < 2:
            raise ValueError('Only one genome is not enough to computed MRCA: {}'.format(genome_set))
        pass

        for g in genome_set:
            if not isinstance(g, genome.Genome):
                raise TypeError("expect subclass obj of '{}', got {}"
                                .format(genome.Genome.__name__,
                                        type(g).__name__))

        children_nodes = set()

        for e in genome_set:
            children_nodes.add(e.taxon)

        common = self.taxonomy.tree.get_common_ancestor(children_nodes)

        return self.get_ancestral_genome(common)

    def _get_mrca_ancestral_genome_using_hog_children(self, hog):
        """
        collect all the genomes insides the hog children and them look for their mrca
        :param hog:
        :return:
        """

        children_genomes = set()

        for child in hog.children:
            children_genomes.add(child.genome)

        return self._get_mrca_ancestral_genome_from_genome_set(children_genomes)

    def _add_missing_taxon(self, child_hog , oldest_hog, missing_taxons):
        """

        :param child_hog: youngest hog
        :param oldest_hog: oldest hog
        :param missing_taxons: array of intermediate taxon sorted from most r4cent to oldest
        :return:
        """
        if not isinstance(child_hog, abstractgene.AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(abstractgene.AbstractGene.__name__,
                                    type(child_hog).__name__))

        if not isinstance(oldest_hog, abstractgene.AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(abstractgene.AbstractGene.__name__,
                                    type(oldest_hog).__name__))

        if oldest_hog == child_hog:
            raise TypeError("Cannot add missing level between an HOG and it self.")

        # the youngest hog is removed from the oldest hog
        oldest_hog.remove_child(child_hog)

        # Then for each intermediate level in between the two hog
        current_child = child_hog
        for tax in missing_taxons:
            # we get the related ancestral genome of this level
            ancestral_genome = self.get_ancestral_genome(tax)

            # we create the related hog and add it to the ancestral genome
            hog = abstractgene.HOG()
            hog.set_genome(ancestral_genome)
            ancestral_genome.add_gene(hog)

            if ancestral_genome.taxon is not current_child.genome.taxon.up:
                raise TypeError("HOG taxon {} is different than child parent taxon {}".format(ancestral_genome.taxon, current_child.genome.taxon.up))

            # we add the child
            hog.add_child(current_child)
            current_child = hog

        oldest_hog.add_child(current_child)

    def _get_oldest_from_genome_pair(self, g1, g2):
        """
        get the oldest genomes from a pair of genomes... I know Adrian, I know!
        :param g1:
        :param g2:
        :return:
        """

        mrca1 = g1.get_common_ancestor(g2)
        if g2 == mrca1:
            return g2

        mrca2 = g2.get_common_ancestor(g1)
        if g1 == mrca2:
            return g1

        return None

    def _get_ancestor_and_descendant(self, genome_set):
        """
        Get from a set of genomes the oldest genomes (or their mrca if oldest not in genomes set) and return the list of descendant (present in the genomes set only)
        :param genome_set:
        :return:
        """
        ancestor = self._get_mrca_ancestral_genome_from_genome_set(genome_set)
        genome_set.discard(ancestor)
        #descendants = genome_set.pop()
        return ancestor, genome_set

    def compare_genomes(self, genomes_set, analysis=None):
        """
        function for level comparaison || Work in progress
        :param genomes_set:
        :param analysis: ADRIAN -> This is not smart
        :return:
        """

        if analysis == "vertical":
            if len(genomes_set) != 2:
                raise TypeError("{} genomes given for vertical HOG mapping, only 2 should be given".format(len(genomes_set)))
            vertical_map = mapper.MapVertical(self)
            vertical_map.add_map(self.get_HOGMap(genomes_set))
            return vertical_map

        elif analysis == "lateral":
            lateral_map = mapper.MapLateral(self)
            anc, desc = self._get_ancestor_and_descendant(copy.copy(genomes_set))
            for g in desc:
                hogmap = mapper.HOGsMap(self, {g, anc})
                lateral_map.add_map(hogmap)
            return lateral_map

        else:
            raise TypeError("Invalid type of genomes comparison")

    def analyze_hog(self):
        pass

    def analyze_ancestral_genomes(self):
        pass
