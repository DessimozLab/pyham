from xml.etree.ElementTree import XMLParser
import ete3
from . import taxonomy as tax
from . import genome
from ham import parsers
import logging
import sys
from . import abstractgene

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

    def get_all_top_level_hogs(self):
        return self.toplevel_hogs

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

    def _get_mrca_ancestral_genome_using_hog_children(self, hog):

            children_genomes = set()
            children_nodes = set()

            for child in hog.children:
                children_genomes.add(child.genome)

            for e in children_genomes:
                children_nodes.add(e.taxon)

            source = children_nodes.pop()
            common = source.get_common_ancestor(children_nodes)

            return self.get_ancestral_genome(common)

    def _add_missing_taxon(self, source_hog , target_hog, missing_taxons):
        """

        :param source_hog: hog
        :param target_hog: hog
        :param missing_taxons: array
        :return:
        """

        target_hog.remove_child(source_hog)

        current_child = source_hog

        for tax in missing_taxons:
            ancestral_genome = self.get_ancestral_genome(tax)
            hog = abstractgene.HOG()
            hog.set_genome(ancestral_genome)
            hog.add_child(current_child)
            current_child = hog

        target_hog.add_child(current_child)
