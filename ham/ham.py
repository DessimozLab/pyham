from xml.etree.ElementTree import XMLParser
import ete3
from . import taxonomy as tax
from . import genome
from ham import parsers


def build_taxonomy_and_ancestral_genomes(newick_str):
    '''
    Take newick tree string as reference to build Taxonomy object and AncestralGenome objects
    :param newick_str: string containing a newick tree
    :return : a Taxonomy object
    '''

    taxonomy = tax.Taxonomy()
    taxonomy.newick_str = newick_str
    taxonomy.tree = ete3.Tree(taxonomy.newick_str, format=1)

    for node in taxonomy.tree.traverse("postorder"):
        if node.is_leaf():
            leaf_genome = genome.ExtantGenome()
            leaf_genome.taxon = node
            node.genome = leaf_genome
            taxonomy.leaves.add(node)

        else:
            internal_genome = genome.AncestralGenome()
            internal_genome.taxon = node
            node.genome = internal_genome
            taxonomy.internal_nodes.add(node)

    return taxonomy


def build_hogs_and_genes(file_object):
    '''
    build AbstractGene.HOG and AbstractGene.Gene using a given orthoXML file object
    :param file_object: orthoXML file object
    :return: a set of AbstractGene.HOG and a set of  AbstractGene.Gene
    '''

    factory = parsers.orthoxmlParser()
    parser = XMLParser(target=factory)

    for line in file_object:
        parser.feed(line)

    return factory.hogs, factory.genes


def resolve_taxonomy_and_hogs():
    pass