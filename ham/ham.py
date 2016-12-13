import ete3
from . import taxonomy as tax

def build_taxonomy_and_ancestral_genomes(newick_str):
    '''
    Take newick tree string as reference to build Taxonomy object and AncestralGenome objects
    :param newick_str: string containing a newick tree
    :return : a Taxonomy object and its related AncestralGenome objects
    '''

    taxonomy = tax.Taxonomy()
    taxonomy.newick_str = newick_str
    taxonomy.tree = ete3.Tree(taxonomy.newick_str, format=1)


    for node in taxonomy.tree.traverse("postorder"):
        # Do some analysis on node
        if node.is_leaf():
            print("\t",node.name)
        else:
            print(node.name)

    ancestral_genomes = None

    return taxonomy, ancestral_genomes


def build_hogs_and_genes(file_object):

    hogs = None
    genes = None

    return hogs, genes


def resolve_taxonomy_and_hogs(taxonomy, hogs, genes):
    pass