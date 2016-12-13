__author__ = 'admin'


class Taxonomy(object):

    def __init__(self):
        self.newick_str = ''  # the original newick string used to build the Ete3 Tree
        self.tree = None  # Ete3 Tree object
        self.internal_nodes = set()  # set of internal node within the self.tree
        self.leaves = set()  # set of leaves within the self.tree
        self.map_name_taxa_node = {}  #
