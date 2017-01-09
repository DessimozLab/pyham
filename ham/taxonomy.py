import ete3
import logging

logger = logging.getLogger(__name__)


def _add_depth(node, depth=0):
    node.add_feature("depth", depth)
    for n in node.get_children():
        _add_depth(n, depth + 1)


class Taxonomy(object):
    def __init__(self, newick_str):
        self.newick_str = newick_str  # the original newick string used to build the Ete3 Tree
        self.tree = ete3.Tree(self.newick_str, format=1)  # Ete3 Tree object
        _add_depth(self.tree.get_tree_root(), depth=0) # Add depth to ete3 node
        self.internal_nodes = set()  # set of internal node within the self.tree
        self.leaves = set()  # set of leaves within the self.tree

    def add_ancestral_genome_to_node(self, node, genome):
        node.add_feature("genome", genome)
        genome.taxon = node
        self.internal_nodes.add(node)

    def add_extant_genome_to_node(self, node, genome):
        node.add_feature("genome", genome)
        genome.taxon = node
        self.leaves.add(node)
