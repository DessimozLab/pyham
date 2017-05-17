from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import str
from future import standard_library
standard_library.install_aliases()
import ete3
import logging
from .genome import ExtantGenome, AncestralGenome, Genome


logger = logging.getLogger(__name__)


class Taxonomy(object):
    """
    Taxonomy is a class to wrap the ete3 Etree used as reference species tree by Ham.
    
    Attributes:
        | newick_str (:obj:`str`): newick tree string used to build the ete3 Etree object.
        | tree (:obj:`ete3 Etree`): species ete3 Etree tree.
        | internal_nodes (:obj:`set`): Set of Etree node that contained a AncestralGenome.
        | leaves (:obj:`set`): Set of Etree node that contained a ExtantGenome.

    """
    def __init__(self, newick_str, use_internal_name=False):
        """
        Args:
            | newick_str (:obj:`str`): Newick str used to build ete3 Etree object.
            | use_internal_name (:obj:`Boolean`, optional): Specify wheter using the given internal node name or use the 
            | concatenatation of the children name. Defaults to False.
        """

        self.newick_str = newick_str
        self.tree = ete3.Tree(self.newick_str, format=1)

        # create internal node name if required
        if use_internal_name is False:
            for node in self.tree.traverse("postorder"):
                if node.is_leaf() is False:
                    self.set_taxon_name(node)

        # check unicity of leaves name.
        self._check_consistency_names()

        # add depth to each node of the tree.
        self._add_depth(self.tree.get_tree_root(), depth=0)

        # tracker for Genome created.
        self.internal_nodes = set()
        self.leaves = set()

    def add_genome_to_node(self, node, genome):
        """  add the given genome to the node attribute "genome".

            Args:
                | node (:obj:`node`): receptor node.
                | genome (:obj:`Genome`): :obj:`Genome` to attach.

        """

        node.add_feature("genome", genome)
        genome.set_taxon(node)

        if isinstance(genome, ExtantGenome):
            self.leaves.add(node)

        elif isinstance(genome, AncestralGenome):
            genome.name = node.name

            self.internal_nodes.add(node)
        else:
            raise TypeError("expect class obj of '{}', got {}".format(type(Genome).__name__,type(genome).__name__))

    def get_path_up(self, lowest_node, ancestor_node):
        """  return the internal node in between two nodes sorted by recentness.

            Args:
                | lowest_node (:obj:`node`): Youngest node.
                | ancestor_node (:obj:`node`): Oldest node.
                
            Returns:
                list of node sorted from most recent to oldest.

        """

        intermediate_level = []

        for tax in lowest_node.iter_ancestors():
            if tax == ancestor_node:
                break
            intermediate_level.append(tax)

        return intermediate_level

    def get_newick_from_tree(self, node):
        """  return the newick tree string (format 8: all names) rooted at the given node.

             Args:
                 node (:obj:`node`): root node.

             Returns:
                 :obj:`str` of the subtree.
         """

        return node.write(format=8, format_root_node=True)

    def set_taxon_name(self, node):
        """  set the node name by concatenation of children name.

             Args:
                 node (:obj:`node`): root node.
        """

        level_name = ""
        for leaf in node:
            level_name += str(leaf.name)
            level_name += "/"

        node.name = str(level_name[:-1])

    def _check_consistency_names(self):

        """  
        Check if leaves names and internal node names are uniques.
        
        """

        leaf_names = []
        int_names = []

        for node in self.tree.traverse():
            if node.is_leaf():
                leaf_names.append(node.name)
            else:
                int_names.append(node.name)


        # check for leaves name
        if len(leaf_names) != len(set(leaf_names)):
            raise KeyError("Leaves names are not unique ! Leaves founded: {}".format(int_names))

        # Check for internal names
        if len(set(int_names)) != len(set(int_names)):
            raise KeyError("Internal Names are not unique. Internal names founded: {}. If you specify use_internal_name=False, please report the bug to us.".format(int_names))

    def _add_depth(self, node, depth=0):
        """  
        Recursive function to add depth to each node of a Etree.
        """
        node.add_feature("depth", depth)
        for n in node.get_children():
            self._add_depth(n, depth + 1)