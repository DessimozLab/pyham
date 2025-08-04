import logging
import warnings
import gzip
from os import PathLike
from typing import Union, Optional
from xml.etree.ElementTree import XMLParser
from ete3 import Tree

from .genome import ExtantGenome, AncestralGenome, Genome

logger = logging.getLogger(__name__)


def create_tree_from_newick(tree: Union[str, PathLike], tree_format: Optional[str] = None, use_internal_name: bool = False, quoted_node_names: bool = True):
    if tree_format is None:
        tree_format = 'newick_string' if isinstance(tree, str) and tree.strip().startswith('(') else 'newick'
    if tree_format == 'newick_string':
        # If tree is a string and not a path, treat it as a newick string
        tree = tree.strip()
        if not tree.endswith(';'):
            tree += ';'
        tree = Tree(tree, format=1, quoted_node_names=quoted_node_names)
    elif tree_format == 'newick':
        # If tree is a path, read the newick file
        tree = Tree(str(tree), format=1, quoted_node_names=quoted_node_names)
    else:
        raise ValueError("tree_format must be 'newick' or 'newick_string'")

    if not use_internal_name:
        # Generate internal node names by concatenating children names
        for node in tree.traverse("postorder"):
            if not node.is_leaf():
                node.name = '/'.join(child.name for child in node.get_children())

    return tree


def create_tree_from_phyloxml(tree: Union[str, PathLike], phyloxml_leaf_name_tag='clade_name', phyloxml_internal_name_tag='clade_name', use_internal_name=True):
    from .parsers import PhyloXMLToETE
    factory = PhyloXMLToETE()
    parser = XMLParser(target=factory)

    open_ = gzip.open if str(tree).endswith('.gz') else open
    with open_(tree, 'rt') as fh:
        for line in fh:
            parser.feed(line)

    def set_leaf_name(node):
        attr = phyloxml_leaf_name_tag.split('_', 1)[-1]
        if not hasattr(node, attr) or len(getattr(node, attr)) == 0:
            raise KeyError(f"Node {node} in the phyloxml file {tree} has no {phyloxml_leaf_name_tag} attribute to populate the species name")
        node.name = getattr(node, attr)

    def set_internal_name(node):
        if not use_internal_name:
            # Generate internal node names by concatenating children names
            node.name = '/'.join(child.name for child in node.get_children())
        else:
            attr = phyloxml_internal_name_tag.split('_', 1)[-1]
            if not hasattr(node, attr) or len(getattr(node, attr)) == 0:
                raise KeyError(f"Node {node} in the phyloxml file {tree} has no {phyloxml_internal_name_tag} attribute to populate the species name")
            node.name = getattr(node, attr)

    for node in factory.root.traverse("postorder"):
        # assign name to extant species
        if node.is_leaf():
            set_leaf_name(node)
        else:
            set_internal_name(node)

    return factory.root


class Taxonomy(object):
    """
    Taxonomy is a class to wrap the ete3 Etree used as reference species tree by Ham.
    
    Attributes:
        | tree_str (:obj:`str`): newick tree string used to build the ete3 Etree object.
        | tree (:obj:`ete3 Etree`): species ete3 Etree tree.
        | internal_nodes (:obj:`set`): Set of Etree node that contained a AncestralGenome.
        | leaves (:obj:`set`): Set of Etree node that contained a ExtantGenome.

    """
    def __init__(self, tree, **kwargs):
        """
        Args:
            | tree_file (:obj:`str`): Path to the file that contained the taxonomy information.
            | tree_format (:obj:`str`): type of inputted tree file. Defaults to newick_string. Can be 'newick', 'phyloxml, 'newick_string'.
            | use_internal_name (:obj:`Boolean`, optional): Specify wheter using the given internal node name or use the 
            | concatenatation of the children name. Defaults to False.
            | quoted_node_names (:obj:'Boolean', optional): Specify whether newick file has quoted node names.
        """
        if not isinstance(tree, Tree):
            # old style API
            warnings.warn("Taxonomy should be constructed from Taxonomy.from_newick or Taxonomy.from_phyloxml", DeprecationWarning)
            tree_format = kwargs.pop('tree_format', 'newick_string')
            if tree_format in ('newick', 'newick_string'):
                quoted_node_names = kwargs.pop('quoted_node_names', True)
                use_internal_name = kwargs.get('use_internal_name', False)
                self.tree = create_tree_from_newick(
                    tree, tree_format=tree_format, use_internal_name=use_internal_name, quoted_node_names=quoted_node_names
                )
            elif tree_format == 'phyloxml':
                phyloxml_leaf_name_tag = kwargs.pop('phyloxml_leaf_name_tag', 'clade_name')
                phyloxml_internal_name_tag = kwargs.pop('phyloxml_internal_name_tag', 'clade_name')
                use_internal_name = kwargs.get('use_internal_name', True)
                self.tree = create_tree_from_phyloxml(
                    tree,
                    phyloxml_leaf_name_tag=phyloxml_leaf_name_tag,
                    phyloxml_internal_name_tag=phyloxml_internal_name_tag,
                    use_internal_name=use_internal_name
                )
        else:
            # new style API: tree is already an Tree object
            self.tree = tree

        # check unicity of leaves name.
        self._check_consistency_names()

        # add depth to each node of the tree.
        self._add_depth(self.tree.get_tree_root(), depth=0)

        # tracker for Genome created.
        self.internal_nodes = set()
        self.leaves = set()

    @classmethod
    def from_newick(cls, tree: Union[str, PathLike], use_internal_name=False, quoted_node_names=True):
        """  Create a Taxonomy object from a newick file or string.

        Args:
            | tree (:obj:`str`): Path to the file that contained the taxonomy information, or a newick string.
            | use_internal_name (:obj:`Boolean`, optional): Specify whether using the given internal node name or use the
            | concatenatation of the children name. Defaults to False.
            | quoted_node_names (:obj:'Boolean', optional): Specify whether the newick tree has quoted node names. Defaults to True.

        Returns:
            obj:`Taxonomy`: Taxonomy object with the ete3 tree.
        """
        tree_inst = create_tree_from_newick(tree, quoted_node_names=quoted_node_names, use_internal_name=use_internal_name)
        return cls(tree_inst)

    @classmethod
    def from_phyloxml(cls, tree: Union[str, PathLike], phyloxml_leaf_name_tag='clade_name', phyloxml_internal_name_tag='clade_name'):
        tree_inst = create_tree_from_phyloxml(tree, phyloxml_leaf_name_tag=phyloxml_leaf_name_tag, phyloxml_internal_name_tag=phyloxml_internal_name_tag)
        return cls(tree_inst)

    @classmethod
    def from_ete_tree(cls, tree: Tree):
        return cls(tree)

    @property
    def tree_str(self):
        return self.tree.write(format=8, format_root_node=True)

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

    def is_child_recursive(self, node, ancestor_node):
        """  Check if the node is a child of the ancestor node.

            Args:
                | node (:obj:`node`): node to check.
                | ancestor_node (:obj:`node`): ancestor node.

            Returns:
                :obj:`bool` True if the node is a child of the ancestor node.
        """
        for tax in node.iter_ancestors():
            if tax == ancestor_node:
                return True
        return False

    def get_newick_from_tree(self, node):
        """  return the newick tree string (format 8: all names) rooted at the given node.

             Args:
                 node (:obj:`node`): root node.

             Returns:
                 :obj:`str` of the subtree.
         """

        return node.write(format=8, format_root_node=True)

    def _check_consistency_names(self):

        """  
        Check if leaves names and internal node names are uniques.
        
        """

        leaf_names = []
        int_names = []

        for node in self.tree.traverse():
            if node.name is None:
                raise KeyError("{} node has no name attribute".format(node))
            if node.is_leaf():
                leaf_names.append(node.name)
            else:
                int_names.append(node.name)

        # check for leaves name
        if len(leaf_names) != len(set(leaf_names)):
            raise KeyError("Leaves names are not unique ! Leaves founded: {}".format(int_names))

        # Check for internal names
        if len(int_names) != len(set(int_names)):
            raise KeyError("Internal Names are not unique. Internal names founded: {}. If you specify use_internal_name=False, please report the bug to us.".format(int_names))

    def _add_depth(self, node, depth=0):
        """  
        Recursive function to add depth to each node of a Etree.
        """
        node.add_feature("depth", depth)
        for n in node.get_children():
            self._add_depth(n, depth + 1)




def build_taxon_node(id, name=None, **kwargs):
    node = Tree(name=name)
    node.add_feature("taxon_id", id)
    return node


