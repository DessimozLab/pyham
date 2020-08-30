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
from six import BytesIO
from ete3 import Phyloxml
import re


logger = logging.getLogger(__name__)


class Taxonomy(object):
    """
    Taxonomy is a class to wrap the ete3 Etree used as reference species tree by Ham.
    
    Attributes:
        | tree_str (:obj:`str`): newick tree string used to build the ete3 Etree object.
        | tree (:obj:`ete3 Etree`): species ete3 Etree tree.
        | internal_nodes (:obj:`set`): Set of Etree node that contained a AncestralGenome.
        | leaves (:obj:`set`): Set of Etree node that contained a ExtantGenome.

    """
    def __init__(self, tree_file, tree_format='newick_string', use_internal_name=False, phyloxml_leaf_name_tag=None, phyloxml_internal_name_tag=None):
        """
        Args:
            | tree_file (:obj:`str`): Path to the file that contained the taxonomy information.
            | tree_format (:obj:`str`): type of inputted tree file. Defaults to newick_string. Can be 'newick', 'phyloxml, 'newick_string'.
            | use_internal_name (:obj:`Boolean`, optional): Specify wheter using the given internal node name or use the 
            | concatenatation of the children name. Defaults to False.
        """

        self.tree_file = tree_file
        self.tree_format = tree_format
        self.use_internal_name = use_internal_name
        self.phyloxml_leaf_name_tag = phyloxml_leaf_name_tag
        self.phyloxml_internal_name_tag = phyloxml_internal_name_tag

        # create tree
        self.tree = self._build_tree(tree_file, tree_format)

        # create internal node name if required
        self._generate_internal_node_name(self.tree)

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

    def _generate_internal_node_name(self, tree):
        if self.use_internal_name is False:
            for node in tree.traverse("postorder"):
                if node.is_leaf() is False:
                    self.set_taxon_name(node)

        if self.tree_format == 'phyloxml':


            #### IMPORTANT ####

            '''
            This code section is an horrible fix due to incompatibility of Phyloxml().export() with python3.
            
            We overwrite the export module (and its functions) here to deal with both byte and unicode.
            
            TODO: this should be replace as ASAP
            '''

            def showIndent(outfile, level):
                for idx in range(level):
                    g = '    '
                    g = g.encode('UTF-8')
                    outfile.write(g)
            namespace_ = 'phy:'
            name_ = 'Phyloxml'
            namespacedef_ = ''
            outfile = BytesIO()
            level = 0
            def export(xxxx, outfile, level, namespace_='phy:', name_='Phyloxml', namespacedef_=''):
                showIndent(outfile, level)
                x = '<%s%s%s' % (namespace_, name_, namespacedef_ and ' ' + namespacedef_ or '',)
                x = x.encode('UTF-8')
                outfile.write(x)
                already_processed = []
                exportAttributes(xxxx, outfile, level, already_processed, namespace_, name_='Phyloxml')
                if hasContent_(xxxx):
                    y = '>\n'
                    y = y.encode('UTF-8')
                    outfile.write(y)
                    exportChildren(xxxx, outfile, level + 1, namespace_, name_)
                    showIndent(outfile, level)
                    z = '</%s%s>\n' % (namespace_, name_)
                    z = z.encode('UTF-8')
                    outfile.write(z)
                else:
                    v = '/>\n'
                    v = v.encode('UTF-8')
                    outfile.write(v)
            def exportAttributes(xxxx, outfile, level, already_processed, namespace_='phy:', name_='Phyloxml'):
                pass
            def exportChildren(xxxx, outfile, level, namespace_='phy:', name_='Phyloxml', fromsubclass_=False):
                for phylogeny_ in xxxx.phylogeny:
                    export(phylogeny_, outfile, level, namespace_, name_='phylogeny')
            def hasContent_(xxxx):
                if hasattr(xxxx, 'phylogeny'):
                    return True
                else:
                    return False

            #### IMPORTANT ####

            # build phyloxml project
            project = Phyloxml()
            project.add_phylogeny(tree)

            # Export phyloxml document
            export(project, outfile, level, namespace_='phy:', name_='Phyloxml', namespacedef_='')
            OUTPUT= outfile


            # Some ad-hoc changes to the phyloxml formatted document to meet the schema definition
            text = OUTPUT.read().decode('UTF-8')
            text = text.replace('phy:', '')
            text = re.sub('branch_length_attr="[^"]+"', "", text)
            header = """
            <phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://www.phyloxml.org"          xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.20/phyloxml.xsd">  
            """
            text = re.sub('<Phyloxml[^>]+>', header, text)
            text = text.replace('Phyloxml', 'phyloxml')
            text = text.replace('\n', '').replace("'", '')

            self.tree_str = text
        else:
            self.tree_str = self.tree.write(format=8, format_root_node=True)

    def _get_name_phyloxml(self, node, phyloxml_species_name_tag):

        if phyloxml_species_name_tag == 'clade_name':
                if node.name != '':
                    return node.name
                else:
                    raise KeyError(
                        "Node {} in the phyloxml file {} have no clade name or phylogeny scientific name to populate the species name".format(
                            node, self.tree_file))

        elif phyloxml_species_name_tag == 'taxonomy_scientific_name':
            if node.phyloxml_clade.taxonomy[0].scientific_name != '':
                return node.phyloxml_clade.taxonomy[0].scientific_name
            else:
                raise KeyError(
                    "Node {} in the phyloxml file {} have no taxonomy scientific name  to populate the species name".format(
                        node, self.tree_file))

        elif phyloxml_species_name_tag == 'taxonomy_code':
            if node.phyloxml_clade.taxonomy[0].code != '':
                return node.phyloxml_clade.taxonomy[0].code
            else:
                raise KeyError(
                    "Node {} in the phyloxml file {} have no taxonomy scientific code  to populate the species name".format(
                        node, self.tree_file))

    def _build_tree(self, tree_file, tree_format):

        if tree_format == 'newick_string':
            self.tree_str = tree_file
            return ete3.Tree(self.tree_str, quoted_node_names=True, format=1)

        elif tree_format == 'newick':
            with open(tree_file, 'r') as nwk_file:
                self.tree_str = nwk_file.read()
            return ete3.Tree(self.tree_str, quoted_node_names=True, format=1)

        elif tree_format == 'phyloxml':
            from ete3 import Phyloxml
            project = Phyloxml()
            project.build_from_file(tree_file)
            self.tree_str = None

            tree = project.get_phylogeny()[0]

            for node in tree.traverse():

                # assign name to extant species
                if node.is_leaf():
                    node.name = self._get_name_phyloxml(node, self.phyloxml_leaf_name_tag)

                # assign name to ancestral species
                elif self.use_internal_name:
                    node.name = self._get_name_phyloxml(node, self.phyloxml_internal_name_tag)

            return tree

    def _check_consistency_names(self):

        """  
        Check if leaves names and internal node names are uniques.
        
        """

        leaf_names = []
        int_names = []

        for node in self.tree.traverse():

            if node.name == None:
                raise KeyError("{} node have no name".format(node))

            if node.is_leaf():
                if self.tree_format == 'phyloxml':
                    nn = self._get_name_phyloxml(node, self.phyloxml_leaf_name_tag)
                    if nn != None:
                        leaf_names.append(nn)
                else:
                    leaf_names.append(node.name)
            else:
                if self.tree_format == 'phyloxml':
                    nn = self._get_name_phyloxml(node,  self.phyloxml_internal_name_tag)
                    if nn != None:
                        int_names.append(nn)
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