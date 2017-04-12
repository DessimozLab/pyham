import json
import re
from string import Template
import collections
import sys
from ham.abstractgene import HOG
import logging
from ete3 import TreeStyle, TextFace, NodeStyle, BarChartFace
import random
logger = logging.getLogger(__name__)

class TreeProfile(object):
    """
    Object that map on a given taxonomy related evolutionary event at each nodes (numbers of ancestral genes,
    duplications, lost, etc...). This can be applied on a single gene family (hog) to represent the evolutionary history
    of those genes or can be used to represent the ancestral states of a multiple genome setup.

    self.ham:  HAM object with related taxonomy object
    self.treemap: ETE3 Tree with required taxonomy and related information at each internal node

     if the TreeProfile is instanciate with only the "ham" argument it will compute the tree profile for the
     whole genomic setup. If the "hog" argument is provided the TreeProfile will be compute for a single hog.
    """
    def __init__(self, ham, hog=None):
        self.ham = ham

        if hog is None:
            self.hog = None
            self.treemap = self.computeTP_whole()

        elif isinstance(hog, HOG):
            self.hog = hog
            self.treemap = self.computeTP_hog(hog)
        else:
            raise KeyError("Invalid argument {} for HOG".format(hog))

    def computeTP_hog(self, hog):
        # deepcopy the required taxonomy using query hog as root level
        treeMap = self.ham.taxonomy.tree.copy(method="newick")

        # create a dictionary that map node with related hogs/genes
        levelGroups = {}

        # add all of subhog to the related level in levelGroups
        for subhog in hog.get_all_descendant_hogs():
            levelGroups.setdefault(subhog.genome.name, []).append(subhog)

        # add empty extant genome to levelGroups
        for extantGenome in treeMap.get_leaves():
            levelGroups[extantGenome.name] = []

        # add genes to related extant genome in levelGroups
        for species, genes in hog.get_all_descendant_genes_clustered_by_species().items():
            levelGroups[species.name] = genes

        # add to each node the number of ancestral|extant genes
        for lvl in treeMap.traverse():
            lvl.add_feature("nbr_genes", len(levelGroups[lvl.name]))
            lvl.add_feature("dupl", None)
            lvl.add_feature("lost", None)
            lvl.add_feature("gain", None)
            lvl.add_feature("single", None)

        return treeMap

    def computeTP_whole(self):

        def _add_annot(node, nbr, dupl, lost, gain, single):
            node.add_feature("nbr_genes", nbr)
            node.add_feature("dupl", dupl)
            node.add_feature("lost", lost)
            node.add_feature("gain", gain)
            node.add_feature("single", single)

        treeMap = self.ham.taxonomy.tree.copy(method="newick")

        for node in treeMap.traverse():
            if node.is_root():
                node_genome = self.ham.get_ancestral_genome_by_name(node.name)
                _add_annot(node, len(node_genome.genes), None, None, None, None)

            else:
                node_genome_up = self.ham.get_ancestral_genome_by_name(node.up.name)

                if node.is_leaf():
                    node_genome = self.ham.get_extant_genome_by_name(name=node.name)
                    nbr = node_genome.get_number_genes(singleton=False)

                else:
                    node_genome = self.ham.get_ancestral_genome_by_name(node.name)
                    nbr = node_genome.get_number_genes()
                hogmap = self.ham.get_HOGMap({node_genome, node_genome_up})
                _add_annot(node, nbr, len(hogmap.DUPLICATE.keys()), len(hogmap.LOSS), len(hogmap.GAIN), len(hogmap.SINGLE.keys()))

        return treeMap

    def export(self, output, layout_function=None, display_internal_histogram=False):

        # maximum number of genes per node in this treeMap
        max_genes = max([d for d in self.treemap.traverse()], key=lambda x:x.nbr_genes).nbr_genes

        def _layout(node):

            _color_scheme = ["#41c1c2","#bdc3c7","#f39c12","#27ae60","#e74c3c"]
            _label = ["Genes","Identicals","Duplicated","Novel","Lost"]

            def _add_face(name_feature, value_feature, cnum=1, pos="branch-right"):
                node.add_face(TextFace("{}: {}".format(name_feature, value_feature)), column=cnum, position=pos)

            def _add_faces(cNbr=1, posNbr="branch-right", cAttr=2, posAtt="branch-right"):

                _add_face("#genes", node.nbr_genes, cnum=cNbr, pos=posNbr)

                if node.single is not None:
                    _add_face("#Identical", node.single, cnum=cAttr, pos=posAtt)

                if node.dupl is not None:
                    _add_face("#Duplicated", node.dupl, cnum=cAttr, pos=posAtt)

                if node.gain is not None:
                    _add_face("#Novel", node.gain, cnum=cAttr, pos=posAtt)

                if node.lost is not None:
                    _add_face("#Lost", node.lost, cnum=cAttr, pos=posAtt)

            if node.is_leaf():
                if display_internal_histogram:
                    node.add_face(BarChartFace([node.nbr_genes,node.single,node.dupl,node.gain,node.lost], deviations=None, width=50, height=25, colors=_color_scheme, labels=_label, min_value=0, max_value=max_genes, label_fsize=6, scale_fsize=6),column=1, position = "branch-right")
                else:
                    _add_faces()

            else:
                if display_internal_histogram:
                    node.add_face(BarChartFace([node.nbr_genes,node.single,node.dupl,node.gain,node.lost], deviations=None, width=50, height=25, colors=_color_scheme, labels=_label, min_value=0, max_value=max_genes, label_fsize=6, scale_fsize=6),column=0, position = "branch-bottom")
                else:
                    _add_faces(cNbr=0, posNbr="branch-top", cAttr=0, posAtt="branch-bottom")

        ts = TreeStyle()

        if layout_function is not None:
            ts.layout_fn = layout_function
        else:
            ts.layout_fn = _layout

        self.treemap.render(output,tree_style=ts)

    def dirty_display(self): # todo to be removed only for dev purposed
        if self.hog:
            att = ["name", "nbr_genes"]
        else:
            att = ["name", "nbr_genes","single", "dupl", "lost", "gain"]
        print(att)
        self.treemap.show()

        #f = open('data_new.txt', 'w')
        #f.write(self.treemap.get_ascii(show_internal=True, compact=True, attributes=att))
        #return self.treemap.get_ascii(show_internal=True, compact=True, attributes=att)

