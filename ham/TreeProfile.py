import json
import re
from string import Template
import collections
from ham.abstractgene import HOG
import logging

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
            self.treemap = self.computeTP_whole()
            self.hog = None
        elif isinstance(hog, HOG):
            self.treemap = self.computeTP_hog(hog)
            self.hog = hog
        else:
            raise KeyError("Invalid argument {} for HOG".format(hog))

    def computeTP_hog(self, hog):
        # deepcopy the required taxonomy using query hog as root level
        treeMap = hog.genome.taxon.copy(method="deepcopy")

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
            lvl.add_feature("dupl", 0)
            lvl.add_feature("lost", 0)
            lvl.add_feature("gain", 0)
            lvl.add_feature("single", 0)

        return treeMap

    def computeTP_whole(self):
        pass

    def write(self, output):

        pass
