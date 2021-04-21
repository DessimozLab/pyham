from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import str
from future import standard_library
standard_library.install_aliases()
from . import abstractgene
import logging
logger = logging.getLogger(__name__)
from collections import defaultdict
import numpy as np


class OrthoXMLParser(object):

    """
    Custom parser for OrthoXML containing HOGs. It can take a FilterParser object to restrict the information to parse.
     
    The parse goes through the whole XML and create on the fly the required Ham objects:
        - In the Xref/header, the parser creates the ExtantGenome and Gene objects.
        - In the Groups section, it creates the HOGs with their hierarchy (parent/children links) and their related
        AncestralGenomes.
    

    Attributes:
        ham_object (:obj:`Ham`): Ham object to feed with created objects.
        filterObj (:obj:`FilterParser`, optional): FilterParser object used to restrict the parsed information. Defaults
        to None.
        
        extant_gene_map (:obj:`dict`): dictionary of gene unique id mapped to their related :obj:`Gene`.
        external_id_mapper (:obj:`dict`): dictionary of xref id mapped to their list of possible unique ids.
        toplevel_hogs (:obj:`dict`): dictionary of top level hog id mapped to their related :obj:`HOG`.
        
        cpt (:obj:`int`): Counter of parsed hogs.
        hog_stack (:obj:`list`): Stack of hogs currently parsed. Reset at each top level hog change.
        current_species (:obj:`ExtantGenome`): Pointer to the current species parsed during xref first step.
        paralog_stack (:obj:`list`): Stack of position in hog stack when paralogy event occured and the related DuplicationNode.
        Reset at each top level hog change.
        in_paralogGroup (:obj:`int`): last position in the stack where a duplication  occured.
        skip_this_hog (:obj:`Boolean`): Boolean to skip the current hog or not (used when filtering option is set).     
    """

    def __init__(self, ham_object, filterObject=None):

        """
        Args:
            ham_object (:obj:`Ham`): Ham object to feed with created objects.
            filterObject (:obj:`FilterParser`, optional): FilterParser object used to restrict the parsed information.
            Defaults to None.
        """

        self.ham_object = ham_object
        self.filterObj = filterObject

        # usefull information
        self.extant_gene_map = {}
        self.external_id_mapper = {}
        self.toplevel_hogs = {}

        # On the fly variable
        self.cpt = 0
        self.hog_stack = []
        self.paralog_stack = []
        self.current_species = None
        self.in_paralogGroup = None
        self.paralogyNode = None
        self.skip_this_hog = False

    def _build_gene(self, attrib):
        gene = abstractgene.Gene(**attrib)
        gene.set_genome(self.current_species)
        self.current_species.add_gene(gene)
        self.extant_gene_map[gene.unique_id] = gene
        for type, Id in attrib.items():
            if type != "id":
                self.external_id_mapper.setdefault(Id, []).append(gene.unique_id)

    def _build_hog(self, attrib):

        hog = abstractgene.HOG(arose_by_duplication=False, **attrib)

        if self.in_paralogGroup == len(self.hog_stack):
            self.paralogyNode.add_child(hog)

        if len(self.hog_stack) > 0:
            self.hog_stack[-1].add_child(hog)

        self.hog_stack.append(hog)

    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = self.ham_object._get_extant_genome_by_name(**attrib)

        elif tag == "{http://orthoXML.org/2011/}gene" and self.filterObj is None:
            self._build_gene(attrib)

        elif tag == "{http://orthoXML.org/2011/}gene" and self.filterObj is not None:
            if attrib["id"] in self.filterObj.geneUniqueId:
                self._build_gene(attrib)

        elif tag == "{http://orthoXML.org/2011/}paralogGroup" and self.skip_this_hog is False:
            dNode = abstractgene.DuplicationNode(self.ham_object)

            self.paralog_stack.append({'depth': len(self.hog_stack), 'node': dNode})

            self.in_paralogGroup = self.paralog_stack[-1]['depth']
            self.paralogyNode = self.paralog_stack[-1]['node']

        elif tag == "{http://orthoXML.org/2011/}geneRef" and self.skip_this_hog is False:

            gene = self.extant_gene_map[attrib['id']]
            self.hog_stack[-1].add_child(gene)

            # if the gene is contained within a paralogousGroup need to update its .arose_by_duplication flag.
            if self.in_paralogGroup == len(self.hog_stack):
                self.paralogyNode.add_child(gene)

        elif tag == "{http://orthoXML.org/2011/}orthologGroup" and self.skip_this_hog is False:
            if len(self.hog_stack) == 0 and self.filterObj is not None:

                if str(attrib["id"]) in self.filterObj.hogsId:
                    self.skip_this_hog = False
                    self._build_hog(attrib)
                else:
                    self.skip_this_hog = True
                    self.hog_stack.append(0)

            else:
                self._build_hog(attrib)



        elif tag == "{http://orthoXML.org/2011/}orthologGroup" and self.skip_this_hog is True:
            self.hog_stack.append(0)

        elif tag == "{http://orthoXML.org/2011/}property" and not self.skip_this_hog:
            self.hog_stack[-1].add_property(attrib["name"], attrib["value"])

        elif tag == "{http://orthoXML.org/2011/}score" and not self.skip_this_hog:
            self.hog_stack[-1].score(attrib['id'], float(attrib['value']))

    def end(self, tag):

        if tag == "{http://orthoXML.org/2011/}species":
            logger.info("Species {} created. ".format(self.current_species.name))
            self.current_species = None

        elif tag == "{http://orthoXML.org/2011/}paralogGroup" and self.skip_this_hog is False:

            ln = self.paralog_stack.pop()

            ln['node'].set_MRCA()

            if len(self.paralog_stack) > 0:
                self.in_paralogGroup = self.paralog_stack[-1]['depth']
                self.paralogyNode = self.paralog_stack[-1]['node']
            else:
                self.in_paralogGroup = None
                self.paralogyNode = None

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":

            # get the latest hog
            hog = self.hog_stack.pop()

            if self.skip_this_hog:
                if len(self.hog_stack) == 0:
                    self.skip_this_hog = False
            else:
                try:
                    # get the ancestral genome related to this hog based on it's children
                    if len(set([child.genome for child in hog.children])) == 1:
                        try:
                            if hog['TaxRange'] == hog.children[0].genome.name:
                                # remove this hog structure.
                                for child in hog.children:
                                    child.parent = hog.parent
                                hog.parent.children.remove(hog)
                                hog.parent.children.extend(hog.children)
                                if hog.arose_by_duplication != False:
                                    dupl_node = hog.arose_by_duplication
                                    dupl_node.remove_child(hog)
                                    for child in hog.children:
                                        dupl_node.add_child(child)
                                    self.paralog_stack[-1]['depth'] -= 1
                                return
                        except KeyError:
                            pass
                        ancestral_genome = self.ham_object._get_ancestral_genome_by_taxon(hog.children[0].genome.taxon.up)
                    else:
                        ancestral_genome = self.ham_object._get_ancestral_genome_by_mrca_of_hog_children_genomes(hog)

                    hog.set_genome(ancestral_genome)
                    ancestral_genome.taxon.genome.add_gene(hog)

                    # get all child clustered by dup if any
                    child_by_duplication = defaultdict(list)
                    for child in hog.children:
                        if child.arose_by_duplication != False:
                            child_by_duplication[child.arose_by_duplication].append(child)

                    # For each duplication
                    for duplication, children in child_by_duplication.items():

                        # add MRCA hog if its missing
                        if duplication.MRCA != hog.genome:

                            # create the MRCA hog
                            mrcahog = abstractgene.HOG()
                            mrcahog.set_genome(duplication.MRCA)
                            duplication.MRCA.add_gene(mrcahog)

                            hog.add_child(mrcahog)

                            # add missing level down (from mrcaHOG to duplicated child)
                            for child in children:
                                hog.remove_child(child)
                                mrcahog.add_child(child)

                                change = self.ham_object.taxonomy.get_path_up(child.genome.taxon, mrcahog.genome.taxon)
                                self.ham_object._add_missing_taxon(child, mrcahog, change)

                            duplication.set_parent(mrcahog)
                            for x in list(duplication.children):
                                duplication.remove_child(x)
                            for y in list(mrcahog.children):
                                duplication.add_child(y)

                        # Otherwise simply add missing taxa between hog and duplicated child
                        else:
                            duplication.set_parent(hog)
                            for child_direct in children:
                                change = self.ham_object.taxonomy.get_path_up(child_direct.genome.taxon, hog.genome.taxon)
                                self.ham_object._add_missing_taxon(child_direct, hog, change)

                                duplication.remove_child(child_direct)
                                duplication.add_child(hog.children[-1])

                    hog_genome = hog.genome
                    change = {} # {child -> [intermediate level]}

                    # for all children of this hog find missing level
                    for child in hog.children:
                        child_genome = child.genome
                        if hog_genome.taxon.depth != child_genome.taxon.depth - 1:
                            change[child] = self.ham_object.taxonomy.get_path_up(child_genome.taxon, hog_genome.taxon)

                    # and add them if required
                    for hog_child, missing in change.items():
                        self.ham_object._add_missing_taxon(hog_child, hog, missing)

                    if len(self.hog_stack) == 0:
                        self.toplevel_hogs[hog.hog_id] = hog
                        self.cpt += 1
                        if self.cpt % 500 == 0:
                            logger.info("{} HOGs parsed. ".format(self.cpt))
                except Exception as e:
                    logger.error("cannot parse HOG '{}': {}".format(hog.hog_id, e))
                    raise

    def data(self, data):
        # Ignore data inside nodes
        pass

    def close(self):
        # Nothing special to do here
        return


class FilterOrthoXMLParser(object):
    """
    Custom OrthoXML parser use to read the orthoxml file and get required information base on specific query set. It use
    the ParserFilter object as input to do the selection.

    The parse goes through the whole XML and create on the fly the required Ham objects:
        - In the Xref/header, the parser creates the ExtantGenome and Gene objects.
        - In the Groups section, it creates the HOGs with their hierarchy (parent/children links) and their related
        AncestralGenomes.


    Attributes:
        filterObj (:obj:`FilterParser`): Filter object with all the filtering information.

        geneUniqueId (:obj:`list`): list of unique gene ids to store for the OrthoXMLParser.
        hogsId (:obj:`list`): list of hogs ids to store for the OrthoXMLParser.

        current_hog (:obj:`int`): Current hog id.
        hog_stack (:obj:`list`): Stack of hogs currently parsed. Reset at each top level hog change.
        hog_generef (:obj:`list`): list of unique gene ids contained into the current parsed hog.
        add_this_hog (:obj:`Boolean`): Boolean to know if we keep the current hog for the OrthoXMLParser.
    """

    def __init__(self, filterO):

        self.filterObj = filterO

        self.geneUniqueId = []
        self.hogsId = []

        self.current_hog = None
        self.hog_stack = []
        self.hog_generef = []
        self.add_this_hog = False

    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}gene":
            if self.filterObj.GeneIntId_filter:
                if attrib['id'] in self.filterObj.GeneIntId_filter:
                    self.geneUniqueId.append(attrib['id'])

            if self.filterObj.GeneExtId_filter:
                for xtid in attrib.values():
                    if xtid in self.filterObj.GeneExtId_filter:
                        self.geneUniqueId.append(attrib['id'])

        elif tag == "{http://orthoXML.org/2011/}geneRef":
            self.hog_generef.append(attrib['id'])
            if attrib['id'] in self.geneUniqueId:
                self.add_this_hog = True

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            if len(self.hog_stack) == 0:
                self.current_hog = attrib["id"]
                if self.current_hog in self.filterObj.HOGId_filter:
                    self.add_this_hog = True

            self.hog_stack.append(1)

    def end(self, tag):

        if tag == "{http://orthoXML.org/2011/}orthologGroup":

            hog = self.hog_stack.pop()

            if len(self.hog_stack) == 0:

                if self.add_this_hog:
                    self.geneUniqueId = self.geneUniqueId + self.hog_generef
                    self.hogsId.append(self.current_hog)

                self.current_hog = None
                self.hog_generef = []
                self.add_this_hog = False

    def data(self, data):
        # Ignore data inside nodes
        pass

    def close(self):
        # Nothing special to do here
        return