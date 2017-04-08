from . import abstractgene
from . import genome
import copy
import logging
logger = logging.getLogger(__name__)

class OrthoXMLParser(object):
    """
    OrthoXML parser use to read the orthoxml file containing the hogs.
    It creates on the fly the gene mapping and the abstractGene.
    """

    def __init__(self, ham_object, filterObject=None):
        self.ham_object = ham_object
        self.filterObj = filterObject

        self.uhi = []

        if filterObject:
            self.uhi = copy.deepcopy(filterObject.geneUniqueId)

        # usefull information
        self.extant_gene_map = {}  # {unique_id -> gene object}
        self.external_id_mapper = {}  # {external_id -> [unique_id]} TODO unittest + more than one id may have the same ext id
        self.toplevel_hogs = {}  # {hog_id -> hog object}


        # On the fly variable
        self.cpt = 0
        self.hog_stack = []
        self.current_species = None  # target the species currently parse
        self.in_paralogGroup = None
        self.skip_this_hog = False


    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}paralogGroup" and self.skip_this_hog is False:
            self.in_paralogGroup = len(self.hog_stack)

        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = self.ham_object.get_extant_genome_by_name(**attrib)

        elif tag == "{http://orthoXML.org/2011/}gene":
            gene = abstractgene.Gene(**attrib)
            gene.set_genome(self.current_species)
            self.current_species.add_gene(gene)
            self.extant_gene_map[gene.unique_id] = gene
            for type, Id in attrib.items():
                if type is not "id":
                    self.external_id_mapper.setdefault(Id,[]).append(gene.unique_id)

        elif tag == "{http://orthoXML.org/2011/}geneRef" and self.skip_this_hog is False:

            gene = self.extant_gene_map[attrib['id']]
            self.hog_stack[-1].add_child(gene)

            # if the gene is contained within a paralogousGroup need to update its .arose_by_duplication flag. TODO unittest
            if self.in_paralogGroup == len(self.hog_stack):
                gene.arose_by_duplication=True

        elif tag == "{http://orthoXML.org/2011/}orthologGroup" and self.skip_this_hog is False:


            # in case we add a filter we have to check if this top level hog if usefull
            if len(self.hog_stack) == 0:
                if self.filterObj is not None:

                    if attrib["id"] in self.filterObj.hogsId:
                        self.skip_this_hog = False
                    else:
                        self.skip_this_hog = True

            if self.skip_this_hog is False:
                if self.in_paralogGroup == len(self.hog_stack):
                    hog = abstractgene.HOG(arose_by_duplication=True,**attrib)
                else:
                    hog = abstractgene.HOG(arose_by_duplication=False,**attrib)

                if len(self.hog_stack) > 0:
                    self.hog_stack[-1].add_child(hog)

                self.hog_stack.append(hog)

        elif tag == "{http://orthoXML.org/2011/}property" and attrib['name'] == "TaxRange":
            #self.hog_stack[-1].set_genome(attrib["value"])
            pass

        elif tag == "{http://orthoXML.org/2011/}score" and self.skip_this_hog is False:
            self.hog_stack[-1].score(attrib['id'], float(attrib['value']))

    def end(self, tag):
        if tag == "{http://orthoXML.org/2011/}species":
            logger.info("Species {} created. ".format(self.current_species.name))
            self.current_species = None

        if tag == "{http://orthoXML.org/2011/}paralogGroup" and self.skip_this_hog is False:
            self.in_paralogGroup = None

        elif tag == "{http://orthoXML.org/2011/}orthologGroup" and self.skip_this_hog is False:

            # get the latest hog
            hog = self.hog_stack.pop()

            # get the ancestral genome related to this hog based on it's children
            ancestral_genome = self.ham_object.get_mrca_ancestral_genome_using_hog_children(hog)
            hog.set_genome(ancestral_genome)
            ancestral_genome.taxon.genome.add_gene(hog)

            hog_genome = hog.genome
            change = {} # {child -> [intermediate level]}

            # for all children of this hog
            for child in hog.children:
                child_genome = child.genome
                if hog_genome.taxon.depth != child_genome.taxon.depth - 1:
                    change[child] = self.ham_object.taxonomy.get_path_up(child_genome.taxon, hog_genome.taxon)

            for hog_child, missing in change.items():
                self.ham_object._add_missing_taxon(hog_child,hog,missing)

            if len(self.hog_stack) == 0:

                self.toplevel_hogs[hog.hog_id] = hog

                self.cpt += 1
                if self.cpt % 500 == 0:
                    logger.info("{} HOGs parsed. ".format(self.cpt))

        elif tag == "{http://orthoXML.org/2011/}orthologGroup" and self.skip_this_hog is True:
            if len(self.hog_stack) == 0:
                self.skip_this_hog = False

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



    self.HOGId_filter = []
    self.GeneExtId_filter = []
    self.GeneIntId_filter = []
    self.taxonomicRange = None # you also need this for the HAM instantiation

    # Information created during buildFilter call that is required for the HOG construction during HAM instantiation.
    self.species = None  # { species -> [geneUniqueId] | "all"}
    self.hogs = None  # [hogId]

    """

    def __init__(self, filterO):

        self.filterObj = filterO

        self.geneUniqueId = []  # [geneUniqueId]
        self.hogsId = []  # [hogId]

        # reset on toplevel end
        self.current_hog = None
        self.hog_stack = []
        self.hog_generef = []
        self.add_this_hog = False


    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}gene":
            if self.filterObj.GeneIntId_filter is not []:
                if attrib['id'] in self.filterObj.GeneIntId_filter:
                    self.geneUniqueId.append(attrib['id'])

            if self.filterObj.GeneExtId_filter is not []:
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
            if attrib["id"] in self.filterObj.HOGId_filter:
                self.add_this_hog = True

            self.hog_stack.append(1)



    def end(self, tag):


        if tag == "{http://orthoXML.org/2011/}orthologGroup":

            # get the latest hog
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