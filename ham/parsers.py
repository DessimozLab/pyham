from . import abstractgene
from . import genome
import logging
logger = logging.getLogger(__name__)

class OrthoXMLParser(object):
    """

    OrthoXML parser use to read the orthoxml file containing the hogs.
    It creates on the fly the gene mapping and the abstractGene.
    """

    def __init__(self, ham_object, hog_filter=None):
        self.extant_gene_map = {} # {unique_id -> [mapped_id]}
        self.current_species = None  # target the species currently parse
        self.hog_stack = []
        self.toplevel_hogs = {}
        if hog_filter is None:
            hog_filter = lambda x: x
        self.filter = hog_filter
        self.ham_object = ham_object
        self.in_paralogGroup = None
        self.cpt = 0

    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}paralogGroup":
            self.in_paralogGroup = len(self.hog_stack)

        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = self.ham_object.get_extant_genome_by_name(**attrib)

        elif tag == "{http://orthoXML.org/2011/}gene":
            gene = abstractgene.Gene(**attrib)
            gene.set_genome(self.current_species)
            self.current_species.add_gene(gene)
            self.extant_gene_map[gene.unique_id] = gene


        elif tag == "{http://orthoXML.org/2011/}geneRef":
            gene = self.extant_gene_map[attrib['id']]
            self.hog_stack[-1].add_child(gene)

            # if the gene is contained within a paralogousGroup need to update its .arose_by_duplication flag. TODO unittest
            if self.in_paralogGroup == len(self.hog_stack):
                gene.arose_by_duplication=True

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
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

        elif tag == "{http://orthoXML.org/2011/}score":
            self.hog_stack[-1].score(attrib['id'], float(attrib['value']))

    def end(self, tag):
        if tag == "{http://orthoXML.org/2011/}species":
            logger.info("Species {} created. ".format(self.current_species.name))
            self.current_species = None

        if tag == "{http://orthoXML.org/2011/}paralogGroup":
            self.in_paralogGroup = None

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":

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

                filter_res = self.filter(hog)
                if filter_res:
                    self.toplevel_hogs[hog.hog_id] = hog
                    ## TODO should we delete the node and all its relatives otherwise ?

                self.cpt += 1
                if self.cpt % 500 == 0:
                    logger.info("{} HOGs parsed. ".format(self.cpt))

    def data(self, data):
        # Ignore data inside nodes
        pass

    def close(self):
        # Nothing special to do here
        return