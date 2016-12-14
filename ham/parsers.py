from . import abstractgene
from . import genome


class OrthoXMLParser(object):
    '''
    OrthoXML parser use to read the orthoxml file containing the hogs.
    It creates on the fly the gene mapping and the abstractGene.
    '''

    def __init__(self, hog_filter=None):
        self.extant_gene_map = {}
        self.current_species = None # target the species currently parse
        self.hog_stack = []
        self.toplevel_hogs = {}
        if hog_filter is None:
            hog_filter = lambda x: x
        self.filter = hog_filter

    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = genome.ExtantGenome(**attrib)

        elif tag == "{http://orthoXML.org/2011/}gene":
            gene = abstractgene.Gene(**attrib)
            self.current_species.add_gene(gene)
            self.extant_gene_map[gene.unique_id] = gene

        elif tag == "{http://orthoXML.org/2011/}geneRef":
            gene = self.extant_gene_map[attrib['id']]
            self.hog_stack[-1].add_child(gene)

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            hog = abstractgene.HOG(**attrib)
            if len(self.hog_stack) > 0:
                self.hog_stack[-1].add_child(hog)
            self.hog_stack.append(hog)

        elif tag == "{http://orthoXML.org/2011/}property" and attrib['name'] == "TaxRange":
            self.hog_stack[-1].set_taxon_range(attrib["value"])

    def end(self, tag):
        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = None

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            hog = self.hog_stack.pop()
            if len(self.hog_stack) == 0:
                filter_res = self.filter(hog)
                if filter_res:
                    self.toplevel_hogs[hog.hog_id] = hog

    def data(self, data):
        # Ignore data inside nodes
        pass

    def close(self):
        # Nothing special to do here
        return