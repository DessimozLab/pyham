from . import abstractGene


class orthoxmlParser(object):
    """
    OrthoXML parser use to read the orthoxml file containing the hogs.
    It creates on the fly the gene mapping and the abstractGene.
    """


    def __init__(self, taxonomy):
        self.current_species = None # target the species currently parse
        self.mapping_id = {}  # Map genes using unique ids (internal to orthoxml) with their external id and their species
        self.mapping_hog = {0:None} # tmp map a depth of an HOG with the last HOG visited at this level (usefull for paralogy).
        self.depth = 0 # current depth parsed
        self.genes = set() # set of all abstractGene.Gene create in this orthoxml
        self.hogs = set() # set of all abstractGene.HOG create in this orthoxml
        self.map_taxon_node = taxonomy.map_name_taxa_node # map a taxon name with a ete3 node
        return

    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = attrib["name"]

        elif tag == "{http://orthoXML.org/2011/}gene":
            self.mapping_id[attrib["id"]] = attrib
            self.mapping_id[attrib["id"]]["species"] = self.current_species

        elif tag == "{http://orthoXML.org/2011/}geneRef":

            mapping_info = self.mapping_id[attrib["id"]]

            extant_gene = abstractGene.Gene()
            extant_gene.unique_id = mapping_info["id"]
            extant_gene.species = mapping_info["species"]
            del mapping_info["id"]
            del mapping_info["species"]
            extant_gene.mapping = mapping_info
            extant_gene.extant_genome =  self.map_taxon_node[extant_gene.species]

            self.mapping_hog[self.depth].children.append(extant_gene)
            self.genes.add(extant_gene)

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":

            self.depth +=1

            current_hog = abstractGene.HOG()
            current_hog.depth = self.depth
            current_hog.parent = self.mapping_hog[self.depth - 1]

            if 'id' in attrib.keys():
                current_hog.hog_id = attrib['id']

            if current_hog.parent != None:
                current_hog.parent.children.append(current_hog)

            self.mapping_hog[self.depth] = current_hog
            self.hogs.add(current_hog)

        elif tag == "{http://orthoXML.org/2011/}property":
            self.mapping_hog[self.depth].taxon = attrib["value"]
            self.mapping_hog[self.depth].ancestral_genome =  self.map_taxon_node[attrib["value"]]

    def end(self, tag):

        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = None

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            self.depth += -1
            if self.depth == 0:
                self.mapping_hog = {0:None}
        pass

    def data(self, data):
        # Ignore data inside nodes
        pass

    def close(self):
        # Nothing special to do here
        return