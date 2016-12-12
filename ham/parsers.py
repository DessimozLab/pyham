from . import abstractGene


class orthoxmlParser(object):

    def __init__(self):
        self.current_species = None
        self.mapping_id = {}
        self.mapping_hog = {0:None}
        self.depth = 0
        return

    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}orthoXML":
            pass

        elif tag == "{http://orthoXML.org/2011/}notes":
            pass

        elif tag == "{http://orthoXML.org/2011/}species":
            self.current_species = attrib["name"]

        elif tag == "{http://orthoXML.org/2011/}database":
            pass

        elif tag == "{http://orthoXML.org/2011/}gene":
            self.mapping_id[attrib["id"]] = attrib
            self.mapping_id[attrib["id"]]["species"] = self.current_species


        elif tag == "{http://orthoXML.org/2011/}groups":
            self.current_species = None

        elif tag == "{http://orthoXML.org/2011/}geneRef":

            extant_gene = abstractGene.Gene()

            info = self.mapping_id[attrib["id"]]
            extant_gene.unique_id = info["id"]
            extant_gene.species = info["species"]
            del info["id"]
            del info["species"]
            extant_gene.mapping = info
            self.mapping_hog[self.depth].children.append(extant_gene)

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":

            self.depth +=1
            current_hog = abstractGene.HOG()
            if 'id' in attrib.keys():
                current_hog.hog_id = attrib['id']
            current_hog.depth = self.depth
            current_hog.parent = self.mapping_hog[self.depth - 1]
            if current_hog.parent != None:
                current_hog.parent.children.append(current_hog)
            self.mapping_hog[self.depth] = current_hog


        elif tag == "{http://orthoXML.org/2011/}property":
            self.mapping_hog[self.depth].taxon = attrib["value"]

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
