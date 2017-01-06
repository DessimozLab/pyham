from . import abstractgene
from . import genome


class OrthoXMLParser(object):
    """

    OrthoXML parser use to read the orthoxml file containing the hogs.
    It creates on the fly the gene mapping and the abstractGene.
    """

    def __init__(self, taxonomy, hog_filter=None):
        self.extant_gene_map = {}
        self.current_species = None # target the species currently parse
        self.hog_stack = []
        self.toplevel_hogs = {}
        if hog_filter is None:
            hog_filter = lambda x: x
        self.filter = hog_filter
        #self.map_taxon_node = taxonomy.map_name_taxa_node
        self.taxonomy = taxonomy

    def start(self, tag, attrib):

        if tag == "{http://orthoXML.org/2011/}species":
            nodes_founded = self.taxonomy.tree.search_nodes(name=attrib['name'])

            if len(nodes_founded) == 1:
                if "genome" in nodes_founded[0].features:
                    self.current_species = nodes_founded[0].genome

                else:
                    self.current_species = genome.ExtantGenome(**attrib)
                    nodes_founded[0].add_feature("genome", self.current_species)
                    self.taxonomy.leaves.add(nodes_founded[0])
                    self.current_species.taxon = nodes_founded[0]
            else:
                print('{} node(s) founded for the species name: {}'.format(len(nodes_founded), self.current_species.name))

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
            #self.hog_stack[-1].set_taxon_range(attrib["value"])
            pass

        elif tag == "{http://orthoXML.org/2011/}score":
            self.hog_stack[-1].score(attrib['id'], float(attrib['value']))

    def end(self, tag):
        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = None

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            hog = self.hog_stack.pop()

            ## Find the taxonomic range based on the MRCA of all the children nodes

            children_genomes = set()
            children_nodes = set()

            print(hog, hog.children)

            for child in hog.children:
                if isinstance(child.taxon, genome.ExtantGenome):
                    children_genomes.add(child.taxon)
                elif isinstance(child.taxon, set):
                    children_genomes.update(child.taxon)

            for e in children_genomes:
                children_nodes.add(e.taxon)

            source = children_nodes.pop()
            common = source.get_common_ancestor(children_nodes)

            if "genome" in common.features:
                hog.set_taxon_range(common.genome)
                common.genome.add_gene(hog)

            else:
                ancestral_genome = genome.AncestralGenome()
                ancestral_genome.taxon = common
                self.taxonomy.internal_nodes.add(common)
                common.add_feature("genome", ancestral_genome)
                hog.set_taxon_range(ancestral_genome)
                ancestral_genome.add_gene(hog)

            ## End

            if len(self.hog_stack) == 0:
                filter_res = self.filter(hog)
                if filter_res:
                    self.toplevel_hogs[hog.hog_id] = hog
                    ## TODO should we delete the node and all its relatives otherwise ?

    def data(self, data):
        # Ignore data inside nodes
        pass

    def close(self):
        # Nothing special to do here
        return