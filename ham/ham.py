from xml.etree.ElementTree import XMLParser
from . import taxonomy as tax
from . import genome
from .TreeProfile import TreeProfile
from ham import parsers
from ham import mapper
import logging
from . import abstractgene
import copy

logger = logging.getLogger(__name__)


def _build_hogs_and_genes(file_object, taxonomy=None, filterObject=None):
    """
    build AbstractGene.HOG and AbstractGene.Gene using a given orthoXML file object
    :param file_object: orthoXML file object
    :param taxonomy: Taxonomy object used to map taxon with ete3 node
    :return: a set of AbstractGene.HOG and a set of  AbstractGene.Gene
    """


    factory = parsers.OrthoXMLParser(taxonomy, filterObject=filterObject)
    parser = XMLParser(target=factory)


    for line in file_object:

        parser.feed(line)

    return factory.toplevel_hogs, factory.extant_gene_map, factory.external_id_mapper

'''
The filter should be abstract from the type of input data (hdf5 or full orthoxml). It should only give a list
of internal/external genes ids, list of toplevel hog ids and a taxonomic range level of interest to slice the HOGS
and/or load everything contained in this level if only things specified is level.

Then in the __init__ of HAM you choose the related applyFilter based on the self.hog_file_type.
'''
class ParserFilter(object): # todo fix problem with str/int query

    def __init__(self):
        """
        ParserFilter Object pre-parse the orthoxml file with a list of queries of interest (TR, hogs, genes) and
        store the minimal information that will be required during the HAM parsing to work on the subdataset of interest.


        """

        # Information provides to select a subdataset of interest during the applyFilter call.
        self.HOGId_filter = []
        self.GeneExtId_filter = []
        self.GeneIntId_filter = []

        # Information created during buildFilter call that is required for the HOG construction during HAM instantiation.
        self.geneUniqueId = None  # [geneUniqueId]
        self.hogsId = None  # [hogId]


    def add_hogs_via_hogId(self, list_id):
        self.HOGId_filter = self.HOGId_filter + list_id

    def add_hogs_via_GeneExtId(self, list_id):
        self.GeneExtId_filter = self.GeneExtId_filter + list_id

    def add_hogs_via_GeneIntId(self, list_id):
        self.GeneIntId_filter = self.GeneIntId_filter + list_id

    def _buildFilter(self, orthoxml_file,  type_hog_file="orthoxml"):

        self.HOGId_filter = set(self.HOGId_filter)
        self.GeneExtId_filter = set(self.GeneExtId_filter)
        self.GeneIntId_filter = set(self.GeneIntId_filter)

        if type_hog_file =="orthoxml":
            factory_filter = parsers.FilterOrthoXMLParser(self)
            parser_filter = XMLParser(target=factory_filter)

            for line in orthoxml_file:
                 parser_filter.feed(line)

            self.geneUniqueId = set(factory_filter.geneUniqueId)
            self.hogsId = set(factory_filter.hogsId)

        elif type_hog_file == "hdf5":
            pass

        else:
            raise TypeError("Invalid type of hog file.")


class HAM(object):
    def __init__(self, newick_str, hog_file, type_hog_file="orthoxml", filterObject=None):

        if newick_str is None or hog_file is None:
            logger.debug('newick_str: {}, hogs_file: {}'.format(newick_str, hog_file))
            raise TypeError('Both newick string or hogs file are required to create HAM object')

        self.newick_str = newick_str
        self.hog_file = hog_file
        self.hog_file_type = type_hog_file
        self.toplevel_hogs = None
        self.extant_gene_map = None
        self.external_id_mapper = None
        self.HOGMaps = {}  # In order to not recompute two time a hog map between two level we are storing them globally
        self.filterObj = filterObject # this can be used in APi query to say if you want something outsite of the loaded information
        self.wholeTreeProfile = None

        self.taxonomy = tax.Taxonomy(self.newick_str)
        logger.info('Build taxonomy: completed.')

        if self.hog_file_type == "orthoxml":

            if self.filterObj is not None: # spend half a day before getting that two parser using the same IOBuffer of with is not working..
                with open(self.hog_file, 'r') as orthoxml_file:
                    self.filterObj._buildFilter(orthoxml_file, self.hog_file_type)
                    logger.info('Filtering Indexing of Orthoxml done: {} top level hogs and {} extant genes will be extract.'.format(len(self.filterObj.hogsId),
                                                                                                len(self.filterObj.geneUniqueId)))

            with open(self.hog_file, 'r') as orthoxml_file:
                self.toplevel_hogs, self.extant_gene_map, self.external_id_mapper = _build_hogs_and_genes(orthoxml_file, self, filterObject=self.filterObj)

            logger.info('Parse Orthoxml: {} top level hogs and {} extant genes extract.'.format(len(self.toplevel_hogs),
                                                                                                len(
                                                                                                    self.extant_gene_map)))

        elif self.hog_file_type == "hdf5":
            # Looping through all orthoXML within the hdf5
            #   for each run self.build_...
            #       update self.toplevel_hogs and self.extant_gene_map for each
            pass
        else:
            raise TypeError("Invalid type of hog file")

        logger.info(
            'Set up HAM analysis: ready to go with {} hogs founded within {} species.'.format(len(self.toplevel_hogs),len(self.taxonomy.leaves)))

    def compare_genomes(self, genomes_set, analysis):
        """
        function for level comparaison || Work in progress
        :param genomes_set:
        :param analysis: ADRIAN -> This is not smart
        :return:
        """

        if analysis == "vertical":
            if len(genomes_set) != 2:
                raise TypeError(
                    "{} genomes given for vertical HOG mapping, only 2 should be given".format(len(genomes_set)))
            vertical_map = mapper.MapVertical(self)
            vertical_map.add_map(self.get_HOGMap(genomes_set))
            return vertical_map

        elif analysis == "lateral":
            if len(genomes_set) < 2:
                raise TypeError(
                    "{} genomes given for lateral HOG mapping, at least 2 should be given".format(len(genomes_set)))
            lateral_map = mapper.MapLateral(self)
            anc, desc = self._get_ancestor_and_descendant(copy.copy(genomes_set))
            for g in desc:
                hogmap = mapper.HOGsMap(self, {g, anc})
                lateral_map.add_map(hogmap)
            return lateral_map

        else:
            raise TypeError("Invalid type of genomes comparison")

    def hogvis(self, hog, outfile=None):  # since i need the taxonomy, etc it's easier to wrap everything here
        """
        :param hog:  HOG object to visualise
        :param outfile: If specify create get_hogvis html file
        :return: the get_hogvis html string but nothing if outfile is specified
        """
        vishtml = hog.get_hogvis(self)

        if outfile is not None:
            with open(outfile, 'w') as fh:
                fh.write(vishtml)
                return

        return vishtml

    def treeProfile(self, hog=None, outfile=None):
        if hog:
            tp = TreeProfile(self, hog=hog)
        else:
            if self.wholeTreeProfile is None:
                self.wholeTreeProfile = TreeProfile(self)
            tp = self.wholeTreeProfile

        if outfile:
            tp.write(outfile)

        return tp

    def show_taxonomy(self):
        print(self.taxonomy.tree.get_ascii())

    # ... QUERY METHODS ... #

    def get_HOGMap(self, genome_pair_set):
        """
        if not already computed, do the hog mapping between a pair of levels.
        :param genome_pair_set:
        :return:
        """
        f = frozenset(genome_pair_set)
        if f in self.HOGMaps.keys():
            return self.HOGMaps[f]
        else:
            m = mapper.HOGsMap(self, genome_pair_set)
            self.HOGMaps[f] = m
            return m

    # ___ gene ___ #

    def get_hog_by_id(self, hog_id):
        """
        return a single HOG object based on a hog_id query
        :param hog_id:
        :return: an HOG object or None if wrong hog_id given
        """
        if hog_id in self.toplevel_hogs.keys():
            return self.toplevel_hogs[hog_id]
        else:
            return None

    def get_hog_by_geneId(self, gene_id):
        """
        get a single HOG object that contained the query gene
        :param gene_id: gene unique id
        :return: an HOG object that matched the query or None
        """
        qgene = self.get_gene_by_id(gene_id)
        if qgene is not None:
            return qgene.get_topLevelHog()
        else:
            return None

    def get_genes_by_external_id(self, external_gene_id): # TODO unittest
        """
        return a list because external id are not unique
        :param external_gene_id:
        :return:
        """
        if external_gene_id in self.external_id_mapper.keys():
            return [self.extant_gene_map[qgene_id] for qgene_id in self.external_id_mapper[external_gene_id]]
        else:
            logger.warning('No extant genes have the external Id {}.'.format(external_gene_id))
            return None

    def get_gene_by_id(self, gene_unique_id):
        """
        get the Gene object that match the query unique id
        :param gene_id:
        :return:
        """
        gene_unique_id = str(gene_unique_id)
        if gene_unique_id in self.extant_gene_map.keys():
            return self.extant_gene_map[gene_unique_id]
        else:
            logger.warning('No extant genes have the unique Id {}.'.format(gene_unique_id))
            return None

    def get_all_top_level_hogs(self):
        return self.toplevel_hogs

    def get_all_extant_genes_dict(self):
        return self.extant_gene_map

    # ___ genome ___ #

    def get_all_extant_genomes(self):
        """
        return all Genome.ExtantGenome
        :return: set of Genome.ExtantGenome
        """
        return set(leaf.genome for leaf in self.taxonomy.leaves)

    def get_all_ancestral_genomes(self):
        """
        return all Genome.AncestralGenome
        :return: set of Genome.AncestralGenome
        """
        return set(internal_node.genome for internal_node in self.taxonomy.internal_nodes)

    def get_ancestral_genome_by_taxon(self, tax_node):

        if "genome" in tax_node.features:
            return tax_node.genome

        else:
            ancestral_genome = genome.AncestralGenome()
            self.taxonomy.add_ancestral_genome_to_node(tax_node, ancestral_genome)

            return ancestral_genome

    def get_ancestral_genome_by_name(self, name):
        """
        :param name: str
        :return:
        """
        nodes_founded = self.taxonomy.tree.search_nodes(name=name)

        if len(nodes_founded) == 1:
            if "genome" in nodes_founded[0].features:
                return nodes_founded[0].genome
        else:
            raise ValueError('{} node(s) founded for the species name: {}'.format(len(nodes_founded), name))

    def get_extant_genome_by_name(self, **kwargs):  # todo: I'm not really happy with this
        """
        :param kwargs: use the name tag to browser a genomes
        :return:
        """

        nodes_founded = self.taxonomy.tree.search_nodes(name=kwargs['name'])

        if len(nodes_founded) == 1:
            if "genome" in nodes_founded[0].features:
                return nodes_founded[0].genome

            else:
                extant_genome = genome.ExtantGenome(**kwargs)
                self.taxonomy.add_extant_genome_to_node(nodes_founded[0], extant_genome)
                return extant_genome
        else:
            raise ValueError('{} node(s) founded for the species name: {}'.format(len(nodes_founded), kwargs['name']))

    def get_mrca_ancestral_genome_from_genome_set(self, genome_set):
        """
        return the MRCA of all the genomes present in genome_set.
        :param genome_set: set of ancestral genomes.
        :return: ancestral genome
        """
        if len(genome_set) < 2:
            raise ValueError('Only one genome is not enough to computed MRCA: {}'.format(genome_set))
        pass

        for g in genome_set:
            if not isinstance(g, genome.Genome):
                raise TypeError("expect subclass obj of '{}', got {}"
                                .format(genome.Genome.__name__,
                                        type(g).__name__))

        children_nodes = set()

        for e in genome_set:
            children_nodes.add(e.taxon)

        common = self.taxonomy.tree.get_common_ancestor(children_nodes)

        return self.get_ancestral_genome_by_taxon(common)

    def get_mrca_ancestral_genome_using_hog_children(self, hog):
        """
        collect all the genomes insides the hog children and them look for their mrca
        :param hog:
        :return:
        """

        children_genomes = set()

        for child in hog.children:
            children_genomes.add(child.genome)

        return self.get_mrca_ancestral_genome_from_genome_set(children_genomes)

    # ... PRIVATE METHODS ... #

    def _add_missing_taxon(self, child_hog, oldest_hog, missing_taxons):
        """

        :param child_hog: youngest hog
        :param oldest_hog: oldest hog
        :param missing_taxons: array of intermediate taxon sorted from most r4cent to oldest
        :return:
        """
        if not isinstance(child_hog, abstractgene.AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(abstractgene.AbstractGene.__name__,
                                    type(child_hog).__name__))

        if not isinstance(oldest_hog, abstractgene.AbstractGene):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(abstractgene.AbstractGene.__name__,
                                    type(oldest_hog).__name__))

        if oldest_hog == child_hog:
            raise TypeError("Cannot add missing level between an HOG and it self.")

        # the youngest hog is removed from the oldest hog
        oldest_hog.remove_child(child_hog)

        # Then for each intermediate level in between the two hog
        current_child = child_hog
        for tax in missing_taxons:

            # we get the related ancestral genome of this level
            ancestral_genome = self.get_ancestral_genome_by_taxon(tax)

            # we create the related hog and add it to the ancestral genome
            hog = abstractgene.HOG()
            hog.set_genome(ancestral_genome)
            ancestral_genome.add_gene(hog)

            if ancestral_genome.taxon is not current_child.genome.taxon.up:
                raise TypeError("HOG taxon {} is different than child parent taxon {}".format(ancestral_genome.taxon,
                                                                                              current_child.genome.taxon.up))

            # we add the child
            hog.add_child(current_child)
            current_child = hog

        oldest_hog.add_child(current_child)

    def _get_oldest_from_genome_pair(self, g1, g2):
        """
        get the oldest genomes from a pair of genomes... I know Adrian, I know!
        :param g1:
        :param g2:
        :return:
        """

        mrca1 = g1.get_common_ancestor(g2)
        if g2 == mrca1:
            return g2

        mrca2 = g2.get_common_ancestor(g1)
        if g1 == mrca2:
            return g1

        return None

    def _get_ancestor_and_descendant(self, genome_set):
        """
        Get from a set of genomes the oldest genomes (or their mrca if oldest not in genomes set) and return the list of descendant (present in the genomes set only)
        :param genome_set:
        :return:
        """
        ancestor = self.get_mrca_ancestral_genome_from_genome_set(genome_set)
        genome_set.discard(ancestor)
        return ancestor, genome_set
