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


def _build_hogs_and_genes(file_object, ham=None, filter_object=None):

    """ This function build from an orthoxml file all data that is required to build HAM object.

        Args:
            file_object (:obj:`FileObject`): File Object of the orthoxml to parse.
            ham (:obj:`str`): :obj:`HAM` used by OrthoXMLParser.
            filter_object (:obj:`ParserFilter`): :obj:`ParserFilter` use by OrthoXMLParser.

        Returns:
            :obj:`set` of top level :obj:`HOG` , :obj:`dict` of unique id with their :obj:`Gene`, :obj:`dict` of
            external id with their :obj:`Gene`.

    """

    factory = parsers.OrthoXMLParser(ham, filterObject=filter_object)
    parser = XMLParser(target=factory)

    for line in file_object:
        parser.feed(line)

    return factory.toplevel_hogs, factory.extant_gene_map, factory.external_id_mapper


def _filter_hogs_and_genes(file_object, filter):

    """ This function collect from an orthoxml file all data that is required to build HAM object based a query set.

        Args:
            file_object (:obj:`FileObject`): File Object of the orthoxml to parse.
            filter (:obj:`str`): :obj:`ParserFilter` used by FilterOrthoXMLParser.

        Returns:
            :obj:`set` of gene unique ids, :obj:`set` of top level hog id.

    """


    factory_filter = parsers.FilterOrthoXMLParser(filter)
    parser_filter = XMLParser(target=factory_filter)

    for line in file_object:
        parser_filter.feed(line)

    return set(factory_filter.geneUniqueId), set(factory_filter.hogsId)


class ParserFilter(object):
    """
    Object containing a list of queries (hogs/genes ids) that will be used by the FilterOrthoXMLParser to collect
    in the orthoxml file the required information for the OrthoXMLParser to only parse data related to this sub-dataset
    of interest.

    The ParserFilter first collect top level HOG ids, unique or external genes ids that are related to the hogs of
    interest (FilterOrthoXMLParser queries) then run the FilterOrthoXMLParser to built the list of all genes and hogs
    (OrthoXMLParser queries) required to be able to work on this subset.

    Attributes:
        HOGId_filter (:obj:`set`): :obj:`set` of HOG ids used by the FilterOrthoXMLParser.
        GeneExtId_filter (:obj:`set`): :obj:`set` of external genes ids used by the FilterOrthoXMLParser.
        GeneIntId_filter (:obj:`set`): :obj:`set` of unique genes ids used by the FilterOrthoXMLParser.

        geneUniqueId (:obj:`set`): :obj:`set` of all required unique gene ids build by the FilterOrthoXMLParser.
        hogsId (:obj:`set`): :obj:`set` of all required top level hog ids build by the FilterOrthoXMLParser.
    """

    def __init__(self):

        # Information used to select a subdataset of interest during the applyFilter call.
        self.HOGId_filter = set()
        self.GeneExtId_filter = set()
        self.GeneIntId_filter = set()

        # Information created during buildFilter call that is required by the main OrthoXMLParser for the HOGs
        # construction during HAM instantiation.
        self.geneUniqueId = None  # [geneUniqueIds]
        self.hogsId = None  # [hogIds]

    def add_hogs_via_hogId(self, list_id):
        self.HOGId_filter = self.HOGId_filter | set(map(lambda x:str(x),list_id))

    def add_hogs_via_GeneExtId(self, list_id):
        self.GeneExtId_filter = self.GeneExtId_filter | set(map(lambda x:str(x),list_id))

    def add_hogs_via_GeneIntId(self, list_id):
        self.GeneIntId_filter = self.GeneIntId_filter | set(map(lambda x:str(x),list_id))

    def buildFilter(self, file_object, type_hog_file="orthoxml"):
        """ This function will use the FilterOrthoXMLParser with the *_filter queries to build geneUniqueId and
        hogsId.

        Args:
            hog_file (:obj:`str`): Path to the file that contained the HOGs information.
            type_hog_file (:obj:`str`):  File type of the hog_file. Can be "orthoxml or "hdf5". Defaults to "orthoxml".
        """

        if type_hog_file == "orthoxml":
            self.geneUniqueId, self.hogsId = _filter_hogs_and_genes(file_object, self)
        elif type_hog_file == "hdf5":
            pass

        else:
            raise TypeError("Invalid type of hog file.")


class HAM(object):
    """
    Attributes:
        hog_file (:obj:`str`): Path to the file that contained the HOGs information.
        hog_file_type (:obj:`str`): File type of the hog_file. Can be "orthoxml or "hdf5". Defaults to "orthoxml".
        top_level_hogs (:obj:`dict`): Dictionary that map hog unique id with its list of related :obj:`HOG`.
        extant_gene_map (:obj:`dict`): Dictionary that map gene unique id with its list of related :obj:`Gene`.
        external_id_mapper (:obj:`dict`): Dictionary that map a gene external id with its list of related :obj:`HOG` or :obj:`Gene`.
        HOGMaps (:obj:`dict`): Dictionary that map a :obj:`frozenset` of a pair of genomes to its :obj:`HOGsMap`.
        filter_obj (:obj:`ParserFilter`): :obj:`ParserFilter` used during the instanciation of HAM. Defaults to None.
        taxonomy: (:obj:`Taxonomy`): :obj:`Taxonomy` build and used by :obj:`HAM` instance.

    """
    def __init__(self, newick_str, hog_file, type_hog_file="orthoxml", filter_object=None):
        """

        Args:
            newick_str (:obj:`str`): Newick str used to build the taxonomy.
            hog_file (:obj:`str`): Path to the file that contained the HOGs information.
            type_hog_file (:obj:`str`, optional): File type of the hog_file. Can be "orthoxml or "hdf5". Defaults
            to "orthoxml".
            filter_object (:obj:`ParserFilter`, optional): :obj:`ParserFilter` used during the instantiation of HAM.
            Defaults to None.
        """

        # HOGs file
        self.hog_file = hog_file
        self.hog_file_type = type_hog_file

        # Filtering
        self.filter_obj = filter_object

        # Taxonomy
        self.taxonomy = tax.Taxonomy(newick_str)
        logger.info('Build taxonomy: completed.')

        # Misc. information
        self.top_level_hogs = None
        self.extant_gene_map = None
        self.external_id_mapper = None
        self.HOGMaps = {}

        # Parsing of data
        if self.hog_file_type == "orthoxml":

            #  If filter_object specified, ham parse a first time to collect required information
            if self.filter_obj is not None:
                with open(self.hog_file, 'r') as orthoxml_file:
                    self.filter_obj.buildFilter(orthoxml_file, self.hog_file_type)
                    logger.info(
                        'Filtering Indexing of Orthoxml done: {} top level hogs and {} extant genes will be extract.'.format(
                            len(self.filter_obj.hogsId),
                            len(self.filter_obj.geneUniqueId)))

            #  This is the actual parser to build HOG/Gene and related Genomes.
            with open(self.hog_file, 'r') as orthoxml_file:
                self.top_level_hogs, self.extant_gene_map, self.external_id_mapper =\
                    _build_hogs_and_genes(orthoxml_file, self, filter_object=self.filter_obj)

            logger.info('Parse Orthoxml: {} top level hogs and {} extant genes extract.'.format(len(self.top_level_hogs),
                                                                                                len(
                                                                                                    self.extant_gene_map)))

        elif self.hog_file_type == "hdf5":
            # Looping through all orthoXML within the hdf5
            #   for each run self.build_...
            #       update self.top_level_hogs and self.extant_gene_map for each
            pass

        else:
            raise TypeError("Invalid type of hog file")

        logger.info(
            'Set up HAM analysis: ready to go with {} hogs founded within {} species.'.format(
                len(self.top_level_hogs),len(self.taxonomy.leaves)))

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

    def treeProfile(self, hog=None, outfile=None, export_with_histogram=True):

        tp = TreeProfile(self, hog=hog)

        if outfile:
            tp.export(outfile, display_internal_histogram=export_with_histogram)

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
        if hog_id in self.top_level_hogs.keys():
            return self.top_level_hogs[hog_id]
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

    def get_genes_by_external_id(self, external_gene_id):  # TODO unittest
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
        return self.top_level_hogs

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
            return g2, g1

        mrca2 = g2.get_common_ancestor(g1)
        if g1 == mrca2:
            return g1, g2

        raise TypeError("The genomes are not in the same lineage: {}".format({g1, g2}))

    def _get_ancestor_and_descendant(self, genome_set):
        """
        Get from a set of genomes the oldest genomes (or their mrca if oldest not in genomes set) and return the list of descendant (present in the genomes set only)
        :param genome_set:
        :return:
        """
        ancestor = self.get_mrca_ancestral_genome_from_genome_set(genome_set)
        genome_set.discard(ancestor)
        return ancestor, genome_set
