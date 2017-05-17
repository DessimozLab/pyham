from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import map
from builtins import open
from builtins import str
from future import standard_library
standard_library.install_aliases()

from xml.etree.ElementTree import XMLParser
from . import taxonomy as tax
from .genome import Genome,AncestralGenome, ExtantGenome
from . import parsers
from . import mapper
from . import abstractgene
from .TreeProfile import TreeProfile
import logging
import copy

logger = logging.getLogger(__name__)


class ParserFilter(object):
    """
    Object containing a list of queries (hogs/genes ids) that will be used by the FilterOrthoXMLParser to collect
    in the orthoxml file the required information for the OrthoXMLParser to only parse data related to this sub-dataset
    of interest.

    The ParserFilter first collect top level HOG ids, unique or external genes ids that are related to the hogs of
    interest (FilterOrthoXMLParser queries) then run the FilterOrthoXMLParser to built the list of all genes and hogs
    (OrthoXMLParser queries) required to be able to work on this subset.

    Attributes:
        | HOGId_filter (:obj:`set`): :obj:`set` of HOG ids used by the FilterOrthoXMLParser.
        | GeneExtId_filter (:obj:`set`): :obj:`set` of external genes ids used by the FilterOrthoXMLParser.
        | GeneIntId_filter (:obj:`set`): :obj:`set` of unique genes ids used by the FilterOrthoXMLParser.

        | geneUniqueId (:obj:`set`): :obj:`set` of all required unique gene ids build by the FilterOrthoXMLParser.
        | hogsId (:obj:`set`): :obj:`set` of all required top level hog ids build by the FilterOrthoXMLParser.
    """

    def __init__(self):

        # Information used to select a subdataset of interest during the applyFilter call.
        self.HOGId_filter = set()
        self.GeneExtId_filter = set()
        self.GeneIntId_filter = set()

        # Information created during buildFilter call that is required by the main OrthoXMLParser for the HOGs
        # construction during Ham instantiation.
        self.geneUniqueId = None  # [geneUniqueIds]
        self.hogsId = None  # [hogIds]

    def add_hogs_via_hogId(self, list_id):
        """
        :param list_id:
        """
        self.HOGId_filter = self.HOGId_filter | set(map(lambda x: str(x), list_id))

    def add_hogs_via_GeneExtId(self, list_id):
        """
                :param list_id:
                """
        self.GeneExtId_filter = self.GeneExtId_filter | set(map(lambda x: str(x), list_id))

    def add_hogs_via_GeneIntId(self, list_id):
        """
                :param list_id:
                """
        self.GeneIntId_filter = self.GeneIntId_filter | set(map(lambda x: str(x), list_id))

    def buildFilter(self, file_object, type_hog_file="orthoxml"):
        """ This function will use the FilterOrthoXMLParser with the *_filter queries to build geneUniqueId and
        hogsId.

        Args:
            | hog_file (:obj:`str`): Path to the file that contained the HOGs information.
            | type_hog_file (:obj:`str`):  File type of the hog_file. Can be "orthoxml or "hdf5". Defaults to "orthoxml".
        """

        if type_hog_file == "orthoxml":
            self.geneUniqueId, self.hogsId = self._filter_hogs_and_genes(file_object)
        elif type_hog_file == "hdf5":
            pass

        else:
            raise TypeError("Invalid type of hog file.")

    def _filter_hogs_and_genes(self, file_object):

        """ This function collect from an orthoxml file all data that is required to build Ham object based this filter
            object.

            Args:
                | file_object (:obj:`FileObject`): File Object of the orthoxml to parse.

            Returns:
                | :obj:`set` of gene unique ids, :obj:`set` of top level hog id.

        """

        factory_filter = parsers.FilterOrthoXMLParser(self)
        parser_filter = XMLParser(target=factory_filter)

        for line in file_object:
            parser_filter.feed(line)

        return set(factory_filter.geneUniqueId), set(factory_filter.hogsId)


class Ham(object):
    """
    
    Attributes:
        | hog_file (:obj:`str`): Path to the file that contained the HOGs information.
        | hog_file_type (:obj:`str`): File type of the hog_file. Can be "orthoxml or "hdf5". Defaults to "orthoxml".
        | top_level_hogs (:obj:`dict`): Dictionary that map hog unique id with its list of related :obj:`pyham.abstractgene.HOG`.
        | extant_gene_map (:obj:`dict`): Dictionary that map gene unique id with its list of related :obj:`pyham.abstractgene.Gene`.
        | external_id_mapper (:obj:`dict`): Dictionary that map a gene external id with its list of related :obj:`pyham.abstractgene.HOG` or :obj:`pyham.abstractGene.gene`.
        | HOGMaps (:obj:`dict`): Dictionary that map a :obj:`frozenset` of a pair of genomes to its :obj:`pyham.mapper.HOGsMap`.
        | filter_obj (:obj:`pyham.pyham.ParserFilter`): :obj:`ParserFilter` used during the instanciation of Ham. Defaults to None.
        | taxonomy: (:obj:`pyham.mapper.Taxonomy`): :obj:`pyham.pyham.Taxonomy` build and used by :obj:`pyham.pyham.Ham` instance.

    """

    def __init__(self, newick_str, hog_file, type_hog_file="orthoxml", filter_object=None, use_internal_name=False):
        """

        Args:
            | newick_str (:obj:`str`): Newick str used to build the taxonomy.
            | hog_file (:obj:`str`): Path to the file that contained the HOGs information.
            | type_hog_file (:obj:`str`, optional): File type of the hog_file. Can be "orthoxml or "hdf5". Defaults to "orthoxml".
            | filter_object (:obj:`pyham.pyham.ParserFilter`, optional): :obj:`pyham.pyham.ParserFilter` used during the instantiation of pyham.pyham.Ham. Defaults to None.
            | use_internal_name (:obj:`Boolean`, optional): Set to decide to use or not the internal naming of the given newick string. This should be set to False when support values are provided in the newick. Defaults to False.
        """

        # HOGs file
        self.hog_file = hog_file
        self.hog_file_type = type_hog_file

        # Filtering
        if isinstance(filter_object, ParserFilter) or filter_object is None:
            self.filter_obj = filter_object
        else:
            raise TypeError("filter_obj should be '{}', got {}"
                            .format(ParserFilter.__name__,
                                    type(filter_object).__name__))

        # Taxonomy
        self.taxonomy = tax.Taxonomy(newick_str, use_internal_name=use_internal_name)
        logger.info('Build taxonomy: completed.')

        # Misc. information
        self.top_level_hogs = None
        self.extant_gene_map = None
        self.external_id_mapper = None
        self.HOGMaps = {}

        # Parsing of data
        if self.hog_file_type == "orthoxml":

            #  If filter_object specified, pyham parse a first time to collect required information
            if self.filter_obj is not None:
                with open(self.hog_file, 'r') as orthoxml_file:
                    self.filter_obj.buildFilter(orthoxml_file, self.hog_file_type)
                    logger.info(
                        'Filtering Indexing of Orthoxml done: {} top level hogs and {} extant genes will be extract.'.format(
                            len(self.filter_obj.hogsId),
                            len(self.filter_obj.geneUniqueId)))

            # This is the actual parser to build HOG/Gene and related Genomes.
            with open(self.hog_file, 'r') as orthoxml_file:
                self.top_level_hogs, self.extant_gene_map, self.external_id_mapper = \
                    self._build_hogs_and_genes(orthoxml_file, filter_object=self.filter_obj)

            logger.info(
                'Parse Orthoxml: {} top level hogs and {} extant genes extract.'.format(len(self.top_level_hogs),
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
            'Set up Ham analysis: ready to go with {} hogs founded within {} species.'.format(
                len(self.top_level_hogs), len(self.taxonomy.leaves)))

    # ... TOOLS ... #

    def compare_genomes_vertically(self, genome1, genome2):

        """
        Function to compute a :obj:`MapVertical` based on the 2 given genomes. The genome order doesn't matter.

        Attributes:
            | genome1 (:obj:`pyham.genome.Genome`): First :obj:`pyham.genome.Genome` to compare.
            | genome2 (:obj:`pyham.genome.Genome`): Second :obj:`pyham.genome.Genome` to compare.

        Returns:
            | :obj:`pyham.mapper.MapVertical`.
        
        Raises:
            | TypeError: if genome1 and genome2 are not :obj:`pyham.genome.Genome`.
        """

        if not isinstance(genome1, Genome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Genome.__name__,
                                    type(genome1).__name__))

        if not isinstance(genome2, Genome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Genome.__name__,
                                    type(genome2).__name__))

        vertical_map = mapper.MapVertical(self)
        vertical_map.add_map(self._get_HOGMap({genome1, genome2}))

        return vertical_map

    def compare_genomes_lateral(self, genome1, genome2):

        """
        Function to compute a :obj:`pyham.mapper.MapLateral` based on the 2 given genomes. The genome order doesn't matter.

        Attributes:
            | genome1 (:obj:`pyham.genome.Genome`): First :obj:`pyham.genome.Genome` to compare.
            | genome2 (:obj:`pyham.genome.Genome`): Second :obj:`pyham.genome.Genome` to compare.

        Returns:
            | :obj:`pyham.mapper.MapLateral`.

        Raises:
            | TypeError: if genome1 and genome2 are not :obj:`pyham.genome.Genome`.
        """

        if not isinstance(genome1, Genome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Genome.__name__,
                                    type(genome1).__name__))

        if not isinstance(genome2, Genome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Genome.__name__,
                                    type(genome2).__name__))

        lateral_map = mapper.MapLateral(self)
        anc, desc = self._get_ancestor_and_descendant(copy.copy({genome1, genome2}))
        for g in desc:
            hogmap = mapper.HOGsMap(self, g, anc)
            lateral_map.add_map(hogmap)

        return lateral_map

    def create_hog_visualisation(self, hog, outfile=None):

        """
        Function to compute a :obj:`pyham.Hogvis`.

        If an outfile is specified, export the :obj:`pyham.Hogvis` as html file.

        Attributes:
            | hog (:obj:`pyham.abstractgene.HOG`): HOG use as template for the :obj:`pyham.Hogvis`.
            | outfile (:obj:`str`, optional): Path to the Hogvis html file.

        Returns:
            | :obj:`pyham.Hogvis` 
        """

        newick_tree =self.taxonomy.get_newick_from_tree(hog.genome.taxon)

        vis = hog.get_hog_vis(newick_tree)

        if outfile is not None:
            with open(outfile, 'w') as fh:
                fh.write(vis.renderHTML)

        return vis

    def create_tree_profile(self, hog=None, outfile=None, export_with_histogram=True):

        """
        Function to compute a :obj:`pyham.TreeProfile`.
        
        If no hog are given the tree profile will be created for the whole Ham setup (all internal nodes with all HOGs).
        Otherwise, the tree profile is build for the specific hog given.
        
        If an outfile is specified, export the create_tree_profile as image into file.

        Attributes:
            | hog (:obj:`pyham.abstractgene.HOG`, optional): HOG use as template for the create_tree_profile.
            | outfile (:obj:`str`, optional): Path to the create_tree_profile output image file. valid extensions are .SVG, .PDF, .PNG.  
            | export_with_histogram (:obj:`Bool`, optional): If True, export image with histogram at each internal node otherwise 
            | display internal node information as text.

        Returns:
            | :obj:`pyham.TreeProfile` 
        """

        tp = TreeProfile(self, hog=hog)

        if outfile:
            tp.export(outfile, display_internal_histogram=export_with_histogram)

        return tp

    def get_ascii_taxonomy(self):
        return self.taxonomy.tree.get_ascii()

    # ... QUERY METHODS ... #

    # ___ Gene ___ #

    def get_gene_by_id(self, gene_unique_id):

        """  Get the :obj:`pyham.abstractgene.Gene` that match the query unique gene Id.

            Args:
                | gene_unique_id (:obj:`str` or :obj:`int`): Unique gene Id.

            Returns:
                :obj:`pyham.abstractgene.Gene`
            
            Raises: 
                KeyError is not `pyham.abstractgene.Gene` match the id.

        """
        gene_unique_id = str(gene_unique_id)

        if gene_unique_id in self.extant_gene_map.keys():
            return self.extant_gene_map[gene_unique_id]

        raise KeyError('Id {} cannot match any Gene unique Id.'.format(gene_unique_id))

    def get_genes_by_external_id(self, external_gene_id):

        """  Get the list of :obj:`pyham.abstractgene.Gene` that match the query external gene Id.

            Args:
                external_gene_id (:obj:`str` or :obj:`int`): External gene Id.

            Returns:
                a list of :obj:`pyham.abstractgene.Gene`
            
            Raises:
                 KeyError if no `pyham.abstractgene.Gene` match id.

        """

        external_gene_id = str(external_gene_id)

        if external_gene_id in self.external_id_mapper.keys():
            return [self.extant_gene_map[qgene_id] for qgene_id in self.external_id_mapper[external_gene_id]]

        raise KeyError('Id {} cannot match any Gene external Id.'.format(external_gene_id))

    def get_list_extant_genes(self):

        """  Get the list of all :obj:`pyham.abstractgene.Gene`.

            Returns:
                a list of :obj:`pyham.abstractgene.Gene`.

        """

        return list(self.extant_gene_map.values())

    def get_dict_extant_genes(self):

        """  Get a dictionary that map all unique gene id with their related :obj:`pyham.abstractgene.Gene`.

            Returns:
                a dictionary mapping unique gene Id (:obj:`str`) with :obj:`pyham.abstractgene.Gene`.

        """

        return self.extant_gene_map

    # ___ HOG ___ #

    def get_hog_by_id(self, hog_id):

        """ Get the top level :obj:`HOG` that match the hog id query.

            Args:
                hog_id (:obj:`str` or :obj:`int`): Top level HOG id.

            Returns:
                :obj:`pyham.abstractgene.HOG`
            
            Raises:
                 KeyError if id match no `pyham.abstractgene.HOG`.

        """

        hog_id = str(hog_id)

        if hog_id in self.top_level_hogs.keys():
            return self.top_level_hogs[hog_id]

        raise KeyError(' Id {} cannot match any HOG Id.'.format(hog_id))

    def get_hog_by_gene(self, gene):

        """  Get the top level :obj:`HOG` that contain the query :obj:`pyham.abstractgene.Gene`. If the :obj:`pyham.abstractgene.Gene` is a singleton it will 
        return itself.

            Args:
                gene (:obj:`pyham.abstractgene.Gene`): :obj:`pyham.abstractgene.Gene` object.

            Returns:
                :obj:`pyham.abstractgene.HOG`
            
            Raises:
                 KeyError is gene is not a :obj:`pyham.abstractgene.Gene`.

        """

        if isinstance(gene, abstractgene.Gene):
            return gene.get_top_level_hog()

        raise KeyError("expect a '{}' as query, got {}".format(abstractgene.Gene, type(gene).__name__))

    def get_list_top_level_hogs(self):

        """  Get the list of all the top level :obj:`pyham.abstractgene.HOG`.

            Returns:
                a list of :obj:`pyham.abstractgene.HOG`.

        """

        return list(self.top_level_hogs.values())

    def get_dict_top_level_hogs(self):

        """  Get a dictionary that map all top level hog id with their related :obj:`pyham.abstractgene.HOG`.

            Returns:
                a dictionary mapping hog Id (:obj:`str`) with :obj:`pyham.abstractgene.HOG`.

        """

        return self.top_level_hogs

    # ___ ExtantGenome ___ #

    def get_list_extant_genomes(self):

        """  
        Get the list of all :obj:`pyham.genome.ExtantGenome` created during the parsing.

            Returns:
                a list of :obj:`pyham.genome.ExtantGenome`.

        """

        return [leaf.genome for leaf in self.taxonomy.leaves]

    def get_extant_genome_by_name(self, name):

        """  
        Get the :obj:`pyham.genome.ExtantGenome` that match the query name.

            Args:
                name (:obj:`str`): Name of the :obj:`pyham.genome.ExtantGenome`.

            Returns:
                :obj:`pyham.genome.ExtantGenome` or raise KeyError

        """

        for taxon in self.taxonomy.leaves:
            if taxon.name == name:
                if "genome" in taxon.features:
                    return taxon.genome

        raise KeyError('No extant genomes match the query name: {}'.format(name))

    # ___ AncestralGenome ___ #

    def get_list_ancestral_genomes(self):

        """  
            Get the list of all :obj:`pyham.genome.AncestralGenome` created during the parsing.

            Returns:
                a list of :obj:`pyham.genome.AncestralGenome`.

        """
        return [internal_node.genome for internal_node in self.taxonomy.internal_nodes]

    def get_ancestral_genome_by_taxon(self, taxon):

        """  
        Get the :obj:`pyham.genome.AncestralGenome` corresponding of the query taxon.

            Args:
                taxon (:obj:`str`): treeNode object of the :obj:`pyham.taxonomy.Taxonomy`.tree object.

            Returns:
                :obj:`pyham.genome.AncestralGenome` or raise KeyError

        """

        if taxon in self.taxonomy.internal_nodes and "genome" in taxon.features:
                return taxon.genome

        raise KeyError("Taxon {} doesn't have a genome attached to it.".format(taxon))

    def get_ancestral_genome_by_name(self, name):

        """  
        Get the :obj:`pyham.genome.AncestralGenome` corresponding of the query name.

            Args:
                name (:obj:`str`): Name of the :obj:`pyham.genome.AncestralGenome`.

            Returns:
                :obj:`pyham.genome.AncestralGenome` or raise KeyError

        """

        for taxon in self.taxonomy.internal_nodes:
            if taxon.name == name:
                if "genome" in taxon.features:
                    return taxon.genome

        raise KeyError('No ancestral genomes match the query name: {}'.format(name))

    def get_ancestral_genome_by_mrca_of_genome_set(self, genome_set):

        """  
        Get the :obj:`pyham.genome.AncestralGenome` corresponding to the MRCA of query genomes.

            Args:
                genome_set (:obj:`set`): Set of :obj:`pyham.genome.AncestralGenome`.

            Returns:
                :obj:`pyham.genome.AncestralGenome` or raise KeyError

        """

        if len(genome_set) < 2:
            raise ValueError('Minimum 2 genomes are required, only {} provided.'.format(len(genome_set)))

        for g in genome_set:
            if not isinstance(g, Genome):
                raise TypeError("expect subclass obj of '{}', got {}"
                                .format(Genome.__name__,
                                        type(g).__name__))

        genome_nodes = set([geno.taxon for geno in genome_set])

        mrca_node = self.taxonomy.tree.get_common_ancestor(genome_nodes)

        return self.get_ancestral_genome_by_taxon(mrca_node)

    # Taxon

    def get_taxon_by_name(self, name):

        """  
        Get the :obj:`ete3.TreeNode` object of the :obj:`pyham.taxonomy.Taxonomy`.tree corresponding of the query name.

            Args:
                name (:obj:`str`): Name of the treeNode.

            Returns:
                :obj:`ete3.TreeNode` or raise KeyError

        """

        nodes_founded = self.taxonomy.tree.search_nodes(name=name)

        if not nodes_founded:
            raise KeyError('No node founded for the species name: {}'.format(name))
        elif len(nodes_founded) == 1:
            return nodes_founded[0]
        else:
            raise KeyError('{} nodes founded for the species name: {}'.format(len(nodes_founded), name))

    # ... PRIVATE METHODS ... #

    def _add_missing_taxon(self, child_hog, oldest_hog, missing_taxons):

        """  
        Add intermediate :obj:`HOG` in between two :obj:`HOG` if their taxon are not direct parent and child in the 
        taxonomy. E.g. if a rodent HOG is connected with a vertebrate HOG it will add an mammal hog in between.

            Args:
                child_hog (:obj:`HOG`): child :obj:`HOG`.
                oldest_hog (:obj:`HOG`): parent :obj:`HOG`.
                missing_taxons (:obj:`HOG`): list of intermediate taxNode between child_hog and oldest_hog sorted 
                from youngest to oldest.

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

        # the youngest hog is removed from the oldest hog children.
        oldest_hog.remove_child(child_hog)

        # Then for each intermediate level in between the two hogs...
        current_child = child_hog
        for tax in missing_taxons:

            # ... we get the related ancestral genome of this level...
            ancestral_genome = self._get_ancestral_genome_by_taxon(tax)

            # ... we create the related hog and add it to the ancestral genome...
            hog = abstractgene.HOG()
            hog.set_genome(ancestral_genome)
            ancestral_genome.add_gene(hog)

            # ... we check if taxon correspond to child parent taxon ...
            if ancestral_genome.taxon is not current_child.genome.taxon.up:
                raise TypeError("HOG taxon {} is different than child parent taxon {}".format(ancestral_genome.taxon,
                                                                                              current_child.genome.taxon.up))

            # ... we add the child if everything is fine.
            hog.add_child(current_child)
            current_child = hog

        oldest_hog.add_child(current_child)

    def _get_oldest_from_genome_pair(self, g1, g2):

        """  
        Get the oldest :obj:`Genome` for a pair of :obj:`Genome`.

            Args:
                g1 (:obj:`Genome`): First :obj:`Genome`.
                g2 (:obj:`Genome`): Second :obj:`Genome`.

            Returns:
                :obj:`Genome`

        """

        mrca = self.taxonomy.tree.get_common_ancestor({g1.taxon,g2.taxon})

        if g1.taxon == mrca:
            return g1, g2
        elif g2.taxon == mrca:
            return g2, g1
        else:
            raise TypeError("The genomes are not in the same lineage: {}".format({g1, g2}))

    def _get_ancestor_and_descendant(self, genome_set):

        """  
        This method fetch from a set of :obj:`Genome`:
            - the oldest :obj:`Genome` from the set (if oldest genome not in set we get their mrca ).
            - the rest of the :obj:`Genome` present in the set.

            Args:
                genome_set (:obj:`set`): A set of :obj:`Genome`.

            Returns:
                :obj:`Genome`, a set of :obj:`Genome`.

        """

        ancestor = self._get_ancestral_genome_by_mrca_of_genome_set(genome_set)
        genome_set.discard(ancestor)
        return ancestor, genome_set

    def _get_HOGMap(self, genome_pair_set):

        """ 
        Get the :obj:`HOGMap` between two genomes.
        
            Args:
                genome_pair_set (:obj:`set`): A set of 2 :obj:`Genome`.

            Returns:
                :obj:`HOGMap`

        """

        f = frozenset(genome_pair_set)

        if f in self.HOGMaps.keys():
            return self.HOGMaps[f]
        else:
            self.HOGMaps[f] = mapper.HOGsMap(self, list(genome_pair_set)[0], list(genome_pair_set)[1])
            return self.HOGMaps[f]

    def _build_hogs_and_genes(self, file_object, filter_object):

        """ This function build from an orthoxml file all data that is required to build this Ham object (using the Ham
        filter object).

            Args:
                file_object (:obj:`FileObject`): File Object of the orthoxml to parse.
                filter_object (:obj:`ParserFilter`): :obj:`ParserFilter` use by OrthoXMLParser.

            Returns:
                :obj:`set` of top level :obj:`HOG` , :obj:`dict` of unique id with their :obj:`Gene`, :obj:`dict` of
                external id with their :obj:`Gene`.

        """

        factory = parsers.OrthoXMLParser(self, filterObject=filter_object)
        parser = XMLParser(target=factory)

        for line in file_object:
            parser.feed(line)

        return factory.toplevel_hogs, factory.extant_gene_map, factory.external_id_mapper

    def _get_extant_genome_by_name(self, **kwargs):

        """ 
        Get the :obj:`ExtantGenome` by name, if not founded in the taxonomy.tree.node.genome then created it.

            Args:
                **kwargs: dictionary of attribute and value required to create the :obj:`ExtantGenome`.

            Returns:
                :obj:`ExtantGenome`

        """

        nodes_founded = self.taxonomy.tree.search_nodes(name=kwargs['name'])

        if len(nodes_founded) == 1:

            node = nodes_founded[0]

            if "genome" in node.features:
                return node.genome

            else:
                extant_genome = ExtantGenome(**kwargs)
                self.taxonomy.add_genome_to_node(node, extant_genome)
                return extant_genome
        else:
            raise KeyError('{} node(s) founded for the species name: {}'.format(len(nodes_founded), kwargs['name']))

    def _get_ancestral_genome_by_taxon(self, tax_node):

        """  
        Get the :obj:`AncestralGenome` corresponding of the query taxon if not founded in the taxonomy.tree
        then created it.

            Args:
                tax_node : treeNode object of the :obj:`Taxonomy`.tree object.

            Returns:
                :obj:`AncestralGenome`

        """

        if "genome" in tax_node.features:
            return tax_node.genome

        else:
            ancestral_genome = AncestralGenome()
            self.taxonomy.add_genome_to_node(tax_node, ancestral_genome)

            return ancestral_genome

    def _get_ancestral_genome_by_mrca_of_hog_children_genomes(self, hog):

        """  
        Get MRCA :obj:`AncestralGenome` of the list of children :obj:`Genome` of the query :obj:`HOG`.
        
            Args:
                hog (:obj:`HOG`): query HOG.
        
            Returns:
                :obj:`AncestralGenome`
        
        """

        children_genomes = set([child.genome for child in hog.children ])

        return self._get_ancestral_genome_by_mrca_of_genome_set(children_genomes)

    def _get_ancestral_genome_by_mrca_of_genome_set(self, genome_set):

        """  
        Get the :obj:`AncestralGenome` corresponding to the MRCA of query genomes.

            Args:
                genome_set (:obj:`set`): Set of :obj:`AncestralGenome`.

            Returns:
                :obj:`AncestralGenome` or raise KeyError

        """

        if len(genome_set) < 2:
            raise ValueError('Minimum 2 genomes are required, only {} provided.'.format(len(genome_set)))

        for g in genome_set:
            if not isinstance(g, Genome):
                raise TypeError("expect subclass obj of '{}', got {}"
                                .format(Genome.__name__,
                                        type(g).__name__))

        genome_nodes = set([gen.taxon for gen in genome_set])

        mrca_node = self.taxonomy.tree.get_common_ancestor(genome_nodes)

        return self._get_ancestral_genome_by_taxon(mrca_node)