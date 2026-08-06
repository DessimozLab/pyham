from . import abstractgene
from .genome import ExtantGenome
from .taxonomy import build_taxon_node, Taxonomy
import logging
logger = logging.getLogger(__name__)
from collections import defaultdict
from tqdm.auto import tqdm
import numpy as np


class OrthoXMLParser:

    """
    Custom parser for OrthoXML containing HOGs. It can take a FilterParser object to restrict the information to parse.
     
    The parse goes through the whole XML and create on the fly the required Ham objects:
        - In the Xref/header, the parser creates the ExtantGenome and Gene objects.
        - the taxonomy group and its taxon elements are parsed to create the taxonomy tree, if no species tree is provided.
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

    def __init__(self, ham_object, taxonomy=None, filterObject=None, id_schema: str = None, with_progress: bool = False):

        """
        Args:
            ham_object (:obj:`Ham`): Ham object to feed with created objects.
            filterObject (:obj:`FilterParser`, optional): FilterParser object used to restrict the parsed information.
            with_progress (:bool:, optional): Whether to display a tqdm progress bar whilst parsing.
            Defaults to None.
        """
        self.ham_object = ham_object
        self.filterObj = filterObject
        self.taxonomy = ham_object.taxonomy

        # usefull information
        self.extant_gene_map = {}
        self.external_id_mapper = defaultdict(list)
        self.toplevel_hogs = {}

        # On the fly variable
        self.cpt = 0
        self.hog_stack = []
        self.paralog_stack = []
        self.taxon_stack = []
        self.current_species = None
        self.in_paralogGroup = None
        self.paralogyNode = None
        self.skip_this_hog = False
        self._genomes_to_add = []

        # save progress bar option
        self.with_progress = with_progress

    def _build_gene(self, attrib):
        gene = abstractgene.Gene(**attrib)
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

    def _update_progress_bar(self, tag):
        if tag == "{http://orthoXML.org/2011/}species":
            if not hasattr(self, 'sp_pbar'):
                self.sp_pbar = tqdm(desc='Parsing Species')
            self.sp_pbar.update()
        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            if not hasattr(self, 'hog_pbar'):
                if hasattr(self, 'sp_pbar'):
                    self.sp_pbar.close()
                    logger.info("Parsing Species completed.")
                    delattr(self, 'sp_pbar')
                self.hog_pbar = tqdm(desc='Parsing HOGs')
            self.hog_pbar.update()
        elif tag == "{http://orthoXML.org/2011/}groups":
            if hasattr(self, 'hog_pbar'):
                self.hog_pbar.close()
                delattr(self, 'hog_pbar')
            logger.info("Parsing completed.")

    def _assign_ancestral_genome(self, hog):
        """This method determines the ancestral genome for a HOG based on its children."""
        child_genomes = [child.genome for child in hog.children]
        tax_node = None
        if hog.taxon_id is not None:
            # If the HOG has a taxon_id, we can directly assign the ancestral genome based on the taxon.
            tax_nodes = list(self.taxonomy.nodes_by_attr(id=hog.taxon_id))
            if len(tax_nodes) == 1:
                tax_node = tax_nodes[0]
                if len(set(child_genomes)) == 1 and child_genomes[0].taxon == tax_node:
                    raise abstractgene.SameLevelHOGError(f"HOG '{hog.hog_id}' at {tax_node.name} has child at same level: {hog.children[0]}.")
        # let's check if we can map it using the TaxRange or taxid property
        if tax_node is None:
            for key in ('TaxRange', 'taxid'):
                try:
                    tax_value = hog[key]
                    tax_node = self.taxonomy.get_node_by_name(tax_value)
                except (KeyError, AttributeError):
                    pass
                else:
                    try:
                        for child_genome in child_genomes:
                            if not self.taxonomy.is_child_recursive(child_genome.taxon, tax_node):
                                raise abstractgene.EvolutionaryConceptError(f"HOG '{hog.hog_id}' at {tax_node.name} has child not below: {child_genome}.")
                        break
                    except abstractgene.EvolutionaryConceptError:
                        if len(set(child_genomes)) == 1 and child_genome.taxon == tax_node:
                            raise abstractgene.SameLevelHOGError(f"HOG '{hog.hog_id}' at {tax_node.name} has child at same level: {child_genome}.")
                        raise

        # if not, let's try with the children genomes
        if tax_node is None:
            if len(set(child_genomes)) == 1:
                try:
                    tax_range = hog['TaxRange']
                    if hog['TaxRange'] == hog.children[0].genome.name:
                        raise abstractgene.SameLevelHOGError(f"HOG '{hog.hog_id}' at {tax_range} has child at same level: {hog.children[0]}.")
                except (KeyError, AttributeError):
                    pass
                tax_node = child_genomes[0].taxon.up
            else:
                tax_node = self.taxonomy.get_mrca_taxnode(*(g.taxon for g in child_genomes))
        genome = self.taxonomy.get_genome_from_taxnode(tax_node)
        genome.add_gene(hog)
        logger.debug(f"Added {hog} to ancestral genome {genome.name}")

    def _create_missing_hogs(self, child: abstractgene.AbstractGene, parent: abstractgene.HOG):
        """Create missing HOGs between a child and a parent genome.
        Returns the oldest created HOG, or the child if no HOG was created."""
        if child == parent:
            raise ValueError("Child and parent genomes are the same, cannot create missing HOGs.")

        missing_taxons = self.taxonomy.get_path_up(child.genome.taxon, parent.genome.taxon)
        if len(missing_taxons) > 0:
            # the youngest hog is removed from the oldest hog's children.
            parent.remove_child(child)
            # Then for each intermediate level in between the two hogs...
            current_child = child
            hog_id = child.hog_id if hasattr(child, 'hog_id') else parent.hog_id
            for tax in missing_taxons:
                # ... we get the related ancestral genome of this level...
                ancestral_genome = self.taxonomy.get_genome_from_taxnode(tax)

                # ... we create the related hog and add it to the ancestral genome...
                hog = abstractgene.HOG(id=hog_id)
                setattr(hog, '_missing_in_xml', parent.genome)
                ancestral_genome.add_gene(hog)

                # ... we check if taxon correspond to child parent taxon ...
                if ancestral_genome.taxon is not current_child.genome.taxon.up:
                    raise TypeError(
                        "HOG taxon {} is different than child parent taxon {}".format(ancestral_genome.taxon,
                                                                                      current_child.genome.taxon.up))
                # ... we add the child if everything is fine.
                hog.add_child(current_child)
                current_child = hog
            parent.add_child(current_child)
        logger.debug(f"Created {len(missing_taxons)} missing HOGs between {child} and {parent}.")
        return current_child if len(missing_taxons) > 0 else child

    def _resolve_missing_levels(self, parent: abstractgene.HOG, children: list):
        """Insert missing intermediate HOG levels between `parent` and `children`.

        If several children are missing the same intermediate taxon (i.e. the
        reference species tree resolves this part of the taxonomy more finely
        than the HOGs were built/augmented at), a single shared ancestral HOG is
        created for all of them instead of one redundant HOG per child. Children
        that are missing no shared level are delegated to `_create_missing_hogs`.
        """
        change = {}
        for child in children:
            if parent.genome.taxon.props["depth"] != child.genome.taxon.props["depth"] - 1:
                change[child] = self.taxonomy.get_path_up(child.genome.taxon, parent.genome.taxon)
        if not change:
            return

        groups = defaultdict(list)
        for child, missing in change.items():
            # missing[-1] is the taxon immediately below `parent`; children sharing
            # any ancestor on the path back to `parent` necessarily share this one too.
            groups[missing[-1]].append(child)

        for tax_node, group_children in groups.items():
            if len(group_children) == 1:
                self._create_missing_hogs(group_children[0], parent)
            else:
                shared_genome = self.taxonomy.get_genome_from_taxnode(tax_node)
                shared_hog = abstractgene.HOG(id=parent.hog_id)
                setattr(shared_hog, "_missing_in_xml", parent.genome)
                shared_genome.add_gene(shared_hog)
                parent.add_child(shared_hog)
                for child in group_children:
                    parent.remove_child(child)
                    shared_hog.add_child(child)
                self._resolve_missing_levels(shared_hog, group_children)

    def start(self, tag, attrib):
        if tag == "{http://orthoXML.org/2011/}species":
            self.current_species = ExtantGenome(**attrib)
            if self.taxonomy is not None:
                taxon = self.taxonomy.get_extant_taxa_by_name(self.current_species.name)
                self.taxonomy.add_genome_to_node(taxon, self.current_species)
            else:
                self._genomes_to_add.append(self.current_species)

        elif tag == "{http://orthoXML.org/2011/}gene" and (self.filterObj is None or attrib["id"] in self.filterObj.geneUniqueId):
            self._build_gene(attrib)

        elif tag == "{http://orthoXML.org/2011/}paralogGroup" and self.skip_this_hog is False:
            cur_depth = len(self.hog_stack)
            if len(self.paralog_stack) > 0 and self.paralog_stack[-1]['depth'] == cur_depth:
                # we have a directly nested paralog group. use the same DuplicationNode
                dNode = self.paralog_stack[-1]['node']
            else:
                dNode = abstractgene.DuplicationNode(self, id=attrib.get('og', None))

            self.paralog_stack.append({'depth': cur_depth, 'node': dNode})

            self.in_paralogGroup = self.paralog_stack[-1]['depth']
            self.paralogyNode = self.paralog_stack[-1]['node']

        elif tag == "{http://orthoXML.org/2011/}geneRef" and self.skip_this_hog is False:

            gene = self.extant_gene_map[attrib['id']]
            if 'LOFT' in attrib:
                gene.set_LOFT(attrib['LOFT'])
            self.hog_stack[-1].add_child(gene)

            # if the gene is contained within a paralogousGroup need to update its .arose_by_duplication flag.
            if self.in_paralogGroup == len(self.hog_stack):
                self.paralogyNode.add_child(gene)

        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            if len(self.hog_stack) == 0:
                # a new roothog
                if self.with_progress:
                    self._update_progress_bar(tag)

                if self.filterObj is not None and not attrib.get("id", None) in self.filterObj.hogsId:
                    self.skip_this_hog = True
            # create HOG and add it to the stack (or a dummy value if we skip this hog)
            if self.skip_this_hog:
                self.hog_stack.append(0)
            else:
                self._build_hog(attrib)

        elif tag == "{http://orthoXML.org/2011/}property" and not self.skip_this_hog:
            self.hog_stack[-1].add_property(attrib["name"], attrib["value"])

        elif tag == "{http://orthoXML.org/2011/}score" and not self.skip_this_hog:
            self.hog_stack[-1].score(attrib['id'], float(attrib['value']))

        elif tag == "{http://orthoXML.org/2011/}taxon":
            taxon = build_taxon_node(**attrib)
            self.taxon_stack.append(taxon)
            if len(self.taxon_stack) > 1:
                self.taxon_stack[-2].add_child(child=taxon)

        elif tag == "{http://orthoXML.org/2011/}groups":
            if self.taxonomy is None:
                raise ValueError("Cannot parse this OrthoXML file without an explicit taxonomy tree. Please provide one.")

    def end(self, tag):

        if tag == "{http://orthoXML.org/2011/}species":
            if self.with_progress:
                self._update_progress_bar(tag)
            logger.info(f"Species {self.current_species.name} created.")
            self.current_species = None

        elif tag == "{http://orthoXML.org/2011/}taxon":
            taxonomy = self.taxon_stack.pop()
            if len(self.taxon_stack) == 0:
                # we are at the root of the taxonomy tree. if no taxonomy is provided upfront, we create it.
                if self.taxonomy is None:
                    self.taxonomy = Taxonomy(taxonomy)
                    # add all extant genomes to the taxonomy tree
                    for genome in self._genomes_to_add:
                        taxon = self.taxonomy.get_extant_taxa_by_name(genome.name)
                        self.taxonomy.add_genome_to_node(taxon, genome)
                    logger.info("Taxonomy tree created from OrthoXML taxonomy. All extant genomes added to the tree.")
                else:
                    logger.info("OrthoXML taxonomy not used, using provided taxonomy tree instead.")

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
                return

            try:
                # assign the ancestral genome to this hog
                self._assign_ancestral_genome(hog)
            except abstractgene.SameLevelHOGError as e:
                # this hog has children at the same level as its parent. Usually an 'augmented' HOG.
                logger.info(str(e))
                # remove this hog structure.
                for child in hog.children:
                    child.parent = hog.parent
                hog.parent.children.remove(hog)
                hog.parent.children.extend(hog.children)
                if hog.arose_by_duplication:
                    dupl_node = hog.arose_by_duplication
                    dupl_node.remove_child(hog)
                    for child in hog.children:
                        dupl_node.add_child(child)
                return

            # get all child clustered by dup if any
            child_by_duplication = defaultdict(list)
            for child in hog.children:
                if child.arose_by_duplication:
                    child_by_duplication[child.arose_by_duplication].append(child)

            # assert that all duplications are have MRCA equal to hog.genome or below
            for duplication in child_by_duplication:
                assert duplication.MRCA == hog.genome or self.taxonomy.is_child_recursive(duplication.MRCA.taxon, hog.genome.taxon), f"Duplication {duplication.MRCA} is not below {hog.genome}"


            # For each duplication
            for duplication, children in child_by_duplication.items():

                # add hog at the level of the duplication MRCA if its missing
                if duplication.MRCA != hog.genome and self.taxonomy.is_child_recursive(duplication.MRCA.taxon, hog.genome.taxon):
                    # create the MRCA hog
                    mrca_genome = duplication.MRCA
                    mrca_hog = abstractgene.HOG(id=hog.hog_id)
                    setattr(mrca_hog, "_missing_in_xml", hog.genome)
                    mrca_genome.add_gene(mrca_hog)

                    # link it to the current hog as parent
                    hog.add_child(mrca_hog)

                    # add missing level down (from mrcaHOG to duplicated child)
                    for child in children:
                        # move child from hog to mrca_hog
                        hog.remove_child(child)
                        mrca_hog.add_child(child)

                    # add missing levels if any, grouping children that are
                    # missing the same intermediate ancestor under one shared HOG
                    self._resolve_missing_levels(mrca_hog, children)

                    duplication.set_parent(mrca_hog)
                    for x in list(duplication.children):
                        duplication.remove_child(x)
                    for y in list(mrca_hog.children):
                        duplication.add_child(y)

                # Otherwise simply add missing taxa between hog and duplicated child
                else:
                    duplication.MRCA = self.taxonomy.get_genome_from_taxnode(hog.genome.taxon)
                    duplication.set_parent(hog)
                    for child_direct in children:
                        new_direct_child = self._create_missing_hogs(child_direct, hog)
                        duplication.remove_child(child_direct)
                        duplication.add_child(new_direct_child)

            # insert any remaining missing intermediate levels, grouping children
            # that are missing the same intermediate ancestor under one shared HOG
            self._resolve_missing_levels(hog, list(hog.children))

            if len(self.hog_stack) == 0:
                hog_id = hog.hog_id.split('_')[0]
                self.toplevel_hogs[hog_id] = hog
                self.cpt += 1
                if self.cpt % 500 == 0:
                    logger.info("{} HOGs parsed. ".format(self.cpt))

        elif tag == "{http://orthoXML.org/2011/}groups":
            if self.with_progress:
                self._update_progress_bar(tag)

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




class PhyloXMLToETE:
    def __init__(self, node_factory=None):
        self.root = None
        self.stack = []
        self.current_tag = None
        self.buffer = ''
        self.in_taxonomy = False
        self.taxonomy_data = None
        self.in_clade = False
        from .taxonomy import Tree
        self.node_factory = node_factory if node_factory else Tree

    def start(self, tag, attrib):
        tag = self._strip_ns(tag)
        self.current_tag = tag

        if tag == "clade":
            dist = attrib.get("branch_length", None)
            node = self.node_factory({"dist": dist})
            if len(self.stack) > 0:
                self.stack[-1].add_child(node)
            else:
                self.root = node
            self.stack.append(node)

        elif tag == "taxonomy":
            self.in_taxonomy = True
            self.taxonomy_data = {}

    def end(self, tag):
        tag = self._strip_ns(tag)

        if tag == "name" and not self.in_taxonomy:
            if len(self.stack) == 0:
                self.buffer = ""
                return
            # name of clade. -> set the name of the last node in the stack
            self.stack[-1].name = self.buffer.strip()

        elif self.in_taxonomy:
            if tag in ('scientific_name', 'id', 'code'):
                self.taxonomy_data[tag] = self.buffer.strip()

            elif tag == "taxonomy":
                node = self.stack[-1]
                sci_name = self.taxonomy_data.get('scientific_name')
                if sci_name:
                    node.add_prop('scientific_name', sci_name)
                taxon_id = self.taxonomy_data.get('id')
                if taxon_id:
                    node.add_prop('taxon_id', taxon_id)
                code = self.taxonomy_data.get('code')
                if code:
                    node.add_prop('code', code)
                self.in_taxonomy = False
                self.taxonomy_data = None

        elif tag == "clade":
            self.stack.pop()

        self.buffer = ""
        self.current_tag = None

    def data(self, data):
        if self.current_tag:
            self.buffer += data

    def close(self):
        return self.root

    def _strip_ns(self, tag):
        return tag.split('}', 1)[-1] if '}' in tag else tag
