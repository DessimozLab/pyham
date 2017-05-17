.. _basics:
The Basics on pyham
=================

pyHam (“HAM: HOG Analysis Method”) is a python library to handle orthoXML containings Hierarchical Orthologous Groups (HOGs). The motivation of pyham is to create an easy python interface to investigate on HOGs and evolutionary information that can be induced from them. pyham provide graphical tools to visualize evolution history of single gene family or comparative genomic setup.



What are hogs ?
###############

Hierarchical Orthologous Groups (HOGs) are defined as set of genes that all descended from a single common ancestral genomes at a specific taxonomic range.

In term of labeled phylogenetic trees, an HOG at a given taxonomic range is corresponding to the subtree rooted by speciation nodes related to this taxon. Regarding the last definition, we can represent gene family evolution (labeled gene tree) by nesting HOGs based on the succession of speciation node and by grouping under paralogous groups HOGs that branch into duplication nodes.

HOGs represent different information depending on the perspective used to analyse them. As a single entity they correspond to gene families evolutionary history while several HOGs all at a specific taxonomic range represents the ancestral genes composition of the related ancestral genomes.

-----------



What is the orthoXML format ?
#############################

The orthoXML format is a XML based format design for phylogenetic purposes. The benefit of orthoXML is to have uniform manner to store orthology related information.

The orthoXML format is composed of two main parts:
    -   **Header**: This part contains information about the species/genes related to this genomic setup. Each species is defined by a (unique) name tag that serves as identifier by pyham to select extant genomes. Each species contain gene elements that are described by an unique orthoxml id tag that serves as identifier by Ham to select extant genes, and other possible id tags (protId, geneId,transcriptId) for cross referencing.
    -   **Groups**: This part contains nested groups representing the HOGs.
        Those nested groups can be divided into orthologGroup and paralogGroup depending on how the children elements of those groups had emerged, respectively here by speciation and duplication events.
        Both groups elements can contained geneRef elements (with orthoXML unique id tag mapping to gene in the header part) or others nested groups.
        The top level groups (corresponding to the root element of an HOG) are defined by an unique id tag.

.. seealso:: Here is the official documentation page http://orthoxml.org/0.3/orthoxml_doc_v0.3.html


-----------


How pyham worked under-the-hood ?
#################################

pyham require as minimal input an orthoxml that contains HOGs and a newick species tree.
**The leaves name in the newick tree should correspond 1:1 with the species name tag use in the header part of the orthoXML.**

The first step is to instantiate a :obj:`pyham.taxonomy.Taxonomy` object by constructing a ete3.Etree object based on the inputted newick tree and assign to each Etree node the related :obj:`pyham.genome.Genome`.

:obj:`pyham.genome.ExtantGenome` are named according to the newick leaves name while :obj:`pyham.genome.AncestralGenome` are either name using the internal node name from the newick tree or the the concatenation of the genome children name.

Then the seconds step is the orthoXML parsing using the :obj:`pyham.parser.OrthoXMLParser` that can be split in 2 steps:
    -   First the header part with the cross reference is parse. For each gene element within species elements, the related :obj:`pyham.abstractgene.ExtantGene` is build using the xref id (id, geneId,protId,transcriptId). **The "id" tag will be the unique identifier of the related genes accross the pyham analysis.**
    -   Second, the groups element part is parse. For each top level groups, we are creating the :obj:`pyham.abstractgene.HOG` for each internal orthologGroups founded insided and connecting them based on the nested groups hierarchy. **The** :obj:`pyham.abstractgene.HOG` **are connected to the related** the :obj:`pyham.genome.AncestralGenome` **they belong to based on the a mrca search of all the species its descendants**  :obj:`pyham.abstractgene.ExtantGene`.
        The hierarchical structure is build by saving the parent/children relations accross :obj:`pyham.abstractgene.HOG` and the :obj:`pyham.abstractgene.ExtantGene`. The paralogGroup information is stored by taging all HOGs emerging through duplication.

.. note:: Ham provide a way to restrict the orthoxml parsing to the information of interest in case of large orthoXML files.
            To proceed the :obj:`pyham.ham.FilterParser` can take as input which HOGs to proccess (based on a gene id, an hog id or a external id) and pre-select for the minimal required information to load for the :obj:`pyham.parser.OrthoXMLParser`.


**Glossary**:
    - **Top level HOG**:  root HOG that have no parent and is direct child of the groups element. This HOG have an unique top level id that act as unique identifier.
    - **Gene Unique id** (protId, geneId, transcriptId): The unique id is orthoXML id founded in the "id" tag in the header part. The others ids are for cross references and are not meant to be unique.
    - **Singleton**: Gene that is present in the orthoxml (in the header part) but belong to any HOGs.
-----------

What are the visualisation tool provide by pyham ?
##################################################

pyham provide two different phylogeny based visualisation tools called: Hogvis and Tree Profile. The goal of those visualisation tool is to synthesise concisely phylogenetic information using different perspective.

Hogvis : tool to visualise how the HOG members genes are clustering based on their ancestral genes membership.

TreeProfile: TreeProfile is a tool to visualise how the genes have evolved in terms of evolutionnary events alonge a phylogenetic tree (duplication, lost, gained).

**If you want to discover more about those visualisation tools, the pyham ipython notebook tutorial provides interactive explanations.**


-----------

