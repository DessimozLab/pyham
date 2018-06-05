.. _basics:
The Basics on pyham
===================

pyHam (“python HOG analysis method”) is a python library to handle orthoXML containings Hierarchical Orthologous Groups (HOGs). The motivation of pyham is to create an easy python interface to investigate on HOGs and evolutionary information that can be induced from them. pyham provide graphical tools to visualize evolution history of single gene family or comparative genomic setup.



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

What are the visualisation tool provide by pyham ?
##################################################

pyham provide two different phylogeny based visualisation tools called: IHAM and Tree Profile. The goal of those visualisation tool is to synthesise concisely phylogenetic information using different perspective.

IHAM : tool to visualise how the HOG members genes are clustering based on their ancestral genes membership.

TreeProfile: TreeProfile is a tool to visualise how the genes have evolved in terms of evolutionnary events along a phylogenetic tree (duplication, lost, gained).

**If you want to discover more about those visualisation tools, the pyham ipython notebook tutorial provides interactive explanations at** https://zoo.cs.ucl.ac.uk/tutorials/tutorial_pyHam_get_started.html .


-----------

