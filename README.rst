pyHam: A Tool to Analyze Hierarchical Orthologous Groups (HOGs)
===============================================================


Motivation
----------
pyHam is a library to facilitate the analysis of hierarchical orthologous groups.
It also contains a few tools to run standard type of analysis.

Currently, Ham is limited to analyze HOGs that are stored in an orthoXML file.
More information on the schema of orthoxml and some examples are
available at http://orthoxml.org.

For extended documentation we refer to the docs folder that contain information
on common use cases and API documentation of the library.

How to cite pyHam
-----------------
If you use pyHam in your work, please consider citing:

*Clément-Marie Train, Miguel Pignatelli, Adrian Altenhoff, Christophe Dessimoz; iHam and pyHam: visualizing and processing hierarchical orthologous groups, Bioinformatics, bty994, https://doi.org/10.1093/bioinformatics/bty994*

Installation
------------
pyHam requires Python 3.9 or newer. Its dependencies (ete4, lxml, requests, tqdm and scipy)
are installed automatically by pip:

.. code-block:: sh

   python -m pip install --upgrade pip
   pip install pyham

To work on pyHam itself, install it in editable mode from a checkout of the git repository:

.. code-block:: sh

   python -m pip install -e ".[dev]"
   pytest

Quick start
-----------

.. code-block:: python

   import pyham

   ham = pyham.Ham("simpleEx.phyloxml", "simpleEx.orthoxml",
                   tree_format="phyloxml", use_internal_name=True)

   # look up a gene family (HOG) and export it as a gene tree
   hog = ham.get_hog_by_id("HOG:0355161")
   print(hog.to_ete_tree().write())

   # compare a species with one of its ancestors
   human = ham.get_extant_genome_by_name("Homo sapiens")
   primates = ham.get_ancestral_genome_by_name("Primates")
   vmap = ham.compare_genomes_vertically(human, primates)
   print(len(vmap.get_retained()), len(vmap.get_duplicated()),
         len(vmap.get_lost()), len(vmap.get_gained()))

Example
-------
The ``example`` folder contains a few ready-to-use scripts that show the main pyHam features.
Run them from inside that folder:

.. code-block:: sh

   cd example
   python run_hog_queries.py
   python run_treeProfile.py
   python run_iHam.py

Getting started
---------------
We create a small introductory blog post about HOGs and pyHam at http://lab.dessimoz.org/blog/2017/06/29/pyham. We highly recommend you to read it before starting using pyHam.

We also create an ipython notebook to help you with basic uses of pyHam API and embedded tools at http://zoo.cs.ucl.ac.uk/tutorials/tutorial_pyHam_get_started.html.

Table of compatibility
----------------------

Support for pyHam by various HOG inference resources.

+-----------------+------------------------------+---------------------------------------+-----------+
| Resource        | Species tree format          | OrthoXML                              | SUPPORT   |
+=================+==============================+=======================================+===========+
||OMA browser|_   | |PhyloXMLo|_ and |Newicko|_  ||All HOGso|_ , or |one HOG at a timeo|_|    YES    |
+-----------------+------------------------------+---------------------------------------+-----------+
||OMA standalone|_| PhyloXML and Newick          | All HOGs                              |    YES    |
+-----------------+------------------------------+---------------------------------------+-----------+
| |FastOMA|_      | Newick                       | All HOGs                              |    YES    |
+-----------------+------------------------------+---------------------------------------+-----------+
| |Ensembl|_      | |Newicke|_                   |    one HOG at a time                  |    YES    |
+-----------------+------------------------------+---------------------------------------+-----------+
| |HieranoidDB|_  | Newick                       |    one HOG at a time                  |    YES    |
+-----------------+------------------------------+---------------------------------------+-----------+

.. |OMA browser| replace:: ``OMA browser``
.. |OMA standalone| replace:: ``OMA standalone``
.. |FastOMA| replace:: ``FastOMA``
.. |Ensembl| replace:: ``Ensembl``
.. |HieranoidDB| replace:: ``HieranoidDB``

.. |PhyloXMLo| replace:: ``PhyloXML``
.. |Newicko| replace:: ``Newick``
.. |PhyloXMLs| replace:: ``PhyloXML``
.. |Newicks| replace:: ``Newick``
.. |Newicke| replace:: ``Newick``

.. |All HOGso| replace:: ``All HOGs``
.. |one HOG at a timeo| replace:: ``one HOG at a time``

.. _OMA browser: https://omabrowser.org
.. _OMA standalone: https://omabrowser.org/standalone/
.. _FastOMA: https://github.com/DessimozLab/FastOMA
.. _Ensembl: https://www.ensembl.org/index.html
.. _HieranoidDB: http://hieranoidb.sbc.su.se/

.. _PhyloXMLo: https://omabrowser.org/All/speciestree.phyloxml
.. _Newicko: https://omabrowser.org/All/speciestree.nwk
.. _Newicke: https://www.ensembl.org/info/about/speciestree.html

.. _All HOGso:  https://omabrowser.org/All/oma-hogs.orthoXML.gz
.. _one HOG at a timeo:  https://omabrowser.org/oma/hogs/


Common input problems
---------------------
Most errors when loading data come from a species tree that doesn't match the OrthoXML file.
To check your inputs before analysing them, run:

.. code-block:: sh

   pyham-validate species_tree.nwk hogs.orthoxml

or, from Python, ``print(pyham.validate_input(tree_file, hog_file, ...))``, which takes the same
arguments as ``pyham.Ham``. It reports every problem it finds at once, each with the affected species,
levels or families and, where it can tell, how to fix it. The most common problems are:

- **Use the species tree the HOGs were inferred with.** pyHam maps each HOG onto the species
  tree, so the tree must contain every species in the OrthoXML file under the same name, and
  its topology must agree with how the HOGs are nested. With OMA standalone and FastOMA, use the
  species tree written alongside the HOGs.
- **use_internal_name=True needs unique internal node names.** Some newick trees label internal
  nodes with numbers (e.g. support values) instead of taxon names. pyHam then reports that the
  internal names are not unique. Pass ``use_internal_name=False`` and pyHam will build names
  from the child species instead.
- **PhyloXML trees need a name field on every node.** Choose the field with
  ``phyloxml_leaf_name_tag`` and ``phyloxml_internal_name_tag``. If the tree has no
  ``<taxonomy>`` elements, use ``clade_name``, or pass ``use_internal_name=False``.
- **Taxonomic conflicts.** If a HOG's taxonomic range (``TaxRange`` or taxid) disagrees with
  the species tree, pyHam raises ``TaxonomicConflictError`` listing every conflicting family it
  found. Pass ``fail_fast=True`` to stop at the first one.
- **Singletons.** Genes that belong to no HOG can be looked up with ``get_gene_by_id``, but
  HOG-based tools such as tree profiles and iHam don't apply to them.

Documentation
-------------
The full documentation of pyHam is available at http://zoo.cs.ucl.ac.uk/doc/pyham/index.html.


