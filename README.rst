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
Ham is written in python3 (also compatible with python2), with little external dependencies, i.e.
currently ete3, lml, six. The setup script should resolve these
dependencies automatically.
Consider using pip to install the package directly from a checked out git repo

.. code-block:: sh

   python -m pip install --upgrade pip
   pip install pyham

Example
-------
We prepare a ready-to-use example (see example folder) with few python scripts to use main pyHam features.
You just have to run the following command in bash:

.. code-block:: sh

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
| |Ensembl|_      | |Newicke|_                   |    one HOG at a time                  |    YES    |
+-----------------+------------------------------+---------------------------------------+-----------+
| |HieranoidDB|_  | |Newickh|_                   |    one HOG at a time                  |    YES    |
+-----------------+------------------------------+---------------------------------------+-----------+

.. |OMA browser| replace:: ``OMA browser``
.. |OMA standalone| replace:: ``OMA standalone``
.. |Ensembl| replace:: ``Ensembl``
.. |HieranoidDB| replace:: ``HieranoidDB``

.. |PhyloXMLo| replace:: ``PhyloXML``
.. |Newicko| replace:: ``Newick``
.. |PhyloXMLs| replace:: ``PhyloXML``
.. |Newicks| replace:: ``Newick``
.. |Newicke| replace:: ``Newick``
.. |Newickh| replace:: ``Newick``

.. |All HOGso| replace:: ``All HOGs``
.. |one HOG at a timeo| replace:: ``one HOG at a time``

.. _OMA browser: https://omabrowser.org
.. _OMA standalone: https://omabrowser.org/standalone/
.. _Ensembl: https://www.ensembl.org/index.html
.. _HieranoidDB: http://hieranoidb.sbc.su.se/

.. _PhyloXMLo: https://omabrowser.org/All/speciestree.phyloxml
.. _Newicko: https://omabrowser.org/All/speciestree.nwk
.. _Newicke: https://www.ensembl.org/info/about/speciestree.html
.. _Newickh: http://hieranoid.sbc.su.se/download/H2/66c.tree

.. _All HOGso:  https://omabrowser.org/All/oma-hogs.orthoXML.gz
.. _one HOG at a timeo:  https://omabrowser.org/oma/hogs/


Documentation
-------------
You can the full documentation of pyHam at http://zoo.cs.ucl.ac.uk/doc/pyham/index.html


