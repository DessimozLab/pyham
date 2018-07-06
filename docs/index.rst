.. Ham documentation master file, created by
   sphinx-quickstart on Mon Dec 12 17:29:31 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to **pyHam**'s documentation !
======================================

.. note:: Documentation for pyham version |release|.

**pyHam** (“**py**\thon **H**\OG **a**\nalysis **m**\ethod”) is a python library for programmatic processing of **H**\ierarchical **O**\rthologous **G**\roups (HOGs) encoded as OrthoXML. The motivation of pyHam is to create an easy python interface to investigate on HOGs (gene families) and the evolutionary information that can be induced from them. pyHam provides interactive visualisation tools to investigate the evolutionary history of single gene families or a whole taxon.

.. toctree::
   :maxdepth: 1

   index
   api

---------------


How to install pyham ?
======================

.. note:: We encourage pyham users to work on virtual environment. If you don't understand the last sentences you should consider reading the following short guide to python virtual environment: http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/

Installing via pip
##################

pyHam can be install via pip using the following command:


    - **First upgrade pip to the latest version**:
        .. code-block:: bash

                python -m pip install --upgrade pip

    - **Then install the pyHam package**:

        .. code-block:: bash

                pip install pyham

        in case you don't have root access (and/or permissions for writing in system directories) you can use the --user flag to install in the local user package:

        .. code-block:: bash

                pip install --user pyham.

Source code
###########


The lastest source code of pyham can be download at https://github.com/DessimozLab/pyham .

---------------

The Basics of pyHam
===================

What are hogs ?
###############

Hierarchical Orthologous Groups (HOGs) are defined as set of genes that all descended from a single common ancestral genomes at a specific taxonomic range.

If you want to know more about HOGs, we made for you a video that explaine what they really are https://www.youtube.com/embed/5p5x5gxzhZA


What is the orthoXML format ?
#############################

The orthoXML format is a XML based format designed to encode gene trees. If you want to know more about orthoXML, we
invite you to visit the official website http://orthoxml.org/0.3/orthoxml_doc_v0.3.html


How to use pyham on my dataset ?
################################

In order to use pyHam you need:
 - A species tree of the phylogeny of interest (Supported format: newick, internal node names are optional, polytomies are accepted)
 - HOGs encoded in an OrthoXML file.

**IMPORTANT**\, species names need to be matching between the species tree and the orthoXML.

We invite you to consult the tutorial at https://zoo.cs.ucl.ac.uk/tutorials/tutorial_pyHam_get_started.html we made to help you set-up your first pyHam instance and also to show the panel
of tools and options pyHam offer


What are the visualisation tool provide by pyham ?
##################################################

pyHam provide two visualisation tools called: iHam and Tree Profile. The goal of those visualisation tools is to synthesise concisely phylogenetic information using different perspective.

iHam: tool to visualise how the HOG members genes are clustered based on their ancestral genes membership.

TreeProfile: TreeProfile can be used to display the proportion of genes that have been retained, lost, duplicated or newly appeared (“gained”) along each branch of the species tree. This can be run for either a
single HOG or a whole taxa.

