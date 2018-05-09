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


Installation
------------
Ham is written in python3 (also compatible with python2), with little external dependencies, i.e.
currently ete3, lml, six. The setup script should resolve these
dependencies automatically.
Consider using pip to install the package directly from a checked out git repo

.. code-block:: sh

   python -m pip install --upgrade pip
   pip install pyham

GETTING STARTED
---------------
We create a small introductory blog post about HOGs and pyHam at http://lab.dessimoz.org/blog/2017/06/29/pyham. We highly recommend you to read it before starting using pyHam.

We also create an ipython notebook to help you with basic uses of pyHam API and embedded tools at http://zoo.cs.ucl.ac.uk/tutorials/tutorial_pyHam_get_started.html.


DOCUMENTATION
-------------
You can the full documentation of pyHam at http://zoo.cs.ucl.ac.uk/doc/pyham/index.html


