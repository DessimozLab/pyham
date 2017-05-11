HAM: A Tool to Analyze Hierarchical Orthologous Groups (HOGs)
=============================================================


Motivation 
----------
HAM is a library to facilitate the analysis of hierarchical orthologous groups.
It also contains a few tools to run standard type of analysis.

Currently, HAM is limited to analyze HOGs that are stored in an orthoXML file.
More information on the schema of orthoxml and some examples are
available at `http://orthoxml.org`_.

For extended documentation we refer to the docs folder that contain information
on common use cases and API documentation of the library.


Installation
------------
HAM is written in python3, with little external dependencies, i.e.
currently ete3, scipy, six. The setup script should resolve these
dependencies automatically. 
Consider using pip to install the package directly from a checked out git repo

.. code-block:: sh

   python -m pip install --upgrade pip
   pip install -e </path/to/HAM-repo/>

