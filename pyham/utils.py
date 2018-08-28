from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import open
from future import standard_library
standard_library.install_aliases()
from . import taxonomy as tax


def get_newick_string(source_path, type="nwk"):
    '''
    Get the newick tree as string from a newick file, a hdf5 file or from the note tag of an orthoXML
    Args:
        | newick_source: path of the file
        | type: type of file where to fetch the newick tree. Available: "nwk", "h5" and "orthoXML"
    '''

    if type == "nwk":
        with open(source_path, 'r') as nwk_file:
            return nwk_file.read()

    elif type == "phylip":
        pass

    elif type == "h5":
        pass

    elif type == 'orthoXML':
        pass


def previsualize_taxonomy(newick_str):
    """
    This function help to previsualyse before running Ham what the topology will looked like and what will be the
    internal node naming.
    Args:
        | tree_str:
    """

    t = tax.Taxonomy(newick_str)

    for node in t.tree.traverse("postorder"):
        t.set_taxon_name(node)

    return t.tree.get_ascii()







