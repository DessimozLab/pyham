
__author__ = 'admin'

def get_newick_string(source_path, type="nwk"):
    '''
    Get the newick tree as string from a newick file, a hdf5 file or from the note tag of an orthoXML
    :param newick_source: path of the file
    :param type: type of file where to fetch the newick tree. Available: "nwk", "h5" and "orthoXML"
    :return:
    '''

    if type == "nwk":
        with open(source_path, 'r') as nwk_file:
            return nwk_file.read()

    elif type == "h5":
        pass

    elif type == 'orthoXML':
        pass






