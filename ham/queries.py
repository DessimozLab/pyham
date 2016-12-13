def get_extant_genomes(taxonomy):
    """
    return all Genome.ExtantGenome attach to a leaf within a taxonomy object
    :param taxonomy: taxonomy.Taxonomy
    :return: set of Genome.ExtantGenome
    """
    return set(leaf.genome for leaf in taxonomy.leaves)


def get_ancestral_genomes(taxonomy):
    """
    return all Genome.AncestralGenome attach to an internal node within a taxonomy object
    :param taxonomy: taxonomy.Taxonomy
    :return: set of Genome.AncestralGenome
    """
    return set(internal_node.genome for internal_node in taxonomy.internal_nodes)


def get_top_level_hogs(hogs):
    """
    return all the top level rooted hogs from set of hogs
    :param hogs: set of abstract.HOG
    :return: set of top level rooted hogs
    """
    return set(hog for hog in hogs if hog.depth == 1)

