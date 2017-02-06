from abc import ABCMeta, abstractmethod


class Analyser(metaclass=ABCMeta):

    """Analyser abstract class for HOGs/Genomes analysis.
    """

    def __init__(self, ham):
        self.ham = ham

class SingleHOGAnalyser(Analyser):

    """SingleHOGAnalyser class for analysis of a single hog (inherit from Analyser).

    Attributes:
        hog    abstract.HOG object of interest
    """

    def __init__(self, hog, hog_viz=True, stats=True, tree_profile=True):
        super(SingleHOGAnalyser, self).__init__()
        self.hog = hog

    def run_hog_viz(self):
        ## call module hog_viz
        pass

    def tree_profile(self):
        ## call module tree_profile
        pass

    def compute_statistics(self):
        ## call methods from HOG object ?
        pass

class AncestralGenomeAnalyser(Analyser):

    """AncestralGenomeAnalyser class for analysis of a single ancestral genome (inherit from Analyser).

    Attributes:
        ancestral_genome    Genome.ancestral_genome object of interest
    """

    def __init__(self, ancestral_genome, clustering=True, stats=True, filter=True):
        super(AncestralGenomeAnalyser, self).__init__()
        self.ancestral_genome = ancestral_genome

    def get_ancestral_genes_clustering(self):
        pass

    def compute_statistics(self):
        ## call methods from AncestralGenome object ?
        pass

    def filter_hogs(self):
        ## Filter all the hogs that fullfil a criterion (no dup, lost or gains at this levels)
        pass

class MultiGenomeAnalyser(Analyser):

    """MultiGenomeAnalyser class for analysis of a multiple ancestral genome (inherit from Analyser).

    Attributes:
        genomes    set of Genome object
    """

    def __init__(self, ancestral_genomes):
        super(MultiGenomeAnalyser, self).__init__()
        self.ancestral_genomes = ancestral_genomes


    def get_hogMap(self):
        ## get the mapping between all hogs
        pass

