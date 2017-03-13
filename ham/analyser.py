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

    def get_genes_history_from_ancestor(self):
        '''
        get which genes came from duplication from parent or have been gains, etc..
        :return:
        '''
        pass


