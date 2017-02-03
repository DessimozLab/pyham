import logging

logger = logging.getLogger(__name__)


class HOGsMap(object):
    '''
    First, the oldest genome (if two genomes in lineage) or the MRCA (if > two genomes) is determined from the given genome_set.
    Once the oldest genome Gi is found, for each "child" genome Gj a "UpMap" is build to find the ancestor Hi of each Hj.
    All Hj in Gj that have no corresponding ancestor Hi in Gi are considered as GAIN.
    Then the "downMap" is reconstructed by regrouping all the Hj bound to a same Hi together.
    '''

    def __init__(self, ham_object, genome_set):
        self.HAM = ham_object
        self.ancestor, self.descendants = self.set_ancestor_and_descendants(genome_set)
        self.upMaps = self.computeUpMaps()
        self.downMap = self.computeDownMap()

    def set_ancestor_and_descendants(self, genome_set):
        ancestor = self.HAM._get_mrca_ancestral_genome_from_genome_set(genome_set)
        descendants = {}  # key|genome -> value|genome
        genome_set.discard(ancestor)
        for genome in genome_set:
            descendants[genome] = genome
        return ancestor, descendants

    def computeDownMap(self):
        downMap = {}  # key|HOGi -> value|{key|genome -> value|HOGj or "None" }
        for hog_ancestor in self.ancestor.genes:
            downMap[hog_ancestor] = {}
            for descandant in self.descendants.values():
                downMap[hog_ancestor][descandant] = [None]
        return downMap

    def computeUpMaps(self):
        upMaps = {}  # key|genome -> value|{key|HOGj -> value|HOGi or "None"}
        for descendant in self.descendants:
            upMaps[descendant] = self.build_UpMap(descendant)
        return upMaps

    def build_UpMap(self, descendant):
        upMap = {}  # value|{key|HOGj -> value|HOGi or "None"
        for hog_source in descendant.genes:
            if hog_source.parent is not None:  # avoid singletons
                hog_target = self.search_ancestor_hog_in_ancestral_genome(hog_source, self.ancestor)
                upMap[hog_source] = hog_target
        return upMap

    def search_ancestor_hog_in_ancestral_genome(self, source_hog, ancestral_genome):
        found = None
        current_hog = source_hog
        while current_hog.parent is not None and found is None:
            current_hog = current_hog.parent
            if current_hog.genome == ancestral_genome:
                found = current_hog
        return found

    def getLost(self):
        pass

    def getGained(self):
        pass

    def conserved(self): # single_copy + duplicated ?
        pass

    def singleCopy(self):
        pass

    def duplicated(self):
        pass


