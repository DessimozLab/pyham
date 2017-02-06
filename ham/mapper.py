import logging

logger = logging.getLogger(__name__)


class HOGsMap(object):
    '''
    Given a set of genomes Gj (and their MRCA Gi --if genomes are on a lineage, the oldest is selected else the MRCA
    is computed from the Taxonomy--), the HOGsMap object HOGs/Genes identify the relations that exist between all Hi
    and Hj respectively the HOGs fron Gi and Gj.

    Those relations are classified into the following categories:
        - SINGLE: if Hi will correspond to a single Hj
        - DUPLICATE: If a duplication occurred in the path between Hi and Hj
        - LOSS: [strict] If Hi is not bound to any Hj, the case where duplication followed by gene loss are not counted
         here # TODO add a distinction between "fully loss" and "partially loss" ?
        - GAIN: if Hj have no corresponding Hi

    First, the oldest genome (if two genomes in lineage) or the MRCA (if > two genomes) is determined from the given genome_set.
    Once the oldest genome Gi is found, for each "child" genome Gj a "UpMap" is build to find the ancestor Hi of each Hj.
    If during the traversal Hj to Hi a duplication occurred, a flag "paralog" is added in the UpMap[Hi|Hj].
    All Hj in Gj that have no corresponding ancestor Hi in Gi are considered as GAIN.

    Then we built the Event Clusters:
        - .LOSS: Hi -> [Gj] | fully lost if duplicated then only one copy lost this will not be count as lost.
        - .GAIN: Gj -> [Hj]
        - .SINGLE: Hi -> {Gj -> Hj}
        - .DUPLICATE: Hi -> {Gj -> [Hj]}

    '''

    def __init__(self, ham_object, genome_set):
        self.HAM = ham_object
        self.ancestor, self.descendants = self.set_ancestor_and_descendants(genome_set)
        self.upMaps = self.computeUpMaps()
        #self.downMap = self.computeDownMap()
        self.LOSS, self.GAIN, self.SINGLE, self.DUPLICATE = self.buildEventClusters()

    def set_ancestor_and_descendants(self, genome_set):
        ancestor = self.HAM._get_mrca_ancestral_genome_from_genome_set(genome_set)
        descendants = {}  # key|genome -> value|genome
        genome_set.discard(ancestor)
        for genome in genome_set:
            descendants[genome] = genome
        return ancestor, descendants

    def computeDownMap(self):  #  Todo Remove if Adrian think buildEventClusters is better
        downMap = {}  # key|HOGi -> value|{key|genome -> value|[HOGj] or [] }

        # build empty shell for downMap
        for hog_ancestor in self.ancestor.genes:
            downMap[hog_ancestor] = {}
            for g in self.descendants.values():
                downMap[hog_ancestor][g] = []

        # feed the downMap
        for descendant, upMap in self.upMaps.items():
            for hog_descendant, hog_target in upMap.items():
                if hog_target[0] is not None:
                    downMap[hog_target[0]][descendant].append(hog_descendant)
        return downMap

    def buildEventClusters(self):
        LOSS = {}  # {Hi -> [Gj]}
        GAIN = {g:[] for g in self.descendants}  # {Gj -> [Hj]}
        SINGLE = {}  # {Hi -> {Gj -> Hj}}
        DUPLICATE = {}  # {Hi -> {Gj -> [Hj]}}

        ancestral_hogs_computed = {g:set() for g in self.descendants}

        # feed the MAPS
        for descendant, upMap in self.upMaps.items():
            for hog_descendant, hog_target in upMap.items():
                if hog_target[0] is None: # If gained
                    GAIN[descendant].append(hog_descendant)
                else:  # If not gained
                    if hog_target[1] is True:  # is duplicated
                        DUPLICATE.setdefault(hog_target[0],{}).setdefault(descendant,[]).append(hog_descendant)
                        ancestral_hogs_computed[descendant].add(hog_target[0])
                    else:  # if single
                        SINGLE.setdefault(hog_target[0],{})[descendant]=hog_descendant ## TODO ad check if already assign
                        ancestral_hogs_computed[descendant].add(hog_target[0])

        # For loss
        for descendant in self.descendants:
            lost = set(self.ancestor.genes) - ancestral_hogs_computed[descendant]
            for lost_hog in lost:
                LOSS.setdefault(lost_hog,[]).append(descendant)

        return LOSS, GAIN, SINGLE, DUPLICATE

    def computeUpMaps(self):
        upMaps = {}  # key|genome -> value|{key|HOGj -> value|HOGi or "None"}
        for descendant in self.descendants:
            upMaps[descendant] = self.build_UpMap(descendant)
        return upMaps

    def build_UpMap(self, descendant):
        upMap = {}  # value|{key|HOGj -> value|[HOGi or "None", paralog{True|False}]
        for hog_source in descendant.genes:
            if hog_source.parent is not None:  # avoid singletons TODO check for if taxon = gain...
                hog_target, paralog = self.search_ancestor_hog_in_ancestral_genome(hog_source, self.ancestor)
                upMap[hog_source] = [hog_target, paralog]
        return upMap

    def search_ancestor_hog_in_ancestral_genome(self, source_hog, ancestral_genome):
        found = None
        current_hog = source_hog
        paralog = current_hog.is_paralog
        while current_hog.parent is not None and found is None:
            current_hog = current_hog.parent
            if current_hog.genome == ancestral_genome:
                return current_hog, paralog
            if current_hog.is_paralog: # make sure been flagged as paralog is immuable
                paralog = True
        return found, paralog




