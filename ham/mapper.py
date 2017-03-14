from abc import ABCMeta, abstractmethod
import copy
import logging

logger = logging.getLogger(__name__)

class HOGsMap(object):
    '''
    Given a set of 2 genomes where Gi denote the ancestor of Gj, the HOGsMap object identify the relations that exist between all Hi
    and Hj respectively the HOGs from Gi and Gj.

    Those relations are classified into the following categories:
        - SINGLE: if Hi will correspond to a single Hj
        - DUPLICATE: If a duplication occurred in the path between Hi and Hj
        - LOSS: If Hi is not bound to any Hj
        //- DUPLICATED_BUT_LOST: --> TODO NOt ure about this
        - GAIN: if Hj have no corresponding Hi

    First, a "UpMap" is build to find the ancestor Hi of each Hj.
    If during the traversal from Hj to Hi a duplication occurred, a flag "paralog" is added in the UpMap[Hi|Hj].
    All Hj in Gj that have no corresponding ancestor Hi in Gi are considered as GAIN.

    Then we built the Event Clusters:
        - .LOSS: [Hi]
        - .GAIN: [Hj]
        - .SINGLE: {Hi -> Hj}
        - .DUPLICATE: {Hi -> [Hj]}

    '''

    def __init__(self, ham_object, genome_set):
        self.HAM = ham_object
        self.verify_genome_set(copy.copy(genome_set))
        self.ancestor, self.descendant = self.set_ancestor_and_descendant(genome_set)
        self.upMap = self.build_UpMap()
        self.LOSS, self.GAIN, self.SINGLE, self.DUPLICATE = self.buildEventClusters()

    def verify_genome_set(self, genome_set):
        if len(genome_set) != 2:
            raise TypeError("Building the HOGMaps required two and only two genomes, invalid genome set: {}".format(genome_set))
        anc, desc = self.HAM._get_ancestor_and_descendant(copy.copy(genome_set))
        x = self.HAM._get_oldest_from_genome_pair(list(genome_set)[0].taxon, list(genome_set)[1].taxon)
        if x is None:
           raise TypeError("The genomes are not in the same lineage: {}".format(genome_set))

    def set_ancestor_and_descendant(self, genome_set):
        ancestor = self.HAM._get_mrca_ancestral_genome_from_genome_set(genome_set)
        genome_set.discard(ancestor)
        descendant = genome_set.pop()
        return ancestor, descendant

    def buildEventClusters(self):
        LOSS = []  # [Hi]
        GAIN = []  # [Hj]
        SINGLE = {}  # {Hi -> Hj}
        DUPLICATE = {}  # {Hi -> [Hj]}

        ancestral_hogs_computed = set()

        # feed the MAPS
        for hog_descendant, hog_target in self.upMap.items():
            if hog_target[0] is None: # If gained
                GAIN.append(hog_descendant)
            else:  # If not gained
                if hog_target[1] is True:  # is duplicated
                    DUPLICATE.setdefault(hog_target[0],[]).append(hog_descendant)
                    ancestral_hogs_computed.add(hog_target[0])
                else:  # if single
                    SINGLE[hog_target[0]] = hog_descendant
                    ancestral_hogs_computed.add(hog_target[0])

        # For loss
        lost = set(self.ancestor.genes) - ancestral_hogs_computed

        for lost_hog in lost:
            LOSS.append(lost_hog)

        return LOSS, GAIN, SINGLE, DUPLICATE

    def build_UpMap(self):
        upMap = {}  # value|{key|HOGj -> value|[HOGi or "None", paralog{True|False}]
        for hog_source in self.descendant.genes:
            if hog_source.parent is not None:  # avoid singletons TODO check for if taxon = gain...
                hog_target, paralog = hog_source.search_ancestor_hog_in_ancestral_genome(self.ancestor)
                upMap[hog_source] = [hog_target, paralog]
        return upMap


class MapResults(metaclass=ABCMeta):
    '''
    The MapResult class manages the HOG mapping between several {extant|ancestral} genomes. There is two type of
    map result:
        - Lineage: for 2 genomes on a same lineages
        - Lateral: for genomes not on the same lineage that are all descendant from the same ancestral genomes.

    The MapResult provides the following methods to compare HOG across genomes:
        - get_loss()
        - get_gain()
        - get_duplicate()
        - get_single()
    '''

    def __init__(self, HAM):
        self.HAM = HAM
        self.ancestor = None

    @abstractmethod
    def add_map(self, map):
        """Add a map to the MapResults object"""
        pass

    @abstractmethod
    def get_lost(self):
        """Return the lost genes"""
        pass

    @abstractmethod
    def get_gained(self):
        """Return the gained genes"""
        pass

    @abstractmethod
    def get_single(self):
        """Return genes that stay single copy"""
        pass

    @abstractmethod
    def get_duplicated(self):
        """Return genes that duplicate"""
        pass

class MapVertical(MapResults):
    '''
    The MapVertical class manages the HOG mapping between two genomes on the same lineages (an ancestral genome Gi
    (with HOG Hi) and its descendant Gj (with HOG Hj)).

    The MapVertical provides the following methods to compare HOG across genomes:
        - get_loss(): return a list of all Hi that are lost in Gj.
        - get_gain(): return a list of all Hj that are novels in Gj and not found in Gi.
        - get_duplicate(): return a dictionary with all Hi map to their descendant Hj.
        - get_single():  return a dictionary with all Hi map to its descendant Hj.
    '''

    def __init__(self, HAM):
        super(MapVertical, self).__init__(HAM)
        self.descendant = None
        self.map = None

    def add_map(self, HogMap):
        """Add a map to the MapResults object"""
        if not isinstance(HogMap, HOGsMap):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(HOGsMap.__name__,
                                    type(HogMap).__name__))
        if self.map != None:
            raise TypeError("MapVertical can only contains one HOGMap object.")
        self.ancestor = HogMap.ancestor
        self.descendant = HogMap.descendant
        self.map = HogMap

    def get_lost(self):
        return self.map.LOSS

    def get_gained(self):
        return self.map.GAIN

    def get_single(self): # TODO return_list_only=False to return list of gene directly ?
        return self.map.SINGLE

    def get_duplicated(self):
        return self.map.DUPLICATE

class MapLateral(MapResults):

    '''
    The MapLateral class manages the HOG mapping between genomes Gx (with HOG Hx) that are not on the same lineage but share an ancestral genome Gi
    (with HOG Hi).

    The MapLateral provides the following methods to compare HOG across genomes:
        - get_loss(): return a dictionary of all Hi that are lost in Gj. {Hi -> [Gx]}
        - get_gain(): return a dictionary of all Gx and their novels HOG Hx with no correspondence in Gi. {Gx -> [Hx]}
        - get_duplicate(): return a dictionary with all Hi map to their descendant Hj clustered by Gj. {Hi -> {Gj -> [Hj]}}
        - get_single():  return a dictionary with all Hi map to its descendant Hx clustered by Gx. {Hx -> {Gx -> Hx}}
    '''


    def __init__(self, HAM):
        super(MapLateral, self).__init__(HAM)
        self.descendants = [] # [Gx]
        self.maps = {} # {Gx -> HOGMap(Gi,Gx)}
        self.LOSS = None  # {Hi -> [Gx]}
        self.GAIN = None # {Gj -> [Hj]}
        self.SINGLE = None # {Hi -> {Gj -> Hj}}
        self.DUPLICATE = None # {Hi -> {Gj -> [Hj]}}

    def add_map(self, HogMap):
        """Add a map to the MapResults object"""
        if not isinstance(HogMap, HOGsMap):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(HOGsMap.__name__,
                                    type(HogMap).__name__))
        if self.maps == {}:
            self.ancestor = HogMap.ancestor
        else:
            if self.ancestor != HogMap.ancestor:
                raise TypeError("You can only add map to MapLateral object that have the same ancestor. The current ancestor is {} and the map ancestor you try to add is {}.".format(self.ancestor, HogMap.ancestor ))
        self.maps[HogMap.descendant] = HogMap
        self.descendants.append(HogMap.descendant)


    def get_lost(self): # Todo only_in=x where x can be descendant genome or "all" ?
        """
        :return:  a dictionary all HOGs Hi and their list of Genome Gj where they are lost.
        """
        if self.LOSS is None:
            loss = {}
            for genome, hmap in self.maps.items():
                for lost_gene in hmap.LOSS:
                    loss.setdefault(lost_gene,[]).append(genome)
            self.LOSS = loss
        return self.LOSS

    def get_gained(self):
        """
        :return: a dictionary of Genome Gj wither their novels genes Hj.
        """
        if self.GAIN is None:
            gain = {}
            for genome, hmap in self.maps.items():
                gain[genome] = hmap.GAIN
            self.GAIN = gain
        return self.GAIN

    def get_single(self): # TODO return_list_only=False to return list of gene directly ?
        """
        :return: a dictionary of HOG Hi with their single Hj clustered by Genome Gj.
        """
        if self.SINGLE is None:
            single = {}
            for genome, hmap in self.maps.items():
                for hi,hj in hmap.SINGLE.items():
                    single.setdefault(hi,{})[genome] = hj
            self.SINGLE = single
        return self.SINGLE

    def get_duplicated(self):
        """
        :return: a dictionary of HOG Hi with their list of single Hj clustered by Genome Gj.
        """
        if self.DUPLICATE is None:
            duplicate = {}
            for genome, hmap in self.maps.items():
                for hi,list_hj in hmap.DUPLICATE.items():
                    duplicate.setdefault(hi,{})[genome] = list_hj
            self.DUPLICATE = duplicate
        return self.DUPLICATE
