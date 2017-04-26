from abc import ABCMeta, abstractmethod
from ham.genome import Genome
import logging

logger = logging.getLogger(__name__)


class HOGsMap(object):
    """
    Class to compare HOG evolutionary relations (duplication, loss, gain and identical) across 2 genomes that are on
    the same lineage.
    
    Let's consider Go the oldest genome and Gy the youngest genomes with their respected HOGs Ho and Hy.
    
    Those relations are classified into the following categories:
        - IDENTICAL: if Ho will correspond to a single Hy.
        - DUPLICATE: If a duplication occurred in Ho in between Go ang Gy.
        - LOSS: If Ho have no representative Hy in Gy.
        - GAIN: if Hy have no ancestor Ho in Go.
        
    First, a "UpMap" is build to find for each Hy its potential ancestor Ho. If during the traversal from Hy to Ho a
    duplication occurred, a flag "paralog" is added in the UpMap[Hy]. Hy with no corresponding Ho are considered as 
    gained meanwhile Ho without corresponding Hy are considered as lost. Ho with only one corresponding Hy (and no
    duplication tag) are clustered as identical otherwhise as duplicated.


    Attributes:
        HAM (:obj:`HAM`): HAM object.
        ancestor (:obj:`Genome`): :obj:`Genome` of Go.
        descendant (:obj:`Genome`): :obj:`Genome` of Gy.
        upMap (:obj:`dict`):  a dictionary that map each Ho with its related Hy (or None if no HOG founded) associated
        a boolean if a duplication occurs or not in between them.
        IDENTICAL (:obj:`dict`): Dictionary that map a Ho with its descendant Hy.
        DUPLICATE (:obj:`dict`): Dictionary that map a Ho with its list of descendants Hy.
        LOSS (:obj:`list`): a list of Ho with no matching in Gy.
        GAIN: (:obj:`list`): a list of Hy with no ancestor in Go.
    
    
    """

    def __init__(self, ham_object, genome1, genome2):
        """
        Args:
            genome1 (:obj:`Genome`): First :obj:`Genome` to compare.
            genome2 (:obj:`Genome`): Second :obj:`Genome` to compare.
            
        """

        if not isinstance(genome1, Genome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Genome.__name__,
                                    type(genome1).__name__))

        if not isinstance(genome2, Genome):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(Genome.__name__,
                                    type(genome2).__name__))

        self.HAM = ham_object
        self.ancestor, self.descendant = self.HAM._get_oldest_from_genome_pair(genome1, genome2)
        self.upMap = self._build_UpMap()
        self.LOSS, self.GAIN, self.IDENTICAL, self.DUPLICATE = self.build_event_clusters()

    def build_event_clusters(self):
        """  
        This method builds all the event clusters based on the UpMap.

        Returns:
            list of loss, list of gain, dict of identical and dict of duplicated.

        """
        GAIN = []  # [Hy]
        IDENTICAL = {}  # {Ho -> Hy}
        DUPLICATE = {}  # {Ho -> [Hy]}
        ancestral_hogs_computed = set()

        for Hy, Ho in self.upMap.items():
            if Ho[0] is None:  # GAIN
                GAIN.append(Hy)
            else:
                if Ho[1]:  # DUPLICATE
                    DUPLICATE.setdefault(Ho[0], []).append(Hy)
                else:  # SINGLE
                    IDENTICAL[Ho[0]] = Hy
                ancestral_hogs_computed.add(Ho[0])

        # LOSS
        LOSS = set(self.ancestor.genes) - ancestral_hogs_computed

        return LOSS, GAIN, IDENTICAL, DUPLICATE

    def _build_UpMap(self):
        """  
        This method construct the UpMap between Ho and Hy.

        Returns:
            A dictionary that map each Ho with its related Hy (or None if no HOG founded) associated a boolean if a
            duplication occurs or not in between them.
        """

        upMap = {}  # {Hy -> [Ho|None, True|False]}

        for Hy in self.descendant.genes:
            if Hy.is_singleton() is False:
                Ho, paralog = Hy.search_ancestor_hog_in_ancestral_genome(self.ancestor)
                upMap[Hy] = [Ho, paralog]

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
        return self.map.IDENTICAL

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
        self.IDENTICAL = None # {Hi -> {Gj -> Hj}}
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
        if self.IDENTICAL is None:
            single = {}
            for genome, hmap in self.maps.items():
                for hi,hj in hmap.IDENTICAL.items():
                    single.setdefault(hi,{})[genome] = hj
            self.IDENTICAL = single
        return self.IDENTICAL

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
