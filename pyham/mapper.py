from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import super
from builtins import object
from future import standard_library
standard_library.install_aliases()
from .genome import Genome
import logging
import abc, six

logger = logging.getLogger(__name__)


class HOGsMap(object):
    """
    Class to compare HOG evolutionary relations (duplication, loss, gain and identical) between an ancestor genome and
    its descendant genome.
    
    Let's consider Go the oldest genome and Gy the youngest genome with their respected HOGs Ho and Hy.
    
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
        | Ham (:obj:`pyham.pyham.Ham`): Ham object.
        | ancestor (:obj:`pyham.genome.Genome`): :obj:`pyham.genome.Genome` of Go.
        | descendant (:obj:`pyham.genome.Genome`): :obj:`pyham.genome.Genome` of Gy.
        | upMap (:obj:`dict`):  a dictionary that map each Ho with its related Hy (or None if no HOG founded) associated a boolean if a duplication occurs or not in between them.
        | IDENTICAL (:obj:`dict`): Dictionary that map a Ho with its descendant Hy.
        | DUPLICATE (:obj:`dict`): Dictionary that map a Ho with its list of descendants Hy.
        | LOSS (:obj:`list`): a list of Ho with no matching in Gy.
        | GAIN: (:obj:`list`): a list of Hy with no ancestor in Go.
    
    
    """

    def __init__(self, ham_object, genome1, genome2):
        """
        Args:
            | genome1 (:obj:`pyham.genome.Genome`): First :obj:`pyham.genome.Genome` to compare.
            | genome2 (:obj:`pyham.genome.Genome`): Second :obj:`pyham.genome.Genome` to compare.
            
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
        self.LOSS, self.GAIN, self.IDENTICAL, self.DUPLICATE = self._build_event_clusters()

    def _build_event_clusters(self):
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
            Ho, paralog = Hy.search_ancestor_hog_in_ancestral_genome(self.ancestor)
            upMap[Hy] = [Ho, paralog]

        return upMap


@six.add_metaclass(abc.ABCMeta)
class MapResults(object):
    """
    Abstract class to map HOGs across multiple genomes through their most recent common ancestral genome. The HOGs are all
    clustered based on their relation with the mrca genome HOGs (duplicated, lost, gained and identical).

    Let's consider Go the mrca genome and Gn the youngest genomes with their respected HOGs Ho and Hn.

    The HOG relations between Ho and Hn are classified into the following categories:
        - IDENTICAL: if Ho will correspond to a single Hn in Gn.
        - DUPLICATE: If a duplication occurred in Ho in between Go ang Gn.
        - LOSS: If Ho have no representative Hn in Gn.
        - GAIN: if Hn have no ancestor Ho in Go.

    There is two type of possible mapping:
        - vertical: compare a genome with its ancestor (that act as "mrca").
        - lateral: compare genomes with their mrca.
    
    The MapResults provides methods to compare HOGs across genomes:
        - get_loss()
        - get_gain()
        - get_duplicate()
        - get_identical()

    Attributes:
        | Ham (:obj:`pyham.pyham.Ham`): Ham object.
        | ancestor (:obj:`pyham.genome.Genome`): :obj:`pyham.genome.Genome` of Go.
    """

    def __init__(self, HAM):
        self.HAM = HAM
        self.ancestor = None

    @abc.abstractmethod
    def add_map(self, map):
        """Add a map to the MapResults object"""
        pass

    @abc.abstractmethod
    def get_lost(self):
        """Return the lost genes"""
        pass

    @abc.abstractmethod
    def get_gained(self):
        """Return the gained genes"""
        pass

    @abc.abstractmethod
    def get_identical(self):
        """Return genes that stay single copy"""
        pass

    @abc.abstractmethod
    def get_duplicated(self):
        """Return genes that duplicate"""
        pass


class MapVertical(MapResults):
    """
    Class to map HOGs between a genome and its ancestor. 

    Let's consider Go the oldest genome and Gn the youngest genome with their respected HOGs Ho and Hn.
    
    Attributes:
        | descendant (:obj:`pyham.genome.Genome`): :obj:`pyham.genome.Genome` of Gn.
        | map (:obj:`pyham.mapper.HOGsMap`): :obj:`pyham.genome.Genome` of Go.
    """

    def __init__(self, ham):
        super(MapVertical, self).__init__(ham)
        self.descendant = None
        self.map = None

    def add_map(self, HogMap):
        """  
        Method to set the HOGsMap.

            Args:
                HogMap (:obj:`pyham.mapper.HOGsMap`): HOGsMap to add.

            Raises:
                TypeError: if HogMap is not :obj:`pyham.mapper.HOGsMap` or if a :obj:`pyham.mapper.HOGsMap` is already set.

        """
        if not isinstance(HogMap, HOGsMap):
            raise TypeError("expect subclass obj of '{}', got {}"
                            .format(HOGsMap.__name__,
                                    type(HogMap).__name__))
        if self.map is not None:
            raise TypeError("MapVertical can only contains one HOGMap object.")
        self.ancestor = HogMap.ancestor
        self.descendant = HogMap.descendant
        self.map = HogMap

    def get_lost(self):
        """
        Getter for lost genes.
        
        Returns:
            List of Ho.
        """

        return self.map.LOSS

    def get_gained(self):
        """
        Getter for gained genes.

        Returns:
            List of Hn.

        """
        return self.map.GAIN

    def get_identical(self):
        """
         Getter for identical genes.

         Returns:
             Dictionary of Ho mapped to their descendant Hy.
         """

        return self.map.IDENTICAL

    def get_duplicated(self):
        """
        Getter for duplicated genes.

        Returns:
            Dictionary of Ho mapped to their list of descendants Hy.
        """

        return self.map.DUPLICATE


class MapLateral(MapResults):
    """
    
    Class to map HOGs between genomes through their ancestor. 

    Let's consider Go the oldest genome and Gn the youngest genomes with their respected HOGs Ho and Hn.

    Attributes:
        | descendant (:obj:`list`): list of Gn.
        | map (:obj:`dict`): dictionary of Gn mapped to a HOGMap(Go,Gn).
        | LOSS (:obj:`dict`): dictionary of Ho mapped to the list of Gn where this HOG is lost.
        | GAIN (:obj:`dict`): dictionary of Gn mapped to their list of gained Hn.
        | IDENTICAL (:obj:`dict`): dictionary of Ho mapped to a dict of Gn mapped to a Hn.
        | DUPLICATE (:obj:`dict`): dictionary of Ho mapped to a dict of Gn mapped to a list of Hn.
    
    """

    def __init__(self, HAM):
        super(MapLateral, self).__init__(HAM)
        self.descendants = []
        self.maps = {}
        self.LOSS = None
        self.GAIN = None
        self.IDENTICAL = None
        self.DUPLICATE = None

    def add_map(self, HogMap):
        """  
        Method to add the HOGsMap.

            Args:
                HogMap (:obj:`pyham.mapper.HOGsMap`): HOGsMap to add.

            Raises:
                TypeError: if HogMap is not :obj:`pyham.mapper.HOGsMap` or if HOGMap ancestor doesn't match the MapLateral one.
        """

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

    def get_lost(self):
        """
        Lazy getter for lost genes.
        
        Returns:
            dictionary of Ho -> [Gn].
        """

        if self.LOSS is None:
            loss = {}
            for genome, hmap in self.maps.items():
                for lost_gene in hmap.LOSS:
                    loss.setdefault(lost_gene, []).append(genome)
            self.LOSS = loss

        return self.LOSS

    def get_gained(self):
        """
        Lazy getter for gained genes.

        Returns:
            dictionary of Gn -> [Hn].
        """
        if self.GAIN is None:
            gain = {}
            for genome, hmap in self.maps.items():
                gain[genome] = hmap.GAIN
            self.GAIN = gain
        return self.GAIN

    def get_identical(self):
        """
        Lazy getter for identical genes.

        Returns:
            dictionary of Hi -> dict of Gn -> Hn.
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
        Lazy getter for duplicated genes.

        Returns:
            dictionary of Hi -> dict of Gn -> list of Hn.
        """
        if self.DUPLICATE is None:
            duplicate = {}
            for genome, hmap in self.maps.items():
                for hi,list_hj in hmap.DUPLICATE.items():
                    duplicate.setdefault(hi,{})[genome] = list_hj
            self.DUPLICATE = duplicate

        return self.DUPLICATE
