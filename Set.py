"""
Implementation of Set
"""
import itertools
import math

class Set(frozenset):
    """
    Definition of a Set

    It's important that Set be a subclass of frozenset, (not set), because:
    1) it makes Set immutable
    2) it allows Set to contains Sets
    """
    def __repr__(self):
        return str(set(self))

    def __str__(self):
        return str(set(self))

    def __add__(self, other):
        """ 
        Returns the union of self and other.
        
        Args:
            other : instance of Set.
        """
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")

        #return Set(self.union(other))
        return Set(self | other)

    union = __add__


    def __sub__(self, other):
        """ 
        Returns the difference of self and other.
        
        Args:
            other : instance of Set.
        """
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")

        #return Set(self.union(other))
        return Set(self.difference(other))

    Difference = __sub__

    def __mul__(self, other):
        """
        Cartesian product of self x other.
        
        Args:
            other : instance of Set.
        """
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        return Set(itertools.product(self,other, repeat=1))
        #return Set((a,b) for a in self for b in other)

    cartesian = __mul__

    def pick(self):
        """
        Return an arbitrary element. (The finite Axiom of Choice is true!)
        """

        if len(self) == 0:
            raise KeyError("This is an empty set")

        for item in self: break
        return item



    def cardinality(self):
        """ 
        Returns the cardinality of self.
        """
        return len(self)


    def is_finite(self):
        """ 
        Checks if self is finite.
        """
        return len(self)<math.inf


    def subsets(self, n=None):
        """ 
        Returns the subsets of self. If n is specified then it returns the
        subgroups of order n.
        
        Args:
            n : integer
        """
        if n != None:
            if not isinstance(n, int) or n < 1:
                raise TypeError("Bad n, it must be greater than 1 ")
            return [Set(i) for i in itertools.combinations(self, n)]
            
        #hasta len+1 porque el mismo conjunto es un subconjunto de él mismo
        else:
            return [Set(i) for p in range(1,len(self)+1) for i in itertools.combinations(self, p)]
 


    def Intersection(self, other):
        """ 
        Returns the intersección of self and other
        
        Args:
            other : instance of Set.
        """
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")

        #return Set([x for x in self if x in other])
        return Set(self & other)




    #Elementos de A y B que no están en ambos a la vez
    def SymmetricDifference(self, other):
        """ 
        Returns the symmetric difference of self and other.
        
        Args:
            other : instance of Set.
        """
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")

        return Set(self ^ other)
