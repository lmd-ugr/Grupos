# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 18:24:47 2020

@author: Alberto
"""

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

    #def __mul__(self, other):
    def cartesian(self, other):
        """Cartesian product"""
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        return Set(itertools.product(self,other, repeat=1))
        #return Set((a,b) for a in self for b in other)

    #Siempre devuelve 1 ¿?
    def pick(self):
        """Return an arbitrary element. (The finite Axiom of Choice is true!)"""
        
        if len(self) == 0:
            raise KeyError("This is an empty set")

        for item in self: break
        return item
    
    
    '''BOOM##################################
    
    ##########################
    Añado funciones, practising
    ##########################
    
    '''


    def cardinality(self):
        return len(self)
    

    def is_finite(self):
        return len(self)<math.inf
    
    
    #Subconjuntos de tamaño n
    def subsetsSize(self,n):
        if not isinstance(n, int) or n < 1:
            raise TypeError("Bad n, it must be greater than 1 ")
        return [Set(i) for i in itertools.combinations(self, n)]
    
    
    #Subconjuntos
    def subsets(self):
        return [Set(i) for p in range(1,len(self)) for i in itertools.combinations(self, p)]
        '''
        Z=[]
        for i in range(2, len(self)):
            for j in itertools.combinations(self, i):
                Z.append(Set(j))
        return Z
        '''
        
    #Funciones elementales y necesarias para hacer el cast Set()
    def Union(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        
        #return Set(self.union(other)) #frozenset, no set
        return Set(self | other)
    
    
    def Intersection(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
            
        #return Set([x for x in self if x in other])
        return Set(self & other)
    
    
    
    def Difference(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
            
        return Set(self - other)
    
    #Elementos de A y B que no están en ambos a la vez
    def SymetricDifference(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
            
        return Set(self ^ other)
        
    


#A=Set({1,2,3})
#print(A.subsets())
#B=Set({3,4})

#print(A.Union(B))
#print(A.SymetricDifference(B))
#print(A.symmetric_difference(B))

#print(C.subsets())

#print(A.__str__)  -------------?
#print(A.__repr__) ---------------?
#Z=[Set({1,2,3,4})]
#Z.append(A)
#print(A.is_finite())
#print(A.subsetsSize(2))
#print(Z)

#print(A.subsets())

#print(A.subsets())

#print(A.cartesian(B))
