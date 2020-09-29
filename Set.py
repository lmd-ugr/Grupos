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

    
    def __sub__(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        
        #return Set(self.union(other))
        return Set(self | other)
    
    #A*B
    def __mul__(self, other):
        """Igual que cartesian"""
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        #return Set((a,b) for a in self for b in other)
        return Set(itertools.product(self,other, repeat=1))

    
    #A+B
    def __add__(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        
        #return Set(self.union(other))
        return Set(self | other)
    
    

    def pick(self):
        """Return an arbitrary element. (The finite Axiom of Choice is true!)"""
        
        if len(self) == 0:
            raise KeyError("This is an empty set")

        for item in self: break
        return item
    
    
    def cardinality(self):
        return len(self)
    

    def is_finite(self):
        return len(self)<math.inf
    
    
    def subsetsSize(self,n):
        if not isinstance(n, int) or n < 1:
            raise TypeError("Bad n, it must be greater than 1 ")
        return [Set(i) for i in itertools.combinations(self, n)]
    
    
    def subsets(self):
        #hasta len+1 porque el mismo conjunto es un subconjunto de él mismo
        return [Set(i) for p in range(1,len(self)+1) for i in itertools.combinations(self, p)]
        '''
        Z=[]
        for i in range(1, len(self)+1):
            for j in itertools.combinations(self, i):
                Z.append(Set(j))
        return Z
        '''
        
    def cartesian(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        #return Set((a,b) for a in self for b in other)
        return Set(itertools.product(self,other, repeat=1))
        
    #Funciones elementales y necesarias para devolver Set()
    def Union(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
        
        #return Set(self.union(other))
        return Set(self | other)
    
    
    def Intersection(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
            
        #return Set([x for x in self if x in other])
        return Set(self & other)
    
    
    
    def Difference(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
            
        return Set(self.difference(other))
                
    
    
    #Elementos de A y B que no están en ambos a la vez
    def SymetricDifference(self, other):
        if not isinstance(other, Set):
            raise TypeError("One of the objects is not a set")
            
        return Set(self ^ other)
        
    

if __name__ == '__main__':
    
    A=Set({1,2,3,4})
    B=Set({1,5})

    C=A*B
    print(1 in A)
    #print(C)
    print(" ")
    #print(C+B)
    
    '''
    for i in A.cartesian(B):
        for j in B.cartesian(A):
            if i==j:
                print(i) #Bien, solamente el (1,1) es igual, (a,b) != (b,a) forall a,b in A,B
    '''   
     
    
    
    #print(A.is_finite())
    #print(A.cardinality())
    #print(C.cardinality())
    
    #print(A.subsets())
    
    
    #print(A.subsets())
    #print(A.Union(B))
    #print(A.SymetricDifference(B))
    #print(A.symmetric_difference(B))
    
    
    
    #P = A.subsets()
    #print(P)
    #print(type(P[3]))
    #print(A.subsets())
    #print(A.cartesian(B))
    