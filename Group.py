# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:18:18 2020
"""

import itertools
import functools
import operator
import math
from Set import Set
from Function import Function
from fractions import gcd
from copy import deepcopy
from sympy.ntheory import factorint,totient
from sympy.utilities.iterables import flatten

class GroupElem:
    """
    Group element definition
    This is mainly syntactic sugar, so you can write stuff like g * h
    instead of group.bin_op(g, h), or group(g, h).
    """

    def __init__(self, elem, group):
        if not isinstance(group, Group):
            raise TypeError(str(group) + " is not a Group")
        if not elem in group.Set:
            raise TypeError( str(elem) + " is not an element of group")
        self.elem = elem
        self.group = group

    def __str__(self):
        return str(self.elem)

    #Se define así para que al llamar al objeto se tome únicamente el objeto
    def __repr__(self):
        return repr(self.elem)

    def __eq__(self, other):
        """
        Two GroupElems are equal if they represent the same element in the same group
        """

        if not isinstance(other, GroupElem):
            return False #raise TypeError("other is not a GroupElem")
        return (self.elem == other.elem) and (self.group.parent==other.group.parent)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.elem)

    def __mul__(self, other):
        """
        If other is a group element, returns self * other.
        If other = n is an int, and self is in an abelian group, returns self**n
        """

        if isinstance(other,Group):
            return Set(self.group.bin_op((self.elem, h)) for h in other.Set)
        #if isinstance(other,set):
        #    return set([self*g for g in other])
        if self.group.is_abelian() and isinstance(other, (int)):
            return self ** other

        if not isinstance(other, GroupElem):
            raise TypeError("other must be a GroupElem, or an int " \
                            "(if self's group is abelian)")
        if not(self.group.parent==other.group.parent):
            raise TypeError("both elements must be in the same group")
        #try:
        return GroupElem(self.group.parent.bin_op((self.elem, other.elem)), \
                             self.group.parent)
        # This can return a TypeError in Funcion.__call__ if self and other
        # belong to different Groups. So we see if we can make sense of this
        # operation the other way around.
        #except TypeError:
        #    return other.__rmul__(self)

    def __rmul__(self, other):
        """
        If other is a group element, returns other * self.
        If other = n is an int, and self is in an abelian group, returns self**n
        """
        if self.group.is_abelian() and isinstance(other, (int)):
            return self ** other

        if not isinstance(other, GroupElem):
            raise TypeError("other must be a GroupElem, or an int " \
                            "(if self's group is abelian)")

        return GroupElem(self.group.bin_op((other.elem, self.elem)), self.group)

    def __add__(self, other):
        """Returns self + other for Abelian groups"""
        if self.group.is_abelian():
            return self * other
        raise TypeError("not an element of an abelian group")

    def __pow__(self, n, modulo=None):
        """
        Returns self**n
        modulo is included as an argument to comply with the API, and ignored
        """
        if not (isinstance(n, int)):
            raise TypeError("n must be an int or a long")

        if n == 0:
            return self.group.e
        elif n < 0:
            return self.group.inverse(self) ** -n
        elif n % 2 == 1:
            return self * (self ** (n - 1))
        else:
            return (self * self) ** (n // 2)
    
    __xor__=__pow__
    
    def __neg__(self):
        """Returns self ** -1 if self is in an abelian group"""
        if not self.group.is_abelian():
            raise TypeError("self must be in an abelian group")
        return self ** (-1)

    def __sub__(self, other):
        """Returns self * (other ** -1) if self is in an abelian group"""
        if not self.group.is_abelian():
            raise TypeError("self must be in an abelian group")
        return self * (other ** -1)


    def conjugate(self,g):
        return g*self*g**-1

    def conjugacy_class(self):
        """
        Returns the conjugacy class of self in self.group
        """
        return Set([g*self*g**-1 for g in self.group])
    
    
    def centralizer(self):
        """
        Returns the centralizer of self, that is, the set of elements that commute with self
        """
        G=self.group
        def prop(g):
            return g*self  ==  self*g
        return G.subgroup_search(prop)

    def order(self):
        if not self in self.group:
            raise ValueError
    
        if self == self.group.identity():
            ind = 1
        else:
            i=1
            x=self
            
            while not x == self.group.identity():
                #aux=x
                x = x*self
                #print(x, "=", aux,"*s",self)
                i = i + 1
            ind = i
        return ind
    
    
    def inverse(self):
        """Returns the inverse of elem"""
        
        if not self in self.group:
            raise TypeError("Element isn't a GroupElem in the Group")
        
        #El elemento a debe ser un entero para poder aplicarle
        #la operación binaria. Lo detecta como groupElem
        for a in self.group.elements():
            if self.group.bin_op((self.elem, a.elem)) == self.group.e.elem:
                return a
        raise RuntimeError("Didn't find an inverse for g")

    

    
'''
def is_associative(Set, bin_op):
    if all(bin_op((a, bin_op((b, c)))) == \
          bin_op((bin_op((a, b)), c)) \
              for a in Set for b in Set for c in Set):
        return True
    else:
        return False
        
def has_inverses(Set, bin_op, e):
        
    # Test for inverses
    for a in G:
        if not any(bin_op((a,  b)) == e for b in Set):
            return False
    return True
    


def has_identity(Set, bin_op, e):
    found_id = False
    for e in Set:
        if all(bin_op((e, a)) == a for a in Set):
            return True
    return False
    
'''

    
class Group:
    

    def __init__(self, G, bin_op, identity=None, parent=None, group_order=None, group_degree=None):
        """Create a group, checking group axioms"""
        
        if not isinstance(G, Set): 
            raise TypeError("G must be a set")
        if not isinstance(bin_op, Function):
            raise TypeError("bin_op must be a function")
        #if bin_op.codomain != G:
         #   raise TypeError("binary operation must have codomain equal to G")
        #if bin_op.domain != G.cartesian(G):
         #   raise TypeError("binary operation must have domain equal to G x G")
        
        # Find the identity
        
        if identity in G:
            e=identity
            #found_id=True
        else:
            found_id = False
            for e in G:
                if all(bin_op((e, a)) == a for a in G):
                    found_id = True
                    break
                
            if not found_id:
                raise ValueError("G doesn't have an identity")
        '''
        found_id = False
        for e in G:
            if all(bin_op((e, a)) == a for a in G):
                found_id = True
                break
        if not found_id:
            print("Not identity")
            raise ValueError("G doesn't have an identity")
        '''
       
        # Test for inverses
        for a in G:
            if not any(bin_op((a,  b)) == e for b in G):
                #print("Not inverses")
                raise ValueError("G doesn't have inverses")
               
        # Test associativity
        if not all(bin_op((a, bin_op((b, c)))) == \
                   bin_op((bin_op((a, b)), c)) \
                   for a in G for b in G for c in G):
            print("Not associative")
            raise ValueError("binary operation is not associative")

        # At this point, we've verified that we have a Group.
        # Now we determine if the Group is abelian:
        '''
        if not(isinstance(abelian,bool)):
            self.abelian = all(bin_op((a, b)) == bin_op((b, a)) \
                               for a in G for b in G)
        else:
            self.abelian=abelian
        '''
            
        self.Set = G
        self.group_elems = Set(GroupElem(g, self) for g in G)
        
        self.e = GroupElem(e, self)
        self.bin_op = bin_op
        if parent==None:
            self.parent=self
        else:
            self.parent=parent
        self.group_gens=list(self.group_elems)
        self.group_order=group_order
        self.group_degree=group_degree



    def __iter__(self):
        """Iterate over the GroupElems in G, returning the identity first"""
        yield self.e
        for g in self.group_elems:
            if g != self.e: yield g

    def __contains__(self, item):
        return item in self.group_elems

    def __hash__(self):
        return hash(self.Set) ^ hash(self.bin_op)

    def __eq__(self, other):
        if not isinstance(other, Group):
            return False

        return id(self) == id(other) or \
               (self.Set == other.Set and self.bin_op == other.bin_op)

    def __call__(self,el):
        return GroupElem(el,self)

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return self.Set.cardinality()

    def __str__(self):
        return "Group with "+str(len(self))+" elements: " + str(self.Set)

    def __repr__(self):
        if self.group_gens!=None:
            gs = "Group( "+str(self.group_gens)+" )"
            if len(gs)>100:
                gs = "Group with "+str(len(self))+" elements"
        else:
            gs = "Group with "+str(len(self))+" elements"
        return gs
    
    def elements(self):
        return self.group_elems
    
    def identity(self):
        return self.e
    
    def inverse(self, g):
        """Returns the inverse of elem"""
        if not g in self.group_elems:
            raise TypeError("g isn't a GroupElem in the Group")
        for a in self:
            if g * a == self.e:
                return a
        raise RuntimeError("Didn't find an inverse for g")

    
    
    def is_abelian(self):
        return all(self.bin_op((a, b)) == self.bin_op((b, a)) \
                               for a in self.Set for b in self.Set)
    
    
    #abcdefghijclme
    #pip install beautifultable
    def Cayley_table(self):
        from beautifultable import BeautifulTable
        head=list(self.Set)
        
        table = BeautifulTable()
        #table.colums_headers = head
        
        table.rows.insert(0, head, "*")
        
        c=1
        for a in self.Set:
           l1 = []
           for b in self.Set:
               l1.append(self.bin_op((a,b)))
           
           
           table.rows.insert(c+1, l1, str(head[c-1]))
           c=c+1
           #print(a, "*", b,"=",self.bin_op((a,b)))
        table.set_style(BeautifulTable.STYLE_BOX)
        return table
        
        
    def cardinality(self):
        return self.Set.cardinality()
    
    def order(self):
        """Return the order of the group.
        """
        if self.group_order != None:
            return self.group_order
        self.group_order = len(self)
        return self.group_order
    
    
    def is_subgroup(self, other):
        #Return True if all elements of self belong to other.
        
        if not isinstance(other, Group):
            #print("No es un grupo")
            return False
        if other.order() % self.order() != 0:
            #print("El subgrupo no tiene la misma cardinalidad")
            return False
        
        if self == other or self.parent==other:
            return True
        
      
           
        if ( self.Set in other.Set.subsets() and \
            all(other.bin_op((a.elem,b.elem)) in self.Set
                for a in self for b in self) and \
            all(other.bin_op((a.elem,b.inverse().elem )) in self.Set
                for a in self for b in self) ):
            return True
        else:
            return False
        
        '''
        
        #¿Incompleto? está comprobando solo que a*b = a*b ???
        #print(self.Set in other.Set.subsets())
        
        return self.Set in other.Set.subsets() and \
            all(self.bin_op((a.elem, b.elem)) == other.bin_op((a.elem, b.elem)) \
                for a in self.group_gens for b in self.group_gens)
        '''
  
    def index(self, other):
        if not isinstance(other,Group):
            raise TypeError(other, " is not a Group")
        return self.order()/other.index()
    
    
    def is_normalSubgroup(self, other):
        if not isinstance(other, Group):
            return False
        
        if self.is_subgroup(other):
            for g in other:
                for h in self:
                    #print(g, "*", h, "*", g.inverse(), " in ", Set(self.group_elems))
                    p = self.bin_op((self.bin_op((g.elem,h.elem)), g.inverse().elem))
                    #p= g.elem*h.elem*(g**-1).elem
                    if not p in self.Set:
                        return False
            return True
            
        
        if self.is_abelian() or other.index(self)==2:
            return True 
        
    '''
    def all_normalSubgroups(self):
        l = self.all_subgroups()
        l2=[]
        for a in l:
            if a.is_normalSubgroup(self):
                l2.append(a)
        return l2
    '''
    
    
    def all_normalSubgroups(self, order=None):
        if order==None:
            return [gr for gr in self.all_subgroups() if gr.is_normalSubgroup(self)]
        
        elif(isinstance(order, int) and order<=self.cardinality()):
            return [gr for gr in self.all_subgroups() if gr.is_normalSubgroup(self) \
                     if gr.cardinality()==order]
        else:
            raise TypeError("Incorrect order value")

            
    
    def all_subgroups(self, order=None):        
        l = [] 
        for a in self.Set.subsets():
            #print(a)
            try:
                #print(a)
                gr = Group(a, self.bin_op, identity=self.e)
                if(gr.is_subgroup(self)):
                    l.append(gr)
            except:
                #print(a, "is not a subgroup")
                pass
        if(order==None or order=="all"):
            return l
        elif(isinstance(order, int) and order<=self.cardinality()):
            #return group of order 'order'
            return [gr for gr in l if gr.cardinality()==order]
        else:
            raise TypeError("Incorrect order value")
    
    
    
    def is_simple(self):
        return len(self.all_subgroups()==2)
    
    
    
    def is_cyclic(self):
        """Checks if self is a cyclic Group"""
        return any(g.order() == self.cardinality() for g in self)
            
            


if __name__ == '__main__':
    
    #Grupo Z_12
    S=Set(range(12))
    b_op12=Function(S*S, S,lambda x: (x[0]+x[1])%12)
    Z12 = Group(S, b_op12)
    
    #Grupo Z_6
    S=Set(range(6))
    b_op6=Function(S*S, S,lambda x: (x[0]+x[1])%6)
    Z6 = Group(S, b_op6)   
    
    
    #Grupo Z_5
    R=Set(range(5))
    b_op5=Function(R*R, R,lambda x: (x[0]+x[1])%5)
    Z5 = Group(R, b_op5)
    
    
    #Grupo trivial 
    C=Set({0})
    b_op0=Function(C*C, C,lambda x: (x[0]+x[1])%6) 
    Z0 = Group(C, b_op0)
        
    #0,3 
    D=Set({0,3})
    b_op3=Function(D*D, D,lambda x: (x[0]+x[1])%6) 
    Z3 = Group(D, b_op3)
    
    
    print(Z12.all_normalSubgroups(13))
    print(" ")
    #print(Z3.is_cyclic())
    #print(Z5.all_subgroups())
    print(" ")
    #print(Z6.all_normalSubgroups2())
    

    
    '''
    print("Subgrupos de Z12: ", Z12.all_subgroups())
    print("Subgrupos de Z6: ", Z6.all_subgroups())
    print("Subgrupos de Z5: ", Z5.all_subgroups())
    
    
    print(Z0.is_subgroup(Z6))
    print(Z0.is_subgroup(Z5))
    print(Z6.Cayley_table())
    print(Z3.is_subgroup(Z6))
    '''
    
    #print(G.all_subgroups())
    #print(G.is_subgroup(H)) #Veamos si H es subgrupo
    #print(G.is_cyclic())
    #print(G.cardinality())
    
    
    
    '''
    #print(elem4.order())
    #print(G.inverse(elem3))
    #print(elem3.inverse())
    
    
    elem2 = GroupElem(2,G)
    elem3 = GroupElem(3,G)
    elem4 = GroupElem(4,G)

    print(elem0.conjugacy_class())
    print(elem1.conjugacy_class())
    print(elem2.conjugacy_class())
    print(elem3.conjugacy_class())
    print(elem4.conjugacy_class())
    '''

    #print(conjugacy_class)
    #print(G.group_elems)
    #print(G)
    #print(len(G)) #It calls __len__() method
    #print(G.is_abelian())
    #print(G.Cayley_table())
    
    
    
    
    
    
    
    
    
    
    
    
    '''
    def _find_identity(self):
    		for element in self.set:
    			if all(element*x == x for x in self.set):
    				return element
    		raise TypeError('The given set has no identity element on the function')



    def Divisores(n):
        return [x for x in range (1,n+1) if n%x==0]
    
    def is_prime(x):
        if x<2:
            return False
        else:
            for n in range(2,x):
                if x%n == 0:
                    return False
            return True
    
    
        
    def PrimosMenores(n):
        return [x for x in range (1,n+1) if is_prime(x)]
    
    c={1,2,3,4}
    
    #unión
    a=c.union({5})
    
    #intersección
    b=c.intersection({2})
    
    v=all(x%2==0 for x in c)
    
    Divisores(6)
    PrimosMenores(8)
    
    isinstance(4,int)
    type(5)
    
    def isnatural(n):
        if not isinstance(n,int):
            return False
        return n>=0
    
    u = set(range(8))
    rl = set((a,b) for a in u for b in u if (a-b)%5 ==0)#
    
    #aquí para ver diagramas https://github.com/pedritomelenas/LMD/blob/master/Relaciones%20y%20Algebras%20de%20Boole/relaciones.py
    '''