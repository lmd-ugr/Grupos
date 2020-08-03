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
from Permutation import permutation




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

    '''
    def __pow__(self, n, modulo):
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
    '''
    
    def __pow__(self, n):
        #Output: self**n
        g = self
        x = self.group.identity()
        if n%2==1:
            x = self*x
        while(n>1):
            g = g*g
            n = n//2
            if n%2 == 1:
                x = x*g
        return x
    
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
        #Returns the conjugacy class of self in self.group
        return Set([g*self*g**-1 for g in self.group])
    
    '''
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
    '''
    
    def power(self, n):
        #Output: self**n
        g = self
        x = self.group.identity()
        if n%2==1:
            x = self*x
        while(n>1):
            g = g*g
            n = n//2
            if n%2 == 1:
                x = x*g
        return x
    
    
    def order(self, n):
        if not self in self.group:
            raise ValueError
        ident = self.group.identity()

        if n==1:
            return 1
        div = [x for x in range(1, n+1) if n%x==0 if is_prime(x)]
    
        for p in div:
            #if(self.power(n//p) == ident):
            if(self**(n//p) == ident):

                return self.order(n//p)
        
        return n
    
        
    
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

    
    
    
    
def is_prime(x):
    if x<2:
        return False
    else:
        for n in range(2,x):
            if x%n == 0:
                return False
        return True

    
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
        
        '''
        for a in self.Set:
            for b in self.Set:
                table[hash_pair(a,b)] = self.bin_op((a,b))
        return table
        '''
        
        
    def cardinality(self):
        return self.Set.cardinality()
    
    def order(self):
        #Return the order of the group.
        
        if self.group_order != None:
            return self.group_order
        self.group_order = len(self)
        return self.group_order
    
    
    def elements_order(self):
        dev = {}
        
        for a in self.group_elems:
            dev[a] = a.order(self.order())
        return dev
    
    
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
        return self.order()/other.order()
    
    
    def is_normalSubgroup(self, other):
        if not isinstance(other, Group):
            return False
        

        if not self.is_subgroup(other):
            return False
        
        if other.is_abelian() or other.index(self)==2:
            return True 
            
        for g in other:
            for h in self:
                #print(g, "*", h, "*", g.inverse(), " in ", Set(self.group_elems))
                #p = self.bin_op((self.bin_op((g.elem,h.elem)), g.inverse().elem))
                p= g.elem*h.elem*(g**-1).elem
                if not p in self.Set:
                    return False
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
            
            
            try:
                #gr = Group(a, self.bin_op, identity=self.e)
                gr = Group(a, self.bin_op)
                #print(gr)
                if(gr.is_subgroup(self)):
                    print(gr)
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
        return len(self.all_subgroups())==2
    
    
    
    def is_cyclic(self):
        """Checks if self is a cyclic Group"""
        return any(g.order() == self.cardinality() for g in self)
            
            

    def gens_cyclic_group(self):
        
        '''
        #phi Euler
        x = self.order()
        
        #if the group is cyclic
        if self.Set == Set(range(x)):
            if x == 1:
                return 1
            else:
                n = [y for y in range(1,x) if math.gcd(x,y)==1]
                return n
        else:
        ''' 
        l = []
        for a in self.group_elems:
            p = Set(a**h for h in range(0,self.order()))
                
            if(p.cardinality() == self.order()):
                print(a,"gens", p)
                l.append(a)
        if len(l) != 0:
            return l
        else:
            for a in self.group_elems:
                for b in self.group_elems:
                    p = Set((a**h) * (b**j)  for h in range(0,self.order()) for j in range(0,self.order()))
                        
                    if(p.cardinality() == self.order()):
                        print(a,",", b, "gens", p)
                        ab = str(a)+" "+str(b)
                        l.append(ab)
            return l
        
    





def SymmetricGroup(n):
    """
    Returns the symmetric group of order n!
    Example:
        >>> S3=SymmetricGroup(3)
        >>> S3.group_elems
        Set([ (2, 3),  (1, 3),  (1, 2),  (1, 3, 2), ( ),  (1, 2, 3)])
    """
    
    G = Set(permutation(list(g)) for g in itertools.permutations(list(range(1,n+1))))
    bin_op = Function(G.cartesian(G), G, lambda x: x[0]*x[1])
    
    if n>2:
        Gr = Group(G, bin_op, identity=permutation(list(range(1,n+1))), 
        group_order=math.factorial(n), group_degree=n)
        Gr.group_gens=[Gr(permutation([tuple(range(1,n+1))])),Gr(permutation((1,2)).extend(n))]
        return Gr
    if n==2:
        Gr = Group(G, bin_op, identity=permutation(list(range(1,3))), 
        group_order=2, group_degree=2)
        Gr.group_gens=[Gr(permutation([tuple(range(1,3))]))]
    if n==1:
        Gr = Group(G, bin_op, identity=permutation(list(range(1,2))), 
        group_order=1, group_degree=1)
        Gr.group_gens=[Gr(permutation([tuple(range(1,2))]))]
    return Gr
    
    
    
def AlternatingGroup(n):
    """
    Returns the alternating group: the subgroup of even permutations of SymmetricGroup(n)
    
    Example:
        >>> A3=AlternatingGroup(3)
        >>> A3<=S3
        True
        >>> A3.is_normal_subgroup(S3)
        True
        >>> Q=S3/A3
        >>> Q.Set
        Set([Set([ (2, 3),  (1, 2),  (1, 3)]), Set([ (1, 2, 3),  (1, 3, 2), ( )])])
    """
    
    #G = Set(permutation(list(g)) for g in itertools.permutations(list(range(1,n+1))) if permutation(list(g)).sign()==1)
    G = Set(permutation(list(g)) for g in itertools.permutations(list(range(1,n+1))) if permutation(list(g)).even_permutation())

    bin_op = Function(G.cartesian(G), G, lambda x: x[0]*x[1])
    if n>2:
        Gr=Group(G, bin_op,identity=permutation(list(range(1,n+1))),
        group_order=math.factorial(n)//2, group_degree=n)
        Gr.group_gens=[Gr.parent(permutation((i,i+1,i+2)).extend(n)) for i in range(1,n-1)]
    if n==2:
        Gr = Group(G, bin_op, identity=permutation(list(range(1,3))),
        group_order=1, group_degree=2)
        Gr.group_gens=[Gr.parent(permutation(list(range(1,3))))]
    return Gr





def CyclicGroup(n, rep="integers"):
    """
    Returns the cylic group of order n

    Args:
        n a positive integer
        rep may be either "integers" and then the output is integers mod n, or "permuations" and the output is the subgroup of S_n generated by the cycle (1..n)

    Example:
        >>> CP=CyclicGroup(3,"permutations")
        >>> CP.Set
        Set([ (1, 2, 3),  (1, 3, 2), ( )])
        >>> C=CyclicGroup(3,"integers")
        >>> C.group_elems
        Set([0, 1, 2])
        >>> CP.is_isomorphic(C)
        True
    """
    if rep=="integers":
        G = Set(range(n))
        bin_op = Function(G.cartesian(G), G, lambda x: (x[0] + x[1]) % n)
        Gr= Group(G, bin_op,identity=0, group_order=n)
        if n==1:
            Gr.group_gens=[Gr(0)]
        else:    
            Gr.group_gens=[Gr(1)]
        return Gr
    if rep=="permutations":
        def rotate_left(x, y):
            if len(x) == 0:
                return []
            y = y % len(x)
            return x[y:] + x[:y]

        def cyclic(n):
            gen = list(range(1,n+1))
            for i in range(n):
                yield permutation(gen)
                gen = rotate_left(gen, 1)
        G=Set(cyclic(n))
        #bin_op = Function(G.cartesian(G), G, lambda x: x[0]*x[1])
        bin_op=SymmetricGroup(n).bin_op.new_domains(G.cartesian(G),G,check_well_defined=False)
        Gr = Group(G, bin_op, identity=permutation(list(range(1,n+1))),
        abelian=True, group_order=n, group_degree=n,parent=SymmetricGroup(n))
        Gr.group_gens=[Gr(permutation([tuple(range(1,n+1))]))]
        return Gr
    raise ValueError("The second argument can be 'integers' or 'permutations'")









    
if __name__ == '__main__':
    
    
    C = CyclicGroup(6)
    print(C.elements_order())
    '''
    print(C)
    print(C.all_subgroups())
    print(C.Cayley_table())
    
    print(C.gens_cyclic_group())
    '''
    
    S = SymmetricGroup(3)
    #A = AlternatingGroup(3)
    print(S.elements_order())
    #print(S.Cayley_table())
    #print(S.is_abelian())
    #print(A.is_abelian())
    #print(S.gens_cyclic_group())
    #print(A.is_normalSubgroup(S))
    
    #print(A)
    #print(S)
    #print(A)
    #print(A.all_normalSubgroups())
    
    
    
    #p=permutation(1,3,2,4)
    #q=permutation(1,2,3,4)
    #print(p*q)   
    
    
    #print(p, "vs" , q)
    
        
    #Grupo Z_12
    S=Set(range(12))
    b_op12=Function(S*S, S,lambda x: (x[0]+x[1])%12)
    Z12 = Group(S, b_op12)
    
    #Grupo Z_6
    S=Set(range(6))
    b_op6=Function(S*S, S,lambda x: (x[0]+x[1])%6)
    Z6 = Group(S, b_op6)   
    
    print(" ")
    #print(Z6.gens_cyclic_group())
    
    #Grupo Z_5
    R=Set(range(5))
    b_op5=Function(R*R, R,lambda x: (x[0]+x[1])%5)
    Z5 = Group(R, b_op5)
    
    
    #Grupo trivial 
    C=Set({0})
    b_op0=Function(C*C, C,lambda x: (x[0]+x[1])%6) 
    Z0 = Group(C, b_op0)
        
    
    #0,6
    D=Set({0,6})
    b_op3=Function(D*D, D, lambda x: (x[0]+x[1])%12) 
    O6 = Group(D, b_op3)
    
    #0,4,8
    D=Set({0,4,8})
    b_op3=Function(D*D, D, lambda x: (x[0]+x[1])%12) 
    O48 = Group(D, b_op3)
    
    
    '''
    print(Z12.all_subgroups())
    print(" ")
    print(O6.gens_cyclic_group())
    print(" ")

    print(O6.is_subgroup(Z12))
    print(O48.is_subgroup(Z12))
    '''
    #print(Z12.elements_order())
    #print(phi_euler(12))
    #print(O.gens_cyclic_group())
    #print(Z12.gens_cyclic_group())

    #print(Z6.Cayley_table())
    
    #print(Z12.all_subgroups())

    #print(Z12.all_normalSubgroups())
    #print(" ")
    #print(Z3.is_cyclic())
    #print(Z5.all_subgroups())
    #print(" ")
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