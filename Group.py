# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:18:18 2020
"""

import itertools
import math
from sympy import pi


from Set import Set
from Function import Function
from Permutation import permutation
from beautifultable import BeautifulTable
from Dihedral import Dihedral
from Quaternion import Quaternion
from Complex import Complex, print_roots
from ToddCoxeter import CosetTable


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
    '''
        
    
    def inverse(self):
        """Returns the inverse of elem"""
        
        
        if not self in self.group:
            print(self, "no está en " , self.group)
            raise TypeError("Element isn't a GroupElem in the Group")
        
        #El elemento a debe ser un entero para poder aplicarle
        #la operación binaria. Lo detecta como groupElem
        for a in self.group.elements():
            if self.group.bin_op((self.elem, a.elem)) == self.group.e.elem:
                return a
        raise RuntimeError("Didn't find an inverse for g")

    
    
    
    


    
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


    def __le__(self, other):
        """Checks if self is a subgroup of other"""
        return self.is_subgroup(other)

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
    
    
    def __truediv__(self, other):
        """ Returns the quotient group self / other """
        if not other.is_normal_subgroup(self):
            raise ValueError("other must be a normal subgroup of self")
        G = Set(Set(self.bin_op((g, h)) for h in other.Set) for g in self.Set)

        def multiply_cosets(x):
            h = x[0].pick()
            return Set(self.bin_op((h, g)) for g in x[1])

        Gr= Group(G, Function(G.cartesian(G), G, multiply_cosets) ,group_order=self.order() // other.order())
        #Gr.group_gens=list(Set(g*other for g in self.group_gens))
        #Gr.group_gens = [Set(self.bin_op((g.elem, h)) for h in other.Set) for g in self.group_gens]

        return Gr
        
    def __div__(self, other):
        """ Returns the quotient group self / other """
        if not other.is_normal_subgroup(self):
            raise ValueError("other must be a normal subgroup of self")
        G = Set(Set(self.bin_op((g, h)) for h in other.Set) for g in self.Set)

        def multiply_cosets(x):
            h = x[0].pick()
            return Set(self.bin_op((h, g)) for g in x[1])

        Gr= Group(G, Function(G.cartesian(G), G, multiply_cosets), group_order=self.order() // other.order())
        #Gr.group_gens=list(Set(g*other for g in self.group_gens))
        #Gr.group_gens = [Set(self.bin_op((g.elem, h)) for h in other.Set) for g in self.group_gens]

        return Gr
    
    
    
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
        table = BeautifulTable()
        head=list(self.Set)
        

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
        
        
    
    def order(self):
        #Return the order of the group.
        
        if self.group_order != None:
            return self.group_order
        self.group_order = len(self)
        return self.group_order
    
    
    def elements_order(self):
        dev = {}
        
        for a in self.group_elems:
            #dev[a] = a.order(self.order())
            dev[a] = a.order()
        return dev
    


  
    def index(self, other):
        if not isinstance(other,Group):
            raise TypeError(other, " is not a Group")
        return self.order()/other.order()
    
    
    
    def is_simple(self):
        return len(self.all_subgroups())==2
    
    
    def is_cyclic(self):
        """Checks if self is a cyclic Group"""
        #return any(g.order(self.cardinality()) == self.cardinality() for g in self)
        return any(g.order() == self.order() for g in self)
    
            
    
    
    def is_subgroup(self, other):
        #Return True if all elements of self belong to other.
        
        if not isinstance(other, Group):
            #print("No es un grupo")
            return False
        if other.order() % self.order() != 0:
            #print("La cardinalidad tiene que dividir al orden del grupo")
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
    def is_subgroup(self, other):
            """Return True if all elements of self belong to other.
            """
            if not isinstance(other, Group):
               return False
            if self == other:
                return True
            if self.parent==other:
                return True
            if other.order() % self.order() != 0:
                return False
            return Set(self.group_gens)<=Set(other.group_elems) and \
                all(self.bin_op((a.elem, b.elem)) == other.bin_op((a.elem, b.elem)) for a in self.group_gens for b in self.group_gens)

        
        #¿Incompleto? está comprobando solo que a*b = a*b ???
        #print(self.Set in other.Set.subsets())
    ''' 
    
    
    
    def is_normal_subgroup(self, other):
            """Checks if self is a normal subgroup of other"""
            if not(self.is_subgroup(other)):
               return False
            if other.is_abelian():
                return True
            if other.index(self)==2:
                return True 
            #gens1 = self.group_gens
            #gens2 = other.group_gens
            gens1 = self.generators()
            #print(gens1)
            gens2 = other.generators()
            #print(gens2)
            for g in gens2:
                for h in gens1:
                    p = self.bin_op((self.bin_op((g,h)), g**-1))

                    #p = g * h * g**-1
                    if p.elem not in self.Set:  #cambio a originales, sin GroupElems
                        print(p, " not in " , other.group_elems)
                        return False
             
            return True
    
    
    
    
    def all_subgroups(self):        
        l = [] 
        for a in self.Set.subsets():
            try:
                gr = Group(a, self.bin_op)
                gr.parent=self
                if(gr.is_subgroup(self)):
                    l.append(gr)
            except:
                pass
        
        orders= [i.order() for i in l]
        return dict(zip(orders,l)) #devolvemos diccionario para mantener estructura

    

    def subgroups_of_orderN(self, order):

        if(isinstance(order, int)  and self.order()%order==0):
            subgr = self.all_subgroups()
            return [subgr[i] for i in subgr.keys() if i==order]
        else:
            raise TypeError("Incorrect order value")
    

    
    def all_normal_subgroups(self):

        subgr = self.all_subgroups()
        good =  [subgr[h] for h in subgr if subgr[h].is_normal_subgroup(self)]
        
        orders= [h.order() for h in good]
        return dict(zip(orders,good)) #devolvemos diccionario para mantener estructura




    def cartesian(self, other):
            
        """
        Returns the cartesian product of the two groups
        """
        if not isinstance(other, Group):
            raise TypeError("other must be a group")

        G = (self.Set).cartesian(other.Set)
        
        bin_op = Function( G.cartesian(G), G, 
            lambda x: (self.bin_op((x[0][0], x[1][0])),  other.bin_op((x[0][1], x[1][1]))) )
        
        #Gr=Group(G, bin_op, identity=(self.e.elem, other.e.elem))
        Gr=Group(G, bin_op)
        
        Gr.group_gens=[Gr((self.e.elem,b.elem))  for b in other.group_gens]+[Gr((a.elem,other.e.elem)) for a in self.group_gens]
        return Gr
    
    direct_product = cartesian

    

    def is_direct_product(self, H, K):
        
        if not isinstance(H, Group) or not isinstance(K, Group):
            raise TypeError("h and k must be groups")
        
        if not (H.intersection(K).is_trivial()):
            return False
        
        if not H.is_normal_subgroup(self) or not K.is_normal_subgroup(self):
            return False
            
        if not (H.direct_product(K)).is_isomorphic(self):
            print(H.direct_product(K), " not = ", self)
            return False
    
        return True

        '''
        for h in H.Set:#other:
            for k in K.Set:
            #for h in self:
                #print(g, "*", h, "*", g.inverse(), " in ", Set(self.group_elems))
                #if not self.bin_op((h,k)) == self.bin_op((k,h)):
                #print(p , " = " , self.bin_op((g.elem,h.elem)), "*", g.inverse().elem() )
                ##if k.elem*h.elem != h.elem*k.elem:
                izq = self.bin_op(( (h[1],k[1]), (h[0],k[0])  ))
                dcha= self.bin_op(( (h[0],k[0]), (h[1],k[1])  ))
                #if not self.bin_op(( (h[0],k[0]), (h[1],k[1]) ))==self.bin_op(( (k[0],h[0]), (k[1],h[1])  )) :
                if not izq == dcha :
                    print("{}*{}={} vs {}".format(h,k,izq,dcha))
                    print("{} != {}".format(izq,dcha))
                    #return False
                #else:
                 #   print("{}*{}={} vs {}".format(h,k,izq,dcha))
                  #  print("{} == {}".format(izq,dcha))
        '''
        
    

    def is_inner_direct_product(self):
 
        subgr = self.all_subgroups().values()
        print(subgr)
        for h in subgr:
            for k in subgr:
                print("comparo: " , h , " con " , k)
                if self.is_direct_product(h, k):
                    return True
        return False
    

    
   
    def generate(self, generators):
            """
            Returns the subgroup of self generated by the list generators of GroupElems
            If any of the items aren't already GroupElems, we will try to convert
            them to GroupElems before continuing.
            generators must be iterable
            """
            idn = self(self.e.elem)
            ord = 1
            element_list = [idn]
            set_element_list = set([idn])
            gens = [g if isinstance(g, GroupElem) else GroupElem(g, self) for g in generators]
            for i in range(len(gens)):
                # D elements of the subgroup G_i generated by gens[:i]
                D = element_list[:]
                N = [idn]
                while N:
                    A = N
                    N = []
                    for a in A:
                        for g in gens[:i + 1]:
                            ag = a*g
                            if ag not in set_element_list:
                            # produce G_i*g
                                for d in D:
                                    ord += 1
                                    ap = d*ag
                                    element_list.append(ap)
                                    set_element_list.add(ap)
                                    N.append(ap)
            G = Set(g.elem for g in element_list)
            #Gr=Group(G, self.bin_op.new_domains(G.cartesian(G), G, check_well_defined=False),
                         #parent=self.parent,check_ass=False,check_inv=False, identity=self.e.elem,group_order=ord)
            Gr=Group(G, self.bin_op.new_domains(G.cartesian(G), G))
                         
            Gr.group_gens=gens
            return Gr
    
    
    
    
    
    def generators(self):
        """
        Returns a list of GroupElems that generate self, with length
        at most log_2(len(self)) + 1
        """
    
        if len(self.group_gens) != len(self):
            return self.group_gens
        result = [self.e.elem]
        H = self.generate(result)

        while len(H) < len(self):
            result.append(next(iter(self.Set - H.Set)))
            H = self.generate(result)

        # The identity is always a redundant generator in nontrivial Groups
        if len(self) != 1:
            result = result[1:]

        self.group_gens= [GroupElem(g, self) for g in result]
        return [GroupElem(g, self) for g in result]
        
        

    def subgroup_search(self, prop):
        """Return a subgroup of all elements satisfying the property."""
        H=set([self(self.e.elem)])
        G=self.group_elems
        F=G-H
        while len(F)>0:
            g=Set(F).pick()
            if prop(g):
                H=set(self.generate(list(H|set([g]))))
                F=F-H
            else:
                Dg=Set([g*h for h in H])
                F=Set(F-Dg)
        G = Set(g.elem for g in H)
        Gr=Group(G, self.bin_op.new_domains(G.cartesian(G), G, check_well_defined=False),
                     parent=self.parent,check_ass=False,check_inv=False, identity=self.e.elem,group_order=len(G))
        return Gr

    def element_search(self, prop):
        """Return an element satisfying the property."""
        H=set([self(self.e.elem)])
        G=self.group_elems
        F=G-H
        while len(F)>0:
            g=Set(F).pick()
            if prop(g):
                return g
            else:
                Dg=Set([g*h for h in H])
                F=Set(F-Dg)
        return g

    def normal_closure(self, other):
        """Return the normal closure of other in self."""
        if not(other<=self):
            raise TypeError("The argument must be subgroup of the given group")
        G=self.group_gens
        U=other.group_gens
        N=[]
        for x in U:
            N.append(x)
        for d in N:
            for g in G:
                c=g*d*g**-1
                if not c in self.generate(N):
                    N.append(c)
        return self.generate(N)
    
    
        
    def Sylow_subgroups(self, p):
        n=len(self)
        if n%p!=0:
           return []
        r=1
        while n%p ==0:
            r=r*p
            n=n/p
        return self.subgroups_of_orderN(r)

    def is_transitive(self):
        f=GroupAction(self,Set(list(range(1,self.group_degree+1))),lambda x,y: x.elem(y))
        return f.is_transitive()

        
        
        
    def find_isomorphism(self, other):
        """
        Returns an isomorphic GroupHomomorphism between self and other,
        or None if self and other are not isomorphic
        Uses Tarjan's algorithm, running in O(n^(log n + O(1))) time, but
        runs a lot faster than that if the group has a small generating set.
        """
        if not isinstance(other, Group):
            raise TypeError("other must be a Group")

        if len(self) != len(other) or self.is_abelian() != other.is_abelian():
            return None

        # Try to match the generators of self with some subset of other
        A = self.generators()
        for B in itertools.permutations(other.group_elems, len(A)):

            func = dict(zip(A, B)) # the mapping
            counterexample = False
            while not counterexample:

                # Loop through the mapped elements so far, trying to extend the
                # mapping or else find a counterexample
                noobs = {}
                for g, h in itertools.product(func, func):
                    if g * h in func:
                        if func[g] * func[h] != func[g * h]:
                            counterexample = True
                            break
                    else:
                        noobs[g * h] = func[g] * func[h]

                # If we've mapped all the elements of self, then it's a
                # homomorphism provided we haven't seen any counterexamples.
                if len(func) == len(self):
                    break

                # Make sure there aren't any collisions before updating
                imagelen = len(set(noobs.values()) | set(func.values()))
                if imagelen != len(noobs) + len(func):
                    counterexample = True
                func.update(noobs)

            if not counterexample:
                return GroupHomomorphism(self, other, lambda x: func[x],check_morphism_axioms=False)

        return None
    
    
    def is_isomorphic(self, other):
        """Checks if self and other are isomorphic"""
        return self.find_isomorphism(other)==True

    def intersection(self, other):
        """
        Computes the intersection of self and other
        """
        if self.parent==None or other.parent==None:
            raise TypeError("self and other must be subgroups of the same Group")
        if self.parent != other.parent:
            raise TypeError("self and other must be subgroups of the same Group")
        common = Set(self.Set & other.Set)
        return Group(common,Function(common.cartesian(common), common, self.bin_op.function), self.parent)



    def conjugacy_classes(self):
        """Compute the set of conjugacy clases of the elements of a group; see conjugacy_class of a group element"""
        cls=set([])
        G=list(self.group_elems)
        while len(G)>0:
            x=G[0]
            H=x.conjugacy_class()
            cls.add(H)
            G=[g for g in G if not(g in H)]
        return cls

    def conjugate_subgroup(self,other,g):
        conj=[g*x*g**-1 for x in other.group_gens]
        return self.generate(conj)

    def conjugacy_class_subgroup(self,other):
        """Computes the set of g*other*g^-1 for all g in self"""
        cls=Set([self.conjugate_subgroup(other,g) for g in self.group_elems])
        return cls

    def conjugacy_classes_subgroups(self):
        """Computes the set of g*H*g^-1 for all H subgroup of self"""
        sbs_d=self.subgroups()
        sbs=[H for j in sbs_d.keys() for H in sbs_d[j]]
        c = set([self.conjugacy_class_subgroup(H) for H in sbs])
        return c 

    def cosets(self,other,side="left"):
        if not(other<=self):
            raise ValueError("The first argument must be a subgroup of the given group")
        if side=="right":
            cls=[set(other.group_elems)]
            D=list((self.group_elems)-(other.group_elems))
            while len(D)>0:
                a=D[0]
                new=set(a*h for h in other)
                cls.append(new)
                D=list(set(D)-new)
            return cls
        if side=="left":
            cls=[set(other.group_elems)]
            D=list((self.group_elems)-(other.group_elems))
            while len(D)>0:
                a=D[0]
                new=set(h*a for h in other)
                cls.append(new)
                D=list(set(D)-new)
            return cls
        raise ValueError("The second argument must be either 'left' or 'right'")




    def normalizer(self,other):
        S=other.group_gens
        H=other.group_elems
        def prop(g):
            return all((g*h*g**-1 in H) for h in S)
        return self.subgroup_search(prop)


    def commutator(self,H,K):
        ggens = H.group_gens
        hgens = K.group_gens
        commutators = []
        for ggen in ggens:
            for hgen in hgens:
                commutator = hgen * ggen * hgen**-1 * ggen**-1
                if commutator not in commutators:
                    commutators.append(commutator)
        res = self.normal_closure(self.generate(commutators))
        return res

    def centralizer(self,other):
        def prop(g):
            return [g*h for h in other.group_gens]  ==  [h*g for h in other.group_gens]
        return self.subgroup_search(prop)


    def center(self):
        """Computes the center of self: the subgroup of element g such that g*h=h*g for all h in G"""
        return self.centralizer(self)


    def derived_subgroup(self):
        """Return the derived subgroup of the group."""
        return self.commutator(self, self)


    def derived_series(self):
        res = [self]
        current = self
        nextgroup = self.derived_subgroup()
        while not current<=nextgroup:
            res.append(nextgroup)
            current = nextgroup
            nextgroup = nextgroup.derived_subgroup()
        return res


    def is_trivial(self):
        return self.order()== 1 and all(g.elem == self.e.elem for g in self)


    def is_soluble(self):
        ds = self.derived_series()
        terminator = ds[len(ds) - 1]
        return terminator.is_trivial()


    def group_lattice(self):
        from graphviz import Digraph
        from IPython.display import display, Image,HTML
        G=Digraph(graph_attr={'rankdir': 'TB'},node_attr={'rank':'source','shape':'plaintext'},format='png')
        sbs_d=self.subgroups()
        sbs=[H for j in sbs_d.keys() for H in sbs_d[j]]
        letters = "ABCDEFGHIJKLNNOPQRSTUVWXYZ"
        toletter = {}
        toelem = {}
        for letter, elem in zip(letters, sbs):
            toletter[elem] = letter
            toelem[letter] = elem
        letters = letters[:len(sbs)]
        table=' <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4"> <TR><TD COLSPAN="2"><B>Subgroups</B></TD></TR>'
        for g in sbs:
            table=table+"<tr><td>"+toletter[g]+"</td><td>"+str(set(g.group_elems))+"</td></tr>"
        table=table+"</TABLE> "
        for g in sbs:
            G.node(repr(g),label=toletter[g],shape='circle',color='Cyan',style='filled')
        edges=[(a,b) for a in sbs for b in sbs if (b <= a) and (a!=b)]
        for a in sbs:
            for b in sbs:
                for c in sbs:
                    if ((a,b) in edges) and ((b,c) in edges and ((a,c))in edges):
                        edges.remove((a,c)) 
        for e in edges:
            G.edge(repr(e[0]),repr(e[1]),color='Red',dir='back')    
        png_str = G.render(cleanup=True)
        display(HTML(table),Image(png_str))

    def CayleyGraph(self):
        from IPython.display import Image,SVG,display
        import graphviz as gv
        G=gv.Digraph(format='png',engine='circo')
        elems=[g for g in self.group_gens]
        vertices=[repr(a) for a in self]
        for v in vertices:
            G.node(v)
        letters = "abcdefghijklmnopqrstuvwxyz"
        colors = ['Red','Blue','Green','Yellow','SkyBlue','PaleGreen','Lime','Tan','YellowGreen',
              'Violet','Tomato','Plum','Magenta','Mocassin','LightPink', 'LightCyan','Khaki',
              'Ivory','Gray','Fuchsia','Coral','Bisque','Aqua','Azure','DeepPink','Gainsboro',
              'HotPink']
        toletter = {}
        toelem = {}
        for letter, elem in zip(letters, elems):
            toletter[elem] = letter
            toelem[letter] = elem
        letters = letters[:len(elems)]
        tocelem = {}
        tocolor = {}
        for elem, color in zip(elems,colors):
            tocolor[elem] = color
            tocelem[color] = elem
        colors = colors[:len(elems)]
        E={}
        for i in range(len(elems)):
            E[i]=[(repr(a),repr(b)) for a in self for b in self if b==a*elems[i]]
            for e in E[i]:
                G.edge(e[0],e[1],toletter[elems[i]],_attributes={'color':tocolor[elems[i]]})
        png_str = G.render(cleanup=True)
        display(Image(data=png_str))


    def automorphism_identity(self):
        return GroupHomomorphism(self, self,lambda x: x,check_morphism_axioms=False)


    def AutomorphismGroup(self):
        """Group of automorphisms of self"""
        dicts=[]
        # Try to match the generators of self with some subset of other
        A = self.group_gens
        for B in itertools.permutations(self.group_elems, len(A)):

            func = dict(zip(A, B)) # the mapping
            counterexample = False
            while not counterexample:

                # Loop through the mapped elements so far, trying to extend the
                # mapping or else find a counterexample
                noobs = {}
                for g, h in itertools.product(func, func):
                    if g * h in func:
                        if func[g] * func[h] != func[g * h]:
                            counterexample = True
                            break
                    else:
                        noobs[g * h] = func[g] * func[h]

                # If we've mapped all the elements of self, then it's a
                # homomorphism provided we haven't seen any counterexamples.
                if len(func) == len(self):
                   break

                # Make sure there aren't any collisions before updating
                imagelen = len(set(noobs.values()) | set(func.values()))
                if imagelen != len(noobs) + len(func):
                    counterexample = True
                func.update(noobs)
            if not counterexample:
                dicts.append(func)
        def fn(d):
            return lambda x:d[x]
        G=Set([GroupHomomorphism(self,self,fn(d)) for d in dicts])
        bin_op= Function(G.cartesian(G), G, lambda x: x[0].homomorphism_compose(x[1]))
        Gr=Group(G, bin_op)
        return Gr

    def automorphism_by_conjugation(self,other):
        return GroupHomomorphism(self, self,lambda x: other*x*other**-1)


    def InnerAutomorphismGroup(self):
        G=Set([self.automorphism_by_conjugation(g) for g in self])
        bin_op= Function(G.cartesian(G), G, lambda x: x[0].homomorphism_compose(x[1]),check_well_defined=False)
        Gr=Group(G, bin_op,check_ass=False,check_inv=False,identity=self.automorphism_identity(), parent=self.AutomorphismGroup())
        return Gr
        
    
    def AllHomomorphisms(self,other):
        """List of homomorphisms of self to other"""
        dicts=[]
        # Try to match the generators of self with some subset of other
        A = self.group_gens
        for B in itertools.permutations(other.group_elems, len(A)):
            func = dict(zip(A, B)) # the mapping
            counterexample = False
            while not counterexample:
                # Loop through the mapped elements so far, trying to extend the
                # mapping or else find a counterexample
                noobs = {}
                for g, h in itertools.product(func, func):
                    if g * h in func:
                        if func[g] * func[h] != func[g * h]:
                            counterexample = True
                            break
                    else:
                        noobs[g * h] = func[g] * func[h]
                # If we've mapped all the elements of self, then it's a
                # homomorphism provided we haven't seen any counterexamples.
                if len(func) == len(self):
                    break
                func.update(noobs)
            if not counterexample:
                dicts.append(func)
        def fn(d):
            return lambda x:d[x]
        return [GroupHomomorphism(self,other,fn(d)) for d in dicts]



    def semidirect_product(self, other,hom):
        """
        Returns the semidirect product of the two groups
        """
        if not isinstance(other, Group):
            raise TypeError("other must be a group")
        bin_op=Function(self.Set.cartesian(other.Set).cartesian(self.Set.cartesian(other.Set)),self.Set.cartesian(other.Set),
                        lambda x:(self.bin_op((x[0][0],hom(other(x[0][1])).elem.function(self(x[1][0])).elem)),\
                                  other.bin_op((x[0][1], x[1][1]))), check_well_defined=False)
            
        Gr=Group((self.Set).cartesian(other.Set), bin_op, group_order=self.group_order*other.group_order)
        #Gr.group_gens=[(self.e.elem,b)  for b in other.group_gens]+[(a,other.e.elem)  for a in self.group_gens]
        return Gr
  

    def is_inner_semidirect_product(self):


        subgr = self.all_subgroups().values()
        print(subgr)
        for h in subgr:
            for k in subgr:
                if h.is_normal_subgroup(self):
                    if h.intersection(k).is_trivial():
                        hom=GroupHomomorphism(k,h.AutomorphismGroup(),lambda x: x,check_morphism_axioms=False)        
                        if h.semidirect_product(k, hom).is_isomorphic(self):
                            return True
        return False





class GroupHomomorphism(Function): #we should add here check_well_defined, and check_group_axioms as options
    """
    The definition of a Group Homomorphism
    A GroupHomomorphism is a Function between Groups that obeys the group
    homomorphism axioms.
    Args:
        domain, codomain and function. domain and codomain are groups, and function is the map
    Example:
        >>> G=CyclicGroup(2)
        >>> H=G.cartesian(G)
        >>> f=GroupHomomorphism(H,G, lambda x:G(x.elem[1]))
        >>> f.kernel()
        Group with 2 elements
        >>> f.is_surjective()
        True
    """

    def __init__(self, domain, codomain, function, check_morphism_axioms=True):
        """Check types and the homomorphism axioms; records the two groups"""

        if not isinstance(domain, Group):
            raise TypeError("domain must be a Group")
        if not isinstance(codomain, Group):
            raise TypeError("codomain must be a Group")
        #if not all(function(elem) in codomain for elem in domain):
        #    raise TypeError("Function returns some value outside of codomain")

        if check_morphism_axioms:
            if not all(function(elem) in codomain for elem in domain):
                raise TypeError("Function returns some value outside of codomain")
            if not all(function(a * b) == function(a) * function(b) \
                       for a in domain for b in domain):
                raise ValueError("function doesn't satisfy the homomorphism axioms")

        self.domain = domain
        self.codomain = codomain
        self.function = function

    def __call__(self,other):
        return self.function(other)

    def __hash__(self):
        return hash(self.domain) + 2 * hash(self.codomain)

    def __eq__(self, other):
        if not isinstance(other, GroupHomomorphism):
            return False

        return self.domain == other.domain and self.codomain==other.codomain and all(self.function(a)==other.function(a) for a in self.domain)

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        if not(self.domain==self.codomain):
            return "Group homomorphism between "+str(self.domain)+" and "+str(self.codomain)
        return "Group endomorphism of "+str(self.domain)

    def __repr__(self):
        if not(self.domain==self.codomain):
            return "Group homomorphism"
        return "Group endomorphism"

    def kernel(self):
        """Returns the kernel of the homomorphism as a Group object"""
        G = Set(g.elem for g in self.domain if self(g) == self.codomain.e)
        return Group(G, self.domain.bin_op.new_domains(G.cartesian(G), G, check_well_defined=False))

        #return Group(G, self.domain.bin_op.new_domains(G.cartesian(G), G, check_well_defined=False),parent=self.domain, check_ass=False,check_inv=False,identity=self.domain.e.elem)

    def image(self):
        """Returns the image of the homomorphism as a Group object"""
        G = Set(g.elem for g in self._image())
        return Group(G, self.codomain.bin_op.new_domains(G.cartesian(G), G, check_well_defined=False))

        #return Group(G, self.codomain.bin_op.new_domains(G.cartesian(G), G, check_well_defined=False),parent=self.codomain, check_ass=False,check_inv=False,identity=self.codomain.e.elem)

    def is_isomorphism(self):
        return self.is_bijective()

    def homomorphism_compose(self,other):
        if not self.domain == other.codomain:
            raise ValueError("codomain of other must match domain of self")
        return GroupHomomorphism(other.domain, self.codomain,lambda x: self(other(x)), check_morphism_axioms=False)

    def automorphism_inverse(self):
        if not self.function.is_bijective():
            raise ValueError("self must be bijective")
        l={}
        for x in self.domain:
            l[x]=self(x)
        inv = {v: k for k, v in l.items()}
        return GroupHomomorphism(self.codomain, self.domain,lambda x: inv[x], check_morphism_axioms=False)






class GroupAction: #we should add here check_well_defined, and check_group_axioms as options
    """
    The definition of a Group Action
    A Group Action is a Function from the cartasian product of a set X and a group G to X, fulfilling some properties
    Example:
        >>> from Group import *
        >>> G=SymmetricGroup(3)
        >>> f=GroupAction(G,Set({1,2,3}),lambda x,y:x.elem(y))
        >>> p=G(permutation(2,3,1))
        >>> f.action(p,3)
        1
        >>> f.orbit(2)
        frozenset({1, 2, 3})
        >>> f.stabilizer(3)
        Group with 2 elements
        >>> list(_)
        [( ),  (1, 2)]
    """

    def __init__(self, group, theset, action, check_action_axioms=True):
        """Check types and the homomorphism axioms; records the two groups"""

        if not isinstance(group, Group):
            raise TypeError("The first argument must be a Group")
        if not isinstance(theset, Set):
            raise TypeError("The secon argument must be a Set")
        # f(g,x) maps to X with g in the group and x in the set
        #if not all(action(g,x) in theset for g in group for x in theset):
        #    raise TypeError("action returns some value outside of codomain")

        if check_action_axioms:
            if not all(action(g,x) in theset for g in group for x in theset):
                raise TypeError("action returns some value outside of codomain")
            #first axiom a*(b*x)=(a b)*x
            if not all(action(a,action(b,x)) == action(a*b,x) \
                       for a in group for b in group for x in theset):
                raise ValueError("action doesn't satisfy the first action axiom")
            #second axiom
            if not all(action(group.e, x)==x for x in theset):
                raise ValueError("action doesn't satisfy the second action axiom")

        self.group = group
        self.set = theset
        self.action = action

    def __str__(self):
        return "Group action"

    def __repr__(self):
        return "Group action from ("+str(self.group)+")x("+str(self.set)+") to "+str(self.set)

    def __call__(self,g,el):
        return self.action(g,el)


    def orbit(self, other):
        if not(other in self.set):
            raise ValueError("other must be in self.set")
        orb=[other]
        S=self.group.group_gens
        for y in orb:
            for a in S:
                if not self.action(a,y) in orb:
                    orb.append(self.action(a,y)) 
        return Set(orb)
        
    def orbits(self):
        lels=list(self.set)
        lorb=[]
        while len(lels)>0:
            el=lels[0]
            orb=self.orbit(el)
            lorb.append(orb)
            lels=[g for g in lels if not(g in orb)]
        return lorb
    
    def stabilizer(self,other):
        def prop(x):
            return self.action(x,other)  ==  other
        return self.group.subgroup_search(prop)

    def is_transitive(self):
        got_orb = False
        for x in self.orbits():
            if len(x) > 1:
                if got_orb:
                    return False
                got_orb = True
        return got_orb




#This is used for todd coxeter, dont delete
def generate(elems):

    oldG = Set(elems)
    
    while True:
        newG = oldG | Set(a * b for a in oldG for b in oldG)
        if oldG == newG: 
            break
        else: 
            oldG = newG
    C = Set(oldG)
    bin_op = Function(C*C, C, lambda x: x[0]*x[1])
    
    G = Group(C, bin_op)
    #G.group_gens = [ GroupElem(g, G) for g in elems ]
    G.group_gens = elems
    return G   
    


















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
        Gr = Group(G, bin_op ,identity=permutation(list(range(1,3)))
        ,group_order=1, group_degree=2)
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
        
        Gr.group_gens = [Gr(0)] if n==1 else [Gr(1)]
        '''
        if n==1:
            Gr.group_gens=[Gr(0)]
        else:    
            Gr.group_gens=[Gr(1)]
        '''
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





def DihedralGroup(n, rep="RS"):
        
    
    """
    Returns the dihedral group of order 2n
    Args:
        n is a positive integer
        rep can be "RS" if we want rotations and symmtries, or "permutations" if we want to see DihedralGroup(n) inside SymmetricGroup(n)
    Example:
        >>> D3=DihedralGroup(3)
        >>> DP3=DihedralGroup(3,"permutations")
        >>> D3.is_isomorphic(DP3)
        True
        >>> D3.Set
        Set(['R0', 'R1', 'R2', 'S2', 'S1', 'S0'])
        >>> DP3.Set
        Set([ (2, 3),  (1, 3),  (1, 2),  (1, 3, 2), ( ),  (1, 2, 3)])
    """
        
    
        
    if rep=="matrix":
        D = Dihedral(n)
        rots = [tuple(x) for x in D.rot]
        refls =[tuple(x) for x in D.refl]

        mlist = rots + refls
            
        G = Set(mlist)
        Gr=Group(G, Function(G.cartesian(G), G, lambda x: D.product_Matrix(x[0],x[1])) ,group_order=2*n)
            
        #Gens:= r, s
        Gr.group_gens=[Gr(rots[1]),Gr(refls[0])]
        return Gr
        
        
    if rep=="RS":
        D = Dihedral(n)        
        G = Set(D.all.keys())
        Gr=Group(G, Function(G.cartesian(G), G, lambda x: D.product_RS(x[0],x[1])))
        
        #Con esto ya funciona is_isomorphic
        Gr.group_gens=[Gr('R1'),Gr('S0')]
        return Gr
        
        
    if rep=="permutations":
        def rotate_left(x, y):
            if len(x) == 0:
                return []
            y = y % len(x)
            return x[y:] + x[:y]
    
        def dihedral(n):
            if n == 1:
                yield permutation([1, 2])
                yield permutation([2, 1])
            elif n == 2:
                yield permutation([1, 2, 3, 4])
                yield permutation([2, 1, 4, 3])
                yield permutation([3, 4, 1, 2])
                yield permutation([4, 3, 2, 1])
            else:
                gen = list(range(1,n+1))
                for i in range(n):
                    yield permutation(gen)
                    yield permutation(gen[::-1])
                    gen = rotate_left(gen, 1)

        G=Set(dihedral(n))
        bin_op = Function(G.cartesian(G), G, lambda x: x[0]*x[1])
        Gr=Group(G, bin_op,identity=permutation(list(range(1,2*n+1))),
                 group_order=2*n,parent=SymmetricGroup(n), group_degree=n)
        Gr.group_gens=[Gr.parent(permutation([1]+list(range(2,n+1))[::-1])),Gr.parent(permutation([tuple(range(1,n+1))]))]
        
        return Gr
        
    raise ValueError("The second argument can be 'matrix' , 'RS' or 'permutations'")
    



def QuaternionGroup(rep="ijk"):
    """
    Example:
        >>> Q = Quaternion.Group(rep="ijk")
        >>> print(Q)
        'Group with 8 elements: { 1,  i,  j,  k,  -k,  -j,  -i,  -1}'
        
        >>> i = Quaternion(letter="i")
        >>> j = Quaternion(letter="j")
        >>> k = Quaternion(letter="j")
            
        >>> print(i*j)
        'k'
        >>> print(i*j*k)
        '-1'
        >>> print(i*i == j*j == k*k == i*j*k == -1)
        'True' 
    """
    if rep=="ijk":
        one = Quaternion(1,0,0,0)
        i   = Quaternion(0,1,0,0)
        j   = Quaternion(0,0,1,0)
        k   = Quaternion(0,0,0,1)
        
        q2=[one, -one, i, -i, j, -j, k, -k]
        
        G=Set(q2)
        Gr=Group(G,Function(G.cartesian(G),G, lambda x: x[0]*x[1]))
        Gr.group_gens = [Gr(i),Gr(j)]
        return Gr
    
    
    if rep=="permutations":
        q1=[permutation([1, 2, 3, 4, 5, 6, 7, 8]), permutation([2, 3, 4, 1, 6, 8, 5, 7]),
            permutation([3, 4, 1, 2, 8, 7, 6, 5]), permutation([4, 1, 2, 3, 7, 5, 8, 6]),
            permutation([5, 7, 8, 6, 3, 2, 4, 1]), permutation([6, 5, 7, 8, 4, 3, 1, 2]),
            permutation([7, 8, 6, 5, 2, 1, 3, 4]), permutation([8, 6, 5, 7, 1, 4, 2, 3])]
        G=Set(q1)
        bin_op = Function(G.cartesian(G), G, lambda x: x[0]*x[1])
        Gr=Group(G, bin_op)
        Gr.group_gens=[Gr(permutation([4, 1, 2, 3, 7, 5, 8, 6])),Gr(permutation([6, 5, 7, 8, 4, 3, 1, 2]))]
        return Gr
    
    raise ValueError("The second argument must be 'ijk' or 'permutations'")






def RootsOfUnitGroup(n):
    '''
    #G = Set(complex( math.cos(2*pi*k/n), math.sin(2*pi*k/n)) for k in range(n))
    G = Set( (np.round(math.cos(2*pi*k/n),3), np.round(math.sin(2*pi*k/n),3)) for k in range(n))
    
    def mul(x):
        #return complex(x[0]*y[0]-x[1]*y[1], x[1]*y[0]+x[0]*y[1])
        return (np.round(x[0][0]*x[1][0]-x[0][1]*x[1][1],3), np.round(x[0][1]*x[1][0]+x[0][0]*x[1][1],3))
        
    bin_op=Function(G.cartesian(G),G,lambda x: mul(x), check_well_defined=False)
    Gr=Group(G,bin_op) 
    
    return Gr
    '''
    
    G = Set(Complex( math.cos(2*pi*k/n), math.sin(2*pi*k/n)) for k in range(n))

    bin_op=Function(G.cartesian(G),G,lambda x: x[0].product(x[1]) ,  check_well_defined=False)
    Gr=Group(G,bin_op) 
    
    Gr.gens_group()
    Gr.gens_group(G.pick())
    
    return Gr



def QuaternionGroup2(rep="ijk"):
    """
    The quaternion group Q2; its elments are 1,i,j,k and their oposite
    Example:
        >>> Q2=QuaternionGroup()
        >>> list(Q2)
        ['1', 'i', 'k', 'j', '-i', '-k', '-j', '-1']
        >>> i=Q2("i")
        >>> j=Q2("j")
        >>> i*j
        'k'
        >>> j*i
        '-k'
    """
    if rep=="ijk":
        q2=[ "1", "-1", "i", "-i", "j", "-j", "k", "-k"]
        table=[[ "1", "-1", "i", "-i", "j", "-j", "k", "-k"],[ "-1", "1", "-i", "i", "-j", "j", "-k", "k"],[ "i", "-i", "-1", "1", "k", "-k", "-j", "j"],[ "-i", "i", "1", "-1", "-k", "k", "j", "-j"],[ "j", "-j", "-k", "k", "-1", "1", "i", "-i"],[ "-j", "j", "k", "-k", "1", "-1", "-i", "i"],[ "k", "-k", "j", "-j", "-i", "i", "-1", "1"],[ "-k", "k", "-j", "j", "i", "-i", "1", "-1"]]
        def product(a,b):
            i=q2.index(a)
            j=q2.index(b)
            return table[i][j]
        G=Set(q2)
        Gr=Group(G,Function(G.cartesian(G),G, lambda x: product(x[0],x[1])))
        Gr.group_gens=[Gr("i"),Gr("j")]
        return Gr
    if rep=="permutations":
        q1=[permutation([1, 2, 3, 4, 5, 6, 7, 8]), permutation([2, 3, 4, 1, 6, 8, 5, 7]),
            permutation([3, 4, 1, 2, 8, 7, 6, 5]), permutation([4, 1, 2, 3, 7, 5, 8, 6]),
            permutation([5, 7, 8, 6, 3, 2, 4, 1]), permutation([6, 5, 7, 8, 4, 3, 1, 2]),
            permutation([7, 8, 6, 5, 2, 1, 3, 4]), permutation([8, 6, 5, 7, 1, 4, 2, 3])]
        G=Set(q1)
        bin_op = Function(G.cartesian(G), G, lambda x: x[0]*x[1])
        Gr=Group(G, bin_op)
        Gr.group_gens=[Gr(permutation([4, 1, 2, 3, 7, 5, 8, 6])),Gr(permutation([6, 5, 7, 8, 4, 3, 1, 2]))]
        return Gr
    raise ValueError("The second argument must be 'ijk' or 'permutations'")



def QuaternionGroupGeneralised(n):

    #Q_n= < a,b | a^{n}=b^{2},a^{2n}=1,b^{-1}ab=a^{-1} >
    genG = ['a','b']
    
    relsG = ['a'*n + 'BB', 'a'*2*n, 'Baba']
    genH = []
    
    C = CosetTable(genG, relsG, genH)
    C.CosetEnumeration()
    gens = C.getGenerators()
    G = generate(gens)
    return G
    



def KleinGroup(rep="integers"):
    """
    The Klein four-group; it can be represented via "integers" as Z_2 x Z_2, and as a subgroup of permutations of A_4
    Example:
        >>> K=KleinGroup()
        >>> KK=KleinGroup("permutations")
        >>> KK.is_isomorphic(K)
        True
        >>> list(K)
        [(0, 0), (0, 1), (1, 0), (1, 1)]
        >>> list(KK)
        [( ),  (1, 4)(2, 3),  (1, 3)(2, 4),  (1, 2)(3, 4)]
    """

    if rep=="integers":
        G=CyclicGroup(2)
        return G.direct_product(G)
    if rep=="permutations":

        G=Set([permutation([1,2,3,4]),permutation((1,2),(3,4)),permutation((1,3),(2,4)),permutation((1,4),(2,3))])
        #bin_op = Function(G.cartesian(G), G, lambda x: x[0]*x[1])
        bin_op=AlternatingGroup(4).bin_op.new_domains(G.cartesian(G),G,check_well_defined=False)
        Gr=Group(G, bin_op, parent=AlternatingGroup(4), group_order=4, group_degree=4)
        Gr.group_gens=[Gr.parent(permutation((1,2),(3,4))), Gr.parent(permutation((1,3),(2,4)))]
        return Gr
    raise ValueError("The second argument can be 'RS' or 'permutations'")


    
if __name__ == '__main__':
    
    
    

    
    '''
    Q = QuaternionGroupGeneralised(2)
    Q2 = QuaternionGroup()
    print(Q)
    print(Q2.is_isomorphic(Q))
    '''
    
    K = KleinGroup()
    Q = KleinGroup("permutations")
    #print(K)
    #print(Q)
    #print(K.is_isomorphic(Q))
    
    
    Q2=QuaternionGroup2()
    #print(list(Q2))
    #i=Q2("i")
    #j=Q2("j")
    
    
    #print(Q2.generate([i,j]))
    
    Q = QuaternionGroup(rep="ijk")
    #Q2 = QuaternionGroup(rep="permutations")
    #print(Q.is_isomorphic(Q2))
    
    '''
    G = RootsOfUnitGroup(5)
    #print_roots(G)
    #print(G.generators())
    print(G.is_abelian())
ithu    print(G.generators())
    
    
    
    Q = QuaternionGroup(rep="ijk")
    Q2 = QuaternionGroup(rep="permutations")
    #print(Q.is_isomorphic(Q2))
    
    print(Q.generators())
    
    
    Dm = DihedralGroup(4, rep='matrix')
    Drs = DihedralGroup(4, rep='RS')
    Dp = DihedralGroup(4, rep='permutations')
    print(Drs.generators())
    #print(Dp.is_isomorphic(Drs))
    #print(Dp.is_isomorphic(Dm))
    #print(Dm.is_isomorphic(Drs))
    '''
    
    
    
    '''
    G= SymmetricGroup(3)
    G2= DihedralGroup(3)
    print(G.is_isomorphic(G2))
    
    
    r1 = Drs('R1')
    s0 = Drs('S1')
    other = Drs.generate([r1, s0])
    print(other)
    
    G3= AlternatingGroup(3)
    #print(G3.is_normalSubgroup(G))    
    #print(H.is_normalSubgroup(G))        
    C = CyclicGroup(6)
    #print(C.elements_order())
    '''
    
    
    '''
    print(C)
    print(C.all_subgroups())
    print(C.Cayley_table())
    
    print(C.gens_cyclic_group())
    '''
    
    S = SymmetricGroup(3)
    #print(S.elements_order())
    #A = AlternatingGroup(3)
    #print(S.elements_order())
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
    Z6_ = Group(S, b_op6)   
    
    
    #Grupo Z_3
    S=Set(range(3))
    b_op3=Function(S*S, S,lambda x: (x[0]+x[1])%3)
    Z3 = Group(S, b_op3)  
    
    #Grupo Z_2
    S=Set(range(2))
    b_op2=Function(S*S, S,lambda x: (x[0]+x[1])%2)
    Z2 = Group(S, b_op2)  
    
    #Grupo Z_1
    S=Set(range(1))
    b_op1=Function(S*S, S,lambda x: (x[0]+x[1])%1)
    Z1 = Group(S, b_op1)  
    
    
    Z6 = Z2.direct_product(Z3)
    
    Z2x1 = Z2.direct_product(Z1)
    Z1x3 = Z1.direct_product(Z3)
    
    '''
    print(Z2x1.is_normalSubgroup(Z6))
    print(Z1x3.is_normalSubgroup(Z6))

    print(Z6.is_direct_product(Z2x1,Z1x3))
    '''
    
    
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