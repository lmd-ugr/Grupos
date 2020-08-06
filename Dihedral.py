#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:18:45 2020

@author: albduranlopez







PROBLEMAS:
-Los elementos no pueden ser listas ni matrices porque no son hashables.
Hay que aplicarle la operación binaria y no compila.

La única representación válida para los elementos que he podido compilar
sin problemas ha sido tuple().

Por otro lado, estaría bien trabajar a nivel de matrices y realizar operaciones
en el grupo Diédrico. Convedría hacer esa modificación?:
    -Trabajar con matrices y a la hora de crear el grupo convertir las matrices a tuplas
    
-Identificar los elementos con RO, R1,...RN, S0, S1,... SN
Para ello podemos usar un diccionario y así quitar rep="RS" y juntar "matrix" y "RS" en una misma
"""



from Set import Set
from Group import Group, SymmetricGroup
from Function import Function
from Permutation import permutation

import math
import numpy as np

class Dihedral():
    
    def __init__(self, n):
    
        if not isinstance(n, int):
            raise TypeError("n is not an integer")
            
        self.n = n
        self.rot =  []
        self.refl = []
        
        def portion(i):
            return (2*math.pi)*i/n
        
        for i in range(n):
            c, s = math.cos(portion(i)), math.sin(portion(i))
            
            #self.rot.append(np.round(  ( (c,-s),(s, c) )  , 4))
            #self.refl.append(np.round( ( (c, s),(s,-c) )  , 4))
            
            #self.rot.append(np.round(np.matrix(  [ [c,-s],[s, c] ]  ), 4))
            #self.refl.append(np.round(np.matrix( [ [c, s],[s,-c] ]  ), 4))
            
            
            self.rot.append(  tuple( np.round([c, -s, s, c],10) )  )  
            self.refl.append( tuple( np.round([c, s, s, -c],10) )  )  
        
        self.all = dict(zip(["R"+str(i) for i in range(0,n)] + ["S"+str(i) for i in range(0,n)], 
                             self.rot+self.refl))
        
        
    def __repr__(self):
        dev = "Rotaciones:"
        for i in self.rot:
            dev = dev+ "\n"+ str(i)+ "\n"
            
        dev = dev+"\nReflexiones"
        for i in self.refl:
            dev = dev+  "\n"+ str(i)+ "\n"
            
        return dev
        #return "Rotaciones: " + str(self.rot) + "\nTraslaciones: "+str(self.refl)
        

    __str__ = __repr__
    
    def __eq__(self, other):
        if not isinstance(other, Dihedral):
            raise TypeError(other, "is not a DIhedral Group")
            
        return self.rot==other.rot and self.refl==other.refl
    
    def __ne__(self, other):
        return not self==other

    def get_rot(self,i):
        return self.rot[i]
    
    def get_refl(self,i):
        return self.refl[i]
    
    def get_all(self):
        return self.all
    
    
    
    def product_RS(self, a, b):
        
        a = self.all[str(a)]
        b = self.all[str(b)]

        if a  in self.rot :
            i = self.rot.index(a)
        else:
            i = self.refl.index(a)

        
        if b in self.rot :
            j = self.rot.index(b)
        else:
            j = self.refl.index(b)
        
        
        keys = list(self.all.keys())
        values = list(self.all.values())

        
        if a in self.rot and b in self.rot:
            this = self.rot[(i+j)%self.n]
        
        if a in self.rot and b in self.refl:
            this = self.refl[(i+j)%self.n]
        
            
        if a in self.refl and b in self.rot:
            this = self.refl[(i-j)%self.n]
           
            
        if a in self.refl and b in self.refl:
            this = self.rot[(i-j)%self.n]
        
        
        return  keys[values.index(this)] 
        
    
    
    def product_Matrix(self, a, b):
        

        if a  in self.rot :
            i = self.rot.index(a)
        else:
            i = self.refl.index(a)

        
        if b in self.rot :
            j = self.rot.index(b)
        else:
            j = self.refl.index(b)
        
        
        if a in self.rot and b in self.rot:
            return self.rot[(i+j)%self.n]
        
        if a in self.rot and b in self.refl:
            return self.refl[(i+j)%self.n]
            
        if a in self.refl and b in self.rot:
            return self.refl[(i-j)%self.n]
            
        if a in self.refl and b in self.refl:
            return self.rot[(i-j)%self.n]
    
    
    
    @classmethod
    def Group(cls, n, rep="RS"):
        
    
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
            
            #Gr.group_gens=[Gr('R1'),Gr('S0')]
            return Gr
        
        
        if rep=="RS":
            D = Dihedral(n)        
            G = Set(D.all.keys())
            Gr=Group(G, Function(G.cartesian(G), G, lambda x: D.product_RS(x[0],x[1])))
            
            Gr.group_gens=['R1','S0']
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
    
    
    
    
    
    




if __name__=="__main__":
    

    G = Dihedral.Group(2,"RS")
    D3 = Dihedral.Group(2,"matrix")    
    D3_1 = Dihedral.Group(2, "permutations")
    

    print("R_i y S_i: ")
    print(G.Cayley_table())
    print("\nForma matricial (tuplas): ")
    print(D3.Cayley_table())
    print("\nCon permutaciones:")
    print(D3_1.Cayley_table())

    
    
    #print(G.all_normalSubgroups())
    #print(D3.all_subgroups())

    
    
    D5 = Dihedral.Group(5,"matrix")
    D3 = Dihedral.Group(3,"matrix")
    S3 = SymmetricGroup(3)
    
    #print(D3.elements_order(), "\n\n", S3.elements_order())
    #print(D3.Cayley_table())
    #print(D5.all_subgroups())
    
    
    
