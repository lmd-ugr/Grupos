# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 07:46:08 2020

@author: Alberto
"""

import numbers
import numpy as np
import math
#from Permutation import permutation
#from Group import Group
#from Set import Set
#from Function import Function

class Quaternion:
    "Quaternions Numbers:  +-{1, i, j, k}"
    
    def __init__(self, real=0, imag_i=0, imag_j=0, imag_k=0, letter=None):
        
        cond = (real==0 and imag_i==0 and imag_j==0 and imag_k==0)
        
        if(letter!=None and cond):
            if letter=="1":
                real=1
            elif letter=="-1":
                real=-1
            elif letter=="i":
                imag_i=1
            elif letter=="-i":  
                imag_i=-1
            elif letter=="j":
                imag_j=1
            elif letter=="-j":
                imag_j=-1
            elif letter=="k":  
                imag_k=1
            elif letter=="-k":
                imag_k=-1
            else:
                raise ValueError("Incorrect letter, not a valid Quaternion")
        
        elif letter!=None and not cond:
            raise ValueError("Not a valid Quaternion")
            
        self.assign_real(real)
        self.assign_imags(imag_i, imag_j, imag_k)
    
    
    
    def __repr__(self):
        
        coords = {"real": self.real, "i": self.imag_i, "j": self.imag_j, "k":self.imag_k}
        
        values = [s for s in coords if coords[s] != 0]
        
        dev = " "
        cont = 0
        for x in values:
            if x=="real":
                dev = dev + str(coords[x]) + str(x)[:-4]
            else:
                sign=""
                if 0< cont <len(values) and coords[x]>0:
                    sign = "+"
                    
                if abs(coords[x])==1:
                    dev = dev + sign+ str(coords[x]).replace(str(1),"") + str(x)
                else:
                    dev = dev + sign+ str(coords[x]) +  str(x)
                    
            cont=cont+1
        return dev
    
    
    __str__ = __repr__
    
    #def __str__ (self):
        #return "({0.real:.2f}{0.imag_i:+.2f}i{0.imag_j:+.2f}j{0.imag_k:+.2f}k)".format(self)

    def assign_real(self, real):
        self.real = real

    def assign_imags(self, imag_i, imag_j, imag_k):
        self.imag_i = imag_i
        self.imag_j = imag_j
        self.imag_k = imag_k

    def get_real(self):
        return self.real

    def get_i(self):
        return self.imag_i
    
    def get_j(self):
        return self.imag_j
    
    def get_k(self):
        return self.imag_k
    
    def __call__(self):
        return self

    def __add__(self, other):
        if not isinstance(other,Quaternion):
            raise TypeError(other, "is not a Quaternion")
            
        return Quaternion(self.real + other.real, self.imag_i + other.imag_i,
                          self.imag_j + other.imag_j, self.imag_k + other.imag_k)

    def __iadd__(self, other):
        self.assign_real(self.real + other.real)
        self.assign_imags(self.imag_i + other.imag_i, self.imag_j + other.imag_j,
                          self.imag_k + other.imag_k)
        return self
    
    def __sub__(self, other):
        if not isinstance(other,Quaternion):
            raise TypeError(other, "is not a Quaternion")
            
        return Quaternion(self.real - other.real, self.imag_i - other.imag_i,
                          self.imag_j - other.imag_j, self.imag_k - other.imag_k)
    
    
    def __eq__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            conditions = (self.real==other and self.imag_i==0 and self.imag_j==0 and self.imag_k==0)
            return True if conditions else False
        
        return (self.real==other.real and self.imag_i==other.imag_i and
                self.imag_j==other.imag_j and self.imag_k==other.imag_k)


    def __mul__(self, other):
        #if not isinstance(other,Quaternion) and not isinstance(other, numbers.Number):
         #   raise TypeError(other, "is not a Quaternion")
            
        if isinstance(other, Quaternion):
            return Quaternion(
                
                self.real*other.real - self.imag_i*other.imag_i - \
                self.imag_j*other.imag_j - self.imag_k*other.imag_k,
                
                self.imag_i*other.real + self.real*other.imag_i  + \
                self.imag_j*other.imag_k - self.imag_k*other.imag_j , 
                
                self.real*other.imag_j - self.imag_i*other.imag_k + \
                self.imag_j*other.real + self.imag_k*other.imag_i, 
                
                self.real*other.imag_k + self.imag_i*other.imag_j - \
                self.imag_j*other.imag_i + self.imag_k*other.real
            )
        elif isinstance(other, numbers.Number):
            return Quaternion(self.real*other, self.imag_i*other,
                              self.imag_j*other, self.imag_k*other)
        else:
            raise TypeError(other, "is not a Quaternion or number")
    
    def __hash__(self):
        return hash(self.real) + 2*hash(self.imag_i)+ \
            3*hash(self.imag_j)+ 4*hash(self.imag_k)
            
            
    def __imul__(self, other):
        if not isinstance(other,Quaternion):
            raise TypeError(other, "is not a Quaternion")
            
        nreal = self.real*other.real - self.imag_i*other.imag_i - \
            self.imag_j*other.imag_j - self.imag_k*other.imag_k
            
        nimag_i = self.imag_i*other.real + self.real*other.imag_i  + \
            self.imag_j*other.imag_k - self.imag_k*other.imag_j
            
        nimag_j = self.real*other.imag_j - self.imag_i*other.imag_k + \
            self.imag_j*other.real + self.imag_k*other.imag_i
        
        nimag_k = self.real*other.imag_k + self.imag_i*other.imag_j - \
            self.imag_j*other.imag_i + self.imag_k*other.real
            
        self.assign_real(nreal)
        self.assign_imags(nimag_i, nimag_j, nimag_k)
        
        return self
    
    
    def __rmul__(self, number):
        if not isinstance(number,numbers.Number):
            raise TypeError(number, "is not a number")
            
        return Quaternion(self.real * number, self.imag_i * number,
                          self.imag_j * number, self.imag_k * number)


    def __neg__(self):
         return Quaternion(-self.real, -self.imag_i, -self.imag_j, -self.imag_k)


    def __div__(self, other):
        return Quaternion(self.real/other, self.imag_i/other, self.imag_j/other, self.imag_k/other)

    __truediv__ = __div__
    
    
    def conjugate(self):
        return Quaternion(self.real, - self.imag_i, - self.imag_j, - self.imag_k)
    
    
    def matrix(self):
        r, i, j, k = self.real, self.imag_i, self.imag_j, self.imag_k 
        
        return np.array([
            [r * r + i * i - j * j - k * k,
                2 * (i * j + r * k),
                2 * (i * k - r * j)],
            [2 * (i * j - r * k),
                r * r - i * i + j * j - k * k,
                2 * (j * k + r * i)],
            [2 * (i * k + r * j),
                2 * (j * k - r * i),
                r * r - i * i - j * j + k * k]
        ])
    
    def norm(self):
        coords=(self.real, self.imag_i, self.imag_j, self.imag_k)
        coords_2= [x*x for x in coords]
        return math.sqrt(sum(coords_2))

    
    def inverse(self):
        return self.conjugate()/(self.norm()**2)


    def trace(self):
        return self+self.conjugate()
    
    '''
    @classmethod        
    def Group(cls, rep="ijk"):
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
            Gr.group_gens=[i,j]
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
    '''

if __name__ == '__main__':
    
    '''
    i2 = Quaternion(0,1,0,0)
    i = Quaternion(letter="i")

    j2 = Quaternion(0,0,1,0)
    j = Quaternion(letter="j")

    k2 = Quaternion(0,0,0,1)
    k = Quaternion(letter="k")
    
    print(i, "*", k, "=", i*k)
    
    
    
    Q = Quaternion.Group(rep="ijk")
    
    #print(Q.all_normalSubgroups())
    #print(Q.gens_group())
    a = Q.group_gens
    b= Q.gens_group()
    
    
    print("##################")
    print(b[0] , a)
    print(type(a[0]), type(b[0]))
    print(type(a), type(b))
    print("{} vs {}".format(a,b))
    #print(Q.elements_order())
    #print(Q.is_abelian())
    #print(Q.Cayley_table())
    '''
    
    
    '''
    q = Quaternion(68,12,3,-9)
    r = Quaternion(-8,-2,2,32)
    print((q*r).conjugate() == r.conjugate()*q.conjugate())
    print((q*r).trace() == (r*q).trace())
    print(round((q*r).norm(), 10) ==  round(q.norm()*r.norm(), 10))
    
    print(( (r*q*r.inverse()).trace()) , " ", q.trace())


    print(i*i == j*j == k*k == i*j*k == -1)
    print(i*j==k, " ", j*k==i, " ", k*i==j)
    print(j*i==-k, " ", k*j==-i, " ", i*k==-j)
  
    
    #print(i.matrix())
    #print(j.matrix())
    #print(k.matrix())
    
    #print(n1.conjugate())
    '''