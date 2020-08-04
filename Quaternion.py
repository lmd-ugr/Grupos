# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 07:46:08 2020

@author: Alberto
"""

import numbers
import numpy as np
import math
from Permutation import permutation
from Group import Group
from Set import Set
from Function import Function

class Quaternion:
    "Quaternions Numbers: i, j, k"
    
    def __init__(self, real=0, imag_i=0, imag_j=0, imag_k=0):
        self.assign_real(real)
        self.assign_imags(imag_i, imag_j, imag_k)
    
    def __repr__(self):
        return "("+ str(self.real) + "," + str(self.imag_i) + "i"
        + "," + str(self.imag_j) + "j" + "," + str(self.imag_k) + "k)"
    
    def __str__ (self):
        if(self.imag_i==0 and self.imag_j==0 and self.imag_k==0):
            return str(self.real)
        
        return "({0.real:.2f}{0.imag_i:+.2f}i{0.imag_j:+.2f}j{0.imag_k:+.2f}k)".format(self)

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
    
    def __eq__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            conditions = (self.real==other and self.imag_i==0 and self.imag_j==0 and self.imag_k==0)
            return True if conditions else False
        
        return (self.real==other.real and self.imag_i==other.imag_i and
                self.imag_j==other.imag_j and self.imag_k==other.imag_k)


    def __mul__(self, other):
        if not isinstance(other,Quaternion):
            raise TypeError(other, "is not a Quaternion")
            
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
        #return math.sqrt(self*self.conjugate())
        coords=(self.real, self.imag_i, self.imag_j, self.imag_k)
        coords_2= [x*x for x in coords]
        return math.sqrt(sum(coords_2))

    
    def inverse(self):
        return self.conjugate()/(self.norm()**2)

    def trace(self):
        return self+self.conjugate()
    
    
    @classmethod
    def Group2(cls):
        print("gola")
        
    def Group(cls, rep="ijk"):
        """
        The quaternion group Q2; its elements are 1,i,j,k and their oposite
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
            one = Quaternion(1,0,0,0)
            i   = Quaternion(0,1,0,0)
            j   = Quaternion(0,0,1,0)
            k   = Quaternion(0,0,0,1)
            
            #q2=[ "1", "-1", "i", "-i", "j", "-j", "k", "-k"]
            q2=[one, -one, i, -i, j, -j, j, -k]
            print(q2)
            table=[[ "1", "-1", "i", "-i", "j", "-j", "k", "-k"],
                   [ "-1", "1", "-i", "i", "-j", "j", "-k", "k"],
                   [ "i", "-i", "-1", "1", "k", "-k", "-j", "j"],
                   [ "-i", "i", "1", "-1", "-k", "k", "j", "-j"],
                   [ "j", "-j", "-k", "k", "-1", "1", "i", "-i"],
                   [ "-j", "j", "k", "-k", "1", "-1", "-i", "i"],
                   [ "k", "-k", "j", "-j", "-i", "i", "-1", "1"],
                   [ "-k", "k", "-j", "j", "i", "-i", "1", "-1"]]
            
            '''
            def product(a,b):
                i=q2.index(a)
                j=q2.index(b)
                return table[i][j]
            '''
            def product(a,b):
                return a*b
            '''
            G=Set(q2)
            Gr=Group(G,Function(G.cartesian(G),G, lambda x: product(x[0],x[1])))
            Gr.group_gens=[Gr("i"),Gr("j")]
            return Gr
            '''
        
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


if __name__ == '__main__':
    
    
    i = Quaternion(0,1,0,0)
    j = Quaternion(0,0,1,0)
    k = Quaternion(0,0,0,1)
    
    q = Quaternion(2,3,1,-2)
    r = Quaternion(-4,5,-8,3)


    q = Quaternion.Group2()
    #print(q.Cayley_table())
    
    '''
    print((q*r).conjugate() == r.conjugate()*q.conjugate())
    print((q*r).trace() == (r*q).trace())
    print(round((q*r).norm(), 10) ==  round(q.norm()*r.norm(), 10))
    print((r*q*r.inverse()).trace() == q.trace())


    print(i*i == j*j == k*k == i*j*k == -1)
    print(i*j==k, " ", j*k==i, " ", k*i==j)
    print(j*i==-k, " ", k*j==-i, " ", i*k==-j)
    '''
    
    #print(i.matrix())
    #print(j.matrix())
    #print(k.matrix())
    
    #print(n1.conjugate())