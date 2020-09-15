#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 21:55:57 2020

@author: albduranlopez
"""
import math

#from Set import Set
#from Function import Function
#from Group import Group

#from sympy import pi
import pylab as plt
import seaborn as sns

class Complex(object):
    def __init__(self, real, imag=0.0):
        #self.real = real
        #self.imag = imag
        self.real = round(real, 3)
        self.imag = round(imag, 3)


    def __add__(self, other):
        return Complex(self.real + other.real,
                       self.imag + other.imag)

    def __sub__(self, other):
        return Complex(self.real - other.real,
                       self.imag - other.imag)

    def __mul__(self, other):
        return Complex(self.real*other.real - self.imag*other.imag,
                       self.imag*other.real + self.real*other.imag)

    def __div__(self, other):
        sr, si, orr, oi = self.real, self.imag, other.real, other.imag # short forms
        r = float(orr**2 + oi**2)
        return Complex((sr*orr+si*oi)/r, (si*orr-sr*oi)/r)

    def __abs__(self):
        return math.sqrt(self.real**2 + self.imag**2)

    def __neg__(self):   # defines -c (c is Complex)
        return Complex(-self.real, -self.imag)

    def __eq__(self, other):
        #return self.real == other.real and self.imag == other.imag
        return abs(self.real-other.real<0.01) and abs(self.imag-other.imag<0.01)
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def product(self,other):
        return Complex(self.real*other.real - self.imag*other.imag,
                       self.imag*other.real + self.real*other.imag)
    '''
    def __str__(self):
        return '(%g, %g)' % (self.real, self.imag)
    '''
    def __pow__(self, power):
        raise NotImplementedError('self**power is not yet impl. for Complex')


    def __hash__(self):
        if not self.imag:
            return hash(self.real)
        return hash((self.real, self.imag))

    def __repr__(self):

        
        if not self.imag:
            return "({})".format(self.real)
        else:
            return "({},{}j)".format(self.real, self.imag)

    def __str__(self):
        if not self.imag:
            return repr(self.real)
        else:
            return "({},{}j)".format(self.real, self.imag)
    
    '''
    def round3(self):
        return Complex(round(self.real, 3), round(self.imag,3))
    
    def __repr__(self):
        #return str(self.real)+","+ str(self.imag)
        return str(self.round3())
    '''
    
    
    '''
    @classmethod
    def RootsOfUnitGroup(cls, n):
        
        
        #G = Set(complex( math.cos(2*pi*k/n), math.sin(2*pi*k/n)) for k in range(n))
        #G = Set( (np.round(math.cos(2*pi*k/n),3), np.round(math.sin(2*pi*k/n),3)) for k in range(n))
    
        #def mul(x):
            #return complex(x[0]*y[0]-x[1]*y[1], x[1]*y[0]+x[0]*y[1])
            #return (np.round(x[0][0]*x[1][0]-x[0][1]*x[1][1],3), np.round(x[0][1]*x[1][0]+x[0][0]*x[1][1],3))
        
        #bin_op=Function(G.cartesian(G),G,lambda x: mul(x), check_well_defined=False)
        #Gr=Group(G,bin_op) 
        #return Gr
        
        
        G = Set(Complex( math.cos(2*pi*k/n), math.sin(2*pi*k/n)) for k in range(n))
        #print(G)
        bin_op=Function(G.cartesian(G),G,lambda x: x[0].product(x[1]) ,  check_well_defined=False)
        Gr=Group(G,bin_op) 
    
        return Gr
    '''


def print_roots(roots):
    
    lroots = roots.Set
        
    n = len(roots.Set)
    colors = sns.color_palette("hls", n)
       
    plt.figure(figsize=(9,9))
    t=0
    arrows, points, text = False, True, True
    for root,c in zip(lroots,colors):
        if arrows:
            plt.arrow(0,0,root.real,root.imag,ec=c,lw=3)
        if points:
            plt.scatter(root.real, root.imag, c=c)
        if text:
            t=t+0.25/n
            #plt.annotate(str(root), (root.real, root.imag))
            plt.text(root.real+.04, root.imag-.025, str(root))
        
        
    plt.xlim(-1.25,1.25+t)
    plt.ylim(-1.25,1.25)
    plt.xlabel('Re(x)')
    plt.ylabel('Im(x)')
    plt.savefig('test/photo.png')
        
    plt.show()
    
    
    

if __name__ == '__main__':
    '''
    G = Complex.RootsOfUnitGroup(5)
    print_roots(G)
    
    print(G.is_abelian())
    '''
    
