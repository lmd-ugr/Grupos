#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 21:55:57 2020

@author: albduranlopez
"""

import math
import pylab as plt
import seaborn as sns
import sympy as sym



class Complex(object):
    def __init__(self, real, imag=0.0):
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

    def __pow__(self, power):
        raise NotImplementedError('self**power is not implemented for Complex')


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
        return "({},{}j)".format(self.real, self.imag)
    
    
    '''
    def __str__(self):
        return '(%g, %g)' % (self.real, self.imag)
    '''






def plot(roots, mode="exp"):
    """ 
    Plot the n roots of unit of group GroupOfUnitGroup
    
    Args: 
        roots: instance of GroupOfUnitGroup
        mode: 
            'exp' if user wants representation with e^{pi i k/n} or
            'real' if user want real representation (a+bj).
    """
            
    lroots = list(roots.Set)
    lreal = [ r.real for r in lroots]
    limag = [ r.imag for r in lroots]
    n,t = len(roots.Set), 0

    colors = sns.color_palette("hls", n)
    plt.axis('square')
    plt.scatter(lreal, limag, c=colors) #plot the roots
    
    k = 0
    for root,c in zip(lroots,colors):
        t=t+0.25/n
        k=k+1

        if mode == "exp" :        
            z = sym.exp(2*sym.pi * sym.I*k/n)
            plt.text(sym.re(z)+.04, sym.im(z)-.025, z)
        elif mode == "binom":
            plt.text(root.real+.04, root.imag-.025, str(root))
        else:
            raise TypeError("mode is not valid.")

    plt.xlim(-1.25,1.25+t)
    plt.ylim(-1.25,1.25)
    plt.xlabel('Re(x)')
    plt.ylabel('Im(x)')
    plt.savefig('test/photo2.png')

    plt.show()

