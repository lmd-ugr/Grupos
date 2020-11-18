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



class Complex:
    def __init__(self, real, imag=0.0):
        self.real = round(real, 3)
        self.imag = round(imag, 3)


    def __add__(self, c2):
        return Complex(self.real + c2.real,
                       self.imag + c2.imag)

    def __sub__(self, c2):
        return Complex(self.real - c2.real,
                       self.imag - c2.imag)

    def __mul__(self, c2):
        return Complex(self.real*c2.real - self.imag*c2.imag,
                       self.imag*c2.real + self.real*c2.imag)

    def __div__(self, c2):
        r = float(c2.real**2 + c2.imag**2)
        return Complex((self.real*c2.real+self.imag*c2.imag)/r, (self.imag*c2.real-self.real*c2.imag)/r)

    def __abs__(self):
        return math.sqrt(self.real**2 + self.imag**2)

    def __neg__(self):
        return Complex(-self.real, -self.imag)

    def __eq__(self, c2):
        return abs(self.real-c2.real<0.01) and abs(self.imag-c2.imag<0.01)

    def __ne__(self, c2):
        return not self.__eq__(c2)

    def product(self,c2):
        return Complex(self.real*c2.real - self.imag*c2.imag,
                       self.imag*c2.real + self.real*c2.imag)

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








def plot(roots, mode="exp"):
    """
    Plot the n roots of unit of group GroupOfUnitGroup

    Args:
        roots: instance of GroupOfUnitGroup
        mode:
            'exp' if user wants representation with e^{pi i k/n} or
            'binom' if user want binomial representation (a+bj).
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
    plt.savefig('roots.png')

    plt.show()
