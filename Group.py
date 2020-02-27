# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:18:18 2020

@author: Alberto
"""


import itertools
from itertools import combinations, chain 

#import functools
#import operator
#import math


#import Set

#from Function import Function
import Function
#from fractions import gcd
#from copy import deepcopy
#from sympy.ntheory import factorint,totient
#from sympy.utilities.iterables import flatten




    

    

















































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


A = set({1,2,3})

#print(3 in A)



