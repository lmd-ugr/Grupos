# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:37:32 2020

@author: Alberto
"""

import functools
import operator
import math

class permutation:
    """
    This is the class of permutations of the set {1..n}

    Attributes:
        tuple: a tuple of n integers that stores the images of 1..n under the permutation
        length: n with the above notation
    """
    def __init__(self, *els):
        """
        Defines a permutation

        Args:
            *els: can be a list of integers or a list of tuples integers or a sequence of integers or a sequence of tuples of integers

        Returns:
            A permutation.
            If a list or sequence of integers is given, then the output is the permutation that sends 1 to the first element in the list or sequence, 2 to the second, and so on.
            If a list or sequence of tuples is given, then they are considered as cycles, and the output is the product of these cycles.

        Example:
            >>> p=permutation([2,1,4,3])
            >>> q=permutation((1,2),(3,4))
            >>> p==q
            True
            >>> p=permutation(2,1,4,3)
            >>> p==q
            True

        """


        def cycle2perm(c):
            """
            Returns a permutation corresponding to the cycle
            given by the tuple c.
            
            Args:
                c : tuple
            """
            m=len(c)
            if len==0: #this will not happen
                return tuple([1])
            n=max(c)
            p=[i+1 for i in range(n)]
            for i in range(1,n+1):
                if (i in c):
                    p[i-1]=c[(c.index(i)+1)%m]
            return permutation(p)


        if len(els)==1:
            t=els[0]
        else:
            t=list(els)
        if not(isinstance(t,list)) and not(isinstance(t,tuple)):
            raise TypeError("expecting a list or sequence of integers or tuples of integers as argument")
        n=len(t)
        if n==0:
            raise TypeError("to avoid ambiguity empty permutations are not allowed")
        if all(isinstance(x,int) for x in t) and isinstance(t,list):
            if set(t)!=set(range(1,n+1)):
                raise TypeError("the input is not a permutation of "+str(range(1,n+1)))
            self.tuple=tuple(t)
            self.length=n
        elif all(isinstance(x,int) for x in t) and isinstance(t,tuple):
            p=cycle2perm(t)
            self.tuple=p.tuple
            self.length=p.length
        elif all(isinstance(x,tuple) for x in t):
            cs=[cycle2perm(c) for c in t]
            p=functools.reduce(operator.mul,cs)
            self.tuple=p.tuple
            self.length=p.length
        else:
            raise TypeError("expecting a list or sequence of integers or tuples of integers as argument")

    def __call__(self,n):
        if 0 < n <= self.degree():
            return self.tuple[n-1]
        else:
            return n


    def __hash__(self):
        """ 
        Returns the hash of self.tuple
        """
        return hash(self.tuple)

    '''
    def __str__(self):

        s=str(list(self.tuple))+" = "

        s2=str(" ")
        for c in self.disjoint_cycles():
            if len(c)>1:
                s2=s2+str(c)
        if s2==str(" "):
            str(list(self.tuple))
            #return s +"( )"
        return s+s2
    '''

    def __repr__(self):
        s2=str("")
        for c in self.disjoint_cycles():
            if len(c)>1:
                s2=s2+str(c)
        if s2==str(""):
            return "()"
        return s2

    __str__ = __repr__

    def __eq__(self, other):
        """
        Tests if the permutations are identical (with the same length).
        """

        if not isinstance(other, permutation):
            raise TypeError("other is not a permutation")
        #if self.order()==other.order():
         #   return True
        return (self.tuple == other.tuple)


    def __ne__(self, other):
        """
        Tests if the permutations are identical  or not (with the same length).
        """

        return not self == other

    def degree(self):
        """ 
        Returns the len of the permutation.
        """
        return len(self.tuple)

    def __mul__(self,other):
        """
        Composition of permutations.

        Example:
            >>> p=permutation((1,3))
            >>> q=permutation([2,1,3])
            >>> p*q
             (1, 2, 3)
        """

        if not(isinstance(other,permutation)):
            raise TypeError("other must also be a permutation")

        l = []
        for i in range(max(self.degree(), other.degree())):
            l.append(self(other(i+1)))

        return permutation(l)

    def inverse(self):
        """
        Inverse of a permutation (as a function).

        Example:
            >>> q=permutation([2,3,1])
            >>> q
             (1, 2, 3)
            >>> q**-1
             (1, 3, 2)
        """
        l=list(self.tuple)
        return permutation([l.index(i)+1 for i in range(1,len(l)+1)])

    def __pow__(self, n):
        """
        Returns the composition of a permutation n times.
        
        Args:
            n : integer
        Example:
            >>> q=permutation([2,3,1])
            >>> q**3
            ( )
            >>> q*q==q**2
            True
        """
        if not (isinstance(n, int)):
            raise TypeError("n must be an int or a long")

        k=self.length
        if n == 0:
            return permutation(list(range(1,k+1)))
        elif n < 0:
            return self.inverse() ** -n
        elif n % 2 == 1:
            return self * (self ** (n - 1))
        else:
            return (self * self) ** (n // 2)

    __xor__=__pow__

    def disjoint_cycles(self):
        """
        Returns a list of disjoint cycles (as tuples) whose product is the given permutation (argument)

        Example:
            >>> p=permutation([2,1,4,3])
            >>> p.disjoint_cycles()
            [(1, 2), (3, 4)]
        """
        l=[]
        p=list(self.tuple)
        els=set(p)
        while len(els)>0:
            e=next(iter(els))
            c=[e]
            while not(p[e-1] in c):
                e=p[e-1]
                c.append(e)
            l.append(tuple(c))
            els=[a for a in els if not(a in c)]
        return l #[c for c in l if len(c.tuple)>1]

    def inversions(self):
        """
        List of inversions of the given permutation p, that is, the set of pairs (i,j) in {1..n}^2 with i<j such that p(i)>p(j)

        Example:
            >>> q=permutation([2,3,1])
            >>> q.inversions()
            [(1, 3), (2, 3)]
        """
        p=list(self.tuple)

        return [tuple([i,j]) for i in range(1,len(p)+1) for j in range(i+1,len(p)+1) if p[i-1]>p[j-1]]

    def sign(self):
        """
        The sign of the permuation, that is, (-1)^i, with i the number of inversions of the permutation
        Example:
            >>> q=permutation([2,3,1])
            >>> q.inversions()
            [(1, 3), (2, 3)]
            >>> q.sign()
            1
        """
        return (-1)**len(self.inversions())

    def order(self):
        """
        The order of the permutation, that is, minimum positive integer such that p**n==identity, with p the argument
        Example:
            >>> p=permutation([4,3,1,2])
            >>> p.order()
            4
        """
        l=[len(c) for c in self.disjoint_cycles()]
        if len(l)==1:
            return l[0]

        return (functools.reduce(operator.mul,l))//(functools.reduce(math.gcd,l))

    def extend(self,n):
        """
        Extends the permutation to a permutation of the set {1..n} leaving the elements above its length untouched
        """
        if not(isinstance(n,int)) or (n<self.length):
            raise ValueError("Either the argument is not an integer or it is less than the length of the permuataion")
        tmp=list(range(1,n+1))
        for i in range(self.length):
            tmp[i]=(self.tuple)[i]
        return permutation(tmp)


    def even_permutation(self):
        """ 
        Checks if self is even permutation.
        """
        dev=0
        for c in self.disjoint_cycles():
            dev= dev + len(c)-1
        return True if dev%2==0 else False


    def odd_permutation(self):
        """ 
        Checks if self is odd permutation.
        """
        return not self.even_permutation()

