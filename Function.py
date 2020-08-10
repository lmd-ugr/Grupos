# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 18:25:42 2020

@author: Alberto


Definition of a function
"""

from Set import Set


class Function:
    
    
    def __init__(self, domain, codomain, function, check_well_defined=True):
        """
        Initialize the function and check that it is well-formed.
        This method can be overwritten by subclasses of Function, so that for
        example GroupHomomorphisms can be between Groups, rather than Sets.
        """
        if not isinstance(domain, Set):
            raise TypeError("Domain must be a Set")
            
        if not isinstance(codomain, Set):
            raise TypeError("Codomain must be a Set")
            
        if check_well_defined:
            #if not all(function(elem) in codomain for elem in domain):
                #raise TypeError("Function returns some value outside of codomain")
            
            for elem in domain:
                #print(elem, "in", codomain)
                if(function(elem) not in codomain):
                    print(function(elem), " not in " , codomain)
                    raise TypeError("Function returns some value outside of codomain")

                    
        self.domain = domain
        self.codomain = codomain
        self.function = function

    
    def __call__(self, elem):
        if elem not in self.domain:
            print(elem, " not in ", self.domain)
            #raise TypeError("Function must be called on elements of the domain")
        return self.function(elem)
    
    
    def __hash__(self):
        """Returns the hash of self"""

        # Need to be a little careful, since self.domain and self.codomain are
        # often the same, and we don't want to cancel out their hashes by xoring
        # them against each other.
        #
        # Also, functions we consider equal, like lambda x: x + 1, and
        # def jim(x): return x + 1, have different hashes, so we can't include
        # the hash of self.function.
        #
        # Finally, we should make sure that if you switch the domain and
        # codomain, the hash will (usually) change, so you can't just add or
        # multiply the hashes together.

        return hash(self.domain) + 2 * hash(self.codomain)


    def __eq__(self, other):
        if not isinstance(other, Function):
            return False

        return id(self) == id(other) or ( \
               self.domain == other.domain and \
               self.codomain == other.codomain and \
               all(self(elem) == other(elem) for elem in self.domain) )

    def __ne__(self, other):
        return not self == other

    def _image(self):
        """The literal image of the function"""
        return Set(self(elem) for elem in self.domain)

    def image(self):
        """
        The API image of the function; can change depending on the subclass.
        For example, GroupHomomorphisms return the image as a Group, not a Set.
        """
        return self._image()

    '''
    #Print
    def __str__(self):
        """Pretty outputing of functions"""

        # Figure out formatting
        maxlen = max(len(str(x)) for x in self.domain) if self.domain else 0
        formatstr1 = "{0:<%d} -> {1}\n" % maxlen
        formatstr2 = "{0:<%d}{1}\n" % (maxlen + 4)
        nothit = self.codomain - self._image()

        return("".join(formatstr1.format(x, self(x)) for x in self.domain) + \
               "".join(formatstr2.format("", y) for y in nothit))
    '''
    
    def __str__(self):
        l1 = []
        for i in self.domain:
            l1.append( "f(" + str(i)+ ")=" + str(self(i)) + "\n")
        
        return ' '.join(l1)
    
    
    def is_surjective(self):
        # Need to make self.domain into a Set, since it might not be in
        # subclasses of Function
        return self._image() == Set(self.codomain)

    def is_injective(self):
        return len(self._image()) == len(self.domain)

    def is_bijective(self):
        return self.is_surjective() and self.is_injective()

    
    def compose(self, other):
        """Returns x -> self(other(x))"""
        if not self.domain == other.codomain:
            raise ValueError("codomain of other must match domain of self")
        return Function(other.domain, self.codomain, lambda x: self(other(x)))
   
    
    def new_domains(self, domain, codomain, check_well_defined=True):
        return Function(domain, codomain, self.function, check_well_defined)




def identity(s):
    if not isinstance(s, Set):
        raise TypeError("s must be a set")
    return Function(s, s, lambda x: x)


def plus1(dom, cod):
    if not isinstance(dom, Set) or not isinstance(cod, Set):
        raise TypeError("both must be a set")
    return Function(dom, cod, lambda x: x+1)
    

    
if __name__ == '__main__':

    F2 = Function(Set({1,2,3}), Set({1,2,3,4}), lambda x: x+1)
    
    print("Imagen: ", F2._image())
    print("Print:\n", F2)
    print(F2.is_bijective())
    
    
    
