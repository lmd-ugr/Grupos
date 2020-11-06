from Set import Set


class Function:
    """
    Definition of a finite function.
    """
    
    def __init__(self, domain, codomain, function, check_well_defined=True):
        """
        Initialize the function and check that it is well-formed.

        This method can be overwritten by subclasses of Function, so that for
        example GroupHomomorphisms can be between Groups, rather than Sets.
        
        Args:
            domain : instance of Set, domain of self.
            codomain : instance of Set, codomain of self.
            function : function of self.
            check_well_defined : boolean.
        """
        if not isinstance(domain, Set):
            raise TypeError("Domain must be a Set")
        if not isinstance(codomain, Set):
            raise TypeError("Codomain must be a Set")
        if check_well_defined:
            if not all(function(elem) in codomain for elem in domain):
                raise TypeError("Function returns some value outside of codomain")

        self.domain = domain
        self.codomain = codomain
        self.function = function

    def __call__(self, elem):
        """ 
        Method to call an instance.
        Args:
            elem : instance of Function.
        """
        if elem not in self.domain:
            raise TypeError("Function must be called on elements of the domain")
        return self.function(elem)

    def __hash__(self):
        """
        Returns the hash of self.
        """

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
        """ 
        Checks if self and other are equals.
        
        Args:
            other : instance of Function.
        """
        if not isinstance(other, Function):
            return False

        return id(self) == id(other) or ( \
               self.domain == other.domain and \
               self.codomain == other.codomain and \
               all(self(elem) == other(elem) for elem in self.domain) )

    def __ne__(self, other):
        """ 
        Checks if self and other are not equals.
        
        Args:
            other : instance of Function.
        """
        return not self == other

    def _image(self):
        """
        The literal image of the function.
        """
        return Set(self(elem) for elem in self.domain)

    def image(self):
        """
        The API image of the function; can change depending on the subclass.
        For example, GroupHomomorphisms return the image as a Group, not a Set.
        """
        return self._image()

    def __str__(self):
        """
        Pretty outputing of functions.
        """

        l1 = []
        for i in self.domain:
            l1.append( "f(" + str(i)+ ")=" + str(self(i)) + "\n")
        
        return ' '.join(l1)
    
    
    def is_surjective(self):
        """
        returns if self is surjective.
        """
        return self._image() == Set(self.codomain)

    def is_injective(self):
        """
        returns if self is injective.
        """
        return len(self._image()) == len(self.domain)

    def is_bijective(self):
        """
        returns if self is bijective.
        """
        return self.is_surjective() and self.is_injective()

    def compose(self, other):
        """
        Returns x -> self(other(x))
        
        Args:
            other : instance of Function.
        """
        if not self.domain == other.codomain:
            raise ValueError("codomain of other must match domain of self")
        return Function(other.domain, self.codomain, lambda x: self(other(x)))

    def new_domains(self, domain, codomain, check_well_defined=True):
        """ 
        Redefines the domain and codomain.
        """
        return Function(domain, codomain, self.function, check_well_defined)


def identity(s):
    """
    Returns the identity function on the set s.
    
    Args: 
        s : instance of Set.
    """
    if not isinstance(s, Set):
        raise TypeError("s must be a set")
    return Function(s, s, lambda x: x)
