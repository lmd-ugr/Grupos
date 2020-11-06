import math
import numpy as np

class Dihedral:
    """
    Stores rotations and reflections to create Dihedral Group.
    """
    def __init__(self, n):
        """
        It creates n rotations and n reflections.
        
        Args:
            n : integer
        
        Example:
            >>>D = Dihedral(2)
            Rotaciones: 
            (1.0, -0.0, 0.0, 1.0)
            (-1.0, -0.0, 0.0, -1.0)
            
            Reflexiones 
            (1.0, 0.0, 0.0, -1.0)
            (-1.0, 0.0, 0.0, 1.0)
        """

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
        
        dev = "Rotaciones: \n"
        for i in self.rot:
            dev = dev+  str(i)+ "\n"

        dev = dev+"\nReflexiones \n" 
        for i in self.refl:
            dev = dev + str(i)+ "\n"

        return dev
        #return "Rotaciones: " + str(self.rot) + "\nTraslaciones: "+str(self.refl)
        

    __str__ = __repr__

    def __eq__(self, other):
        """ 
        Checks if self is equal to other.
        
        Args:
            other : Dihedral instance.
        """
        if not isinstance(other, Dihedral):
            raise TypeError(other, "is not a Dihedral Group")

        return self.rot==other.rot and self.refl==other.refl

    def __ne__(self, other):
        """ 
        Checks if self is equal or not to other.
        
        Args:
            other : Dihedral instance.
        """
        return not self==other


    def get_rot(self,i):
        return self.rot[i]

    def get_rots(self):
        return self.rot
    
    
    def get_refls(self):
        return self.refls

    def get_refl(self,i):
        return self.refl[i]

    def get_all(self):
        return self.all



    def product_RS(self, a, b):
        """
        Product if representation is RS
        
        Args:
            a : rotation or reflection
            b : rotation or reflection
        """

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
        """
        Product if representation is matrix.
        
        Args:
            a : rotation or reflection
            b : rotation or reflection
        """
        
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

