#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from beautifultable import BeautifulTable
from IPython.display import display, Image
import networkx as nx
from Permutation import permutation



gens = "aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ"


def ConvertToWord(w, s, N):
    for item in s:
        elem = gens.find(item)
        if elem >= N:
            return False
        w.append(elem)
    return True



def inv(x):
    return x+1 if x%2==0 else x-1


def replace_all(text, rep):
    for i, j in rep.items():
        text = text.replace(i,j)
    return text



class EquivalenceClass(object):
    """
    It controls the equivalence class of coset/vertexs
    We represent each vertex as the representant of it equivalence class.
    If a coincidence is found, we need to delete the coset with
    the higher value.
    """

    def __init__(self):
        self.p = [0]


    def __getitem__(self, index):
        return self.p[index]

    def __repr__(self):
        dev = " "
        for i in range(0,len(self.p)):
            dev = dev + str(self.p[i]) + " "
        return dev


    __str__ = __repr__

    def add(self):
        self.p.append(len(self.p))


    def rep(self, k):
        """
        k : coset (int)

        Returns the representant of the equivalence class of the coset K
        """

        l = k
        r = self.p[l]
        while r != l:
            l = r
            r = self.p[l]
        m = k
        r = self.p[m]
        while r != l:
            self.p[m] = l
            m=r
            r = self.p[m]
        return l

    def merge(self, k, l):
        """
        k : coset (int)
        l : coset (int)

        Returns the representent of both cosets or false
        if the cosets are equals.
        """

        v = self.rep(k)
        c = self.rep(l)

        if v != c:
            m = min(v,c)
            u = max(v,c)
            self.p[u]=m
            return True, u

        return False, c



class CosetTable(object):
    """
    This class stores a coset table that is equivalent to the
    Schreier graph that reflects the action of G on G/H.
    """


    def __init__(self, *arguments):
        """
        *arguments: admits multiples ways to read a group G given
            as generators and relators and a subgroup H given as generators.
        """

        l = len(arguments)
        gen = arguments[0][0] if l==1 else arguments[0]
        rel = arguments[0][1] if l==1 else arguments[1]
        gen_H = arguments[0][2] if l==1 else arguments[2]

        #variables to control memory
        self.M = 1E6
        self.n = 0
        self.cosets=[]

        self.table = None
        self.q = []
        self.relator = []
        self.generator_of_H = []

        self.NGENS = 2*len(gen) #Generators and inverses
        self.p = EquivalenceClass()


        self.tab = []
        self.tab.append([-1 for i in range(0, self.NGENS)])

        for i in range(0, len(gen_H)):
            w = []
            ConvertToWord(w, gen_H[i], self.NGENS)
            self.generator_of_H.append(w)

        for i in range(0, len(rel)):
            w = []
            s = rel[i]
            ConvertToWord(w,s, self.NGENS)
            self.relator.append(w)



    def finalCosets(self):
        """
        returns the index of [G:H], the number of living cosets.
        """
        cont=0
        for k in range(0,len(self.tab)):
            if self.isAlive(k):
                cont=cont+1
        return cont



    def usedCosets(self):
        """
        returns the total number of cosets used in the execution.
        """
        return len(self.tab)


    def isAlive(self, k):
        """
        k : coset (int)

        checks if coset k is alive, i.e., p[k]=k.
        """
        return self.p[k]==k


    def isDefined(self, k, x):
        """
        k : coset (int)
        x : element of A:= X^{+1} \cup X^{-1}

        returns True if x^k is defined.
        """
        return self.tab[k][x]>=0


    def undefine(self,k, x):
        """
        k : coset (int)
        x : element of A:= X^{+1} \cup X^{-1}

        undefine x^k.
        """
        self.tab[k][x] = -1


    def coset_table(self):
        """
        returns the coset table.
        """
        return self.table if self.table != None else self.tab


    def define(self, k, x):
        """
        Args:
            k : coset (int)
            x : generator of X or inverse.

        New definition, k^x = new available natural.
        """
        if self.n==self.M:
            raise MemoryError("Not enough memory.")
        self.n=self.n+1

        l = len(self.tab)
        self.tab[k][x]=l
        r = [-1 for i in range(0, self.NGENS)]
        r[inv(x)] = k

        self.tab.append(r)
        self.p.add()



    def merge(self,k, l):
        """
        Args:
            k : coset (int)
            l : coset (int)

        Returns the representent of both cosets or nothing
        if the cosets are equals.
        """

        var, u = self.p.merge(k,l)
        if var:
            self.q.append(u)  #q is queue


    def coincidence(self, a, b):
        """
        Args:
            a : coset (int)
            b : coset (int)

        It processes a coincidence between coset a and b
        """

        self.merge(a,b)

        while not len(self.q)==0: #mientras que no sea empty
            y = self.q[0]   #q.front()
            self.q.pop(0)   #q.pop

            for x in range(0, self.NGENS):
                if not self.isDefined(y,x):
                    continue

                d = self.tab[y][x]
                self.undefine(d,inv(x))
                m = self.p.rep(y)
                v = self.p.rep(d)

                if self.isDefined(m,x):
                    self.merge(v, self.tab[m][x])
                elif self.isDefined(v,inv(x)):
                    self.merge(m, self.tab[v][inv(x)])
                else:
                    self.tab[m][x] = v
                    self.tab[v][inv(x)] = m


    def scan_and_fill(self, a, w):
        """
        Args:
            a : coset (int)
            w : word (list)

        Scans the word w over the coset a.
        This method defines new cosets if needed.
        """

        i, j = 0, len(w)-1
        f, b= a, a

        while True:
            while i<=j and self.isDefined(f, w[i]):
                f=self.tab[f][w[i]]
                i = i+1

            if i>j:
                if f != b:
                    self.coincidence(f,b)
                return

            while j>=i and self.isDefined(b, inv(w[j])):
                b = self.tab[b][inv(w[j])]
                j = j-1

            if j<i:
                self.coincidence(f,b)
                return

            elif j==i:
                self.tab[f][w[i]]=b
                self.tab[b][inv(w[i])]=f
                return
            else:
                self.define(f,w[i])



    def CosetEnumeration(self):
        """
        HLT Method for Todd Coxeter algorithm.
        """

        for i in self.generator_of_H:
            self.scan_and_fill(0, i)


        k=0
        while k < len(self.tab):

            for item in self.relator:
                if not self.isAlive(k):
                    continue
                self.scan_and_fill(k, item)


            if(self.isAlive(k)):
                for x in range(0, self.NGENS):
                    if not self.isDefined(k,x):
                        self.define(k,x)


            k=k+1

        self.pretty_print()




    def pretty_print(self):
        """
        returns the Coset table
        """

        table = BeautifulTable()

        head= [gens[i] for i in range(0, self.NGENS)]
        table.rows.insert(0, head, "C")

        inserted=0
        for a in range(0, self.usedCosets()):
            if inserted==self.finalCosets():
                break

            l1=[]
            if self.isAlive(a):
                self.cosets.append(a)
                for b in range(0,len(self.tab[a])):
                   l1.append(self.tab[a][b]+1)

                table.rows.insert(a+1, l1, str(a+1))
                inserted=inserted+1

        table.set_style(BeautifulTable.STYLE_BOX)
        self.table = table
        return table



    def schreier_graph(self, notes=True):
        """
        Args:
            notes : boolean

        Return the Schreier graph of self.
        """

        if self.table == None:
            raise ValueError("Debe llamar primero a la función Todd Cox")

        l = [ list(self.table[i]) for i in range(1,len(self.cosets)+1)]
        vertexs = [ self.cosets[i]+1 for i in range(len(self.cosets))]
        cols = [[row[i] for row in l] for i in range (0,len(l[0]),2) ]
        colors = ['red','blue','green','black','orange','grey','purple','orange'] #max 8 gens

        Grafo = nx.MultiDiGraph(engine='circo')
        Grafo.graph['node']={'shape':'circle'}

        for i in vertexs:
            Grafo.add_node(i)
        for idx, i in enumerate(vertexs):
            for j in range(len(l[0])//2):
                Grafo.add_edge(i, cols[j][idx], color= colors[j], label=" "+ gens[2*j])
                if notes:
                    print("Arrow from", i, " to ", cols[j][idx], " coloured of ", colors[j])

        G = nx.nx_agraph.to_agraph(Grafo)
        G.layout('circo') #circo
        G.draw('T&C_SchreierGraph.png')
        display(Image('T&C_SchreierGraph.png'))



    def follow(self, c, d):
        """
        It follows the path of the cosets given and returns a permutation.

        Args:
            c : coset (int)
            d : coset (int)
        """
        c = self.p.rep(c)
        k = self.tab[c]

        return self.p.rep(k[d])



    def getGenerators(self):
        """
        Return the Schreier generators after executing the Todd Coxeter algorithm.
        We can use these generators to obtain a permutation representation of G.
        """

        perms=[]
        for g in range(self.NGENS//2):
            l=[]
            for i, c in enumerate(self.cosets):
                l.append(self.cosets.index(self.follow(c, 2*g)))
            perms.append(l)

        for i in range(len(perms)):
            for j in range(len(perms[i])):
                perms[i][j] = perms[i][j]+1

        gens=[]
        for i in range(self.NGENS//2):
            gens.append(permutation(perms[i]))

        return gens


def readGroup(file):
    """
    Read a group G given as generators and relators and a subgroup H given as generators.
    Args:
        file: a path

    Example:
        >>> f = readGroup("Groups/S3.txt")
        >>> print(f)
        (['a', 'b'], ['aaa', 'bb', 'abAAB'], [])
    """

    try:
        f = open(file, "r")
        rep = {",": " ", "=": " ", "1": " ", "{": " ", "}": " "}
        gen = replace_all(f.readline(), rep).split()
        rels = replace_all(f.readline(), rep).split()
        genH = replace_all(f.readline(), rep).split()

        f.close()
    except IOError:
        print("Could not read file ", file)

    return gen, rels, genH
