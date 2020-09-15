#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Permutation import permutation
from Set import Set
from Function import Function
from Group import Group, GroupElem, GroupHomomorphism
from beautifultable import BeautifulTable

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



class ClassEquivalence(object):
    
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
        v = self.rep(k)
        c = self.rep(l)
        
        if v != c:
            m = min(v,c)
            u = max(v,c)
            self.p[u]=m
            return True, u
        
        return False, c
        


class CosetTable(object):
    

    def __init__(self, total, gen, rel, gen_H):
        
        #variables to control memory
        self.M = 1E8
        self.n = 0
        self.cosets=[]
        
        self.table = None
        self.q = []
        self.relator = []
        self.generator_of_H = []
        
        ############################
        
        self.NGENS = total #Generators and inverses
        self.p = ClassEquivalence()
        
        
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
        cont=0
        for k in range(0,len(self.tab)):
            if self.isAlive(k):
                cont=cont+1
        return cont
       
    
    
    def usedCosets(self):
        return len(self.tab)
    
    def isAlive(self, k):
        return self.p[k]==k
    
    def isDefined(self, k, x):
        return self.tab[k][x]>=0
    
    def undefine(self,k, x):
        self.tab[k][x] = -1
    
    def define(self, a, x):
        
        if self.n==self.M:
            raise MemoryError("Not enough memory.")
        self.n=self.n+1
        
        l = len(self.tab)
        self.tab[a][x]=l
        r = [-1 for i in range(0, self.NGENS)]
        r[inv(x)] = a
        
        
        self.tab.append(r)
        self.p.add()
        
 
    
    def merge(self,k, l):
        var, u = self.p.merge(k,l)
        if var:
            self.q.append(u)  #q is queue

      
    def coincidence(self, a, b):
        
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
        



    def pretty_print(self):
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



    def pretty_print2(self):
        print(self.tab)
        
        
    def schreier_graph(self):
        if self.table == None:
            raise ValueError("Debe llamar primero a la función Todd Cox")
        
        cols = [self.table.columns[i] for i in range(0, self.getnlive()) if i%2==0]
        table=[]
        for i in range(1, len(cols[0])):
            l = []
            for j in range(0, len(cols)):
                l.append(cols[j][i])
            table.append(l)
        
        '''
        import networkx as nx
        import matplotlib.pyplot as plt
        import graphviz as gv
        G = gv.Digraph(format='png', engine='circo')
        
        for i in range(0, self.getnlive()):
            G.node(i)
        '''
            
        
        print(table)
        '''
        A = nx.nx_agraph.tp_agraph(G)
        A.layout('dot')
        A.draw('salida.png')
        '''    
    
    
    
    def find(self,c):
        c2 = self.p[c]
        if c == c2:
            return c
        else:
            c2 = self.find(c2)
            self.p[c] = c2
            return c2
    
    def follow(self, c, d):
        c = self.find(c)
        ns = self.tab[c]
        
        return self.find(ns[d])
            



def generate2(gens):
        
    idn = permutation([1])
    
    set_element_list = set([idn])
    dev = []
    for i in range(len(gens)):
        D = list(set_element_list)
        N = [idn]
        while N:
            A = N
            N = []
            for a in A:
                for g in gens[:i + 1]:
                    ag = a*g
                    #print(ag , " =| " , a , " * " , g)
                    if ag not in set_element_list:
                        #print(ag , " not in " , set_element_list)
                        # produce G_i*g
                        for d in D:
                            ap = d*ag
                            dev.append(ap)
                            set_element_list.add(ap)
                            N.append(ap)
   
    S = Set(dev) if len(dev) != 0 else Set({1})
    bin_op = Function(S*S, S, lambda x: x[0]*x[1])
    G = Group(S, bin_op)
    G.group_gens=gens
    return G



def generate(elems):

    oldG = Set(elems)
    
    while True:
        newG = oldG | Set(a * b for a in oldG for b in oldG)
        if oldG == newG: 
            break
        else: 
            oldG = newG
    C = Set(oldG)
    bin_op = Function(C*C, C, lambda x: x[0]*x[1])
    
    G = Group(C, bin_op)
    G.group_gens = [ GroupElem(g, G) for g in elems ]
   
    return G   
    

if __name__ == "__main__":
    
    file = "Groups/D10.txt"
    try:
        
        f = open(file, "r")
        rep = {",": " ", "=": " ", "1": " ", "{": " ", "}": " "}
        gen = replace_all(f.readline(), rep).split()
        rels = replace_all(f.readline(), rep).split()
        genH = replace_all(f.readline(), rep).split()
        
        nueva = CosetTable(2*len(gen), gen, rels, genH)
        nueva.CosetEnumeration()
        
        big = CosetTable(2*len(gen), gen, rels, [])
        big.CosetEnumeration()

        print("The order of G is {}".format(big.finalCosets()))
        print("The index [G:H] is {}".format(nueva.finalCosets()))
        print("The order of H is {}".format(big.finalCosets()//nueva.finalCosets()))
        print("Number of used cosets: " , nueva.usedCosets())
        
       
        table = nueva.pretty_print()
        if nueva.finalCosets() <= 25:
            pass
            print(table)
           
        '''
        En nueva.cosets almacenamos el índice de los vértices que están vivos.
        Estos no tienen por qué ser consecutimos, puede ser por ejemplo:
            0,1,4,5,6,8,9
        
        Se tiene la estructura del grafo de Schreier y a partir de las aristas
        del grafo es posible construir una permutación
        
        En nueva.cosets se almacenan las clases que son válidas!
        recorro cada una de ellas llamando a la función follow.
        Esta función toma un vértice c y encuentra al vecino en la dirección 2*g
        '''
        perms=[]
        for g in range(len(gen)):
            l=[]
            for i, c in enumerate(nueva.cosets):
                l.append(nueva.cosets.index(nueva.follow(c, 2*g)))
            perms.append(l)
             
        #perms = [[nueva.cosets.index(nueva.follow(c, 2*g)) for i, c in enumerate(nueva.cosets)] for g in range(len(gen))]
 
        for i in range(len(perms)):
            for j in range(len(perms[i])):
                perms[i][j] = perms[i][j]+1
                
        print("\nGenerators of G:")
        gens=[]
        for i in range(len(gen)):
            gens.append(permutation(perms[i]))
            print("g{} = {}".format(i, permutation(perms[i])))
        
        
        print("\n")
        G = generate(gens)
        
        print("Orden del grupo: ", G.order())
        #print(G.Cayley_table())
        #print(G.elements_order())


        #S = SymmetricGroup(3)
        #print(S.is_isomorphic(S))
        '''
        I = S.find_isomorphism(G)
        if I==None:
            print("No son isomorfos")
        elif I.is_isomorphism():
            print("Son isomorfos")
        '''
        f.close()
        
    except IOError:
        print("Could not read file ", file)