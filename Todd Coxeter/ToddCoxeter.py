#!/usr/bin/env python3
# -*- coding: utf-8 -*-



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
        from beautifultable import BeautifulTable
        head= [gens[i] for i in range(0, self.NGENS)]
        
        table = BeautifulTable()
        
        table.rows.insert(0, head, "Coset")
        
        #for a in range(0, len(self.tab)):
        yea=False
        less=False
        for a in range(0, self.finalCosets()):
           l1 = []
           
           
           if a==self.finalCosets()-1:
               yea = any(i==-1 for i in self.tab[a])
               
           for b in range(0,len(self.tab[a])):
               if yea:
                   a=a+1
                   yea=False
                   less=True
                   
                   
               l1.append(self.tab[a][b]+1)

           if less:   
               table.rows.insert(a+1, l1, str(a))
           else:
               table.rows.insert(a+1, l1, str(a+1))

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
            


          

    

if __name__ == "__main__":
    
    file = "Groups/Q2.txt"
    try:
        f = open(file, "r")
        rep = {",": " ", "=": " ", "1": " ", "{": " ", "}": " "}
        gen = replace_all(f.readline(), rep).split()
        rels = replace_all(f.readline(), rep).split()
        genH = replace_all(f.readline(), rep).split()

    
        nueva = CosetTable(len(gen)*2, gen, rels, genH)
        nueva.CosetEnumeration()
        print("The index of H in G is {}".format(nueva.finalCosets()))
        print("Number of used cosets: " , nueva.usedCosets())
        
        if nueva.finalCosets() <= 25:
            print(nueva.pretty_print())
            #print(nueva.schreier_graph())
        
        f.close()
    except IOError:
        print("Could not read file ", file)
    



