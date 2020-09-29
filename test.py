#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 16:25:37 2020

@author: albduranlopez
"""


from Group import *

Q = QuaternionGroup(rep="ijk")
sub_Q = Q.all_subgroups()


i = Quaternion(0,1,0,0)
j = Quaternion(0,0,1,0)


G = Q.generate([i])
G2 = Q.generate([j])

    
#print(G.is_normal_subgroup(Q))

#print(Q.all_subgroups())
#print(Q.subgroups_of_orderN(8))
#print(Q.all_normal_subgroups())

D = DihedralGroup(2, "RS")
r = D('R1')
R = D.generate([r])

#print(D.all_normal_subgroups())
#print((D/R).identity())



J = CyclicGroup(3)
print(J.is_inner_direct_product())
#print(D.is_inner_direct_product())
#print(G.is_inner_semidirect_product())
#print(G.is_inner_direct_product())












