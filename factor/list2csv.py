#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 12:39:38 2020

@author: andreweckford
"""


def list2csv(l):
    for i in range(0,len(l)-1):
        print(str(l[i])+',',end='')
    print(str(l[-1]))
    
def printP(P):
    numStates = P.shape[0]
    for i in range(0,numStates):
        list2csv(P[i,:])
        
def listOfLists2csv(l):
    for i in range(0,len(l)):
        list2csv(l[i])