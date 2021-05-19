#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:21:58 2020

@author: andreweckford
"""

import numpy as np
from .myTools import randIntNonUniform,getSteadyStateDist

class Simulator:
    
    def __init__(self,receptor):
        self.P0 = receptor.P0
        self.statemap = receptor.statemap
        self.px = None
    
    # returns two lists, states and ion channels, given input sequence
    def getReceptorState(self,x):
        n = len(x)
        r = np.zeros(n)
        c = np.zeros(n)
        initState = randIntNonUniform(getSteadyStateDist(self.P0))
        for i in range(0,n):
            if (i > 0):
                r[i] = randIntNonUniform(self.P0[int(r[i-1])])
            else:
                r[i] = randIntNonUniform(self.P0[int(initState)]) 
             
            c[i] = self.statemap[int(r[i])]
            
        return r,c
    
    # returns three lists, states, ion channels, and noise-corrupted ion channels
    # noise variance is sigma2 -- this is scaled to the current, which is assumed to be unit
    def getReceptorStateNoisy(self,x,sigma2):
        
        r,c = self.getReceptorState(x)
        y = np.sqrt(sigma2)*np.random.randn(len(x)) + c
        
        return r,c,y
