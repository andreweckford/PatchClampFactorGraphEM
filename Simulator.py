#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:21:58 2020

@author: andreweckford
"""

import numpy as np
from myTools import randIntNonUniform,getSteadyStateDist

class Simulator:
    
    def __init__(self,receptor,px):
        self.P0 = receptor.P0
        self.P1 = receptor.P1
        self.statemap = receptor.statemap
        self.px = px
    
    # returns two lists, states and ion channels, given input sequence
    def getReceptorState(self,x):
        n = len(x)
        r = np.zeros(n)
        c = np.zeros(n)
        initState = randIntNonUniform(getSteadyStateDist(self.px[0] * self.P0 + self.px[1] * self.P1))
        for i in range(0,n):
            if int(x[i]) == 0:
                if (i > 0):
                    r[i] = randIntNonUniform(self.P0[int(r[i-1])])
                else:
                    r[i] = randIntNonUniform(self.P0[int(initState)]) 
            else:
                if (i > 0):
                    r[i] = randIntNonUniform(self.P1[int(r[i-1])])
                else:
                    r[i] = randIntNonUniform(self.P1[int(initState)])
             
            c[i] = self.statemap[int(r[i])]
            
        return r,c