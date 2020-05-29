#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 09:46:58 2020

@author: andreweckford
"""

import numpy as np
from myTools import randIntNonUniform, getSteadyStateDist

# assumed binary
class IIDChannelInputs:
    
    def __init__(self,px):
        self.px = px
        
    def getChannelInputs(self,n):
        r = np.zeros(n)
        for i in range(0,n):
            r[i] = randIntNonUniform(self.px)
            
        return r

class FakeACHSystem:
    
    def __init__(self):
        self.P0 = np.array(
            [[0.8, 0.1, 0.,  0.1, 0.],
              [0.1, 0.8, 0.1, 0.,  0.],
              [0.,  0.1, 0.8, 0.1, 0.],
              [0.1, 0.,  0.1, 0.7, 0.1],
              [0.,  0.,  0.,  0.1, 0.9]])
        self.P1 = np.array(
            [[0.1, 0.8, 0.,  0.1, 0.],
              [0.1, 0.8, 0.1, 0.,  0.],
              [0.,  0.1, 0.8, 0.1, 0.],
              [0.1, 0.,  0.5, 0.3, 0.1],
              [0.,  0.,  0.,  0.6, 0.4]])
        self.statemap = [1,1,0,0,0]


class ACHSystem:
    
    def __init__(self):
        dt = 1e-5
        self.P0 = np.array(
            [[-3050.,50.,0.,3000.,0.],
             [0.66,-500.66,500,0.,0.],
             [0.,15000.,-19000,4000.,0.],
             [15.,0.,50.,-2065.,2000.],
             [0.,0.,0.,10.,-10.]])
        self.P0 = np.eye(5) + dt*self.P0
        self.P1 = np.array(
            [[-8000.,5000.,0.,3000.,0.],
             [0.66,-500.66,500,0.,0.],
             [0.,15000.,-19000,4000.,0.],
             [15.,0.,5000.,-7015.,2000.],
             [0.,0.,0.,1000.,-1000.]])
        self.P1 = np.eye(5) + dt*self.P1
        self.statemap = [1,1,0,0,0]

class ThreeStateChain:
    
    def __init__(self):
        self.P0 = np.array(
            [[0.9,0.1,0.],
             [0.1,0.8,0.1],
             [0.,0.1,0.9]])
        self.P1 = np.array(
            [[0.5,0.5,0.],
             [0.1,0.8,0.1],
             [0.,0.1,0.9]])
        self.statemap = [0,0,1]
        
        
class Receptor:
    
    def __init__(self,receptorParams,px):
        self.P0 = receptorParams.P0
        self.P1 = receptorParams.P1
        self.statemap = receptorParams.statemap
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
        
if (__name__ == "__main__"):
    xx = IIDChannelInputs(np.array([0.9,0.1]))
    ss = ACHSystem(xx.px)
        
    x = xx.getChannelInputs(1000)
    r,c = ss.getReceptorState(x)
    print(x)
    print(r)
    print(c)
        
        
        
        
        
        
