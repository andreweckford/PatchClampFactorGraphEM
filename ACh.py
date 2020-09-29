#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 14:32:10 2020

@author: andreweckford
"""

import numpy as np

class Receptor:
    
    def __init__(self):
        
        # State numbering: (labels from Colquhoun and Hawkes, 1982)
        # 0 = AR
        # 1 = A2R
        # 2 = A2T
        # 3 = AT
        # 4 = T
        
        #maps = [[[0],[1],[2],[3],[4]],[[0,1],[2,3,4]]]
        
        bigDelta = 0.05 # discrete timestep in ms
        
        # concentrations for low (0) and high (1) input
        c = np.array([2e-5,1e-5])
    
        # parameter values taken from Colquhoun and Hawkes, 1982
        # rates at the reference are in seconds, these are in ms
        kPlus1 = 50000 * c # per mol
        kMinus1 = 2 # insensitive to concentration
        kPlus2 = 500000 * c # per mol
        kMinus2 = 2 # insensitive to concentration
        beta1 = 0.015 # insensitive
        alpha1 = 3 # insensitive
        beta2 = 15 # insensitive
        alpha2 = 0.5 # insensitive
        kStarPlus2 = 50000 * c # per mol
        kStarMinus2 = 0.00033 # insensitive
                
        # xA_array = np.array([1e-7,1e-5])
        
        i = 0
        R0 = np.array([
            [-1*(alpha1+kStarPlus2[i]),kStarPlus2[i],0,alpha1,0],
            [2*kStarMinus2,-1*(alpha2 + 2*kStarMinus2),alpha2,0,0],
            [0,beta2,-1*(beta2+2*kMinus2),2*kMinus2,0],
            [beta1,0,kPlus2[i],-1*(beta1+kPlus2[i]+kMinus1),kMinus1],
            [0,0,0,2*kPlus1[i],-2*kPlus1[i]]
        ])
        

        i = 1
        R1 = np.array([
            [-1*(alpha1+kStarPlus2[i]),kStarPlus2[i],0,alpha1,0],
            [2*kStarMinus2,-1*(alpha2 + 2*kStarMinus2),alpha2,0,0],
            [0,beta2,-1*(beta2+2*kMinus2),2*kMinus2,0],
            [beta1,0,kPlus2[i],-1*(beta1+kPlus2[i]+kMinus1),kMinus1],
            [0,0,0,2*kPlus1[i],-2*kPlus1[i]]
        ])
        
        self.P0 = np.eye(5) + bigDelta * R0
        self.P1 = np.eye(5) + bigDelta * R1
        
        self.Pmask = np.array([
            [1.,1.,0.,1.,0.],
            [0.,1.,1.,0.,0.],
            [0.,1.,1.,1.,0.],
            [1.,0.,1.,1.,1.],
            [0.,0.,0.,1.,1.]
        ])
        
        self.statemap = np.array([0.,0.,1.,1.,1.])