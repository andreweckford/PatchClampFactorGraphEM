#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 14:32:10 2020

@author: andreweckford
"""

import numpy as np

class Receptor:
    
    def __init__(self,c0,wt=False):

        # this is the Zhang model for acetylcholine with 4 states
        # grabbed it from:
        # J. Chen et al., "Activation kinetics of recombinant mouse 
        # nicotinic acetylcholine receptors: Mutations of \alpha-subunit
        # tyrosine 190 affect both binding and gating," Biophys. J., 1995.
        
        # State labels:
        # 0 = C
        # 1 = AC
        # 2 = A2C
        # 3 = A2O
        # Blocked state is ignored

        # not a great model for WT ... huge mismatch in parameter values

        bigDelta = 0.05 # discrete timestep in ms

        # concentrations for low (0) and high (1) input
        c = np.array([c0,1e-5])

        # rates at the reference are in seconds^-1, these are in ms^-1
        if (wt is True):
            kPlus1 = 0.021 # per micromol
            kMinus1 = 0.650
            kPlus2 = 0.04 # per micromol
            kMinus2 = 25
            alpha1 = 60
            beta1 = 0.24
        else:
            # mutants are \alpha Y190F mutant given at above reference
            kPlus1 = 0.00049
            kMinus1 = 0.366
            kPlus2 = 0.002
            kMinus2 = 2.43
            alpha1 = 0.150
            beta1= 0.5
            

        # xA_array = np.array([1e-7,1e-5])

        i = 0
        R0 = np.array([
            [-1*kPlus1*c[i],kPlus1*c[i],0.,0.],
            [kMinus1,-1*(kMinus1 + kPlus2*c[i]),kPlus2*c[i],0.],
            [0.,kMinus2,-1*(kMinus2 + alpha1),alpha1],
            [0.,0.,-1*beta1,beta1]
        ])        

        i = 1
        R1 = np.array([
            [-1*kPlus1*c[i],kPlus1*c[i],0.,0.],
            [kMinus1,-1*(kMinus1 + kPlus2*c[i]),kPlus2*c[i],0.],
            [0.,kMinus2,-1*(kMinus2 + alpha1),alpha1],
            [0.,0.,-1*beta1,beta1]
        ])
        
        self.P0 = np.eye(4) + bigDelta * R0
        self.P1 = np.eye(4) + bigDelta * R1
        
        self.Pmask = np.array([
            [1.,1.,0.,0.],
            [1.,1.,1.,0.],
            [0.,1.,1.,1.],
            [0.,0.,1.,1.]
        ])
        
        self.statemap = np.array([0.,0.,0.,1.])
