#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:43:29 2020

@author: andreweckford
"""

import numpy as np
from scipy.linalg import expm

# Edit this file to change the CFTR parameters.
# P0 and P1 are expressed as transition probability matrices:
#  - each entry is the probability of the given transition
#    i.e. P0[3,4] is Pr(next state is O1 given current state is C2)
#  - rows must sum to 1 by the rules of probability

class Receptor:
    
    def __init__(self,C1aExitProb=0.1,TwoConductanceStates=False):

        # matrix of rates
        R = np.array(
            [[-9., 9.,    0.,    0.,    0.,   0.,    0.],     # C1a
             [5.,  -12.7, 7.7,   0.,    0.,   0.,    0.],     # C1b
             [0.,  5.8,   -10.7, 4.9,   0.,   0.,    0.],     # C2
             [0.,  0.,    10.,   -17.1, 7.1,  0.,    0.],     # O1
             [0.,  0.,    0.,    0.,    -3.,  3.,    0.],     # O2
             [0.,  0.,    0.,    0.,    7.,   -13.,  6.],     # C3
             [1.7, 0.,    0.,    0.,    0.,   12.8,  -14.5]]) # C4
        
        dt = 0.01

        # we specify C1aExitProbability to stay compatible with earlier versions
        # but what it really means is the out rate from C1a normalized by dt
        R[0,0] = -1*C1aExitProb/dt
        R[0,1] = C1aExitProb/dt        

        #self.P0 = np.eye(7) + R*dt
        self.P0 = expm(dt*R)

        # C1aExitProb is [ATP] * 9000
        self.P1 = None # self.P1 used to have meaning but is now unused
        
        self.Pmask = np.array(
            [[1., 1., 0.,  0.,  0.,  0.,  0.],   # C1a
              [1., 1., 1., 0.,  0.,  0.,  0.],   # C1b
              [0.,  1., 1., 1., 0.,  0.,  0.],   # C2
              [0.,  0.,  1., 1., 1., 0.,  0.],   # O1
              [0.,  0.,  0.,  0.,  1., 1., 0.],   # O2
              [0.,  0.,  0.,  0.,  1., 1., 1.],  # C3
              [1., 0.,  0.,  0.,  0.,  1., 1.]]) # C4)
        
        if (TwoConductanceStates is True):        
            self.statemap = [0.,0.,0.,1.,2.,0.,0.]
        else:
            self.statemap = [0.,0.,0.,1.,1.,0.,0.]
