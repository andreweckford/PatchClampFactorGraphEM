#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:43:29 2020

@author: andreweckford
"""

import numpy as np

# Edit this file to change the CFTR parameters.
# P0 and P1 are expressed as transition probability matrices:
#  - each entry is the probability of the given transition
#    i.e. P0[3,4] is Pr(next state is O1 given current state is C2)
#  - rows must sum to 1 by the rules of probability

class Receptor:
    
    def __init__(self,C1aExitProb=0.1,TwoConductanceStates=False):

        self.P0 = np.array(
            [[1.-C1aExitProb, C1aExitProb, 0.,  0.,  0.,  0.,  0.],   # C1a
              [0.1, 0.8, 0.1, 0.,  0.,  0.,  0.],   # C1b
              [0.,  0.1, 0.8, 0.1, 0.,  0.,  0.],   # C2
              [0.,  0.,  0.1, 0.8, 0.1, 0.,  0.],   # O1
              [0.,  0.,  0.,  0.,  0.9, 0.1, 0.],   # O2
              [0.,  0.,  0.,  0.,  0.1, 0.8, 0.1],  # C3
              [0.1, 0.,  0.,  0.,  0.,  0.1, 0.8]]) # C4
        
        self.P1 = np.array(
            [[0.9, 0.1, 0.,  0.,  0.,  0.,  0.],   # C1a
              [0.1, 0.8, 0.1, 0.,  0.,  0.,  0.],   # C1b
              [0.,  0.1, 0.8, 0.1, 0.,  0.,  0.],   # C2
              [0.,  0.,  0.1, 0.8, 0.1, 0.,  0.],   # O1
              [0.,  0.,  0.,  0.,  0.9, 0.1, 0.],   # O2
              [0.,  0.,  0.,  0.,  0.1, 0.8, 0.1],  # C3
              [0.1, 0.,  0.,  0.,  0.,  0.1, 0.8]]) # C4
        
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
