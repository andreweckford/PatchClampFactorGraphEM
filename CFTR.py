#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:43:29 2020

@author: andreweckford
"""

import numpy as np

class Receptor:
    
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
