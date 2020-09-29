#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 08:29:17 2020

@author: andreweckford
"""

import numpy as np

numTimeInstants_v = [20000]
maxEMIterations_v = [400]
confidence_v = [0.8]
C1aXP_v = np.arange(0.05,0.2,0.002)
numRuns = 4
num_cores = 4
coreID = 3
confidenceDigits = 1
C1aXPDigits = 3

counter = 0


# string formatting - floats, returns string with given number of 
# digits after the decimal
# f should be less than 1
def mysf(f,digits):
    r = str(int(np.round(f*10**digits)))
    
    # pad with leading zeros
    while (len(r) < digits):
        r = '0'+r
        
    return r
    


for numTimeInstants in numTimeInstants_v:
    for maxEMIterations in maxEMIterations_v:
        for confidence in confidence_v:
            for C1aXP in C1aXP_v:
                for run in range(0,numRuns):
                
                    if (counter % num_cores == coreID):
                        scriptString = 'python3 main.py ' + str(numTimeInstants)
                        scriptString += ' ' + str(maxEMIterations) 
                        scriptString += ' 0.' + mysf(confidence,confidenceDigits)
                        scriptString += ' 0.' + mysf(C1aXP,C1aXPDigits) 
                        scriptString += ' -l > result' + str(run) 
                        scriptString += '_' + str(numTimeInstants)
                        scriptString += '_' + str(maxEMIterations) 
                        scriptString += '_0' + mysf(confidence,confidenceDigits)
                        scriptString += '_0' + mysf(C1aXP,C1aXPDigits) + '_l.csv'
                        
                        print(scriptString)

                    counter += 1
