#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 08:59:31 2020

@author: andreweckford
"""

print(__name__)
if __package__ is None:
    print("It's none!")
else:
    print("It's not none!")
    
import sys
import csv
from .. import argvHandler as a

# read from a CSV file, if filename is none read from stdin
# the file should be structured as follows:
# - (if readParams == True) ... first 2 lines are parameters
# - (if readEstimates == True) ... next several lines are data
# -- first line in the data is the true state values
# -- next line is the state estimate with known parameters
# -- next line is the state estimate with known parameters, with non-confident estimates blanked out (with -1)
# -- next: if lastOnly is True, two more lines with the EM estimates; if lastOnly is False, two lines for each EM iteration
# - (if readP == True) ... the matrix of transition probability estimates
def parseResults(filename = None, readParams = True, readEstimates = True, readP = True):
    
    if filename is None:
        r = csv.reader(sys.stdin.readlines())
    else:
        with open(filename) as csvfile:
            r = csv.reader(csvfile)

    # create a big list of everything in the csv
    l = []
    for i in r:
        l.append(list(i))

    print(l)

    #params = a.clp.paramsFromCSV(l[0],l[1])
    
    
    
    # # 
    
    # m = None
    # for i in r:
    #     z = np.array([float(k) for k in list(i)])
    #     if m is None:
    #         m = z
    #     else:
    #         m = np.vstack((m,z))

    # em_iter = int((m.shape[0]-3)/2)
    # numTimeInstants = int(m.shape[1])

    # kp_err = 0
    # em_err = np.zeros(em_iter)
    # em_conf_err = np.zeros(em_iter)
    # for i in range(0,m.shape[1]):
    #     if m[0,i] != m[1,i]:
    #         kp_err += 1
    #     for k in range(0,em_iter):
    #         if (m[0,i] != m[2*k+3,i]):
    #             em_err[k] += 1

    # r_em[q] += em_err[-1]
    # r_kp[q] += kp_err
    
