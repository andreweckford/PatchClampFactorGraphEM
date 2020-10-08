#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 08:59:31 2020

@author: andreweckford
"""

import sys
import csv
import numpy as np
from factor.mainArgvHandler import clp

class Results:
    
    def __init__(self):
        self.params = None
        self.states = None
        self.kp = None
        self.kp_conf = None
        self.emEstimates = None
        self.emEstimates_conf = None
        self.last = None
        self.last_conf = None
        self.P = None

    # read from a CSV file, if filename is none read from stdin
    # the file should be structured as follows:
    # - (if readParams == True) ... first 2 lines are parameters
    # - (if readEstimates == True) ... next several lines are data
    # -- first line in the data is the true state values
    # -- next line is the state estimate with known parameters
    # -- next line is the state estimate with known parameters, with non-confident estimates blanked out (with -1)
    # -- next: if lastOnly is True, two more lines with the EM estimates; if lastOnly is False, two lines for each EM iteration
    # - (if readP == True) ... the matrix of transition probability estimates
    def parseResults(self,filename = None, reader = None, readParams = True, readEstimates = True, readP = True):
        
        if reader is not None:
            r = reader
        elif filename is not None:
            with open(filename) as csvfile:
                r = csv.reader(csvfile)
        else:
            r = csv.reader(sys.stdin.readlines())
    
        # create a big list of everything in the csv
        l = []
        for i in r:
            l.append(list(i))
    
        currentIndex = 0
    
        # what do we do if readParams is not true?
        # use the default values I guess?
        if (readParams is True):
            self.params = clp.paramsFromCSV(l[0],l[1])
            currentIndex = 2
        else:
            self.params = dict(clp.initParams) # remember to use the copy constructor
            
        if (readEstimates is True):
            
            # known states
            self.states = l[currentIndex]
            
            # estimates with known parameters
            self.kp = l[currentIndex + 1]
            
            # estimates with known parameters, non-confident estimates blanked out
            self.kp_conf = l[currentIndex + 2]
            
            currentIndex += 3
            
            self.emEstimates = []
            self.emEstimates_conf = []
            
            # parsed parameters are strings here
            # in the future, maybe cast to boolean?
            if (self.params["lastOnly"] == "True"):
                numEstimates = 1
            else:
                numEstimates = self.params["maxEMIterations"]
    
            for i in range(0,numEstimates):
                self.emEstimates.append(l[currentIndex])
                self.emEstimates_conf.append(l[currentIndex + 1])
                currentIndex += 2
    
            # the last EM estimates are particularly important ... store here
            self.last = l[currentIndex - 2]
            self.last_conf = l[currentIndex - 1]            
        
        if (readP is True):
            # what's the receptor?
            q = clp.createReceptor(self.params["receptor"],self.params["receptorParameter"])
            numStates = len(q.statemap)
            for i in range(0,numStates):
                newLine = np.array([float(j) for j in l[currentIndex]])
                if (i == 0):
                    self.P = newLine
                else:
                    self.P = np.vstack((self.P,newLine))
                currentIndex += 1
    
    def __errorsHelper(self,v):
        e = 0
        for i in range(0,len(v)):
            if (v[i] != self.states[i]):
                e += 1
                
        return e
    
    # only count when the state vector is in onlyCount
    # onlyCount should be a *list*, even if it has only one element
    def __errorsHelperOnlyCount(self,v,onlyCount):
        e = 0
        n = 0
        for i in range(0,len(v)):
            if (self.states[i] in onlyCount):
                n += 1
                if (v[i] != self.states[i]):
                    e += 1
                    
        return np.array([e,n])
    
    # don't count when the miss vector is in missCount (e.g., -1)
    # missCount should be a *list*, even if it has only one element
    def __errorsHelperMissCount(self,v,missCount):
        e = 0
        n = 0
        for i in range(0,len(v)):
            if (v[i] not in missCount):
                n += 1
                if (v[i] != self.states[i]):
                    e += 1
                    
        return np.array([e,n])
                
    def kpErrors(self):
        return self.__errorsHelper(self.kp)
    
    def emErrors(self):
        r = np.zeros(len(self.emEstimates))
        for i in range(0,len(self.emEstimates)):
            r[i] = self.__errorsHelper(self.emEstimates[i])
        return r
    
    # conf count is a 2-element array, first element is errors, second is count of confident estimates
    def kpConfErrors(self):
        return self.__errorsHelperMissCount(self.kp_conf,['-1.0'])

    def emConfErrors(self):
        n = []
        for i in range(0,len(self.emEstimates)):
            n.append(self.__errorsHelperMissCount(self.emEstimates[i],'-1.0'))
        return n

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
    
