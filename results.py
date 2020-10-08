#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:03:39 2020

@author: andreweckford
"""

from resultHandling.resultHandler import Results
from resultHandling.resultArgvHandler import ResultCLP
from factor.list2csv import list2csv,printP,listOfLists2csv
import sys


def main():  
    
    
    # command line parameters are handled here
    # Note: Probably not a good idea to have more than one of -i, -l, -n active at once, but we don't check
    params = ResultCLP.argvHandler(sys.argv)
    rv = [] # list of result objects
    kp = []
    em = []
    emLast = []
    pList = []
    
    # if readFiles is true, we read a list of files from stdin
    if (params["readFiles"] is True):
        f = sys.stdin.readlines()
        for fn in f:
            # strip trailing whitespace and \n
            fn = fn.rstrip()
            rv.append(Results())
            rv[-1].parseResults(filename=fn)
            
    # if not, we read one set of results directly from stdin
    else:
        rv.append(Results())
        rv[-1].parseResults() # no parameters -- read from stdin    
    
    for r in rv:
        
        #r.parseResults() # no parameters -- read from stdin
    
        if (params["confidence"] is True):
    
            if (params["kpErrors"] is True):
                kp.append(r.kpConfErrors())
    
            if (params["emErrors"] is True):
                em.append(list(r.emConfErrors()))
            
            if (params["emLastErrors"] is True):
                emLast.append(r.emConfErrors()[-1])
    
        else:
            
            if (params["kpErrors"] is True):
                kp.append(r.kpErrors())
                
            if (params["emErrors"] is True):
                em.append(list(r.emErrors()))
            
            if (params["emLastErrors"] is True):
                emLast.append(r.emErrors()[-1])
                
        if (params["displayP"] is True):
            pList.append(r.P)
                
    if (params["displayP"] is True):
        for P in pList:
            printP(P)

    if (params["kpErrors"] is True):
        list2csv(kp)
            
    if (params["emErrors"] is True):
        listOfLists2csv(em) 
        
    if (params["emLastErrors"] is True):
        list2csv(emLast)

    

if __name__ == "__main__":
    main()