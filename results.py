#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:03:39 2020

@author: andreweckford
"""

from resultHandling.resultHandler import Results
from resultHandling.resultArgvHandler import ResultCLP
from factor.list2csv import list2csv,printP,listOfLists2csv
import numpy as np
import sys
from bpr import bpr_md_new,bpr_fa_new,openToCloseIntervals,countClosings


def main():  
    
    
    # command line parameters are handled here
    # Note: Probably not a good idea to have more than one of -i, -l, -n active at once, but we don't check
    params = ResultCLP.argvHandler(sys.argv)
    rv = [] # list of result objects
    kp = []
    em = []
    cc = []
    emLast = []
    pList = []
    
    # if readFiles is true, we read a list of files from stdin
    if (params["readFiles"] is True):
        f = sys.stdin.readlines()
        for fn in f:
            # strip trailing whitespace and \n
            fn = fn.rstrip()
            rv.append(Results())
            rv[-1].parseResults(filename=fn,rawData=params["rawData"])
            
    # if not, we read one set of results directly from stdin
    else:
        rv.append(Results())
        rv[-1].parseResults(rawData=params["rawData"]) # no parameters -- read from stdin
    
    for r in rv:
        
        # note, in these lines we cast to np.array to make the output similar to the confidence outputs
        # because we use similar scripts handle them when we generate plots
        # and we didn't notice beforehand that the confidence outputs are np.array
        # ... the code is acquiring a lot of special cases, one of these days we should refactor
    
        if (params["transitionInitial"] is not None) and (params["transitionFinal"] is not None):
            if (params["confidence"] is True):
                foo = r.transition(r.emEstimates_conf[-1],params["transitionInitial"],params["transitionFinal"])
            else:
                foo = r.transition(r.emEstimates[-1],params["transitionInitial"],params["transitionFinal"])
            if len(foo) > 0:
                list2csv(foo)
            sys.exit()
        
        elif (params["openIntervals"] is True):
            list2csv(openToCloseIntervals(r,openings=True))
            sys.exit()
            
        elif (params["closedIntervals"] is True):
            list2csv(openToCloseIntervals(r,openings=False))
            sys.exit()

        elif (params["countClosings"] is True):
            cc.append(np.array(countClosings(r.states)))
        
        elif (params["permissiveMD"] is True):
            
            if (params["kpErrors"] is True):
                kp.append(np.array(bpr_md_new(r,conf=params["confidence"],getEmResults=False)))

            if (params["emErrors"] is True):
                print("Invalid")

            if (params["emLastErrors"] is True):
                emLast.append(np.array(bpr_md_new(r,conf=params["confidence"],getEmResults=True)))
        
        elif (params["permissiveFA"] is True):

            if (params["kpErrors"] is True):
                kp.append(np.array(bpr_fa_new(r,conf=params["confidence"],getEmResults=False)))

            if (params["emErrors"] is True):
                print("Invalid")

            if (params["emLastErrors"] is True):
                emLast.append(np.array(bpr_fa_new(r,conf=params["confidence"],getEmResults=True)))

        elif (params["confidence"] is True):
    
            if (params["kpErrors"] is True):
                kp.append(r.kpConfErrors())
    
            if (params["emErrors"] is True):
                em.append(list(r.emConfErrors()))
            
            if (params["emLastErrors"] is True):
                emLast.append(r.emConfErrors()[-1])
    
        else:
            
            if (params["errorsForTrue"] is None):
            
                if (params["kpErrors"] is True):
                    kp.append(r.kpErrors())
                    
                if (params["emErrors"] is True):
                    em.append(list(r.emErrors()))
                
                if (params["emLastErrors"] is True):
                    emLast.append(r.emErrors()[-1])
                    
            else:
                
                if (params["kpErrors"] is True):
                    kp.append(r.kpStateErrors(params['errorsForTrue']))
                    
                if (params["emErrors"] is True):
                    print("Invalid")
                
                if (params["emLastErrors"] is True):
                    emLast.append(r.emStateErrors(params['errorsForTrue'])[-1])                
                
        if (params["displayP"] is True):
            pList.append(r.P)

    for r in rv:
        
        if (params["suppressParams"] is False):
            #list2csv(list(r.params.keys()))
            list2csv(list(r.params.values()))
                
    if (params["displayP"] is True):
        for P in pList:
            printP(P)

    if (params["countClosings"]) is True:
        list2csv(cc)

    if (params["kpErrors"] is True):
        if (params["permissiveMD"] is True) or (params["permissiveFA"] is True):
            list2csv(kp)
        else:
            list2csv(kp)
            
    if (params["emErrors"] is True):
        listOfLists2csv(em) 
        
    if (params["emLastErrors"] is True):
        if (params["permissiveMD"] is True) or (params["permissiveFA"] is True):
            list2csv(emLast)
        else:
            list2csv(emLast)

if __name__ == "__main__":
    main()
