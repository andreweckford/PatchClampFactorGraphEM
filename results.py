#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:03:39 2020

@author: andreweckford
"""

from resultHandling.resultHandler import Results
from resultHandling.resultArgvHandler import ResultCLP
import sys


def main():  
    params = ResultCLP.argvHandler(sys.argv)

    r = Results()
    r.parseResults() # no parameters -- read from stdin


    if (params["confidence"] is True):

        if (params["kpErrors"] is True):
            print(r.kpConfErrors())

        if (params["emErrors"] is True):
            print(r.emConfErrors())
        
        if (params["emLastErrors"] is True):
            print(r.emConfErrors()[-1])

    else:
        
        if (params["kpErrors"] is True):
            print(r.kpErrors())
            
        if (params["emErrors"] is True):
            print(r.emErrors())
        
        if (params["emLastErrors"] is True):
            print(r.emErrors()[-1])
            
    if (params["displayP"] is True):
        print(r.P)
    
    

if __name__ == "__main__":
    main()