#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 10:40:38 2020

@author: andreweckford
"""

from receptors import CFTR,ACh

class clp:
    
    # entries of the form: dictionary name, flag name, bare flag?, type
    # bare flags do not take any additional values, type is ignored
    # e.g. -l is a bare flag, -n=10000 is not

    flags = [
        ['numTimeInstants','-n',False,'int'],
        ['maxEMIterations','-i',False,'int'],
        ['confidence','-c',False,'float'],
        ['receptorParameter','-p',False,'string'],
        ['lastOnly','-l',True,''],
        ['receptor','-r',False,'string']
    ]


    # this is where the initial parameters are defined
    # dictionary keys are given above in flags
    initParams = {
        flags[0][0] : 1000,
        flags[1][0] : 10,
        flags[2][0] : 0.8,
        flags[3][0] : "0.1",
        flags[4][0] : False,
        flags[5][0] : "CFTR",
        "printP" : False,
        "validArgv" : True
    }
    
    # handle command line parameters
    def argvHandler(argv):
        
        params = dict(clp.initParams) # copy constructor
        
        # argv[0] is not useful
        for i in range(1,len(argv)):
            # scroll through all the possible flags looking for a match
            for j in range(0,len(clp.flags)):
                # match the short or long flag name
                if argv[i].startswith('--'+clp.flags[j][0]) or argv[i].startswith(clp.flags[j][1]):
                    # with a match, check to see if it's a bare flag
                    if (clp.flags[j][2] is True):
                        # if bare, just set the parameter to True
                        params[clp.flags[j][0]] = True
                    else:                        
                        # if not bare, extract the parameter
                        foo = argv[i].split('=')
                        if (len(foo) == 2):
                            # find the correct data type
                            if (clp.flags[j][3] == 'int'):
                                params[clp.flags[j][0]] = int(foo[1])
                            elif (clp.flags[j][3] == 'float'):
                                params[clp.flags[j][0]] = float(foo[1])
                            else:
                                params[clp.flags[j][0]] = foo[1]
        
        return params

    # creates and returns a receptor object with the given parameters
    # rString and rParams are given by the command line parameters
    # each possible receptor in the receptors directory should have an entry
    def createReceptor(rString,rParams):
        
        if (rString == "CFTR"):
            # parameter for CFTR is transition probability 1->2
            r = CFTR.Receptor(float(rParams))
            return r
        
        if (rString == "ACh"):
            # parameter for ACh is concentration in the P0 state
            r = ACh.Receptor(float(rParams))
            return r
                    
        return None
            