#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 17:06:47 2020

@author: andreweckford
"""

import sys

class ResultCLP:
    
    # entries of the form: dictionary name, flag name, bare flag?, type
    # bare flags do not take any additional values, type is ignored
    # e.g. -l is a bare flag, -n=10000 is not

    # no short flag should be the prefix of any other short flag!
    flags = [
        ['kpErrors','-n',True,''],
        ['emErrors','-i',True,''],
        ['emLastErrors','-l',True,''],
        ['displayP','-p',True,''],
        ['confidence','-c',True,''],
        ['readFiles','-f',True,''],
        ['suppressParams','-s',True,''],
        ['errorsForTrue','-t',False,'int'],
        ['permissiveMD','-m',True,''],
        ['permissiveFA','-a',True,''],
        ['openIntervals','-ov',True,''],
        ['closedIntervals','-zv',True,''],
        ['transitionInitial','-qi',False,'state'],
        ['transitionFinal','-qf',False,'state']
    ]

    # this is where the initial parameters are defined
    # dictionary keys are given above in flags
    initParams = {
        flags[0][0] : False,
        flags[1][0] : False,
        flags[2][0] : False,
        flags[3][0] : False,
        flags[4][0] : False,
        flags[5][0] : False,
        flags[6][0] : False,
        flags[7][0] : None,
        flags[8][0] : False,
        flags[9][0] : False,
        flags[10][0] : False,
        flags[11][0] : False,
        flags[12][0] : None,
        flags[13][0] : None
    }
        
    # handle command line parameters
    def argvHandler(argv):
        
        # improve on help later
        if (argv[1] == '--help') or (argv[1] == '-h'):
            print(ResultCLP.flags)
            sys.exit()
        
        params = dict(ResultCLP.initParams) # copy constructor
        
        # argv[0] is not useful
        for i in range(1,len(argv)):
            # scroll through all the possible flags looking for a match
            for j in range(0,len(ResultCLP.flags)):
                # match the short or long flag name
                if argv[i].startswith('--'+ResultCLP.flags[j][0]) or argv[i].startswith(ResultCLP.flags[j][1]):
                    # with a match, check to see if it's a bare flag
                    if (ResultCLP.flags[j][2] is True):
                        # if bare, just set the parameter to True
                        params[ResultCLP.flags[j][0]] = True
                    else:                        
                        # if not bare, extract the parameter
                        foo = argv[i].split('=')
                        if (len(foo) == 2):
                            # find the correct data type - some of these are unused in this package
                            if (ResultCLP.flags[j][3] == 'int'):
                                params[ResultCLP.flags[j][0]] = int(foo[1])
                            elif (ResultCLP.flags[j][3] == 'float'):
                                params[ResultCLP.flags[j][0]] = float(foo[1])
                            elif (ResultCLP.flags[j][3] == 'state'):
                                # open and close methods are expecting states of the form '2.0' for state 2
                                params[ResultCLP.flags[j][0]] = foo[1] + '.0'
                            else:
                                params[ResultCLP.flags[j][0]] = foo[1]
        
        return params
