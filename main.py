#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:27:26 2020

@author: andreweckford
"""


import numpy as np
import SumProduct as sp
from myTools import getSteadyStateDist
from System import Receptor, IIDChannelInputs, ThreeStateChain
import matplotlib.pyplot as plt

# Channelrhodopsin-like graph
#P0 = np.array([[0.9,0.1,0.],[0.,0.9,0.1],[0.1,0.,0.9]])
#P1 = np.array([[0.1,0.9,0.],[0.,0.9,0.1],[0.1,0.,0.9]])
#statemap = [1,0,0]

# probability of each input
px = np.array([0.999,0.001])

# ACh-like graph
c = IIDChannelInputs(px)
a = Receptor(ThreeStateChain(),px)

P0 = a.P0 
P1 = a.P1
statemap = a.statemap


# sequence of channel inputs (binary: 0,1)
#inputs = c.getChannelInputs(100)
inputs = np.zeros(102)

# sequence of channel openings
#r,openings = a.getReceptorState(inputs)

openings = np.zeros(102)
openings[0] = 1.
openings[-1] = 1.

# number of states
numStates = len(statemap)
numTimeInstants = len(openings)

P = [P0,P1]

intrvl = [10,20,30,40,50,60,70,80,90,100]

for z in range(0,len(intrvl)):
    
    inputs = np.zeros(intrvl[z]+2)
    openings = np.zeros(intrvl[z]+2)
    openings[0] = 1.
    openings[-1] = 1.
    
    numTimeInstants = len(openings)

    # create list of required nodes
    v = []
    s = []
    c = []
    
    for i in range(0,numTimeInstants):
        v.append(sp.StateNode())
        s.append(sp.MarkovFactorNode(P[int(inputs[i])]))
        c.append(sp.IonChannelNode(openings[i],statemap))
     
    initRightMessage = getSteadyStateDist(px[0]*P[0]+px[1]*P[1])
    initLeftMessage = np.ones(numStates)
    
    # rightward messages
    v[0].setChannelMessage(c[0].message())
    v[0].setRightInMessage(initRightMessage)
    #print(v[0].rightOutMessage())
    
    for i in range(1,numTimeInstants):
        s[i].setRightInMessage(v[i-1].rightOutMessage())
        v[i].setRightInMessage(s[i].rightOutMessage())
        v[i].setChannelMessage(c[i].message())
        #print(v[i].rightOutMessage())
    
    
    v[-1].setLeftInMessage(initLeftMessage)
    for i in range(numTimeInstants-2,-1,-1):
        s[i].setLeftInMessage(v[i+1].leftOutMessage())
        v[i].setLeftInMessage(s[i].leftOutMessage())
        #print(v[i].leftOutMessage())
        
        
    means = np.array([])
    states = np.array([1.,2.,3.])
    for i in range(1,numTimeInstants-1):
        #print(np.sum(v[i].aPosteriori()*states))
        means = np.append(means,np.sum(v[i].aPosteriori()*states))
        
    #for i in range(0,len(means)-1):
    #    print(str(means[i])+',',end='')
    #print(means[-1])
        
    
    plt.figure(1)
    plt.plot(np.arange(-intrvl[z]/2+1,intrvl[z]/2+1),means)
    plt.xlabel('Time from midpoint of Y=0 interval')
    plt.ylabel('Expected value of state number')
    plt.title('Y=0 Inteval Lengths from 10 to 100')
    plt.axis([-50,50,1.,2.])
    plt.grid(True)
    plt.savefig('result.pdf')

