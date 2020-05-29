#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:27:26 2020

@author: andreweckford
"""


import numpy as np
import SumProduct as sp
from myTools import getSteadyStateDist
from CFTR import Receptor
from Simulator import Simulator
import matplotlib.pyplot as plt

# Channelrhodopsin-like graph
#P0 = np.array([[0.9,0.1,0.],[0.,0.9,0.1],[0.1,0.,0.9]])
#P1 = np.array([[0.1,0.9,0.],[0.,0.9,0.1],[0.1,0.,0.9]])
#statemap = [1,0,0]

# length of the simulation, in time instants
numTimeInstants = 100

# ACh-like graph
cftrModel = Receptor() # CFTR parameter object

# get the parameters out from the model
P0 = cftrModel.P0 
P1 = cftrModel.P1
statemap = cftrModel.statemap

# Normally we can have a varying input, here with binary inputs "0" and "1" ...
# ... for this example we only use input "0"
# ... so the state transition matrix is given entirely by P0
# ... and P1 is irrelevant
px = np.array([1.,0.])

# simulator object
sim = Simulator(cftrModel,px)

# following the above, the inputs are all "0"
inputs = np.zeros(numTimeInstants)

# sequence of states and channel openings from the simulator
states,ionChannels = sim.getReceptorState(inputs)

# number of states
numStates = len(statemap)
numTimeInstants = len(ionChannels)

P = [P0,P1]

# here's where the sum-product algoirthm is implemented

# create list of required nodes
v = []
s = []
c = []

for i in range(0,numTimeInstants):
    v.append(sp.StateNode())
    s.append(sp.MarkovFactorNode(P[int(inputs[i])]))
    c.append(sp.IonChannelNode(ionChannels[i],statemap))
 
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
    
for i in range(0,numTimeInstants):
    if i == 0:
        aPosterioriProbs = np.array([v[i].aPosteriori()])
    else:
        aPosterioriProbs = np.vstack((aPosterioriProbs,v[i].aPosteriori()))
      
plt.figure(1)
plt.plot(ionChannels)
for i in range(0,4):
    plt.plot(aPosterioriProbs[:,i])

    

