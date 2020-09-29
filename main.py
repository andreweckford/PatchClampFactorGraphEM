#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:27:26 2020

@author: andreweckford
"""

import numpy as np
import SumProductNodes as sp
from myTools import getSteadyStateDist
#from CFTR import Receptor
from receptor.ACh import Receptor
from Simulator import Simulator
#import matplotlib.pyplot as plt
import sys

# The default value below allows you to simply run this script.
# If you want to change the value, better to execute it from the console,
# for example if you want a length of 1000:
# 
# from main import main
# main(1000)

def main(numTimeInstants = 1000, 
         maxEMIterations = 10, 
         confidence = 0.8, 
         C1aXP=0.1, 
         lastOnly=False, 
         inputData=None,
         printP=False):    

    if (len(sys.argv) == 5):
        numTimeInstants = int(sys.argv[1])
        maxEMIterations = int(sys.argv[2])
        confidence = float(sys.argv[3])
        C1aXP = float(sys.argv[4])
        
    if (len(sys.argv) == 6):
        numTimeInstants = int(sys.argv[1])
        maxEMIterations = int(sys.argv[2])
        confidence = float(sys.argv[3])
        C1aXP = float(sys.argv[4])
        if (sys.argv[5] == '-l'):
            lastOnly = True
            
    # if we have input data, i.e., inputData is not None,
    # numTimeInstants need not be specified and C1aXP is irrelevant
    # also, the simulation is irrelevant, but we will keep a dummy simulation
    # so that the data processing doesn't need to be modified
    
    if inputData is not None:
        numTimeInstants = len(inputData)
    
    # uncomment for ACh-like graph
    receptorModel = Receptor()
    
    # uncomment for CFTR-like graph
    #receptorModel = Receptor(C1aExitProb=C1aXP) # CFTR parameter object
        
    # get the parameters out from the model
    P0 = receptorModel.P0 
    P1 = receptorModel.P1
    statemap = receptorModel.statemap
    
    # Normally we can have a varying input, here with binary inputs "0" and "1" ...
    # ... for this example we only use input "0"
    # ... so the state transition matrix is given entirely by P0
    # ... and P1 is irrelevant
    px = np.array([1.,0.])
    
    # simulator object
    sim = Simulator(receptorModel,px)
    
    # following the above, the inputs are all "0"
    inputs = np.zeros(numTimeInstants)
    
    # sequence of states and channel openings from the simulator
    if inputData is not None:
        ionChannels = inputData
        states = np.zeros(len(inputData))
    else:
        states,ionChannels = sim.getReceptorState(inputs)
    
    list2csv(states)
    
    # number of states
    numStates = len(statemap)
    #numTimeInstants = len(ionChannels)
    
    P = [P0,P1]
    
    # here's where the sum-product algoirthm is implemented        

    #np.set_printoptions(precision=3,suppress=True)
    
    # initial estimate of P

    # this gives the initial estimate a diagonal bias
    dt = 1.
    P = (1-dt)*np.eye(numStates) + dt*receptorModel.Pmask
    #P = P0
    
    # force the initial estimate to be a probability, if it's not already
    for i in range(0,numStates):
        P[i,:] = P[i,:] / np.sum(P[i,:])
    
    # get the state estimates when the system parameters are perfectly known
    vvv = eStep(P0,ionChannels,statemap)
    list2csv(stateGuesses(vvv))
    list2csv(confidentStateGuesses(vvv,confidence))
    
    for emIter in range(0,maxEMIterations):
            
        # fix later ... we changed P from [P0,P1] to just P ...
        # now we are just ignoring px
        initRightMessage = getSteadyStateDist(P)
        initLeftMessage = np.ones(numStates)
        
        # E-step

        # create list of required nodes
        v = []
        s = []
        c = []
        
        for i in range(0,numTimeInstants):
            # v[i] refers to state variable s_i
            v.append(sp.StateNode())
            c.append(sp.IonChannelNode(ionChannels[i],statemap))
    
        for i in range(0,numTimeInstants-1):
            # s[i] refers to factor p(s_{i+1} | s_i)
            # note that there is one less factor node than variable node
            # because the factor nodes appear between the variables
            s.append(sp.MarkovFactorNode(P))
            # NOTE: For now it doesn't matter because we effectively use P0 only,
            # but we will have to check whether inputs[i] is correctly aligned
            # with this factor node (my guess: it's not, should be inputs[i+1])

        
        # rightward messages
        # initial, for the first variable and state
        v[0].setChannelMessage(c[0].message())
        v[0].setRightInMessage(initRightMessage)        
        s[0].setRightInMessage(v[0].rightOutMessage())
        
        for i in range(1,numTimeInstants-1):
            v[i].setChannelMessage(c[i].message())
            v[i].setRightInMessage(s[i-1].rightOutMessage())
            s[i].setRightInMessage(v[i].rightOutMessage())
            #print(v[i].rightOutMessage())
            
        # finally ... the last rightward messages
        v[numTimeInstants-1].setChannelMessage(c[numTimeInstants-1].message())
        v[numTimeInstants-1].setRightInMessage(s[numTimeInstants-2].rightOutMessage())
        
        
        # leftward messages
        # initial, from the left
        # no need to set the channel messages
        v[numTimeInstants-1].setLeftInMessage(initLeftMessage)
        for i in range(numTimeInstants-2,-1,-1):
            s[i].setLeftInMessage(v[i+1].leftOutMessage())
            v[i].setLeftInMessage(s[i].leftOutMessage())
            #print(v[i].leftOutMessage())
            
        # M-step
            
        # numTimeInstants-1 because that is the number of initialized factors
        Q = np.zeros((numStates,numStates))
        for i in range(1,numTimeInstants-1):
            #print(Q)
            #print(s[i])
            #print(s[i].aPosteriori())
            Q += s[i].aPosteriori()
            
        for i in range(0,numStates):
            Q[i,:] = Q[i,:] / np.sum(Q[i,:])
        
        P = Q    
                
        if ((lastOnly is False) or (emIter == maxEMIterations-1)):
            list2csv(stateGuesses(v))
            list2csv(confidentStateGuesses(v,confidence))
            
    if printP is True:
        print(P)

def eStep(P,ionChannels,statemap):
    
    numStates = len(statemap)
    numTimeInstants = len(ionChannels)
    # fix later ... we changed P from [P0,P1] to just P ...
    # now we are just ignoring px
    initRightMessage = getSteadyStateDist(P)
    initLeftMessage = np.ones(numStates)
    
    # E-step

    # create list of required nodes
    v = []
    s = []
    c = []
    
    for i in range(0,numTimeInstants):
        # v[i] refers to state variable s_i
        v.append(sp.StateNode())
        c.append(sp.IonChannelNode(ionChannels[i],statemap))

    for i in range(0,numTimeInstants-1):
        # s[i] refers to factor p(s_{i+1} | s_i)
        # note that there is one less factor node than variable node
        # because the factor nodes appear between the variables
        s.append(sp.MarkovFactorNode(P))
        # NOTE: For now it doesn't matter because we effectively use P0 only,
        # but we will have to check whether inputs[i] is correctly aligned
        # with this factor node (my guess: it's not, should be inputs[i+1])

    
    # rightward messages
    # initial, for the first variable and state
    v[0].setChannelMessage(c[0].message())
    v[0].setRightInMessage(initRightMessage)        
    s[0].setRightInMessage(v[0].rightOutMessage())
    
    for i in range(1,numTimeInstants-1):
        v[i].setChannelMessage(c[i].message())
        v[i].setRightInMessage(s[i-1].rightOutMessage())
        s[i].setRightInMessage(v[i].rightOutMessage())
        #print(v[i].rightOutMessage())
        
    # finally ... the last rightward messages
    v[numTimeInstants-1].setChannelMessage(c[numTimeInstants-1].message())
    v[numTimeInstants-1].setRightInMessage(s[numTimeInstants-2].rightOutMessage())
    
    
    # leftward messages
    # initial, from the left
    # no need to set the channel messages
    v[numTimeInstants-1].setLeftInMessage(initLeftMessage)
    for i in range(numTimeInstants-2,-1,-1):
        s[i].setLeftInMessage(v[i+1].leftOutMessage())
        v[i].setLeftInMessage(s[i].leftOutMessage())
        #print(v[i].leftOutMessage())
        
    return v
    

def stateGuesses(v):
    numTimeInstants = len(v)
    
    for i in range(0,numTimeInstants):
        if i == 0:
            aPosterioriProbs = np.array([v[i].aPosteriori()])
        else:
            aPosterioriProbs = np.vstack((aPosterioriProbs,v[i].aPosteriori()))
          
    # highest probability state
    hpState = np.zeros(numTimeInstants)
    for i in range(0,numTimeInstants):
        hpState[i] = np.argmax(aPosterioriProbs[i,:])
        
    return hpState

def confidentStateGuesses(v,conf):
    numTimeInstants = len(v)
    
    for i in range(0,numTimeInstants):
        if i == 0:
            aPosterioriProbs = np.array([v[i].aPosteriori()])
        else:
            aPosterioriProbs = np.vstack((aPosterioriProbs,v[i].aPosteriori()))
          
    # highest probability state
    hpState = np.zeros(numTimeInstants)
    for i in range(0,numTimeInstants):
        hpState[i] = np.argmax(aPosterioriProbs[i,:])
        if aPosterioriProbs[i,int(hpState[i])] < conf:
            hpState[i] = -1
        
    return hpState    

def errorCount(states,hpState,countMinusOne=True):
    
    numTimeInstants = len(states)

    stateErrorCount = 0
    for i in range(0,numTimeInstants):
        if hpState[i] != states[i]:
            if (countMinusOne is True) or ((countMinusOne is False) and (hpState[i] != -1)):
                stateErrorCount += 1
        
    return stateErrorCount

def missedDetectionCount(states,hpState,target):

    numTimeInstants = len(states)

    numState = 0 # counts the number of times the target state occurs
    numMissed = 0 # counts the number of times the target state was missed
    for i in range(0,numTimeInstants):
        if states[i] == target:
            numState += 1
            if hpState[i] != target:
                numMissed += 1

    return [numMissed,numState]        


def falseAlarmCount(states,hpState,target):
    numTimeInstants = len(states)

    numGuess = 0 # counts the number of times the target state was guessed
    numFalse = 0 # counts the number of times the target state was guessed wrong
    for i in range(0,numTimeInstants):
        if hpState[i] == target:
            numGuess += 1
            if states[i] != target:
                numFalse += 1

    return [numFalse,numGuess]        

def noDecisionCount(chpState):
    r = 0
    for i in range(0,len(chpState)):
        if chpState[i] == -1:
            r += 1
            
    return r


def plots():
    pass

# attempt to calculate the actual likelihood function
# this will work with only a single matrix P
def likelihood(P,ionChannels,statemap):
    
    numTimeInstants = len(ionChannels)
    numStates = len(statemap)
    
    v = []
    c = []
    s = []
    
    # backwards messages only
    for i in range(0,numTimeInstants):
        # v[i] refers to state variable s_i
        v.append(sp.StateNode(normalize=False))
        c.append(sp.IonChannelNode(ionChannels[i],statemap))

    for i in range(0,numTimeInstants-1):
        # s[i] refers to factor p(s_{i+1} | s_i)
        # note that there is one less factor node than variable node
        # because the factor nodes appear between the variables
        s.append(sp.MarkovFactorNode(P,normalize=False))

    initRightMessage = getSteadyStateDist(P)
    initLeftMessage = np.ones(numStates)
    
    for i in range(0,numTimeInstants):
        v[i].setChannelMessage(c[i].message())

    v[numTimeInstants-1].setLeftInMessage(initLeftMessage)
    qqq = 0
    for i in range(numTimeInstants-2,-1,-1):
        l1 = v[i+1].leftOutMessage()
        qqq += np.log(np.sum(l1))
        s[i].setLeftInMessage(l1/np.sum(l1))
        l2 = s[i].leftOutMessage()
        qqq += np.log(np.sum(l2))
        v[i].setLeftInMessage(l2/np.sum(l2))    
        
    #print("L: "+str(qqq))
    #print("LL: "+str(qqq + np.log(np.sum(initRightMessage * v[0].leftOutMessage()))))
    return qqq + np.log(np.sum(initRightMessage * v[0].leftOutMessage()))

def qratio(a,b):
    if (b == 0):
        return -1.
    
    return a/b

def list2csv(l):
    for i in range(0,len(l)-1):
        print(str(l[i])+',',end='')
    print(str(l[-1]))
    
    
if (__name__ == '__main__'):
    
    #d = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0])
    #main(inputData = d,maxEMIterations=100,lastOnly=True,printP=True)
    main()
    

