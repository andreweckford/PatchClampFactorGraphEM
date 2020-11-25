#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:27:26 2020

@author: andreweckford
"""

import numpy as np
import factor.SumProductNodes as sp
from factor.myTools import getSteadyStateDist
from factor.mainArgvHandler import clp
from factor.Simulator import Simulator
import sys
from factor.list2csv import list2csv,printP

def main(inputData = None):    

    # default parameters are in clp    
    
    params = clp.argvHandler(sys.argv)
    
    # exit if the parameters are invalid ... this is not fully implemented yet
    if (params["validArgv"] is False):
        print(params)
        sys.exit()
    
        
    # check to see if we print the parameters, if not then print them
    if (params["suppressParameters"] is False):
        clp.csvParams(params)    
    
    # if we have input data, i.e., inputData is not None,
    # numTimeInstants need not be specified and C1aXP is irrelevant
    # also, the simulation is irrelevant, but we will keep a dummy simulation
    # so that the data processing doesn't need to be modified
    
    if inputData is not None:
        params["numTimeInstants"] = len(inputData)
    
    # obtain the correct receptor model
    receptorModel = clp.createReceptor(params["receptor"],params["receptorParameter"])
        
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
    inputs = np.zeros(params["numTimeInstants"])
    
    # sequence of states and channel openings from the simulator
    if inputData is not None:
        ionChannels = inputData
        states = np.zeros(len(inputData))
    else:
        if (params["noiseVariance"] == 0.):
            # no noise
            states,ionChannels = sim.getReceptorState(inputs)
        else:
            # if there is nonzero noise, the observations are still in ionChannels,
            # and the true ion channel state is in trueIonChannels
            states,trueIonChannels,ionChannels = sim.getReceptorStateNoisy(inputs,params["noiseVariance"])

    # don't print the states if we are suppressing the estimates    
    if (params["suppressEstimates"] is False):
        list2csv(states)
    
    # number of states
    numStates = len(statemap)
    #numTimeInstants = len(ionChannels)
    
    P = [P0,P1]
    
    # here's where the sum-product algoirthm is implemented        

    
    # initial estimate of P

    # this gives the initial estimate a diagonal bias
    dt = 1. - params["diagonalBias"]
    P = (1-dt)*np.eye(numStates) + dt*receptorModel.Pmask
    #P = P0
    #P = (1-dt)*np.eye(numStates) + dt*np.ones((numStates,numStates))/numStates
    
    # force the initial estimate to be a probability, if it's not already
    for i in range(0,numStates):
        P[i,:] = P[i,:] / np.sum(P[i,:])
        
    # initial estimate of sigma2
    # equal to the true variance if we are not estimating
    if params["estimateNoiseVariance"] is True:
        sigma2 = params["initialNoiseVarianceEstimate"]
    else:
        sigma2 = params["noiseVariance"]
        
    
    # get the state estimates when the system parameters are perfectly known
    # don't print these estimates if we are suppressing the estimates    
    vvv = eStep(P0,ionChannels,statemap,params["noiseVariance"])
    if (params["suppressEstimates"] is False):
        list2csv(stateGuesses(vvv))
        list2csv(confidentStateGuesses(vvv,params["confidence"]))
    
    for emIter in range(0,params["maxEMIterations"]):
            
        # fix later ... we changed P from [P0,P1] to just P ...
        # now we are just ignoring px
        initRightMessage = getSteadyStateDist(P)
        initLeftMessage = np.ones(numStates)
        
        # E-step

        # create list of required nodes
        v = []
        s = []
        c = []
        
        for i in range(0,params["numTimeInstants"]):
            # v[i] refers to state variable s_i
            v.append(sp.StateNode())
            if (params["noiseVariance"] == 0.):
                # no noise case
                c.append(sp.IonChannelNode(ionChannels[i],statemap))
            else:
                c.append(sp.IonChannelNodeNoisy(ionChannels[i],statemap,sigma2))
    
        for i in range(0,params["numTimeInstants"]-1):
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
        
        for i in range(1,params["numTimeInstants"]-1):
            v[i].setChannelMessage(c[i].message())
            v[i].setRightInMessage(s[i-1].rightOutMessage())
            s[i].setRightInMessage(v[i].rightOutMessage())
            #print(v[i].rightOutMessage())
            
        # finally ... the last rightward messages
        v[params["numTimeInstants"]-1].setChannelMessage(c[params["numTimeInstants"]-1].message())
        v[params["numTimeInstants"]-1].setRightInMessage(s[params["numTimeInstants"]-2].rightOutMessage())
        
        
        # leftward messages
        # initial, from the left
        # no need to set the channel messages
        v[params["numTimeInstants"]-1].setLeftInMessage(initLeftMessage)
        for i in range(params["numTimeInstants"]-2,-1,-1):
            s[i].setLeftInMessage(v[i+1].leftOutMessage())
            v[i].setLeftInMessage(s[i].leftOutMessage())
            #print(v[i].leftOutMessage())
            
        # M-step
        
        # transition probability matrix estimate
        # numTimeInstants-1 because that is the number of initialized factors
        Q = np.zeros((numStates,numStates))
        for i in range(1,params["numTimeInstants"]-1):
            #print(Q)
            #print(s[i])
            #print(s[i].aPosteriori())
            Q += s[i].aPosteriori()
            
        for i in range(0,numStates):
            Q[i,:] = Q[i,:] / np.sum(Q[i,:])

        P = Q    
        
        # noise estimate
        # we don't estimate the noise if the noise variance is zero, obviously
        if (params["estimateNoiseVariance"] is True) and (params["noiseVariance"] > 0):
            Q_sigma2 = estimateNoiseVariance(c,v,ionChannels)
            #print(Q_sigma2)
            sigma2 = Q_sigma2
        
        # print the estimates
        # don't print if (a) we're not on the last estimate and lastOnly is true, or
        # (b) suppressEstimates is True
        if ((params["lastOnly"] is False) or (emIter == params["maxEMIterations"]-1)) and (params["suppressEstimates"] is False):
            list2csv(stateGuesses(v))
            list2csv(confidentStateGuesses(v,params["confidence"]))

    # after the estimates, print the estimated matrix if the flag is set
    if params["suppressP"] is False:
        printP(P)

# do the M step estimation of the noise variance
def estimateNoiseVariance(c,v,ionChannels):
    z2Sum = 0 # sum of squares of channel observations
    yzSum = 0 # sum of channel obs * E[y_i]
    y2Sum = 0 # sum of E[y_i^2]
    
    # calculate the sums
    for i in range(0,len(c)):
        z2Sum += np.power(ionChannels[i],2)
        
        # posterior distribution of ion channel state y given everything known
        foo = c[i].emPosteriori(v[i].aPosterioriNoChannel())
        
        # since y \in {0,1}, E[y_i] = E[y_i^2] = foo[1]
        yzSum += foo[1] * ionChannels[i]
        y2Sum += np.power(foo[1],2)
    
    q = (z2Sum - 2*yzSum + y2Sum)/len(c)
    return q
    

def eStep(P,ionChannels,statemap,sigma2=None):
    
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
        if (sigma2 is None) or (sigma2 == 0.):
            # no noise case
            c.append(sp.IonChannelNode(ionChannels[i],statemap))
        else:
            # case with noise
            c.append(sp.IonChannelNodeNoisy(ionChannels[i],statemap,sigma2))

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


    
    
if (__name__ == '__main__'):
    
    #d = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0])
    #main(inputData = d,maxEMIterations=100,lastOnly=True,printP=True)
    main()
    

