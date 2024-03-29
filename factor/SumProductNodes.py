#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:25:42 2020

@author: andreweckford
"""

# Note: "left" and "right" messages are named for their direction of travel
# "left" messages are from the future (on the right) to the past (on the left)
# "right" messages are from the past to the future

import numpy as np
import time

class IonChannelNode:
    
    # statemap is the map of kinetic states to ion channel states
    # for example, 5 states, [0,1,0,0,1]: states 0,2,3 are closed, 1,4 are open
    # the ion channel may be nonbinary, such as two conduction states numbered 1 and 2
    def __init__(self,o,statemap):
        self.o = o
        self.statemap = statemap
        
    def message(self):
        r = np.zeros(len(self.statemap))
        for i in range(0,len(self.statemap)):
            if self.o == self.statemap[i]:
                r[i] = 1.
            else:
                r[i] = 0.
        
        return r

class IonChannelNodeNoisy:
    
    def __init__(self,o,statemap,currentEstimate,sigma2):
        self.o = o
        self.statemap = statemap
        self.sigma2 = sigma2
        self.currentEstimate = currentEstimate
        
    def message(self):
        r = np.zeros(len(self.statemap))
        for i in range(0,len(self.statemap)):
            r[i] = 1/np.sqrt(2*np.pi*self.sigma2)
            if (self.statemap[i] == 0):
                r[i] *= np.exp(-1/(2*self.sigma2) * np.power(self.o - self.currentEstimate[0],2))
            elif (self.statemap[i] == 1):
                r[i] *= np.exp(-1/(2*self.sigma2) * np.power(self.o - self.currentEstimate[1],2))
            else:
                print("Boo!")
        return r
    
    def emPosteriori(self,stateExtrinsicMessage):
        # stateExtrinsicMessage carries information from the rest of the graph
        # not including the channel information from this node
        
        # extrinsic message for y
        # first element is y = 0, second element is y = 1
        emsg_y = np.zeros(2)
        
        for i in range(0,len(self.statemap)):
            if self.statemap[i] == 0:
                emsg_y[0] += stateExtrinsicMessage[i]
            else:
                emsg_y[1] += stateExtrinsicMessage[i]
                
        # channel message for y
        #cmsg_y = np.exp(-1/(2*self.sigma2) * np.power(self.o - np.array([0.,1.]),2))
        #print(self.currentEstimate)
        cmsg_y = np.exp(-1/(2*self.sigma2) * np.power(self.o - self.currentEstimate,2))
        cmsg_y *= 1/np.sqrt(2*np.pi*self.sigma2)
        
        return (emsg_y * cmsg_y)/np.sum(emsg_y * cmsg_y)

    
class StateNode:
    
    # left and right messages are np.array vectors (row vectors)
    def __init__(self,normalize=True):
        self.channelMessage = None
        self.rightInMessage = None
        self.leftInMessage = None
        self.normalize = normalize
    
    # inbound message from channel
    def setChannelMessage(self,channelMessage):
        self.channelMessage = channelMessage
        
    # inbound message towards the right (from the past)
    def setRightInMessage(self,rightInMessage):
        self.rightInMessage = rightInMessage
    
    # inbound message towards the left (from the future)
    def setLeftInMessage(self,leftInMessage):
        self.leftInMessage = leftInMessage
    
    def leftOutMessage(self):
        if self.channelMessage is None:
            return None
        
        if self.leftInMessage is None:
            return None
        
        r = self.leftInMessage * self.channelMessage
        if self.normalize is True:
            return r / np.sum(r)
        
        return r

    def rightOutMessage(self):
        if self.channelMessage is None:
            return None
        
        if self.rightInMessage is None:
            return None
        
        r = self.rightInMessage * self.channelMessage
        if self.normalize is True:
            return r / np.sum(r)
        
        return r
    
    def aPosteriori(self):
        if self.channelMessage is None:
            return None
        
        if self.rightInMessage is None:
            return None
        
        if self.leftInMessage is None:
            return None
        
        r = self.rightInMessage * self.channelMessage * self.leftInMessage
        if self.normalize is True:
            return r / np.sum(r)
        
        return r
    
    def aPosterioriNoChannel(self):
        if self.rightInMessage is None:
            return None
        
        if self.leftInMessage is None:
            return None
        
        r = self.rightInMessage * self.leftInMessage
        if self.normalize is True:
            return r / np.sum(r)
        
        return r
        

class MarkovFactorNode:
    
    def __init__(self,P,normalize=True):
        self.rightInMessage = None
        self.leftInMessage = None
        self.P = P
        self.normalize = normalize
    
    def updateP(self,P):
        self.P = P
    
    # set inbound message from the right
    def setRightInMessage(self,rightInMessage):
        self.rightInMessage = rightInMessage
        
    # set inbound message from the left
    def setLeftInMessage(self,leftInMessage):
        self.leftInMessage = leftInMessage
        
    def leftOutMessage(self):
        if self.leftInMessage is None:
            return None
        
        return self.P @ self.leftInMessage.T
    
    def rightOutMessage(self):
        if self.rightInMessage is None:
            return None
        
        return self.rightInMessage @ self.P        
    
    def aPosteriori(self):
        if self.rightInMessage is None:
            return None
        
        if self.leftInMessage is None:
            return None
        
        r = np.diag(self.rightInMessage) @ self.P @ np.diag(self.leftInMessage)
        if self.normalize is False:
            return r
        
        return r/np.sum(r)
        
