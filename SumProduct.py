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

class IonChannelNode:
    
    # o = 1 if open, 0 if closed
    # statemap is the map of kinetic states to ion channel states
    # for example, 5 states, [0,1,0,0,1]: states 0,2,3 are closed, 1,4 are open
    def __init__(self,o,statemap):
        self.o = o
        self.statemap = statemap
        
    def message(self):
        if self.o == 0:
            return 1-np.array(self.statemap)
        
        return np.array(self.statemap)
    
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
    
class MarkovFactorNode:
    
    def __init__(self,P):
        self.rightInMessage = None
        self.leftInMessage = None
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