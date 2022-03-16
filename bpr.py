import numpy as np
from resultHandling.resultHandler import Results

# What do we want to do here?
#
# Given a transition, e.g. 2-3 close-open:
# - number of transitions in the ground truth (looking for missed detections)
# -- number of those correctly detected
# -- number of those mistaken for the other kind of transition (here 5-4)
# -- number of those mistaken for close-close, open-close, or open-open
# -- same as above but at high confidence
# - number of transitions claimed in the estimated sequence (looking for false alarms)
# -- number of those correctly detected
# -- number of those mistaken for the other kind of transition (here 5-4)
# -- number of those mistaken for close-close, open-close, or open-open
# -- same as above but at high confidence

# r is a results object
def bpr_md(r,doOpenStates=True,conf=False,getEmResults=True):


    # count number of transitions in the ground truth sequence
    #for i in range(0,len(r.kp)):
    #    print(str(i) + ' ' + r.kp[i])
    
    openTransition_v = [['2.0','3.0'],['5.0','4.0']]
    closeTransition_v = [['3.0','2.0'],['4.0','5.0']]
    openStates = ['3.0','4.0']
    closedStates = ['0.0','1.0','2.0','5.0','6.0']
    open_gt = []
    open_kp = []
    open_em = []

    states = r.states

    if conf is True:
        kp = r.kp_conf
        em = r.emEstimates_conf[-1]
    else:
        kp = r.kp
        em = r.emEstimates[-1]

    if doOpenStates is True:
        iterate_v = openTransition_v
    else:
        iterate_v = closeTransition_v
        
    for ot in iterate_v:
        # transitions in the ground truth sequence
        open_gt.append(r.transition(states,ot[0],ot[1]))
        
        # transitions in the estimated sequence with known probabilities
        open_kp.append(r.transition(kp,ot[0],ot[1]))
        
        # transitions in the EM estimated sequence
        open_em.append(r.transition(em,ot[0],ot[1]))
        
    #print(open_gt)
    #print(open_kp)
    #print(open_em)
    
    # count the number correctly detected
    gtLen = [len(open_gt[0]),len(open_gt[1])]
    kpCorrectlyDetected = [0]*len(open_gt)
    emCorrectlyDetected = [0]*len(open_gt)
    
    for i in range(0,len(open_gt)):
        for j in open_gt[i]:
            if j in open_kp[i]:
                kpCorrectlyDetected[i] += 1
            if j in open_em[i]:
                emCorrectlyDetected[i] += 1
    
    #print(gtLen)
    #print(kpCorrectlyDetected)
    #print(emCorrectlyDetected)

    # count the number mistaken for hte other transition
    kpOtherTransition = [0]*len(open_gt)
    emOtherTransition = [0]*len(open_gt)
    
    # did we mistake one type of opening for another?
    for i in range(0,len(open_gt)):
        for j in open_gt[i]:
            for k in range(0,len(open_gt)):
                if (i != k):
                    if j in open_kp[k]:
                        kpOtherTransition[i] += 1
                    if j in open_em[k]:
                        emOtherTransition[i] += 1
                     
    #print(kpOtherTransition)
    #print(emOtherTransition)
    
    # count the number where the ground truth has a transition, but the estimate stays open or closed
    kpStaysOpen = [0]*len(open_gt)
    emStaysOpen = [0]*len(open_gt)
    kpStaysClosed = [0]*len(open_gt)
    emStaysClosed = [0]*len(open_gt)
   
    for i in range(0,len(open_gt)):
        for j in open_gt[i]:
            if (kp[j-1] in openStates) and (kp[j] in openStates):
                kpStaysOpen[i] += 1
            if (kp[j-1] in closedStates) and (kp[j] in closedStates):
                kpStaysClosed[i] += 1

        for j in open_gt[i]:
            if (em[j-1] in openStates) and (em[j] in openStates):
                emStaysOpen[i] += 1
            if (em[j-1] in closedStates) and (em[j] in closedStates):
                emStaysClosed[i] += 1

    #print(kpStaysOpen)
    #print(emStaysOpen)
    #print(kpStaysClosed)
    #print(emStaysClosed)
    
    # count the number where the ground truth has a transition, but the estimate has the opposite transition
    kpWrongWay = [0]*len(open_gt)
    emWrongWay = [0]*len(open_gt)
    
    # did we mistake one type of opening for another?
    for i in range(0,len(open_gt)):
    
        if doOpenStates is True:
            for j in open_gt[i]:
                if (kp[j] in closedStates) and (kp[j-1] in openStates):
                    kpWrongWay += 1
            for j in open_gt[i]:
                if (em[j] in closedStates) and (em[j-1] in openStates):
                    emWrongWay += 1
        else:
            for j in open_gt[i]:
                if (kp[j] in openStates) and (kp[j-1] in closedStates):
                    kpWrongWay += 1
            for j in open_gt[i]:
                if (em[j] in openStates) and (em[j-1] in closedStates):
                    emWrongWay += 1


                     
    #print(kpWrongWay)
    #print(emWrongWay)
    
    # assemble results
    # total, correct, other transition, stays open, stays closed, wrong way
    kpResult = []
    emResult = []
    for i in range(0,len(open_gt)):
        kpResult.append([len(open_gt[i]),
                         kpCorrectlyDetected[i],
                         kpOtherTransition[i],
                         kpStaysOpen[i],
                         kpStaysClosed[i],
                         kpWrongWay[i]])
        emResult.append([len(open_gt[i]),
                         emCorrectlyDetected[i],
                         emOtherTransition[i],
                         emStaysOpen[i],
                         emStaysClosed[i],
                         emWrongWay[i]])

    if getEmResults is True:
        foo = [i for s in emResult for i in s]
    else:
        foo = [i for s in kpResult for i in s]
        
    return foo

# false alarms are different enough to warrant a new method
def bpr_fa(r,doOpenStates=False,conf=False,getEmResults=True):


    # count number of transitions in the ground truth sequence
    #for i in range(0,len(r.kp)):
    #    print(str(i) + ' ' + r.kp[i])
    
    openTransition_v = [['2.0','3.0'],['5.0','4.0']]
    closeTransition_v = [['3.0','2.0'],['4.0','5.0']]
    openStates = ['3.0','4.0']
    closedStates = ['0.0','1.0','2.0','5.0','6.0']
    open_gt = []
    open_kp = []
    open_em = []

    states = r.states

    if conf is True:
        kp = r.kp_conf
        em = r.emEstimates_conf[-1]
    else:
        kp = r.kp
        em = r.emEstimates[-1]

    if doOpenStates is True:
        iterate_v = openTransition_v
    else:
        iterate_v = closeTransition_v
        
    for ot in iterate_v:
        # transitions in the ground truth sequence
        open_gt.append(r.transition(states,ot[0],ot[1]))
        
        # transitions in the estimated sequence with known probabilities
        open_kp.append(r.transition(kp,ot[0],ot[1]))
        
        # transitions in the EM estimated sequence
        open_em.append(r.transition(em,ot[0],ot[1]))
        
    #print(open_gt)
    #print(open_kp)
    #print(open_em)
    
    # count the number correctly detected
    gtLen = [len(open_gt[0]),len(open_gt[1])]
    kpCorrectlyDetected = [0]*len(open_gt)
    emCorrectlyDetected = [0]*len(open_gt)
    
    for i in range(0,len(open_gt)):
        for j in open_kp[i]:
            if j in open_gt[i]:
                kpCorrectlyDetected[i] += 1
        for j in open_em[i]:
            if j in open_gt[i]:
                emCorrectlyDetected[i] += 1
    
    #print(gtLen)
    #print(kpCorrectlyDetected)
    #print(emCorrectlyDetected)

    # count the number mistaken for hte other transition
    kpOtherTransition = [0]*len(open_gt)
    emOtherTransition = [0]*len(open_gt)
    
    # did we mistake one type of opening for another?
    for i in range(0,len(open_gt)):
        for j in open_kp[i]:
            for k in range(0,len(open_gt)):
                if (i != k) and (j in open_gt[k]):
                    kpOtherTransition[i] += 1
                    
        for j in open_em[i]:
            for k in range(0,len(open_gt)):
                if (i != k) and (j in open_gt[k]):
                    emOtherTransition[i] += 1
                     
    #print(kpOtherTransition)
    #print(emOtherTransition)
    
    # count the number where the ground truth has a transition, but the estimate stays open or closed
    kpStaysOpen = [0]*len(open_gt)
    emStaysOpen = [0]*len(open_gt)
    kpStaysClosed = [0]*len(open_gt)
    emStaysClosed = [0]*len(open_gt)
   
    for i in range(0,len(open_gt)):
        for j in open_kp[i]:
            if (states[j-1] in openStates) and (states[j] in openStates):
                kpStaysOpen[i] += 1
            if (states[j-1] in closedStates) and (states[j] in closedStates):
                kpStaysClosed[i] += 1

        for j in open_em[i]:
            if (states[j-1] in openStates) and (states[j] in openStates):
                emStaysOpen[i] += 1
            if (states[j-1] in closedStates) and (states[j] in closedStates):
                emStaysClosed[i] += 1

    #print(kpStaysOpen)
    #print(emStaysOpen)
    #print(kpStaysClosed)
    #print(emStaysClosed)
    
    # count the number where the ground truth has a transition, but the estimate has the opposite transition
    kpWrongWay = [0]*len(open_gt)
    emWrongWay = [0]*len(open_gt)
    
    # did we mistake one type of opening for another?
    for i in range(0,len(open_gt)):
    
        for j in open_kp[i]:
            if doOpenStates is True:
                if (states[j] in closedStates) and (kp[j-1] in openStates):
                    kpWrongWay += 1
            else:
                if (states[j] in openStates) and (kp[j-1] in closedStates):
                    kpWrongWay += 1
                    
        for j in open_em[i]:
            if doOpenStates is True:
                if (states[j] in closedStates) and (kp[j-1] in openStates):
                    emWrongWay += 1
            else:
                if (states[j] in openStates) and (kp[j-1] in closedStates):
                    emWrongWay += 1

                     
    #print(kpWrongWay)
    #print(emWrongWay)
    
    # assemble results
    # total, correct, other transition, stays open, stays closed, wrong way
    kpResult = []
    emResult = []
    for i in range(0,len(open_gt)):
        kpResult.append([len(open_kp[i]),
                         kpCorrectlyDetected[i],
                         kpOtherTransition[i],
                         kpStaysOpen[i],
                         kpStaysClosed[i],
                         kpWrongWay[i]])
        emResult.append([len(open_em[i]),
                         emCorrectlyDetected[i],
                         emOtherTransition[i],
                         emStaysOpen[i],
                         emStaysClosed[i],
                         emWrongWay[i]])

    if getEmResults is True:
        foo = [i for s in emResult for i in s]
    else:
        foo = [i for s in kpResult for i in s]
        
    return foo

