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

def transition(v,firstState,lastState):
  # find a specific transition
  currentLocation = 1
  transitionIndices = []

  while (currentLocation < len(v)):

    if (v[currentLocation-1] == firstState) and (v[currentLocation] == lastState):
      transitionIndices.append(currentLocation)

    currentLocation += 1

  return transitionIndices


def isNonPermissiveClosing(v,index):

    # v[index] should be '5.0'
    if v[index] != '5.0':
        return False

    n = nextOpen(v,index)

    if n == np.Inf:
        return False

    if v[n] == '3.0':
        return True

    return False

# state is a list, we return the index of the next appearance of a state in the list
def getIndexOfNextState(v,start,state):
  currentLocation = start+1
  
  while currentLocation < len(v):
    if (v[currentLocation] not in state):
      currentLocation += 1
    else:
      return currentLocation
    
  return np.Inf # this makes it easier to compare

def nextOpen(v,index):
    return getIndexOfNextState(v,index,['3.0','4.0'])

def countClosings(v):
    
    # 3 to 2 cannot possibly be nonpermissive
    c32 = transition(v,'3.0','2.0')
    c45 = transition(v,'4.0','5.0')
    
    npc45 = 0
    for i in c45:
        if isNonPermissiveClosing(v,i):
            npc45 += 1
    
    return [len(c32)+len(c45),npc45]

def bpr_md_new(r,conf=False,getEmResults=True):
    
    t = transition(r.states,'4.0','5.0')
    
    md_n = 0
    md = 0
    mdc_n = 0
    mdc = 0
    
    # remember that emEstimates is nominally a list with one entry per iteration
    # but when looking for permissive closures, we are only concerned with the last EM estimate
    if getEmResults is True:
        comparison = r.emEstimates[-1]
        comparison_conf = r.emEstimates_conf[-1]
    else:
        comparison = r.kp
        comparison_conf = r.kp_conf
    
    # fix variable names
    for i in range(0,len(t)):
        if (isNonPermissiveClosing(r.states,t[i])):
            md_n += 1
            if isNonPermissiveClosing(comparison,t[i]):
                md += 1
            qqq = nextOpen(comparison,t[i])
            if (qqq != np.Inf):
                if comparison_conf[qqq] != '-1.0' and comparison_conf[t[i]] != '-1.0':
                    if isNonPermissiveClosing(comparison,t[i]):
                        mdc += 1

    mdc_n = md_n

    if conf is True:
        return [mdc_n,mdc]
        
    return [md_n,md]

def bpr_fa_new(r,conf=False,getEmResults=True):

    fa_n = 0
    fa = 0
    fac_n = 0
    fac = 0

    # remember that emEstimates is nominally a list with one entry per iteration
    # but when looking for permissive closures, we are only concerned with the last EM estimate
    if getEmResults is True:
        comparison = r.emEstimates[-1]
        comparison_conf = r.emEstimates_conf[-1]
    else:
        comparison = r.kp
        comparison_conf = r.kp_conf

    t = transition(comparison,'4.0','5.0')

    for i in range(0,len(t)):
        if (isNonPermissiveClosing(comparison,t[i])):
            fa_n += 1
            if isNonPermissiveClosing(r.states,t[i]):
                fa += 1
            qqq = nextOpen(comparison,t[i])
            if qqq != np.Inf:
                if comparison_conf[qqq] != '-1.0' and comparison_conf[t[i]] != '-1.0':
                    fac_n += 1
                    if isNonPermissiveClosing(r.states,t[i]):
                        fac += 1
    
    if conf is True:
        return [fac_n,fac]
        
    return [fa_n,fa]

def openToCloseIntervals(r,openings=True):

    openStates = ['3.0','4.0']
    closedStates = ['0.0','1.0','2.0','5.0','6.0']
    
    index = 0
    intervalCounter = 0
    openIntervals = []
    closedIntervals = []
    
    if r.emEstimates[-1][0] in openStates:
        indexOpen = True
    else:
        indexOpen = False
        
    for i in range(1,len(r.emEstimates[-1])):
        if indexOpen is True:
            if r.emEstimates[-1][i] in openStates:
                intervalCounter += 1
            else:
                indexOpen = False
                openIntervals.append(intervalCounter)
                intervalCounter = 0
        else:
            if r.emEstimates[-1][i] in closedStates:
                intervalCounter += 1
            else:
                indexOpen = True
                closedIntervals.append(intervalCounter)
                intervalCounter = 0
    
    if openings is True:
        return openIntervals
        
    return closedIntervals

#def bpr_md(r,doOpenStates=True,conf=False,getEmResults=True):
#
#
#    # count number of transitions in the ground truth sequence
#    #for i in range(0,len(r.kp)):
#    #    print(str(i) + ' ' + r.kp[i])
#
#    openTransition_v = [['2.0','3.0'],['5.0','4.0']]
#    closeTransition_v = [['3.0','2.0'],['4.0','5.0']]
#    openStates = ['3.0','4.0']
#    closedStates = ['0.0','1.0','2.0','5.0','6.0']
#    open_gt = []
#    open_kp = []
#    open_em = []
#
#    states = r.states
#
#    if conf is True:
#        kp = r.kp_conf
#        em = r.emEstimates_conf[-1]
#    else:
#        kp = r.kp
#        em = r.emEstimates[-1]
#
#    # this code was initially written with channel openings in mind
#    # so some of the variable names refer only to openings, e.g. open_kp, open_em
#    # however, if doOpenStates is false, we handle channel closings rather than openings
#    if doOpenStates is True:
#        iterate_v = openTransition_v
#    else:
#        iterate_v = closeTransition_v
#
#    for ot in iterate_v:
#        # transitions in the ground truth sequence
#        open_gt.append(r.transition(states,ot[0],ot[1]))
#
#        # transitions in the estimated sequence with known probabilities
#        open_kp.append(r.transition(kp,ot[0],ot[1]))
#
#        # transitions in the EM estimated sequence
#        open_em.append(r.transition(em,ot[0],ot[1]))
#
#    #print(open_gt)
#    #print(open_kp)
#    #print(open_em)
#
#    # count the number correctly detected
#    gtLen = [len(open_gt[0]),len(open_gt[1])]
#    kpCorrectlyDetected = [0]*len(open_gt)
#    emCorrectlyDetected = [0]*len(open_gt)
#
#    for i in range(0,len(open_gt)):
#        for j in open_gt[i]:
#            if j in open_kp[i]:
#                kpCorrectlyDetected[i] += 1
#            if j in open_em[i]:
#                emCorrectlyDetected[i] += 1
#
#    #print(gtLen)
#    #print(kpCorrectlyDetected)
#    #print(emCorrectlyDetected)
#
#    # count the number mistaken for hte other transition
#    kpOtherTransition = [0]*len(open_gt)
#    emOtherTransition = [0]*len(open_gt)
#
#    # did we mistake one type of opening for another?
#    for i in range(0,len(open_gt)):
#        for j in open_gt[i]:
#            for k in range(0,len(open_gt)):
#                if (i != k):
#                    if j in open_kp[k]:
#                        kpOtherTransition[i] += 1
#                    if j in open_em[k]:
#                        emOtherTransition[i] += 1
#
#    #print(kpOtherTransition)
#    #print(emOtherTransition)
#
#    # count the number where the ground truth has a transition, but the estimate stays open or closed
#    kpStaysOpen = [0]*len(open_gt)
#    emStaysOpen = [0]*len(open_gt)
#    kpStaysClosed = [0]*len(open_gt)
#    emStaysClosed = [0]*len(open_gt)
#
#    for i in range(0,len(open_gt)):
#        for j in open_gt[i]:
#            if (kp[j-1] in openStates) and (kp[j] in openStates):
#                kpStaysOpen[i] += 1
#            if (kp[j-1] in closedStates) and (kp[j] in closedStates):
#                kpStaysClosed[i] += 1
#
#        for j in open_gt[i]:
#            if (em[j-1] in openStates) and (em[j] in openStates):
#                emStaysOpen[i] += 1
#            if (em[j-1] in closedStates) and (em[j] in closedStates):
#                emStaysClosed[i] += 1
#
#    #print(kpStaysOpen)
#    #print(emStaysOpen)
#    #print(kpStaysClosed)
#    #print(emStaysClosed)
#
#    # count the number where the ground truth has a transition, but the estimate has the opposite transition
#    kpWrongWay = [0]*len(open_gt)
#    emWrongWay = [0]*len(open_gt)
#
#    # did we mistake one type of opening for another?
#    for i in range(0,len(open_gt)):
#
#        if doOpenStates is True:
#            for j in open_gt[i]:
#                if (kp[j] in closedStates) and (kp[j-1] in openStates):
#                    kpWrongWay[i] += 1
#            for j in open_gt[i]:
#                if (em[j] in closedStates) and (em[j-1] in openStates):
#                    emWrongWay[i] += 1
#        else:
#            for j in open_gt[i]:
#                if (kp[j] in openStates) and (kp[j-1] in closedStates):
#                    kpWrongWay[i] += 1
#            for j in open_gt[i]:
#                if (em[j] in openStates) and (em[j-1] in closedStates):
#                    emWrongWay[i] += 1
#
#
#
#    #print(kpWrongWay)
#    #print(emWrongWay)
#
#    # assemble results
#    # total, correct, other transition, stays open, stays closed, wrong way
#    kpResult = []
#    emResult = []
#
#    # for each type of transition in open_gt, give a list of results
#    for i in range(0,len(open_gt)):
#        kpResult.append([len(open_gt[i]),
#                         kpCorrectlyDetected[i],
#                         kpOtherTransition[i],
#                         kpStaysOpen[i],
#                         kpStaysClosed[i],
#                         kpWrongWay[i]])
#        emResult.append([len(open_gt[i]),
#                         emCorrectlyDetected[i],
#                         emOtherTransition[i],
#                         emStaysOpen[i],
#                         emStaysClosed[i],
#                         emWrongWay[i]])
#
#    # the above is a list of lists, the code below flattens the list of lists into a single list
#    if getEmResults is True:
#        foo = [i for s in emResult for i in s]
#    else:
#        foo = [i for s in kpResult for i in s]
#
#    return foo
#
## false alarms are different enough to warrant a new method
#def bpr_fa(r,doOpenStates=False,conf=False,getEmResults=True):
#
#
#    # count number of transitions in the ground truth sequence
#    #for i in range(0,len(r.kp)):
#    #    print(str(i) + ' ' + r.kp[i])
#
#    openTransition_v = [['2.0','3.0'],['5.0','4.0']]
#    closeTransition_v = [['3.0','2.0'],['4.0','5.0']]
#    openStates = ['3.0','4.0']
#    closedStates = ['0.0','1.0','2.0','5.0','6.0']
#    open_gt = []
#    open_kp = []
#    open_em = []
#
#    states = r.states
#
#    if conf is True:
#        kp = r.kp_conf
#        em = r.emEstimates_conf[-1]
#    else:
#        kp = r.kp
#        em = r.emEstimates[-1]
#
#    # this code was initially written with channel openings in mind
#    # so some of the variable names refer only to openings, e.g. open_kp, open_em
#    # however, if doOpenStates is false, we handle channel closings rather than openings
#    if doOpenStates is True:
#        iterate_v = openTransition_v
#    else:
#        iterate_v = closeTransition_v
#
#    for ot in iterate_v:
#        # transitions in the ground truth sequence
#        open_gt.append(r.transition(states,ot[0],ot[1]))
#
#        # transitions in the estimated sequence with known probabilities
#        open_kp.append(r.transition(kp,ot[0],ot[1]))
#
#        # transitions in the EM estimated sequence
#        open_em.append(r.transition(em,ot[0],ot[1]))
#
#    #print(open_gt)
#    #print(open_kp)
#    #print(open_em)
#
#    # count the number correctly detected
#    gtLen = [len(open_gt[0]),len(open_gt[1])]
#    kpCorrectlyDetected = [0]*len(open_gt)
#    emCorrectlyDetected = [0]*len(open_gt)
#
#    for i in range(0,len(open_gt)):
#        for j in open_kp[i]:
#            if j in open_gt[i]:
#                kpCorrectlyDetected[i] += 1
#        for j in open_em[i]:
#            if j in open_gt[i]:
#                emCorrectlyDetected[i] += 1
#
#    #print(gtLen)
#    #print(kpCorrectlyDetected)
#    #print(emCorrectlyDetected)
#
#    # count the number mistaken for hte other transition
#    kpOtherTransition = [0]*len(open_gt)
#    emOtherTransition = [0]*len(open_gt)
#
#    # did we mistake one type of opening for another?
#    for i in range(0,len(open_gt)):
#        for j in open_kp[i]:
#            for k in range(0,len(open_gt)):
#                if (i != k) and (j in open_gt[k]):
#                    kpOtherTransition[i] += 1
#
#        for j in open_em[i]:
#            for k in range(0,len(open_gt)):
#                if (i != k) and (j in open_gt[k]):
#                    emOtherTransition[i] += 1
#
#    #print(kpOtherTransition)
#    #print(emOtherTransition)
#
#    # count the number where the ground truth has a transition, but the estimate stays open or closed
#    kpStaysOpen = [0]*len(open_gt)
#    emStaysOpen = [0]*len(open_gt)
#    kpStaysClosed = [0]*len(open_gt)
#    emStaysClosed = [0]*len(open_gt)
#
#    for i in range(0,len(open_gt)):
#        for j in open_kp[i]:
#            if (states[j-1] in openStates) and (states[j] in openStates):
#                kpStaysOpen[i] += 1
#            if (states[j-1] in closedStates) and (states[j] in closedStates):
#                kpStaysClosed[i] += 1
#
#        for j in open_em[i]:
#            if (states[j-1] in openStates) and (states[j] in openStates):
#                emStaysOpen[i] += 1
#            if (states[j-1] in closedStates) and (states[j] in closedStates):
#                emStaysClosed[i] += 1
#
#    #print(kpStaysOpen)
#    #print(emStaysOpen)
#    #print(kpStaysClosed)
#    #print(emStaysClosed)
#
#    # count the number where the ground truth has a transition, but the estimate has the opposite transition
#    kpWrongWay = [0]*len(open_gt)
#    emWrongWay = [0]*len(open_gt)
#
#    # did we mistake one type of opening for another?
#    for i in range(0,len(open_gt)):
#
#        for j in open_kp[i]:
#            if doOpenStates is True:
#                if (states[j] in closedStates) and (kp[j-1] in openStates):
#                    kpWrongWay[i] += 1
#            else:
#                if (states[j] in openStates) and (kp[j-1] in closedStates):
#                    kpWrongWay[i] += 1
#
#        for j in open_em[i]:
#            if doOpenStates is True:
#                if (states[j] in closedStates) and (kp[j-1] in openStates):
#                    emWrongWay[i] += 1
#            else:
#                if (states[j] in openStates) and (kp[j-1] in closedStates):
#                    emWrongWay[i] += 1
#
#
#    #print(kpWrongWay)
#    #print(emWrongWay)
#
#    # assemble results
#    # total, correct, other transition, stays open, stays closed, wrong way
#    kpResult = []
#    emResult = []
#    for i in range(0,len(open_gt)):
#        kpResult.append([len(open_kp[i]),
#                         kpCorrectlyDetected[i],
#                         kpOtherTransition[i],
#                         kpStaysOpen[i],
#                         kpStaysClosed[i],
#                         kpWrongWay[i]])
#        emResult.append([len(open_em[i]),
#                         emCorrectlyDetected[i],
#                         emOtherTransition[i],
#                         emStaysOpen[i],
#                         emStaysClosed[i],
#                         emWrongWay[i]])
#
#    if getEmResults is True:
#        foo = [i for s in emResult for i in s]
#    else:
#        foo = [i for s in kpResult for i in s]
#
#    return foo

def openToCloseIntervals(r,openings=True):

    openStates = ['3.0','4.0']
    closedStates = ['0.0','1.0','2.0','5.0','6.0']
    
    index = 0
    intervalCounter = 0
    openIntervals = []
    closedIntervals = []
    
    if r.emEstimates[-1][0] in openStates:
        indexOpen = True
    else:
        indexOpen = False
        
    for i in range(1,len(r.emEstimates[-1])):
        if indexOpen is True:
            if r.emEstimates[-1][i] in openStates:
                intervalCounter += 1
            else:
                indexOpen = False
                openIntervals.append(intervalCounter)
                intervalCounter = 0
        else:
            if r.emEstimates[-1][i] in closedStates:
                intervalCounter += 1
            else:
                indexOpen = True
                closedIntervals.append(intervalCounter)
                intervalCounter = 0
    
    if openings is True:
        return openIntervals
        
    return closedIntervals
    
    

