import numpy as np
import scipy as sp

# Given column vectors Z, returns (possibly fewer) column vectors
# that form an orthonormal basis that spans Z
# using the Gram-Schmidt process
#
def getOrthonormalBasis(Z):
    
    R = np.zeros(np.shape(Z))
    tolerance = 1e-12
    vectorIndex = 0
    
    for i in range(0,np.shape(Z)[1]):
        v = Z[:,i]
        # remove components from v using the already obtained basis vectors
        for j in range(0,vectorIndex):
            v = v - np.sum(R[:,j]*v)*R[:,j]
            
        # there might be nothing left!
        # only add the vector if the residual norm is greater than the tolerance
        if (np.sum(v*v) >= tolerance):
            v = v / np.sqrt(np.sum(v*v))
            R[:,vectorIndex] = v
            vectorIndex += 1
        
    return R[:,0:vectorIndex]


# If P is a square row-stochastic matrix, returns the 
# normalized (sum 1) eigenvector associated with the 
# unit eigenvalue.
#
# For Markov chains with state transition probability 
# matrix P, this vector gives the steady state 
# distribution over the states.
#
# If P does not have a unit eigenvalue, returns None.
# However, sensible results are not guaranteed for
# arbitrary matrices.

def getSteadyStateDist(P):
    
    tolerance = 1e-10
    
    # need the left eigenvectors
    [u,v] = np.linalg.eig(np.transpose(P))
    v = np.transpose(v)
    
    index = 0
    for i in u:
        if np.abs(i - 1.) < tolerance:
            return np.real(v[index,:] / np.sum(v[index,:]))
        index += 1
        
    return None 

# Generates a binomial random number where
# - n is the number of trials, and
# - p is the success probability in any trial

def randBinomial(n,p):
    r = np.random.rand()
    
    # in future make this more efficient by precalculating the distribution
    # and allowing multiple outputs
    
    # the cdf ... we are trying to find the largest k so that p_comp < r
    k = 0
    p_comp = sp.misc.comb(n,k)*np.power(p,k)*np.power(1-p,n-k)
    
    while r > p_comp:
        k += 1
        p_comp += sp.misc.comb(n,k)*np.power(p,k)*np.power(1-p,n-k)
        
    return k

# p is a vector of probabilities; p must sum to 1
# returns a discrete random variable supported in {0,1,...,len(p)-1}
# with the distribution given by {p[0],p[1],...,p[len(p)-1]}
def randIntNonUniform(p):
    q = np.random.rand()
    r = 0
    for i in range(0,len(p)):
        r += p[i]
        if (q < r):
            return i

    # should never happen
    return -1

# the partial entropy function
# (some would define it as the negative of this)
def phi(p):
    if p == 0:
        return 0
    return p * np.log(p)

# m is a Markov chain transition probability matrix (row stochastic)
# returns the entropy rate of the Markov chain defined by m
# p allows you to provide an alternative steady state distribution
def markovChainEntropyRate(m,p=None):
    n = np.size(m,0) # number of rows
    
    # steady state distribution
    if (p is None):
        pi_ssd = getSteadyStateDist(m)
    else:
        pi_ssd = p
    
    # result
    r = 0
    
    # entropy for each row
    for i in range(0,n):
        row_ent = 0
        for j in range(0,n):
            row_ent -= phi(m[i,j])
        r += row_ent*pi_ssd[i]

    return r
        
