## cdi
## SXD 0908, 0914, 1029, 1124, 1201, 0219
'''
CDI15 utilities
'''

#random tools
import random as rd
import numpy as np
import copy as cp

''' Probability '''

# Checks whether W is a positive list with at least 2 elements
def weight_list(W):
    if len(W)<2: return False
    for w in W:
        if w <= 0: return False
    return True

# Returns the probability distribution of a weight distribution
def normalize(W):
    if not weight_list(W): 
        return 'cdi/normalize: input is not a weight distribution'
    S = sum(W)
    return [w/S for w in W]
    
# Intersection union and difference of two lists, ranges or tuples
def intersect(X,Y): return set(X) & set(Y)
def union(X,Y): return set(X) | set(Y)
def difference(X,Y): return set(X)-set(Y)

#A=range(5); B=list(range(7,2,-1)); C={3,4,7,11,15}

# Checks whether X is a subset of U
def subset(X,U): 
    if intersect(X,U)==set(X): return True
    else: return False

# Probability of X with respect to list of weights
# X has to be a subset of range(1,len(w)+1).
def probability(X,w):
    if not weight_list(w): return 'cdi/probability: Wrong weight list'
    U = range(1,len(w)+1)
    if not subset(X,U): return 'cdi/probability: Wrong event data'
    p = 0
    for k in X:
        p += w[k-1]
    return p/sum(w)

# Probability of Y conditioned to the occurrence of X
# with with respect to w.
def conditional_probability(Y,X,w): 
    if list(X) == []: 
        return 'cdi/conditional_probability: Wrong data for conditioning event'
    return probability(intersect(X,Y),w)/probability(X,w)

# Independence of X and Y for probabilities defined by w.   
def independent(X,Y,w): 
    x = probability(X,w)
    y = probability(Y,w)
    xy = probability(intersect(X,Y),w)
    if xy == x*y: return True
    else: return False


# Interaction index of X and Y for probabilities defined by w.
def interaction(X,Y,w):
    if list(X) == [] or list(Y) == []: 
        return 'cdi/interaction: events have to be non empty'
    return probability(intersect(X,Y),w)/probability(X,w)/probability(Y,w)
    
# Bayes formula
def bayes(Y,X,w): return probability(Y,w)*interaction(X,Y,w)

# Average of a list of quantities E weighted with respect
# to given list of weights w.
def weighted_average(E,w):
    if not weight_list(w): 
        return 'cdi/weighted_average: Wrong weight list'
    if len(E) != len(w): 
        return 'cdi/weighted_average: unequal lengths'
    a = 0
    for k in range(0,len(E)):
        a += w[k]*E[k]
    return a/sum(w)

# average(X): X can be a numerical list, tuple or set.
def average(E): 
    if len(E)==0: 
        return 'cdi/average: E has to be non-empty'
    return sum(E)/len(E) 

# With N independent repeats, probability that an event
# of probability p occurs at least once  
def at_least_once_in(N,p): return 1-(1-p)**N
at_least_once=at_least_once_in



def rd_int(a=0, b=1):
    return rd.randint(a,b)
    
def ri(n): 
    if n == 1:
        return rd_int(0,9)
    else:
        return rd_int(10**(n-1),10**n)

def rd_float(a=0, b=1):
    return rd.uniform(a,b)

def rd_gauss(mu = 0, sigma = 1):
    return rd.gauss(mu,sigma)

#def rd_choice(X,k = 1):
#    Y = []
#    for _ in range(k):
#        Y +=[rd.choice(X)]
#    return Y

def rd_choice(A, m=1): 
    if isinstance(A,int): return rd_choice(list(range(A)),m)
    if len(A)==0: return None
    if len(A)<m: return "rd_choice: not enough elements to choose from"
    X = list(cp.deepcopy(A))
    Y = []
    for j in range(m):
        t = rd_int(0,len(X)-1)
        x = X[t]
        Y += [x]
        X = X[:t] + X[t+1:]
    return Y
rd_selection = rd_sample = rd_choice

#def rd_sample(X,k = 1):
#    return rd.sample(X,k)

def rd_message(A,N = 1):
    Y = ''
    for _ in range(N):
        Y +=rd.choice(A)
    return Y


''' Information theory '''

L2 = lambda x: np.log(x)/np.log(2)

# Entropy of a weight distribution
def entropy(W):
    S = sum(W)
    if S == 0: return 'entropy(W): W is not a weight distribution' 
    H = 0
    for w in W:
        if w < 0: return 'entropy(W): W is not a weight distribution'
        if w == 0: continue
        else: H -= w * L2(w)
    H += S * L2(S)
    return H/S
'''Remark:
This function yields the entropy of the probability distribution
p_j = W[j]/S, that is, sum_j -p_j * log2(p_j). The proof is a
simple computation.
'''

def joint_entropy(P): return entropy(unfold(P))

# Marginals

def row_marginal(P): return [sum(p) for p in P]

def col(k,P): return [P[j][k] for j in range(len(P))]

def col_marginal(P): return [sum(col(k,P)) for k in range(len(P[0]))]


'''Prefix codes'''

# ells2ens(L) transforms a list of positive lenghts 
# to a table with their frequencies.

def ells2ens(L): 
    return {j:L.count(j) for j in range(1,1+max(L))}


# ens2ells(N), where N = {1:n1,...,l:nl} 
# with the nj non negative integers, 
# yields the list [l1,...,ln] in which 
# j in [1,...,l] is repeated nj times.

def ens2ells(N):
    L = []
    for j in N:
        L += N[j] * [j]
    return L


# cartesian(X,C) yields the list of strigns
# obtained by concatenating, for x in X and c in C,
# str(x) and str(c), with the convention that 
# if C (or X) is empty, then we get the list of
# str(x) (respectively str(c)). 

def cartesian(X,C):
    if not C:
        return sorted({str(x) for x in X})
    if not X:
        return sorted({str(c) for c in C})
    return sorted({str(x)+str(c) for x in X for c in C})


# Define a function that constructs, given r and a feasible list
# of lengths, a prefix encoding with those lengths

# to check whether Kraft's inequality holds
def kraft(L, r=2): return sum(r**(-l) for l in L) <= 1

def ells2code(L,r=2):
    sorted(L)
    if not kraft(L,r): 
        return "ells2code: the list does not satisfy Kraft's inequality" 
    return( make_code( list(ells2ens(L).values()), r ) )

def make_code(N,r=2):
    #sorted(N)
    n = sum(N)
    l = len(N)
    A = range(n)
    C = range(r)
    T = dict()
    X = []  # list of length j-1 codes 
    Y = []  # list of length j-1 words eligible for constructing length j codes
    sj=0
    for j in range(l):
        nj = N[j]
        if nj==0:
           Y = cartesian(Y,C)
           continue 
        Y = cartesian(Y,C) 
        #print(Y)
        X = Y[:nj]
        Y = list(set(Y)-set(X))
        T.update({i:X[i-sj] for i in range(sj,sj+nj)})
        sj = sj+nj
    return T


''' Arithmetic coding '''

# binary numbers
def dec2bin(x, nb = 58):
    if (x<0) | (x>1): 
        return 'dec2bin: was expecting a number in [0,1]'
    if x == 1: return nb*'1'
    xb = ''
    for j in range(nb):
        x = 2*x
        if x < 1: xb += '0'
        else: 
            xb += '1'
            x -= 1
    return xb

def bin2dec(xb):
    x =0.0; j = 1
    for b in xb:
        b = int(b)
        if b == 0: pass
        else: x += b/2**j
        j += 1
    return x

# Segment defined by a binary string
def segment(x):
    l = bin2dec(x)
    u = 2**(-len(x))
    h = l + u
    return (l,h)




''' Filters '''

def filterx(h,x):
    m = len(h); n = len(x)
    y = []
    for k in range(m+n-1):
        a = max(0,k-n+1); b = min(m,k+1)
        s = sum(h[j]*x[k-j] for j in range(a,b))
        y += [s]
    for i in range(n,n+m-1):
        y[i%n] += y[i]
    return y[:n]

def up_sample(a):
    x = []
    for t in a:
        x += [t,0]
    return x

''' Signals '''

# sample(h,N, a=0, b=1): if h is a function defined on [a,b] 
# and N is a positive integer, it returns the list of 
# samples h(a+j*s) of h for j = 0,...,N, s = (b-a)/N
def sample(h,N,a=0,b=1):
    s = (b-a)/N
    return [h(a+j*s) for j in range(N+1)]

# energy(X): X can be a numerical list, tuple or set.
def energy(f): return sum(x*x for x in f)

def energy_map(x):
    y = []
    s = 0
    for t in x:
        s += t**2
        y += [s]
    return y

# If f is a sequence, decimate(f) 
# returns [f[1],f[3],...]. Instead of starting
# at the index 1, we can provide a starting value
# as the second argument
def decimate(f,start=1): return [f[j] for j in range(start, len(f),2)]

# It returns h reversed and with an alternative change of sign.
# [a,b,c] => [c,-b,a]. 
def dual(h): 
    s = 1
    hd = []
    for x in reversed(h):
        hd += [s*x]
        s = -s
    return hd


''' Miscellaneous '''

# a as percent of b>0
def percent(a,b): return 100*a/b

# 
def unfold(P):
    X = []
    for p in P:
        X += p
    return X

# Rounding  
def is_sequence(arg):
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))

# To round a float to n decimal places. It also works
# for a list or a matrix of floats.
def roundint(f,n):
    if isinstance(f,int): return float(f)
    if isinstance(f,float): return float("%.*f" % (n,f))
    if is_sequence(f): return [roundint(x,n) for x in f]  
 