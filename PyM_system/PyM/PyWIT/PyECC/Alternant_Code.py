# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 15:34:31 2017

@author: admin_local
"""

#from PyECC_globals import *
import PyM.PyWIT.PyECC.PyECC_globals as PGl
from PyM.PyWIT.PyECC.pfactor import *

def vandermonde(a,r, i = 0):
    a = create_vector(a)
    n = a.size()
    v = a**i
    V = create_matrix(a.domain(),r,n)
    V[0]=v
    
    for k in range(1,r):
        v = a*v
        V[k]=v
    return V
vandermonde_matrix = vandermonde

def h_(C):
    return C._h

def a_(C):
    return C._alpha

def r_(C):
    return C._r

def b_(C):
    return C._beta

def G_(C):
    return C._G

def H_(C):
    return C._H

def set_G_(C,G):
    return C.set_G(G)

class Alternant_Code():
    """Alternant code class"""
    def __init__(self, h,alpha,r, K = ''):
        self._h = create_vector(h)
        self._alpha = create_vector(alpha)
        self._r = r
        self._construct_code(K)

    def _construct_code(self,K):
        if (self._alpha.is_invertible()):
            self._beta = 1/self._alpha
        else:
            self._beta = None
        self._G = None
        V = vandermonde(self._alpha,self._r)
        for i in range(ncols(V)):
            for j in range(len(V)):
                V[j,i] = V[j,i] * self._h[i]
        self._H = V
        if K != '':
            self._K = K
        else:
            D = self._H.domain().subdomain()
            if D == D.base_domain():
                self._K = D
            else:
                self._K = D.subdomain()
    
    def __str__(self):
        res = "Alternant Code "
        if PGl.TYPES == 1:
            res += "with alpha = "+self._alpha.__str__()+" h = "+self._h.__str__()+" and r = "+self._r.__str__()
        return res
    
    def __repr__(self):
        return self.__str__()
    
    def __hash__(self):
        return PGl.HAC + hash(self._h) + hash(self._alpha) + r
    
    def h_(self):
        return self._h
    
    def a_(self):
        return self._alpha
    
    def r_(self):
        return self._r
    
    def b_(self):
        return self._beta
    
    def G_(self):
        return self._G
    
    def H_(self):
        return self._H
    
    def K_(self):
        return self._K
    
    def set_G(self,G):
        self._G = G.cp()
#
AC = alternant_code = Alternant_Code

#class Reed_Salomon_Code(Alternant_Code):
#    """RSC as a alternant code class"""
#    def __init__(self,alpha,k):
#        self._alpha = create_vector(alpha)
#        self._r = self._alpha.size() - k
#        aux = []
#        for i in range(self._alpha.size()):
#            val = 1
#            for j in range(self._alpha.size()):
#                if i == j:
#                    continue
#                val = val*(self._alpha[j]-self._alpha[i])
#            aux.append(val)
#        h = create_vector(aux)
#        self._h = 1/h
#        self._construct_code()
#    def __str__(self):
#        return "Reed Salomon Code as an Alternant Code with alpha = "+self._alpha.__str__()+" h = "+self._h.__str__()+" and r = "+self._r.__str__()
#
    
    
def Reed_Salomon_Code(a,k):
    n = len(a)
    a = create_vector(a)
    h = []
    for i in range(n):
        hi = 1
        for j in range(n):
            if j == i: continue
            hi = hi*(a[j]-a[i])
        h += [hi]
    h = create_vector(h)
    h = 1/h
    C = AC(h,a,n-k)
    C._G = vandermonde(a,k)
    return C

RS = Reed_Salomon_Code

def BCH_Code(alpha,d,l=1,K=''):
    print("kaki-1")
    n = order(alpha)
    print("kaki0")
    h = geometric_series(alpha**l,n)
    print("kaki1")
    a = geometric_series(alpha, n)
    print("kaki2")
    return AC(h,a,d-1,K)
# alias    
BCH = BCH_Code

#class Goppa_Classic_Code(Alternant_Code):
#    """Goppa classic code as a alternant code class"""
#    def __init__(self,g,alpha):
#        self._r = g.degree()
#        self._alpha = create_vector(alpha)
#        aux = []
#        for i in range(self._alpha.size()):
#            aux.append(1/g.evaluate(self._alpha[i]))
#        self._h = create_vector(aux)
#        self._construct_code()
#    
#    def __str__(self):
#        return "Goppa Classic Code as an Alternant Code with alpha = "+self._alpha.__str__()+" h = "+self._h.__str__()+" and r = "+self._r.__str__()
#

def Goppa_Classic_Code(g,alpha):
        alpha = create_vector(alpha)
        aux = []
        for i in range(alpha.size()):
            aux.append((1/g.evaluate(alpha[i]))._value)
        h = create_vector(aux)
        return AC(h,alpha,g.degree())

Goppa=goppa=classical_goppa_code=Goppa_Classic_Code
                  
def GRS(h,a,k):
    n = len(h)
    if len(a)!=n: return 'GRS error: length mismatch'
    C = AC(h,a,n-k)
    V = alternant_matrix(h,a,n-1)
    h1 = left_kernel(V)[0]
    G = alternant_matrix(h1,a,k)
    C._G = G
    return C            
            
def PRS(K,k):
    n = cardinal(K)-1
    alpha = primitive_root(K) 
    C = BCH(alpha,n-k+1)
    a = a_(C)
    C._G = vandermonde(a,k)
    return C            
            