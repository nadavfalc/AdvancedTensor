# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:44:57 2018

@author: UPC-ESAII
"""

from PyM.PyWIT.PyECC.Poly import *
from PyM.PyWIT.PyECC.Matrix import *
from PyM.PyWIT.PyECC.Fq_type import *

import time

def field(x): return x.min_domain()

def pick_element(k,A,variable = ['x']):
    if is_Zn(A):
        n = cardinal(A)
        return (k%n)>>A
    B = extension_base(A)   # A = B[x] = B[X]/(f), f = X**r+···
#    x = generator(A)
#    r = degree(A)            # cardinal(A) = cardinal(B)**r
    b = cardinal(B)
    K = convert(k,b)
    S = []
    for i in range(len(K)):
        S += [pick_element(K[i],B)]
    #D = [pick_element(d,B) for d in D]
    if isinstance(A,Poly):
        return A.element_list(S,variable)
    return A.element(S)#hohner(D,x)
#
element = pick = pick_element    

def pick_nonzero_element(k,A):
    if is_Zn(A):
        n = cardinal(A)
        return (k%n)>>A
    B = extension_base(A)   # A = B[x] = B[X]/(f), f = X**r+···
#    x = generator(A)
#    r = degree(A)            # cardinal(A) = cardinal(B)**r
    b = cardinal(B)
    K = convert(k,b)
    S = []
    for i in range(1,len(K)):
        S += [pick_element(K[i],B)]
    #D = [pick_element(d,B) for d in D]
    return A.element(S)#hohner(D,x)
#
nonzero_element = pick_nonzero = pick_nonzero_element  

def cardinal(A): return A.cardinal()

def characteristic(F): return F._characteristic

def extension_base(A): return A.subdomain()

def K_(f): 
    if isinstance(f,Base_type.Base_type):
        return f.K_()
    if isinstance(f,Base.Base):
        return f.K_()
    if isinstance(f,int):
        return ZZ()
    if isinstance(f,float):
        return QQ(ZZ())
    if isinstance(f,Symbol_type):
        return ZZ()
#    if isinstance(f,Alternant_Code):
#        return f.K_()
    if isinstance(f,list):
        return K_(vector(f))
    return "Type not found"

# cal definir is_Zn(A), true quan A=Zn(n) per algun n, False altrament
def is_Zn(A): return isinstance(A,Zn)#==Zn

def md_power(f,n,g): 
    return g.domain()._md_pow_(f,n,g)

def md_double_power(f,a,b,g): 
    return g.domain()._md_double_pow_(f,a,b,g)

def md_mul(f,g,h): 
    return h.domain()._md_mul_(f,g,h)      
md_mult = md_mul   
         
def convert(x,b):
    if x<b:
        return [x]
    else:
        return convert(x//b,b) + [x%b]
    
    
def is_irreducible(f,K=''):
    # error if K is not a field and f not univariate polynomial
    # if coefficient_ring(f) != K: f = reduce(f,K)
    # (in case f is given with integer coefficients)
    
    if isinstance(f,list):
        K = K_(f)
        P = Poly(K)
        aux = []
        x = Symbol_type('x')
        for i in range(len(f)):
            if f[i] != 0:
                aux.append(Monomial_type(f[i]>>K,[x],[len(f) - 1 - i],P))
        if len(aux) == 0:
            ff = Monomial_type(0>>K,[x],[0],P)
        elif len(aux) == 1:
            ff = aux[0]
        else:
            ff = Poly_type(aux,[x],P,aux[0]._degree, aux[0]._totIndex)
        return is_irreducible(ff,K)    
    if isinstance(f,Monomial_type):
        if f._totIndex < 2:
            return True
        else:
            return False
    if isinstance(f,Symbol_type):
        return True
    if K == '': 
        K = K_(f)
    if (f.domain().subdomain().is_sub_domain(K)) and f.domain().subdomain()!=K:
        for i in range(len(f._value)):
            if not K.is_element(f._value[i]._value.transform_to_min_domain()):
                raise TypeError("This polynomial is not defined in this domain")
        aux = []
        aind = []
        Dom = Poly(K)
        if f._varOrder == None:
            for i in range(len(f._value)):
                aux.append(Monomial_type(f._value[i]._value, f._variables[:],f._value[i]._index[:],Dom))
            f = Poly_type(aux,f._variables[:],Dom,f._degree, f._mdegree)
        else:
            for i in range(len(f._value)):
                aux.append(Monomial_type(f._value[i]._value, f._variables[:],f._value[i]._index[:],Dom,f._varOrder[:]))
        f = Poly_type(aux,f._variables[:],Dom,f._degree, f._mdegree,f._ord[:], f._varOrder[:])
    if not K.is_sub_domain(f._structure.subdomain()):
        raise TypeError("This polynomial cannot be computed in this domain")
    if K.is_sub_domain(f.K_()) and K!=f.K_():
        D = Poly(K)
        if isinstance(f,Monomial_type):
            if f._varOrder == None:
                s = f._value >> K
                if s != 0:
                    f = Monomial_type(s, f._variables[:],f._index[:],D)
                else:
                    f = Monomial_type(s, f._variables[:],[0]*(len(f._variables)),D)
            else:
                if s != 0:
                    f = Monomial_type(f._value >> K, f._variables[:],f._index[:],D,f._varOrder[:])
                else:
                    f = Monomial_type(s, f._variables[:],[0]*(len(f._variables)),D)
        elif isinstance(f,Poly_type):
            aux = []
            if f._varOrder == None:
                for i in range(len(f._value)):
                    s = f._value[i]._value >> K
                    if s != 0:
                        aux.append(Monomial_type(s, f._variables[:],f._value[i]._index[:],D))
                if len(aux) == 0:
                    f = Monomial_type(0>>K, f._variables[:],[0]*(len(f._variables)),D)   
                elif len(aux) == 1:
                    f = aux[0]
                else:
                    f = Poly_type(aux,f._variables[:],D,f._degree,f._mdegree)
            else:
                for i in range(len(f._value)):
                    s = f._value[i]._value >> K
                    if s != 0:
                        aux.append(Monomial_type(s, f._variables[:],f._value[i]._index[:],D,f._varOrder[:]))
                if len(aux) == 0:
                    f = Monomial_type(0>>K, f._variables[:],[0]*(len(f._variables)),D)   
                elif len(aux) == 1:
                    f = aux[0]
                else:
                    f = Poly_type(aux,f._variables[:],D,f._degree,f._mdegree,f._ord[:],f._varOrd[:])
                    
    if len(f._variables) > 1:
        raise NotImplemented("Only implemented for univariate polynomials")
    if K.characteristic() == 0:
        faux = f
        while True:
            lc = faux.constant_coeff()>>ZZ()
            if isinstance(lc,int):
                break
            elif isinstance(lc,Monomial_type):
                if isinstance(lc._value,int):
                    lc = lc._value
                    break
            faux = lc
        
        D = divisors(lc)
        for d in D:
            if f.eval_poly(d) == 0:
                return False
        return True
    
    n =  f.degree()#n = degree(f)
    D = f.domain()
    if n == 0: return False
    if n == 1: return True
    X = f.variable()>>f.domain()  #X = variable(f)
    q = K.cardinal()  #q = cardinal(K)
    
    P = list(set(prime_factors(n)))
    h = X;
    nant = 0
    for i in range(len(P)):
        nact = n//P[len(P)-1 - i]
        #r = D.variable_d(q**(n//p))
        h = md_power(h,q**(nact-nant),f)
        a1 = D.gcd(h-X,f)
        if a1 != 1:
            return False
        nant = nact
    h = md_power(h,q**(n-nact),f)
    if h == X:
        return True
    return False