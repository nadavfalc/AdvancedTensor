# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 17:47:22 2020

@author: narcis
"""

from PyM.PyWIT.wit_mor import *

def chi(X,F,T='T'):
    if dim(X) == 0:
        return rk(F)
    return integral(X,HRR(X,F,T))

def HRR(X,F,T='T'):
    n = dim(X)
    if todd_class(X) == null_vector():
        X._todd_class = todd_vector(tan_bundle(X))
    if n == 0:
        return rk(F)
    else:
        P = collect(chern_character(F,T)*Td(X,T),T)
        return coeffs_inc(P)[n]