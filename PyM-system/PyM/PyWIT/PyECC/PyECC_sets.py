# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 17:11:17 2018

@author: narcis
"""

import PyM.PyWIT.PyECC.PyECC_globals as PGl

def PyECC_set(variable, value):
    if hasattr(PGl, variable):
        setattr(PGl,variable,value)
    else:
        print('PyECC set error: Variable not found')

def PyECC_get(variable):
    if hasattr(PGl, variable):
        return getattr(PGl,variable) 
    else:
        print('PyECC get error: Variable not found')
        
        
