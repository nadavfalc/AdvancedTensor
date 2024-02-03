# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:36:59 2016

@author: admin_local
"""

from PyM.PyWIT.PyECC.arithmetic.ifactor import *
import PyM.PyWIT.PyECC.Base_type as Base_type
import PyM.PyWIT.PyECC.PyECC_globals as PGl


class Zn_type(Base_type.Base_type):
    def __init__(self, value, domain):
        self._structure = domain
        self._value = value % domain._characteristic
    
    def _len_value_(self):
        return 1
    
    def _sign_value_(self):
        return 0
    
    def __repr__(self):
        return self.__str__() 
        
    def __str__(self):
        if PGl.TYPES == 1:
            return self._print_() + " :: " + self._structure._name
        else:
            return self._print_()
        
    def _print_(self):
        return self._value.__str__()
    
    def _inverse(self):
        return self._structure._inv(self)
    
    def inv(self):
        return self._inverse()
        
    def is_zero(self):
        if self._value == 0:
            return True
        return False
        
    def is_invertible(self):
        return self._structure.is_invertible(self)
    
    def is_primitive(self):
        return self._structure.is_primitive(self)
    
    def __hash__(self):
        return hash(self._value)
    
    def __neg__(self):
        return self._structure._opposite_(self)
        
    def opposite(self):
        return self._structure._opposite_(self)
    
    def negation(self):
        return self._structure._opposite_(self)
                    
    def _characteristic_(self):
        return self._structure._characteristic
    
    def order(self):
        return self._structure._order_(self)
        
    def change_characteristic(self,n):
        if self.domain().characteristic() == n:
            return self
        if n == 0:
            return self._value
        if (self.domain().characteristic() % n) == 0:
            D = self.domain().change_characteristic(n)
            return Zn_type(self._value,D)
        raise TypeError("Wrong characteristic")
    
    def transform_to_ZZ(self):
        return self._value
        
    def K_(self):
        return self.domain().K_()
       
    def degree(self):
        return 0    
    
    def cp(self):
        a = Zn_type(self._value,self._structure.cp())
        a._is_inv = self._is_inv
        return a
