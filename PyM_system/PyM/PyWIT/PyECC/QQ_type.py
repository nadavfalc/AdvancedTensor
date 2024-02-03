# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:43:25 2017

@author: admin_local
"""

import PyM.PyWIT.PyECC.Base_type as Base_type
import PyM.PyWIT.PyECC.PyECC_globals as PGl

import math
from fractions import Fraction


class QQ_type(Base_type.Base_type):
    def __init__(self, domain, numerator, denominator = 1):
        self._structure = domain
        if numerator == 0:
            if denominator != 0:
                self._value = [0>>self.domain().subdomain(), 1>>self.domain().subdomain()]
        if denominator == 0:
            raise TypeError("Denominator value must be non zero")
        self._value = [numerator>>self.domain().subdomain(), denominator>>self.domain().subdomain()]
        self.domain().reduce(self)
    
    def cp(self):
        a = QQ_type(self._structure.cp(), self._value[0], self._value[1])
        a._is_inv = self._is_inv
        return a
    
    def is_zero(self):
        if self._value[0] == 0:
            return True
        return False
        
    def _len_value_(self):
        if self._value[1] != 1:
            return len(self._value)
        else:
            if isinstance(self._value[0],int):
                return 1
            return self._value[0]._len_value_()
    
    def _sign_value_(self):
        if self._value[1] == 1:
            if isinstance(self._value[0],int):
                if self._value[0] >= 0:
                    return 0
                else:
                    return 1
            if isinstance(self.domain().subdomain().characteristic() == 0):
                if (self._len_value_() == 1):
                    lc = self._value[0].leading_coeff()
                    if isinstance(lc,int):
                        if lc >= 0:
                            return 0
                        else:
                            return 1
                    else:
                        return lc._sign_value_()
                else:
                    return 0
            else:
                return 0
        else:
            return 0
    
    def __round__(self, d = None):
        x = self._value[0]/self._value[1]
        return round(x,d)

    def __le__(self,b):
        try:
            if isinstance(b,int):
                a = self._value[0] / self._value[1]
                return a <= b
            if isinstance(b,float):
                a = self._value[0] / self._value[1]
                return a <= b
            if isinstance(b,QQ_type):
                a = self._value[0] / self._value[1]
                bb = b._value[0] / b._value[1]
                return a <= bb
            else:
                return b.__ge__(self)
        except:
            return 'comparison <=: No possible comparison'
    
    def __lt__(self,b):
        try:
            if isinstance(b,int):
                a = self._value[0] / self._value[1]
                return a < b
            if isinstance(b,float):
                a = self._value[0] / self._value[1]
                return a < b
            if isinstance(b,QQ_type):
                a = self._value[0] / self._value[1]
                bb = b._value[0] / b._value[1]
                return a < bb
            else:
                return b.__gt__(self)
        except:
            return 'comparison <: No possible comparison'
    
    def __ge__(self,b):
        try:
            if isinstance(b,int):
                a = self._value[0] / self._value[1]
                return a >= b
            if isinstance(b,float):
                a = self._value[0] / self._value[1]
                return a >= b
            if isinstance(b,QQ_type):
                a = self._value[0] / self._value[1]
                bb = b._value[0] / b._value[1]
                return a >= bb
            else:
                return b.__le__(self)
        except:
            return 'comparison >=: No possible comparison'
            
    def __gt__(self,b):
        try:
            if isinstance(b,int):
                a = self._value[0] / self._value[1]
                return a > b
            if isinstance(b,float):
                a = self._value[0] / self._value[1]
                return a > b
            if isinstance(b,QQ_type):
                a = self._value[0] / self._value[1]
                bb = b._value[0] / b._value[1]
                return a > bb
            else:
                return b.__lt__(self)
        except:
            return 'comparison >: No possible comparison'     
    
    def __float__(self):
        if isinstance(self._value[0],int):
            return self._value[0]/self._value[1]
        else:
            a = self.transform_to_min_domain()
            if isinstance(a._value,int):
                return float(a._value)
            if isinstance(a._value[0],int):
                return a._value[0]/a._value[1]
            raise TypeError("This number cannot be converted into float")
    
    def sqrt(self):
        return math.sqrt(self.__float__())
    
    def __sqrt__(self):
        a = self.__float__()
        return math.sqrt(a)
    
    def __abs__(self):
        a = self._value[0] / self._value[1]
        if a < 0:
            return QQ_type(self.domain(),-self._value[0],self._value[1])
        else:
            return self.cp()
       
    def _print_(self):
        if self._value[1] == 1:
            if isinstance(self._value[0],int):
                return self._value[0].__str__()
            return self._value[0]._print_()
        if isinstance(self._value[0],int):
            if PGl.Q_DEC == -1:
                return self._value[0].__str__() +'/'+self._value[1].__str__()
            else:
                a = self._value[0] / self._value[1]
                res = round(a,PGl.Q_DEC)
                return res.__str__()
        res = ''
        flagP = False
        if self._value[0]._len_value_()>1:
            res += '('
            flagP = True
        res += self._value[0]._print_()
        if flagP:
            res += ')'
            flagP = False
        res += '/'
        if self._value[1]._len_value_()>1:
            res += '('
            flagP = True
        aux = self._value[1]._print_()
        if not flagP and aux.count('*')>0:
            res += '('
            flagP = True
        res += self._value[1]._print_()
        if flagP:
            res += ')'
        return res
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._print_() + ' :: ' + self._structure._print_() 
        else:
            return self._print_()
        
    def __repr__(self):
        return self.__str__()
    
    def __hash__(self):
        return hash(self._value[0]/self._value[1])
        
    def __neg__(self):
        return self._structure._opposite_(self)
    
    def opposite(self):
        return self._structure._opposite_(self)
    
    def negation(self):
        return self._structure._opposite_(self)
        
    def numerator(self):
        return self.domain()._numerator_(self)

    def denominator(self):
        return self.domain()._denominator_(self)
    
    def min_domain(self):
        if self._value[1] == 1:
            if isinstance(self._value[0],int):
                return self.domain().subdomain()
            return self._value[0].min_domain()
        return self.domain()
    
    def transform_to_min_domain(self):
        if self._value[1] == 1:
            if isinstance(self._value[0],int):
                return self._value[0]
            return self._value[0].transform_to_min_domain()
        return self
    
    def change_characteristic(self,n):
        if self.domain().characteristic() == n:
            return self
        if (n == 0):
            return self.transform_to_ZZ()
        cself = self.domain().characteristic()
        if ((cself%n) == 0):
            D = self.domain().change_characteristic(n)
            aux = []
            for a in self._value:
                if isinstance(a,int):
                    A = self.domain().subdomain()
                    aux.append(a>>A)
                else:
                    aux.append(a.change_characteristic(n))
            return QQ_type(D,self._value[0],self._value[1])
        raise TypeError("Wrong characteristic")
            
    def transform_to_ZZ(self):
        if self.domain().characteristic() == 0:
            return self
        Dom = self.domain().transform_to_ZZ()
        aux = []
        for a in self._value:
            aux.append(a.transform_to_ZZ())
        return QQ_type(Dom,self._value[0],self._value[1])
    
    def K_(self):
        return self.domain().K_()
        
    def is_invertible(self):
        return self._structure.is_invertible(self)
            
