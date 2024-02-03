# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:32:50 2017

@author: admin_local
"""

from PyM.PyWIT.PyECC.Power_Series_type import *
from PyM.PyWIT.PyECC.Poly import *
import PyM.PyWIT.PyECC.PyECC_globals as PGl


def const_zero_funct(x):
    if x == 0:
        return 1
    return 0

class Power_Series(Base.Base):
    def __init__(self, base, symbol = 'x', name = ''):
        self._ring = True
        self._group = True
        self._power_structure = True
        self._characteristic = base.characteristic()
        self._cardinal = 'Infinity'
        self._symbol = symbol
        self._subDomain = base
        self._base_Domain = base.base_domain()
        if name == '':
            self._name = 'PS[' + base._print_() + ']'
        else:
            self._name = name
    
    def cp(self):
        return Power_Series(self._subDomain, self._symbol, self._name)
            
    def characteristic(self):
        return self._characteristic
    
    def cardinal(self):
        return self._cardinal
    
    def element(self, func, point = 0, n = 3):
        D = Power_Series(self.subdomain(),self._symbol,self._name)
        return Power_Series_type(D,func,point,n)
        
    def _print_(self):
        return self._name
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._name + ' :: Ring'
        else:
            return self._name
    
    def __repr__(self):
        return self.__str__()
    
    def is_field(self):
        return False
    
    def _ZZ_(self):
        return 
    
    def is_invertible(self, a):
        return False
        
    def _inv_(self, a):
        raise NotImplementedError
    
    def gcd(self,a,b):
        raise NotImplementedError
    
    def _div_(self, a, b):
        raise NotImplementedError
            
        
    def _add_(self, a, b):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux+b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a+aux
        if (self.is_sub_element(a))and(not self.is_element(a)):
            f = b._value
            g = lambda x: a + f(x)
            return self.element(g, b._point, b._nPrint)
        if (self.is_sub_element(b))and(not self.is_element(b)):
            f = a._value
            g = lambda x: f(x) + b 
            return self.element(g, a._point, a._nPrint)
        f = a._value
        g = b._value
        if a._point != b._point:
            raise TypeError("Points of the series must be the same")
        h = lambda x: f(x) + g(x)
        n = max(a._nPrint,b._nPrint)
        return self.element(h, a._point, n)
        
    def _sub_(self, a, b):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux-b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a-aux
        if (self.is_sub_element(a))and(not self.is_element(a)):
            f = b._value
            g = lambda x: a - f(x)
            return self.element(g, b._point, b._nPrint)
        if (self.is_sub_element(b))and(not self.is_element(b)):
            f = a._value
            g = lambda x: f(x) - b 
            return self.element(g, a._point, a._nPrint)
        f = a._value
        g = b._value
        if a._point != b._point:
            raise TypeError("Points of the series must be the same")
        h = lambda x: f(x) - g(x)
        n = max(a._nPrint,b._nPrint)
        return self.element(h, a._point, n)      
        
    def _mul_(self, a, b):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux*b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a*aux
        if (self.is_sub_element(a))and(not self.is_element(a)):
            f = b._value
            g = lambda x: a*f(x)
            return self.element(g, b._point, b._nPrint)
        if (self.is_sub_element(b))and(not self.is_element(b)):
            f = a._value
            g = lambda x: f(x)*b 
            return self.element(g, a._point, a._nPrint)
        raise NotImplementedError
    
    def _shift_ps_(self,a,d):
        f = a._value
        g = lambda x: f(x+d)
        return self.element(g, a._point, a._nPrint)   
        
    def _der_(self,a):
        f = a._value
        g = lambda x: f(x+1)*(x+1)
        return self.element(g, a._point, a._nPrint)    
    
    def _opposite_(self, a):
        f = a._value
        g = lambda x: -f(x)
        return self.element(g, a._point, a._nPrint)
    
    def __hash__(self):
        return PGl.HPS + hash(self.subdomain())
    
    def __eq__(self, other):
        if isinstance(other, Base.Base):
            if other._subDomain == self._subDomain:
                if other._power_structure:
                    return True
        return False
        
    def __neq__(self, other):
        return not self__eq__(other)
        
    def _eq_(self, a, b, n = ''):
        if isinstance(a,float):
            return False
        if isinstance(b,float):
            return False
        if a.domain() == b.domain():
            if a._point == b._point:
                if n == '':
                    n = 10*max(abs(a._nPrint),abs(b._nPrint))
                for i in range(-n,n+1):
                    if a(i)!= b(i):
                        return False
                return True
        return False
    
    def __rrshift__(self, a):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return self.__rrshift__(aux)
        if (self.is_sub_element(a))and(not self.is_element(a)):
            f = const_zero_funct
            g = lambda x: a*f(x)
            return self.element(g)
        if (self.is_element(a)):
            f = a._value
            point = a._point
            n = a._nPrint
            return self.element(f,point,n)
        lam = lambda x:x
        if type(a) == type(lam):
            f = a
            return self.element(f)
        if isinstance(a,tuple):
            if len(a) == 2:
                if type(a[0]) == type(lam):
                    f = a[0]
                    return self.element(f,a[1])
                else:
                    raise TypeError("Wrong parameter")
            elif len(a) == 3:
                if type(a[0]) == type(lam):
                    f = a[0]
                    return self.element(f,a[1],a[2])
                else:
                    raise TypeError("Wrong parameter")
            else:
                raise TypeError("Wrong number of parameters")
        if isinstance(a,list):
            if len(a) == 2:
                if type(a[0]) == type(lam):
                    f = a[0]
                    return self.element(f,a[1])
                else:
                    raise TypeError("Wrong parameter")
            elif len(a) == 3:
                if type(a[0]) == type(lam):
                    f = a[0]
                    return self.element(f,a[1],a[2])
                else:
                    raise TypeError("Wrong parameter")
            else:
                raise TypeError("Wrong number of parameters")
        else:
            raise TypeError("This element cannot be converted into a Power Series")
    
    def _create_Power_Series_Structure_(self,D, symbol = ''):
        if symbol != '':
            Dom = Power_Series(D,symbol)
        else:
            Dom = Power_Series(D,self._symbol)
        return Dom
    
    
    def _rd(self,r):
        raise NotImplementedError
    
    def _rd_nonzero(self,r): 
        raise NotImplementedError
         
    def change_characteristic(self,n):    
        if self.characteristic() == n:
            return self
        if (n == 0):
            return self.transform_to_ZZ()
        cself = self.characteristic()
        if ((cself%n) == 0):
            SDom = self.subdomain().change_characteristic(n)
            return Power_Series(SDom,self._symbol)
        raise TypeError("Wrong characteristic")        
    
    def transform_to_ZZ(self):
        if self.characteristic() == 0:
            return self
        SDom = self.subdomain().transform_to_ZZ()
        return Power_Series(SDom,self._symbol)
    
    def set_symbol(self, symbol):
        self._symbol = symbol
    
    def get_symbol(self):
        return self._symbol
    
    def K_(self):
        return self.subdomain().cp()
        
    def _serie_to_poly_(self,a,n0,n1 = ''):
        if n1 == '':
            n1 = abs(n0) 
            n0 = 0
        if n1 == 0:
            raise TypeError("range is not valid")
        if n0 < 0:
            n0 = 0
        D = Poly(self.subdomain())
        x = Symbol_type(self._symbol)
        aux = a.to_list(n0,n1)
        aux.reverse()
        M = []
        indexs = [[i] for i in range(n1-1,n0-1,-1)]
        for i in range(len(aux)):
            Mi = Monomial_type(aux[i],[x],indexs[i])
            if Mi._value != 0:
                M.append(Mi)
        if len(M) == 0:
            return Monomial_type(0>>self._subDomain,[x],[0],D)
        elif len(M) == 1.:
            return M[0]
        else:
            return Poly_type(M,[x],D,M[0]._degree,M[0]._degree)