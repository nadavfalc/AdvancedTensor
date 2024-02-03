# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:37:55 2017

@author: admin_local
"""

import PyM.PyWIT.PyECC.Base as Base
from PyM.PyWIT.PyECC.QQ_type import *
import PyM.PyWIT.PyECC.PyECC_globals as PGl

from fractions import Fraction


class QQ(Base.Base):
    def __init__(self, base, name = ''):
        self._ring = True
        self._group = True
        self._fraction_structure = True
        self._characteristic = base.characteristic()
        if isinstance(base, QQ):
            self._cardinal = base._cardinal
            self._subDomain = base.subdomain()
            self._base_Domain = base.base_domain()
            self._field = 1
            if name == '':
                self._name = base._name
            else:
                self._name = name
        else:
            if base._cardinal == 'Infinity':
                self._cardinal = 'Infinity'
            else:
                self._cardinal = base._cardinal*(base._cardinal - 1)
            self._subDomain = base
            self._base_Domain = base.base_domain()
            self._field = 1
            if name == '':
                self._name = 'Q['+base._print_()+']'
            else:
                self._name = name
    
    def cp(self):
        return QQ(self._subDomain, self._name)
    
    def characteristic(self):
        return self._characteristic
    
    def cardinal(self):
        return self._cardinal
    
    def element(self, numerator, denominator = 1):
        if self.subdomain().is_field():
            aux = QQ_type(self,numerator,denominator)
            if denominator != 1:
                aux._value[0] = aux._value[0]/aux._value[1]
            return aux._value[0]
        return QQ_type(self,numerator,denominator)
    
    def reduce(self,a):
        if self.subdomain().is_field():
            a._value[0] = a._value[0]/a._value[1]
            a._value[1] = 1>>self.subdomain()
        else:
            try:
                d = self.subdomain().gcd(a._value[0],a._value[1])
            except:
                d = 1
            a._value[0] = a._value[0]//d
            a._value[1] = a._value[1]//d
            if isinstance(a._value[0],float):
                a._value[0] = int(a._value[0])
                a._value[1] = int(a._value[1])
            if self.subdomain().characteristic() == 0:
                if a._value[1] == -1:
                    a._value[0] = -a._value[0]
                    a._value[1] = 1>>self.subdomain()
                elif a._value[0] == 0:
                    a._value[0] = 0>>self.subdomain()
                    a._value[1] = 1>>self.subdomain()
                if isinstance(a._value[0],int):
                    if a._value[1] < 0:
                        a._value[0] = - a._value[0]
                        a._value[1] = - a._value[1]
                elif isinstance(a._value[1], Base_type.Base_type):
                    if a._value[1].domain().is_invertible(a._value[1]):
                        a._value[0] = a._value[0]/a._value[1]
                        a._value[1] = 1>>self.subdomain()
                    if hasattr(a._value[1],'_index'):
                        if isinstance(a._value[1].leading_coeff(),int):
                            if a._value[1].leading_coeff() < 0 :
                                a._value[0] = -a._value[0]
                                a._value[1] = -a._value[1]
                        elif a._value[1] == (-1)>>self.subdomain():
                            a._value[0] = -a._value[0]
                            a._value[1] = 1>>self.subdomain()
        
    def _print_(self):
        return self._name
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._name + ' :: Field'
        return self._name
    
    def __repr__(self):
        return self.__str__()
    
    def is_field(self):
        if self._field == 1:
            return True
        return False
    
    def is_invertible(self, a):
        if isinstance(a,float):
            aux = a>>self
            return self.is_invertible(aux)
        if isinstance(a, int):
            aux = self.element(a)
            return self.is_invertible(aux)
        if isinstance(a, Base_type.Base_type):
            if (self.is_sub_element(a)):
                return not a.is_zero()
        return False
        
    def _inv_(self, a):
        if not self.is_invertible(a>>self):
            raise ZeroDivisionError("This number is not invertible")
        if (self.is_sub_element(a))and(not self.is_element(a)): 
            return self.element(1,a) 
        return self.element(a._value[1],a._value[0])    
        
    def _add_(self, a, b):
        if isinstance(a,float):
             l = Fraction(a).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._add_(aux,b)
        if isinstance(b,float):
             l = Fraction(b).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._add_(a,aux)
#        if (self.subdomain()._cardinal == 'Infinity') and (self.subdomain() == self.base_domain()):
#            if isinstance(a,int):
#                res = a + b._value[0]/b._value[1]
#                res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#                return QQ_type(self,res.numerator,res.denominator)
#            if isinstance(b,int):
#                res = a._value[0]/a._value[1] + b 
#                res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#                return QQ_type(self,res.numerator,res.denominator)
#            res = a._value[0]/a._value[1] + b._value[0]/b._value[1]
#            res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#            return QQ_type(self,res.numerator,res.denominator)
            
        if (self.is_sub_element(a))and(not self.is_element(a)):
            num = a*b._value[1]+b._value[0]
            denum = b._value[1]
            return self.element(num,denum)
        if (self.is_sub_element(b))and(not self.is_element(b)):
            num = b*a._value[1]+a._value[0]
            denum = a._value[1]
            return self.element(num,denum)
        num = a._value[0]*b._value[1]+b._value[0]*a._value[1]
        denum = a._value[1]*b._value[1]
        return self.element(num,denum)            
        
    def _sub_(self, a, b):
        if isinstance(a,float):
             l = Fraction(a).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._sub_(aux,b)
        if isinstance(b,float):
             l = Fraction(b).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._sub_(a,aux)  
#        if (self.subdomain()._cardinal == 'Infinity') and (self.subdomain() == self.base_domain()):
#            if isinstance(a,int):
#                res = a - b._value[0]/b._value[1]
#                res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#                return QQ_type(self,res.numerator,res.denominator)
#            if isinstance(b,int):
#                res = a._value[0]/a._value[1] - b 
#                res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#                return QQ_type(self,res.numerator,res.denominator)
#            res = a._value[0]/a._value[1] - b._value[0]/b._value[1]
#            res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#            return QQ_type(self,res.numerator,res.denominator)
        if (self.is_sub_element(a))and(not self.is_element(a)):
            num = a*b._value[1]-b._value[0]
            denum = b._value[1]
            return self.element(num,denum)
        if (self.is_sub_element(b))and(not self.is_element(b)):
            num = a._value[0]-b*a._value[1]
            denum = a._value[1]
            return self.element(num,denum)
        num = a._value[0]*b._value[1]-b._value[0]*a._value[1]
        denum = a._value[1]*b._value[1]
        return self.element(num,denum)      
        
    def _mul_(self, a, b):
        if isinstance(a,float):
             l = Fraction(a).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._mul_(aux,b)
        if isinstance(b,float):
             l = Fraction(b).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._mul_(a,aux)
#        if (self.subdomain()._cardinal == 'Infinity') and (self.subdomain() == self.base_domain()):
#            if isinstance(a,int):
#                res = a * b._value[0]/b._value[1]
#                res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#                return QQ_type(self,res.numerator,res.denominator)
#            if isinstance(b,int):
#                res = a._value[0]/a._value[1] * b 
#                res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#                return QQ_type(self,res.numerator,res.denominator)
#            res = a._value[0]/a._value[1] * b._value[0]/b._value[1]
#            res = Fraction(res).limit_denominator(PGl.Q_MAXDENOM)
#            return QQ_type(self,res.numerator,res.denominator)
        if (self.is_sub_element(a))and(not self.is_element(a)):
            return self.element(a*b._value[0],b._value[1])
        if (self.is_sub_element(b))and(not self.is_element(b)):#and(not self._subDomain.is_element(b)):
            return self.element(a._value[0]*b,a._value[1])
        return self.element(a._value[0]*b._value[0],a._value[1]*b._value[1])
    
    def _div_(self,a,b):
        """operates two elements using this domain operation /"""
        if isinstance(a,float):
             l = Fraction(a).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._div_(aux,b)
        if isinstance(b,float):
             l = Fraction(b).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._div_(a,aux)  
        return a*self._inv_(b)
    
    def _opposite_(self, a):
        return self.element(-a._value[0],a._value[1])
    
    def __hash__(self):
        return PGl.HQ
    
    def __eq__(self, other):
        if isinstance(other, Base.Base):
            if other._subDomain == self._subDomain:
                if other._fraction_structure:
                    return True
        return False
        
    def __neq__(self, other):
        return not self__eq__(other)
        
    def _eq_(self, a, b):
        if isinstance(a,float):
             l = Fraction(a).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._eq_(aux,b)
        if isinstance(b,float):
             l = Fraction(b).limit_denominator()
             aux = self.element(l.numerator,l.denominator)
             return self._eq_(a,aux)  
        if self.is_element(b):
            if a._value[0] == b._value[0]:
                if a._value[1] == b._value[1]:
                    return True
            return False
        else:
            aux = b>>self
            return self._eq_(a,aux)
    
    def __rrshift__(self, a):
        if isinstance(a,float):
            l = Fraction(a).limit_denominator()
            return self.element(l.numerator>>self.subdomain(),l.denominator>>self.subdomain())
        if (self.is_sub_element(a))and(not self.is_element(a)):
            return self.element(a,1>>self.subdomain())
        if (self.is_element(a)):
            return a
        if isinstance(a, list):
            if len(a)!= 2:
                raise TypeError("a must have 2 elements")
            return self.element(a[0],a[1])
    
    def _numerator_(self,f):
        return f._value[0]
    
    def _denominator_(self,f):
        return f._value[1]
        
    def _rd(self,r=''): # random polynomial of degree r
        A = self.subdomain()
        if A._poly_structure:
            nm = A.rd(r)
            dm = A.rd(r)
            return self.element(nm,dm)
        nm = A.rd()
        dm = A.rd()
        return self.element(nm,dm)
    
    def _rd_nonzero(self,r): # random polynomial of degree r
        A = self.subdomain()
        A = self.subdomain()
        if A._poly_structure:
            nm = A._rd_nonzero(r)
            dm = A._rd_nonzero(r)
            return self.element(nm,dm)
        nm = A._rd_nonzero()
        dm = A._rd_nonzero()
        return self.element(nm,dm)
        
    def change_characteristic(self,n):    
        if self.characteristic() == n:
            return self
        if (n == 0):
            return self.transform_to_ZZ()
        cself = self.characteristic()
        if ((cself%n) == 0):
            SDom = self.subdomain().change_characteristic(n)
            return QQ(SDom)
        raise TypeError("Wrong characteristic")        
    
    def transform_to_ZZ(self):
        if self.characteristic() == 0:
            return self
        SDom = self.subdomain().transform_to_ZZ()
        return QQ(SDom)
    
    def is_sub_element(self, a):
        if isinstance(a,float) and self.subdomain().characteristic() == 0:
            return True
        return super(self.__class__, self).is_sub_element(a)
    
    def K_(self):
        return self.cp()
    