# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:26:55 2017

@author: admin_local
"""

import PyM.PyWIT.PyECC.Base_type as Base_type
import PyM.PyWIT.PyECC.PyECC_globals as PGl


def change_mod(x,n,D=''):
    if isinstance(x,int):
        y = x>>D.subdomain()
    else:
        y = x.change_characteristic(n)
    return y

class Power_Series_type(Base_type.Base_type):
    """Power Serie Element class. """
        
    def __init__(self, domain, function, point = 0, n = 3):
        self._structure = domain
        if type(function) != type(lambda x:x):
            raise TypeError("Function must be of function Type")
        self._value = lambda x: function(x) >> self.domain().subdomain()
        self._point = point>>self.domain().subdomain()
        self._nPrint = n
    
    def cp(self):
        return Power_Series_type(self._structure.cp(), self._value, self._point, self._nPrint)
    
    def is_zero(self,a = '', b = ''):
        if a == '':
            a = -10*self._nPrint + 1
            b = 10*self._nPrint
        elif b == '':
            b = a
            a = -a+1
        aux = self.to_list(a,b)
        for i in range(len(aux)):
            if aux[i] != (0>>self._structure.subdomain()):
                return False
        return True
    
    def to_list(self,a='',b=''):
        if a == '':
            a = 0
            b = self._nPrint
        elif b == '':
            b = a
            a = 0
        aux = []
        for i in range(a,b):
            aux.append(self._value(i))
        return aux
        
    def _len_value_(self):
        return self._nPrint + 1

    def _sign_value_(self):
        return 0
        
    def _print_(self,a='',b=''):
        if a == '':
            a = 0
            b = self._nPrint
        elif b == '':
            b = a
            a = 0
        res = ''
        aux2 = (1>>self._structure._subDomain)
        flaginici = True
        if self._point != 0:
            if isinstance(self._point,int):
                if self._point > 0:
                    symb1 = '(' +  self._structure._symbol + '-' + self._point.__str__() + ')'
                    symb2 = self._structure._symbol + '-' + self._point.__str__()
                else:
                    symb1 = '(' +  self._structure._symbol + '+' + abs(self._point).__str__() + ')'
                    symb2 = self._structure._symbol + '+' + abs(self._point).__str__()
            elif isinstance(self._point,Base_type.Base_type):
                aux = -self._point.min_domain()
                if isinstance(aux,int):
                    if aux < 0:
                        symb1 = '(' +  self._structure._symbol + '-' + abs(aux).__str__() + ')'
                        symb2 = self._structure._symbol + '-' + abs(aux).__str__()
                    else:
                        symb1 = '(' +  self._structure._symbol + '+' + aux.__str__() + ')'
                        symb2 = self._structure._symbol + '+' + aux.__str__()
                else:
                    symb1 = '(' + self._structure._symbol + '+' + aux._print_() + ')'
                    symb2 = self._structure._symbol + '+' + aux._print_()
        else:
            symb1 = self._structure._symbol
            symb2 = self._structure._symbol
        for i in range(a,b):
            aux = self._value(i)
            if aux == 0:
                continue
            if abs(i)>1:
                if isinstance(aux,int):
                    if flaginici:
                        flaginici = False
                    elif aux>0:
                        res += '+ '
                    if i>0:
                        if aux == 1:
                            res += symb1 + '**' + i.__str__() + ' '
                        elif aux == -1:
                            res += '- ' + symb1 + '**' + i.__str__() + ' '
                        elif aux != 0:
                            res += aux.__str__() + '*' + symb1 + '**' + i.__str__() + ' '
                    else:
                        if aux == 1:
                            res += symb1 + '**' + '(' + i.__str__() + ') '
                        elif aux == -1:
                            res += '- ' + symb1 + '**' + '(' + i.__str__() + ') '
                        elif aux != 0:
                            res += aux.__str__() + '*' + symb1 + '**' + '(' + i.__str__() + ') '
                elif aux == aux2:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    if i > 0:
                        res += symb1 + '**' + i.__str__() + ' '
                    else:
                        res += symb1 + '**(' + i.__str__() + ') '
                elif aux != 0:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    flagPar = False
                    if aux._len_value_() > 1:
                        res += '('
                        flagPar = True
                        
                    res += aux._print_()
                    if flagPar:
                        res += ')'
                    if i > 0:
                        res += '*' + symb1 + '**' + i.__str__() + ' '
                    else: 
                        res += '*' + symb1 + '**(' + i.__str__() + ') '
            elif abs(i) == 1:
                if isinstance(aux,int):
                    if flaginici:
                        flaginici = False
                    elif aux>0:
                        res += '+ '
                    if aux == 1:
                        if i == 1:
                            res += symb2
                        else:
                            res += symb1
                    elif aux == -1:
                        res += '- ' + symb1 
                    elif aux != 0:
                        res += aux.__str__() + '*' + symb1
                    if i == -1:
                        res += '**(-1) '
                    else:
                        res += ' '
                elif aux == aux2:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    if i == -1:
                        res += symb1 + '**(-1) '
                    else:
                        res += symb2 + ' '
                else:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    flagPar = False
                    if aux._len_value_() > 1:
                        res += '('
                        flagPar = True
                    res += aux._print_()
                    if flagPar:
                        res += ')'
                    res += '*' + '*' + symb1  
                    if i == -1:
                        res += '**(-1) '
                    else:
                        res += ' '
            else:
                if flaginici:
                    flaginici = False
                elif isinstance(aux,int):
                    if aux>0:
                        res += '+ '
                else:
                    res += '+ '
                if isinstance(aux,int):
                    if aux != 0:
                        res += aux.__str__() + ' '
                else:
                    if aux != 0:
                        res += aux._print_() + ' '
        if len(res)==0:
            res = '0 '
        if abs(b) > 1:
            res += '+ O(' + symb1 + '**'
            if b>0:
                res += b.__str__() +')'
            else:
                res += '(' + b.__str__() +'))'
        elif b == 1:
            res += '+ O(' + symb2 + ')'
        elif b == -1:
            res += '+ O(' + symb1 + '**(-1))'
        else:
            res += '+ O(' + symb1 + '**(0))'
        return res
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._print_() + ' :: PS[' + self._structure.subdomain()._print_() + ']' 
        else:
            return self._print_()
        
    def __repr__(self):
        return self.__str__()
    
    def __getitem__(self,key):
        if isinstance(key,slice):
            a0 = key.start
            if a0 == None:
                a0 = 0
            a1 = key.stop
            if a1 == None:
                a1 = self._nPrint
            a2 = key.step
            if a2 == None:
                a2 = 1
            indexs = list(range(a0,a1,a2))
        if isinstance(key,int):
            indexs = [key]
        if isinstance(key,list):
            indexs = key
        indexs.sort()
        aux = []
        if len(indexs) == 1:
            return self._value(indexs[0])
        for i in indexs:
            aux.append(self._value(i))
        return aux
        
    def __call__(self,n,m = ''):
        return self.domain()._serie_to_poly_(self,n,m)
    
    def __hash__(self):
        return hash(self.domain()._serie_to_poly_(self,0,10))
        
    def __neg__(self):
        return self._structure._opposite_(self)
    
    def opposite(self):
        return self._structure._opposite_(self)
    
    def negation(self):
        return self._structure._opposite_(self)
    
    def is_invertible(self):
        return self._structure.is_invertible(self)
    
    def min_domain(self):
        return self.domain()
    
    def transform_to_min_domain(self):
        D = self.min_domain()
        f = self._value
        return Power_Series_type(D,f,self._point,self._nPrint)

    def change_characteristic(self,n):
        D = self.min_domain()
        f = self._value
        if self.domain().characteristic() == n:
            return Power_Series_type(D,f,self._point,self._nPrint)
        else:
            D = D.change_characteristic(n)
            aux = f(0)
            if isinstance(aux,int):
                g = lambda x: change_mod(f(x),n,D)
            else:
                g = lambda x: f(x).change_characteristic(n)
        if isinstance(self._point,int):
            point = self._point>>D.subdomain()
        else:
            point = self._point.change_characteristic(n)
        return Power_Series_type(D,g,point,self._nPrint)
    
    def der(self):
        return self.domain()._der_(self)
            
    def transform_to_ZZ(self):
        if self.domain().characteristic() == 0:
            return self
        return self.change_characteristic(0)
    
    def K_(self):
        return self.domain().K_()
        
    def __add__(self, other):
        if isinstance(other,Base_type.Base_type):
            if other.domain()._power_structure:
                if not (self._structure.is_sub_domain(other._structure)):
                    if not (other._structure.is_sub_domain(self._structure)):
                        if (self._structure.subdomain().is_sub_domain(other._structure.subdomain())):
                            f = other._value
                            g = lambda x: f(x)>>self._structure.subdomain()
                            point = other._point >> self._structure.subdomain()
                            aux = Power_Series_type(self.domain(),g,point,other._nPrint)
                            return self._structure._add_(self, aux)
                        if (other._structure.subdomain().is_sub_domain(self._structure.subdomain())):
                            f = self._value
                            g = lambda x: f(x)>>other._structure.subdomain()
                            point = self._point >> other._structure.subdomain()
                            aux = Power_Series_type(other.domain(),g,point,self._nPrint)
                            return other._structure._add_(aux, other)
            elif not (self._structure.is_sub_domain(other._structure)):
                if not (other._structure.is_sub_domain(self._structure)):
                    if (other._structure.is_sub_domain(self._structure.subdomain())):
                        D = self.domain()._create_Power_Series_Structure_(other.domain())
                        f = self._value
                        g = lambda x: f(x)>>D.subdomain()
                        point = self._point >> D.subdomain()
                        aux = Power_Series_type(D,g,point,self._nPrint)
                        return aux.domain()._add_(aux,other)
        return Base_type.Base_type.__add__(self,other)
    
    def __sub__(self, other):
        if isinstance(other,Base_type.Base_type):
            if other.domain()._power_structure:
                if not (self._structure.is_sub_domain(other._structure)):
                    if not (other._structure.is_sub_domain(self._structure)):
                        if (self._structure.subdomain().is_sub_domain(other._structure.subdomain())):
                            f = other._value
                            g = lambda x: f(x)>>self._structure.subdomain()
                            point = other._point >> self._structure.subdomain()
                            aux = Power_Series_type(self.domain(),g,point,other._nPrint)  
                            return self._structure._sub_(self, aux)
                        if (other._structure.subdomain().is_sub_domain(self._structure.subdomain())):
                            f = self._value
                            g = lambda x: f(x)>>other._structure.subdomain()
                            point = self._point >> other._structure.subdomain()
                            aux = Power_Series_type(other.domain(),g,point,self._nPrint)
                            return other._structure._sub_(aux, other)
            elif not (self._structure.is_sub_domain(other._structure)):
                if not (other._structure.is_sub_domain(self._structure)):
                    if (other._structure.is_sub_domain(self._structure.subdomain())):
                        D = self.domain()._create_Power_Series_Structure_(other.domain())
                        f = self._value
                        g = lambda x: f(x)>>D.subdomain()
                        point = self._point >> D.subdomain()
                        aux = Power_Series_type(D,g,point,self._nPrint)
                        return aux.domain()._sub_(aux,other)
        return Base_type.Base_type.__sub__(self,other)
    
    def __mul__(self, other):
        if isinstance(other,Base_type.Base_type):
            if other.domain()._power_structure:
                if not (self._structure.is_sub_domain(other._structure)):
                    if not (other._structure.is_sub_domain(self._structure)):
                        if (self._structure.subdomain().is_sub_domain(other._structure.subdomain())):
                            f = other._value
                            g = lambda x: f(x)>>self._structure.subdomain()
                            point = other._point >> self._structure.subdomain()
                            aux = Power_Series_type(self.domain(),g,point,other._nPrint)   
                            return self._structure._mul_(self, aux)
                        if (other._structure.subdomain().is_sub_domain(self._structure.subdomain())):
                            f = self._value
                            g = lambda x: f(x)>>other._structure.subdomain()
                            point = self._point >> other._structure.subdomain()
                            aux = Power_Series_type(other.domain(),g,point,self._nPrint)
                            return other._structure._mul_(aux, other)
            elif not (self._structure.is_sub_domain(other._structure)):
                if not (other._structure.is_sub_domain(self._structure)):
                    if (other._structure.is_sub_domain(self._structure.subdomain())):
                        D = self.domain()._create_Power_Series_Structure_(other.domain())
                        f = self._value
                        g = lambda x: f(x)>>D.subdomain()
                        point = self._point >> D.subdomain()
                        aux = Power_Series_type(D,g,point,self._nPrint)
                        return aux.domain()._mul_(aux,other)
        return Base_type.Base_type.__mul__(self,other)

    def __radd__(self, other):
        if isinstance(other,Base_type.Base_type):
            if other.domain()._power_structure:
                if not (self._structure.is_sub_domain(other._structure)):
                    if not (other._structure.is_sub_domain(self._structure)):
                        if (self._structure.subdomain().is_sub_domain(other._structure.subdomain())):
                            f = other._value
                            g = lambda x: f(x)>>self._structure.subdomain()
                            point = other._point >> self._structure.subdomain()
                            aux = Power_Series_type(self.domain(),g,point,other._nPrint)
                            return self._structure._add_(aux, self)
                        if (other._structure.subdomain().is_sub_domain(self._structure.subdomain())):
                            f = self._value
                            g = lambda x: f(x)>>other._structure.subdomain()
                            point = self._point >> other._structure.subdomain()
                            aux = Power_Series_type(other.domain(),g,point,self._nPrint)
                            return other._structure._add_(other, aux)
            elif not (self._structure.is_sub_domain(other._structure)):
                if not (other._structure.is_sub_domain(self._structure)):
                    if (other._structure.is_sub_domain(self._structure.subdomain())):
                        D = self.domain()._create_Power_Series_Structure_(other.domain())
                        f = self._value
                        g = lambda x: f(x)>>D.subdomain()
                        point = self._point >> D.subdomain()
                        aux = Power_Series_type(D,g,point,self._nPrint)
                        return aux.domain()._add_(other, aux)
        return Base_type.Base_type.__radd__(self,other)
        
    def __rsub__(self, other):
        if isinstance(other,Base_type.Base_type):
            if other.domain()._power_structure:
                if not (self._structure.is_sub_domain(other._structure)):
                    if not (other._structure.is_sub_domain(self._structure)):
                        if (self._structure.subdomain().is_sub_domain(other._structure.subdomain())):
                            f = other._value
                            g = lambda x: f(x)>>self._structure.subdomain()
                            point = other._point >> self._structure.subdomain()
                            aux = Power_Series_type(self.domain(),g,point,other._nPrint)  
                            return self._structure._sub_(aux ,self)
                        if (other._structure.subdomain().is_sub_domain(self._structure.subdomain())):
                            f = self._value
                            g = lambda x: f(x)>>other._structure.subdomain()
                            point = self._point >> other._structure.subdomain()
                            aux = Power_Series_type(other.domain(),g,point,self._nPrint)
                            return other._structure._sub_(other, aux)
            elif not (self._structure.is_sub_domain(other._structure)):
                if not (other._structure.is_sub_domain(self._structure)):
                    if (other._structure.is_sub_domain(self._structure.subdomain())):
                        D = self.domain()._create_Power_Series_Structure_(other.domain())
                        f = self._value
                        g = lambda x: f(x)>>D.subdomain()
                        point = self._point >> D.subdomain()
                        aux = Power_Series_type(D,g,point,self._nPrint)
                        return aux.domain()._sub_(other, aux)
        return Base_type.Base_type.__rsub__(self,other)
    
    def __rmul__(self, other):
        if isinstance(other,Base_type.Base_type):
            if other.domain()._power_structure:
                if not (self._structure.is_sub_domain(other._structure)):
                    if not (other._structure.is_sub_domain(self._structure)):
                        if (self._structure.subdomain().is_sub_domain(other._structure.subdomain())):
                            f = other._value
                            g = lambda x: f(x)>>self._structure.subdomain()
                            point = other._point >> self._structure.subdomain()
                            aux = Power_Series_type(self.domain(),g,point,other._nPrint)   
                            return self._structure._mul_(aux, self)
                        if (other._structure.subdomain().is_sub_domain(self._structure.subdomain())):
                            f = self._value
                            g = lambda x: f(x)>>other._structure.subdomain()
                            point = self._point >> other._structure.subdomain()
                            aux = Power_Series_type(other.domain(),g,point,self._nPrint)
                            return other._structure._mul_(other, aux)
            elif not (self._structure.is_sub_domain(other._structure)):
                if not (other._structure.is_sub_domain(self._structure)):
                    if (other._structure.is_sub_domain(self._structure.subdomain())):
                        D = self.domain()._create_Power_Series_Structure_(other.domain())
                        f = self._value
                        g = lambda x: f(x)>>D.subdomain()
                        point = self._point >> D.subdomain()
                        aux = Power_Series_type(D,g,point,self._nPrint)
                        return aux.domain()._mul_(other, aux)
        return Base_type.Base_type.__rmul__(self,other)