# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 17:24:37 2016

@author: admin_local
"""

import PyM.PyWIT.PyECC.Base_type as Base_type
import PyM.PyWIT.PyECC.PyECC_globals as PGl

import time


class Fq_type(Base_type.Base_type):
    def __init__(self, value, domain, index = ''):
        self._structure = domain
        self._value = []
        self._index = []
        fIni = False
        auxZ = (0>>self._structure._subDomain)
        if index != '':
            if len(index) == len(value):
                for i in range(len(index)):
                    s = value[i]>>self._structure._subDomain
                    if s != auxZ:
                        self._value.append(s)
                        self._index.append(index[i])
                if len(self._value) == 0:
                    self._value.append(auxZ)
                    self._index.append(0)
            else:
                raise TypeError('Value and index must have the same lenght')
        else:
            if isinstance(value,list):
                for i in range(len(value)):
                    #Falta mirar si cada element de value està al domain
                    if fIni == False:
                        if value[i]==auxZ:
                            continue
                    fIni = True
                    e = value[i]>>self._structure._subDomain
                    if e != auxZ:
                        self._value.append(e)
                        self._index.append(len(value)-i-1)
            elif isinstance(value,int):
                self._value = [value>>self._structure._subDomain()]
                self._index = [0]
            elif isinstance(value, Base_type.Base_type):
                if self._structure.is_element(value):
                    self._value = value._value
                    self._index = self._index
                elif self._structure.is_sub_element(value):
                    self._value = [value>>self._structure._base_field]
                    self._index = [0]
                else:
                    raise TypeError("Parameter Value cannot be converted to this space")
            else:
                raise TypeError("Parameter Value cannot be converted to Fq_type")    
            if len(self._value)==0:
                self._value.append(auxZ)
                self._index = [0]
        (a,b) = self.domain()._reduct_(self._value[:],self._index[:])
        self._value = a
        self._index = b
        self._degree = self._index[0]

    def cp(self):
        return Fq_type(self._value[:],self._structure.cp(), self._index[:])
    
    def is_zero(self):
        if len(self._value) == 1:
            if self._value[1] == (0>>self._structure._sudDomain):
                return True
        return False
    
    def _len_value_(self):
        if len(self._value)>1:
            return len(self._value)
        else:
            if isinstance(self._value[0],int):
                return 1
            return self._value[0]._len_value_()
    
    def _sign_value_(self):
        return 0
        
    def _characteristic_(self):
        return self._structure._characteristic
            
    def _print_(self):
        res = ''
        aux2 = (1>>self._structure._subDomain)
        flaginici = True
        for kk in range(len(self._index)):
            if PGl.POL_DES_SORT:
                i = kk
            else:
                i = len(self._index) - kk - 1
            if self._index[i]>1:
                if isinstance(self._value[i],int):
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    if self._value[i] == 1:
                        res += self._structure._symbol + '**' + self._index[i].__str__() + ' '
                    else:
                        res += self._value[i].__str__() + '*' + self._structure._symbol + '**' + self._index[i].__str__() + ' '
                elif self._value[i] == aux2:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    res += self._structure._symbol + '**' + self._index[i].__str__() + ' '
                else:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    flagPar = False
                    if self._value[i]._len_value_() > 1:
                        res += '('
                        flagPar = True
                        
                    res += self._value[i]._print_()
                    if flagPar:
                        res += ')'
                    res += '*' + self._structure._symbol + '**' + self._index[i].__str__() + ' '
            elif self._index[i] == 1:
                if isinstance(self._value[i],int):
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    if self._value[i] == 1:
                        res += self._structure._symbol + ' '
                    else:
                        res += self._value[i].__str__() + '*' + self._structure._symbol + ' '
                elif self._value[i] == aux2:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    res += self._structure._symbol + ' '
                else:
                    if flaginici:
                        flaginici = False
                    else:
                        res += '+ '
                    flagPar = False
                    if self._value[i]._len_value_() > 1:
                        res += '('
                        flagPar = True
                    res += self._value[i]._print_()
                    if flagPar:
                        res += ')'
                    res += '*' + self._structure._symbol + ' '
            else:
                if flaginici:
                    flaginici = False
                else:
                    res += '+ '
                if isinstance(self._value[i],int):
                    res += self._value[i].__str__() + ' '
                else:
                    res += self._value[i]._print_() + ' '
        if len(res)>0:
            if res[len(res)-1]==' ':
                res=res[:len(res)-1]
        return res       
        
    def __str__(self):
        if PGl.TYPES == 1:
            return self._print_() + ' :: ' + self._structure._print_()
        else:
            return self._print_()
        
    def __repr__(self):
        return self.__str__()
        
    def reverse(self):
        aux = []
        aindex = []
        for i in range(len(self._value)-1,-1,-1):
            aux.append(self._value[i])
            aindex.append(self._degree - self._index[i])
        return Fq_type(aux,self._structure,aindex)
    
    def __hash__(self):
        f = self.transform_to_min_domain()
        if isinstance(f,int):
            return hash(f)
        elif f.domain() != self.domain():
            return hash(f)
        aux = 0
        for i in range(len(self._index)):
            aux = aux + hash(self._value[i]*PGl.HFQT**self._index[i])
            
        return hash(self._structure)*self.degree() + aux
            
    def __neg__(self):
        return self._structure._opposite_(self)
    
    def opposite(self):
        return self._structure._opposite_(self)
    
    def negation(self):
        return self._structure._opposite_(self)    

    def degree(self):
        return self._degree
        
    def leading_coeff(self):
        return self._value[0]
    
    def constant_coeff(self):
        c = self._index[-1]
        if c == 0: return self._value[-1]
        return 0>>self.domain().subdomain()
    
    def leading_term(self):
        return Fq_type(self._value[0],self._structure, self._index[0])
        
    def element_to_list(self):
        auxZero = 0>>self._structure._subDomain
        aux = [auxZero]*(self._degree + 1)
        for i in range(len(self._value)):
            aux[len(aux) - self._index[i] - 1] = self._value[i]
        return aux
    
    def to_list(self):
        auxZero = 0>>self._structure._subDomain
        aux = [auxZero]*(self._degree + 1)
        for i in range(len(self._value)):
            aux[self._index[i]] = self._value[i]
        return aux
    
    def is_invertible(self):
        return self._structure.is_invertible(self)
   
    def is_primitive(self):
        return self._structure.is_primitive(self)
    
    def order(self):
        return self._structure._order_(self)
        
    def min_domain(self):
        if self.degree() == 0:
            return self._value[0].min_domain()
        return self.domain()
    
    def transform_to_min_domain(self):
        if self.degree() == 0:
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
            return Fq_type(aux,D,self._index)
        raise TypeError("Wrong characteristic")
    
    def transform_to_ZZ(self):
        if self.domain().characteristic() == 0:
            return self
        Dom = self.domain().transform_to_ZZ()
        aux = []
        for a in self._value:
            aux.append(a.transform_to_ZZ())
        return Fq_type(aux,Dom,self._index)
    
    def K_(self):
        return self.domain().K_()
        
        