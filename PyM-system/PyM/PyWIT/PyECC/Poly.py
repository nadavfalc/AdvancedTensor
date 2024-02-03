# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:21:20 2020

@author: narcis
"""
from PyM.PyWIT.PyECC.QQ import *
from PyM.PyWIT.PyECC.arithmetic.ZZ import * 
import PyM.PyWIT.PyECC.PyECC_globals as PGl
from PyM.PyWIT.PyECC.Zn import *

def _order_symbols_(x,y):
    #TO DO:versió provisional
    if x._value < y._value:
        return -1
    elif x._value > y._value:
        return 1
    else:
        return 0
    #

class Symbol_type():
    def __init__(self, symb):
        self._value = symb
        
    def __str__(self):
        res = self._print_()
        if PGl.TYPES == 1:
            res += ' : symbol'
        return res
        
    def __repr__(self):
        return self.__str__()
    
    def cp(self):
        return Symbol_type(self._value)
    
    def __hash__(self):
        return hash(1)*PGl.HPOLV*PGl.HPOLT
    
    def __eq__(self,x):
        if isinstance(x,Symbol_type):
            if x._value == self._value:
                return True
        elif isinstance(x, Monomial_type):
            return x.__eq__(self)
        elif isinstance(x, Poly_type):
            return False
        return False
    
    def __req__(self, other):
        return self.__eq__(other)
    
    def __ne__(self, other):
        return not self.__eq__(other)

    def __rne__(self, other):
        return self.__ne__(other)
    
    def _sign_value_(self):
        return 0
    
    def is_zero(self):
        return False
    
    def _print_(self):
        return self._value
    
    def __neg__(self):
        return Monomial_type(-1,[self],[1], Poly(ZZ()))
    
    def opposite(self):
        return self.__neg__()
    
    def negation(self):
        return self.__neg__()
    
    def degree(self,symb = None):
        if symb == None:
            return 1
        else:
            x = _get_symbol_(symb)
        if x == self:
            return 1
        else:
            return 0
        
    def mdegree(self):
        return 1
    
    def total_degree(self):
        return 1
    
    def _to_monomial_(self):
        return Monomial_type(1,[self],[1], Poly(ZZ()))
    
    def reverse(self):
        return 1
    
    def reverse_m(self,m,var = None):
        if var == None:
            var = self.cp()
        x = _get_symbol_(var)
        if self == x:
            if m == 0:
                raise TypeError('m must be at least 1')
            return Monomial_type(1, [self],[m - 1],Poly(ZZ()))
        else:
            [variables,pos] = _add_symbol_([self.cp()],x)
            if pos == 0:
                aindex = [m,1]
            else:
                aindex = [1,m]
            return Monomial_type(1, variables,aindex,Poly(ZZ()))
    
    def leading_coeff(self, var = None):
        if var == None:
            Monomial_type(1,[],[],Poly(ZZ()))
        x = _get_symbol_(var)
        if x == self:
            Monomial_type(1,[],[],Poly(ZZ()))
        else:
            Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
    
    def leading_term(self, var = None):
        Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
    
    def trailing_coeff(self,  var = None):
        if var == None:
            Monomial_type(1,[],[],Poly(ZZ()))
        x = _get_symbol_(var)
        if x == self:
            Monomial_type(1,[],[],Poly(ZZ()))
        else:
            Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
    
    def trailing_term(self, var = ''):    
        Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
    
    def constant_coeff(self, var = None):
        if var == None:
            return Monomial_type(0,[],[],Poly(ZZ()))
        else:
            x = _get_symbol_(var)
            if x == self:
                return Monomial_type(0,[],[],Poly(ZZ()))
            else:
                return Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
        
    def constant_term(self, var = None):
        if var == None:
            return Monomial_type(0,[self.cp()],[0],Poly(ZZ()))
        else:
            x = _get_symbol_(var)
            if x == self:
                return Monomial_type(0,[self.cp()],[0],Poly(ZZ()))
            else:
                return Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
    
    def zero_degree_term(self):
        return Monomial_type(0,[self.cp()],[0],Poly(ZZ()))
            
    def trailing_degree(self, var = None):
        if var == None:
            return 1
        else:
            x = _get_symbol_(var)
            if x == self:
                return 1
            else:
                return 0
    
    def der(self,symb = None):
        if symb == None:
            return Monomial_type(1,[self.cp()],[0],Poly(ZZ()))
        else:
            x = _get_symbol_(symb)
            if x == self:
                return Monomial_type(1,[self.cp()],[0],Poly(ZZ()))
            else:
                return Monomial_type(0,[self.cp()],[0],Poly(ZZ()))
    
    def collect(self,var):
        h = self._to_monomial_()
        return h.collect(var)
    
    def eval_symbol(self, value, variables = None):
        if variables == None:
            return value
        else:
            if isinstance(variables, list):
                for i in range(len(variables)):
                    x = _get_symbol_(variables[i])
                    if x == self:
                        return value[i]
                return self.cp()
            x = _get_symbol_(variables)
            if x == self:
                return value
            else:
                return self.cp()
    
    def evaluate(self,value,variables = None):
        return self.eval_symbol(value,variables)
    
    def to_list(self):
        return self.coeffs()
            
    def coeffs(self, var = None):
        if var == None:
            return [Monomial_type(1,[],[],Poly(ZZ())), Monomial_type(0,[],[],Poly(ZZ()))]
        else:
            x = _get_symbol_(var)
            if x == self:
                return [Monomial_type(1,[],[],Poly(ZZ())), Monomial_type(0,[],[],Poly(ZZ()))]
            else:
                return [Monomial_type(1,[self.cp()],[1], Poly(ZZ()))]
    
    def min_domain(self):
        return Poly(ZZ())
    
    def transform_to_min_domain(self):
        return Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
    
    def K_(self):
        return ZZ()
    
    def __add__(self,other):
        if isinstance(other, Symbol_type):
            return Poly(ZZ())._add_(self,other)
        elif isinstance(other,int):
            aux = [Monomial_type(1,[self.cp()],[1], Poly(ZZ())),Monomial_type(other,[self.cp()],[0], Poly(ZZ()))]
            return Poly_type(aux,[self.cp()],Poly(ZZ()),1,1)
        elif isinstance(other,float):
            aux = [Monomial_type(1>>QQ(ZZ()),[self.cp()],[1], Poly(QQ(ZZ()))),Monomial_type(other>>QQ(ZZ()),[self.cp()],[0], Poly(QQ(ZZ())))]
            return Poly_type(aux,[self.cp()],Poly(QQ(ZZ())),1,1)
        elif isinstance(other,Base_type.Base_type):
            s = self._proj_to_domain_(other.domain())
            return s.domain()._add_(s,other)
        try:
            return other.__radd__(self)
        except:
            s = Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
            return s.__add__(other)
        
    def __sub__(self,other):
        if isinstance(other, Symbol_type):
            return Poly(ZZ())._sub_(self,other)
        elif isinstance(other,int):
            aux = [Monomial_type(1,[self.cp()],[1], Poly(ZZ())),Monomial_type(-other,[self.cp()],[0], Poly(ZZ()))]
            return Poly_type(aux,[self.cp()],Poly(ZZ()),1,1)
        elif isinstance(other,float):
            aux = [Monomial_type(1>>QQ(ZZ()),[self.cp()],[1], Poly(QQ(ZZ()))),Monomial_type((-other)>>QQ(ZZ()),[self.cp()],[0], Poly(QQ(ZZ())))]
            return Poly_type(aux,[self.cp()],Poly(QQ(ZZ())),1,1)
        elif isinstance(other,Base_type.Base_type):
            s = self._proj_to_domain_(other.domain())
            return s.domain()._sub_(s,other)
        try:
            return other.__rsub__(self)
        except:
            s = Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
            return s.__sub__(other)
        
    def __mul__(self,other):
        if isinstance(other, Symbol_type):
            return Poly(ZZ())._mul_(self,other)
        elif isinstance(other,int):
            return Monomial_type(other,[self.cp()],[1], Poly(ZZ()))
        elif isinstance(other,float):
            return Monomial_type(other>>QQ(ZZ()),[self.cp()],[1], Poly(QQ(ZZ())))
        elif isinstance(other,Base_type.Base_type):
            s = self._proj_to_domain_(other.domain())
            return s.domain()._mul_(s,other)
        try:
            return other.__rmul__(self)
        except:
            s = Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
            return s.__mul__(other)
        
    def __radd__(self,other):
        if isinstance(other, Symbol_type):
            return Poly(ZZ())._add_(other,self)
        elif isinstance(other,int):
            aux = [Monomial_type(1,[self.cp()],[1], Poly(ZZ())),Monomial_type(other,[self.cp()],[0], Poly(ZZ()))]
            return Poly_type(aux,[self.cp()],Poly(ZZ()),1,1)
        elif isinstance(other,float):
            aux = [Monomial_type(1>>QQ(ZZ()),[self.cp()],[1], Poly(QQ(ZZ()))),Monomial_type(other>>QQ(ZZ()),[self.cp()],[0], Poly(QQ(ZZ())))]
            return Poly_type(aux,[self.cp()],Poly(QQ(ZZ())),1,1)
        elif isinstance(other,Base_type.Base_type):
            s = self._proj_to_domain_(other.domain())
            return s.domain()._add_(other,s)
        try:
            return other.__add__(self)
        except:
            s = Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
            return s.__radd__(other)
        
    def __rsub__(self,other):
        if isinstance(other, Symbol_type):
            return Poly(ZZ())._sub_(other,self)
        elif isinstance(other,int):
            aux = [Monomial_type(-1,[self.cp()],[1], Poly(ZZ())),Monomial_type(other,[self.cp()],[0], Poly(ZZ()))]
            return Poly_type(aux,[self.cp()],Poly(ZZ()),1,1)
        elif isinstance(other,float):
            aux = [Monomial_type((-1)>>QQ(ZZ()),[self.cp()],[1], Poly(QQ(ZZ()))),Monomial_type(other>>QQ(ZZ()),[self.cp()],[0], Poly(QQ(ZZ())))]
            return Poly_type(aux,[self.cp()],Poly(QQ(ZZ())),1,1)
        elif isinstance(other,Base_type.Base_type):
            s = self._proj_to_domain_(other.domain())
            return s.domain()._sub_(other,s)
        try:
            return other.__sub__(self)
        except:
            s = Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
            return s.__rsub__(other)
        
    def __rmul__(self,other):
        if isinstance(other, Symbol_type):
            return Poly(ZZ())._mul_(other,self)
        elif isinstance(other,int):
            return Monomial_type(other,[self.cp()],[1], Poly(ZZ()))
        elif isinstance(other,float):
            return Monomial_type(other>>QQ(ZZ()),[self.cp()],[1], Poly(QQ(ZZ())))
        elif isinstance(other,Base_type.Base_type):
            s = self._proj_to_domain_(other.domain())
            return s.domain()._mul_(other,s)
        try:
            return other.__mul__(self)
        except:
            s = Monomial_type(1,[self.cp()],[1],Poly(ZZ()))
            return s.__rmul__(other)
    
    def __pow__(self,other):
        if isinstance(other, int):
            if other < 0:
                raise TypeError("A symbol cannot be invertible")
            elif other == 0:
                return 1
            else:
                return Monomial_type(1,[self.cp()],[other],Poly(ZZ()))
    
    def _proj_to_domain_(self,D):
        if D._poly_structure:
            return Monomial_type(1>>D.subdomain(), [self.cp()], [1], D)
        else:
            return Monomial_type(1>>D,[self.cp()],[1],Poly(D))
        
    def _push_to_union(self,A):
        return self._proj_to_domain_(A)
        
###Monomial functions and class     
    
def _get_symbol_(var):
    if isinstance(var,str):
        return Symbol_type(var)
    elif isinstance(var,Symbol_type):
        return var
    elif isinstance(var,Monomial_type):
        if len(var._variables) == 1:
            return var._variables[0]
        else:
            ind = 0
            indCount = 0
            for i in range(len(var._index)):
                if var._index[i] > 0:
                    indCount += 1
                    ind = i
            if indCount == 1:
                return var._variables[ind]
            else:
                raise TypeError('var must be a symbol or a univariate monomial or polynomial')
    elif isinstance(var, Poly_type):
        if len(var._variables) == 1:
            return var._variables[0]
        else:
            raise TypeError('var must be a symbol or a univariate monomial or polynomial')
    elif isinstance(var,list):
        if len(var) == 1:
            return _get_symbol_(var[0])
        else:
            raise TypeError('var must be a symbol or a univariate monomial or polynomial')
    else:
        raise TypeError('var must be a symbol or a univariate monomial or polynomial')

def _add_symbol_(L,var):
    x = _get_symbol_(var)
    if x in L:
        return L,L.index(x)
    else:
        for i in range(len(L)):
            if _order_symbols_(x,L[i]) < 0:
                LL = L[:i] + [x] + L[i:]
                return LL,i
        return L+[x],len(L)

# 0 if equal order, -1 if b > a and 1 if a > b
def _order_monomials_(a,b):
    if a._totIndex == 0:
        if b._totIndex == 0:
            return 0
        else:
            return -1
    if b._totIndex == 0:
        return 1
    ##TO DO: versió provisional
    ia = len(a._variables) - 1
    ib = len(b._variables) - 1
    while min(ia,ib) > -1:
        if a._index[ia] == 0:
            ia -= 1
        elif b._index[ib] == 0:
            ib -= 1
        else:
            fab = _order_symbols_(a._variable[ia], b._variable[ib])
            if fab == -1:
                return -1
            if fab == 1:
                return 1
            else:
                if a._index[ia] < b._index[ib]:
                    return -1
                elif a._index[ia] > b._index[ib]:
                    return 1
                else:
                    ia -= 1
                    ib -= 1
    if ia == -1:
        if ib == -1:
            return 0
        else:
            if sum(b._index[:ib+1]) > 0:
                return -1
            else:
                return 0
    else:
        if sum(a._index[:ia+1]) > 0:
            return 1
        else:
            return 0
    ##    

def _order_monomials_same_variables_(a,b):
    if a._totIndex == 0:
        if b._totIndex == 0:
            return 0
        else:
            return -1
    if b._totIndex == 0:
        return 1
    for i in range(len(a._variables)-1,-1,-1):
        if a._index[i] < b._index[i]:
            return -1
        if b._index[i] < a._index[i]:
            return 1
    return 0    

def _find_index_monomials(index,M,indinici):
    pos = len(M)
    fpos = False
    for i in range(indinici,len(M)):
        fend = 0
        for j in range(len(index)-1,-1,-1):
            if index[j] > M[i]._index[j]:
                pos = i
                fend = 1
                break
            elif index[j] < M[i]._index[j]:
                fend = -1
                break
        if fend == 0:
            pos = i
            fpos = True
            break
        if fend == 1:
            break
    return pos,fpos
                   
        

class Monomial_type(Base_type.Base_type):
    def __init__(self, coeff, variables, index, domain, varOrder = None):
        self._structure = domain
        self._value = coeff
        self._varOrder = varOrder
        if len(variables) == 0:
            self._index = []
            self._variables = []
            self._degree = 0
            self._totIndex = 0
        else:
            self._index = index
            self._variables = variables
            self._degree = index[-1]
            self._totIndex = sum(index)
        
    def cp(self):
        if self._varOrder == None:
            return Monomial_type(self._value,self._variables[:], self._index[:], self._structure)
        else:
            return Monomial_type(self._value,self._variables[:], self._index[:], self._structure, self._varOrder[:])
    
    def is_zero(self):
        if self._value == 0:
            return True
        else:
            return False
        
    def _sign_value_(self):
        if (self.domain().characteristic() == 0):
            if isinstance(self._value,int):
                if self._value >= 0:
                    return 0
                else:
                    return 1
            elif len(self._value) > 1:
                return 0
            else:
                return self._value._sign_value_()
        else:
            return 0
        
    def __hash__(self):
        f = self.transform_to_min_domain()
        if not isinstance(f,Monomial_type):
            return hash(f)
        
        aux = hash(self._value)
        for j in range(len(self._variables)):
            aux *= PGl.HPOLV**self._index[j]
        aux *= PGl.HPOLT
        return hash(self._structure)*self.mdegree() + aux   
    
    def _len_value_(self):
        if isinstance(self._value,int):
            return 1
        return self._value._len_value_()
        
    
    def _set_var_order(self, varOrd):
        if len(varOrd) == self._variables:
            return Monomial_type(self._value,self._variables[:],self._index[:],self.domain(),varOrd)
        else:
            raise TypeError("Wrong nuber of variables")
        
    
    def _get_var_order(self):
        if self._varOrder == None:
            return []
        else:
            return self._varOrder[:]

    def _erase_var_order(self):
        return Monomial_type(self._value,self._variables[:],self._index[:],self.domain())
    
    def _print_(self):
        res = ''
        
        if len(self._variables) == 0:
            if isinstance(self._value, int):
                return self._value.__str__()
            else:
                return self._value._print_()
        if self._varOrder == None:
            if isinstance(self._value, int):
                if self._totIndex == 0:
                    return self._value.__str__()
                elif self._value == 1:
                    for i in range(len(self._variables)):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
                elif self._value == -1:
                    res = '-'
                    for i in range(len(self._variables)):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
                else:
                    res = self._value.__str__() + '*'
                    for i in range(len(self._variables)):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
            else:
                if self._totIndex == 0:
                    return self._value._print_()
                elif self._value == 1:
                    for i in range(len(self._variables)):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
                elif self._value == -1 and self._structure._characteristic == 0:
                    res = '-'
                    for i in range(len(self._variables)):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
                else:
                    if self._value._len_value_() > 1:
                        res = '(' + self._value._print_() + ')*'
                    else:
                        res = self._value._print_() + '*'
                    for i in range(len(self._variables)):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
        else:
            if isinstance(self._value, int):
                if self._totIndex == 0:
                    return self._value.__str__()
                elif self._value == 1:
                    for i in range(len(self._variables)):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
                elif self._value == -1:
                    res = '-'
                    for i in range(len(self._variables)):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
                else:
                    res = self._value.__str__() + '*'
                    for i in range(len(self._variables)):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
            else:
                if self._totIndex == 0:
                    return self._value._print_()
                elif self._value == 1:
                    for i in range(len(self._variables)):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
                elif self._value == -1 and self._structure._characteristic == 0:
                    res = '-'
                    for i in range(len(self._variables)):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
                else:
                    if self._value._len_value_() > 1:
                        res = '(' + self._value._print_() + ')*'
                    else:
                        res = self._value._print_() + '*'
                    for i in range(len(self._variables)):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
    
    def _print_coeff_(self):
        res = ''
        
        if len(self._variables) == 0:
            if isinstance(self._value, int):
                return self._value.__str__()
            else:
                return self._value._print_()
        if self._varOrder == None:
            if isinstance(self._value, int):
                if self._totIndex == 0:
                    return self._value.__str__()
                elif self._value == 1:
                    for i in range(len(self._variables)-1):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
                elif self._value == -1:
                    res = '-'
                    for i in range(len(self._variables)-1):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    if len(res) == 1:
                        return res
                    return res[:-1]
                else:
                    res = self._value.__str__() + '*'
                    for i in range(len(self._variables)-1):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
            else:
                if self._totIndex == 0:
                    return self._value._print_()
                elif self._value == 1:
                    for i in range(len(self._variables)-1):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
                elif self._value == -1 and self._structure._characteristic == 0:
                    res = '-'
                    for i in range(len(self._variables)-1):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    if len(res) == 1:
                        return res
                    return res[:-1]
                else:
                    if self._value._len_value_() == 1:
                        res = self._value._print_() + '*'
                    else:
                        if self._totIndex - self._index[-1] > 0:
                            res = '(' + self._value._print_() + ')*'
                        else:
                            res = self._value._print_() + '*'
                    for i in range(len(self._variables)-1):
                        if self._index[i] > 1:
                            res += self._variables[i]._print_() + '**' + self._index[i].__str__() + '*'
                        elif self._index[i] == 1:
                            res += self._variables[i]._print_() + '*'
                    return res[:-1]
        else:
            if isinstance(self._value, int):
                if self._totIndex == 0:
                    return self._value.__str__()
                elif self._value == 1:
                    for i in range(len(self._variables)-1):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
                elif self._value == -1:
                    res = '-'
                    for i in range(len(self._variables)-1):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    if len(res) == 1:
                        return res
                    return res[:-1]
                else:
                    res = self._value.__str__() + '*'
                    for i in range(len(self._variables)-1):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
            else:
                if self._totIndex == 0:
                    return self._value._print_()
                elif self._value == 1:
                    for i in range(len(self._variables)-1):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
                elif self._value == -1 and self._structure._characteristic == 0:
                    res = '-'
                    for i in range(len(self._variables)-1):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    if len(res) == 1:
                        return res
                    return res[:-1]
                else:
                    if self._value._len_value_() == 1:
                        res = self._value._print_() + '*'
                    else:
                        if self._totIndex - self._index[self._varOrder[-1]] > 0:
                            res = '(' + self._value._print_() + ')*'
                        else:
                            res = self._value._print_() + '*'
                    for i in range(len(self._variables)-1):
                        if self._index[self._varOrder[i]] > 1:
                            res += self._variables[self._varOrder[i]]._print_() + '**' + self._index[self._varOrder[i]].__str__() + '*'
                        elif self._index[self._varOrder[i]] == 1:
                            res += self._variables[self._varOrder[i]]._print_() + '*'
                    return res[:-1]
        
    def __str__(self):
        if PGl.TYPES == 1:
            if len(self._variables) == 0:
                return self._print_() + ' :: ' + self._structure._print_()                    
            if self.domain()._printdefault:
                res = self._print_() + ' :: ' + self._structure.subdomain()._print_()
                if self._varOrder == None:
                    for i in range(len(self._variables)):
                        res += '[' + self._variables[i]._print_() + ']'
                    return res
                else:
                    for i in range(len(self._variables)):
                        res += '[' + self._variables[self._varOrder[i]]._print_() + ']'
                    return res
            else:
                return self._print_() + ' :: ' + self._structure._print_()
        else:
            return self._print_()
        
    def __repr__(self):
        return self.__str__()
    
    def reverse(self,var=None):
        if var == None:
            if len(self._variables) == 0:
                return self.cp()
            else:
                if self._totIndex == 0:
                    return self.cp()
                if self._varOrder == None:
                    for i in range(len(self._index),-1,-1):
                        ind = 0
                        if self._index[i] > 0:
                            ind = i + 1
                    return Monomial_type(self._value,self._variables[:ind],self._index[:ind],self._structure)
                else:
                    fbo = [1]*len(self._variables)
                    for i in range(len(self._index),-1,-1):
                        ind = 0
                        if self._index[self._varOrder[i]] > 0:
                            ind = i + 1
                        else:
                            fbo[self._varOrder[i]] = 0
                    aindex = []
                    avars = []
                    for i in range(len(self._variables)):
                        if fbo[i] == 1:
                            aindex.append(self._index[i])
                            avars.append(self._variables[i])
                    return Monomial_type(self._value, avars, aindex, self._structure, self._varOrder[:])
        else:
            x = _get_symbol_(var)
            m = self.degree(x)
            return self.reverse_m(m,x)
    
    def reverse_m(self,m,var = None):
        if var == None:
            x = self.variable()
            if x == []:
                raise TypeError('This polynomial has no variable')
            return self.reverse_m(m,x)
        else:
            x = _get_symbol_(var)
            if len(self._variables) == 0:
                return Monomial_type(self._value,[x],[m],self._structure)
            if m < self.degree(x):
                raise TypeError('m must be equal or larger than the degree of the monomial')
            if x in self._variables:
                aindex = self._index[:]
                pos = aindex.index(x)
                aindex[pos] = m - aindex[pos]
                if self._varOrder == None:
                    return Monomial_type(self._value, self._variables[:],aindex,self._structure)
                else:
                    return Monomial_type(self._value, self._variables[:],aindex,self._structure,self._varOrder[:])
            else:
                [variables, pos] = _add_symbol_(self._variables[:],x)
                aindex = self._index[:pos]+[m]+self._index[pos:]
                return Monomial_type(self._value, variables,aindex,self._structure)
            
    def trunc(self,k,var=None):
        if k == 0:
            return Monomial_type(0>>self._structure.subdomain(), self._variables[:],[0]*len(self._variables),self._structure)
        if self._totIndex == 0:
            return self.cp()
        if var == None:
            if self.degree() < k:
                return self.cp()
            else:
                return Monomial_type(0>>self._structure.subdomain(), self._variables[:],[0]*len(self._variables),self._structure)
        else:
            x = _get_symbol_(var)
            if x in self._variables:
                if self.degree(x) < k:
                    return self.cp()
                else:
                    return Monomial_type(0>>self._structure.subdomain(), self._variables[:],[0]*len(self._variables),self._structure)   
            else:
                [avariables,pos] = _add_symbol_(self._variables[:],x)
                if k > 0:
                    return Monomial_type(0>>self._structure.subdomain(), avariables,[0]*len(avariables),self._structure) 
                else:
                    return Monomial_type(self._value,avariables, self._index[:pos]+[0]+self._index[pos:],self._structure)
                
    def trunc_total_degree(self, d):
        if self._totIndex < d:
            return self.cp()
        else:
            return Monomial_type(0>>self._structure.subdomain(), self._variables,[0]*len(self._variables),self._structure)
    
    def get_subpolynomial_degree(self,d, var = None):
        return self.trunc(d + 1 ,var) - self.trunc(d,var)
    
    def get_subpolynomial_total_degree(self,d):
        return self.trunc_total_degree(d + 1) - self.trunc_total_degree(d)

    def truncated_inverse(self,k):
        return self._structure._trunc_inv_2_(self,k)       
        
    def __neg__(self):
        if self._varOrder == None:
            return Monomial_type(-self._value, self._variables[:],self._index[:],self._structure)
        else:
            return Monomial_type(-self._value, self._variables[:],self._index[:],self._structure,self._varOrder[:])
    
    def opposite(self):
        return self.__neg__()
    
    def negation(self):
        return self.__neg__()
    
    def degree(self, symb = None):
        if symb == None:
            if len(self._variables) == 0:
                return 0
            if self._varOrder == None:
                return self._index[-1]
            else:
                return self._index[self._varOrder[-1]]
        else:
            x = _get_symbol_(symb)
        if x in self._variables:
            ind = self._variables.index(x)
            return self._index[ind]
        else:
            return 0
    
    def mdegree(self):
        return self._totIndex
    
    def total_degree(self):
        return self._totIndex
    
    def leading_coeff(self, var = None):
        if var == None:
            if self._varOrder == None:
                return Monomial_type(self._value, self._variables[:-1], self._index[:-1], self._structure)
            else:
                varOrder = []
                v = self._varOrder[-1]
                for i in range(len(self._varOrder) - 1):
                    if self._varOrder[i] > v:
                        varOrder.append(self._varOrder[i] - 1)
                    else:
                        varOrder.append(self._varOrder[i])
                return Monomial_type(self._value, self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1]+1:], self._index[:self._varOrder[-1]] + self._index[self._varOrder[-1] + 1:], self._structure,varOrder)
        x = _get_symbol_(var)
        if x in self._variables:
            ind = self._variables.index(x)
            if self._varOrder == None:
                return Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure)
            else:
                varOrder = []
                for i in range(len(self._varOrder)):
                    if self._varOrder[i] > ind:
                        varOrder.append(self._varOrder[i] - 1)
                    elif self._varOrder[i] < ind:
                        varOrder.append(self._varOrder[i])
                return Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure, varOrder)  
        else:
            if self._varOrder == None:
                return self.cp()
            else:
                return Monomial_type(self._value, self._variables[:], self._index[:], self._structure)
        
    
    def leading_term(self, var = None):
        if var == None:
            return self.cp()
        else:
            x = _get_symbol_(var)
            if x in self._variables:
                return self.cp()
            else:
                [avariables, pos] = _add_symbol_(self._variables[:], x)
                return Monomial_type(self._value, avariables, self._index[:pos] + [0] + self._index[pos:], self._structure)
    
    def trailing_coeff(self,  var = None):
        if var == None:
            if self._varOrder == None:
                return Monomial_type(self._value, self._variables[:-1], self._index[:-1], self._structure)
            else:
                varOrder = []
                v = self._varOrder[-1]
                for i in range(len(self._varOrder) - 1):
                    if self._varOrder[i] > v:
                        varOrder.append(self._varOrder[i] - 1)
                    else:
                        varOrder.append(self._varOrder[i])
                return Monomial_type(self._value, self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1]+1:], self._index[:self._varOrder[-1]] + self._index[self._varOrder[-1] + 1:], self._structure, varOrder)
        x = _get_symbol_(var)
        if x in self._variables:
            ind = self._variables.index(x)
            if self._varOrder == None:
                return Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure)
            else:
                varOrder = []
                for i in range(len(self._varOrder)):
                    if self._varOrder[i] > ind:
                        varOrder.append(self._varOrder[i] - 1)
                    elif self._varOrder[i] < ind:
                        varOrder.append(self._varOrder[i])
                return Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure, varOrder)  
        else:
            if self._varOrder == None:
                return self.cp()
            else:
                return Monomial_type(self._value, self._variables[:], self._index[:], self._structure)
    
    def trailing_term(self, var = None):    
        if var == None:
            return self.cp()
        else:
            x = _get_symbol_(var)
            if x in self._variables:
                return self.cp()
            else:
                [avariables, pos] = _add_symbol_(self._variables[:], x)
                return Monomial_type(self._value, avariables, self._index[:pos] + [0] + self._index[pos:], self._structure)

    
    def constant_coeff(self, var = None):
        if var == None:
            if len(self._variables) == 0:
                return Monomial_type(self._value, self._variables[:], self._index[:], self._structure)
            if self._varOrder == None:
                if self._index[-1] == 0:
                    return Monomial_type(self._value, self._variables[:-1], self._index[:-1], self._structure)
                else:
                    return Monomial_type(0 >> self._structure.subdomain(), self._variables[:-1],[0]*(len(self._variables)-1), self._structure)
            else:
                if self._index[self._varOrder[-1]] == 0:
                    varOrder = []
                    v = self._varOrder[-1]
                    for i in range(len(self._varOrder) - 1):
                        if self._varOrder[i] > v:
                            varOrder.append(self._varOrder[i] - 1)
                        else:
                            varOrder.append(self._varOrder[i])
                    return Monomial_type(self._value, self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1]+1:], self._index[:self._varOrder[-1]] + self._index[self._varOrder[-1] + 1:], self._structure, varOrder)
                else:
                    return Monomial_type(0 >> self._structure.subdomain(), self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1]+1:],[0]*(len(self._variables)-1), self._structure)
        else:
            x = _get_symbol_(var)
            if x in self._variables:
                ind = self._variables.index(x)
                if self._index[ind] == 0:
                    if self._varOrder == None:
                        return Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure)
                    else:
                        varOrder = []
                        for i in range(len(self._varOrder)):
                            if self._varOrder[i] > ind:
                                varOrder.append(self._varOrder[i] - 1)
                            elif self._varOrder[i] < ind:
                                varOrder.append(self._varOrder[i])
                        return Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure, varOrder)  
                else:
                    return Monomial_type(0 >> self._structure.subdomain(), self._variables[:ind]+self._variables[ind+1:], [0]*(len(self._variables)-1), self._structure)
            else:
                return Monomial_type(self._value, self._variables[:], self._index[:], self._structure)
        
    def constant_term(self, var = None):
        if var == None:
            if len(self._variables) == 0:
                return Monomial_type(self._value, self._variables[:], self._index[:], self._structure)
            if self._varOrder == None:
                if self._index[-1] == 0:
                    return Monomial_type(self._value, self._variables[:], self._index[:], self._structure)
                else:
                    return Monomial_type(0>>self._structure.subdomain(), [0]*len(self._variables), self._index[:], self._structure)
            else:
                if self._index[self._varOrder[-1]] == 0:
                    return self.cp()
                else:
                    return Monomial_type(0>>self._structure.subdomain(), [0]*len(self._variables), self._index[:], self._structure)

        else:
            x = _get_symbol_(var)
            if x in self._variables:
                ind = self._variables.index(x)
                if self._index[ind] == 0:
                    return self.cp()
                else:
                    return Monomial_type(0>>self._structure.subdomain(), self._variables[:], [0]*len(self._variables), self._structure)
            else:
                [avariables, pos] = _add_symbol_(self._variables[:], x)
                return Monomial_type(self._value, avariables, self._index[:pos] + [0] + self._index[pos:], self._structure)

        
    def trailing_degree(self, var = None):
        if var == None:
            if len(self._variables) == 0:
                return 0
            if self._varOrder == None:
                return self._index[-1]
            else:
                return self._index[self._varOrder[-1]]
        else:
            x = _get_symbol_(var)
            if x in self._variables:
                ind = self._variables.index(x)
                return self._index[ind]
            else:
                return 0
            
    def zero_degree_term(self):
        if self._totIndex == 0:
            return self.cp()
        else:
            return Monomial_type(0>>self.domain().subdomain(),self._variables[:],[0]*len(self._variables),self.domain())
    
    def der(self,symb = None):
        if symb == None:
            if len(self._variables) == 0:
                return Monomial_type(0>>self._structure.subdomain(), self._variables[:], [0]*len(self._variables), self._structure)
            if self._varOrder == None:
                if self._index[-1] == 0:
                    return Monomial_type(0>>self._structure.subdomain(), self._variables[:], [0]*len(self._variables), self._structure)
                else:
                    return Monomial_type(self._index[-1]*self._value, self._variables[:], self._index[:-1]+[self._index[-1]-1], self._structure)
            else:
                if self._index[self._varOrder[-1]] == 0:
                    return Monomial_type(0>>self._structure.subdomain(), self._variables[:], [0]*len(self._variables), self._structure)
                else:
                    return Monomial_type(self._index[self._varOrder[-1]]*self._value, self._variables[:], self._index[:self._varOrder[-1]]+[self._index[self._varOrder[-1]]-1] + self._index[self._varOrder[-1]:], self._structure, self._varOrder[:])
        else:
            x = _get_symbol_(symb)
            if x in self._variables:
                ind = self._variables.index(x)
                if self._index[ind] == 0:
                    return Monomial_type(0>>self._structure.subdomain(), self._variables[:], [0]*len(self._variables), self._structure)
                else:
                    if self._varOrder == None:
                        return Monomial_type(self._index[ind]*self._value, self._variables[:], self._index[:ind]+[self._index[ind]-1]+self._index[ind+1:], self._structure)
                    else:
                        return Monomial_type(self._index[ind]*self._value, self._variables[:], self._index[:ind]+[self._index[ind]-1]+self._index[ind+1:], self._structure, self._varOrder[:])
            else:
                [avariables, pos] = _add_symbol_(self._variables[:], x)
                return Monomial_type(0>>self._structure.subdomain(), avariables, [0]*len(avariables), self._structure)
            
    def __pow__(self,other):
        if isinstance(other, int):
            if other < 0:
                raise TypeError("A monomial cannot be invertible")
            elif other == 0:
                return Monomial_type(1>>self.domain().subdomain(),self._variables[:],[0],self.domain())
            else:
                if self._varOrder == None:
                    return Monomial_type(self._value**other,self._variables[:],[k*other for k in self._index],self.domain())
                else:
                    return Monomial_type(self._value**other,self._variables[:],[k*other for k in self._index],self.domain(), self._varOrder[:])
                    
    def to_dict(self):
        T = {}
        if len(self._variables) == 0:
            T[0] = self._value[0]
            return T
        T[tuple(self._index)] = self._value
        return T            
        
    def skeleton(self):
        if len(self._variables) == 0:
            return self._value, self._totIndex
        return self._value, self._index    
    
    def eval_monomial(self, value, variables = None):
        if len(self._variables) == 0:
            return self.cp()
        if variables == None:
            if self._varOrder == None:
                return value**self._index[-1]*Monomial_type(self._value, self._variables[:-1], self._index[:-1], self._structure)
            else:
                varOrder = []
                v = self._varOrder[-1]
                for i in range(len(self._varOrder) - 1):
                    if self._varOrder[i] > v:
                        varOrder.append(self._varOrder[i] - 1)
                    else:
                        varOrder.append(self._varOrder[i])
                return value**self._index[self._varOrder[-1]]*Monomial_type(self._value, self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1] + 1:], self._index[:self._varOrder[-1]] + self._index[self._varOrder[-1] + 1:], self._structure, varOrder)
        else:
            if isinstance(variables, list):
                res = 1>>self._structure.subdomain()
                fvar = [0]*len(self._variables)
                for i in range(len(variables)):
                    x = _get_symbol_(variables[i])
                    if x in self._variables:
                        ind = self._variables.index(x)
                        fvar[ind] = 1
                        res *= value[i]**self._index[ind]
                aindex = []
                avariables = []
                varOrder = []
                for i in range(len(self._variables)):
                    if fvar[i] == 0:
                        aindex.append(self._index[i])
                        avariables.append(self._variables[i])
                    if self._varOrder != None:
                        if fvar[self._varOrder[i]] != 1:
                            varOrder.append(self._varOrder[i] - sum(fvar[:self._varOrder[i]]))
                if len(aindex) == 0:
                    return Monomial_type(res*self._value,[],[], self._structure)
                else:
                    if self._varOrder == None:
                        return res*Monomial_type(self._value, avariables, aindex, self._structure)
                    else:
                        return res*Monomial_type(self._value, avariables, aindex, self._structure, varOrder)
            else:    
                x = _get_symbol_(variables)
                if x in self._variables:
                    ind = self._variables.index(x)
                    if self._varOrder == None:
                        return value**self._index[ind]*Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure)
                    else:
                        varOrder = []
                        for i in range(len(self._varOrder)):
                            if self._varOrder[i] > ind:
                                varOrder.append(self._varOrder[i] - 1)
                            elif self._varOrder[i] < ind:
                                varOrder.append(self._varOrder[i])
                        return value**self._index[ind]*Monomial_type(self._value, self._variables[:ind]+self._variables[ind+1:], self._index[:ind]+self._index[ind+1:], self._structure, varOrder)
                else:
                    return self.cp()
    
    def evaluate(self,value,variables = None):
        return self.eval_monomial(value,variables)
    
    def coeffs(self, var = None):
        if len(self._variables) == 0:
            return [Monomial_type(self._value, [], [], self._structure)]
        d = self.degree(var)
        if var == None:
            if self._varOrder == None:
                cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:-1], [0]*(len(self._variables)-1), self._structure)]*d
            else:
                cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1] + 1], [0]*(len(self._variables)-1), self._structure, self._varOrder[:-1])]*d
            cf = [self.leading_coeff()]+cf
            return cf
        else:
            x = _get_symbol_(var)
            if x in self._variables:
                ind = self._variables.index(x)
                if self._varOrder == None:
                    cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:ind] + self._variables[ind+1:], [0]*(len(self._variables)-1), self._structure)]*d
                else:
                    varOrder = []
                    for i in range(len(self._varOrder)):
                        if self._varOrder[i] > ind:
                            varOrder.append(self._varOrder[i] - 1)
                        elif self._varOrder[i] < ind:
                            varOrder.append(self._varOrder[i])
                    cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:ind] + self._variables[ind+1:], [0]*(len(self._variables)-1), self._structure, varOrder)]*d
                cf = [self.leading_coeff(x)]+cf
                return cf
            else:
                if self._varOrder == None:
                    return [self.cp()]
                else:
                    return [Monomial_type(self._value,self._variables[:], self._index[:],self._structure)]
    
    def to_list(self):
        return self.coeffs()
    
    def min_domain(self):
        if self._totIndex == 0:
            if isinstance(self._value, int):
                return ZZ()
            else:
                return self._value.min_domain()
        else:
            return self._structure
        
    def transform_to_min_domain(self):
        if self._totIndex == 0:
            if isinstance(self._value, int):
                return self._value
            else:
                return self._value.transform_to_min_domain()
        else:
            avariables = []
            aindex = []
            fvar = [0]*len(self._variables)
            for i in range(len(self._variables)):
                if self._index[i] != 0:
                    avariables.append(self._variables[i])
                    aindex.append(self._index[i])
                else:
                    if self._varOrder != None:
                        fvar[i] = 1
            if self._varOrder == None:
                return Monomial_type(self._value,avariables,aindex,self._structure)
            else:
                varOrder = []
                for i in range(len(self._variables)):
                    if fvar[self._varOrder[i]] != 1:
                        varOrder.append(self._varOrder[i] - sum(fvar[:self._varOrder[i]]))
                return Monomial_type(self._value,avariables,aindex,self._structure, varOrder)
    
    def change_characteristic(self,n):
        if n == self._structure.characteristic():
            return self.cp()
        cself = self.domain().characteristic()
        if (cself%n) == 0:
            D = self._structure.change_characteristic(n)
            if isinstance(self._value, int):
                if self._varOrder == None:
                    s = self._value >> D.subdomain()
                    if s != 0:
                        return Monomial_type(s, self._variables[:], self._index[:], D)
                    else:
                        return Monomial_type(s, self._variables[:], [0]*len(self._variables), D)
                else:
                    s = self._value >> D.subdomain()
                    if s != 0:
                        return Monomial_type(s, self._variables[:], self._index[:], D, self._varOrder[:])
                    else:
                        return Monomial_type(s, self._variables[:], [0]*len(self._variables), D)
            else:
                if self._varOrder == None:
                    s = self._value.change_characteristic(n)
                    if s != 0:
                        return Monomial_type(s, self._variables[:], self._index[:], D)
                    else:
                        return Monomial_type(s, self._variables[:], [0]*len(self._variables), D)
                else:
                    s = self._value.change_characteristic(n)
                    if s != 0:
                        return Monomial_type(s, self._variables[:], self._index[:], D, self._varOrder[:])
                    else:
                        return Monomial_type(s, self._variables[:], [0]*len(self._variables), D)
        raise TypeError('Wrong characteristic') 
        
    def transform_to_ZZ(self):
        if self.domain().characteristic() == 0:
            return self.cp()
        else:
            D = self.domain().transform_to_ZZ()
            if self._varOrder == None:
                return Monomial_type(self._value.transform_to_ZZ(), self._variables[:], self._index[:], D)
            else:
                return Monomial_type(self._value.transform_to_ZZ(), self._variables[:], self._index[:], D, self._varOrder[:])
    
    def is_invertible(self):
        return self._structure.is_invertible(self)
    
    def K_(self):
        return self.domain().K_()
               
    def variable(self):
        if len(self._variables) == 0:
            return self._variables
        if self._varOrder == None:
            return self._variables[-1]
        else:
            return self._variables[self._varOrder[-1]]
    
    def variables(self):
        return self._variables[:]
    
    def nvariables(self):
        return len(self._variables)
    
    def add_variable(self,v):
        [avariables, pos] = _add_symbol_(self._variables[:],v)
        if len(avariables) == len(self._variables):
            return self.cp()
        else:
            return Monomial_type(self._value, avariables, self._index[:pos]+[0]+self._index[pos:], self._structure)
    
    def add_variables(self,V):
        p = self.cp()
        for v in V:
            p = p.add_variable(v)
        return p
        
    def variables_name(self):
        return [i._value for i in self._variables]
    
    def collect(self,var):
        x = _get_symbol_(var)
        if x == self.variable():
            return self.cp()
        if x == self._variables[-1]:
            return Monomial_type(self._value,self._variables[:], self._index[:],self._structure)
        if x in self._variables:
            pos = self._variables.index(x)
            varOrd = list(range(len(self._variables)))
            varOrd = varOrd[:pos] + varOrd[pos+1:] + [varOrd[pos]]
            return Monomial_type(self._value,self._variables[:], self._index[:],self._structure,varOrd)
        else:
            [variables,pos] = _add_symbol_(self._variables[:],x)
            aindex = self._index[:pos] + [0] + self._index[pos:]
            varOrd = list(range(len(variables)))
            varOrd = varOrd[:pos] + varOrd[pos+1:] + [varOrd[pos]]
            return Monomial_type(self._value,variables, aindex,self._structure,varOrd)
    
    def __floordiv__(self, other):
        if isinstance(other, Symbol_type):
            return self.__floordiv__(Monomial_type(1>>self.K_(),[other],[1], self._structure))
        if isinstance(other,Base_type.Base_type):
            if other.domain()._poly_structure:
                if self.domain() == other.domain():
                    [d,r] = self._structure._div_res_(self,other)
                    return d
                if self.domain().subdomain().is_sub_domain(other.domain().subdomain()):
                    [d,r] = self._structure._div_res(self,other >> self._structure)
                    return d
                if other.domain().subdomain().is_sub_domain(self.domain().subdomain()):
                    [d,r] = other._structure._div_res(self >> other._structure,other)
                    return d 
                else:
                    try:
                        D = (self.K_().union_domain(other.K_()))[0]
                        PD = Poly(D)
                        [d,r] = PD._div_res_(self >> PD,other >> PD)
                        return r
                    except:
                        try:
                            return other.__rtruediv__(self)
                        except:
                            raise TypeError("I do not know how to opperate this elements")
            elif self.domain().subdomain().is_sub_domain(other.domain()):
                aux = Monomial_type(other >> self.K_(), self._variables[:],[0]*len(self._variables), self._structure)
                [d,r] = self._structure._div_res(self,aux)
                return d
            if other.domain().is_sub_domain(self.domain().subdomain()):
                D = self.domain()._create_Poly_Structure_(other.domain())
                aux = Monomial_type(self._value >> other.K_(),self._varables[:],self._index[:],D)
                [d,r] = D._structure._div_res(aux,other>>D)
                return d 
            else:
                try:
                    D = (self.K_().union_domain(other.K_()))[0]
                    PD = Poly(D)
                    [d,r] = PD._div_res_(self >> PD,other >> PD)
                    return r
                except:
                    try:
                        return other.__rtruediv__(self)
                    except:
                        raise TypeError("I do not know how to opperate this elements")
        if isinstance(other, int):
            [d,r] = self._structure._div_res_(self,other)
            return d
        if isinstance(other, float):
            return self.__floordiv__((1>>self.domain())*other)
        try:
            return other.__rtruediv__(self)
        except:
            raise TypeError("I do not know how to opperate this elements")
    
    def __mod__(self,other):
        return self.residual_poly(other)
    
    def residual_poly(self, other):
        if isinstance(other, int):
            [d,r] = self._structure._div_res_(self,other)
            return r
        if isinstance(other, float):
            return self.residual_poly((1>>self.domain())*other)
        if isinstance(other, Symbol_type):
            return self.residual_poly(Monomial_type(1>>self._K(),[other],[1],self._structure))
        if isinstance(other, Monomial_type) or isinstance(other, Poly_type):
            if self._structure == other._structure:
                [d,r] = self._structure._div_res_(self,other)
                return r
            elif self.K_().is_sub_domain(other.K_()):
                [d,r] = self._structure._div_res_(self,other >> self._structure)
                return r
            elif other.K_().is_sub_domain(self.K_()):
                [d,r] = other._structure._div_res_(self >> other._structure,other)
                return r
            else:
                D = (self.K_().union_domain(other.K_()))[0]
                PD = Poly(D)
                [d,r] = PD._div_res_(self >> PD,other >> PD)
                return r
        if isinstance(other, Base_type.Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                [d,r] = self._structure._div_res_(self, Monomial_type(other>>self.K_(),self._variables[:],[0]*len(self._variables), self._structure))
                return r
            if other._structure.is_sub_domain(self._structure):
                [d,r] = other._structure._div_res_(self,other)
                return r
            if other.domain().is_sub_domain(self._structure.subdomain()):
                D = self.domain()._create_Poly_Structure_(other.domain())
                [d,r] = D._div_res_(self >> D,other>>D)
                return r
            else:
                D = (self.K_().union_domain(other.K_()))[0]
                PD = Poly(D)
                [d,r] = PD._div_res_(self >> PD,other >> PD)
                return r
        try:
            return other.__rtruediv__(self)
        except:
            raise TypeError("I do not know how to opperate this elements")          
        

    
### Polynomial class
def _own_order_monomials_(a,b):
    if a._totIndex == 0:
        if b._totIndex > 0:
            return -1
        else:
            return 0
    if b._totIndex == 0:
        return 1
    for i in range(len(a._varOrder - 1), -1 , -1):
        if a._index[a._varOrder[i]] < b._index[b._varOrder[i]]:
            return -1
        elif a._index[a._varOrder[i]] > b._index[b._varOrder[i]]:
            return 1
    return 0
            
def _get_own_order_monomials_(M):
    i = 0
    V = list(range(len(M)))
    while i < (len(M) - 1):
        if _own_order_monomials_(M[i], M[i+1]) > 0:
            v = V[i]
            V[i] = V[i+1]
            V[i+1] = v
            m = M[i]
            M[i] = M[i+1]
            M[i+1] = m
        i += 1
            
    return V
    

class Poly_type(Base_type.Base_type):
    def __init__(self,monomials,variables,domain,degree, total_degree,order = None, varOrder = None):
        self._structure = domain
        self._value = monomials
        self._variables = variables
        self._degree = degree
        self._mdegree = total_degree
        self._ord = order
        self._varOrder = varOrder
        
    def cp(self):
        if self._varOrder == None:
            return Poly_type(self._value[:], self._variables[:], self.domain(), self._degree, self._mdegree)
        else:
            return Poly_type(self._value[:], self._variables[:], self.domain(), self._degree, self._mdegree, self._ord[:], self._varOrder[:])
    
    def is_zero(self):
        if len(self._value) == 1:
            return self._value[0].is_zero()
        return False

    def _sign_value_(self):
        if (self.domain().subdomain().characteristic() == 0):
            lc = self.leading_coeff()
            if isinstance(lc,int):
                if lc >= 0:
                    return 0
                else:
                    return 1
            elif len(lc) > 1:
                return 0
            else:
                return lc._sign_value_()
        else:
            return 0
        
    def __hash__(self):
        aux = 0
        f = self.transform_to_min_domain()
        if not isinstance(f,Poly_type):
            return hash(f)
        for i in range(len(self._value)):
            aux += hash(self._value[i])
        return aux
        
    def _len_value_(self):
        if len(self._value)>1:
            return len(self._value)
        else:
            return self._value[0]._len_value_()
        
    def _print_(self):
        if PGl.POL_DES_SORT:
            if self._varOrder == None:
                value = self._value[:]
            else:
                value = [self._value[i] for i in self._ord]
        else:
            if self._varOrder == None:
                value = self._value[::-1]
            else:
                value = [self._value[self._varOrder[i]] for i in range(len(self._ord)-1,-1,-1)]
        
        if PGl.POL_EXP:
            res = ''
            for i in range(len(value)):
                aux = value[i]._print_()
                if i == 0:
                    res += aux
                elif aux[0] == '-':
                    res = res + ' - ' + aux[1:]
                else:
                    res = res + ' + ' + aux
            return res
        else:
            res = ''
            flagIni = True
            i = 0
            
            while(i < len(value)):
                nTerms = 1
                d = value[i].degree()
                jind = len(value)
                resaux = value[i]._print_coeff_()
                if len(resaux) == 0:
                    fcanviN = False
                elif resaux[0] == '-':
                    fcanviN = True
                else:
                    fcanviN = False
                poscanviN = []
                for j in range(i+1,len(value)):
                    dj = value[j].degree()
                    if dj == d:
                        nTerms += 1
                        rTerm = value[j]._print_coeff_()
                        if len(rTerm) == 0:
                            fcanviN = False
                            resaux += ' + 1'
                        elif rTerm == '-':
                            poscanviN.append(len(resaux) + 1)
                            resaux += ' - 1'
                        if rTerm[0] == '-':
                            poscanviN.append(len(resaux) + 1)
                            resaux += ' - ' + rTerm[1:]
                        else:
                            fcanviN = False
                            resaux += ' + ' + rTerm
                    else:
                        jind = j - 1
                        break
                if nTerms > 1:
                    if resaux[0] == ' ':
                        resaux = '1' + resaux
                    elif resaux[:2] == '- ':
                        resaux = '-1' + resaux[1:]
                        if fcanviN:
                            for j in range(len(poscanviN)):
                                poscanviN[j] +=1
                if self._varOrder == None:
                    if d > 1:
                        if nTerms > 1:
                            if fcanviN:
                                resaux = list(resaux)
                                for j in poscanviN:
                                    resaux[j] = '+'
                                resaux = ''.join(resaux)
                                resaux = resaux[1:]
                                if flagIni:
                                    flagIni = False
                                    res += '-(' + resaux + ')*' + self._variables[-1]._print_() + '**' + d.__str__()
                                else:
                                    res += ' - (' + resaux + ')*' + self._variables[-1]._print_() + '**' + d.__str__()
                            else:
                                if flagIni:
                                    flagIni = False
                                    res += '(' + resaux + ')*' + self._variables[-1]._print_() + '**' + d.__str__()
                                else:
                                    res += ' + (' + resaux + ')*' + self._variables[-1]._print_() + '**' + d.__str__()
                        else:
                            flagEx = False
                            if not isinstance(value[i]._value,int):
                                if value[i]._value != 1:
                                    if not ((value[i]._value == 1)and(self._structure._characteristic == 0)):
                                        if value[i]._value._len_value_() > 1:
                                            flagEx = True
                                            if flagIni:
                                                flagIni = False
                                                res += '(' + resaux + ')*' + self._variables[-1]._print_() + '**' + d.__str__()
                                            else:
                                                res += ' + (' + resaux + ')*' + self._variables[-1]._print_() + '**' + d.__str__()
                            if not flagEx:
                                if flagIni:
                                    flagIni = False
                                    if len(resaux) == 0:
                                        res += self._variables[-1]._print_() + '**' + d.__str__()
                                    elif resaux == '-':
                                        res += resaux + self._variables[-1]._print_() + '**' + d.__str__()
                                    else:
                                        res += resaux + '*' + self._variables[-1]._print_() + '**' + d.__str__()
                                else:
                                    if len(resaux) == 0:
                                        res += ' + ' + self._variables[-1]._print_() + '**' + d.__str__()
                                    elif resaux == '-':
                                        res += ' - ' + self._variables[-1]._print_() + '**' + d.__str__()
                                    elif resaux[0] == '-':
                                        res += ' - ' + resaux[1:] + '*' + self._variables[-1]._print_() + '**' + d.__str__()
                                    else:
                                        res += ' + ' + resaux + '*' + self._variables[-1]._print_() + '**' + d.__str__()
                    elif d > 0:
                        if nTerms > 1:
                            if fcanviN:
                                resaux = list(resaux)
                                for j in poscanviN:
                                    resaux[j] = '+'
                                resaux = ''.join(resaux)
                                resaux = resaux[1:]
                                if flagIni:
                                    flagIni = False
                                    res += '-(' + resaux + ')*' + self._variables[-1]._print_() 
                                else:
                                    res += ' - (' + resaux + ')*' + self._variables[-1]._print_() 
                            else:
                                if flagIni:
                                    flagIni = False
                                    res += '(' + resaux + ')*' + self._variables[-1]._print_() 
                                else:
                                    res += ' + (' + resaux + ')*' + self._variables[-1]._print_() 
                        else:
                            flagEx = False
                            if not isinstance(value[i]._value,int):
                                if value[i]._value != 1:
                                    if not ((value[i]._value == 1)and(self._structure._characteristic == 0)):
                                        if value[i]._value._len_value_() > 1:
                                            flagEx = True
                                            if flagIni:
                                                flagIni = False
                                                res += '(' + resaux + ')*' + self._variables[-1]._print_() 
                                            else:
                                                res += ' + (' + resaux + ')*' + self._variables[-1]._print_() 
                            if not flagEx:
                                if flagIni:
                                    flagIni = False
                                    if len(resaux) == 0:
                                        res += self._variables[-1]._print_() 
                                    elif resaux == '-':
                                        res += resaux + self._variables[-1]._print_()
                                    else:
                                        res += resaux + '*' + self._variables[-1]._print_() 
                                else:
                                    if len(resaux) == 0:
                                        res += ' + ' + self._variables[-1]._print_() 
                                    elif resaux == '-':
                                        res += ' - ' + self._variables[-1]._print_() 
                                    elif resaux[0] == '-':
                                        res += ' - ' + resaux[1:] + '*' + self._variables[-1]._print_() 
                                    else:
                                        res += ' + ' + resaux + '*' + self._variables[-1]._print_() 
                    else:
                        if flagIni:
                            flagIni = False
                            res += resaux  
                        else:
                            if resaux[0] == '-':
                                res += ' - ' + resaux[1:] 
                            else:
                                res += ' + ' + resaux 
                else:
                    if d > 1:
                        if nTerms > 1:
                            if fcanviN:
                                resaux = list(resaux)
                                for j in poscanviN:
                                    resaux[j] = '+'
                                resaux = ''.join(resaux)
                                resaux = resaux[1:]
                                if flagIni:
                                    flagIni = False
                                    res += '-(' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                else:
                                    res += ' - (' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                            else:
                                if flagIni:
                                    flagIni = False
                                    res += '(' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                else:
                                    res += ' + (' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                        else:
                            flagEx = False
                            if not isinstance(value[i]._value,int):
                                if value[i]._value != 1:
                                    if not ((value[i]._value == 1)and(self._structure._characteristic == 0)):
                                        if value[i]._value._len_value_() > 1:
                                            flagEx = True
                                            if flagIni:
                                                flagIni = False
                                                res += '(' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                            else:
                                                res += ' + (' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                            if not flagEx:
                                if flagIni:
                                    flagIni = False
                                    if len(resaux) == 0:
                                        res += self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                    elif resaux == '-':
                                        res += resaux + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                    else:
                                        res += resaux + '*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                else:
                                    if len(resaux) == 0:
                                        res += ' + ' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                    elif resaux == '-':
                                        res += ' - ' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                    elif resaux[0] == '-':
                                        res += ' - ' + resaux[1:] + '*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                                    else:
                                        res += ' + ' + resaux + '*' + self._variables[self._varOrder[-1]]._print_() + '**' + d.__str__()
                    elif d > 0:
                        if nTerms > 1:
                            if fcanviN:
                                resaux = list(resaux)
                                for j in poscanviN:
                                    resaux[j] = '+'
                                resaux = ''.join(resaux)
                                resaux = resaux[1:]
                                if flagIni:
                                    flagIni = False
                                    res += '-(' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() 
                                else:
                                    res += ' - (' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() 
                            else:
                                if flagIni:
                                    flagIni = False
                                    res += '(' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() 
                                else:
                                    res += ' + (' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() 
                        else:
                            flagEx = False
                            if not isinstance(value[i]._value,int):
                                if value[i]._value != 1:
                                    if not ((value[i]._value == 1)and(self._structure._characteristic == 0)):
                                        if value[i]._value._len_value_() > 1:
                                            flagEx = True
                                            if flagIni:
                                                flagIni = False
                                                res += '(' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_() 
                                            else:
                                                res += ' + (' + resaux + ')*' + self._variables[self._varOrder[-1]]._print_()
                            if not flagEx:
                                if flagIni:
                                    flagIni = False
                                    if len(resaux) == 0:
                                        res += self._variables[self._varOrder[-1]]._print_() 
                                    elif resaux == '-':
                                        res += resaux + self._variables[self._varOrder[-1]]._print_() 
                                    else:
                                        res += resaux + '*' + self._variables[self._varOrder[-1]]._print_() 
                                else:
                                    if len(resaux) == 0:
                                        res += ' + ' + self._variables[self._varOrder[-1]]._print_() 
                                    elif resaux == '-':
                                        res += ' - ' + self._variables[self._varOrder[-1]]._print_() 
                                    elif resaux[0] == '-':
                                        res += ' - ' + resaux[1:] + '*' + self._variables[self._varOrder[-1]]._print_() 
                                    else:
                                        res += ' + ' + resaux + '*' + self._variables[self._varOrder[-1]]._print_() 
                    else:
                        if flagIni:
                            flagIni = False
                            res += resaux  
                        else:
                            if resaux[0] == '-':
                                res += ' - ' + resaux[1:] 
                            else:
                                res += ' + ' + resaux
                i = jind + 1
            if len(res) == 0:
                return '0'
            else:
                return res
            
    def __str__(self):
        if PGl.TYPES == 1:
            if len(self._variables) == 0:
                return self._print_() + ' :: ' + self._structure._print_()                    
            if self.domain()._printdefault:
                res = self._print_() + ' :: ' + self._structure.subdomain()._print_()
                if self._varOrder == None:
                    for i in range(len(self._variables)):
                        res += '[' + self._variables[i]._print_() + ']'
                    return res
                else:
                    for i in range(len(self._variables)):
                        res += '[' + self._variables[self._varOrder[i]]._print_() + ']'
                    return res
            else:
                return self._print_() + ' :: ' + self._structure._print_()
        else:
            return self._print_()
        
    def __repr__(self):
        return self.__str__()

    def reverse(self,var=None):
        if len(self._variables) == 0:
            return self.cp()
        if var == None:
            d = self.degree()
            return self.reverse_m(d)
        else:
            x = _get_symbol_(var)
            d = self.degree(x)
            return self.reverse_m(d,x)
                
    def reverse_m(self,m,var=None):
        if var == None:
            if len(self._variables) == 0:
                raise TypeError('This polynomial has no variable')
            if m < self.degree():
                 raise TypeError('m must be equal or larger than the degree of the polynomial')
            mtot = self._mdegree - m
            i = 0
            aux = []
            x = self.variable()
            while(i < len(self._value)):
                aux2 = [self._value[i].reverse_m(m, x)]
                d = self._value[i].degree()
                if mtot < aux2[-1]._totIndex:
                    mtot = aux2[-1]._totIndex
                jind = len(self._value)
                for j in range(i+1,len(self._value)):
                    dj = self._value[j].degree()
                    if dj == d:
                        aux2.append(self._value[j].reverse_m(m, x))
                        if mtot < aux2[-1]._totIndex:
                            mtot = aux2[-1]._totIndex
                    else:
                        jind = j
                        break
                aux = aux2 + aux
                i = jind
            if self._ord == None:
                return Poly_type(aux,self._variables[:],self._structure, aux[0].degree(),mtot)
            else:
                aord = _get_own_order_monomials_(aux[:])
                return Poly_type(aux,self._variables[:],self._structure, aux[0]._degree, mtot, aord, self._varOrder[:])
        else:
            x = _get_symbol_(var)
            if len(self._variables) == 0:
                return self._value[0].reverse_m(m,x)
            dx = self.degree(x)
            if m < dx:
                 raise TypeError('m must be equal or larger than the degree of the polynomial')
            if len(self._value) == 1:
                return self._value[0].reverse_m(m,x)
            else:
                mtot = self._mdegree - m
                i = 0
                aux = []
                while(i < len(self._value)):
                    aux2 = [self._value[i].reverse_m(m, x)]
                    d = self._value[i].degree(x)
                    if mtot < aux2[-1]._totIndex:
                        mtot = aux2[-1]._totIndex
                    jind = len(self._value)
                    for j in range(i+1,len(self._value)):
                        dj = self._value[j].degree(x)
                        if dj == d:
                            aux2.append(self._value[j].reverse_m(m, x))
                            if mtot < aux2[-1]._totIndex:
                                mtot = aux2[-1]._totIndex
                        else:
                            jind = j
                            break
                    aux = aux2 + aux
                    i = jind
                if self._ord == None:
                    return Poly_type(aux,self._variables[:],self._structure, aux[0].degree(),mtot)
                else:
                    aord = _get_own_order_monomials_(aux[:])
                    return Poly_type(aux,self._variables[:],self._structure, aux[0]._degree, mtot, aord, self._varOrder[:])
       
    def trunc(self,k,var=None):
        if k < 1:
            return Monomial_type(0>>self._structure.subdomain(),self._variables[:], [0]*len(self._variables), self._structure)
        if (var == None):
            
            
            if self._ord == None:
                if self._value[-1].degree() >= k:
                    return Monomial_type(0>>self._structure.subdomain(),self._variables[:], [0]*len(self._variables), self._structure)
                aux = []
                pos = -1
                for i in range(len(self._value)):
                    if self._value[i].degree() < k:
                        pos = i
                        break
                    
                d = self._value[pos]._degree
                md = self._value[pos].mdegree()
                aux.append(self._value[pos].cp())
                for i in range(pos + 1, len(self._value)):
                    aux.append(self._value[i].cp())
                    if md < self._value[i].mdegree():
                        md = self._value[i].mdegree()
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,self._variables[:],self._structure, d,md)
            else:
                if self._value[self._ord[-1]].degree() >= k:
                    return Monomial_type(0>>self._structure.subdomain(),self._variables[:], [0]*len(self._variables), self._structure)
                aux = []
                v = [0]*len(self._ord)
                d = 0
                md = 0
                for i in range(len(self._value)):
                    if self._value[i].degree() < k:
                        aux.append(self._value[i].cp())
                        di = self._value[i]._degree
                        if di > d:
                            d = di
                        mdi = self._value[i].mdegree()
                        if md < mdi:
                            md = mdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,self._variables[:],self._structure,d,md,aord, self._varOrder[:])  
        else:
            x = _get_symbol_(var)
            if self.trailing_degree(x) >= k:
                 return Monomial_type(0>>self._structure.subdomain(),self._variables[:], [0]*len(self._variables), self._structure)
            if self._ord == None:
                aux = []
                d = 0
                md = 0
                for i in range(len(self._value)):
                    if self._value[i].degree(x) < k:
                        aux.append(self._value[i].cp())
                        di = self._value[i]._degree
                        if di > d:
                            d = di
                        mdi = self._value[i].mdegree()
                        if md < mdi:
                            md = mdi
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,self._variables[:],self._structure,d,md)
            else:
                aux = []
                v = [0]*len(self._ord)
                d = 0
                md = 0
                for i in range(len(self._value)):
                    if self._value[i].degree(x) < k:
                        aux.append(self._value[i].cp())
                        di = self._value[i]._degree
                        if di > d:
                            d = di
                        mdi = self._value[i].mdegree()
                        if md < mdi:
                            md = mdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,self._variables[:],self._structure,d,md,aord, self._varOrder[:])
            
    def trunc_total_degree(self, d):
        aux = []
        vord = [0]*len(self._value)
        for i in self._value:
            if i._totIndex < d:
                aux.append(i.cp())
            elif self._ord != None:
                vord[i] = 1
        if self._ord == None:
            if len(aux) == 0:
                return Monomial_type(0>>self._structure.subdomain(), self._variables[:], [0]*len(self._variables), self._structure)
            elif len(aux) == 1:
                return aux[0]
            else:
                if self._ord == None:
                    return Poly_type(aux,self._variables[:], self._structure, aux[0]._index[-1], d)
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if vord[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(vord[:self._ord[i]]))
                    return Poly_type(aux,self._variables[:], self._structure, aux[0]._index[-1], d, aord, self._varOrder[:])
    
    def get_subpolynomial_degree(self,d, var = None):
        return self.trunc(d + 1 ,var) - self.trunc(d,var)
    
    def get_subpolynomial_total_degree(self,d):
        return self.trunc_total_degree(d + 1) - self.trunc_total_degree(d)
        
    def truncated_inverse(self,k):
        return self._structure._trunc_inv_2_(self,k)       
        
    def __neg__(self):
        return self._structure._opposite_(self)
    
    def opposite(self):
        return self._structure._opposite_(self)
    
    def negation(self):
        return self._structure._opposite_(self)        
        
    def degree(self, symb = None):
        if len(self._variables) == 0:
            return 0
        if symb == None:
            if self._varOrder == None:
                return self._degree
            else:
                return self._value[self._ord[0]].degree()
        else:
            x = _get_symbol_(symb)
            if x == self.variable():
                return self.degree()
            d = self._value[0].degree(x)
            for i in range(1,len(self._value)):
                di = self._value[i].degree(x)
                if d < di:
                    d = di
            return d
    
    def mdegree(self):
        return self._mdegree
    
    def total_degree(self):
        return self._mdegree
        
    def leading_coeff(self, var = None):
        if var == None:
            if self._ord == None:
                d = self._value[0].degree()
                aux = [self._value[0].leading_coeff()]
                dn = aux[0].degree()
                td = aux[0].mdegree()
                for i in range(1,len(self._value)):
                    di = self._value[i].degree()
                    if di == d:
                        aux.append(self._value[i].leading_coeff())
                        tdi = aux[-1].mdegree()
                        if tdi > td:
                            td = tdi
                    else:
                        break
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,self._variables[:-1],self._structure,dn,td)
            else:
                aux = []
                d = self._value[self._ord[0]].degree()
                v = [0]*len(self._value)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree()
                    if di == d:
                        aux.append(self._value[i].leading_coeff())
                        tdi = aux[-1].mdegree()
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,aux.variables[:],self._structure,dn,td,aord, aux[0]._varOrder[:])
        else:
            x = _get_symbol_(var)
            if x == self.variable():
                return self.leading_coeff()
            if self._ord == None:
                d = self.degree(x)
                dn = 0
                td = 0
                aux = []
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].leading_coeff(x))
                        tdi = aux[-1].mdegree()
                        if tdi > td:
                            td = tdi
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
            else:
                aux = []
                d = self.degree(x)
                v = [0]*len(self._value)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].leading_coeff(x))
                        tdi = aux[-1].mdegree()
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    if aux[0]._varOrder == None:
                        return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td,aord, aux[0]._varOrder[:])
    
    def leading_term(self, var = None):
        if var == None:
            if self._ord == None:
                d = self._value[0].degree()
                aux = [self._value[0].leading_term()]
                td = aux[0].mdegree()
                for i in range(1,len(self._value)):
                    di = self._value[i].degree()
                    if di == d:
                        aux.append(self._value[i].leading_term())
                        tdi = aux[-1].mdegree()
                        if tdi > td:
                            td = tdi
                    else:
                        break
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,self._variables[:],self._structure,d,td)
            else:
                aux = []
                d = self._value[self._ord[0]].degree()
                v = [0]*len(self._value)
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree()
                    if di == d:
                        aux.append(self._value[i].leading_term())
                        tdi = aux[-1].mdegree()
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,self._variables[:],self._structure,d,td,aord, self._varOrder[:])
        else:
            x = _get_symbol_(var)
            if x == self.variable():
                return self.leading_term()
            if self._ord == None:
                aux = []
                d = self.degree(x)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].leading_term(x))
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                        tdi = aux[-1].mdegree()
                        if tdi > td:
                            td = tdi
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
            else:
                aux = []
                d = self.degree(x)
                v = [0]*len(self._value)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].leading_term(x))
                        tdi = aux[-1].mdegree()
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    if aux[0]._varOrder == None:
                        return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td,aord, self._varOrder[:])
    
    def trailing_coeff(self, var = None):
        if var == None:
            if self._ord == None:
                d = self._value[-1].degree()
                aux = [self._value[-1].trailing_coeff()]
                td = aux[0].mdegree()
                for i in range(len(self._value)-2,-1,-1):
                    di = self._value[i].degree()
                    if di == d:
                        aux = [self._value[i].trailing_coeff()] + aux
                        tdi = aux[0].mdegree()
                        if tdi > td:
                            td = tdi
                    else:
                        break
                dn = aux[0].degree()
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,self._variables[:-1],self._structure,dn,td)
            else:
                aux = []
                d = self._value[self._ord[-1]].degree()
                v = [0]*len(self._value)
                dn = 0
                td = 0
                for i in range(len(self._value)-1,-1,-1):
                    di = self._value[i].degree()
                    if di == d:
                        aux = [self._value[i].trailing_coeff()] + aux
                        tdi = aux[0].mdegree()
                        dni = aux[0].degree()
                        if dni > dn:
                            dn = dni
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td,aord, aux[0]._varOrder[:])
        else:
            x = _get_symbol_(var)
            if x == self.variable():
                return self.trailing_coeff()
            if self._ord == None:
                aux = []
                d = self.trailing_degree(x)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].trailing_coeff(x))
                        tdi = aux[-1].mdegree()
                        if tdi > td:
                            td = tdi
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
            else:
                aux = []
                d = self.trailing_degree(x)
                v = [0]*len(self._value)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].trailing_coeff(x))
                        tdi = aux[-1].mdegree()
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    if aux[0]._varOrder == None:
                        return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td,aord, aux[0]._varOrder[:])

    def trailing_term(self, var = None):
        if var == None:
            if self._ord == None:
                d = self._value[-1].degree()
                aux = [self._value[-1].trailing_term()]
                td = aux[0].mdegree()
                for i in range(len(self._value)-2,-1,-1):
                    di = self._value[i].degree()
                    if di == d:
                        aux = [self._value[i].trailing_term()] + aux
                        tdi = aux[0].mdegree()
                        if tdi > td:
                            td = tdi
                    else:
                        break
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,self._variables[:],self._structure,d,td)
            else:
                aux = []
                d = self._value[self._ord[-1]].degree()
                v = [0]*len(self._value)
                dn = 0
                td = 0
                for i in range(len(self._value)-1,-1,-1):
                    di = self._value[i].degree()
                    if di == d:
                        aux = [self._value[i].trailing_term()] + aux
                        tdi = aux[0].mdegree()
                        dni = aux[0].degree()
                        if dni > dn:
                            dn = dni
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,self._variables[:],self._structure,dn,td,aord, self._varOrder[:])
        else:
            x = _get_symbol_(var)
            if x == self.variable():
                return self.trailing_term()
            if self._ord == None:
                aux = []
                d = self.trailing_degree(x)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].trailing_term(x))
                        tdi = aux[-1].mdegree()
                        if tdi > td:
                            td = tdi
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                if len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,self._variables[:],self._structure,dn,td)
            else:
                aux = []
                d = self.trailing_degree(x)
                v = [0]*len(self._value)
                dn = 0
                td = 0
                for i in range(len(self._value)):
                    di = self._value[i].degree(x)
                    if di == d:
                        aux.append(self._value[i].trailing_term(x))
                        tdi = aux[-1].mdegree()
                        dni = aux[-1].degree()
                        if dni > dn:
                            dn = dni
                        if tdi > td:
                            td = tdi
                    else:
                        v[i] = 1
                if len(aux) == 1:
                    return aux[0]
                else:
                    if aux[0]._varOrder == None:
                        return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
                    aord = []
                    for i in range(len(self._ord)):
                        if v[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                    return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td,aord, aux[0]._varOrder[:])
    
    def constant_coeff(self, var = None):
        if var == None:
            if self._ord == None:
                d = self._value[-1].degree()
                if d == 0:
                    aux = [self._value[-1].constant_coeff()]
                    td = aux[0].mdegree()
                    for i in range(len(self._value)-2,-1,-1):
                        di = self._value[i].degree()
                        if di == 0:
                            aux = [self._value[i].constant_coeff()] + aux
                            tdi = aux[0].mdegree()
                            if tdi > td:
                                td = tdi
                        else:
                            break
                    dn = aux[0].degree()
                    if len(aux) == 1:
                        return aux[0]
                    else:
                        return Poly_type(aux,self._variables[:-1],self._structure,dn,td)
                else:
                    return Monomial_type(0>>self._structure.subdomain(),self._variables[:-1], [0]*(len(self._variables)-1), self._structure)
            else:
                d = self._value[self._ord[-1]].degree()
                if d == 0:
                    aux = []
                    v = [0]*len(self._value)
                    dn = 0
                    td = 0
                    for i in range(len(self._value)-1,-1,-1):
                        di = self._value[i].degree()
                        if di == 0:
                            aux = [self._value[i].trailing_coeff()] + aux
                            tdi = aux[0].mdegree()
                            dni = aux[0]._degree
                            if dni > dn:
                                dn = dni
                            if tdi > td:
                                td = tdi
                        else:
                            v[i] = 1
                    if len(aux) == 1:
                        return aux[0]
                    else:
                        aord = []
                        for i in range(len(self._ord)):
                            if v[self._ord[i]] == 0:
                                aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                        return Poly_type(aux,self._variables[:-1],self._structure,dn,td,aord, aux[0]._varOrder[:])
                else:
                    return Monomial_type(0>>self._structure.subdomain(),self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1]+1:], [0]*(len(self._variables)-1), self._structure)
        else:
            x = _get_symbol_(var)
            if x == self.variable():
                return self.constant_coeff()
            if self._ord == None:
                aux = []
                d = self.trailing_degree(x)
                if d == 0:
                    dn = 0
                    td = 0
                    for i in range(len(self._value)):
                        di = self._value[i].degree(x)
                        if di == 0:
                            aux.append(self._value[i].trailing_coeff(x))
                            tdi = aux[-1].mdegree()
                            if tdi > td:
                                td = tdi
                            dni = aux[-1].degree()
                            if dni > dn:
                                dn = dni
                    if len(aux) == 1:
                        return aux[0]
                    else:
                        return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
                else:
                    if x in self._variables:
                        ind = self._variables.index(x)
                        return Monomial_type(0>>self._structure.subdomain(),self._variables[:ind] + self._variables[ind+1:],[0]*(len(self._variables)-1), self._structure)
                    else:
                        return Monomial_type(0>>self._structure.subdomain(),self._variables[:],[0]*len(self._variables), self._structure)
            else:
                d = self.trailing_degree(x)
                if d == 0:
                    aux = []
                    v = [0]*len(self._value)
                    dn = 0
                    td = 0
                    for i in range(len(self._value)):
                        di = self._value[i].degree(x)
                        if di == 0:
                            aux.append(self._value[i].trailing_coeff(x))
                            tdi = aux[-1].mdegree()
                            dni = aux[-1]._degree
                            if dni > dn:
                                dn = dni
                            if tdi > td:
                                td = tdi
                        else:
                            v[i] = 1
                    if len(aux) == 1:
                        return aux[0]
                    else:
                        if aux[0]._varOrder == None:
                            return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td)
                        aord = []
                        for i in range(len(self._ord)):
                            if v[self._ord[i]] == 0:
                                aord.append(self._ord[i] - sum(v[:self._ord[i]]))
                        return Poly_type(aux,aux[0]._variables[:],self._structure,dn,td,aord, aux[0]._varOrder[:])
                else:
                    return Monomial_type(0>>self._structure.subdomain(),self._variables[:],[0]*len(self._variables), self._structure)


    def constant_term(self, var = None):
        return self.trunc(1,var)
    
    def trailing_degree(self, var = None):
        if var == None:
            return self._value[-1].trailing_degree()
        else:
            x = _get_symbol_(var)
            if x == self.variable():
                return self._value[-1].trailing_degree()
            d = self._degree
            for i in self._value:
                di = i.trailing_degree(x)
                if di < d:
                    d = di
            return d
        
    def zero_degree_term(self):
        if self._value[-1]._totIndex == 0:
            return self._value[-1].cp()
        else:
            return Monomial_type(0>>self._structure.subdomain(), self._variables[:], [0]*len(self._variables), self._structure)
    
    def der(self, var = None):
        if var == None:
            aux = []
            vord = [0]*len(self._value)
            for i in range(len(self._value)):
                deri = self._value[i].der()
                if not deri.is_zero():
                    aux.append(deri)
                elif self._ord != None:
                    vord[i] = 1
            if len(aux) == 0:
                return Monomial_type(0>>self._structure.subdomain(),self._variables[:],[0]*len(self._variables), self._structure)
            elif len(aux) == 1:
                return aux[0]
            else:
                if self._ord == None:
                    return Poly_type(aux,aux[0]._variables[:],self._structure,self._degree - 1, self._mdegree - 1)
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if vord[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(vord[:self._ord[i]]))
                    return Poly_type(aux,aux[0]._variables[:],self._structure,self._degree - 1, self._mdegree - 1,aord, self._varOrder[:])
        else:
            aux = []
            x = _get_symbol_(var)
            if x == self.variable():
                return self.der()
            vord = [0]*len(self._value)
            for i in range(len(self._value)):
                deri = self._value[i].der(x)
                if not deri.is_zero():
                    aux.append(deri)
                elif self._ord != None:
                    vord[i] = 1
            if len(aux) == 0:
                return Monomial_type(0>>self._structure.subdomain(),self._variables[:],[0]*len(self._variables), self._structure)
            elif len(aux) == 1:
                return aux[0]
            else:
                if self._ord == None:
                    return Poly_type(aux,aux[0]._variables[:],self._structure,self._degree - 1, self._mdegree - 1)
                else:
                    aord = []
                    for i in range(len(self._ord)):
                        if vord[self._ord[i]] == 0:
                            aord.append(self._ord[i] - sum(vord[:self._ord[i]]))
                    return Poly_type(aux,aux[0]._variables[:],self._structure,self._degree - 1, self._mdegree - 1,aord, self._varOrder[:])

    def to_dict(self):
        T = {}
        for i in self._value:
            T[tuple(i._index)] = i._value
        return T
    
    def skeleton(self):
        aval = []
        aindex = []
        for i in self._value:
            aval.append(i._value)
            aindex.append(i._index)
        return aval, aindex
        
    def eval_poly(self, value, variables = None):
        res = self._value[0].eval_monomial(value,variables)
        for i in range(1,len(self._value)):
            res +=  self._value[i].eval_monomial(value,variables)
        return res
    
    def evaluate(self,value,variables = None):
        return self.eval_poly(value,variables)
    
    def coeffs(self, var = None):
        d = self.degree(var)
        if var == None:
            if self._varOrder == None:
                cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:-1], [0]*(len(self._variables)-1), self._structure)]*(d + 1)
            else:
                cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:self._varOrder[-1]] + self._variables[self._varOrder[-1] + 1], [0]*(len(self._variables)-1), self._structure, self._varOrder[:-1])]*(d + 1)
            for i in self._value:
                di = i.degree()
                cf[d - di] = cf[d - di] + i.leading_coeff()
            return cf
        else:
            x = _get_symbol_(var)
            if x in self._variables:
                ind = self._variables.index(x)
                if self._varOrder == None:
                    cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:ind] + self._variables[ind+1:], [0]*(len(self._variables)-1), self._structure)]*(d + 1)
                else:
                    varOrder = []
                    for i in range(len(self._varOrder)):
                        if self._varOrder[i] > ind:
                            varOrder.append(self._varOrder[i] - 1)
                        elif self._varOrder[i] < ind:
                            varOrder.append(self._varOrder[i])
                    cf = [Monomial_type(0>>self._structure.subdomain(), self._variables[:ind] + self._variables[ind+1:], [0]*(len(self._variables)-1), self._structure, varOrder)]*(d + 1)
                for i in self._value:
                    di = i.degree(x)
                    cf[d - di] = cf[d - di] + i.leading_coeff(x)
                return cf
            else:
                if self._varOrder == None:
                    return [self.cp()]
                else:
                    return [Monomial_type(self._value,self._variables[:], self._index[:],self._structure)]
    
    def to_list(self):
        return self.coeffs()
    
    def min_domain(self):
        return self._structure
        
    def transform_to_min_domain(self):
        fvar = [1]*len(self._variables)
        for i in self._value:
            for j in range(len(self._variables)):
                if i._index[j] > 0:
                    fvar[j] = 0
        avariables = []
        aux = []
        varOrder = []
        for i in range(len(self._variables)):
            if fvar[i] == 0:
                avariables.append(self._variables[i])
                if self._ord != None:
                    varOrder.append(self._varOrder[i] - fvar[:self._varOrder[i]])
        for i in self._value:
            aindex = []
            for j in range(len(self._variables)):
                if fvar[j] == 0:
                    aindex.append(i._index[j])
            if self._ord == None:
                aux.append(Monomial_type(i._value,avariables[:],aindex,self._structure))
            else:
                aux.append(Monomial_type(i._value,avariables[:],aindex,self._structure, varOrder))
        if self._ord == None:
            return Poly_type(aux,avariables[:],self._structure,self._degree, self._mdegree)
        else:
            return Poly_type(aux,avariables[:],self._structure,self._degree, self._mdegree, self._ord[:], varOrder)
    
    def change_characteristic(self,n):
        if n == self._structure.characteristic():
            return self.cp()
        cself = self.domain().characteristic()
        if (cself%n) == 0:
            D = self._structure.change_characteristic(n)
            aux = []
            md = 0
            for i in self._value:
                aux.append(i.change_characteristic(n))
                if aux[-1]._value == 0:
                    aux = aux[:-1]
                else:
                    if aux[-1]._totIndex > md:
                        md = aux[-1]._totIndex
            if len(aux) == 0:
                return Monomial_type(0>>D.subdomain(),self._variables[:],self._index[:],D)
            elif len(aux) == 1:
                return aux[0]
            else:
                if self._ord == None:
                    return Poly_type(aux,self._variables[:],D,self._degree, md)
                else:
                    aord = _get_own_order_monomials_(aux[:])
                    return Poly_type(aux,self._variables[:],D,self._degree, md, aord, self._varOrder[:])
        raise TypeError('Wrong characteristic') 
        
    def transform_to_ZZ(self):
        if self.domain().characteristic() == 0:
            return self.cp()
        else:
            D = self.domain().transform_to_ZZ()
            aux = []
            for i in self._value:
                aux.append(i.transform_to_ZZ())
            if self._ord == None:
                return Poly_type(aux,self._variables[:],D,self._degree, self._mdegree)
            else:
                return Poly_type(aux,self._variables[:],D,self._degree, self._mdegree, self._ord[:], self._varOrder[:])
    
    def is_invertible(self):
        return self._structure.is_invertible(self)
    
    def K_(self):
        return self.domain().K_()
                   
    def variable(self):
        if self._varOrder == None:
            return self._variables[-1]
        else:
            return self._variables[self._varOrder[-1]]
    
    def variables(self):
        return self._variables[:]
    
    def nvariables(self):
        return len(self._variables)
    
    def add_variable(self,v):
        [avariables, pos] = _add_symbol_(self._variables[:],v)
        if len(avariables) == len(self._variables):
            return self.cp()
        else:
            aux = []
            for i in self._value:
                aux.append(Monomial_type(i._value, avariables[:], i._index[:pos] + [0] + i._index[pos:], self._structure))
            return Poly_type(aux, avariables, self._structure, self._degree, self._mdegree)
    
    def add_variables(self,V):
        p = self.cp()
        for v in V:
            p = p.add_variable(v)
        return p
        
    def variables_name(self):
        return [i._value for i in self._variables]
    
    def collect(self,var):
        x = _get_symbol_(var)
        if x == self.variable():
            return self.cp()
        aux = []
        for i in self._value:
            aux.append(i.collect(var))
        if aux[0]._varOrder == None:
            return Poly_type(aux, aux[0]._variables[:], self._structure, self._degree, self._mdegree)
        else:   
            aord = _get_own_order_monomials_(aux[:])
            return Poly_type(aux, aux[0]._variables[:], self._structure, self._degree, self._mdegree, aord, aux[0]._varOrder[:])
        
    def __pow__(self,other):
        if isinstance(other, int):
            if other < 0:
                raise TypeError("A polynomial cannot be invertible")
            elif other == 0:
                if self!=0:
                    if self._varOrder == None:
                        return Monomial_type(1>>self.domain().subdomain(),self._variables[:],[0]*len(self._variables),self.domain())
                    else:
                        return Monomial_type(1>>self.domain().subdomain(),self._variables[:],[0]*len(self._variables),self.domain(),self._varOrder[:])
                else:
                    raise TypeError("0**0 indetermination")
            else:
                if self._varOrder == None:
                    aux1 = Monomial_type(1>>self.domain().subdomain(),self._variables[:],[0]*len(self._variables),self.domain())
                else:
                    aux1 = Monomial_type(1>>self.domain().subdomain(),self._variables[:],[0]*len(self._variables),self.domain(),self._varOrder[:])
                aux2 = self.cp()
                while other>1:
                    if other%2==0:
                        other//=2
                        aux2 = self.domain()._mul_core_(aux2,aux2)
                    else:
                        aux1 = self.domain()._mul_core_(aux2,aux1)
                        aux2 = self.domain()._mul_core_(aux2,aux2)
                        other = (other-1)//2
                return self.domain()._mul_core_(aux1,aux2)
        else:
            raise TypeError("This element cannot be operated")
      
        
    def __floordiv__(self, other):
        if isinstance(other, Symbol_type):
            return self.__floordiv__(Monomial_type(1>>self.K_(),[other],[1], self._structure))
        if isinstance(other,Base_type.Base_type):
            if other.domain()._poly_structure:
                if self.domain() == other.domain():
                    [d,r] = self._structure._div_res_(self,other)
                    return d
                if self.domain().subdomain().is_sub_domain(other.domain().subdomain()):
                    [d,r] = self._structure._div_res(self,other >> self._structure)
                    return d
                if other.domain().subdomain().is_sub_domain(self.domain().subdomain()):
                    [d,r] = other._structure._div_res(self >> other._structure,other)
                    return d 
                else:
                    try:
                        D = (self.K_().union_domain(other.K_()))[0]
                        PD = Poly(D)
                        [d,r] = PD._div_res_(self >> PD,other >> PD)
                        return r
                    except:
                        try:
                            return other.__rtruediv__(self)
                        except:
                            raise TypeError("I do not know how to opperate this elements")
            elif self.domain().subdomain().is_sub_domain(other.domain()):
                aux = Monomial_type(other >> self.K_(), self._variables[:],[0]*len(self._variables), self._structure)
                [d,r] = self._structure._div_res(self,aux)
                return d
            if other.domain().is_sub_domain(self.domain().subdomain()):
                D = self.domain()._create_Poly_Structure_(other.domain())
                aux = []
                for i in self._value:
                    aux.append(Monomial_type(i._value >> D.K_(),i._variables[:],i._index[:],D))
                aux = Poly_type(aux,self._variables[:],D,self._degree, self._mdegree)
                [d,r] = other._structure._div_res(aux,other>>D)
                return d 
            else:
                try:
                    D = (self.K_().union_domain(other.K_()))[0]
                    PD = Poly(D)
                    [d,r] = PD._div_res_(self >> PD,other >> PD)
                    return r
                except:
                    try:
                        return other.__rtruediv__(self)
                    except:
                        raise TypeError("I do not know how to opperate this elements")
        if isinstance(other, int):
            [d,r] = self._structure._div_res_(self,other)
            return d
        if isinstance(other, float):
            return self.__floordiv__((1>>self.domain())*other)
        try:
            return other.__rtruediv__(self)
        except:
            raise TypeError("I do not know how to opperate this elements")
    
    def __mod__(self,other):
        return self.residual_poly(other)
    
    def residual_poly(self, other):
        if isinstance(other, int):
            [d,r] = self._structure._div_res_(self,other)
            return r
        if isinstance(other, float):
            return self.residual_poly((1>>self.domain())*other)
        if isinstance(other, Symbol_type):
            return self.residual_poly(Monomial_type(1>>self._K(),[other],[1],self._structure))
        if isinstance(other, Monomial_type) or isinstance(other, Poly_type):
            if self._structure == other._structure:
                [d,r] = self._structure._div_res_(self,other)
                return r
            elif self.K_().is_sub_domain(other.K_()):
                [d,r] = self._structure._div_res_(self,other >> self._structure)
                return r
            elif other.K_().is_sub_domain(self.K_()):
                [d,r] = other._structure._div_res_(self >> other._structure,other)
                return r
            else:
                D = (self.K_().union_domain(other.K_()))[0]
                PD = Poly(D)
                [d,r] = PD._div_res_(self >> PD,other >> PD)
                return r
        if isinstance(other, Base_type.Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                [d,r] = self._structure._div_res_(self, Monomial_type(other>>self.K_(),self._variables[:],[0]*len(self._variables), self._structure))
                return r
            if other._structure.is_sub_domain(self._structure):
                [d,r] = other._structure._div_res_(self,other)
                return r
            if other.domain().is_sub_domain(self._structure.subdomain()):
                D = self.domain()._create_Poly_Structure_(other.domain())
                [d,r] = D._div_res_(self >> D,other>>D)
                return r
            else:
                D = (self.K_().union_domain(other.K_()))[0]
                PD = Poly(D)
                [d,r] = PD._div_res_(self >> PD,other >> PD)
                return r
        try:
            return other.__rtruediv__(self)
        except:
            raise TypeError("I do not know how to opperate this elements")   
    
    def _erase_var_order(self):
        aux = []
        for i in self._value:
            aux.append(i._erase_var_order())
        return Poly_type(aux,self._variables[:],self._structure, self._degree, self._mdegree)
    
    def _set_var_order(self, varOrd):
        if len(varOrd) == self._variables:
            aux = []
            for i in self._value:
                aux.append(Monomial_type(i._value,i._variables[:],i._index[:],i.domain(),varOrd))
            aord = _get_own_order_monomials_(aux[:])
            return Poly_type(aux,self._variables[:],self.domain(),self._degree, self._mdegree,aord,varOrd)
        else:
            raise TypeError("Wrong nuber of variables")
        
    
    def _get_var_order(self):
        if self._varOrder == None:
            return []
        else:
            return self._varOrder[:]
     
##Poly domain class and functions
def _transform_to_global_variable_order(p):
    if isinstance(p, Monomial_type):
        return Monomial_type(p._value, p._variables[:], p._index[:], p.domain())
    elif isinstance(p, Poly_type):
        aux = []
        for m in p._value:
            aux.append(Monomial_type(m._value, m._variables[:], m._index[:], m.domain()))
        return Poly_type(aux,p._variables[:],p.domain(), p._degree, p._mdegree)
    else:
        return p
                    
class Poly(Base.Base):
    def __init__(self, base, name = ''):
        self._ring = True
        self._group = True
        self._field = -1
        self._poly_structure = True
        self._characteristic = base.characteristic()
        self._cardinal = 'Infinity'
        self._subDomain = base
        if base._poly_structure:
            self._subDomain = base.subdomain()
        self._base_Domain = base.base_domain()
        if name == '':
            self._name = 'Pol['+self._subDomain._print_()+']'
            self._printdefault = True
        else:
            self._name = name
            self._printdefault = False
    
    def cp(self):
        return Poly(self._subDomain, self._name)     
    
    def variable(self, symb = 'x'):
        return _get_symbol_(symb)
    
    def variables(self, symb = ['x']):
        V = []
        for i in range(len(symb)):
            V.append(_get_symbol_(symb[i]))
        return V
    
    def characteristic(self):
        return self._characteristic
    
    def cardinal(self):
        return self._cardinal
        
    def reduce(self,value):
        if not isinstance(value, list):
            raise("Value must be a list of elements of " + self._subDomain._name)
        res = []        
        for i in value:
            res.append(i>>self.subdomain())
        return res
    
    def element_monomial(self,value,variables,index,varOrder = None):
        return Monomial_type(value >> self.subdomain(), variables, index, self, varOrder)
        
    def element(self,monomials,variables,domain,degree, total_degree,order = None, varOrder = None):
        if len(monomials) == 1:
            if len(variables) == 0:
                return Monomial_type(monomials[0]._value >> self.subdomain(), [], [], self, varOrder)
            elif len(variables) == 1:
                if total_degree == 1:
                    return variables[0].cp()
                else:
                    return Monomial_type(monomials[0]._value >> self.subdomain(), variables[:], monomials[0]._index[:], self, varOrder)
            else:
                return Monomial_type(monomials[0]._value >> self.subdomain(), variables[:], monomials[0]._index[:], self, varOrder)     
        else:
            aux = []
            for i in range(len(monomials)):
                aux.append(Monomial_type(monomials[i]._value >> self.subdomain(), variables[:], monomials[i]._index[:], self, varOrder))
            return Poly_type(aux,variables,domain,degree,total_degree, order, varOrder)
        
    def element_list(self,values,variable):
        if len(values) == 1:
            x = Symbol_type(variable)
            return Monomial_type(values[0]>>self.subdomain(),[x],[0],self)
        else:
            aux = []
            x = Symbol_type(variable)
            for i in range(len(values)):
                v = values[i]>>self.subdomain()
                if v != 0:
                    aux.append(Monomial_type(v,[x],[len(values)-1-i],self))
            if len(aux) == 0:
                return Monomial_type(0>>self.subdomain(),[x],[0],self)
            elif len(aux) == 1:
                return aux[0]
            else:
                return Poly_type(aux,[x],self,aux[0]._degree,aux[0]._degree)
        
    def _print_(self):
        return self._name
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._name + ' :: Ring'
        else:
            return self._name
    
    def __repr__(self):
        return self.__str__()
    
    def _transform_to_same_variables_varOrder_(self,a,b):
        if isinstance(a, Poly_type):
            if isinstance(b, Poly_type):
                if a.domain() == b.domain():
                    if a._variables == b._variables:
                        if a._varOrder == b._varOrder:
                            return a.cp(),b.cp()
                        else:
                            p0 = _transform_to_global_variable_order(a)
                            p1 = _transform_to_global_variable_order(b)
                            return p0,p1
                    else:
                        return a.add_variables(b._variables[:]), b.add_variables(a._variables[:])
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
            elif isinstance(b, Monomial_type):
                if a.domain() == b.domain():
                    if a._variables == b._variables:
                        if a._varOrder == b._varOrder:
                            return a.cp(),b.cp()
                        else:
                            p0 = _transform_to_global_variable_order(a)
                            p1 = _transform_to_global_variable_order(b)
                            return p0,p1
                    else:
                        return a.add_variables(b._variables[:]), b.add_variables(a._variables[:])
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
            elif isinstance(b, Symbol_type):
                if a._variables == [b]:
                    return a.cp(), Monomial_type(1>>self.subdomain(),[b],[1],self, None)
                else:
                    p0 = _transform_to_global_variable_order(a)
                    p0 = p0.add_variable(b)
                    p1 = Monomial_type(1>>self.subdomain(),[b],[1],self, None)
                    p1.add_variables(a._variables)
                    return p0,p1
            elif isinstance(b,int):
                if a._varOrder == None:
                    return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                else:
                    return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self,a._varOrder[:])
            elif isinstance(b,float):
                return self._transform_to_same_variables_varOrder_(a,b>>QQ(ZZ()))
            else:
                if self.subdomain().is_sub_domain(b.domain()) and a.domain() == self:
                    if a._varOrder == None:
                        return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                    else:
                        return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self,a._varOrder[:])
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
        elif isinstance(a, Monomial_type):
            if isinstance(b, Poly_type):
                if a.domain() == b.domain():
                    if a._variables == b._variables:
                        if a._varOrder == b._varOrder:
                            return a.cp(),b.cp()
                        else:
                            p0 = _transform_to_global_variable_order(a)
                            p1 = _transform_to_global_variable_order(b)
                            return p0,p1
                    else:
                        return a.add_variables(b._variables[:]), b.add_variables(a._variables[:])
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
            elif isinstance(b, Monomial_type):
                if a.domain() == b.domain():
                    if a._variables == b._variables:
                        if a._varOrder == b._varOrder:
                            return a.cp(),b.cp()
                        else:
                            p0 = _transform_to_global_variable_order(a)
                            p1 = _transform_to_global_variable_order(b)
                            return p0,p1
                    else:
                        return a.add_variables(b._variables[:]), b.add_variables(a._variables[:])
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
            elif isinstance(b, Symbol_type):
                if a._variables == [b]:
                    return a.cp(), Monomial_type(1>>self.subdomain(),[b],[1],self, None)
                else:
                    p0 = _transform_to_global_variable_order(a)
                    p0 = p0.add_variable(b)
                    p1 = Monomial_type(1>>self.subdomain(),[b],[1],self, None)
                    p1.add_variables(a._variables)
                    return p0,p1
            elif isinstance(b,int):
                if a._varOrder == None:
                    return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                else:
                    return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self,a._varOrder[:])
            elif isinstance(b,float):
                return self._transform_to_same_variables_varOrder_(a,b>>QQ(ZZ()))
            else:
                if self.subdomain().is_sub_domain(b.domain()) and self == a.domain():
                    if a._varOrder == None:
                        return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                    else:
                        return a.cp(), Monomial_type(b>>self.subdomain(),a._variables[:],[0]*len(a._variables),self,a._varOrder[:])
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
        elif isinstance(a, Symbol_type):
            if isinstance(b, Poly_type):
                if b._variables == [a]:
                    return Monomial_type(1>>self.subdomain(),[a],[1],self, None), b.cp()
                else:
                    p0 = Monomial_type(1>>self.subdomain(),[a],[1],self, None)
                    p0.add_variables(b._variables)
                    p1 = _transform_to_global_variable_order(b)
                    p1 = p1.add_variable(a)
                    return p0,p1
            elif isinstance(b, Monomial_type):
                if b._variables == [a]:
                    return Monomial_type(1>>self.subdomain(),[a],[1],self, None), b.cp()
                else:
                    p0 = Monomial_type(1>>self.subdomain(),[a],[1],self, None)
                    p0.add_variables(b._variables)
                    p1 = _transform_to_global_variable_order(b)
                    p1 = p1.add_variable(a)
                    return p0,p1
            elif isinstance(b, Symbol_type):
                if a == b:
                    return Monomial_type(1>>self.subdomain(),[a],[1],self), Monomial_type(1>>self.subdomain(),[b],[1],self)
                else:
                    p0 = Monomial_type(1>>self.subdomain(),[a],[1],self)
                    p1 = Monomial_type(1>>self.subdomain(),[b],[1],self)
                    p0 = p0.add_variable(b)
                    p1 = p1.add_variable(a)
                    return p0,p1
            elif isinstance(b,int):
                return Monomial_type(1>>self.subdomain(),[a],[1],self), Monomial_type(b>>self.subdomain(),[a],[0],self)
            elif isinstance(b,float):
                return self._transform_to_same_variables_varOrder_(a,b>>QQ(ZZ()))
            else:
                if self==b.domain():
                    return Monomial_type(1>>self.subdomain(),[a],[1],self), Monomial_type(b>>self.subdomain(),[a],[0],self)
                else:
                    D = Poly(b.domain())
                    return Monomial_type(1>>D.subdomain(),[a],[1],D), Monomial_type(b>>D.subdomain(),[a],[0],D)
        elif isinstance(a,int):
            if isinstance(b, Poly_type):
                if b.domain() == self:
                    if b._varOrder == None:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self), b.cp()
                    else:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self,b._varOrder[:]), b.cp()
                else:
                    if b._varOrder == None:
                        return Monomial_type(a>>b.domain().subdomain(),b._variables[:],[0]*len(b._variables),b.domain()), b.cp()
                    else:
                        return Monomial_type(a>>b.domain().subdomain(),b._variables[:],[0]*len(b._variables),b.domain(),b._varOrder[:]), b.cp()
            elif isinstance(b, Monomial_type):
                if b.domain() == self:
                    if b._varOrder == None:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self), b.cp()
                    else:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self,b._varOrder[:]), b.cp()
                else:
                    if b._varOrder == None:
                        return Monomial_type(a>>b.domain().subdomain(),b._variables[:],[0]*len(b._variables),b.domain()), b.cp()
                    else:
                        return Monomial_type(a>>b.domain().subdomain(),b._variables[:],[0]*len(b._variables),b.domain(),b._varOrder[:]), b.cp()
            elif isinstance(b, Symbol_type):
                if self.subdomain() == ZZ():
                    return Monomial_type(a>>self.subdomain(),[b],[0],self), Monomial_type(1>>self.subdomain(),[b],[1],self)
                else:
                    D = Poly(ZZ())
                    return Monomial_type(a,[b],[0],D), Monomial_type(1,[b],[1],D)
            elif isinstance(b,int):
                D = Poly(ZZ())
                return Monomial_type(a,[],[],D), Monomial_type(b,[],[],D)
            elif isinstance(b,float):
                D = Poly(QQ(ZZ()))
                return Monomial_type(a>>D.subdomain(),[],[],D), Monomial_type(b>>D.subdomain(),[],[],D)
            else:
                D = Poly(b.domain())
                return Monomial_type(a>>D.subdomain(),[],[],D), Monomial_type(b,[],[],D)
        elif isinstance(a,float):
            return self._transform_to_same_variables_varOrder_(a>>QQ(ZZ()),b)
        else:
            if isinstance(b, Poly_type):
                if b.domain() == self and self.subdomain().is_sub_domain(a.domain()):
                    if b._varOrder == None:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self), b.cp()
                    else:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self,b._varOrder[:]), b.cp()
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
            elif isinstance(b, Monomial_type):
                if b.domain() == self and self.subdomain().is_sub_domain(a.domain()):
                    if b._varOrder == None:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self), b.cp()
                    else:
                        return Monomial_type(a>>self.subdomain(),b._variables[:],[0]*len(b._variables),self,b._varOrder[:]), b.cp()
                else:
                    p0 = a._push_to_union(b.domain())
                    p1 = b._push_to_union(a.domain())
                    return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
            elif isinstance(b, Symbol_type):
                if self.is_sub_domain(a.domain()):
                    return Monomial_type(a>>self.subdomain(),[b],[0],self), Monomial_type(1>>self.subdomain(),[b],[1],self)
                else:
                    D = Poly(a.domain())
                    return Monomial_type(a,[b],[0],D), Monomial_type(1>>D,[b],[1],D)
            elif isinstance(b,int):
                D = Poly(a.domain())
                return Monomial_type(a,[],[],D), Monomial_type(b>>D,[],[],D)
            elif isinstance(b,float):
                return self._transform_to_same_variables_varOrder_(a,b>>QQ(ZZ()))
            else:
                p0 = a._push_to_union(b.domain())
                p1 = b._push_to_union(a.domain())
                return p0.domain()._transform_to_same_variables_varOrder_(p0,p1)
            
    def is_field(self):
        return False
    
    def is_invertible(self, a):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return self.is_invertible(aux)
        if isinstance(a, int):
            return self._subDomain.is_invertible(a>>self._subDomain)
        if isinstance(a, Poly_type):
            return False
        if isinstance(a, Monomial_type):
            if a.mdegree() > 0:
                return False
            else:
                return self._subDomain.is_invertible(a)
        if isinstance(a, Symbol_type):
            return False
        if isinstance(a,Base_type.Base_type):
            if self.is_sub_element(a):
                return self._subDomain.is_invertible(a>>self._subDomain)
        return False
        
    def _inv_(self, a):
        if not self.is_invertible(a):
            return 'inv error: This number is not invertible'
        if (self.is_sub_element(a))and(not self.is_element(a)): 
            aux = a>>self
            return self._inv_(aux)
        if a._varOrder == None:
            return Monomial_type(self._subDomain._inv_(a._value), [], [], self) 
        else:
            return Monomial_type(self._subDomain._inv_(a._value), a._variables[:], a._index[:], self, a._varOrder[:])
    
    def is_leading_coeff_inv(self, a):
        aux = a.leading_coeff()
        return self.is_invertible(aux)
    
    def gcd(self,a,b):
        r0,r1 = self._transform_to_same_variables_varOrder_(a,b)
        ##TO DO: Pensar què fer en el cas varOrder != None
        r0 = r0._erase_var_order()
        r1 = r1._erase_var_order()
        ##
        if len(r0._variables) > 1:
            raise NotImplementedError
        else:
            if r0.domain().subdomain().is_field:
                [a1,a2,a3]=r0.domain()._extended_euclidean_algorithm_core_(r0,r1)
                return a1
            else:
                raise NotImplementedError
                
        
    def extended_euclidean_algorithm(self, a, b):
        r0,r1 = self._transform_to_same_variables_varOrder_(a,b)
        ##TO DO: Pensar què fer en el cas varOrder != None
        r0 = r0._erase_var_order()
        r1 = r1._erase_var_order()
        if len(r0._variables) > 1:
            return 'extended_euclidean_algorithm error: Not possible to compute the Bézout of multivariate polynomials'
        ##
        return r0.domain()._extended_euclidean_algorithm_core_(r0,r1)
    
    def _extended_euclidean_algorithm_core_(self,a,b):
        if a._varOrder == None:
            zero = Monomial_type(0>>self.subdomain(),a._variables[:], [0]*len(a._variables), self)
            one = Monomial_type(1>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
            
            r0 = a.cp()
            r1 = b.cp()
            
            s0 = one.cp()
            s1 = zero.cp()
            t0 = zero.cp()
            t1 = one.cp()
            
            while(r1 != zero):
                [dv, rs] = self._div_res_core_(r0, r1)
                s2 = s1
                t2 = t1
                r0 = r1
                r1 = rs
                s1 = self._sub_core_(s0, self._mul_core_(dv,s1))
                t1 = self._sub_core_(t0, self._mul_core_(dv,t1))
                s0 = s2
                t0 = t2
            if isinstance(r0,Monomial_type):
                lc = r0._value
            elif isinstance(r0,Poly_type):
                lc = r0._value[0]._value
            if lc != 1:
                t0=t0/lc
                s0=s0/lc
                r0=r0/lc
            r0 = r0
            s0 = s0
            t0 = t0
            return [r0, s0, t0]
                
        else:
            raise NotImplementedError 

    def _add_(self,a,b):
        r0,r1 = self._transform_to_same_variables_varOrder_(a,b)
        return r0.domain()._add_core_(r0,r1)

    def _add_core_(self,a,b):
        if isinstance(a,Monomial_type):
            if isinstance(b,Monomial_type):
                if a == 0:
                    return b.cp()
                if b == 0:
                    return a.cp()
                if a._varOrder == None:
                    o = _order_monomials_same_variables_(a,b)
                    if o == 0:
                        s = a._value+b._value
                        if s == 0:
                            return Monomial_type(s, a._variables[:], [0]*len(a._variables), a.domain())
                        else:
                            return Monomial_type(s, a._variables[:], a._index[:], a.domain())
                    elif o == -1:
                        aux = [b.cp(), a.cp()]
                        return Poly_type(aux,a._variables[:],a.domain(),b._degree, max(a._totIndex,b._totIndex))
                    else:
                        aux = [a.cp(), b.cp()]
                        return Poly_type(aux,a._variables[:],a.domain(),a._degree, max(a._totIndex,b._totIndex))
                else:
                    o = _own_order_monomials_(a,b)
                    if o == 0:
                        s = a._value+b._value
                        if s == 0:
                            return Monomial_type(s, a._variables[:], [0]*len(a._variables), a.domain())
                        else:
                            return Monomial_type(s, a._variables[:], a._index[:], a.domain(), a.domain(),a._varOrder[:])
                    elif o == -1:
                        aux = [b.cp(), a.cp()]
                        ordr = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],a.domain(),b._degree, max(a._totIndex,b._totIndex),ordr,a._varOrder[:])
                    else:
                        aux = [a.cp(), b.cp()]
                        ordr = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],a.domain(),a._degree, max(a._totIndex,b._totIndex),ordr,a._varOrder[:])
            elif isinstance(b, Poly_type):
                if a == 0:
                    return b.cp()
                if a._varOrder == None:
                    for i in range(len(b._value)):
                        o = _order_monomials_same_variables_(a,b._value[i])
                        if o == 0:
                            s = a._value+b._value[i]._value
                            if s == 0:
                                if len(b._value) == 2:
                                    if i == 0:
                                        return b._value[1].cp()
                                    else:
                                        return b._value[0].cp()
                                else:
                                    aux = b._value[:i]+b._value[i+1:]
                                    if b._value[i]._totIndex == b._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == b._mdegree:
                                                md = b._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = b._mdegree
                                    return Poly_type(aux,b._variables[:],b.domain(),aux[0]._degree,md)
                            else:
                                aux = Monomial_type(s,a._variables[:],a._index[:],a.domain())
                                return Poly_type(b._value[:i] + [aux] + b._value[i+1:],b._variables[:],b.domain(),b._degree,b._mdegree)
                        elif o == 1:
                            return Poly_type(b._value[:i] + [a] + b._value[i:], b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree))
                    return Poly_type(b._value[:] + [a], b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree))
                else:
                    for i in range(len(b._value)):
                        o = _order_monomials_same_variables_(a,b._value[i])
                        if o == 0:
                            s = a._value+b._value[i]._value
                            if s == 0:
                                if len(b._value) == 2:
                                    if i == 0:
                                        return b._value[1].cp()
                                    else:
                                        return b._value[0].cp()
                                else:
                                    aux = b._value[:i]+b._value[i+1:]
                                    if b._value[i]._totIndex == b._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == b._mdegree:
                                                md = b._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = b._mdegree
                                    vord = b._ord[:]
                                    nbo = len(b._ord)
                                    dbo = 0
                                    j = 0
                                    while(j < nbo - dbo):
                                        if vord[j] == i:
                                            vord = vord[:j] + vord[j+1:]
                                            j -= 1
                                        elif i < vord[j]:
                                            vord[j] -= 1
                                        j += 1
                                    return Poly_type(aux,b._variables[:],b.domain(),aux[0]._degree,md,vord,b._varOrder[:])
                            else:
                                aux = Monomial_type(s,a._variables[:],a._index[:],a.domain())
                                return Poly_type(b._value[:i] + [aux] + b._value[i+1:],b._variables[:],b.domain(),b._degree,b._mdegree, b._ord[:], b._varOrder[:])
                        elif o == 1:
                            aux = b._value[:i] + [a] + b._value[i:]
                            vord = _get_own_order_monomials_(aux[:])
                            return Poly_type(aux, b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree),vord,b._varOrder[:])
                    aux = b._value[:] + [a]
                    vord = _get_own_order_monomials_(aux[:])
                    return Poly_type(aux, b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree), vord, b._varOrder[:])
        elif isinstance(a,Poly_type):
            if isinstance(b, Monomial_type):
                if b == 0:
                    return a.cp()
                if a._varOrder == None:
                    for i in range(len(a._value)):
                        o = _order_monomials_same_variables_(b,a._value[i])
                        if o == 0:
                            s = a._value[i]._value+b._value
                            if s == 0:
                                if len(a._value) == 2:
                                    if i == 0:
                                        return a._value[1].cp()
                                    else:
                                        return a._value[0].cp()
                                else:
                                    aux = a._value[:i]+a._value[i+1:]
                                    if a._value[i]._totIndex == a._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == b._mdegree:
                                                md = b._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = b._totIndex
                                    return Poly_type(aux,a._variables[:],a.domain(),aux[0]._degree,md)
                            else:
                                aux = Monomial_type(s,b._variables[:],b._index[:],b.domain())
                                return Poly_type(a._value[:i] + [aux] + a._value[i+1:],a._variables[:],a.domain(),a._degree,a._mdegree)
                        elif o == 1:
                            return Poly_type(a._value[:i] + [b] + a._value[i:], a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree))
                    return Poly_type(a._value[:] + [b], a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree))
                else:
                    for i in range(len(a._value)):
                        o = _order_monomials_same_variables_(b,a._value[i])
                        if o == 0:
                            s = a._value[i]._value+b._value
                            if s == 0:
                                if len(a._value) == 2:
                                    if i == 0:
                                        return a._value[1].cp()
                                    else:
                                        return a._value[0].cp()
                                else:
                                    aux = a._value[:i]+a._value[i+1:]
                                    if a._value[i]._totIndex == a._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == a._mdegree:
                                                md = a._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = a._totIndex
                                    vord = a._ord[:]
                                    nao = len(a._ord)
                                    dao = 0
                                    j = 0
                                    while(j < nao - dao):
                                        if vord[j] == i:
                                            vord = vord[:j] + vord[j+1:]
                                            j -= 1
                                        elif i < vord[j]:
                                            vord[j] -= 1
                                        j += 1
                                    return Poly_type(aux,a._variables[:],a.domain(),aux[0]._degree,md,vord,a._varOrder[:])
                            else:
                                aux = Monomial_type(s,b._variables[:],b._index[:],b.domain())
                                return Poly_type(a._value[:i] + [aux] + a._value[i+1:],a._variables[:],a.domain(),a._degree,a._mdegree, a._ord[:], a._varOrder[:])
                        elif o == 1:
                            aux = a._value[:i] + [b] + a._value[i:]
                            vord = _get_own_order_monomials_(aux[:])
                            return Poly_type(aux, a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree),vord,a._varOrder[:])
                    aux = a._value[:] + [b]
                    vord = _get_own_order_monomials_(aux[:])
                    return Poly_type(aux, a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree), vord, a._varOrder[:])
            elif isinstance(b,Poly_type):
                if a._varOrder == None:
                    aux = a._value[:]
                    i = 0
                    j = 0
                    n = len(aux)
                    af = 0
                    while(i < n + af):
                        o = _order_monomials_same_variables_(b._value[j], aux[i])
                        if o == 0:
                            s = aux[i]._value + b._value[j]._value
                            if s == 0:
                                aux = aux[:i] + aux[i+1:]
                                af -= 1
                            else:
                                aux = aux[:i] + [Monomial_type(s,a._variables[:],b._value[j]._index[:],a.domain())] + aux[i+1:]
                                i += 1
                            j += 1
                        elif o == 1:
                            aux = aux[:i] + [b._value[j]] + aux[i:]
                            j += 1
                            i += 1
                            af += 1
                        else:
                            i += 1
                        if (j >= len(b._value)):
                            break
                    if j < len(b._value):
                        aux = aux + b._value[j:]
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(), a._variables[:],[0]*len(a._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        return Poly_type(aux,a._variables[:],self,max(a._degree,b._degree),max(a._mdegree,b._mdegree))
                else:
                    aux = a._value[:]
                    i = 0
                    j = 0
                    n = len(aux)
                    af = 0
                    while(i < n + af):
                        o = _order_monomials_same_variables_(b._value[j], aux[i])
                        if o == 0:
                            s = aux[i]._value + b._value[j]._value
                            if s == 0:
                                aux = aux[:i] + aux[i+1:]
                                af -= 1
                            else:
                                aux = aux[:i] + [Monomial_type(s,a._variables[:],b._value[j]._index[:],a.domain(),a._varOrder[:])] + aux[i+1:]
                                i += 1
                            j += 1
                        elif o == 1:
                            aux = aux[:i] + [b._value[j]] + aux[i:]
                            j += 1
                            i += 1
                            af += 1
                        else:
                            i += 1
                        if (j >= len(b._value)):
                            break
                    if j < len(b._value):
                        aux = aux + b._value[j:]
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(), a._variables[:],[0]*len(a._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        vord = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],self,max(a._degree,b._degree),max(a._mdegree,b._mdegree),vord,a._varOrder[:])
                    
    def _sub_(self,a,b):
        r0,r1 = self._transform_to_same_variables_varOrder_(a,b)
        return r0.domain()._sub_core_(r0,r1)
    
    def _sub_core_(self,a,b):
        if isinstance(a,Monomial_type):
            if isinstance(b,Monomial_type):
                if a == 0:
                    return -b
                if b == 0:
                    return a.cp()
                if a._varOrder == None:
                    o = _order_monomials_same_variables_(a,b)
                    if o == 0:
                        s = a._value-b._value
                        if s == 0:
                            return Monomial_type(s, a._variables[:], [0]*len(a._variables), a.domain())
                        else:
                            return Monomial_type(s, a._variables[:], a._index[:], a.domain())
                    elif o == -1:
                        aux = [Monomial_type(-b._value,b._variables[:],b._index[:],b.domain()), a.cp()]
                        return Poly_type(aux,a._variables[:],a.domain(),b._degree, max(a._totIndex,b._totIndex))
                    else:
                        aux = [a.cp(), Monomial_type(-b._value,b._variables[:],b._index[:],b.domain())]
                        return Poly_type(aux,a._variables[:],a.domain(),a._degree, max(a._totIndex,b._totIndex))
                else:
                    o = _own_order_monomials_(a,b)
                    if o == 0:
                        s = a._value-b._value
                        if s == 0:
                            return Monomial_type(s, a._variables[:], [0]*len(a._variables), a.domain())
                        else:
                            return Monomial_type(s, a._variables[:], a._index[:], a.domain(), a.domain(),a._varOrder[:])
                    elif o == -1:
                        aux = [Monomial_type(-b._value,b._variables[:],b._index[:],b.domain(),b._varOrder[:]), a.cp()]
                        ordr = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],a.domain(),b._degree, max(a._totIndex,b._totIndex),ordr,a._varOrder[:])
                    else:
                        aux = [a.cp(), Monomial_type(-b._value,b._variables[:],b._index[:],b.domain(), b._varOrder[:])]
                        ordr = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],a.domain(),a._degree, max(a._totIndex,b._totIndex),ordr,a._varOrder[:])
            elif isinstance(b, Poly_type):
                if a == 0:
                    return -b
                if a._varOrder == None:
                    bm = [Monomial_type(-i._value, i._variables, i._index[:],i.domain()) for i in b._value]
                    for i in range(len(bm)):
                        o = _order_monomials_same_variables_(a,bm[i])
                        if o == 0:
                            s = a._value+bm[i]._value
                            if s == 0:
                                if len(bm) == 2:
                                    if i == 0:
                                        return bm[1]
                                    else:
                                        return bm[0]
                                else:
                                    aux = bm[:i]+bm[i+1:]
                                    if bm[i]._totIndex == b._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == b._mdegree:
                                                md = b._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = b._mdegree
                                    return Poly_type(aux,b._variables[:],b.domain(),aux[0]._degree,md)
                            else:
                                aux = Monomial_type(s,a._variables[:],a._index[:],a.domain())
                                return Poly_type(bm[:i] + [aux] + bm[i+1:],b._variables[:],b.domain(),b._degree,b._mdegree)
                        elif o == 1:
                            return Poly_type(bm[:i] + [a] + bm[i:], b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree))
                    return Poly_type(bm[:] + [a], b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree))
                else:
                    bm = [Monomial_type(-i._value, i._variables, i._index[:],i.domain()) for i in b._value]
                    for i in range(len(bm)):
                        o = _order_monomials_same_variables_(a,bm[i])
                        if o == 0:
                            s = a._value+bm[i]._value
                            if s == 0:
                                if len(bm) == 2:
                                    if i == 0:
                                        return bm[1]
                                    else:
                                        return bm[0]
                                else:
                                    aux = bm[:i]+bm[i+1:]
                                    if bm[i]._totIndex == b._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == b._mdegree:
                                                md = b._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = b._mdegree
                                    vord = b._ord[:]
                                    nbo = len(b._ord)
                                    dbo = 0
                                    j = 0
                                    while(j < nbo - dbo):
                                        if vord[j] == i:
                                            vord = vord[:j] + vord[j+1:]
                                            j -= 1
                                        elif i < vord[j]:
                                            vord[j] -= 1
                                        j += 1
                                    return Poly_type(aux,b._variables[:],b.domain(),aux[0]._degree,md,vord,b._varOrder[:])
                            else:
                                aux = Monomial_type(s,a._variables[:],a._index[:],a.domain())
                                return Poly_type(bm[:i] + [aux] + bm[i+1:],b._variables[:],b.domain(),b._degree,b._mdegree, b._ord[:], b._varOrder[:])
                        elif o == 1:
                            aux = bm[:i] + [a] + bm[i:]
                            vord = _get_own_order_monomials_(aux[:])
                            return Poly_type(aux, b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree),vord,b._varOrder[:])
                    aux = bm[:] + [a]
                    vord = _get_own_order_monomials_(aux[:])
                    return Poly_type(aux, b._variables[:],b.domain(),max(a._degree,b._degree), max(a._totIndex,b._mdegree), vord, b._varOrder[:])
        elif isinstance(a,Poly_type):
            if isinstance(b, Monomial_type):
                if b == 0:
                    return a.cp()
                if a._varOrder == None:
                    for i in range(len(a._value)):
                        o = _order_monomials_same_variables_(b,a._value[i])
                        if o == 0:
                            s = a._value[i]._value-b._value
                            if s == 0:
                                if len(a._value) == 2:
                                    if i == 0:
                                        return a._value[1].cp()
                                    else:
                                        return a._value[0].cp()
                                else:
                                    aux = a._value[:i]+a._value[i+1:]
                                    if a._value[i]._totIndex == a._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == a._mdegree:
                                                md = a._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = a._mdegree
                                    return Poly_type(aux,a._variables[:],a.domain(),aux[0]._degree,md)
                            else:
                                aux = Monomial_type(s,b._variables[:],b._index[:],b.domain())
                                return Poly_type(a._value[:i] + [aux] + a._value[i+1:],a._variables[:],a.domain(),a._degree,a._mdegree)
                        elif o == 1:
                            return Poly_type(a._value[:i] + [Monomial_type(-b._value,b._variables[:],b._index[:],b.domain())] + a._value[i:], a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree))
                    return Poly_type(a._value[:] + [Monomial_type(-b._value,b._variables[:],b._index[:],b.domain())], a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree))
                else:
                    for i in range(len(a._value)):
                        o = _order_monomials_same_variables_(b,a._value[i])
                        if o == 0:
                            s = a._value[i]._value-b._value
                            if s == 0:
                                if len(a._value) == 2:
                                    if i == 0:
                                        return a._value[1].cp()
                                    else:
                                        return a._value[0].cp()
                                else:
                                    aux = a._value[:i]+a._value[i+1:]
                                    if a._value[i]._totIndex == a._mdegree:
                                        md = 0
                                        for j in aux:
                                            if j._totIndex == a._mdegree:
                                                md = a._mdegree
                                                break
                                            elif j._totIndex > md:
                                                md = j._totIndex
                                    else:
                                        md = a._mdegree
                                    vord = a._ord[:]
                                    nao = len(a._ord)
                                    dao = 0
                                    j = 0
                                    while(j < nao - dao):
                                        if vord[j] == i:
                                            vord = vord[:j] + vord[j+1:]
                                            j -= 1
                                        elif i < vord[j]:
                                            vord[j] -= 1
                                        j += 1
                                    return Poly_type(aux,a._variables[:],a.domain(),aux[0]._degree,md,vord,a._varOrder[:])
                            else:
                                aux = Monomial_type(s,b._variables[:],b._index[:],b.domain())
                                return Poly_type(a._value[:i] + [aux] + a._value[i+1:],a._variables[:],a.domain(),a._degree,a._mdegree, a._ord[:], a._varOrder[:])
                        elif o == 1:
                            aux = a._value[:i] + [Monomial_type(-b._value,b._variables[:],b._index[:],b.domain(),b._varOrder[:])] + a._value[i:]
                            vord = _get_own_order_monomials_(aux[:])
                            return Poly_type(aux, a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree),vord,a._varOrder[:])
                    aux = a._value[:] + [Monomial_type(-b._value,b._variables[:],b._index[:],b.domain(),b._varOrder[:])]
                    vord = _get_own_order_monomials_(aux[:])
                    return Poly_type(aux, a._variables[:],a.domain(),max(a._degree,b._degree), max(b._totIndex,a._mdegree), vord, a._varOrder[:])
            elif isinstance(b,Poly_type):
                if a._varOrder == None:
                    bm = [Monomial_type(-i._value, i._variables, i._index[:],i.domain()) for i in b._value]
                    aux = a._value[:]
                    i = 0
                    j = 0
                    n = len(aux)
                    af = 0
                    while(i < n + af):
                        o = _order_monomials_same_variables_(bm[j], aux[i])
                        if o == 0:
                            s = aux[i]._value + bm[j]._value
                            if s == 0:
                                aux = aux[:i] + aux[i+1:]
                                af -= 1
                            else:
                                aux = aux[:i] + [Monomial_type(s,a._variables[:],b._value[j]._index[:],a.domain())] + aux[i+1:]
                                i += 1
                            j += 1
                        elif o == 1:
                            aux = aux[:i] + [bm[j]] + aux[i:]
                            j += 1
                            i += 1
                            af += 1
                        else:
                            i += 1
                        if (j >= len(bm)):
                            break
                    if j < len(bm):
                        aux = aux + bm[j:]
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(), a._variables[:],[0]*len(a._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        return Poly_type(aux,a._variables[:],self,max(a._degree,b._degree),max(a._mdegree,b._mdegree))
                else:
                    aux = a._value[:]
                    bm = [Monomial_type(-i._value, i._variables, i._index[:],i.domain()) for i in b._value]
                    i = 0
                    j = 0
                    n = len(aux)
                    af = 0
                    while(i < n + af):
                        o = _order_monomials_same_variables_(bm[j], aux[i])
                        if o == 0:
                            s = aux[i]._value + bm[j]._value
                            if s == 0:
                                aux = aux[:i] + aux[i+1:]
                                af -= 1
                            else:
                                aux = aux[:i] + [Monomial_type(s,a._variables[:],b._value[j]._index[:],a.domain(), a._varOrder[:])] + aux[i+1:]
                                i += 1
                            j += 1
                        elif o == 1:
                            aux = aux[:i] + [bm[j]] + aux[i:]
                            j += 1
                            i += 1
                            af += 1
                        else:
                            i += 1
                        if (j >= len(bm)):
                            break
                    if j < len(bm):
                        aux = aux + bm[j:]
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(), a._variables[:],[0]*len(a._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        vord = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],self,max(a._degree,b._degree),max(a._mdegree,b._mdegree),vord,a._varOrder[:])
                    
    def _mul_(self,a,b):
        r0,r1 = self._transform_to_same_variables_varOrder_(a,b)
        return r0.domain()._mul_core_(r0,r1)
                    
    def _mul_core_(self,a,b):
        if isinstance(a,Monomial_type):
            if isinstance(b,Monomial_type):
                if a == 0:
                    return a.cp()
                if b == 0:
                    return b.cp()
                if a._varOrder == None:
                    s = a._value*b._value
                    if s == 0:
                        return Monomial_type(s,a._variables[:], [0]*len(a._variables), self)
                    else:
                        index = [a._index[i] + b._index[i] for i in range(len(a._index))]
                        return Monomial_type(s,a._variables[:], index, self)
                else:
                    s = a._value*b._value
                    if s == 0:
                        return Monomial_type(s,a._variables[:], [0]*len(a._variables), self)
                    else:
                        index = [a._index[i] + b._index[i] for i in range(len(a._index))]
                        return Monomial_type(s,a._variables[:], index, self, a._varOrder[:])
            elif isinstance(b,Poly_type):
                if a == 0:
                    return a.cp()
                if a._varOrder == None:
                    aux = []
                    md = 0
                    for i in range(len(b._value)):
                        s = a._value*b._value[i]._value
                        if s != 0:
                            index = [a._index[j] + b._value[i]._index[j] for j in range(len(a._index))]
                            aux.append(Monomial_type(s,a._variables[:], index, self))
                            if aux[-1]._totIndex > md:
                                md = aux[-1]._totIndex
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(),a._variables[:], [0]*len(a._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        return Poly_type(aux,a._variables[:],self,aux[0]._degree, md)
                else:
                    aux = []
                    md = 0
                    for i in range(len(b._value)):
                        s = a._value*b._value[i]._value
                        if s != 0:
                            index = [a._index[j] + b._value[i]._index[j] for j in range(len(a._index))]
                            aux.append(Monomial_type(s,a._variables[:], index, self,a._varOrder[:]))
                            if aux[-1]._totIndex > md:
                                md = aux[-1]._totIndex
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(),a._variables[:], [0]*len(a._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        vord = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],self,aux[0]._degree, md, vord, b._varOrder[:])
        elif isinstance(a,Poly_type):
            if isinstance(b,Monomial_type):
                if b == 0:
                    return b.cp()
                if a._varOrder == None:
                    aux = []
                    md = 0
                    for i in range(len(a._value)):
                        s = a._value[i]._value*b._value
                        if s != 0:
                            index = [a._value[i]._index[j] + b._index[j] for j in range(len(b._index))]
                            aux.append(Monomial_type(s,b._variables[:], index, self))
                            if aux[-1]._totIndex > md:
                                md = aux[-1]._totIndex
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(),b._variables[:], [0]*len(b._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        return Poly_type(aux,a._variables[:],self,aux[0]._degree, md)
                else:
                    aux = []
                    md = 0
                    for i in range(len(a._value)):
                        s = a._value[i]._value*b._value
                        if s != 0:
                            index = [a._value[i]._index[j] + b._index[j] for j in range(len(b._index))]
                            aux.append(Monomial_type(s,b._variables[:], index, self,b._varOrder[:]))
                            if aux[-1]._totIndex > md:
                                md = aux[-1]._totIndex
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(),a._variables[:], [0]*len(a._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        vord = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],self,aux[0]._degree, md, vord, a._varOrder[:])
            elif isinstance(b,Poly_type):
                if a._varOrder == None:
                    aux = []
                    for i in range(len(a._value)):
                        ind = 0
                        for j in range(len(b._value)):
                            index = [a._value[i]._index[k]+b._value[j]._index[k] for k in range(len(a._value[i]._index))]
                            pos,fo = _find_index_monomials(index,aux,ind)
                            if fo:
                                s = aux[pos]._value + a._value[i]._value*b._value[j]._value
                                if s == 0:
                                    aux = aux[:pos] + aux[pos+1:]
                                else:
                                    aux[pos] = Monomial_type(s,a._variables[:],aux[pos]._index,self)
                            else:
                                s = a._value[i]._value*b._value[j]._value
                                if s != 0:
                                    aux = aux[:pos] + [Monomial_type(s,a._variables[:],index,self)] + aux[pos:]
                            ind = pos
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(),b._variables[:], [0]*len(b._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        l = a._mdegree + b._mdegree
                        md = 0
                        for i in aux:
                            if l == i._totIndex:
                                md = l
                                break
                            elif md < i._totIndex:
                                md = i._totIndex
                        return Poly_type(aux,a._variables[:],self,aux[0]._degree, md)
                else:
                    aux = []
                    for i in range(len(a._value)):
                        ind = 0
                        for j in range(len(b._value)):
                            index = [a._value[i]._index[k]+b._value[j]._index[k] for k in range(len(a._value[i]._index))]
                            pos,fo = _find_index_monomials(index,aux,ind)
                            if fo:
                                s = aux[pos]._value + a._value[i]._value*b._value[j]._value
                                if s == 0:
                                    aux = aux[:pos] + aux[pos+1:]
                                else:
                                    aux[pos] = Monomial_type(s,a._variables[:],aux[pos]._index,self,a._varOrder[:])
                            else:
                                s = a._value[i]._value*b._value[j]._value
                                if s != 0:
                                    aux = aux[:pos] + [Monomial_type(s,a._variables[:],index,self,a._varOrder[:])] + aux[pos:]
                            ind = pos
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(),b._variables[:], [0]*len(b._variables), self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        l = a._mdegree + b._mdegree
                        md = 0
                        for i in aux:
                            if l == i._totIndex:
                                md = l
                                break
                            elif md < i._totIndex:
                                md = i._totIndex
                        vord = _get_own_order_monomials_(aux[:])
                        return Poly_type(aux,a._variables[:],self,aux[0]._degree, md, a._varOrder[:],vord)
                    
                    
    def _div_res_(self, a, b):
        r0,r1 = self._transform_to_same_variables_varOrder_(a,b)
        ##TO DO: Pensar què fer en el cas varOrder != None
        r0 = r0._erase_var_order()
        r1 = r1._erase_var_order()
        ##
        return r0.domain()._div_res_core_(r0,r1)
    
    def _div_res_core_(self,a,b):
        if isinstance(a,Monomial_type):
            if isinstance(b,Monomial_type):
                if a._varOrder == None:
                    index = []
                    fdiv = True
                    for i in range(len(a._index)):
                        if b._index[i] <= a._index[i]:
                            index.append(a._index[i] - b._index[i])
                        else:
                            fdiv = False
                            break
                    if fdiv:
                        [d,r] = self.subdomain()._div_res_(a._value,b._value)
                        if r == 0:
                            return Monomial_type(d,a._variables[:],index,self), Monomial_type(r,a._variables[:],[0]*len(a._variables),self)
                        else:
                            D = QQ(self.subdomain())
                            MD = Poly(D)
                            v = d + QQ_type(D,r,b._value)
                            return Monomial_type(v,a._variables[:],index,MD), Monomial_type(0>>D,a._variables[:],[0]*len(a._variables),MD)
                    else:
                        return Monomial_type(0>>self.subdomain(),a._variables[:],[0]*len(a._variables),self), a.cp()
                else:
                    raise NotImplementedError
            elif isinstance(b,Poly_type):
                if a._varOrder == None:
                    q = Monomial_type(0>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                    p = a.cp()
                    r = q.cp()
                    ltb = b.leading_term()
                    fsalt = False
                    while(p != 0):
                        ltp = p.leading_term()
                        [d,ri] = self._div_res_core_(ltp,ltb)
                        if ri == 0:
                            if d.domain() != self:
                                q = q + d
                                p = p - d*b
                                fsalt = True
                            elif fsalt:
                                q = q + d
                                p = p - d*b
                            else:
                                q = self._add_core_(q,d)
                                p = self._sub_core_(p,self._mul_core_(d,b))
                        elif fsalt:
                            r = r + ltp
                            p = p - ltp
                        else:
                            r = self._add_core_(r,ltp)
                            p = self._sub_core_(p,ltp)
                    return q,r
                else:
                    raise NotImplementedError
        elif isinstance(a,Poly_type):
            if isinstance(b,Monomial_type):
                if a._varOrder == None:
                    q = Monomial_type(0>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                    p = a.cp()
                    r = q.cp()
                    ltb = b.leading_term()
                    fsalt = False
                    while(p != 0):
                        ltp = p.leading_term()
                        [d,ri] = self._div_res_core_(ltp,ltb)
                        if ri == 0:
                            if d.domain() != self:
                                q = q + d
                                p = p - d*b
                                fsalt = True
                            elif fsalt:
                                q = q + d
                                p = p - d*b
                            else:
                                q = self._add_core_(q,d)
                                p = self._sub_core_(p,ltp)
                        elif fsalt:
                            r = r + ltp
                            p = p - ltp
                        else:
                            r = self._add_core_(r,ltp)
                            p = self._sub_core_(p,ltp)
                    return q,r
                else:
                    raise NotImplementedError
            elif isinstance(b,Poly_type):
                if a._varOrder == None:
                    q = Monomial_type(0>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                    p = a.cp()
                    r = q.cp()
                    ltb = b.leading_term()
                    fsalt = False
                    while(p != 0):
                        ltp = p.leading_term()
                        [d,ri] = self._div_res_core_(ltp,ltb)
                        if ri == 0:
                            if d.domain() != self:
                                q = q + d
                                p = p - d*b
                                fsalt = True
                            elif fsalt:
                                q = q + d
                                p = p - d*b
                            else:
                                q = self._add_core_(q,d)
                                p = self._sub_core_(p,self._mul_core_(d,b))
                        elif fsalt:
                            r = r + ltp
                            p = p - ltp
                        else:
                            r = self._add_core_(r,ltp)
                            p = self._sub_core_(p,ltp)
                    return q,r
                else:
                    raise NotImplementedError

    def _div_(self, a, b):
        r0,r1 = self._transform_to_same_variables_varOrder_(a,b)
        ##TO DO: Pensar què fer en el cas varOrder != None
        r0 = r0._erase_var_order()
        r1 = r1._erase_var_order()
        ##
        return r0.domain()._div_core_(r0,r1)
    
    def _div_core_(self,a, b):
        if b == 0:
            raise TypeError("Polinomial mod must be different of zero")
        [dv, rs] = self._div_res_core_(a,b)
        if rs == 0:
            return dv
        else:
            D = QQ(dv.domain())
            return QQ_type(D,a,b)
        
    def _shift_pol_(self,pol,d):
        if d == 0:
            return pol.cp()
        elif d>0:
            if len(pol._variables) == 0:
                raise TypeError('Impossible to shift a costant without variable')
            if pol.nvariables() > 1:
                if isinstance(pol,Monomial_type):
                    if pol._varOrder == None:
                        aindex = pol._index[:]
                        aindex[-1] += d
                        return Monomial_type(pol._value, pol._variables[:],aindex, self)
                    else:
                        aindex = pol._index[:]
                        aindex[pol._varOrder[-1]] += d
                        return Monomial_type(pol._value, pol._variables[:],aindex, self, pol._varOrder[:])
                elif isinstance(pol,Poly_type):
                    if pol._varOrder == None:
                        aux = []
                        for i in pol._value:
                            aindex = pol._index[:]
                            aindex[-1] += d
                            aux.append(Monomial_type(pol._value, pol._variables[:],aindex, self))
                        return Poly_type(aux,pol._variables[:],self,pol._degree + d, pol._mdegree + d)
                    else:
                        aux = []
                        for i in pol._value:
                            aindex = pol._index[:]
                            aindex[pol._varOrder[-1]] += d
                            aux.append(Monomial_type(pol._value, pol._variables[:],aindex, self, pol._varOrder[:]))
                        return Poly_type(aux,pol._variables[:],self,pol._degree + d, pol._mdegree + d, pol._ord[:], pol._varOrder[:])
            if pol.nvariables() == 1:
                if isinstance(pol,Monomial_type):
                    return Monomial_type(pol._value, pol._variables[:],[pol._index[0] + d], self)
                elif isinstance(pol,Poly_type):
                    aux = []
                    for i in pol._value:
                        aux.append(Monomial_type(i._value, i._variables[:],[i._index[0] + d], self))
                    return Poly_type(aux,pol._variables[:],self,pol._degree + d, pol._mdegree + d)
        else:
            if pol.trailing_term().degree() < d:
                raise TypeError('d must be at least the trailing degree')
            if len(pol._variables) == 0:
                raise TypeError('Impossible to shift a costant without variable')
            if pol.nvariables() > 1:
                if isinstance(pol,Monomial_type):
                    if pol._varOrder == None:
                        aindex = pol._index[:]
                        aindex[-1] -= d
                        return Monomial_type(pol._value, pol._variables[:],aindex, self)
                    else:
                        aindex = pol._index[:]
                        aindex[pol._varOrder[-1]] -= d
                        return Monomial_type(pol._value, pol._variables[:],aindex, self, pol._varOrder[:])
                elif isinstance(pol,Poly_type):
                    if pol._varOrder == None:
                        aux = []
                        for i in pol._value:
                            aindex = pol._index[:]
                            aindex[-1] -= d
                            aux.append(Monomial_type(pol._value, pol._variables[:],aindex, self))
                        return Poly_type(aux,pol._variables[:],self,pol._degree - d, pol._mdegree - d)
                    else:
                        aux = []
                        for i in pol._value:
                            aindex = pol._index[:]
                            aindex[pol._varOrder[-1]] -= d
                            aux.append(Monomial_type(pol._value, pol._variables[:],aindex, self, pol._varOrder[:]))
                        return Poly_type(aux,pol._variables[:],self,pol._degree - d, pol._mdegree - d, pol._ord[:], pol._varOrder[:])
            if pol.nvariables() == 1:
                if isinstance(pol,Monomial_type):
                    return Monomial_type(pol._value, pol._variables[:],[pol._index[0] - d], self)
                elif isinstance(pol,Poly_type):
                    aux = []
                    for i in pol._value:
                        aux.append(Monomial_type(i._value, i._variables[:],[i._index[0] - d], self))
                    return Poly_type(aux,pol._variables[:],self,pol._degree - d, pol._mdegree - d)
                
    def __hash__(self):
        return PGl.HPOL + hash(self.subdomain())
    
    def __eq__(self, other):
        if isinstance(other, Base.Base):
            if other._poly_structure:
                if other._subDomain == self._subDomain:
                    return True
        return False
        
    def __neq__(self, other):
        return not self.__eq__(other)
      
    def _eq_(self, a, b):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return self._eq_(aux,b)
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return self._eq_(a,aux)
        
        if self.is_element(a):
            if not self.is_element(b):
                if isinstance(a, Poly_type):
                    return False
                if len(a._variables) == 0:
                    return a._value == b
                else:
                    if a._totIndex > 0:
                        return False
                    else:
                        return a._value == b
            else:
                if isinstance(a, Poly_type):
                    if isinstance(b, Poly_type):
                        if len(a._value) != len(b._value):
                            return False
                        for i in range(len(a._value)):
                            if a._value[i]!=b._value[i]:
                                return False
                        return True
                    else:
                        return False
                elif isinstance(a, Monomial_type):
                    if isinstance(b,Poly_type):
                        return False
                    elif isinstance(b, Monomial_type):
                        if a._value != b._value:
                            return False
                        ia = 0
                        ib = 0
                        while ((ia < len(a._index))and(ib < len(b._index))):
                            if a._index[ia] == 0:
                                ia += 1
                                continue
                            if b._index[ib] == 0:
                                ib += 1
                                continue
                            if a._variables[ia] != b._variables[ib]:
                                return False
                            if a._index[ia] != b._index[ib]:
                                return False
                            ia += 1
                            ib += 1
                        if ia == len(a._index):
                            if ib == len(b._index):
                                return True
                            else:
                                for i in range(ib, len(b._index)):
                                    if b._index[i] > 0:
                                        return False
                                return True
                        else:
                            for i in range(ia, len(a._index)):
                                if a._index[i] > 0:
                                    return False
                            return True
                    elif isinstance(b,Symbol_type):
                        return self._eq_(a, Monomial_type(1>>self.subdomain(),[b],[1],self))
                elif isinstance(a,Symbol_type):
                    if isinstance(b,Poly_type):
                        return False
                    elif isinstance(b, Monomial_type):
                        return self._eq_(Monomial_type(1>>self.subdomain(),[a],[1],self), b)
                    elif isinstance(b, Symbol_type):
                        return a.__eq__(b)
            
    def __rrshift__(self, value, symb = None):
        if isinstance(value,float):
            aux = value>>QQ(self.base_domain())
            D = self.union_domain(aux.domain())[0]
            return D.__rrshift__(aux)
        if isinstance(value, Symbol_type):
            return Monomial_type(1>>self.subdomain(),[value],[1],self)
        if (self.is_sub_element(value))and(not self.is_element(value)):
            if symb == None:
                return Monomial_type(value >> self.subdomain(), [],[], self)
            elif isinstance(symb, list):
                avariables = []
                for s in symb:
                    avariables.append(_get_symbol_(s))
                return Monomial_type(value >> self.subdomain(), avariables,[0]*len(avariables), self)
        if (self.is_element(value)):
            if symb == None:
                return value.cp()
            else:
                aux = value.cp()
                for s in symb:
                    x = _get_symbol_(s)
                    if x in aux._variables:
                        continue
                    else:
                       aux = aux.add_variable(x)
                return aux
        if value.domain()._poly_structure:
            if value.domain().subdomain().is_sub_domain(self.subdomain()):
                a = value.transform_to_min_domain()
                if a.domain() == self:
                    return a
                elif self.subdomain().is_sub_domain(a.K_()):
                    if isinstance(a,Monomial_type):
                        s = a._value >> self.subdomain()
                        if s != 0:
                            return Monomial_type(s,a._variables[:],a._index[:],self)
                        else:
                            return Monomial_type(s,a._variables[:],[0]*len(a._variables),self)
                    else:
                        aux = []
                        md = 0
                        for i in a._value:
                            s = i._value >> self.subdomain()
                            if s != 0:
                                aux.append(Monomial_type(s,a._variables[:],i._index[:],self))
                                if aux[-1]._totIndex > md:
                                    md = aux[-1]._totIndex
                        if len(aux) == 0:
                            return Monomial_type(0>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                        elif len(aux) == 1:
                            return aux[0]
                        else:
                            return Poly_type(aux,a._variables[:],self,aux[0]._degree,md)
                elif a.K_().is_sub_domain(self.subdomain()):
                    return a
                else:
                    raise TypeError("This element cannot be converted")
            elif self.subdomain().is_sub_domain(value.K_()):
                if isinstance(value,Monomial_type):
                    s = value._value >> self.subdomain()
                    if s != 0:
                        return Monomial_type(s,value._variables[:],value._index[:],self)
                    else:
                        return Monomial_type(s,value._variables[:],[0]*len(value._variables),self)
                else:
                    aux = []
                    md = 0
                    for i in value._value:
                        s = i._value >> self.subdomain()
                        if s != 0:
                            aux.append(Monomial_type(s,value._variables[:],i._index[:],self))
                            if aux[-1]._totIndex > md:
                                md = aux[-1]._totIndex
                    if len(aux) == 0:
                        return Monomial_type(0>>self.subdomain(),value._variables[:],[0]*len(value._variables),self)
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        return Poly_type(aux,value._variables[:],self,aux[0]._degree,md)
            else:
                raise TypeError("This element cannot be converted")
        else:
            a = value.transform_to_min_domain()
            if self.subdomain().is_sub_domain(a.K_()):
                return Monomial_type(a>>self._subdomain(),[],[],self)
            elif a.subdomain().is_sub_domain(self.subdomain()):
                return Monomial_type(a,[],[],Poly(a.domain()))
            else:
                raise TypeError("This element cannot be converted")
                    
    
    def _create_Poly_Structure_(self,D):
        return Poly(D)
    
    def variable_d(self,n, symb = 'x'):
        x = _get_symbol_(symb)
        return Monomial_type(1>>self.subdomain(),[x],[n],self)
    
    def _leading_term_(self,f):
        return f.leading_term()
    
    def _leading_coef_(self,f):
        return f.leading_coeff()
        
    def _trunc_inv_2_(self,f,k):
        if (self.is_sub_element(f)) and (not self.is_element(f)):
            return Monomial_type((1/f) >> self.subdomain(), [], [], self)
        if len(f._variables) == 0:
            return (1/f)
        if len(f._variables) == 1:
            if isinstance(f, Monomial_type):
                if f._totIndex > 0:
                    raise TypeError('Polynomial must have independent term')
                else:
                    return Monomial_type(1/f._value, f._variables[:], f._index[:], self)
            elif isinstance(f,Symbol_type):
                raise TypeError('Polynomial must have independent term')
            else:
                x = f._variables[0]
                ind = 0
                if f._varOrder == None:
                    f._varOrder = None
                else:
                    varOrder = f._varOrder[:]
                if f._value[-1]._index[0] != 0:
                    raise TypeError('Polynomial must have independent term')
                g = (1/f._value[-1]._value) >> self
        else:
            if f._varOrder == None:
                x = f._variables[-1]
                ind = -1
                varOrder = None
            else:
                x = f._variables[f._varOrder[-1]]
                ind = f._varOrder[-1]
                varOrder = f._varOrder[:]
            if f._value[-1]._totIndex != 0:
                    raise TypeError('Polynomial must have independent term')
            g = (1/f._value[-1]._value) >> self

        r = math.ceil(math.log2(k))
        t = 2**r
        dos = Monomial_type(2>>self.subdomain(),f._variables[:],[0]*len(f._variables[:]),self,varOrder)
        for j in range(r):
            t = t/2
            e = Monomial_type(1>>self.subdomain(), f._variables[:], [0]*(ind) + [math.ceil(k/t)] + [0]*(len(f._variables) - ind - 1), self, varOrder)
            g2 = self._md_square_core_(g,e)
            g = self._md_mul_core_(dos,g,e)
            g2 = self._md_mul_core_(f,g2,e)
            g = self._sub_core_(g,g2) 
        return g   
        
    def _trunc_inv_(self,f,k):
        if (self.is_sub_element(f)) and (not self.is_element(f)):
            return Monomial_type((1/f) >> self.subdomain(), [], [], self)
        if len(f._variables) == 0:
            return (1/f)
        if len(f._variables) == 1:
            if isinstance(f, Monomial_type):
                if f._totIndex > 0:
                    raise TypeError('Polynomial must have independent term')
                else:
                    return Monomial_type(1/f._value, f._variables[:], f._index[:], self)
            elif isinstance(f,Symbol_type):
                raise TypeError('Polynomial must have independent term')
            else:
                x = f._variables[0]
                ind = 0
                if f._varOrder == None:
                    f._varOrder = None
                else:
                    varOrder = f._varOrder[:]
                if f._value[-1]._index[0] != 0:
                    raise TypeError('Polynomial must have independent term')
                g = (1/f._value[-1]._value) >> self
        else:
            if f._varOrder == None:
                x = f._variables[-1]
                ind = -1
                varOrder = None
            else:
                x = f._variables[f._varOrder[-1]]
                ind = f._varOrder[-1]
                varOrder = f._varOrder[:]
            if f._value[-1]._totIndex != 0:
                    raise TypeError('Polynomial must have independent term')
            g = (1/f._value[-1]._value) >> self
        
        r = math.ceil(math.log2(k))
        t = 1
        dos = Monomial_type(2>>self.subdomain(),f._variables[:],[0]*len(f._variables[:]),self,varOrder)
        for j in range(r):
            t = 2*t
            e = Monomial_type(1>>self.subdomain(), f._variables[:], [0]*(ind) + [t] + [0]*(len(f._variables) - ind - 1), self, varOrder)
            g2 = self._md_mul_core_(g,f,e)
            g = self._md_mul_core_(dos,g,e)
            g2 = self._md_mul_core_(g,g2,e)
            g = self._sub_core_(g,g2) 
        return g
    
    def _fast_quo_rem_(self,f0,f1):
        r0,r1 = self._transform_to_same_variables_varOrder_(f0,f1)
        ##TO DO: Pensar què fer en el cas varOrder != None
        r0 = r0._erase_var_order()
        r1 = r1._erase_var_order()
        ##
        return r0.domain()._fast_quo_rem_core_(r0,r1)
        
        
    def _fast_quo_rem_core_(self, a, b):
        if isinstance(a,Monomial_type):
            return self._div_res_core_(a,b)
        if isinstance(b, Monomial_type):
            return self._div_res_core_(a,b)
        
        if a._varOrder == None:
            m = a.degree() - b.degree()
            h = b.reverse()
            q = a.reverse()*h.truncated_inverse(m+1)
            
            q = q.trunc(m+1)
            q = q.reverse_m(m)
            return [q, self._sub_core_(a, self._mul_core_(q,b))]  
        else:
            raise NotImplementedError
          
    def _rd(self,r, symb = 'x'): # random polynomial of degree r
        x = _get_symbol_(symb)
        A = self.subdomain()
        M = []
        for k in range(r,-1,-1):
            a = A._rd()
            if a != 0:
                M.append(Monomial_type(a,[x],[k],self))
        if len(M) == 0:
            return Monomial_type(0>>self.subdomain(),[x],[0],self)
        elif len(M) == 1:
            return M[0]
        else:
            return Poly_type(M,[x],self,M[0]._degree, M[0]._degree)
    
    def _rd_nonzero(self,r, symb = 'x'): # random polynomial of degree r
        x = _get_symbol_(symb)
        A = self.subdomain()
        q = A.cardinal()
        M = []
        l = rd_int(0,r)
        for k in range(r,-1,-1):
            if k == l:
                a = A._rd_nonzero()
            else:
                a = A._rd()
            
            if a != 0:
                M.append(Monomial_type(a,[x],[k],self))
        if len(M) == 1:
            return M[0]
        else:
            return Poly_type(M,[x],self,M[0]._degree, M[0]._degree)
         
    def change_characteristic(self,n):    
        if self.characteristic() == n:
            return self.cp()
        if (n == 0):
            return self.transform_to_ZZ()
        cself = self.characteristic()
        if ((cself%n) == 0):
            SDom = self.subdomain().change_characteristic(n)
            return Poly(SDom)
        raise TypeError("Wrong characteristic")        
    
    def transform_to_ZZ(self):
        if self.characteristic() == 0:
            return self
        SDom = self.subdomain().transform_to_ZZ()
        return Poly(SDom)
    
        
    def K_(self):
        return self.subdomain().K_()
        
            
    def _get_subpolynomial_degree_(self,f,d, var = None):
        return f.get_subpolynomial_degree(d,var)
    
    def _eval_(self,pol,a,variables=None):
        if isinstance(pol,Symbol_type):
            return pol.eval_symbol(a,variables)
        elif isinstance(pol,Monomial_type):
            return pol.eval_monomial(a,variables)
        elif isinstance(pol,Poly_type):
            return pol.eval_poly(a,variables)
        else:
            return pol>>self
    
    def _opposite_(self, a):
        if isinstance(a,Symbol_type):
            return a.opposite()
        elif isinstance(a,Monomial_type):
            return a.opposite()
        elif isinstance(a, Poly_type):
            M = []
            if a._varOrder == None:
                varOrder = None
                aorder = None
            else:
                varOrder = a._varOrder[:]
                aorder = a._ord[:]
            for m in a._value:
                M.append(Monomial_type(-m._value,m._variables[:],m._index[:],self, varOrder))
            return Poly_type(M,a._variables[:],self,a._degree,a._mdegree,aorder,varOrder)
        else:
            return -a
    
    def _der_(self,a, symb = None):
        if isinstance(a,Symbol_type):
            return a.der(symb)
        elif isinstance(a,Monomial_type):
            return a.der(symb)
        elif isinstance(a,Poly_type):
            return a.der(symb)
        else:
            if symb == None:
                return Monomial_type(0>>self.subdomain(),[],[],self)
            else:
                x = _get_symbol_(symb)
                return Monomial_type(0>>self.subdomain(),[x],[0],self)
    
#    def _reduce_aux_(self,coef,degf,degg,g,gt):
#        s = floor(degf/degg)
#        r = degf % degg
#        if len(g._variables) == 1:
#            f = Poly_type([coef],self,g._variables,[[r]], g._pri_var, goodValues=1)
#        else:
#            f = Poly_type([coef],self,g._variables,[[0]*(g._pri_var-1) + [r] + [0]*(len(g._variables) - g._pri_var)],g._pri_var,goodValues=1)
#        aux = g**s#self._md_pow_core_(g,s,gt)
#        f = aux*f#self._md_mul_core_(aux,f,gt)
#        return f

    def _md_reduce_core_(self,f,g):
        if isinstance(f,Monomial_type):
            if isinstance(g,Monomial_type):
                if len(f._variables) == 0:
                    raise TypeError('The degree of the modul must be at least 1')
                elif len(f._variables) == 1:
                    if f._degree < g._degree:
                        return f.cp()
                    else:
                        return Monomial_type(0>>self.subdomain(),f._variables[:],[0],self)
                else:
                    try:
                        x = _get_symbol_(g)
                    except:
                        for i in range(len(f._index)):
                            if f._index[i] < g._index[i]:
                                return f.cp()
                        return Monomial_type(0>>self.subdomain(),f._variables[:],[0]*len(f._variables),self)
                    else:
                        if f.degree(x) < g.degree(x):
                            return f.cp()
                        else:
                            return Monomial_type(0>>self.subdomain(),f._variables[:],[0],self)
            elif (g,Poly_type):
                if len(f._variables) == 0:
                    raise TypeError('The degree of the modul must be at least 1')
                elif len(f._variables) == 1:
                    if f._degree < g._degree:
                        return f.cp()
                    else:
                        s = floor(f._degree/g._degree)
                        r = f._degree % g._degree
                        if s < 2:
                            if len(g._value) == 2:
                                if isinstance(f._value,int):
                                    cf =-g._value[1]._value*f._value / g._value[0]._value >> QQ(ZZ()) 
                                    if cf._value[1] == 1:
                                        ff = Monomial_type(cf._value[0], f._variables[:],[g._value[1]._index[0] + r],self)
                                    else:
                                        raise TypeError("leading coefficient is not invertible")
                                else:
                                    cf =-g._value[1]._value*f._value / g._value[0]._value 
                                    ff = Monomial_type(cf, f._variables[:],[g._value[1]._index[0] + r],self)
                            else:
                                aux = []
                                for i in range(1,len(g._value)):
                                    aux.append(Monomial_type(-g._value[i]._value*f._value/g._value[0]._value, f._variables[:],[g._value[i]._index[0] + r],self))
                                if g._value[0]._totIndex < g._mdegree:
                                    md = g._degree
                                else:
                                    md = 0
                                    for a in aux:
                                        if a._totIndex == g._degree:
                                            md = g._degree
                                            break
                                        elif a._totIndex > md:
                                            md = a._totIndex
                                ff = Poly_type(aux,f._variables[:], self, aux[0]._degree, md)
                            return self._md_reduce_core_(ff,g)
                        else:
                            aux = []
                            for i in range(1,len(g._value)):
                                aux.append(Monomial_type(-g._value[i]._value/g._value[0]._value, f._variables[:],g._index[:],self))
                            if len(aux) == 1:
                                resmod = aux[0]
                            else:
                                if g._value[0]._totIndex < g._mdegree:
                                    md = g._degree
                                else:
                                    md = 0
                                    for a in aux:
                                        if a._totIndex == g._degree:
                                            md = g._degree
                                            break
                                        elif a._totIndex > md:
                                            md = a._totIndex
                                resmod = Poly_type(aux,f._variables[:],self,aux[0]._degree,md)
                            ff = self._md_pow_(resmod,s,g)
                            return self._md_reduce_core_(self._mul_core_(ff,Monomial_type(f._value,f._variables[:],[r],self)),g)
                else:
                    try:
                        x = _get_symbol_(g)
                    except:
                        mind = f._mdegree
                        for i in range(len(f._index)):
                            if f._index[i] < g._value[0]._index[i]:
                                return f.cp()
                            a = f._index[i] // g._value[0]._index[i]
                            if a < mind:
                                mind = a
                        
                        aux = []
                        for i in range(1,len(g._value)):
                            aux.append(Monomial_type(-g._value[i]._value/g._value[0]._value, f._variables[:],g._index[:],self))
                        if len(aux) == 1:
                            resmod = aux[0]
                        else:
                            if g._value[0]._totIndex < g._mdegree:
                                md = g._degree
                            else:
                                md = 0
                                for a in aux:
                                    if a._totIndex == g._degree:
                                        md = g._degree
                                        break
                                    elif a._totIndex > md:
                                        md = a._totIndex
                            resmod = Poly_type(aux,f._variables[:],self,aux[0]._degree,md)
                        ff = self._md_pow_(resmod,mind,g)
                        aindex = []
                        for i in range(len(f._index)):
                            aindex.append(f._index[i] - g._index[i]*mind)
                        return self._md_reduce_core_(self._mul_core_(ff,Monomial_type(f._value,f._variables[:],aindex,self)),g)                                
                    else:
                        degf = f.degree(x)
                        degg = g.degree(x)
                        indx = f._variables.index(x)
                        if f._varOrder == None:
                            varOrder = None
                        else:
                            varOrder = f._varOrder[:]
                        if degf < degg:
                            return f.cp()
                        else:
                            s = floor(degf/degg)
                            r = degf % degg
                            if s < 2:
                                if len(g._value) == 2:
                                    ff = Monomial_type(-g._value[1]._value*f._value / g._value[0]._value, f._variables[:],g._value[1]._index[:indx] + [g._value[1]._index[indx] + r] + g._value[1]._index[indx+1:],self,varOrder)
                                else:
                                    aux = []
                                    for i in range(1,len(g._value)):
                                        aux.append(Monomial_type(-g._value[i]._value*f._value/g._value[0], f._variables[:],g._value[1]._index[:indx] + [g._value[1]._index[indx] + r] + g._value[1]._index[indx+1:],self,varOrder))
                                    if g._value[0]._totIndex < g._mdegree:
                                        md = g._degree
                                    else:
                                        md = 0
                                        for a in aux:
                                            if a._totIndex == g._degree:
                                                md = g._degree
                                                break
                                            elif a._totIndex > md:
                                                md = a._totIndex
                                        aord = []
                                        for i in range(len(g._ord)):
                                            if g._ord[i] > 0:
                                                aord.append(g._ord[i] - 1)
                                    ff = Poly_type(aux,f._variables[:], self, aux[0]._degree, md,aord,varOrder)
                                return self._md_reduce_core_(ff,g)
                            else:
                                aux = []
                                for i in range(1,len(g._value)):
                                    aux.append(Monomial_type(-g._value[i]._value/g._value[0]._value, f._variables[:],g._index[:],self,varOrder))
                                if len(aux) == 1:
                                    resmod = aux[0]
                                else:
                                    if g._value[0]._totIndex < g._mdegree:
                                        md = g._degree
                                    else:
                                        md = 0
                                        for a in aux:
                                            if a._totIndex == g._degree:
                                                md = g._degree
                                                break
                                            elif a._totIndex > md:
                                                md = a._totIndex
                                        aord = []
                                        for i in range(len(g._ord)):
                                            if g._ord[i] > 0:
                                                aord.append(g._ord[i] - 1)
                                    resmod = Poly_type(aux,f._variables[:],self,aux[0]._degree,md,aord,varOrder)
                                ff = self._md_pow_(resmod,s,g)
                                fnx = Monomial_type(f._value,f._variables[:],f._index[:indx]  + [r] + f._index[indx+1:])
                                return self._md_reduce_core_(self._mul_core_(ff,fnx),g)
        elif isinstance(f,Poly_type):
            if isinstance(g,Monomial_type):
                if len(f._variables) == 0:
                    raise TypeError('The degree of the modul must be at least 1')
                elif len(f._variables) == 1:
                    if f._degree < g._degree:
                        return f.cp()
                    else:
                        return Monomial_type(0>>self.subdomain(),f._variables[:],[0],self)
                else:
                    try:
                        x = _get_symbol_(g)
                    except:
                        for i in range(len(f._index)):
                            if f._index[i] < g._index[i]:
                                return f.cp()
                        return Monomial_type(0>>self.subdomain(),f._variables[:],[0]*len(f._variables),self)
                    else:
                        if f.degree(x) < g.degree(x):
                            return f.cp()
                        else:
                            return Monomial_type(0>>self.subdomain(),f._variables[:],[0],self)
            elif isinstance(g,Poly_type):
                if len(f._variables) == 0:
                    raise TypeError('The degree of the modul must be at least 1')
                elif len(f._variables) == 1:
                    if f._degree < g._degree:
                        return f.cp()
                    else:
                        if len(f._value) == 2:
                            fnl = f._value[1].cp()
                        else:
                            if f._varOrder == None:
                                fnl = Poly_type(f._value[1:],f._variables[:],self,f._value[1]._degree,f._value[1]._degree)
                            else:
                                aord = []
                                for i in range(len(f._ord)):
                                    if f._ord[i]>0:
                                        aord.append(f._ord[i] - 1)
                                fnl = Poly_type(f._value[1:],f._variables[:],self,f._value[1]._degree,f._value[1]._degree,aord,f._varOrder[:])
                        return self._add_core_(self._md_reduce_core_(fnl,g),self._md_reduce_core_(f._value[0],g))
                else:
                    try:
                        x = _get_symbol_(g)
                    except:
                        if len(f._value) == 2:
                            fnl = f._value[1].cp()
                        else:
                            if f._varOrder == None:
                                if f._value[0]._totIndex < f._mdegree:
                                    md = f._mdegree
                                else:
                                    md = 0
                                    for i in range(1,len(f._value)):
                                        if f._value[i]._totIndex == f._mdegree:
                                            md = f._mdegree
                                            break
                                        elif f._value[i]._totIndex > md:
                                            md = f._value[i]._totIndex
                                fnl = Poly_type(f._value[1:],f._variables[:],self,f._value[1]._degree,md)
                            else:
                                aord = []
                                for i in range(len(f._ord)):
                                    if f._ord[i]>0:
                                        aord.append(f._ord[i] - 1)
                                fnl = Poly_type(f._value[1:],f._variables[:],self,f._value[1]._degree,md,aord,f._varOrder[:])
                            return self._add_core_(self._md_reduce_core_(fnl,g),self._md_reduce_core_(f._value[0],g))
                    else:
                        if f.degree(x) < g.degree(x):
                            return f.cp()
                        else:
                            if len(f._value) == 2:
                                fnl = f._value[1].cp()
                            else:
                                if f._varOrder == None:
                                    if f._value[0]._totIndex < f._mdegree:
                                        md = f._mdegree
                                    else:
                                        md = 0
                                        for i in range(1,len(f._value)):
                                            if f._value[i]._totIndex == f._mdegree:
                                                md = f._mdegree
                                                break
                                            elif f._value[i]._totIndex > md:
                                                md = f._value[i]._totIndex
                                    fnl = Poly_type(f._value[:1],f._variables[:],self,f._value[1]._degree,md)
                                else:
                                    aord = []
                                    for i in range(len(f._ord)):
                                        if f._ord[i]>0:
                                            aord.append(f._ord[i] - 1)
                                    fnl = Poly_type(f._value[:1],f._variables[:],self,f._value[1]._degree,md,aord,f._varOrder[:])
                                return self._add_core_(self._md_reduce_core_(fnl,g),self._md_reduce_core_(f._value[0],g))
    
    def _md_reduce_(self,f,g):
        '''return f mod g'''
        r0,r1 = self._transform_to_same_variables_varOrder_(f,g)
        return r0.domain()._md_reduce_core_(r0,r1)     
            
    def _md_mul_core_(self,f1,f2,g):
        if isinstance(g,Monomial_type):
            if len(g._variables) == 0:
                raise TypeError('The degree of the modul must be at least 1')
            elif len(g._variables) == 1:
                if isinstance(f1,Monomial_type):
                    if isinstance(f2,Monomial_type):
                        if f1._degree + f2._degree >= g._degree:
                            return Monomial_type(0>>self.subdomain(),g._variables[:],[0],self) 
                        if f1._varOrder == None:
                            s = f1._value*f2._value
                            if s == 0:
                                return Monomial_type(s,g._variables[:], [0]*len(g._variables), self)
                            else:
                                return Monomial_type(s,g._variables[:], [f1._index[0] + f2._index[0]], self)
                        else:
                            s = f1._value*f2._value
                            if s == 0:
                                return Monomial_type(s,g._variables[:], [0]*len(g._variables), self)
                            else:
                                return Monomial_type(s,g._variables[:], [f1._index[0] + f2._index[0]], self, g._varOrder[:])
                    elif isinstance(f2,Poly_type):
                        if g._varOrder == None:
                            aux = []
                            md = 0
                            for i in range(len(f2._value)):
                                if f1._degree + f2._value[i]._degree >= g._degree:
                                    continue
                                s = f1._value*f2._value[i]._value
                                if s != 0:
                                    index = [f1._degree + f2._value[i]._degree]
                                    aux.append(Monomial_type(s,g._variables[:], index, self))
                                    if aux[-1]._totIndex > md:
                                        md = aux[-1]._totIndex
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, md)
                        else:
                            aux = []
                            md = 0
                            for i in range(len(f2._value)):
                                if f1._degree + f2._value[i]._degree >= g._degree:
                                    continue
                                s = f1._value*f2._value[i]._value
                                if s != 0:
                                    index = [f1._degree + f2._value[i]._degree]
                                    aux.append(Monomial_type(s,g._variables[:], index, self,g._varOrder[:]))
                                    if aux[-1]._totIndex > md:
                                        md = aux[-1]._totIndex
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                vord = _get_own_order_monomials_(aux[:])
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, md, vord, g._varOrder[:])
                elif isinstance(f1,Poly_type):
                    if isinstance(f2,Monomial_type):
                        if g._varOrder == None:
                            aux = []
                            md = 0
                            for i in range(len(f1._value)):
                                if f1._value[i]._degree + f2._degree >= g._degree:
                                    continue
                                s = f1._value[i]._value*f2._value
                                if s != 0:
                                    index = [f1._value[i]._degree + f2._degree]
                                    aux.append(Monomial_type(s,g._variables[:], index, self))
                                    if aux[-1]._totIndex > md:
                                        md = aux[-1]._totIndex
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, md)
                        else:
                            aux = []
                            md = 0
                            for i in range(len(f1._value)):
                                if f1._value[i]._degree + f2._degree >= g._degree:
                                    continue
                                s = f1._value[i]._value*f2._value
                                if s != 0:
                                    index = [f1._value[i]._degree + f2._degree]
                                    aux.append(Monomial_type(s,g._variables[:], index, self,g._varOrder[:]))
                                    if aux[-1]._totIndex > md:
                                        md = aux[-1]._totIndex
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                vord = _get_own_order_monomials_(aux[:])
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, md, vord, g._varOrder[:])
                    elif isinstance(f2,Poly_type):
                        if g._varOrder == None:
                            aux = []
                            for i in range(len(f1._value)):
                                ind = 0
                                for j in range(len(f2._value)):
                                    if f1._value[i]._degree + f2._value[j]._degree >= g._degree:
                                        continue
                                    index = [f1._value[i]._degree + f2._value[j]._degree]
                                    pos,fo = _find_index_monomials(index,aux,ind)
                                    if fo:
                                        s = aux[pos]._value + f1._value[i]._value*f2._value[j]._value
                                        if s == 0:
                                            aux = aux[:pos] + aux[pos+1:]
                                        else:
                                            aux[pos] = Monomial_type(s,g._variables[:],aux[pos]._index,self)
                                    else:
                                        s = f1._value[i]._value*f2._value[j]._value
                                        if s != 0:
                                            aux = aux[:pos] + [Monomial_type(s,g._variables[:],index,self)] + aux[pos:]
                                    ind = pos
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                l = f1._mdegree + f2._mdegree
                                md = 0
                                for i in aux:
                                    if l == i._totIndex:
                                        md = l
                                        break
                                    elif md < i._totIndex:
                                        md = i._totIndex
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, md)
                        else:
                            aux = []
                            for i in range(len(f1._value)):
                                ind = 0
                                for j in range(len(f2._value)):
                                    if f1._value[i]._degree + f2._value[j]._degree >= g._degree:
                                        continue
                                    index = [f1._value[i]._degree + f2._value[j]._degree]
                                    pos,fo = _find_index_monomials(index,aux,ind)
                                    if fo:
                                        s = aux[pos]._value + f1._value[i]._value*f2._value[j]._value
                                        if s == 0:
                                            aux = aux[:pos] + aux[pos+1:]
                                        else:
                                            aux[pos] = Monomial_type(s,g._variables[:],aux[pos]._index,self,g._varOrder[:])
                                    else:
                                        s = f1._value[i]._value*f2._value[j]._value
                                        if s != 0:
                                            aux = aux[:pos] + [Monomial_type(s,g._variables[:],index,self,g._varOrder[:])] + aux[pos:]
                                    ind = pos
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                l = f1._mdegree + f2._mdegree
                                md = 0
                                for i in aux:
                                    if l == i._totIndex:
                                        md = l
                                        break
                                    elif md < i._totIndex:
                                        md = i._totIndex
                                vord = _get_own_order_monomials_(aux[:])
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, md, vord, g._varOrder[:])
            else:
                try:
                    x = _get_symbol_(g)
                except:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f1,f2),g)
                else:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f1,f2),g)
        elif isinstance(g,Poly_type):
            if len(g._variables) == 0:
                raise TypeError('The degree of the modul must be at least 1')
            elif len(g._variables) == 1:
                #TO DO: fer-ho bé
                return self._md_reduce_core_(self._mul_core_(f1,f2),g)
            else:
                try:
                    x = _get_symbol_(g)
                except:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f1,f2),g)
                else:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f1,f2),g)
        
    def _md_mul_(self,a,b,c):
        f1,f2 = self._transform_to_same_variables_varOrder_(a,b)
        f1,g = f1.domain()._transform_to_same_variables_varOrder_(f1,c)
        f2,g = f1.domain()._transform_to_same_variables_varOrder_(f2,g)
        return f1.domain()._md_mul_core_(f1,f2,g)
    
    def _md_square_core_(self,f,g):
        if isinstance(g,Monomial_type):
            if len(g._variables) == 0:
                raise TypeError('The degree of the modul must be at least 1')
            elif len(g._variables) == 1:
                if isinstance(f,Monomial_type):
                    if 2*f._degree >= g._degree:
                        return Monomial_type(0>>self.subdomain(),g._variables[:],[0],self) 
                    if g._varOrder == None:
                        s = f._value*f._value
                        if s == 0:
                            return Monomial_type(s,g._variables[:], [0]*len(g._variables), self)
                        else:
                            return self._md_reduce_core_(Monomial_type(s,g._variables[:], [f._index[0] + f._index[0]], self),g)
                    else:
                        s = f._value*f._value
                        if s == 0:
                            return Monomial_type(s,g._variables[:], [0]*len(g._variables), self)
                        else:
                            return Monomial_type(s,g._variables[:], [f._index[0] + f._index[0]], self, g._varOrder[:])
                elif isinstance(f,Poly_type):
                    if g._varOrder == None:
                        if self.characteristic() == 2:
                            aux = []
                            for i in range(len(f._value)):
                                if 2*f._value[i]._degree >= g._degree:
                                    continue
                                s = f._value[i]._value * f._value[i]._value
                                if s != 0:
                                    aux.append(Monomial_type(s,g._variables[:],[2*f._value[i]._degree],self))
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, aux[0]._degree)
                        else:
                            dos = 2 >> self.subdomain()
                            aux = []
                            for i in range(len(f._value)):
                                ind = 0
                                for j in range(i,len(f._value)):
                                    if f._value[i]._degree + f._value[j]._degree >= g._degree:
                                        continue
                                    index = [f._value[i]._degree + f._value[j]._degree]
                                    pos,fo = _find_index_monomials(index,aux,ind)
                                    if fo:
                                        if i == j:
                                            s = aux[pos] + f._value[i]._value*f._value[j]._value
                                        else:
                                            s = aux[pos] + dos*f._value[i]._value*f._value[j]._value
                                        if s == 0:
                                            aux = aux[:pos] + aux[pos:]
                                        else:
                                            aux[pos] = Monomial_type(s,g._variables[:],aux[pos]._index,self)
                                    else:
                                        if i == j:
                                            s = f._value[i]._value*f._value[j]._value
                                        else:
                                            s = dos*f._value[i]._value*f._value[j]._value
                                        if s != 0:
                                            aux = aux[:pos] + [Monomial_type(s,g._variables[:],index,self)]
                                    ind = pos
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, aux[0]._degree)
                    else:
                        if self.characteristic() == 2:
                            aux = []
                            for i in range(len(f._value)):
                                if 2*f._value[i]._degree >= g._degree:
                                    continue
                                s = f._value[i]._value * f._value[i]._value
                                if s != 0:
                                    aux.append(Monomial_type(s,g._variables[:],[2*f._value[i]._degree],self))
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, aux[0]._degree,f._ord[:],g._varOrder[:])
                        else:
                            dos = 2 >> self.subdomain()
                            aux = []
                            for i in range(len(f._value)):
                                ind = 0
                                for j in range(len(i,f._value)):
                                    if f._value[i]._degree + f._value[j]._degree >= g._degree:
                                        continue
                                    index = [f._value[i]._degree + f._value[j]._degree]
                                    pos,fo = _find_index_monomials(index,aux,ind)
                                    if fo:
                                        if i == j:
                                            s = aux[pos] + f._value[i]._value*f._value[j]._value
                                        else:
                                            s = aux[pos] + dos*f._value[i]._value*f._value[j]._value
                                        if s == 0:
                                            aux = aux[:pos] + aux[pos:]
                                        else:
                                            aux[pos] = Monomial_type(s,g._variables[:],aux[pos]._index,self,g._varOrder[:])
                                    else:
                                        if i == j:
                                            s = f._value[i]._value*f._value[j]._value
                                        else:
                                            s = dos*f._value[i]._value*f._value[j]._value
                                        if s != 0:
                                            aux = aux[:pos] + [Monomial_type(s,g._variables[:],index,self,g._varOrder[:])]
                                    ind = pos
                            if len(aux) == 0:
                                return Monomial_type(0>>self.subdomain(),g._variables[:], [0]*len(g._variables), self)
                            elif len(aux) == 1:
                                return aux[0]
                            else:
                                l = 2*f._mdegree
                                md = 0
                                for i in aux:
                                    if l == i._totIndex:
                                        md = l
                                        break
                                    elif md < i._totIndex:
                                        md = i._totIndex
                                vord = _get_own_order_monomials_(aux[:])
                                return Poly_type(aux,g._variables[:],self,aux[0]._degree, md, vord, g._varOrder[:])
            else:
                try:
                    x = _get_symbol_(g)
                except:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f,f),g)
                else:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f,f),g)
        elif isinstance(g,Poly_type):
            if len(g._variables) == 0:
                raise TypeError('The degree of the modul must be at least 1')
            elif len(g._variables) == 1:
                #TO DO: fer-ho bé
                return self._md_reduce_core_(self._mul_core_(f,f),g)
            else:
                try:
                    x = _get_symbol_(g)
                except:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f,f),g)
                else:
                    #TO DO: fer-ho bé
                    return self._md_reduce_core_(self._mul_core_(f,f),g)
            
    def _md_square_(self,a,c):
        f,g = self._transform_to_same_variables_varOrder_(a,c)
        return f1.domain()._md_square_core_(f,g)

    
    def _md_double_pow_(self,f,a,b,g):
        r0,r1 = self._transform_to_same_variables_varOrder_(f,g)
        return r0.domain()._md_double_pow_core_(r0,a,b,r1)
    
    def _md_double_pow_core_(self,f,a,b,g):
        if f.degree() >= g.degree():
            f = self._md_reduce_core_(f,g)
        for i in range(b):
            f = self._md_pow_core_(f,a,g)
        return f
        
    def _md_pow_core_(self,a,b,c):
        if isinstance(b,int):
            if b<0:
                return self._md_pow_core_(1/a,-b,c)
            if b == 0:
                if a!=0:
                    if c._varOrder == None:
                        return Monomial_type(1>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
                    else:
                        return Monomial_type(1>>self.subdomain(),a._variables[:],[0]*len(a._variables),self,c._varOrder[:])
                else:
                    raise TypeError("0**0 indetermination")
            aux2 = a.cp()
            if c._varOrder == None:
                aux1 = Monomial_type(1>>self.subdomain(),a._variables[:],[0]*len(a._variables),self)
            else:
                aux1 = Monomial_type(1>>self.subdomain(),a._variables[:],[0]*len(a._variables),self,c._varOrder[:])
            
            while b>1:
                if b%2==0:
                    b//=2
                    aux2 = self._md_square_core_(aux2,c)
                else:
                    aux1 = self._md_mul_core_(aux2,aux1,c)
                    aux2 = self._md_square_core_(aux2,c)
                    b = (b-1)//2
            return self._md_mul_core_(aux1,aux2,c)
        else:
            raise TypeError("This element cannot be operated")
        
    def _md_pow_(self,a,b,c):
        f,g = self._transform_to_same_variables_varOrder_(a,c)
        return g.domain()._md_pow_core_(f,b,g)
    
    def _get_own_order_monomials_(self,L):
        return _get_own_order_monomials_(L)
    
    def is_element(self, a):
        """returns if a is un element of THIS space"""
        if isinstance(a, Base_type.Base_type):
            if a._structure == self:
                return True
        if isinstance(a,Symbol_type):
            return True
        return False
            
            
            
            
            
            
            
            
            
            