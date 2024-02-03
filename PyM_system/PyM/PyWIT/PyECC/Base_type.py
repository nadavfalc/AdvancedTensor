# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 16:42:49 2016

@author: admin_local
"""
import copy
import time

##
class Base_type(object):
    """Represents an abstract Element. """
    
    _value = None
    _structure = None
    _is_inv = -1
    
    
    
    def __init__(self):
        raise NotImplementedError
    
    def __add__(self, other):
        if isinstance(other, float):
            return self._structure._add_(self, other)
        if isinstance(other, int):
            return self._structure._add_(self, other)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._add_(self, other)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._add_(self, other)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._add_(aux1,aux2)                    
        try:
            return other.__radd__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
    
    def __sub__(self, other):
        if isinstance(other, float):
            return self._structure._sub_(self, other)
        if isinstance(other, int):
            return self._structure._sub_(self, other)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._sub_(self, other)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._sub_(self, other)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._sub_(aux1,aux2)
        try:
            return other.__rsub__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
    
    def __mul__(self, other):
        if isinstance(other, float):
            return self._structure._mul_(self, other)
        if isinstance(other, int):
            return self._structure._mul_(self, other)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._mul_(self, other)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._mul_(self, other)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._mul_(aux1,aux2)
        try:
            return other.__rmul__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
        
    def __pow__(self,other):
        if isinstance(other, int):
            if other<0:
                return (1/self)**(-other)
            if other == 0:
                return 1>>self._structure
            aux1 = 1>>self._structure
            aux2 = self
            while other>1:
                if other%2==0:
                    other //= 2
                    aux2 = aux2*aux2
                else:
                    aux1 = aux2*aux1
                    aux2 = aux2*aux2
                    other = (other-1)//2
            return aux1 * aux2
        else:
            raise TypeError("This element cannot be operated")
        
    def __truediv__(self, other):
        if isinstance(other, float):
            return self._structure._div_(self, other)
        if isinstance(other, int):
            return self._structure._div_(self, other)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._div_(self, other)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._div_(self, other)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._div_(aux1,aux2)
        try:
            return other.__rtruediv__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
        
    def __floordiv__(self, other):
        return self.__truediv__(other)
        
    def __radd__(self, other):
        if isinstance(other, float):
            return self._structure._add_(other, self)
        if isinstance(other, int):
            return self._structure._add_(other, self)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._add_(other, self)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._add_(other, self)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._add_(aux2,aux1)
        try:
            return other.__add__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
    
    def __rsub__(self, other):
        if isinstance(other, float):
            return self._structure._sub_(other, self)
        if isinstance(other, int):
            return self._structure._sub_(other, self)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._sub_(other, self)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._sub_(other, self)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._sub_(aux2,aux1)
        try:
            return other.__sub__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
    
    def __rmul__(self,other):
        if isinstance(other, float):
            return self._structure._mul_(other, self)
        if isinstance(other, int):
            return self._structure._mul_(other, self)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._mul_(other, self)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._mul_(other, self)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._mul_(aux2,aux1)
        try:
            return other.__mul__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
        
    def __rpow__(self,other):
        raise TypeError("This element cannot be operated")
        
    def __rtruediv__(self, other):
        if isinstance(other, float):
            return self._structure._div_(other, self)
        if isinstance(other, int):
            return self._structure._div_(other, self)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._div_(other, self)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._div_(other, self)
            else:
                aux1 = self._push_to_union(other.domain())
                aux2 = other._push_to_union(self.domain())
                return aux1.domain()._div_(aux2,aux1)
        try:
            return other.__truediv__(self)
        except:
            raise TypeError("I do not know how to opperate these elements")
        
    def __rfloordiv__(self, other):
        self.__rtruediv__(other)
            
    def _characteristic_(self):
        """returns the characteristic of its domain"""
        return self.domain().characteristic()
        
    def _cardinal_(self):
        """Returns the cardinal of its domain. """
        raise NotImplementedError('cardinal()')
        
    def __eq__(self, other):
        if isinstance(other,float):
            return self._structure._eq_(self, other)
        if isinstance(other, int):
            return self._structure._eq_(self, other)
        if isinstance(other, Base_type):
            if (self._structure.is_sub_domain(other._structure)):
                return self._structure._eq_(self, other)
            if (other._structure.is_sub_domain(self._structure)):
                return other._structure._eq_(other, self)
            else:
                return False
        return False
           
    
    def __req__(self, other):
        return self.__eq__(other)
    
    def __ne__(self, other):
        return not self.__eq__(other)

    def __rne__(self, other):
        return self.__ne__(other)
    
    def __neg__(self):
        raise NotImplementedError
        
    def domain(self):
        """returns the domain of the element"""
        return self._structure

    def _sign_value_(self):
        raise NotImplementedError
    
    def _len_value_(self):
        raise NotImplementedError
    
    def min_domain(self):
        """returns the minim domain where the element is define"""
        return self.domain()
    
    def transform_to_min_domain(self):
        """Transforms the element to itself in the minim domain where it is define"""
        return self
    
    def change_characteristic(self,n):
        """changes the characteristic of the element"""
        raise NotImplementedError
    
    def transform_to_ZZ(self):
        """Tranforms the element to have ZZ as base domain"""
        raise NotImplementedError
    
    def set_name(self,name):
        """changes the name of the element"""
        self._name = name
    
    def get_name(self):
        """returns the name of the element"""
        return self._name
    
    def cp(self):
        raise NotImplementedError
    
    def _push_to_subunion(self,A,D,S,KS,KA):
        if D.is_sub_domain(self.domain()):
            return self>>D
        if S[0][0] == 1:
            a = self._push_to_subunion(A,D.subdomain(),S[1:],KS,KA)
            return a>>D
        else:
            if S[0][1] == 1:
                aux = []
                if isinstance(self._value,int):
                    if self._varOrder == None:
                        M = D.element_monomial(self._value>>KA, self._variables[:], self._index[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0>>KA,self._variables[:],[0]*len(self._variables))
                    else:
                        M = D.element_monomial(self._value>>KA, self._variables[:], self._index[:],self._varOrder[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0>>KA,self._variables[:],[0]*len(self._variables))
                elif isinstance(self._value,list):
                    md = 0
                    if self._varOrder == None:
                        for i in self._value:
                            if isinstance(i._value,int):
                                M = D.element_monomial(self._value>>KA, self._variables[:], self._index[:])
                                if M._value != 0:
                                    aux.append(M)
                                    if md < M._totIndex:
                                        md = M._totIndex
                            else:
                                M = D.element_monomial(self._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA), self._variables[:], self._index[:])
                                if M._value != 0:
                                    aux.append(M)
                                    if md < M._totIndex:
                                        md = M._totIndex
                        if len(aux) == 0:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
                        elif len(aux) == 1:
                            return aux[0]
                        else:
                            return D.element(aux,self._variables[:],D,aux[0]._degree, md)
                    else:
                        for i in self._value:
                            if isinstance(i._value,int):
                                M = D.element_monomial(self._value>>KA, self._variables[:], self._index[:],self._varOrder[:])
                                if M._value != 0:
                                    aux.append(M)
                                    if md < M._totIndex:
                                        md = M._totIndex
                            else:
                                M = D.element_monomial(self._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA), self._variables[:], self._index[:],self._varOrder[:])
                                if M._value != 0:
                                    aux.append(M)
                                    if md < M._totIndex:
                                        md = M._totIndex
                        if len(aux) == 0:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
                        elif len(aux) == 1:
                            return aux[0]
                        else:
                            vord = D._get_own_order_monomials_(aux)
                            return D.element(aux,self._variables[:],D,aux[0]._degree, md,vord,self._varOrder[:])
                else:
                    if self._varOrder == None:
                        M = D.element_monomial(self._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA), self._variables[:], self._index[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0>>D.subdomain(),self._variables[:],[0]*len(self._variables))
                    else:
                        M = D.element_monomial(self._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA), self._variables[:], self._index[:],self._varOrder[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0>>D.subdomain(),self._variables[:],[0]*len(self._variables))
            elif S[0][1] == 2:
                f = self._value
                if isinstance(self.domain().subdomain().element(0),int):
                    func = lambda x: f(x) >> KA
                    point = self._point >> KA
                else:
                    func = lambda x: f(x)._push_to_union(A)
                    point = self._point._push_to_union(A)
                return D.element(func,point,self._nPrint)
            elif S[0][1] == 3:
                if isinstance(self._value[0],list):
                    if isinstance(self._value[0][0],int):
                        aux = []
                        for i in range(len(self._value)):
                            aux2 = []
                            for j in self._value[i]:
                                aux2.append(j>>KA)
                            aux.append(aux2)
                    else:
                        aux = []
                        for i in range(len(self._value)):
                            aux2 = []
                            for j in self._value[i]:
                                aux2.append(j._push_to_subunion(A,D.subdomain(),S[1:],KS,KA))
                            aux.append(aux2)
                    return D.element(aux)
                else:
                    if isinstance(self._value[0],int):
                        aux = []
                        for i in self._value:
                            aux.append(i>>KA)
                    else:
                        aux = []
                        for i in self._value:
                            aux.append(i._push_to_subunion(A,D.subdomain(),S[1:],KS,KA))
                    if self._direction:
                        return D.Vector_type(aux)
                    else:
                        return D.Vector_type(aux,'c')
    
    def _push_to_union(self,A):
        if A == self.domain():
            return self
        if A.is_sub_domain(self.domain()):
            return self>>A
        if self.domain().is_sub_domain(A):
            if A.is_sub_domain(self.min_domain()):
                return (self.transform_to_min_domain) >> A
            else:
                return self>>self.domain()
        [D,S,KS,KA] = self.domain().union_domain(A)
        S.reverse()
        if S[0][0] == 1:
            a = self._push_to_subunion(A,D.subdomain(),S[1:],KS,KA)
            return a>>D
        else:
            if S[0][1] == 1:
                aux = []
                if isinstance(self._value,int):
                    if self._varOrder == None:
                        M = D.element_monomial(self._value, self._variables[:], self._index[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
                    else:
                        M = D.element_monomial(self._value, self._variables[:], self._index[:],self._varOrder[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
                elif isinstance(self._value,list):
                    aux =[]
                    md = 0
                    if self._varOrder == None:
                        for i in self._value:
                            if isinstance(i._value,int):
                                M = D.element_monomial(i._value,i._variables[:],i._index[:])
                            else:
                                M = D.element_monomial(i._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA),i._variables[:],i._index[:])
                            if M._value != 0:
                                aux.append(M)
                                if md < aux[-1]._totIndex:
                                    md = aux[-1]._totIndex
                        if len(aux) == 0:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
                        elif len(aux) == 1:
                            return aux[0]
                        else:
                            return D.element(aux,self._variables[:],D,aux[0]._degree, md)
                    else:
                        for i in self._value:
                            if isinstance(i._value,int):
                                M = D.element_monomial(i._value,i._variables[:],i._index[:], i._varOrder[:])
                            else:
                                M = D.element_monomial(i._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA),i._variables[:],i._index[:], i._varOrder[:])
                            if M._value != 0:
                                aux.append(M)
                                if md < aux[-1]._totIndex:
                                    md = aux[-1]._totIndex
                        if len(aux) == 0:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
                        elif len(aux) == 1:
                            return aux[0]
                        else:
                            vord = D._get_own_order_monomials_(aux)
                            return D.element(aux,self._variables[:],D,aux[0]._degree, md,vord,self._varOrder[:])
                else:
                    if self._varOrder == None:
                        M = D.element_monomial(self._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA), self._variables[:], self._index[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
                    else:
                        M = D.element_monomial(self._value._push_to_subunion(A,D.subdomain(),S[1:],KS,KA), self._variables[:], self._index[:],self._varOrder[:])
                        if M._value != 0:
                            return M
                        else:
                            return D.element_monomial(0,self._variables[:],[0]*len(self._variables))
            elif S[0][1] == 2:
                f = self._value
                if isinstance(self.domain().subdomain().element(0),int):
                    func = lambda x: f(x) >> KA
                    point = self._point >> KA
                else:
                    func = lambda x: f(x)._push_to_union(A)
                    point = self._point._push_to_union(A)
                return D.element(func,point,self._nPrint)
            elif S[0][1] == 3:
                if isinstance(self._value[0],list):
                    if isinstance(self._value[0][0],int):
                        aux = []
                        for i in range(len(self._value)):
                            aux2 = []
                            for j in self._value[i]:
                                aux2.append(j>>KA)
                            aux.append(aux2)
                    else:
                        aux = []
                        for i in range(len(self._value)):
                            aux2 = []
                            for j in self._value[i]:
                                aux2.append(j._push_to_subunion(A,D.subdomain(),S[1:],KS,KA))
                            aux.append(aux2)
                    return D.element(aux)
                else:
                    if isinstance(self._value[0],int):
                        aux = []
                        for i in self._value:
                            aux.append(i>>KA)
                    else:
                        aux = []
                        for i in self._value:
                            aux.append(i._push_to_subunion(A,D.subdomain(),S[1:],KS,KA))
                    if self._direction:
                        return D.Vector_type(aux)
                    else:
                        return D.Vector_type(aux,'c')
        
    #Ficar les que surgeixin.