
from PyM.PyWIT.PyECC.Zn_type import *
from PyM.PyWIT.PyECC.QQ import *
import PyM.PyWIT.PyECC.PyECC_globals as PGl

import math


class Zn(Base.Base):
    def __init__(self, base = 2, name = ''):
        self._ring = True
        self._group = True
        self._characteristic = base
        self._cardinal = base
        if name == '':
            if self._characteristic < 100:
                self._name = 'Z'+self._characteristic.__str__()
            else:
                self._name = 'Zn'
        else:
            self._name = name
    
    def cp(self):
        return Zn(self._characteristic, self._name)
    
    def element(self,a):
        if isinstance(a,Zn_type):
            if a.domain() == self:
                return a
        return Zn_type(a,self)
            
    def __repr__(self):
        if PGl.TYPES == 1:
            if(self.is_field()):
                return self._name + " :: Field"
            else:
                return self._name + " :: Ring"
        else:
             return self._name   
    
    def __str__(self):
        return self.__repr__()
    
    def _print_(self):
        return self._name
    
    def is_field(self):
        if (self._field == 0):
            aux = rabin_miller(self._characteristic)
            if aux:
                self._field = 1
            else:
                self._field = -1
        
        if self._field == 1:
            return True
        else:
            return False
    
    def _add_(self, a, b):
        if isinstance(a,int):
            res = a + b._value
            return Zn_type(res,self)
        
        if isinstance(b, int):
            res = a._value + b
            return Zn_type(res,self)
            
        if a._structure._characteristic != b._structure._characteristic:
            raise TypeError("Parameters a and b have different Base")
        
        res = a._value + b._value
        return Zn_type(res,self)        

    def _sub_(self, a, b):
        if isinstance(a,int):
            res = a - b._value
            return Zn_type(res,self)
        
        if isinstance(b, int):
            res = a._value - b
            return Zn_type(res,self)
            
        if a._structure._characteristic != b._structure._characteristic:
            raise TypeError("Parameters a and b have different Base")
        
        res = a._value - b._value
        return Zn_type(res,self)     
        
    def _mul_(self, a, b):
        if isinstance(a,int):
            res = a * b._value
            return Zn_type(res,self)
            
        if isinstance(b,int):
            res = a._value * b
            return Zn_type(res,self)
            
            
        if a._structure._characteristic != b._structure._characteristic:
            raise TypeError("Parameters a and b have different Base")
        
        res = a._value * b._value
        return Zn_type(res,self)        
    
    def is_invertible(self, a):
        if isinstance(a, int):
            return self.is_invertible(a>>self)
        if isinstance(a, Base_type.Base_type):
            if a._is_inv == -1:
                if igcd(a._value,self._characteristic) == 1:
                    a._is_inv = 1
                    return True
                else:
                    a._is_inv = 0
                    return False
            elif a._is_inv == 0:
                return False
            return True
        return False
        
    def _inv_(self, a):
        if not self.is_invertible(a):
            raise ZeroDivisionError("This number is not invertible")
        if (self.is_sub_element(a))and(not self.is_element(a)): 
            return self._inv_(a>>self)
        a1,a2,d=extended_euclidean_algorithm(a._value, self._characteristic)
        return Zn_type(a1,self)
    
    def _div_res_(self,a,b):
        dv = a//b
        return dv,a-dv*b
    
    def _div_(self, a, b):
        if isinstance(a,int):
            res = a * (self._inv_(b))._value
            return Zn_type(res,self)
            
        if isinstance(b,int):
            aux = Zn_type(b,self)
            res = a._value * (self._inv_(aux))._value
            return Zn_type(res,self)
            
            
        if a._structure._characteristic != b._structure._characteristic:
            raise TypeError("Parameters a and b have different Base")      
        
        res = a._value * (self._inv_(b))._value
        return Zn_type(res,self)
    
    def _opposite_(self, a):
        res = self._characteristic-a._value;
        return Zn_type(res,self)
        
    def __rrshift__(self, a):
        if isinstance(a, int):
            res = a % self._characteristic
            return Zn_type(res, self)
        
        if isinstance(a, Zn_type):
            if a._characteristic_()==self._characteristic:
                return a
            else:
                raise TypeError("I do not know how to project into this domain")
        if isinstance(a, list):
            res = []
            for i in a:
                if isinstance(i, int):
                    aux = a % self._characteristic
                    res.append(Zn_type(aux, self))
                if isinstance(i, Zn_type):
                    if i._characteristic_()==self._characteristic:
                        res.append(i)
                else:
                    raise TypeError("At least, one of the elements of the list cannot be converted")
            return res
        if isinstance(a,Base_type.Base_type):
            if a._characteristic_() == self._characteristic:
                return a.cp()
            elif a._characteristic_() == 0:
                return a.change_characteristic(self._characteristic)
            elif a._characteristic_() % self._characteristic == 0:
                return a.change_characteristic(self._characteristic)
            else:
                raise TypeError("This element cannot be converted")
            
        raise TypeError("This element cannot be converted")
    
    def __hash__(self):
        return PGl.HZN + self._characteristic
    
    def __eq__(self, other):
        if isinstance(other, Zn):
            if other.characteristic() == self.characteristic():
                return True
        return False
    
    def __ne__(self, other):
        return not self.__eq__(other)
        
    def __req__(self, other):
        return self.__eq__(other)
    
    def __rne__(self, other):
        return self.__ne__(other)
        
    def _order_(self, other):
        if not other.is_invertible():
            raise ZeroDivisionError("This number is not invertible")
        odr = order(other._value, self._characteristic)
        return odr
        
    def is_primitive(self,other):
        if other.is_invertible():
            if other.order() == (self.cardinal() - 1):
                return True
        return False
            
    def get_primitive_element(self):
        ordT = self.cardinal()-1
        for a in range(ordT):
            s = (a+1)>>self
            if s.order() == ordT:
                return s
    
    def characteristic(self):
        return self._characteristic
    
    def cardinal(self):
        return self._cardinal
        
    def _eq_(self, a, b):
        if isinstance(b,Zn_type):
            if a._characteristic_()!=b._characteristic_():
                return False
            if a._value == b._value:
                return True
            return False
        if isinstance(b,int):
            aux = b>>self
            return self._eq_(a,aux)
    
    def _ne_(self, a, b):
        return not self._eq_(a,b)
    
    def is_sub_element(self, a):
        if isinstance(a, int):
            return True
        if isinstance(a, Zn_type):
            if a._characteristic_()==self._characteristic:
                return True
        return False
    
    def is_sub_domain(self,other):
        if isinstance(other,Zn):
            if other._characteristic == self._characteristic:
                return True
        if isinstance(other,ZZ):
            return True
        return False            

    def is_element(self, a):
        if isinstance(a, Zn_type):
            if a._characteristic_()==self._characteristic:
                return True
        return False
                
    def _rd(self): # random polynomial of degree r
        q = self.cardinal()
        k = rd_int(0,q-1)
        #x = A._rd()
        return k>>self
    
    def _rd_nonzero(self): # random polynomial of degree r
        q = self.cardinal()
        k = rd_int(1,q-1)
        #x = A._rd()
        return k>>self
    
    def change_characteristic(self,n):
        if self.characteristic() == n:
            return self
        if (n == 0):
            return ZZ()
        cself = self.characteristic()
        if ((cself%n) == 0):
            return Zn(n)
        raise TypeError("Wrong characteristic")
    
    def transform_to_ZZ(self):
        return ZZ()
    
    def K_(self):
        return self.cp()
    
    def blow(self,a,K):
        if K == self:
            return [a]
        else:
            raise TypeError("K parameter is not a subdomain")
    
    def degree(self):
        return 0    

    def subdomain(self):
        return self
    


class ZZ(Base.Base):
    """ZZ domain """
    
    def __init__(self,name = ''):
        self._ring = True
        self._group = True
        self._field = -1
        self._characteristic = 0
        self._cardinal = 'Infinity'
        if name == '':
            self._name = 'Z'
        else:
            self._name = name
    
    def cp(self):
        return self
    
    def subdomain(self):
        return self
    
    def __repr__(self):
        return self._name + " :: Ring"
    
    def __str__(self):
        return self.__repr__()
        
    def _print_(self):
        return self._name
    
    def is_field(self):
        return False 
    
    def element(self,a):
        return self.transform_to_ZZ(a)
    
    def transform_to_ZZ(self,a=''):
        if a == '':
            return self
        if isinstance(a,int):
            return a
        if isinstance(a,Base_type.Base_type):
            return a.transform_to_ZZ()
        if a == self:
            return self
        if isinstance(a,Base.Base):
            return a.transform_to_ZZ()
        if isinstance(a,float):
            return int(a)
        if isinstance(a,list):
            aux = []
            for i in a:
                aux.append(self.transform_to_ZZ(i))
            return aux
    
    def _div_(self,a,b):
        if b!=1 or b!=-1:
            D = QQ(self)
            return QQ_type(D,a,b)
        else:
            return a//b
                
    def __rrshift__(self, a):
        return self.transform_to_ZZ(a)
            
    def characteristic(self):
        return self._characteristic
        
    def cardinal(self):
        return self._cardinal
    
    def __hash__(self):
        return PGl.HZN + 10
        
    def __eq__(self, other):
        if isinstance(other,Base.Base):
            if other.characteristic() == 0:
                if other._poly_structure == False:
                    if other._matrix_structure == False:
                        if not other._fraction_structure:
                            if not other._power_structure:
                                return True
        return False
        
    def __neq__(self, other):
        return not self.__eq__(other)
    
    def is_sub_element(self, a):
        if isinstance(a,int):
            return True
        return False
    
    def base_domain(self):
        return self
        
    def is_sub_domain(self,other):
        return self.__eq__(other)
    
    def is_element(self, a):
        if isinstance(a, int):
            return True
        return False
    
    def is_invertible(self, a):
        if a == 1:
            return True
        if a == -1:
            return True
        return False
    
    def _inv_(self, a):
        if a == 1:
            return 1
        if a == -1:
            return -1
        raise TypeError("This element is not invertible")
    #Ficar les que surgeixin.    
        
    def change_characteristic(self,n):
        if n == 0:
            return self
        else:
            return Zn(n)
    
    def cp(self):
        return ZZ(self._name)
    
    def K_(self):
        return self.cp()

    def blow(self,a,K):
        if K == self:
            return [a]
        else:
            raise TypeError("K parameter is not a subdomain")
    
    def gcd(self,a,b):
        return igcd(a,b)
    
    def can_be_compared(self,other):
        return True
    
    def _rd(self):
        p = rd_float(0, 1)
        l = math.ceil(-math.log2(p))
        s = 2*rd_int()-1
        return s*ri(l)

    def _div_res_(self,a,b):
        dv = a//b
        return dv,a-dv*b