import PyM.PyWIT.PyECC.pfactorFunctions as pfactorFunctions
from PyM.PyWIT.PyECC.Fq_type import *
from PyM.PyWIT.PyECC.Matrix import *
from PyM.PyWIT.PyECC.Poly import *
import PyM.PyWIT.PyECC.PyECC_globals as PGl

import time


class Fq(Base.Base):
    """Structure of a non prime finite field"""
    def __init__(self, subDomain, poly,  symbol = 'x', name = ''):
        self._ring = True
        self._group = True
        self._characteristic = subDomain.characteristic()
        self._poly = []
        self._symbol = symbol
        self._subDomain = subDomain
        self._base_Domain = self._subDomain.base_domain()
        self._ext_Domain = True
        self.__create_poly__(poly,subDomain)
        if self._subDomain.is_field():
            if pfactorFunctions.is_irreducible(self._poly,subDomain):
                self._field = 1
            else:
                self._field = -1
        else:
            self._field = -1
        if subDomain.characteristic() == 0:
            self._cardinal = 'Infinity'
        else:
            self._cardinal = subDomain.cardinal()**(self._poly.degree())
        if name == '':
            self._name = subDomain._print_()+'['+symbol+']'
        else:
            self._name = name
    
    def cp(self):
        return Fq(self._subDomain, self._poly.cp(), self._symbol, self._name)
            
    def degree(self):
        """returns the degree of the extension"""
        return self._poly.degree()
            
    def variable(self):
        """returns the basic element of the extension"""
        return Fq_type([1],self,[1])            
    
    def __create_poly__(self, poly, domain):
        """creates the polynomial that defines the extension"""
        self._poly = []
        if isinstance(poly,list):
            if len(poly)<3:
                raise TypeError("Polynomial degree must be at least 2")
            if poly[0]!=1:
                raise TypeError("Polynomial must be monic")
            aux = []
            aindex = []
            x = Symbol_type('x')
            P = Poly(self.subdomain())
            for i in range(len(poly)):
                s = poly[i]>>self._subDomain
                if s!= 0:
                    aux.append(Monomial_type(s,[x],[len(poly)-i-1],P))
            if len(aux) == 0:
                self._poly = Monomial_type(0>>self.subdomain(),[x],[0],P)
            elif len(aux) == 1:
                self._poly = aux[0]
            else:
                self._poly = Poly_type(aux,[x],P,aux[0]._degree,aux[0]._degree)
        elif isinstance(poly,Poly_type):
            if poly._structure._subDomain == self._subDomain:
                if len(poly._variables) > 1:
                    raise TypeError("Polynomial must be univariate")
                if poly.degree()<2:
                    raise TypeError("Polynomial degree must be at least 2")
                if poly.leading_coeff()!=1:
                    raise TypeError("Polynomial must be monic")
                self._poly = poly
            elif domain.is_sub_domain(poly.K_()):
                if len(poly._variables) > 1:
                    raise TypeError("Polynomial must be univariate")
                if poly.degree()<2:
                    raise TypeError("Polynomial degree must be at least 2")
                if poly.leading_coeff()!=1:
                    raise TypeError("Polynomial must be monic")
                z = 0>>domain
                self._poly = poly + z
            else:
                raise TypeError("This polynomial space is not valid")
        else:
            raise TypeError("Polynomial must be a list of at least 2 coefficients with the first element 1")
            
    def _reduct_(self,value,index):
        """returns the polynomial to this space"""
        if index[0] < self.degree():
            return (value,index)
        aux = []
        for i in range(len(value)):
            if value[i]!=0:
                aux.append(Monomial_type(value[i], self._poly._variables[:],[index[i]],self._poly.domain()))
        if len(aux) == 0:
            return ([0>>self.subdomain()],[0])
        elif len(aux) == 1:
            aux = aux[0]
        else:
            aux = Poly_type(aux,self._poly._variables[:],self._poly.domain(),aux[0]._degree, aux[0]._degree)
        res = self._poly.domain()._md_reduce_core_(aux,self._poly)
        if isinstance(res,Monomial_type):
            return ([res._value],res._index)
        else:
            vals = []
            inds = []
            for i in res._value:
                vals.append(i._value)
                inds.append(i._index[0])
        return (vals,inds)                        
    
    def element(self, lst, ind = ''):
        """returns the element defined by lst and ind"""
        aux =[]
        if ind == '':
            if isinstance(lst, list):
                for i in lst:
                    if (not self._subDomain.is_sub_element(i)):
                        raise TypeError("This list cannot be converted to Fq_type")
                    aux.append(i>>self._subDomain)
            elif isinstance(lst,Base_type):
                if self.is_element(lst):
                    return lst
                elif self.is_sub_element(lst):
                    aux.append(lst>>self._subDomain)
                else:
                    raise TypeError("This element cannot be converted to Fq_type")
            return Fq_type(aux, self)
        else:
            return Fq_type(lst,self,ind)

    def _print_(self):
        return self._name
    
    def __str__(self):
        if PGl.TYPES == 1:
            if(self.is_field()):
                
                return self._name + " :: Field"
            else:
                return self._name + " :: Ring"
        else:
            return self._name
    
    def __repr__(self):
        return self.__str__()
    
    def is_field(self):
        """returns if this space is a field"""
        if self._field == 1:
            return True
        return False
    
    def __rrshift__(self, a):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return self.__rrshift__(aux)
        if (self.is_element(a)):
            return a
        if (self.is_sub_element(a)):
            aux = []
            aux.append(a>>self._subDomain)
            aindex = [0]
            return Fq_type(aux,self,aindex)
        if isinstance(a,Symbol_type):
            return Monomial_type(1>>self,[a.cp()],[1],Poly(self))
        if isinstance(a, Base_type.Base_type):
            if a.K_() == self:
                return a.cp()
            if a.K_().is_sub_domain(self):
                b = a.transform_to_min_domain()
                K = b.K_()
                if K == self:
                    return b
                elif K.is_sub_domain(self):
                    return b
                else:
                    return b >> self
            elif self.is_sub_domain(a.K_()):
                if isinstance(a,Monomial_type):
                    if a._varOrder == None:
                        s = a._value >> self
                        if s != 0:
                            return Monomial_type(s, a._variables[:], a._index[:],Poly(self))
                        else:
                            return Monomial_type(s, a._variables[:], [0]*len(a._variables),Poly(self))
                    else:
                        s = a._value >> self
                        if s != 0:
                            return Monomial_type(s, a._variables[:], a._index[:],Poly(self),a._varOrder[:])
                        else:
                            return Monomial_type(s, a._variables[:], [0]*len(a._variables),Poly(self))
                        
                if isinstance(a,Poly_type):
                    aux = []
                    md = 0
                    for i in a._value:
                        if a._varOrder == None:
                            s = i._value >> self
                            if s != 0:
                                aux.append(Monomial_type(s, i._variables[:], i._index[:],Poly(self)))
                                if aux[-1]._totIndex > md:
                                    md = aux[-1]._totIndex
                        else:
                            s = i._value >> self
                            if s != 0:
                                aux.append(Monomial_type(s, i._variables[:], i._index[:],Poly(self), a._varOrder[:]))
                                if aux[-1]._totIndex > md:
                                    md = aux[-1]._totIndex
                    if len(aux) == 0:
                        return Monomial_type(0>>self, a._variables[:], [0]*len(a._variables),Poly(self))
                    elif len(aux) == 1:
                        return aux[0]
                    else:
                        if a._varOrder == None:
                            return Poly_type(aux,a._variables[:],Poly(self),aux[0]._degree,md)
                        else:
                            aord = _get_own_order_monomials_(aux[:])
                            return Poly_type(aux,a._variables[:],Poly(self),aux[0]._degree,md,aord,a._varOrder[:])
                if isinstance(a,Vector_type):
                    e = a._value >> self
                    if len(e) == 0:
                        if a._direction:
                            return Vector_type(Matrix(self),[],'r')
                        else:
                            return Vector_type(Matrix(self),[],'c')
                    D = e[0].domain()
                    for i in range(1,len(e)):
                        D = D.union_domain(e[i].domain())
                    aux = []
                    for b in e:
                        aux.append(b>>D)
                    if a._direction:
                        return Vector_type(Matrix(D),aux,'r')
                    else:
                        return Vector_type(Matrix(D),aux,'c')
                if isinstance(a,Matrix_type):
                    e = [k >> self for k in a._value]
                    if len(e) == 0:
                        return Matrix_type(Matrix(self),[[]])
                    D = e[0][0].domain()
                    for i in range(1,len(e[0])):
                        D = D.union_domain(e[0][i].domain())
                    for i in range(1,len(e)):
                        for j  in range(len(e[i])):
                            D = D.union_domain(e[i][j].domain())
                    aux = []
                    for i in range(len(e)):
                        aux2 = []
                        for j in range(len(e[i])):
                            aux2.append(e[i][j]>>D)
                        aux.append(aux2)
                    return Matrix_type(Matrix(D),aux)
                
                    
        if isinstance(a, list):
            aux = []
            for i in a:
                aux.append(i>>self)
            return aux
        else:
            raise TypeError("This element cannot be converted into this domain")

    def variable_d(self,n):
        """returns the monomial x**n"""
        aux = [1]
        aindex = [n]
        return Fq_type(aux,self,aindex)
    
    def _leading_term_(self,f):
        aux = [f._value[0]]
        aindex = [f._index[0]]
        return Fq_type(aux,self,aindex)
        
    def _leading_coef_(self,f):
        return f._value[0]
    
    def _add_(self, a, b):
        """operates two elements using this domain operation +"""
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux+b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a+aux
        azero = 0>>self._subDomain
        if (self.is_sub_element(a))and(not self.is_element(a)):
            aux = []
            aindex = []
            for i in range(len(b._value)):
                aux.append((b._value[i])>>self._subDomain)
                aindex.append((b._index[i]))
            if aindex[len(aindex)-1] == 0:
                aux[len(aux)-1] = (a + aux[len(aux)-1] )>>self._subDomain
                if aux[len(aux)-1]==azero:
                    aux = aux[:len(aindex)-1]
                    aindex = aindex[:len(aindex)-1]
            else:
                s = a>>self._subDomain
                if s != azero:
                    aux.append(s)
                    aindex.append(0)
            if (len(aux) == 0):
                aux = [0]
                aindex = [0]
            return Fq_type(aux,self,aindex)
        if (self.is_sub_element(b))and(not self.is_element(b)):
            aux = []
            aindex = []
            for i in range(len(a._value)):
                aux.append((a._value[i])>>self._subDomain)
                aindex.append((a._index[i]))
            if aindex[len(aindex)-1] == 0:
                aux[len(aux)-1] = (aux[len(aux)-1] + b)>>self._subDomain
                if aux[len(aux)-1]==azero:
                    aux = aux[:len(aindex)-1]
                    aindex = aindex[:len(aindex)-1]
            else:
                s = b>>self._subDomain
                if s != azero:
                    aux.append(s)
                    aindex.append(0)
            if (len(aux) == 0):
                aux = [0]
                aindex = [0]
            return Fq_type(aux,self,aindex)
        
        aux = []
        aindex = []
        ia = 0
        ib = 0
        
        while True:
            if a._index[ia]>b._index[ib]:
                aux.append(a._value[ia])
                aindex.append(a._index[ia])
                ia = ia + 1
                if ia >= len(a._index):
                    ia = -1
                    break
            elif a._index[ia]<b._index[ib]:
                aux.append(b._value[ib])
                aindex.append(b._index[ib])
                ib = ib + 1
                if ib >= len(b._index):
                    ib = -1
                    break
            else:
                s = a._value[ia] + b._value[ib]
                if s!= azero:
                    aux.append(s)
                    aindex.append(a._index[ia])
                ia = ia + 1
                ib = ib + 1
                if ia >= len(a._index):
                    ia = -1
                    break
                if ib >= len(b._index):
                    ib = -1
                    break
        if ia == -1:
            if ib != -1:
                for i in range(ib,len(b._index)):
                    aux.append(b._value[i])
                    aindex.append(b._index[i])
        if ib == -1:
            if ia != -1:
                for i in range(ia,len(a._index)):
                    aux.append(a._value[i])
                    aindex.append(a._index[i])
        if (len(aux) == 0):
            aux = [0]
            aindex = [0]
        return Fq_type(aux,self,aindex)
        
            
    def _sub_(self, a, b):
        """operates two elements using this domain operation -"""
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux-b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a-aux
        azero = 0>>self._subDomain
        if (self.is_sub_element(a))and(not self.is_element(a)):
            aux = []
            aindex = []
            for i in range(len(b._value)):
                aux.append((-b._value[i])>>self._subDomain)
                aindex.append((b._index[i]))
            if aindex[len(aindex)-1] == 0:
                aux[len(aux)-1] = (a + aux[len(aux)-1] )>>self._subDomain
                if aux[len(aux)-1]==azero:
                    aux = aux[:len(aindex)-1]
                    aindex = aindex[:len(aindex)-1]
            else:
                s = a>>self._subDomain
                if s != azero:
                    aux.append(s)
                    aindex.append(0)
            if (len(aux) == 0):
                aux = [0]
                aindex = [0]
            return Fq_type(aux,self,aindex)
        if (self.is_sub_element(b))and(not self.is_element(b)):
            aux = []
            aindex = []
            for i in range(len(a._value)):
                aux.append((a._value[i])>>self._subDomain)
                aindex.append((a._index[i]))
            if aindex[len(aindex)-1] == 0:
                aux[len(aux)-1] = (aux[len(aux)-1] - b)>>self._subDomain
                if aux[len(aux)-1]==azero:
                    aux = aux[:len(aindex)-1]
                    aindex = aindex[:len(aindex)-1]
            else:
                s = -b>>self._subDomain
                if s != azero:
                    aux.append(s)
                    aindex.append(0)
            if (len(aux) == 0):
                aux = [0]
                aindex = [0]
            return Fq_type(aux,self,aindex)
        aux = []
        aindex = []
        ia = 0
        ib = 0
        
        while True:
            if a._index[ia]>b._index[ib]:
                aux.append(a._value[ia])
                aindex.append(a._index[ia])
                ia = ia + 1
                if ia >= len(a._index):
                    ia = -1
                    break
            elif a._index[ia]<b._index[ib]:
                aux.append(-b._value[ib])
                aindex.append(b._index[ib])
                ib = ib + 1
                if ib >= len(b._index):
                    ib = -1
                    break
            else:
                s = a._value[ia] - b._value[ib]
                if s != azero:
                    aux.append(s)
                    aindex.append(a._index[ia])
                ia = ia + 1
                ib = ib + 1
                if ia >= len(a._index):
                    ia = -1
                    break
                if ib >= len(b._index):
                    ib = -1
                    break
        if ia == -1:
            if ib != -1:
                for i in range(ib,len(b._index)):
                    aux.append(-b._value[i])
                    aindex.append(b._index[i])
        if ib == -1:
            if ia != -1:
                for i in range(ia,len(a._index)):
                    aux.append(a._value[i])
                    aindex.append(a._index[i])
        if (len(aux) == 0):
            aux = [0]
            aindex = [0]
        return Fq_type(aux,self,aindex)
        
    def _mul_(self, a, b):
        """operates two elements using this domain operation *"""
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux*b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a*aux
        azero = 0>>self._subDomain
        if (self.is_sub_element(a))and(not self.is_element(a)):
            aux = []
            aindex = []
            for i in range(len(b._value)):
                s = (a*b._value[i])>>self._subDomain
                if s!= azero:
                    aux.append(s)
                    aindex.append(b._index[i])
            if (len(aux) == 0):
                aux = [0]
                aindex = [0]
            return Fq_type(aux,self,aindex)
        if (self.is_sub_element(b))and(not self.is_element(b)):#and(not self._subDomain.is_element(b)):
            aux = []
            aindex = []
            for i in range(len(a._value)):
                s = (a._value[i]*b)>>self._subDomain
                if s != azero:
                    aux.append(s)
                    aindex.append(a._index[i])
            if (len(aux) == 0):
                aux = [0]
                aindex = [0]
            return Fq_type(aux,self,aindex)
        aux = []
        aindex = []
        for i in range(len(a._value)):
            posini = 0
            for j in range(len(b._value)):
                inij = a._index[i] + b._index[j]
                trobat = False
                pos = 0
                for k in range(posini,len(aux)):
                    if aindex[k] < inij:
                        pos = k
                        break
                    elif aindex[k] == inij:
                        pos = k
                        trobat = True
                        break
                    elif k >= len(aux) - 1:
                        pos = len(aux)
                        break
                posini = pos
                if trobat:
                    aux[pos] = aux[pos] + a._value[i] * b._value[j]
                    if aux[pos] == azero:
                        aux.pop(pos)
                        aindex.pop(pos)
                        posini = posini - 1
                else:
                    s = a._value[i] * b._value[j]
                    if s != azero:
                        aux.insert(pos,s)
                        aindex.insert(pos,inij)
        if (len(aux) == 0):
            aux = [0]
            aindex = [0]
        return Fq_type(aux,self,aindex)        

    def _opposite_(self, a):
        aux = []
        for i in a._value:
            aux.append((-i)>>self._subDomain)
        return Fq_type(aux,self,a._index)
    
    def __hash__(self):
        return PGl.HFQ + hash(self._subDomain) + hash(self._poly)
    
    def __eq__(self, other):
        if isinstance(other, Base.Base):
            if other._ext_Domain:
                if other._subDomain == self._subDomain:
                    if self._poly.degree() == other._poly.degree():
                        return self._poly == other._poly
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
        if self.is_element(b):
            if (len(a._value)!= len(b._value)):
                return False
            for i in range(len(a._value)):
                if (a._value[i]!=b._value[i]):
                    return False
                if (a._index[i]!=b._index[i]):
                    return False
            return True
        else:
            aux = b>>self
            return self._eq_(a,aux)
            
    def characteristic(self):
        """Return the characteristic of this domain. """
        return self._characteristic

    def cardinal(self):
        """Return the cardinal of this domain. """
        return self._cardinal
        
    def _rd(self,r=''): 
        """random element"""
        if r == '':
            q = self.cardinal()
            k = rd_int(0,q-1)
            #x = A._rd()
            return pfactorFunctions.pick_element(k,self) 
        else:
            A = self.subdomain()
            q = A.cardinal()
            C = []
            for k in range(r,-1,-1):
                x = A._rd()
                C += [x]
            return C>>self
    
    def _rd_nonzero(self,r=''): 
        """random element"""
        if r == '':
            q = self.cardinal()
            k = rd_int(1,q-1)
            #x = A._rd()
            return pfactorFunctions.pick_element(k,self) 
        else:
            A = self.subdomain()
            q = A.cardinal()
            k = rd_int(1,q**(r+1)-1)
            #x = A._rd()
            return pfactorFunctions.pick_element(k,self)
         
    def _inv_(self,a):
        """returns the inverse of the element a in this space"""
        P = self._poly._structure
        aux = a>>self
        vals = []
        for i in range(len(a._value)):
            vals.append(Monomial_type(a._value[i],self._poly._variables[:],[a._index[i]],self._poly.domain()))
        if len(vals) == 0:
            g = Monomial_type(0>>self.subdomain(),self._poly._variables[:],[0],self._poly.domain())
        elif len(vals) == 1:
            g = vals[0]
        else:
            g = Poly_type(vals,self._poly._variables[:],self._poly.domain(),vals[0]._degree,vals[0]._degree)
        [a1,a2,a3] = P.extended_euclidean_algorithm(self._poly,g)
        if a1 != 1:
            a._is_inv == 0
            raise TypeError("Not Invertible")
        else:
            a._is_inv = 1
            if len(a3._variables) == 0:
                return Fq_type(a3._value,self,[0])
            if isinstance(a3,Monomial_type):
                return Fq_type([a3._value],self,a3._index)
            else:
                vals = []
                inds = []
                for i in range(len(a3._value)):
                    vals.append(a3._value[i]._value)
                    inds.append(a3._value[i]._index[0])
            return Fq_type(vals,self,inds)
    
    def is_invertible(self,a):
        """returns if a is invertible in this space"""
        if isinstance(a,float):
            return self.is_invertible(a>>self)
        if isinstance(a,int):
            return self.is_invertible(a>>self)
        if isinstance(a, Base_type.Base_type):
            if a._is_inv == -1:
                try:
                    self._inv_(a)
                except TypeError:
                    a._is_inv = 0
                    return False
                a._is_inv = 1
                return True
            elif a._is_inv == 0:
                return False
            return True
        return False
            
    def _div_(self,a,b):
        """operates two elements using this domain operation /"""
        return a*self._inv_(b)
    
    def _div_res_(self,a,b):
        dv = a//b
        return dv,a-dv*b
    
    def _order_(self,other):
        """returns the order of the element"""
        if not other.is_invertible():
            raise ZeroDivisionError("This number is not invertible")
        D = divisors(self.cardinal()-1)
        one = 1>>self
        for d in D:
            if other**d == one:
                return d
    
    def is_primitive(self,other):
        """returns if the element is primitive or not"""
        if other.is_invertible():
            one = 1>>self
            ordT = self.cardinal() - 1
            P = prime_factors(ordT)
            for p in P:
                if (other **(ordT//p)) == one:
                    return False
            return True
        return False
                
    def get_primitive_element(self):
        """returns a primitive element of the extension"""        
        a = self.subdomain().cardinal()
        while True:
            s = pfactorFunctions.pick_element(a,self)
            if s.is_primitive():
                return s
            a = a + 1
    
    def get_polynomial(self):
        """returns the polynomial that defines the extension"""
        return self._poly.cp()
        
    def change_characteristic(self,n):
        """changes the characteristic of this space to n and transforms all objects that define this space"""
        if self.characteristic() == n:
            return self
        if (n == 0):
            return self.transform_to_ZZ()
        cself = self.characteristic()
        if ((cself%n) == 0):
            SDom = self.subdomain().change_characteristic(n)
            poly = self._poly.change_characteristic(n)
            return Fq(SDom, poly)
        raise TypeError("Wrong characteristic")
    
    def transform_to_ZZ(self):
        """changes the characteristic of this space to 0 and transforms all objects that define this space"""
        if self.characteristic() == 0:
            return self
        poly = self._poly.transform_to_ZZ()
        SDom = self.subdomain().transform_to_ZZ()
        return Fq(SDom,poly)
    
    def set_symbol(self,symbol):
        """changes the symbol of the space"""
        if self._name == (self.subdomain()._print_()+'['+self._symbol+']'):
            self._name = self.subdomain()._print_()+'['+symbol+']'
        self._symbol = symbol
    
    def get_symbol(self):
        """returns the symbol of the space"""
        return self._symbol
    
    def K_(self):
        return self.cp()
    
    def blow(self,a,K):
        """Returns a list corresponding to a from self into a K vector (K\subset self)"""
        if not self.is_sub_domain(K):
            raise TypeError("Parameter K must be a subdomain")
        if K == self:
            return [a]
        res = a.to_list()
        if len(res)<self.degree():
            res = res+[0>>self.subdomain()]*(self.degree() - len(res))
        ni = len(res)-1
        for i in range(ni,-1,-1):
            res = res[0:i] + self.subdomain().blow(res[i],K) + res[i+1:]
        return res    
        
            
        
        
        