# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:08:09 2017

@author: admin_local
"""
from PyM.PyWIT.PyECC.Matrix_type import *
from PyM.PyWIT.PyECC.Zn import *
import PyM.PyWIT.PyECC.PyECC_globals as PGl


#Moure a pfactor o a algun altre lloc (crida abançada vector)        

class Matrix(Base.Base):
    def __init__(self, base, name = ''):
        self._ring = False
        self._group = False
        self._field = -1
        self._matrix_structure = True
        self._characteristic = base.characteristic()
        self._cardinal = 'Infinity'
        self._subDomain = base
        self._base_Domain = base.base_domain()
        if name == '':
            self._name = 'Mat['+base._print_()+']'
        else:
            self._name = name
    
    def cp(self):
        return Matrix(self._subDomain, self._name)
            
    def eye(self,n,m=''):
        if m == '':
            m = n
        if m>n:
            nd = n
        else:
            nd = m
        aux = [[(0>>self.subdomain())]*m for i in range(n)]
        for i in range(nd):
            aux[i][i] = 1>>self.subdomain()
        D = Matrix(self.subdomain(),self._name)
        return Matrix_type(D,aux)
    
    def e_i(self,k,n,d = 'r'):
        aux = [0>>self.subdomain()]*n
        aux[k] = 1>>self.subdomain()
        D = Matrix(self.subdomain(),self._name)
        if d == 'r':
            return Vector_type(D,aux)
        else:
            return Vector_type(D,aux,'c')
        
    def characteristic(self):
        return self._characteristic
    
    def cardinal(self):
        return self._cardinal
        
    def element(self, lst, ind = ''):
        if isinstance(lst,Matrix_type):
            if lst.domain() == self:
                return lst.cp()
        if isinstance(lst,Vector_type):
            if lst.domain() == self:
                return lst
        D = Matrix(self.subdomain(),self._name)
        return Matrix_type(D,lst)
    
    def Vector_type(self,lst,direction = 'r'):
        D = Matrix(self.subdomain(),self._name)
        return Vector_type(D,lst,direction)
        
    def _print_(self):
        return self._name
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._name + ' :: Array'
        else:
            return self._name
    
    def __repr__(self):
        return self.__str__()
    
    def is_field(self):
        return False
    
    def _create_Matrix_(self,D):
        Dom = Matrix(D)
        return Dom
    
    def _ZZ_(self):
        return ZZ()
        
    def _div_(self, a, b):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux/b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a/aux
            
        #Treure els elements que de veritat són del subdomini
        if isinstance(a,Matrix_type):
            if a.ncols() == 1 and a.nrows() == 1:
                return self._div_(a[0][0],b) 
        if isinstance(b,Matrix_type):
            if b.ncols() == 1 and b.nrows() == 1:
                return self._div_(a,b[0][0]) 
        if isinstance(a,Vector_type):
            if a.size()==1:
                return self._div_(a[0],b) 
        if isinstance(b,Vector_type):
            if b.size()==1:
                return self._div_(a,b[0]) 
        if self.subdomain().is_field():
            D = self
        else:
            aux = QQ(self.subdomain())
            D = Matrix(aux)
        if self.is_invertible(b,D):
            return a*self._inv_(b)
        else:
            raise TypeError("Elment not invertible")
        
    def _opposite_(self, a):
        if isinstance(a,Matrix_type):
            aux = [a._value[i][:] for i in range(len(a._value))]
            for i in range(a._rows):
                for j in range(a._cols):
                    aux[i][j] = -aux[i][j]
            D = Matrix(self.subdomain(),self._name)
            return Matrix_type(D,aux)
        if isinstance(a,Vector_type):
            aux2 = a._value[:]
            aux = [-i for i in aux2]
            D = Matrix(self.subdomain(),self._name)
            if a._direction:
                return Vector_type(D,aux,'r')
            else:
                return Vector_type(D,aux,'c')
                
           
    def __hash__(self):
        return PGl.HMAT + hash(self.subdomain())
   
    def __eq__(self, other):
        if isinstance(other, Base.Base):
            if other._subDomain == self._subDomain:
                if other._matrix_structure:
                    return True
        return False
        
    def __neq__(self, other):
        return not self__eq__(other)
         
    def _eq_(self, a, b):
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return self._eq_(aux,b)
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return self._eq_(a,aux)
        if self.is_element(b):
            if isinstance(a,Vector_type):
                if isinstance(b,Vector_type):
                    if a._value == b._value:
                        return True
                else:
                    if len(a._value) == 0:
                        if len(b._value)*len(b._value[0]) == 0:
                            return True
                    elif a._direction:
                        if b._rows == 1:
                            if b._value[0] == a._value:
                                return True
                    else:
                        if b._cols == 1:
                            for i in range(b._rows):
                                if b._value[0][i] != a._value[i]:
                                    return False
                            return True
            else:
                if isinstance(b,Matrix_type):
                    if len(a._value) == len(b._value):
                        if a._rows == b._rows:
                            if a._rows*a._cols == 0:
                                return True
                            else:
                                for i in range(a._rows):
                                    if a._value[i] != b._value[i]:
                                        return False
                                return True
                else:
                    if len(b._value) == 0:
                        if len(a._value)*len(a._value[0]) == 0:
                            return True
                    elif b._direction:
                        if a._rows == 1:
                            if a._value[0] == b._value:
                                return True
                    else:
                        if a._cols == 1:
                            for i in range(a._rows):
                                if a._value[0][i] != b._value[i]:
                                    return False
                            return True            
                        
                        
            if a._value == b._value:
                if isinstance(a,Matrix_type):
                    return True
                if isinstance(b,Matrix_type):
                    return True
                if a._direction == b._direction:
                    return True
            return False
        else:
            try:
                aux = b>>self
                return self._eq_(a,aux)
            except:
                return False
  
    def __rrshift__(self, a):
        if isinstance(a,list):
            D = Matrix(self.subdomain(),self._name)
            return Matrix_type(D,a)
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return self.__rrshift__(aux)
        if isinstance(a,int):
            D = Matrix(self.subdomain(),self._name)
            return Matrix_type(D,[a])
        if isinstance(a,Matrix_type):
            if a.domain() == self:
                return a.cp()
            if self.is_sub_domain(a.domain().subdomain()):
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,a._value)
            return a.cp()
        if isinstance(a,Vector_type):
            if a.domain() == self:
                return a.cp()
            if self.is_sub_domain(a.domain().subdomain()):
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,a._value)
                else:
                    return Vector_type(D,a._value,'c')
            return a.cp()
        if isinstance(a,Base_type.Base_type):
            D = Matrix(self.subdomain(),self._name)
            return Matrix_type(D,[a])
    
        
    def is_invertible(self, a, D = ''):
        if D == '':
            D = self
        if D == self:
            if isinstance(a,Vector_type):
                for i in a._value:
                    if isinstance(i,int):
                        if i*i != 1:
                            return False
                    elif not i.is_invertible():
                        return False
                return True
            if isinstance(a,Matrix_type):
                s = self._det_(a)
                if isinstance(s,int):
                    if s == 1 or s == -1:
                        return True
                    else:
                        return False
                else:
                    return s.is_invertible()
            if isinstance(a,float):
                aux = a>>QQ(self.base_domain())
                return self.is_invertible(aux)
            if isinstance(a,int):
                if a == 1 or a== -1:
                    return True
            else:
                return (a>>self.subdomain()).is_invertible()
        else:
            aux = a>>D
            return D.is_invertible(a>>D)
        
    def _add_(self, a, b):
        """operates two elements using this domain operation +"""
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux+b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a+aux
            
        #Treure els elements que de veritat són del subdomini
        if isinstance(a,Matrix_type):
            if a.ncols() == 1 and a.nrows() == 1:
                return self._add_(a[0][0],b) 
        if isinstance(b,Matrix_type):
            if b.ncols() == 1 and b.nrows() == 1:
                return self._add_(a,b[0][0]) 
        if isinstance(a,Vector_type):
            if a.size()==1:
                return self._add_(a[0],b) 
        if isinstance(b,Vector_type):
            if b.size()==1:
                return self._add_(a,b[0])               
            
            
        if isinstance(a,Matrix_type):
            if isinstance(b,Matrix_type):
                if a._rows == b._rows:
                    if a._cols == b._cols:
                        aux = []
                        for i in range(a._rows):
                            aux2 = []
                            for j in range(a._cols):
                                aux2.append(a._value[i][j]+b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                raise TypeError("Matrixes have different sizes")
            elif isinstance(b,Vector_type):
                if b._direction:
                    if b.size() == a._cols:
                        aux = []
                        for i in range(a._rows):
                            aux2=[]
                            for j in range(a._cols):
                                aux2.append(a._value[i][j]+b._value[j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
                if not b._direction:
                    if b.size() == a._rows:
                        aux = []
                        for i in range(a._rows):
                            aux2=[]
                            for j in range(a._cols):
                                aux2.append(a._value[i][j]+b._value[i])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
            elif isinstance(b,Base_type.Base_type):
                aux = []
                for i in range(a._rows):
                    aux2 = []
                    for j in range(a._cols):
                        aux2.append(a._value[i][j]+b)
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
            elif isinstance(b,int):
                aux = []
                for i in range(a._rows):
                    aux2 = []
                    for j in range(a._cols):
                        aux2.append(a._value[i][j]+b)
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
        if isinstance(a,Vector_type):
            if isinstance(b,Matrix_type):
                if a._direction:
                    if a.size() == b._cols:
                        aux = []
                        for i in range(b._rows):
                            aux2 = []
                            for j in range(b._cols):
                                aux2.append(a._value[j]+b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Vector and Matrix have different sizes")
                if not a._direction:
                    if a.size() == b._rows:
                        aux = []
                        for i in range(b._rows):
                            aux2 = []
                            for j in range(b._cols):
                                aux2.append(a._value[i]+b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                raise TypeError("Vector and Matrix have different sizes")
            elif isinstance(b,Vector_type):
                if b.size() == a.size():
                    aux = []
                    for i in range(a.size()):
                        aux.append(a._value[i]+b._value[i])
                    D = Matrix(self.subdomain(),self._name)
                    if a._direction == b._direction:
                        if a._direction:
                            return Vector_type(D,aux,'r')
                        else:
                            return Vector_type(D,aux,'c')
                    return Vector_type(D,aux,'r')
                raise TypeError("Vectors size not match")
            elif isinstance(b,Base_type.Base_type):
                aux = []
                for i in range(a.size()):
                    aux.append(a._value[i]+b)
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
            elif isinstance(b,int):
                aux = []
                for i in range(a.size()):
                    aux.append(a._value[i]+b)
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
                    
        if isinstance(b,Matrix_type):
            if isinstance(a,Vector_type):
                if a._direction:
                    if a.size() == b._cols:
                        aux = []
                        for i in range(b._rows):
                            aux2=[]
                            for j in range(b._cols):
                                aux2.append(a._value[j]+b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
                if not a._direction:
                    if a.size() == b._rows:
                        aux = []
                        for i in range(b._rows):
                            aux2=[]
                            for j in range(b._cols):
                                aux2.append(a._value[i]+b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
            elif isinstance(a,Base_type.Base_type):
                aux = []
                for i in range(b._rows):
                    aux2 = []
                    for j in range(b._cols):
                        aux2.append(a+b._value[i][j])
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
            elif isinstance(a,int):
                aux = []
                for i in range(b._rows):
                    aux2 = []
                    for j in range(b._cols):
                        aux2.append(a+b._value[i][j])
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
        if isinstance(b,Vector_type):
            if isinstance(a,Base_type.Base_type):
                aux = []
                for i in range(b.size()):
                    aux.append(a+b._value[i])
                D = Matrix(self.subdomain(),self._name)
                if b._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
            elif isinstance(a,int):
                aux = []
                for i in range(b.size()):
                    aux.append(a+b._value[i])
                D = Matrix(self.subdomain(),self._name)
                if b._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
        else:
            return a+b
    
    def _sub_(self, a, b):
        """operates two elements using this domain operation -"""
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux-b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a-aux
            
        #Treure els elements que de veritat són del subdomini
        if isinstance(a,Matrix_type):
            if a.ncols() == 1 and a.nrows() == 1:
                return self._sub_(a[0][0],b) 
        if isinstance(b,Matrix_type):
            if b.ncols() == 1 and b.nrows() == 1:
                return self._sub_(a,b[0][0]) 
        if isinstance(a,Vector_type):
            if a.size()==1:
                return self._sub_(a[0],b) 
        if isinstance(b,Vector_type):
            if b.size()==1:
                return self._sub_(a,b[0])             
            
            
        if isinstance(a,Matrix_type):
            if isinstance(b,Matrix_type):
                if a._rows == b._rows:
                    if a._cols == b._cols:
                        aux = []
                        for i in range(a._rows):
                            aux2 = []
                            for j in range(a._cols):
                                aux2.append(a._value[i][j]-b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                raise TypeError("Matrixes have different sizes")
            elif isinstance(b,Vector_type):
                if b._direction:
                    if b.size() == a._cols:
                        aux = []
                        for i in range(a._rows):
                            aux2=[]
                            for j in range(a._cols):
                                aux2.append(a._value[i][j]-b._value[j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
                if not b._direction:
                    if b.size() == a._rows:
                        aux = []
                        for i in range(a._rows):
                            aux2=[]
                            for j in range(a._cols):
                                aux2.append(a._value[i][j]-b._value[i])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
            elif isinstance(b,Base_type.Base_type):
                aux = []
                for i in range(a._rows):
                    aux2 = []
                    for j in range(a._cols):
                        aux2.append(a._value[i][j]-b)
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
            elif isinstance(b,int):
                aux = []
                for i in range(a._rows):
                    aux2 = []
                    for j in range(a._cols):
                        aux2.append(a._value[i][j]-b)
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
        if isinstance(a,Vector_type):
            if isinstance(b,Matrix_type):
                if a._direction:
                    if a.size() == b._cols:
                        aux = []
                        for i in range(b._rows):
                            aux2 = []
                            for j in range(b._cols):
                                aux2.append(a._value[j]-b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                raise TypeError("Vector and Matrix have different sizes")
                if not a._direction:
                    if a.size() == b._rows:
                        aux = []
                        for i in range(b._rows):
                            aux2 = []
                            for j in range(b._cols):
                                aux2.append(a._value[i]-b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                raise TypeError("Vector and Matrix have different sizes")
            elif isinstance(b,Vector_type):
                if b.size() == a.size():
                    aux = []
                    for i in range(a.size()):
                        aux.append(a._value[i]-b._value[i])
                    D = Matrix(self.subdomain(),self._name)
                    if a._direction == b._direction:
                        if a._direction:
                            return Vector_type(D,aux,'r')
                        else:
                            return Vector_type(D,aux,'c')
                    return Vector_type(D,aux,'r')
                raise TypeError("Vectors size not match")
            elif isinstance(b,Base_type.Base_type):
                aux = []
                for i in range(a.size()):
                    aux.append(a._value[i]-b)
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
            elif isinstance(b,int):
                aux = []
                for i in range(a.size()):
                    aux.append(a._value[i]-b)
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
                    
        if isinstance(b,Matrix_type):
            if isinstance(a,Vector_type):
                if a._direction:
                    if a.size() == b._cols:
                        aux = []
                        for i in range(b._rows):
                            aux2=[]
                            for j in range(b._cols):
                                aux2.append(a._value[j]-b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
                if not a._direction:
                    if a.size() == b._rows:
                        aux = []
                        for i in range(b._rows):
                            aux2=[]
                            for j in range(b._cols):
                                aux2.append(a._value[i]-b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
            elif isinstance(a,Base_type.Base_type):
                aux = []
                for i in range(b._rows):
                    aux2 = []
                    for j in range(b._cols):
                        aux2.append(a-b._value[i][j])
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
            elif isinstance(a,int):
                aux = []
                for i in range(b._rows):
                    aux2 = []
                    for j in range(b._cols):
                        aux2.append(a-b._value[i][j])
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
        if isinstance(b,Vector_type):
            if isinstance(a,Base_type.Base_type):
                aux = []
                for i in range(b.size()):
                    aux.append(a-b._value[i])
                D = Matrix(self.subdomain(),self._name)
                if b._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
            elif isinstance(a,int):
                aux = []
                for i in range(b.size()):
                    aux.append(a-b._value[i])
                D = Matrix(self.subdomain(),self._name)
                if b._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
        else:
            return a - b
    
    def _mul_(self, a, b):
        """operates two elements using this domain operation *"""
        if isinstance(a,float):
            aux = a>>QQ(self.base_domain())
            return aux*b
        if isinstance(b,float):
            aux = b>>QQ(self.base_domain())
            return a*aux
        
        #Treure els elements que de veritat són del subdomini
        if isinstance(a,Matrix_type):
            if a.ncols() == 1 and a.nrows() == 1:
                return self._mul_(a[0][0],b) 
        if isinstance(b,Matrix_type):
            if b.ncols() == 1 and b.nrows() == 1:
                return self._mul_(a,b[0][0]) 
        if isinstance(a,Vector_type):
            if a.size()==1:
                return self._mul_(a[0],b) 
        if isinstance(b,Vector_type):
            if b.size()==1:
                return self._mul_(a,b[0]) 
            
        if isinstance(a,Matrix_type):
            if isinstance(b,Matrix_type):
                if a._cols == b._rows:
                    aux = []
                    zeroval = 0>>self.subdomain()
                    for i in range(a._rows):
                        aux2 = []
                        for j in range(b._cols):
                            sumval = zeroval
                            for k in range(a._cols):
                                sumval = sumval + a._value[i][k]*b._value[k][j]
                            aux2.append(sumval)
                        aux.append(aux2)
                    D = Matrix(self.subdomain(),self._name)
                    return Matrix_type(D,aux)
                raise TypeError("Matrixes sizes does not match")
            elif isinstance(b,Vector_type):
                if b._direction:
                    if a._cols == 1:
                        aux = []
                        for i in range(a._rows):
                            aux2=[]
                            for j in range(b.size()):
                                aux2.append(a._value[i][0]*b._value[j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif b.size() == a._cols:
                        aux = []
                        for i in range(a._rows):
                            aux2=0>>self.subdomain()
                            for j in range(a._cols):
                                aux2 = aux2 + a._value[i][j]*b._value[j]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux,'c')
                    elif b.size() == a._rows:
                        return b*a
                    raise TypeError("Matrix and vector size not match")
                if not b._direction:
                    if b.size() == a._cols:
                        aux = []
                        for i in range(a._rows):
                            aux2=0>>self.subdomain()
                            for j in range(a._cols):
                                aux2 = aux2 + a._value[i][j]*b._value[j]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux,'c')
                    elif b.size() == a._rows:
                        return b*a
                    raise TypeError("Matrix and vector size not match")
            elif isinstance(b,Base_type.Base_type):
                aux = []
                for i in range(a._rows):
                    aux2 = []
                    for j in range(a._cols):
                        aux2.append(a._value[i][j]*b)
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
            elif isinstance(b,int):
                aux = []
                for i in range(a._rows):
                    aux2 = []
                    for j in range(a._cols):
                        aux2.append(a._value[i][j]*b)
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
        if isinstance(a,Vector_type):
            if isinstance(b,Matrix_type):
                if a._direction:
                    if a.size() == b._rows:
                        aux = []
                        for i in range(b._cols):
                            aux2 = 0>>self.subdomain()
                            for j in range(b._rows):
                                aux2 = aux2 + a._value[j]*b._value[j][i]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux)
                    elif a.size() == b._cols:
                        return b*a.transpose()
                    raise TypeError("Vector and Matrix have different sizes")
                if not a._direction:
                    if 1 == b._rows:
                        aux = []
                        for i in range(a.size()):
                            aux2 = []
                            for j in range(b._cols):
                                aux2.append(a._value[i]*b._value[0][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif a.size() == b._rows:
                        return a.transpose()*b
                    raise TypeError("Vector and Matrix have different sizes")
            elif isinstance(b,Vector_type):
                if b.size() == a.size():
                    aux = []
                    for i in range(a.size()):
                        aux.append(a._value[i]*b._value[i])
                    D = Matrix(self.subdomain(),self._name)
                    if a._direction == b._direction:
                        if a._direction:
                            return Vector_type(D,aux,'r')
                        else:
                            return Vector_type(D,aux,'c')
                    return Vector_type(D,aux,'r')
                raise TypeError("Vectors size not match")
            elif isinstance(b,Base_type.Base_type):
                aux = []
                for i in range(a.size()):
                    aux.append(a._value[i]*b)
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
            elif isinstance(b,int):
                aux = []
                for i in range(a.size()):
                    aux.append(a._value[i]*b)
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
                    
        if isinstance(b,Matrix_type):
            if isinstance(a,Vector_type):
                if a._direction:
                    if a.size() == b._cols:
                        aux = []
                        for i in range(b._rows):
                            aux2=[]
                            for j in range(b._cols):
                                aux2.append(a._value[j]*b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
                if not a._direction:
                    if a.size() == b._rows:
                        aux = []
                        for i in range(b._rows):
                            aux2=[]
                            for j in range(b._cols):
                                aux2.append(a._value[i]*b._value[i][j])
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    raise TypeError("Matrix and vector size not match")
            elif isinstance(a,Base_type.Base_type):
                aux = []
                for i in range(b._rows):
                    aux2 = []
                    for j in range(b._cols):
                        aux2.append(a*b._value[i][j])
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
            elif isinstance(a,int):
                aux = []
                for i in range(b._rows):
                    aux2 = []
                    for j in range(b._cols):
                        aux2.append(a*b._value[i][j])
                    aux.append(aux2)
                D = Matrix(self.subdomain(),self._name)
                return Matrix_type(D,aux)
        if isinstance(b,Vector_type):
            if isinstance(a,Base_type.Base_type):
                aux = []
                for i in range(b.size()):
                    aux.append(a*b._value[i])
                D = Matrix(self.subdomain(),self._name)
                if b._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
            elif isinstance(a,int):
                aux = []
                for i in range(b.size()):
                    aux.append(a*b._value[i])
                D = Matrix(self.subdomain(),self._name)
                if b._direction:
                    return Vector_type(D,aux,'r')
                else:
                    return Vector_type(D,aux,'c')
        else:
            return a*b

    def _inv_(self, a):
        if self.subdomain().is_field():
            D = self
        else:
            aux = QQ(self.subdomain())
            D = Matrix(aux)
        if not self.is_invertible(a,D):
            raise ZeroDivisionError("This element is not invertible")
        flagD = 0
        if isinstance(a,Vector_type):
            aux = []
            for i in range(len(a._value)):
                if D == self:
                    aux.append(1/a._value[i])
                else:
                    aux2 = a._value[i]>>D.subdomain()
                    aux2 = 1/aux2
                    if flagD>0:
                        aux.append(aux2)
                    elif aux2._value[1] == 1:
                        aux.append(aux2._value[0])
                    else:
                        aux.append(aux2)
                        flagD = i
            if flagD == 0:
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux)
                else:
                    return Vector_type(D,aux,'c')
#            for i in range(flagD):
#                aux[i] = aux[i]>>D.subdomain()
            if a._direction:
                return Vector_type(D,aux)
            else:
                return Vector_type(D,aux,'c')    
        if isinstance(a,Matrix_type):
            m = a.ncols()
            Dom = Matrix(self.subdomain(),self._name)
            Ia = Matrix_type(Dom,m,m)
            self.gaussian_rows_elimination(a)
            for i in range(m):
                x = a.P*self.e_i(i,m,'c')
                for j in range(1,m):
                    for k in range(j):
                        x[j] = x[j] - a.L[j,k]*x[k]
                        if a.L.domain() == D and flagD == 0:
                            flagD = j
#                if flagD>0:
#                    for i in range(flagD):
#                        x[j] = x[j]>>D.subdomain()
                if D == self:
                    x[m-1] = x[m-1]/a.U[m-1,m-1]
                else:
                    aux2 = a.U[m-1,m-1]>>D.subdomain()
                    aux2 = 1/aux2
                    if flagD>0 or flagD <0:
                        x[m-1] = x[m-1]*aux2
                    elif aux2._value[1] == 1:
                        x[m-1] = x[m-1]*aux2._value[0]
                    else:
                        x[m-1] = x[m-1]*aux2
                        flagD = m-1
                for j in range(m-2,-1,-1):
                    for k in range(j+1,m):
                        x[j] = x[j] - a.U[j,k]*x[k]
                    if D == self:
                        x[j] = x[j]/a.U[j,j]
                    else:
                        aux2 = a.U[j,j]>>D.subdomain()
                        aux2 = 1/aux2
                        if flagD>0 or flagD<0:
                           x[j] = x[j]*aux2
                        elif aux2._value[1] == 1:
                            x[j] = x[j]*aux2._value[0]
                        else:
                            x[j] = x[j]*aux2
                            flagD = j
                Ia[:,i] = x
            return Ia
        else:
            return D.subdomain()._inv_(a)
    
    def _rd(self,n,m): # random matrix of size nxm
        A = self.subdomain()
        q = A.cardinal()
        aux = [[]]
        for h in range(n):
            C = []
            for k in range(m):
                x = A._rd()
                C += [x]
            aux.append(C)
        D = Matrix(self.subdomain(),self._name)
        return Matrix_type(D,aux)
       
    def change_characteristic(self,n):    
        if self.characteristic() == n:
            return self
        if (n == 0):
            return self.transform_to_ZZ()
        cself = self.characteristic()
        if ((cself%n) == 0):
            SDom = self.subdomain().change_characteristic(n)
            return Matrix(SDom)
        raise TypeError("Wrong characteristic")        
    
    def transform_to_ZZ(self):
        if self.characteristic() == 0:
            return self
        SDom = self.subdomain().transform_to_ZZ()
        return Matrix(SDom)
    
    def K_(self):
        return self.subdomain().K_()
    
    def _ext_vector_(self,a,b): 
        if isinstance(a,float):
            DD = QQ(self.base_domain())
            aux = a>>DD
            if self.is_sub_domain(DD):
                return self._ext_vector_(aux,b)
            DD = QQ(self.subdomain())
            if DD == self.subdomain():
                aux = (aux._value[0]>>self.subdomain())/(aux._value[1]>>self.subdomain())
                return self._ext_vector_(aux,b)
            else:
                D = Matrix(DD)
                if isinstance(b,Vector_type):
                    if b._direction:
                        aux2 = D.Vector_type(b._value)
                    else:
                        aux2 = D.Vector_type(b._value,'c')
                if isinstance(b,Matrix_type):
                    aux2 = D.element(b._value)
                return D._ext_vector_(aux,aux2)
        if isinstance(b,float):
            DD = QQ(self.base_domain())
            aux = b>>DD
            if self.is_sub_domain(DD):
                return self._ext_vector_(a,aux)
            DD = QQ(self.subdomain())
            if DD == self.subdomain():
                aux = (aux._value[0]>>self.subdomain())/(aux._value[1]>>self.subdomain())
                return self._ext_vector_(a, aux)
            else:
                D = Matrix(DD)
                if isinstance(a,Vector_type):
                    if a._direction:
                        aux2 = D.Vector_type(a._value)
                    else:
                        aux2 = D.Vector_type(a._value,'c')
                if isinstance(a,Matrix_type):
                    aux2 = D.element(a._value)
                return D._ext_vector_(aux2,aux)
        if isinstance(a,Vector_type):
            if isinstance(b,list):
                if isinstance(b[0],list):
                    if len(b[0])!=1:
                        raise TypeError("list not valid to extend the vector")
                    if isinstance(b[0][0],int):
                        D = ZZ()
                    elif isinstance(b[0][0],float):
                        D = QQ()
                    else:
                        D = b[0][0].domain()
                    aux = a._value
                    for i in range(len(b)):
                        if len(b[0])!= len(b[i]):
                            raise TypeError("list not valid to extend the vector")
                        aux.append(b[i][0])
                        if isinstance(b[i][0],int):
                            continue
                        if isinstance(b[i][0],float):
                            D = D.union_domain(QQ())[0]
                        else:
                            D = D.union_domain(b[i][0].domain())[0]
                    if D == self:
                        D = Matrix(self._subdomain(),self._name)
                    else:
                        D = Matrix(D)
                    if a._direction:
                        return Vector_type(D,aux)
                    else:
                        return Vector_type(D,aux,'c')
                else:
                    aux = a._value
                    if isinstance(b[0],int):
                        D = ZZ()
                    elif isinstance(b[0],float):
                        D = QQ()
                    else:
                        D = b[0].domain()
                    aux = a._value
                    for i in range(len(b)):
                        aux.append(b[i])
                        if isinstance(b[i],int):
                            continue
                        if isinstance(b[0],float):
                            D = D.union_domain(QQ())[0]
                        else:
                            D = D.union_domain(b[i].domain())[0]
                    D = Matrix(self.subdomain(),self._name)
                    if a._direction:
                        return Vector_type(D,aux)
                    else:
                        return Vector_type(D,aux,'c')
            if isinstance(b,Vector_type):
                aux = a._value + b._value
                if a.domain != b.domain():
                    Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                else:
                    Dun = a.domain().subdomain()
                D = Matrix(Dun[0])   
                if a._direction:
                    return Vector_type(D,aux)
                else:
                    return Vector_type(D,aux,'c')
            if isinstance(b,Matrix_type):
                if b._rows == 1:
                    aux = a._value + b._value[0]
                    if a.domain != b.domain():
                        Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                    else:
                        Dun = a.domain().subdomain()
                    D = Matrix(Dun[0])
                    if a._direction:
                        return Vector_type(D,aux)
                    else:
                        return Vector_type(D,aux,'c')
                elif b._cols == 1:
                    aux = a._value
                    for i in range(b._rows):
                        aux.append(b[0][i])
                    if a.domain != b.domain():
                        Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                    else:
                        Dun = a.domain().subdomain()
                    D = Matrix(Dun[0])
                    if a._direction:
                        return Vector_type(D,aux)
                    else:
                        return Vector_type(D,aux,'c')
                else:
                    raise TypeError("Sizes not valid")
            if isinstance(b,int):
                aux = a._value + [b>>self.subdomain()]
                D = Matrix(self.subdomain(),self._name)
                if a._direction:
                    return Vector_type(D,aux)
                else:
                    return Vector_type(D,aux,'c')
            if isinstance(b,Base_type.Base_type):
                aux = a._value + [b>>self.subdomain()]
                if a.domain != b.domain():
                    Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                else:
                    Dun = a.domain().subdomain()
                D = Matrix(Dun[0])
                if a._direction:
                    return Vector_type(D,aux)
                else:
                    return Vector_type(D,aux,'c')
            else:
                raise TypeError("Wrong Parameter")
        else:
            if isinstance(b,Vector_type):
                D = Matrix(self.subdomain(),self._name)
                if b._direction:
                    A = Vector_type(D,a)
                else:
                    A = Vector_type(D,a,'c')
                return self._ext_vector_(A,b)
            else:
                D = Matrix(self.subdomain(),self._name)
                A = Vector_type(D,a)
            return self._ext_vector_(A,b)   
            
    def _ext_with_columns_(self,a,b):
        if isinstance(a,float):
            DD = QQ(self.base_domain())
            aux = a>>DD
            if self.is_sub_domain(DD):
                return self._ext_with_columns_(aux,b)
            DD = QQ(self.subdomain())
            if DD == self.subdomain():
                aux = (aux._value[0]>>self.subdomain())/(aux._value[1]>>self.subdomain())
                return self._ext_with_columns_(aux,b)
            else:
                D = Matrix(DD)
                if isinstance(b,Vector_type):
                    if b._direction:
                        aux2 = D.Vector_type(b._value)
                    else:
                        aux2 = D.Vector_type(b._value,'c')
                if isinstance(b,Matrix_type):
                    aux2 = D.element(b._value)
                return D._ext_with_columns_(aux,aux2)
        if isinstance(b,float):
            DD = QQ(self.base_domain())
            aux = b>>DD
            if self.is_sub_domain(DD):
                return self._ext_with_columns_(a,aux)
            DD = QQ(self.subdomain())
            if DD == self.subdomain():
                aux = (aux._value[0]>>self.subdomain())/(aux._value[1]>>self.subdomain())
                return self._ext_with_columns_(a, aux)
            else:
                D = Matrix(DD)
                if isinstance(a,Vector_type):
                    if a._direction:
                        aux2 = D.Vector_type(a._value)
                    else:
                        aux2 = D.Vector_type(a._value,'c')
                if isinstance(a,Matrix_type):
                    aux2 = D.element(a._value)
                return D._ext_with_columns_(aux2,aux)
        if isinstance(a,Vector_type):
            if isinstance(b,list):
                if isinstance(b[0],list):
                    if a._direction:
                        if len(b)>1:
                            raise TypeError("Wrong sizes")
                        aux = a._value + b[0]
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux)
                    else:
                        if len(b)!=a.size():
                            raise TypeError("Wrong list size")
                        aux = []
                        for i in range(len(b)):
                            if len(b[i])>len(b[0]):
                                raise TypeError("Wrong list sizes")
                            aux2 = [a._value[i]] + b[i]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        if len(aux)>1:
                            return Matrix_type(D,aux)
                        else:
                            return Vector_type(D,aux)
                else:
                    D = Matrix(self.subdomain(),self._name)
                    if a._direction:
                        aux = a._value + b
                        return Vector_type(D,aux)
                    else:
                        if a.size()==1:
                            aux = a._value + b
                            return Vector_type(D,aux)
                        elif a.size() == len(b):
                            aux = []
                            for i in range(len(b)):
                                aux2 = [a._value[i]] + [b[i]]
                                aux.append(aux2)
                            return Matrix_type(D,aux)
                        else:
                            raise TypeError("Wrong list size")
            if isinstance(b,Vector_type):
                if a._direction:
                    if b._direction:
                        aux = a._value + b._value
                        if a.domain() == b.domain():
                            D = Matrix(self.subdomain(),self._name)
                        else:
                            Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                            D = Matrix(Dun[0])
                        return Vector_type(D,aux)
                    else:
                        if len(b)==1:
                            aux = a._value + b._value
                            if a.domain() == b.domain():
                                D = Matrix(self.subdomain(),self._name)
                            else:
                                Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                                D = Matrix(Dun[0])
                            return Vector_type(D,aux)
                        else:
                            raise TypeError("Wrong sizes")
                else:
                    if b._direction:
                        if a.size() == 1:
                            aux = a._value + b._value
                            if a.domain() == b.domain():
                                D = Matrix(self.subdomain(),self._name)
                            else:
                                Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                                D = Matrix(Dun[0])
                            return Vector_type(D,aux)
                        else:
                            raise TypeError("Wrong sizes")
                    else:
                        if b.size() == a.size():
                            aux = []
                            for i in range(len(a)):
                                aux2 = [a._value[i]] + [b._value[i]]
                                aux.append(aux2)
                            if a.domain() == b.domain():
                                D = Matrix(self.subdomain(),self._name)
                            else:
                                Dun = a.domain().subdomain().union_domain(b.domain().subdomain())
                                D = Matrix(Dun[0])
                            return Matrix_type(D,aux)
                        else:
                            raise TypeError("Vectors' sizes does not match")
            if isinstance(b,Matrix_type):
                if a._direction:
                    if b._rows == 1:
                        aux = a._value + b._value[0]
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D, aux)
                    elif a.size() == b._rows:
                        aux = []
                        for i in range(len(a)):
                            aux2 = [a._value[i]] + b._value[i]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")
                else:
                    if a.size() == b._rows:
                        aux = []
                        for i in range(len(a)):
                            aux2 = [a._value[i]] + b._value[i]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif b._rows == 1:
                        aux = a._value + b._value[0]
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D, aux)
                    else:
                        raise TypeError("Wrong sizes")
            if isinstance(b,int):
                if a._direction:
                    aux = a._value + [b>>self.subdomain()]
                    D = Matrix(self.subdomain(),self._name)
                    return Vector_type(D,aux)
                else:
                    if a.size() == 1:
                        aux = a._value + [b>>self.subdomain()]
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")
            if isinstance(b,Base_type.Base_type):
                if a._direction:
                    aux = a._value + [b>>self.subdomain()]
                    D = Matrix(self.subdomain(),self._name)
                    return Vector_type(D,aux)
                else:
                    if a.size() == 1:
                        aux = a._value + [b>>self.subdomain()]
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")                     
        if isinstance(a,Matrix_type):
            if isinstance(b,list):
                if isinstance(b[0],list):
                    if a._rows == len(b):
                        aux = []
                        for i in range(a._rows):
                            aux2 = a._value[i] + b[i]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                else:
                    if a._rows == 1:
                        aux = [a._value[0] + b]
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif a._rows == len(b):
                        aux = []
                        for i in range(a._rows):
                            aux2 = a._value[i] + [b[i]]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")
            if isinstance(b,Vector_type):
                if b._direction:
                    if a._rows == 1:
                        aux = [a._value[0] + b._value]
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif a._rows == b.size():
                        aux = []
                        for i in range(a._rows):
                            aux2 = a._value[i] + [b._value[i]]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")
                else:
                    if b.size() == a._rows:
                        aux = []
                        for i in range(a._rows):
                            aux2 = a._value[i] + [b._value[i]]
                            aux.append(aux2)
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif a._rows == 1:
                        aux = [a._value[0] + b._value]
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Vectors' sizes does not match")
            if isinstance(b,Matrix_type):
                if a._rows == b._rows:
                    aux = []
                    for i in range(a._rows):
                        aux2 = a._value[i] + b._value[i]
                        aux.append(aux2)
                    D = Matrix(self.subdomain(),self._name)
                    return Matrix_type(D, aux)
                else:
                    raise TypeError("Wrong sizes")
            if isinstance(b,int):
                if a._rows == 1:
                    aux = [a._value[0] + [b>>self.subdomain()]]
                    D = Matrix(self.subdomain(),self._name)
                    return Matrix_type(D,aux)
                else:
                    raise TypeError("Wrong sizes")
            if isinstance(b,Base_type.Base_type):
                if a._rows == 1:
                    aux = [a._value[0] + [b>>self.subdomain()]]
                    D = Matrix(self.subdomain(),self._name)
                    return Matrix_type(D,aux)
                else:
                    raise TypeError("Wrong sizes")
        else:
            D = Matrix(self.subdomain(),self._name)
            A = Matrix_type(D,a)
            return self._ext_with_columns_(A,b)
                        
    def _ext_with_rows_(self,a,b):
        if isinstance(a,float):
            DD = QQ(self.base_domain())
            aux = a>>DD
            if self.is_sub_domain(DD):
                return self._ext_with_rows_(aux,b)
            DD = QQ(self.subdomain())
            if DD == self.subdomain():
                aux = (aux._value[0]>>self.subdomain())/(aux._value[1]>>self.subdomain())
                return self._ext_with_rows_(aux,b)
            else:
                D = Matrix(DD)
                if isinstance(b,Vector_type):
                    if b._direction:
                        aux2 = D.Vector_type(b._value)
                    else:
                        aux2 = D.Vector_type(b._value,'c')
                if isinstance(b,Matrix_type):
                    aux2 = D.element(b._value)
                return D._ext_with_rows_(aux,aux2)
        if isinstance(b,float):
            DD = QQ(self.base_domain())
            aux = b>>DD
            if self.is_sub_domain(DD):
                return self._ext_with_rows_(a,aux)
            DD = QQ(self.subdomain())
            if DD == self.subdomain():
                aux = (aux._value[0]>>self.subdomain())/(aux._value[1]>>self.subdomain())
                return self._ext_with_rows_(a, aux)
            else:
                D = Matrix(DD)
                if isinstance(a,Vector_type):
                    if a._direction:
                        aux2 = D.Vector_type(a._value)
                    else:
                        aux2 = D.Vector_type(a._value,'c')
                if isinstance(a,Matrix_type):
                    aux2 = D.element(a._value)
                return D._ext_with_rows_(aux2,aux)
        if isinstance(a,Vector_type):
            if isinstance(b,list):
                if isinstance(b[0],list):
                    if a._direction:
                        if len(b[0])!= a.size():
                            raise TypeError("Wrong sizes")
                        aux = [a._value]
                        for i in range(len(b)):
                            if len(b[0])!=len(b[i]):
                                raise TypeError("Wrong list sizes")
                            aux.append(b[i])
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        if len(b[0])!=1:
                            raise TypeError("Wrong list size")
                        aux = a._value
                        for i in range(len(b)):
                            if len(b[i])>len(b[0]):
                                raise TypeError("Wrong list sizes")
                                aux.append(b[i])
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux,'c')
                else:
                    if a._direction:
                        if a.size() == len(b):
                            aux = [a._value]
                            aux.append(b)
                            D = Matrix(self.subdomain(),self._name)
                            return Matrix_type(D, aux)
                    else:
                        aux = a._value + b
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux,'c')
            if isinstance(b,Vector_type):
                if a._direction:
                    if b._direction:
                        if b.size() == a.size():
                            aux = [a._value]
                            aux.append(b._value)
                            D = Matrix(self.subdomain(),self._name)
                            return Matrix_type(D,aux)
                        elif a.size() == 1:
                            aux = a._value + b._value
                            D = Matrix(self.subdomain(),self._name)
                            return Vector_type(D,aux,'c')
                        else:
                            raise TypeError("Vectors' sizes does not match")
                    else:
                        if a.size() == 1:
                            aux = a._value + b._value
                            D = Matrix(self.subdomain(),self._name)
                            return Vector_type(D,aux,'c')
                        elif a.size() == b.size():
                            aux = [a._value]
                            aux.append(b._value)
                            D = Matrix(self.subdomain(),self._name)
                            return Matrix_type(D,aux)
                        else:
                            raise TypeError("Vectors' sizes does not match")
                else:
                    aux = a._value + b._value
                    D = Matrix(self.subdomain(),self._name)
                    return Vector_type(D,aux,'c')
            if isinstance(b,Matrix_type):
                if not a._direction:
                    if b._cols == 1:
                        aux = a._value
                        for i in range(a._rows):
                            aux.append(b._value[i][0])
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D, aux,'c')
                    elif a.size() == b._cols:
                        aux = [a._value]+b._value
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")
                else:
                    if a.size() == b._cols:
                        aux = [a._value]+b._value
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif b._cols == 1:
                        aux = a._value
                        for i in range(a._rows):
                            aux.append(b._value[i][0])
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D, aux,'c')
                    else:
                        raise TypeError("Wrong sizes")
            if isinstance(b,int):
                if not a._direction:
                    aux = a._value + [b>>self.subdomain()]
                    D = Matrix(self.subdomain(),self._name)
                    return Vector_type(D,aux,'c')
                else:
                    if a.size() == 1:
                        aux = a._value + [b>>self.subdomain()]
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux,'c')
                    else:
                        raise TypeError("Wrong sizes")
            if isinstance(b,Base_type.Base_type):
                if not a._direction:
                    aux = a._value + [b>>self.subdomain()]
                    D = Matrix(self.subdomain(),self._name)
                    return Vector_type(D,aux,'c')
                else:
                    if a.size() == 1:
                        aux = a._value + [b>>self.subdomain()]
                        D = Matrix(self.subdomain(),self._name)
                        return Vector_type(D,aux,'c')
                    else:
                        raise TypeError("Wrong sizes")                     
        if isinstance(a,Matrix_type):
            if isinstance(b,list):
                if isinstance(b[0],list):
                    if a._cols == len(b[0]):
                        aux = a._value
                        for i in range(len(b)):
                            if len(b[0])!=len(b[i]):
                                raise TypeError("Wrong list size")
                            aux.append(b[i])
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                else:
                    if a._cols == len(b):
                        aux = a._value + [b]
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif a._cols == 1:
                        aux = a._value
                        for i in range(len(b)):
                            aux.append([b[i]])
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")
            if isinstance(b,Vector_type):
                if b._direction:
                    if a._cols == b.size():
                        aux = a._value + [b._value]
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif a._cols == 1:
                        aux = a._value
                        for i in range(b.size()):
                            aux.append([b._value[i]])
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Wrong sizes")
                else:
                    if a._cols == 1:
                        aux = a._value
                        for i in range(b.size()):
                            aux.append([b._value[i]])
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    elif b.size() == a._cols:
                        aux = a._value + [b._value]
                        D = Matrix(self.subdomain(),self._name)
                        return Matrix_type(D,aux)
                    else:
                        raise TypeError("Vectors' sizes does not match")
            if isinstance(b,Matrix_type):
                if a._cols == b._cols:
                    aux = a._value + b._value
                    D = Matrix(self.subdomain(),self._name)
                    return Matrix_type(D, aux)
                else:
                    raise TypeError("Wrong Matrix sizes")
            if isinstance(b,int):
                if a._cols == 1:
                    aux = a._value + [[b>>self.subdomain()]]
                    D = Matrix(self.subdomain(),self._name)
                    return Matrix_type(D,aux)
                else:
                    raise TypeError("Wrong sizes")
            if isinstance(b,Base_type.Base_type):
                if a._cols == 1:
                    aux = a._value + [[b>>self.subdomain()]]
                    D = Matrix(self.subdomain(),self._name)
                    return Matrix_type(D,aux)
                else:
                    raise TypeError("Wrong sizes")                          
        else:
            D = Matrix(self.subdomain(),self._name)
            A = Matrix_type(D,a)
            return self._ext_with_rows_(A,b)

    def gaussian_rows_elimination(self,A):
        if not isinstance(A,Matrix_type):
            raise TypeError("A must be a matrix")
        if A.L != None:
            return [A.L, A.U, A.P]
        else:
            if self.subdomain().is_field():
                D = self
            else:
                aux = QQ(self.subdomain())
                D = Matrix(aux)
            U = A.cp()
            L = self.eye(A.nrows())
            P = self.eye(A.nrows())
            if (A.nrows()>A.ncols()):
                n = A.ncols()
            else:
                n = A.nrows()
            flagD = False
            for k in range(n):
                s = U[k,k]
                i = k;
                restartF = False
                while True:
                    if isinstance(s,int):
                        if s != 0:
                            break
                    else:
                        if s!=0:
                            break
                    i = i + 1
                    if i == A.nrows():
                        restartF = True
                        break
                    s = U[i,k]
                if restartF:
                    continue
                if i != k:
                    A.signPLU = - A.signPLU
                    aux = U.row(i)
                    U[i] = U.row(k)
                    U[k] = aux
                    aux = L[k,0:k]
                    L[k,0:k] = L[i,0:k]
                    L[i,0:k] = aux
                    aux = P.row(i)
                    P[i] = P.row(k)
                    P[k] = aux
                for j in range(k+1,A.nrows()):
                    if D == self:
                        L[j,k] = U[j,k]/U[k,k]
                    else:
                        aux2 = U[k,k]>>D.subdomain()
                        aux2 = 1/aux2
                        if aux2.domain()._fraction_structure:
                            if aux2._value[1] == 1:
                                L[j,k] = U[j,k]*aux2._value[0]
                            else:
                                L[j,k] = U[j,k]*aux2
                        else:
                            L[j,k] = U[j,k]*aux2
                    U[j,:] = U[j,:] - L[j,k]*U[k,:]
            A.L = L
            A.U = U
            A.P = P
            return [A.L, A.U, A.P]
    
    def _det_(self,A):
        if A.ncols() != A.nrows():
            raise TypeError("Matrix must be square")
        if self.subdomain() == ZZ():
            B = Matrix_type(Matrix(QQ(ZZ())),[A._value[i][:] for i in range(len(A._value))])
            s = B.det()
            if s._value[1] == 1:
                return s._value[0]
            return s
        if not self.subdomain().is_field():
            if A.L == None:
                self.pivot_matrix(A)
                divD = 1
                n = A.ncols()
                for i in range(n-2):
                    if A.Piv[i,i]==0:
                        return A.Piv[i,i]
                    divD = divD*(A.Piv[i,i]**(n-2-i))
                if A.Piv[n-2,n-2]==0:
                    return A.Piv[n-2,n-2]
                return A.signPiv*A.Piv[n-1,n-1]//divD
            else:
                val = 1>>self.subdomain()
                for i in range(A.ncols()):
                    val = val * A.U[i,i]
                if val.domain() != self.subdomain():
                    val = val._value[0]
                return A.signPLU*val
        else:
            self.gaussian_rows_elimination(A)
            val = 1>>self.subdomain()
            for i in range(A.ncols()):
                val = val * A.U[i,i]
            if val.domain() != self.subdomain():
                val = val._value[0]
            return A.signPLU*val
    
    def _trace_(self,A):
        if A.nrows() > A.ncols():
            n = A.ncols()
        else:
            n = A.nrows()
        val = 0>>self.subdomain()
        for i in range(n):
            val = val + A[i,i]
        return val
    
    def _rank_(self,A):
        if A.nrows() > A.ncols():
            A2 = A.cp()
            n = A.ncols()
        else:
            A2 = A.transpose()
            n = A.nrows()
        self.gaussian_rows_elimination(A2)
        rang = 0
        for i in range(n):
            if A2.U[i,i]!= 0:
                rang = rang + 1
        return rang
    
    def rank_block(self,A):
        S = self.rank_block_ext(A)
        return S[0]
        
    def rank_block_ext(self,A):
        if A.nrows() > A.ncols():
            A2 = A.cp()
            n = A.ncols()
            ftrans = False
        else:
            A2 = A.transpose()
            n = A.nrows()
            ftrans = True
        self.gaussian_rows_elimination(A2)
        R = []
        for i in range(n):
            if A2.U[i,i]!= 0:
                R.append(i)
        S = []
        PP = A2.P.transpose()
        for i in R:
            S.append(PP[:,i].to_list().index(1>>self.subdomain()))
        A3 = A2[S,R]
        if ftrans:
            return (A3.transpose(),R,S)
        else:
            return (A3,S,R)
                
    def prune(self,A):
        if A.nrows() < A.ncols():
            A2 = A.cp()
            n = A.nrows()#A.ncols()
            ftrans = False
        else:
            A2 = A.transpose()
            n = A.ncols()#A.nrows()
            ftrans = True
        self.gaussian_rows_elimination(A2)
        R = []
        for i in range(n):
            if A2.U[i,i]!= 0:
                R.append(i)
        S = []
        PP = A2.P.transpose()
        for i in R:
            S.append(PP[:,i].to_list().index(1>>self.subdomain()))
        if ftrans:
             A3 = A2[:,S]
             return A3.transpose()
        else:
            A3 = A2[S,:]
            return A3

    def blow(self,A,K):
        if not self.subdomain().is_sub_domain(K):
            raise TypeError("K must be a subdomain")
        if K == self.subdomain():
            return A.cp()
        s = self.subdomain().blow(A[0,0],K)
        ns = len(s)
        DomM = Matrix(K)
        R = Matrix_type(DomM,ns*A.nrows(),A.ncols())
        for i in range(A.nrows()):
            for j in range(A.ncols()):
                s = self.subdomain().blow(A[i,j],K)
                R[slice(ns*i,ns*(i+1)),j]=s
        return R
    
    def kernel(self,A):
        if A.domain().subdomain() == ZZ():
            return self.kernelZZ(A)
        B = self._ext_with_rows_(A,self.eye(A.ncols()))
        B = B.transpose()
        B.lu()
        S = B.P.transpose()*B.U
        S = S.transpose()
        T = []
        for i in range(S.ncols()):
            nflag = True
            for j in range(A.nrows()):
                if S[j,i]!= 0:
                    nflag = False
                    break
            if nflag:
                T += [i]
        if len(T) == 0:
            return self.element([0>>self.subdomain()]*A.nrows()).transpose()
        return S[A.nrows():,T]
        
    def kernelZZ(self,A):
        nr = A.nrows()
        U = self._ext_with_rows_(A,self.eye(A.ncols()))
        n = U.ncols()
        nrz = 0
        ncolnzero = 0
        for k in range(nr):
            i = k-nrz
            s = U[k,i]
            restartF = False
            while True:
                if s != 0:
                    ncolnzero = ncolnzero + 1
                    break
                i = i + 1
                if i == n:
                    restartF = True
                    break
                s = U[k,i]
            if restartF:
                nrz = nrz + 1
                continue
            if i != k-nrz:
                aux = U.col(i)
                U[:,i] = U.col(k-nrz)
                U[:,k-nrz] = aux
            for j in range(k+1-nrz,n):
                if U[k,j]!= 0:
                    auxgcd = ilcm(U[k,j],U[k,k-nrz])
                    U[:,j] = U[:,j]*(auxgcd//U[k,j])-U[:,k-nrz]*(auxgcd//U[k,k-nrz])
        if ncolnzero>= U.ncols():
            return A[:,0]-A[:,0]
        return U[nr:,ncolnzero:]
    
    def solve_pivot_matrix_system(self,A,z):
        if isinstance(z,Matrix_type):
            z = z.to_vector()
        if not z._direction:
            b = z
        else:
            b = z.transpose()
        if A.nrows() != len(b):
            raise TypeError("A and z must have the same number of rows")
        
        x = b
        x[len(b)-1] /= A[len(b) - 1, len(b) - 1]        
        for i in range(len(b) - 2, -1, -1):
            for j in range(i+1,A.ncols()):
                b[i] -= A[i,j]*b[j]
            b[i]= b[i]/A[i,i]
        return b
    
    def solve_linear_system(self,A,z):
        if isinstance(z,Matrix_type):
            z = z.to_vector()
        if not isinstance(z,Vector_type):
            if isinstance(z,Base_type.Base_type):
                z = Vector_type(self.cp(),[z],'c')
            elif isinstance(z,int):
                z = Vector_type(self.cp(),[z],'c')
        if not z._direction:
            b = z
        else:
            b = z.transpose()
        B = self._ext_with_columns_(A,b)
        rA = A.rank()
        rB = B.rank()
        if rA != rB:
            return []
        if rA == A.ncols():
            [B,fA,cA] = self.rank_block_ext(A)
            B.lu()
            x = (1/B.L)*B.P*b[fA]
            if type(x) != Vector_type:
                x = Vector_type(self.cp(),[x],'c')
            #x = y.to_list()
            if x.domain().subdomain().is_field():
                D = x.domain()
            else:
                aux = QQ(x.domain().subdomain())
                D = Matrix(aux)
            
            if D == b.domain():
                x[len(x)-1] = x[len(x)-1]/B.U[len(x)-1,len(x)-1]
            else:
                aux2 = B.U[len(x)-1,len(x)-1]>>D.subdomain()
                aux2 = 1/aux2
                if aux2._value[1] == 1:
                    x[len(x)-1] = x[len(x)-1]*aux2._value[0]
                else:
                    x[len(x)-1] = x[len(x)-1]*aux2
            flagD = 0
            for i in range(len(x)-2,-1,-1):
                for j in range(i+1,len(x)):
                    if D == b.domain():
                        x[i] = x[i] - B.U[i,j]*x[j]
                    else:
                        if flagD>0:
                            x[i] = x[i] - B.U[i,j]*x[j]
                        elif D == B.U.domain():
                            if B.U[i,j]._value[1] == 1:
                                x[i] = x[i] - B.U[i,j]._value[0]*x[j]
                            else:
                                x[i] = x[i] - B.U[i,j]*x[j]
                                flagD = i
                        else:
                            x[i] = x[i] - B.U[i,j]*x[j]
                if D == b.domain():
                    if B.U.domain() == D:
                        x[i]=x[i]/B.U[i,i]
                    else:
                        aux2 = B.U[i,i]>>D.subdomain()
                        aux2 = 1/aux2
                        x[i]=x[i]*aux2
                else:
                    if B.U.domain() == D:
                        if flagD>0:
                            x[i]=x[i]/B.U[i,i]
                        elif B.U[i,i]._value[0] == 1:
                            x[i]=x[i]*B.U[i,i]._value[1]
                        else:
                            flagD = i
                            x[i]=x[i]/B.U[i,i]
                    else:
                        aux2 = B.U[i,i]>>D.subdomain()
                        aux2 = 1/aux2
                        if flagD>0:
                            x[i]=x[i]*aux2
                        elif aux2._value[1] == 1:
                            x[i]=x[i]*aux2._value[0]
                        else:
                            flagD = i
                            x[i]=x[i]*aux2
            if flagD>0:
                 for i in range(flagD):
                     x[i] = x[i]>>D.subdomain()
            return x.transform_to_min_domain()>>self
#            if self.characteristic()!=0:
#                D = Matrix(self.subdomain(),self._name)
#                xx = Vector_type(D,x,'c')
#                return xx
#            else:
#                flagInt = True
#                for i in range(len(x)):
#                    if isinstance(x[i],float):
#                        if int(x[i]) == x[i]:
#                            x[i] = int(x[i])
#                        else:
#                            flagInt = False
#                            break
#                if flagInt:
#                    D = Matrix(self.subdomain(),self._name)
#                    xx = Vector_type(D,x,'c')
#                    return xx
#            return x
        else:
            H = self.kernel(A)
            (B,rows,cols) = self.rank_block_ext(A)
            x = self.solve_linear_system(B,b[rows])
            D = x.domain()._create_Matrix_(x.domain().subdomain())
            xout = Vector_type(D,A.ncols(),'c')#[0>>self.subdomain()]*A.ncols()
            for i in range(len(cols)):
                xout[cols[i]] = x[i]
            return (xout,H)
#            print(xout)
#            if self.characteristic()!=0:
#                D = Matrix(self.subdomain(),self._name)
#                xx = Vector_type(D,xout,'c')
#                return (xx,H)
#            else: 
#                flagInt = True
#                for i in range(len(xout)):
#                    if isinstance(xout[i],float):
#                        if int(xout[i]) == xout[i]:
#                            xout[i] = int(xout[i])
#                        else:
#                            flagInt = False
#                            break
#                if flagInt:
#                    D = Matrix(self.subdomain(),self._name)
#                    xx = Vector_type(D,xout,'c')
#                    return (xx,H)
#            return (xout,H)

    def pivot_matrix(self,A):
        if not isinstance(A,Matrix_type):
            raise TypeError("A must be a matrix")
        if A.Piv != None:
            return [A.Piv, A.PivP]
        else:
            U = A.cp()
            P = self.eye(A.nrows())
            if (A.nrows()>A.ncols()):
                n = A.ncols()
            else:
                n = A.nrows()
            for k in range(n):
                s = U[k,k]
                i = k;
                restartF = False
                while True:
                    if isinstance(s,int):
                        if s != 0:
                            break
                    else:
                        if s != 0:#s.is_invertible():
                            break
                    i = i + 1
                    if i == A.nrows():
                        restartF = True
                        break
                    s = U[i,k]
                if restartF:
                    continue
                if i != k:
                    A.signPiv = - A.signPiv
                    aux = U.row(i)
                    U[i] = U.row(k)
                    U[k] = aux
                    aux = P.row(i)
                    P[i] = P.row(k)
                    P[k] = aux
#                for j in range(k-1,-1,-1):
#                    if self.characteristic() == 0:
#                        if U[k,k] == 1 or U[k,k] == -1:
#                            U[j] = U[j] - (U[j,k]//U[k,k])*U[k]
#                        else:
#                            U[j] = U[j]*U[k,k] - U[j,k]*U[k]
#                    else:
#                        if U[k,k].is_invertible():
#                            U[j] = U[j] - (U[j,k]/U[k,k])*U[k]
#                        else:
#                            U[j] = U[j]*U[k,k] - U[j,k]*U[k]
                for j in range(k+1,A.nrows()):
                    valU = U[j,k]
                    for ll in range(k,A.ncols()):
                        U[j,ll] = U[j,ll]*U[k,k] - U[k,ll]*valU
#                    if self.characteristic() == 0:
#                        U[j,k:] = U[j,k:]*U[k,k] - U[k,k:]*U[j,k]
#                    else:
#                        U[j,k:] = U[j,k:]*U[k,k] - U[k,k:]*U[j,k]
            A.Piv = U
            A.PivP = P
            return [A.Piv, A.PivP]

  
                        
                        
                        
                        
     