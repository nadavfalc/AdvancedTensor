# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 11:08:33 2017

@author: admin_local
"""

import PyM.PyWIT.PyECC.Base_type as Base_type
import PyM.PyWIT.PyECC.PyECC_globals as PGl

import copy
    

class Vector_type(Base_type.Base_type):
    """Vector Element class. """
        
    def __init__(self, domain, rows = '', cols = ''):
        '''With 1 parameter:
               - Creates a 0 element vector.
           With 2 parameters:
               - Domain, n: it creates a zero row-vector of size elements
               - Domain, list: it inizelize a row vector with the element list.
           With 3 parameters:
               - Domain, n,1 or (1,n) it creates a zero row(column)- vector of size n
               - Domain, list of elements and r/c initelize a row/column vector with the list.        
        '''
        self._structure = domain
        if rows == '':
            if isinstance(domain,Vector_type):
                self._structure = domain.domain()
                self._direction = domain._direction
                self._value = domain._value[:]
            elif isinstance(domain, Matrix_type):
                if domain._rows == 1:
                    self._structure = domain.domain()
                    self._value = []
                    for i in range(domain._cols):
                        self._value.append(domain._value[0][i])
                    self._direction = True
                elif domain._cols == 1:
                    self._structure = domain.domain()
                    self._value = []
                    for i in range(domain._rows):
                        self._value.append(domain._value[i][0])
                    self._direction = False
                else:
                    raise TypeError("This matrix cannot be converted to a vector")
            else:
                self._value = []#[0>>domain.subdomain()]
                self._direction = True
        elif cols == '':
            if isinstance(rows,list):
                self._value = rows[:]
                self._direction = True
                for i in range(len(self._value)):
                    if self._structure.subdomain().is_sub_element(self._value[i]):
                        self._value[i] = self._value[i]>> self.domain().subdomain()
                    else:
                        self._value[i] = self._value[i]._push_to_union(self.domain().subdomain())
            elif isinstance(rows,int):
                self._value = [0>>domain.subdomain()]*rows
                self._direction  = True
            elif isinstance(rows,Base_type.Base_type):
                if domain.subdomain().is_sub_element(rows):
                    self._value = [rows>>domain.subdomain()]
                    self._direction = True
                elif self.domain() == rows.domain():
                    self._value = rows._value
                    self._direction = rows._direction
                elif rows.domain()._Matrix and rows.domain().subdomain().is_sub_domain(domain.subdomain()):
                    self._structure = rows.domain()
                    self._value = rows._value[:]
                    self._direction = rows._direction
                else:
                    self._value = [rows._push_to_union(self.domain().subdomain())]
                    self._direction = True
            else:
                raise TypeError("Data not valid")
        else:
            if isinstance(cols,int):
                if cols == 1:
                    self._value = [0>>domain.subdomain()]*rows
                    self._direction  = True
                elif rows == 1:
                    self._value = [0>>domain.subdomain()]*cols
                    self._direction  = False
                else:
                    raise TypeError("Wrong parameters")
            elif cols == 'r':
                if isinstance(rows,int):
                    self._value = [0>>domain.subdomain()]*rows
                    self._direction  = True
                elif isinstance(rows,list):
                    self._value = []
                    self._direction = True
                    for i in range(len(rows)):
                        if self._structure.subdomain().is_sub_element(rows[i]):
                            self._value.append(rows[i]>> self.domain().subdomain())
                        else:
                            raise TypeError("One of the elements does not match with the domain")
                else:
                    raise TypeError("Wrong parameters")
            elif cols == 'c':
                if isinstance(rows,int):
                    self._value = [0>>domain.subdomain()]*rows
                    self._direction  = False
                elif isinstance(rows,list):
                    self._value = []
                    self._direction = False
                    for i in range(len(rows)):
                        if self._structure.subdomain().is_sub_element(rows[i]):
                            self._value.append(rows[i]>> self.domain().subdomain())
                        else:
                            self._value.append(rows[i]._push_to_union(self.domain().subdomain()))
                else:
                    raise TypeError("Wrong parameters")
            else:
                raise TypeError("Wrong parameters")
    
    def cp(self):
        if self._direction:
            return Vector_type(self.domain(), self._value[:], 'r')
        else:
            return Vector_type(self.domain(), self._value[:], 'c')
                
    def size(self):
        return len(self._value)
    
    def __len__(self):
        return self.size()
    
    def is_zero(self):
        if len(self._value) == 0:
            return False
        for i in range(len(self._value)):
            if self._value[i] != (0>>self._structure.subdomain()):
                return False
        return True
    
    def to_Matrix(self):
        return Matrix_type(self)
    
    def to_Diagonal_Matrix(self):
        aux = [[(0>>self.domain().subdomain())]*self.size() for i in range(self.size())]
        for i in range(self.size()):
            aux[i][i] = self._value[i]
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        return Matrix_type(D,aux)
    
    def to_list(self):
        return self._value[:]
        
    def _len_value_(self):
        return 1
    
    def _sign_value_(self):
        return 0
        
    def _print_(self):
        res = '['
        for i in range(len(self._value)):
            if isinstance(self._value[i],int):
                res = res + self._value[i].__str__()
            else:
                flagPar = False
#                if self._value[i]._len_value_() > 1:
#                    res += '('
#                    flagPar = True
                res += self._value[i]._print_()
                if flagPar:
                    res += ')'
            if self._direction:
                res += ', '
            else:
                res += ',\n '
        while(res[len(res)-1]!=','):
            if len(res) < 2:
                break
            res = res[:len(res)-1]
        if len(res)> 1:
            res=res[:len(res)-1]
        res += ']'
        return res
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._print_() + ' :: Vector[' + self._structure.subdomain()._print_() + ']' 
        else:
            return self._print_()
        
    def __repr__(self):
        return '\n' + self.__str__()
    
    def transpose(self):
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        if self._direction:
            return Vector_type(D,self._value,'c')
        else:
            return Vector_type(D,self._value,'r')
    
    def sub_vector(self,begin,end):
        aux = self._value[begin:end]
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        if self._direction:
            return Vector_type(D,aux,'r')
        else:
            return Vector_type(D,aux,'c')
    
    def __getitem__(self,key):
        fvec = False
        if isinstance(key,slice):
            a0 = key.start
            if a0 == None:
                a0 = 0
                fvec = True
            a1 = key.stop
            if a1 == None:
                a1 = self.size()
                fvec = True
            a2 = key.step
            if a2 == None:
                a2 = 1
            indexs = list(range(a0,a1,a2))
        if isinstance(key,int):
            indexs = [key]
        if isinstance(key,list):
            indexs = key
        #indexs.sort()
        aux = []
        if len(indexs) == 1:
            if fvec:
                aux = [copy.deepcopy(self._value.__getitem__(indexs[0]))]
                D = self.domain()._create_Matrix_(self.domain().subdomain())
                if self._direction:
                    return Vector_type(D,aux)
                else:
                    return Vector_type(D,aux, 'c')
            else:
                return copy.deepcopy(self._value.__getitem__(indexs[0]))
        for i in indexs:
            if i >= self.size():
                aux.append(0>>self.domain().subdomain())
            elif i < 0:
                aux.append(0>>self.domain().subdomain())
            else:
                aux.append(self._value[i])
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        if self._direction:
            return Vector_type(D,aux)
        else:
            return Vector_type(D,aux, 'c')
        
    
    def __setitem__(self,key,value): 
        if isinstance(key,slice):
            a0 = key.start
            if a0 == None:
                a0 = 0
            a1 = key.stop
            if a1 == None:
                a1 = self.size()
            a2 = key.step
            if a2 == None:
                a2 = 1
            indexs = list(range(a0,a1,a2))
        if isinstance(key,int):
            indexs = [key]
        if isinstance(key,list):
            indexs = key
        #indexs.sort()
        flagRep = False
        if isinstance(value,list):
            if len(value)!= len(indexs):
                raise TypeError("Value size does not match")
            for i in range(len(indexs)):
                if isinstance(value[i],float):
                    value[i] = value[i]*(1>>self.domain())
                if self.domain().is_sub_element(value[i]):
                    self._value[indexs[i]] = value[i]>>self.domain().subdomain()
                else:
                    flagRep = True
                    self.domain().__init__((self.domain().union_domain(value[i].domain())[0]).subdomain())
                    self._value[indexs[i]] = value[i].cp()._push_to_union(self.domain().subdomain())
        elif isinstance(value,Vector_type):
            if len(indexs)!=value.size():
                raise TypeError("Sizes do not match")
            if value.domain() != self.domain():
                aux = value.to_list()
                self.__setitem__(key,aux)
                return
            for i in range(len(indexs)):
                self._value[indexs[i]] = value._value[i]>>self.domain().subdomain()
        elif isinstance(value,Matrix_type):
            vv = value._toVector()
            self.__setitem__(key,vv)
        elif isinstance(value,int):
            for i in indexs:
                self._value[i] = value>>self.domain().subdomain()
        elif isinstance(value,float):
            v = value*(1>>self.domain())
            if self.domain().is_sub_element(v):
                for i in indexs:
                    self._value[i] = v>>self.domain().subdomain()
            else:
                flagRep = True
                self.domain().__init__((self.domain().union_domain(v.domain())[0]).subdomain())
                for i in indexs:
                    self._value[i] = v._push_to_union(self.domain().subdomain())
        elif isinstance(value,Base_type.Base_type):
            if self.domain().is_sub_element(value):
                for i in indexs:
                    self._value[i] = value>>self.domain().subdomain()
            elif value.domain().is_sub_domain(self.domain().subdomain()):
                flagRep = True
                self.domain().__init__((self.domain().union_domain(value.domain())[0]).subdomain())
                for i in indexs:
                    self._value[i] = value.cp()._push_to_union(self.domain().subdomain())
        else:
            raise TypeError("Value not valid")
        if flagRep:
            for i in range(len(self._value)):
                if self.domain().is_sub_element(self._value[i]):
                    self._value[i] = self._value[i]>>self.domain().subdomain()
                else:
                    self._value[i] = self._value[i]._push_to_union(self.domain().subdomain())
    
    def __hash__(self):
        f = self.transform_to_min_domain()
        if isinstance(f,int):
            return hash(f)
        elif f.domain() != self.domain():
            return hash(f)
        
        aux = 0
        for i in range(self.size()):
            if self._direction:
                aux += hash(self._value[i])*PGl.HMATR**i
            else:
                aux *= hash(self._value[i])*PGl.HMATC**i
        return aux
        
    
    def __neg__(self):
        return self._structure._opposite_(self)
    
    def opposite(self):
        return self._structure._opposite_(self)
    
    def negation(self):
        return self._structure._opposite_(self)
    
    def is_invertible(self):
        return self._structure.is_invertible(self)
    
    def min_domain(self):
        if len(self._value) == 0:
            if isinstance(self._value[0],int):
                return self.domain()._ZZ_()
            return self._value[0].min_domain()
        if isinstance(self._value[0],int):
            D = self.domain()._ZZ_()
        else:
            D = self._value[0].min_domain()
        if D == self.domain().subdomain():
            return self.domain()._create_Matrix_(self.domain().subdomain())
        for i in range(1,len(self)):
            if isinstance(self._value[i],int):
                continue
            DS = self._value[i].min_domain()
            if DS == self.domain().subdomain():
                return self.domain()._create_Matrix_(self.domain().subdomain())
            if not D.is_sub_domain(DS):
                D = DS
        return self.domain()._create_Matrix_(D)
    
    def transform_to_min_domain(self):
        D = self.min_domain()
        if D == self.domain():
            if self._direction:
                return Vector_type(D,self._value,'r')
            else:
                return Vector_type(D,self._value,'c')
        if self.domain().is_sub_domain(D):
            if isinstance(self._value[0],int):
                return self._value[0]
            return self._value[0].transform_to_min_domain()
        aux = []
        for i in range(len(self)):
            if isinstance(self._value[i],int):
                aux.append(self._value[i]>>D.subdomain())
            else:
                aux.append(self._value[i].transform_to_min_domain()>>D.subdomain())
        if self._direction:
            return Vector_type(D,aux)
        else:
            return Vector_type(D,aux,'c')

    def change_characteristic(self,n):
        aux =[]
        for i in range(len(self._value)):
            aux.append(self._value[i].change_characteristic(n))
        if self._direction:
            return Vector_type(self.domain().change_characteristic(n),aux,'r') 
        else:
            return Vector_type(self.domain().change_characteristic(n),aux,'c')             
            
    def transform_to_ZZ(self):
        if self.domain().characteristic() == 0:
            return self
        return self.change_characteristic(0)
    
    def K_(self):
        return self.domain().K_()

    def __pow__(self,other):
        if isinstance(other, int):
            if other<0:
                return (1/self)**(-other)
            if other == 0:
                D = self.domain()._create_Matrix_(self.domain().subdomain())
                if self._direction:
                    return Vector_type(D,[1>>self.domain().subdomain()]*self.size())
                else:
                    return Vector_type(D,[1>>self.domain().subdomain()]*self.size(),'c')
            D = self.domain()._create_Matrix_(self.domain().subdomain())
            if self._direction:
                aux1 = Vector_type(D,[1>>self.domain().subdomain()]*self.size())
            else:
                aux1 = Vector_type(D,[1>>self.domain().subdomain()]*self.size(),'c')

            aux2 = self
            while other>1:
                if other%2==0:
                    other /= 2
                    aux2 = aux2*aux2
                else:
                    aux1 = aux2*aux1
                    aux2 = aux2*aux2
                    other = (other-1)/2
            return aux1 * aux2
        else:
            raise TypeError("This element cannot be operated")
    

class Matrix_type(Base_type.Base_type):
    """Matrix Element class. """
        
    def __init__(self, domain, rows = '', cols = ''):
        '''With 3 parameters, it creates a zero rowsxcols matrix of domain
           With 2 parameters, it inizializates matrix with rows
        '''
        self.L = None
        self.P = None
        self.U = None
        self.Piv = None
        self.PivP = None
        self.signPLU = 1
        self.signPiv = 1
        self._structure = domain
        if rows == '':
            if isinstance(domain,Matrix_type):
                self._structure = domain.domain()
                self._value = [domain._value[i][:] for i in range(len(domain._value))]
                self._rows = domain._rows
                self._cols = domain._cols
            elif isinstance(domain,Vector_type):
                self._structure = domain.domain()
                if domain._direction:
                    self._value = [domain._value[:]]
                    self._rows = 1
                    self._cols = domain.size()
                else:
                    self._value = [[domain._value[0]]]
                    for i in range(1,domain.size()):
                        self._value.append([domain._value[i]])
                    self._rows = domain.size()
                    self._cols = 1
            else:
                self._value = [[]]#[[0>>domain.subdomain()]]
                self._rows = 0
                self._cols = 0
        elif cols == '':
            if isinstance(rows,list):
                if len(rows) > 0:
                    if isinstance(rows[0],list):
                        self._value =[rows[i][:] for i in range(len(rows))]
                    else:
                        self._value = [rows[:]]
                    for i in range(len(self._value)):
                        if len(self._value[i])!=len(self._value[0]):
                            raise TypeError("List of lists always has to have the same lenght")
                        for j in range(len(self._value[0])):
                            if self._structure.subdomain().is_sub_element(self._value[i][j]):
                                self._value[i][j] = self._value[i][j]>> self.domain().subdomain()
                            else:
                                self._value[i][j] = self._value[i][j]._push_to_union(self.domain().subdomain())
                else:
                    self._value = [[]]
            elif isinstance(rows,int):
                self._value = [[rows>>domain.subdomain()]]
            elif isinstance(rows,Base_type.Base_type):
                if domain.subdomain().is_sub_element(rows):
                    self._value[[rows>>domain.subdomain()]]
                else:
                    self._value[[rows._push_to_union(self.domain().subdomain())]]
            else:
                raise TypeError("Data not valid")
            if self._value == [[]]:
                self._rows = 0
                self._cols = 0
            else:
                self._rows = len(self._value)
                self._cols = len(self._value[0])            
        else:
            if cols == 0 or rows == 0:
                self._value = [[]]
                self._rows = 0
                self._cols = 0
            else:
                self._value = [[(0>>domain.subdomain())]*cols for i in range(rows)]
                self._rows = rows
                self._cols = cols
    
    def cp(self):
        return Matrix_type(self)
        
    def is_zero(self):
        if self._cols == 0 or self._rows == 0:
            return False
        for i in range(self._rows):
            for j in range(self._cols):
                if self._value[i][j] != (0>>self._structure.subdomain()):
                    return False
        return True
        
    def _len_value_(self):
        return 1
    
    def _sign_value_(self):
        return 0
    
    def __len__(self):
        return self.nrows()
        
    def _print_(self):
        res = '['
        for i in range(self._rows):
            res +='['
            for j in range(self._cols):
                if isinstance(self._value[i][j],int):
                    res = res + self._value[i][j].__str__()
                else:
                    flagPar = False
#                    if self._value[i][j]._len_value_() > 1:
#                        res += '('
#                        flagPar = True
                    res += self._value[i][j]._print_()
                    if flagPar:
                        res += ')'
                res += '\t'
            res=res[:len(res)-1]
            res += ']\n '
        if len(res) > 1:
            res=res[:len(res)-2]
        res += ']'
        return res
    
    def __str__(self):
        if PGl.TYPES == 1:
            return self._print_() + ' :: Matrix[' + self._structure.subdomain()._print_() + ']' 
        else:
            return self._print_()
        
    def __repr__(self):
        return '\n' +self.__str__()
    
    def transpose(self):
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        return Matrix_type(D,list(map(list,zip(*self._value))))
    
    def det(self):
        return self._structure._det_(self)
    
    def trace(self):
        return self._structure._trace_(self)
    
    def is_invertible(self):
        return self._structure.is_invertible(self)
        
    def _toVector(self):
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        if self._rows == 1:
            return Vector_type(D,self._value[0],'r')
        elif self._cols == 1:
            aux = []
            for i in range(self._rows):
                aux.append(self._value[i][0])
            return Vector_type(D,aux,'c')
        else:
            raise TypeError("It is not possible to transform to a Vector")
    
    def ncols(self):
        return self._cols
    
    def nrows(self):
        return self._rows
       
    def col(self,n):
        A = self.transpose()
        A = [A._value[n]]
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        return (Matrix_type(D,A)).transpose()._toVector()
    
    
    def row(self,n):
        A = self._value[n]
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        return Vector_type(D,A)
    
    def pivot_matrix(self):
        return self.domain().pivot_matrix(self)
    
    def sub_matrix(self,rows=[],cols=[]):
        aux = []
        if rows == []:
            rows = range(self._rows)
        if isinstance(rows,int):
            rows = [rows]
        if cols == []:
            cols = range(self._cols)
        if isinstance(cols,int):
            cols = [cols]
        for i in rows:
            aux2 = []
            if i>=self._rows:
                for j in cols:
                    aux2.append(0>>self.domain().subdomain())
            elif i < 0:
                for j in cols:
                    aux2.append(0>>self.domain().subdomain())
            else:
                for j in cols:
                    if j>=self._cols:
                        aux2.append(0>>self.domain().subdomain())
                    elif j < 0:
                        aux2.append(0>>self.domain().subdomain())
                    else:
                        aux2.append(self._value[i][j])
            aux.append(aux2)
        D = self.domain()._create_Matrix_(self.domain().subdomain())
        return Matrix_type(D,aux)
    
    def __getitem__(self,key):
        if isinstance(key,tuple):
            if len(key) == 2:
                if isinstance(key[0],slice):
                    a0 = key[0].start
                    if a0 == None:
                        a0 = 0
                    a1 = key[0].stop
                    if a1 == None:
                        a1 = self._rows
                    a2 = key[0].step
                    if a2 == None:
                        a2 = 1
                    rows = list(range(a0,a1,a2))
                if isinstance(key[1],slice):
                    a0 = key[1].start
                    if a0 == None:
                        a0 = 0
                    a1 = key[1].stop
                    if a1 == None:
                        a1 = self._cols
                    a2 = key[1].step
                    if a2 == None:
                        a2 = 1
                    cols = list(range(a0,a1,a2))
                if isinstance(key[0],int):
                    rows = key[0]
                if isinstance(key[1],int):
                    cols = key[1]
                if isinstance(key[0], list):
                    rows = key[0]
                if isinstance(key[1], list):
                    cols = key[1]
                #if isinstance(rows,list):
                    #rows.sort()
                #if isinstance(cols, list):
                    #cols.sort()
                if isinstance(rows,list) and isinstance(cols,list):
                    return self.sub_matrix(rows,cols)
                if isinstance(rows,list) and isinstance(cols, int):
                    if len(rows) == self._rows:
                        return self.col(cols)
                    else:
                        aux = []
                        for i in rows:
                            aux.append(self._value[i][cols])
                        D = self.domain()._create_Matrix_(self.domain().subdomain())
                        return Vector_type(D,aux,'c')
                if isinstance(rows, int) and isinstance(cols, list):
                    if len(cols) == self._cols:
                        return self.row(rows)
                    else:
                        aux = []
                        for i in cols:
                            aux.append(self._value[rows][i])
                        D = self.domain()._create_Matrix_(self.domain().subdomain())
                        return Vector_type(D,aux)
                return self._value[key[0]][key[1]]
            else:
                raise TypeError("Too much dimensions")
        else:
            A = self._value.__getitem__(key)
            if isinstance(A,list):
                D = self.domain()._create_Matrix_(self.domain().subdomain())
                return Vector_type(D,A)
            else:
                return A
        
    
    def __setitem__(self,key,value):
        self.L = None
        self.P = None
        self.U = None
        self.Piv = None
        self.PivP = None
        self.signPLU = 1
        self.signPiv = 1
        flagRep = False
        if isinstance(key,tuple):
            if len(key) == 2:
                if isinstance(key[0],slice):
                    a0 = key[0].start
                    if a0 == None:
                        a0 = 0
                    a1 = key[0].stop
                    if a1 == None:
                        a1 = self._rows 
                    a2 = key[0].step
                    if a2 == None:
                        a2 = 1
                    rows = list(range(a0,a1,a2))
                if isinstance(key[1],slice):
                    a0 = key[1].start
                    if a0 == None:
                        a0 = 0
                    a1 = key[1].stop
                    if a1 == None:
                        a1 = self._cols
                    a2 = key[1].step
                    if a2 == None:
                        a2 = 1
                    cols = list(range(a0,a1,a2))
                if isinstance(key[0],int):
                    rows = [key[0]]
                if isinstance(key[1],int):
                    cols = [key[1]]
                if isinstance(key[0], list):
                    rows = key[0]
                if isinstance(key[1], list):
                    cols = key[1]
                if isinstance(value,list):
                    if isinstance(value[0],list):
                        if len(value[0]) != len(cols) or len(value)!= len(rows):
                            raise TypeError("Value size does not match")
                        for i in range(len(rows)):
                            for j in range(len(cols)):
                                if isinstance(value[i][j],float):
                                    value[i][j] = value[i][j]*(1>>self.domain())
                                if self.domain().is_sub_element(value[i][j]):
                                    self._value[rows[i]][cols[j]] = value[i][j]>>self.domain().subdomain()
                                else:
                                    flagRep = True
                                    self.domain().__init__((self.domain().union_domain(value[i][j].domain())[0]).subdomain())
                                    self._value[rows[i]][cols[j]] = value[i][j]._push_to_union(self.domain().subdomain())
                                    
#                                elif value[i][j].domain().is_sub_domain(self.domain().subdomain()):
#                                    flagRep = True
#                                    self.domain().__init__(value[i][j].domain())
#                                    self._value[rows[i]][cols[j]] = value[i][j]>>self.domain().subdomain()
                    else:
                        if len(value)!= len(cols)*len(rows):
                            raise TypeError("Value size does not match")
                        for i in range(len(rows)):
                            for j in range(len(cols)):
                                if isinstance(value[i*len(cols) + j],float):
                                    value[i*len(cols) + j] = value[i*len(cols) + j]*(1>>self.domain().subdomain())
                                if self.domain().is_sub_element(value[i]):
                                    self._value[rows[i]][cols[j]] = value[i*len(cols) + j]>>self.domain().subdomain()
                                else:
                                    flagRep = True
                                    self.domain().__init__((self.domain().union_domain(value[i*len(cols) + j].domain())[0]).subdomain())
                                    self._value[rows[i]][cols[j]] = value[i*len(cols) + j]._push_to_union(self.domain().subdomain())

#                                elif value[i][j].domain().is_sub_domain(self.domain().subdomain()):
#                                    flagRep = True
#                                    self.domain().__init__(value[i].domain())
#                                    self._value[rows[i]][cols[j]] = value[i*len(cols) + j]>>self.domain().subdomain()
                elif isinstance(value,Vector_type):
                    if len(rows)>1 and len(cols)>1:
                        raise TypeError("Sizes do not match")
                    if len(rows)*len(cols)!=value.size():
                        raise TypeError("Sizes do not match")
                    if value.domain() != self.domain():
                        aux = value.to_list()
                        self.__setitem__(key,aux)
                        return
                    if len(rows)==1:
                        for i in range(len(cols)):
                            self._value[rows[0]][cols[i]] = value._value[i]>>self.domain().subdomain()
                    else:
                        for i in range(len(rows)):
                            self._value[rows[i]][cols[0]] = value._value[i]>>self.domain().subdomain()
                elif isinstance(value,Matrix_type):
                    if len(rows)!=value._rows or len(cols)!=value._cols:
                        raise TypeError("Sizes do not match")
                    if value.domain() != self.domain():
                        aux = value.to_list()
                        self.__setitem__(key,aux)
                        return
                    for i in range(len(rows)):
                        for j in range(len(cols)):
                            self._value[rows[i]][cols[j]] = value._value[i][j]>>self.domain().subdomain()
                elif isinstance(value,int):
                    for i in rows:
                        for j in cols:
                            self._value[i][j] = value>>self.domain().subdomain()
                elif isinstance(value,float):
                    v= value*(1>>self.domain())
                    if self.domain().is_sub_element(v):
                        for i in rows:
                            for j in cols:
                                self._value[i][j] = v>>self.domain().subdomain()
                    else:
                        flagRep = True
                        self.domain().__init__((self.domain().union_domain(v.domain())[0]).subdomain())
                        for i in rows:
                            for j in cols:
                                self._value[i][j] = v._push_to_union(self.domain().subdomain())

                elif isinstance(value,Base_type.Base_type):
                    if not self.domain().is_sub_domain(value.domain()):
                        flagRep = True
                        self.domain().__init__((self.domain().union_domain(value.domain())[0]).subdomain())
                    for i in rows:
                        for j in cols:
                            self._value[i][j] = value._push_to_union(self.domain().subdomain())
                else:
                    raise TypeError("Value not valid")
            else:
                raise TypeError("Too much dimensions")
        else:
            if isinstance(value,list):
                if len(value) != self._cols:
                    raise TypeError("Wrong number of elements")
                aux = []
                for i in range(self._cols):
                    if isinstance(value[i],float):
                        value[i]= value[i]*(1>>self.domain())
                    if self.domain().is_sub_element(value[i]):
                        aux.append(value[i]>>self.domain().subdomain())
                    else:
                        flagRep = True
                        self.domain().__init__((self.domain().union_domain(value[i].domain())[0]).subdomain())
                        aux.append(value[i]._push_to_union(self.domain().subdomain()))
                self._value.__setitem__(key,aux)
            elif isinstance(value,Vector_type):
                if value.size() != self._cols:
                    raise TypeError("Wrong number of elements")
                if value.domain() != self.domain():
                    aux = value.to_list()
                    self.__setitem__(key,aux)
                    return
                self._value.__setitem__(key,value._value)
            elif isinstance(value,Matrix_type):
                if value._rows == 1:
                    if value.domain() != self.domain():
                        aux = value.to_list()
                        self.__setitem__(key,aux[0])
                        return
                    self._value.__setitem__(key,value._value[0])
                else:
                    raise TypeError("Size not valid")
            elif isinstance(value,int):
                for i in range(self._cols):
                    self._value[key][i] = value>>self.domain().subdomain()
            elif isinstance(value,float):
                v= value*(1>>self.domain())
                if self.domain().is_sub_element(v):
                    for i in range(self._rows):
                        self._value[key][i] = v>>self.domain().subdomain()
                else:
                    flagRep = True
                    self.domain().__init__((self.domain().union_domain(v.domain())[0]).subdomain())
                    for i in range(self._rows):
                        self._value[key][i] = v._push_to_union(self.domain().subdomain())

            elif isinstance(value,Base_type.Base_type):
                if not self.domain().is_sub_domain(value.domain()):
                    flagRep = True
                    self.domain().__init__((self.domain().union_domain(value.domain())[0]).subdomain())
                for i in range(self._cols):
                    self._value[key][i] = value._push_to_union(self.domain().subdomain())
        if flagRep:
            for i in range(self._rows):
                for j in range(self._cols):
                    if self.domain().is_sub_element(self._value[i][j]):
                        self._value[i][j] = self._value[i][j]>>self.domain().subdomain()
                    else:
                        self._value[i][j] = self._value[i][j]._push_to_union(self.domain().subdomain())
    
    def find_index_col(self,a):
        if isinstance(a,list):
            D = self.domain()._create_Matrix_(self.domain().subdomain())
            A = Vector_type(D,a,'c')
        if isinstance(a,Vector_type):
            A = a
        if isinstance(a,Matrix_type):
            if a._cols == 1:
                A = a._toVector()
            elif a._rows == 1:
                A = a._toVector()
            else:
                return -1
        else:
            D = self.domain()._create_Matrix_(self.domain().subdomain())
            A = Vector_type(D,a)
        for i in range(self._cols):
            if A == self.col(i):
                return i
        return -1    
    
    def find_index_row(self,a):
        if isinstance(a,list):
            D = self.domain()._create_Matrix_(self.domain().subdomain())
            A = Vector_type(D,a)
        if isinstance(a,Vector_type):
            A = a
        if isinstance(a,Matrix_type):
            if a._cols == 1:
                A = a._toVector()
            elif a._rows == 1:
                A = a._toVector()
            else:
                return -1
        else:
            D = self.domain()._create_Matrix_(self.domain().subdomain())
            A = Vector_type(D,a)
        for i in range(self._rows):
            if A == self.row(i):
                return i
        return -1    
    
    def __hash__(self):
        f = self.transform_to_min_domain()
        if isinstance(f,int):
            return hash(f)
        elif f.domain() != self.domain():
            return hash(f)
        
        aux = 0
        for i in range(self.nrows()):
            for j in range(self.ncols()):
                aux *= hash(self._value[i][j])*PGl.HMATC**i*PGl.HMATC**j
        return aux
                    
    def __neg__(self):
        return self._structure._opposite_(self)
    
    def opposite(self):
        return self._structure._opposite_(self)
    
    def negation(self):
        return self._structure._opposite_(self)
    
    def to_list(self):
        return self._value
    
    def min_domain(self):
        if self._rows == 1 and self._cols == 1:
            return self._value[0][0].min_domain()
        else:
            if isinstance(self._value[0][0],int):
                D = self.domain()._ZZ_()
            else:
                D = self._value[0][0].min_domain()
            if D == self.domain().subdomain():
                return self.domain()
            for i in range(self.nrows()):
                for j in range(self.ncols()):
                    if isinstance(self._value[i][j],int):
                        continue
                    if self._value[i][j].min_domain() == self.domain().subdomain():
                        D = self.domain()._create_Matrix_(self.domian().subdomain())
                        return D
                    if self._value[i][j].min_domain().is_sub_domain(D):
                        D = self._value[i][j].min_domain()
            return self.domain()._create_Matrix_(D)
    
    def transform_to_min_domain(self):
        DM = self.min_domain()
        if self.domain() != DM:
            if self.domain().is_sub_domain(DM):
                if isinstance(self._value[0][0],int):
                    return self._value[0][0]
                return self._value[0][0].transform_to_min_domain()
            else:
                aux = []
                for i in range(self.nrows()):
                    aux2 = []
                    for j in range(self.ncols()):
                        if isinstance(self._value[i][j],int):
                            aux2.append(self._value[i][j]>>DM.subdomain())
                        else:
                            aux2.append((self._value[i][j].transform_to_min_domain())>>DM.subdomain())
                    aux.append(aux2)
                return Matrix_type(DM,aux)
        elif self._rows == 1:
            return self._toVector()
        elif self._cols == 1:
            return self._toVector()
        else:
            return self.cp()

    def change_characteristic(self,n):
        aux =[]
        for i in range(self._rows):
            aux2 = []
            for j in range(self._cols):
                aux2.append(self._value[i][j].change_characteristic(n))
            aux.append(aux2)
        return Matrix_type(self.domain().change_characteristic(n),aux)                
            
    def transform_to_ZZ(self):
        if self.domain().characteristic() == 0:
            return self.cp()
        return self.change_characteristic(0)
    
    def rank(self):
        return self.domain()._rank_(self)
    
    def lu(self):
        return self.domain().gaussian_rows_elimination(self)
                    
    def __pow__(self,other):
        if isinstance(other, int):
            if self._cols != self._rows:
                raise TypeError("Only square matrices can be operated by **")
            if other<0:
                return (1/self)**(-other)
            if other == 0:
                D = self.domain()._create_Matrix_(self.domain().subdomain())
                return (Vector_type(D,[1>>self.domain().subdomain()]*self._rows)).to_Diagonal_Matrix()
            D = self.domain()._create_Matrix_(self.domain().subdomain())
            aux1 = (Vector_type(D,[1>>self.domain().subdomain()]*self._rows)).to_Diagonal_Matrix()

            aux2 = self
            while other>1:
                if other%2==0:
                    other /= 2
                    aux2 = aux2*aux2
                else:
                    aux1 = aux2*aux1
                    aux2 = aux2*aux2
                    other = (other-1)/2
            return aux1 * aux2
        else:
            raise TypeError("This element cannot be operated")
    
    def rank_block(self):
        return self.domain().rank_block(self)
    
    def prune(self):
        return self.domain().prune(self)
    
    def characteristic(self):
        return self.domain().characteristic()
    
    def K_(self):
        return self.domain().K_()

