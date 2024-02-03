## pfactor: factorizing polynomials

from PyM.PyWIT.PyECC.Fq import *
from PyM.PyWIT.PyECC.Power_Series import *
from PyM.PyWIT.PyECC.pfactorFunctions import *


'''
Theorem: Let K be a finite fiel of cardinal q, f in K[X], and n = degree(f). 
Then f is irreducible iff
1) f divides X**(q**n), and
2) gcd(X**(q**(n//p)), f) = 1 for all prime divisors p of n.

See vzGathen-Gerhard-2003, pp 396-397. 
'''

def univariate_polynomial_ring(K,symbol = 'X', name = ''):
    P = Poly(K,name)
    x = P.variable([symbol]);
    return [P,x]
    
def create_vector(elem,size='',direction=''):
    '''creates a Vector_type
         1 parameter:
             - elem == list Transforms to row vector
             - elem == domain returns 1x1 zero row vector
             - elem == Base_type transform to a 1x1 row vector
         2 parameters:
             - elem == domain && size == list --> same as 1 parameter with elem = size>>elem
             - size = int --> same as 1 parameter but with elem = [elem]*direction
             - size = 'r' --> same as 1 parameter
             - size = 'c' --> same as 1 parameter but with column vector
         3 parameters:
             - direction = 'r' --> same as 2 parameters
             - direction = 'c' --> same as 2 parameters but column vector
    '''
    if size == '':
        if isinstance(elem,list):
            if len(elem)==0:
                return Vector_type(Matrix(ZZ()),[])
            if isinstance(elem[0],int):
                domain = ZZ()
            elif isinstance(elem[0],float):
                domain = QQ(ZZ())
            elif isinstance(elem[0],Base_type.Base_type):
                domain = elem[0].domain()
            elif isinstance(elem[0], Symbol_type):
                domain = Poly(ZZ())
            else:
                raise TypeError("Some elements of the list are not valid")
            for i in range(1,len(elem)):
                if isinstance(elem[i],int):
                    domaux = ZZ()
                elif isinstance(elem[i],float):
                    domaux = QQ(ZZ())
                elif isinstance(elem[i], Symbol_type):
                    domaux = Poly(ZZ())
                else:
                    domaux = elem[i].domain()
                if domaux!=domain:
                    if domain.is_sub_domain(domaux):
                        continue
                    elif domaux.is_sub_domain(domain):
                        domain = domaux
                    elif domain == ZZ() and domaux != ZZ():
                        domain = domaux
                    elif domaux == ZZ():
                        continue
                    else:
                        domain = domain.union_domain(domaux)[0]
            DomM = Matrix(domain)
            return Vector_type(DomM,elem)
        elif isinstance(elem,Vector_type):
            return Vector_type(elem)
        elif isinstance(elem,Matrix_type):
            return Vector_type(elem)
        elif isinstance(elem,Base_type.Base_type):
            DomM = Matrix(elem.domain())
            return Vector_type(DomM,[elem])
        elif isinstance(elem, int):
            DomM = Matrix(ZZ())
            return Vector_type(DomM, elem)
        elif isinstance(elem, float):
            DomM = Matrix(QQ(ZZ()))
            return Vector_type(DomM,elem)
        elif isinstance(elem,Symbol_type):
            DomM = Matrix(Poly(ZZ()))
            return Vector_type(DomM,elem._to_monomial_())
        elif isinstance(elem,Base.Base):
            if elem._matrix_structure:
                return Vector_type(elem,[])
            else:
                DomM = Matrix(elem)
                return Vector_type(elem,[])
        else:
            raise TypeError("Parameter elem not valid")
    elif direction == '':
        if size == 'r':
            return create_vector(elem)
        elif size == 'c':
            if isinstance(elem,list):
                if len(elem) == 0:
                    return Vector_type(Matrix(ZZ()),[],'c')
                if isinstance(elem[0],int):
                    domain = ZZ()
                elif isinstance(elem[0],float):
                    domain = QQ(ZZ())
                elif isinstance(elem[0],Base_type.Base_type):
                    domain = elem[0].domain()
                elif isinstance(elem[0], Symbol_type):
                    domain = Poly(ZZ())
                else:
                    raise TypeError("Some elements of the list are not valids")
                for i in range(1,len(elem)):
                    if isinstance(elem[i],int):
                        domaux = ZZ()
                    elif isinstance(elem[i],float):
                        domaux = QQ(ZZ())
                    elif isinstance(elem[i], Symbol_type):
                        domaux = Poly(ZZ())
                    else:
                        domaux = elem[i].domain()
                    if domaux!=domain:
                        if domain.is_sub_domain(domaux):
                            continue
                        elif domaux.is_sub_domain(domain):
                            domain = domaux
                        else:
                            domain = domain.union_domain(domaux)[0]
                DomM = Matrix(domain)
                return Vector_type(DomM,elem,'c')
            elif isinstance(elem,Vector_type):
                return Vector_type(elem.domain(),elem._value,'c')
            elif isinstance(elem,Matrix_type):
                return Vector_type(elem)
            elif isinstance(elem,Base_type.Base_type):
                DomM = Matrix(elem.domain())
                return Vector_type(DomM,[elem],'c')
            elif isinstance(elem, int):
                DomM = Matrix(ZZ())
                return Vector_type(DomM, [elem],'c')
            elif isinstance(elem, float):
                DomM = Matrix(QQ(ZZ()))
                return Vector_type(DomM, [elem],'c')
            elif isinstance(elem,Symbol_type):
                DomM = Matrix(Poly(ZZ()))
                return Vector_type(DomM, [elem._to_monomial_()],'c')
            elif isinstance(elem,Base.Base):
                if elem._matrix_structure:
                    return Vector_type(elem,[],'c')
                else:
                    DomM = Matrix(elem)
                    return Vector_type(DomM,[],'c')
            else:
                raise TypeError("Parameter elem not valid")
        elif isinstance(size,int):
            if isinstance(elem,Base.Base):
                if elem._matrix_structure:
                    return Vector_type(elem,size,1)
                else:
                    DomM = Matrix(elem)
                    return Vector_type(DomM,size,1)
            elif isinstance(elem,Vector_type):
                return Vector_type(elem.domain(),elem._value*size)
            elif isinstance(elem,Matrix_type):
                return Vector_type(elem)
            elif isinstance(elem,Base_type.Base_type):
                DomM = Matrix(elem.domain())
                return Vector_type(DomM,[elem]*size)
            elif isinstance(elem, int):
                DomM = Matrix(ZZ())
                return Vector_type(DomM,[elem]*size) 
            elif isinstance(elem, float):
                DomM = Matrix(QQ(ZZ()))
                return Vector_type(DomM,[elem]*size) 
            elif isinstance(elem, Symbol_type):
                DomM = Matrix(Poly(ZZ()))
                return Vector_type(DomM,[elem._to_monomial_()]*size) 
            elif isinstance(elem,list):
                return create_vector(elem*size)
            else:
                raise TypeError("Wrong parameters")
        elif isinstance(size,list):
            if isinstance(elem,Base.Base):
                aux = []
                for i in range(len(size)):
                    aux.append(size[i]>>elem)
                return create_vector(aux)
            else:
                raise TypeError("Wrong parameters")
        else:
            raise TypeError("Wrong parameters")
    elif direction == 'r':
        return create_vector(elem,size)
    elif direction == 'c':
        if isinstance(size,list):
            if len(size) == 0:
                    return Vector_type(Matrix(ZZ()),[],'c')
            if isinstance(elem,Base.Base):
                aux = []
                for i in range(len(size)):
                    aux.append(size[i]>>elem)
                return create_vector(aux,'c')
            else:
                raise TypeError("Wrong parameters")
        if isinstance(elem,Base.Base):
            if elem._matrix_structure:
                return Vector_type(elem,1,size)
            else:
                DomM = Matrix(elem)
                return Vector_type(elem,1,size)
        elif isinstance(elem,Vector_type):
            return Vector_type(elem.domain(),elem._value*size,'c')
        elif isinstance(elem,Matrix_type):
            return Vector_type(elem)
        elif isinstance(elem,Base_type.Base_type):
            DomM = Matrix(elem.domain())
            return Vector_type(DomM,[elem]*size,'c')
        elif isinstance(elem, int):
            DomM = Matrix(ZZ())
            return Vector_type(DomM,[elem]*size,'c') 
        elif isinstance(elem, float):
            DomM = Matrix(QQ(ZZ()))
            return Vector_type(DomM,[elem]*size,'c') 
        elif isinstance(elem, Symbol_type):
            DomM = Matrix(Poly(ZZ()))
            return Vector_type(DomM,[elem._to_monomial_()]*size,'c') 
        elif isinstance(elem,list):
            if isinstance(elem[0],int):
                domain = ZZ()
            elif isinstance(elem[0],float):
                domain = QQ(ZZ())
            elif isinstance(elem[0],Base_type.Base_type):
                domain = elem[0].domain()
            elif isinstance(elem[0], Symbol_type):
                domain = Poly(ZZ())
            else:
                raise TypeError("Some elements of the list are not valids")
            for i in range(1,len(elem)):
                if isinstance(elem[i],int):
                    domaux = ZZ()
                elif isinstance(elem[i],float):
                    domaux = QQ(ZZ())
                elif isinstance(elem[i], Symbol_type):
                    domaux = Poly(ZZ())
                else:
                    domaux = elem[i].domain()
                if domaux!=domain:
                    if domain.is_sub_domain(domaux):
                        continue
                    elif domaux.is_sub_domain(domain):
                        domain = domaux
                    else:
                        raise TypeError("Some elements domains cannot be compared")
            DomM = Matrix(domain)
            return Vector_type(DomM,elem*size,'c')
        else:
            raise TypeError("Wrong parameters")
            
def create_matrix(elem,rows='',cols='',s = ''):
    '''creates a Matrix_type
         1 parameter:
             - elem == list Transforms to matrix
             - elem == domain returns 1x1 zero matrix
             - elem == Base_type transform to a 1x1 matrix
         2 parameters:
             - elem == domain && rows == list --> same as 1 parameter with elem = rows>>elem
             - elem = Base_type & rows = int --> same as 1 parameter but with elem = [[elem]*rows]*rows
             - elem = Base & rows = int -> rowsxrows zero matrix
         3 parameters:
             - elem == domain:
                 - rows & cols == int, zero>>elem rowsxcols matrix
             - elem == Base_type:
                 - rows & cols == int, returns a rowsxcols matrix with each element = elem
             - elem == list:
                 - rows & cols == int, returns a reshape of the list to a rowsxcols matrix
         4 parametres:
             - elem == domain, rows == list, cols & s == int, returns a reshape of the list>>domain to a colsxs matrix
    '''
    if rows == '':
        if isinstance(elem,list):
            if len(elem) == 0:
                return Matrix_type(Matrix(ZZ()),[[]])
            if isinstance(elem[0],Vector_type):
                aux = elem[0].cp()
                for i in range(1,len(elem)):
                    aux = stack(aux,elem[i])
                return aux
            if not isinstance(elem[0],list):
                return create_matrix([elem])
            if isinstance(elem[0][0],int):
                domain = ZZ()
            elif isinstance(elem[0][0],float):
                domain = QQ(ZZ())
            elif isinstance(elem[0][0],Symbol_type):
                domain = Poly(ZZ())
            elif isinstance(elem[0][0],Base_type.Base_type):
                domain = elem[0][0].domain() #FALTA list,list
            else:
                raise TypeError("Some elements of the list are not valid")
            for i in range(len(elem)):
                if len(elem[i])!=len(elem[0]):
                    raise TypeError("Bad list")
                for j in range(len(elem[0])):
                    if isinstance(elem[i][j],int):
                        domaux = ZZ()
                    elif isinstance(elem[i][j],float):
                        domaux = QQ(ZZ())
                    elif isinstance(elem[i][j],Symbol_type):
                        domaux = Poly(ZZ())
                    else:
                        domaux = elem[i][j].domain()
                    if domaux!=domain:
                        if domain.is_sub_domain(domaux):
                            continue
                        elif domaux.is_sub_domain(domain):
                            domain = domaux
                        else:
                            domain = domain.union_domain(domaux)[0]
            DomM = Matrix(domain)
            return Matrix_type(DomM,elem)
        elif isinstance(elem,Vector_type):
            return Matrix_type(elem)
        elif isinstance(elem,Matrix_type):
            return elem
        elif isinstance(elem,Base_type.Base_type):
            DomM = Matrix(elem.domain())
            return Matrix_type(DomM,[[elem]])
        elif isinstance(elem, int):
            DomM = Matrix(ZZ())
            return Matrix_type(DomM, [[elem]])
        elif isinstance(elem, float):
            DomM = Matrix(QQ(ZZ()))
            return Matrix_type(DomM, [[elem]])
        elif isinstance(elem,Symbol_type):
            DomM = Matrix(Poly(ZZ()))
            return Matrix_type(DomM, [[elem._to_monomial_()]])
        elif isinstance(elem,Base.Base):
            if elem._matrix_structure:
                return Matrix_type(elem,1,1)
            else:
                DomM = Matrix(elem)
                return Matrix_type(DomM,1,1)
        else:
            raise TypeError("Parameter elem not valid")
    elif cols == '':
        if isinstance(rows,list):
            if isinstance(elem,Base.Base):
                if elem._matrix_structure:
                    return Matrix_type(elem,rows)
                else:
                    DomM = Matrix(elem)
                    if len(rows) == 0:
                        return Matrix_type(DomM,[[]])
                    if isinstance(rows[0],list):
                        aux = []
                        for i in range(len(rows)):
                            aux2 = []
                            if len(rows[i])!= len(rows[0]):
                                raise TypeError("The list is not valid")
                            for j in range(len(rows[0])):
                                aux2.append(rows[i][j]>>elem)
                            aux.append(aux2)
                        return Matrix_type(DomM,aux)
                    else:
                        return Matrix_type(DomM,[rows])
            else:
                raise TypeError("Parameter elem not valid")
        elif isinstance(rows,int):
            if isinstance(elem,Base.Base):
                if elem._matrix_structure:
                    return Matrix_type(elem,rows)
                else:
                    DomM = Matrix(elem)
                    return Matrix_type(DomM,rows,rows)
            elif isinstance(elem,Vector_type):
                return Matrix_type(elem.domain(),elem._value*rows)
            elif isinstance(elem,Matrix_type):
                return Matrix_type(elem)
            elif isinstance(elem,Base_type.Base_type):
                DomM = Matrix(elem.domain())
                return Matrix_type(DomM,[[elem]*rows for i in range(rows)])
            elif isinstance(elem, int):
                DomM = Matrix(ZZ())
                return Matrix_type(DomM,[[elem>>ZZ()]*rows for i in range(rows)]) 
            elif isinstance(elem, float):
                DomM = Matrix(QQ(ZZ()))
                return Matrix_type(DomM,[[elem>>QQ(ZZ())]*rows for i in range(rows)]) 
            elif isinstance(elem, Symbol_type):
                DomM = Matrix(Poly(ZZ()))
                return Matrix_type(DomM,[[elem._to_monomial_()]*rows for i in range(rows)]) 
            else:
                raise TypeError("Wrong parameters")
        else:
            raise TypeError("Wrong parameters")
    elif s == '':
        if isinstance(cols, int):
            if isinstance(rows, int):
                if isinstance(elem, Base.Base):
                    if elem._matrix_structure:
                        return Matrix_type(elem,rows,cols)
                    else:
                        domM = Matrix(elem)
                        return Matrix_type(domM,rows,cols)
                elif isinstance(elem, Matrix_type):
                    return elem
                elif isinstance(elem, Vector_type):
                    raise TypeError("Wrong elem")
                elif isinstance(elem, Base_type.Base_type):
                    DomM = Matrix(elem.domain())
                    return Matrix_type(DomM,[[elem]*cols for i in range(rows)])
                elif isinstance(elem, int):
                    DomM = Matrix(ZZ())
                    return Matrix_type(DomM,[[elem>>ZZ()]*cols for i in range(rows)]) 
                elif isinstance(elem, float):
                    DomM = Matrix(QQ(ZZ()))
                    return Matrix_type(DomM,[[elem>>QQ(ZZ())]*cols for i in range(rows)]) 
                elif isinstance(elem, Symbol_type):
                    DomM = Matrix(Poly(ZZ()))
                    return Matrix_type(DomM,[[elem._to_monomial_()]*cols for i in range(rows)])
                elif isinstance(elem, list):
                    if len(elem)!=rows*cols:
                        raise TypeError("List and rows and cols does not match")
                    aux = []
                    for i in range(rows):
                        aux2 = []
                        for j in range(cols):
                            aux2.append(elem[i*cols + j])
                        aux.append(aux2)
                    return create_matrix(aux)
                else:
                    raise TypeError("Wrong parameters")
            else:
                raise TypeError("Wrong parameters")
    elif isinstance(s,int):
        if isinstance(cols, int):
            if isinstance(rows,list):
                if isinstance(elem, Base.Base):
                    aux = []
                    for i in range(cols):
                        aux2 = []
                        for j in range(s):
                            aux2.append(rows[i*s + j])
                        aux.append(aux2)
                    if elem._matrix_structure:
                        return Matrix_type(elem,aux)
                    else:
                        domM = Matrix(elem)
                        return Matrix_type(domM,aux)
        raise TypeError("Wrong parameters")
    else:
        raise TypeError("Wrong parameters")
        
        
        
def max_domain(*args):
    if len(args) == 1:
        if isinstance(args[0],list):
            aux = domain(args[0][0])
            for i in range(1,len(args[0])):
                try:
                    aux = union_domain(aux,domain(args[0][i]))
                except:
                    return 'Error: max domain does not exist'
            return aux
        else:
            return domain(args[0])
    else:
        aux = domain(args[0])
        for i in range(1,len(args)):
            try:
                aux = union_domain(aux,domain(args[i]))
            except:
                return 'Error: max domain does not exist'
        return aux
    
    
def union_domain(A,B):
    return (A.union_domain(B))[0]    
    

def append_vector(*args):
    if len(args)< 1:
        return 'append_vector error: Parameters missing'
    if isinstance(args[0], list):
        return append_vector(args[0])
    else:
        if len(args) == 1:
            return create_vector(args[0])
        else:
            LD = []
            aux = []
            for i in range(len(args)):
                if isinstance(args[i], int):
                    LD.append(ZZ())
                    aux.append(args[i])
                elif isinstance(args[i], float):
                    LD.append(QQ(ZZ()))
                    aux.append(args[i]>>QQ(ZZ()))
                elif isinstance(args[i], Symbol_type):
                    LD.append(Poly(ZZ()))
                    aux.append(args[i]._to_monomial_())
                else:
                    LD.append(args[i].domain())
                    aux.append(args[i])
            D = max_domain(LD)
            if D._matrix_structure:
                D = D
            else:
                D = Matrix(D)
            if isinstance(aux[0], Vector_type):
                S = aux[0]
            elif isinstance(aux[0], Matrix_type):
                S = aux[0]
            else:
                S = create_vector(D.subdomain(),[aux[0]])
            flagNull = False
            if isinstance(S,Matrix_type):
                if S._rows == 0 or S._cols == 0:
                    flagNull = True
            if isinstance(S, Vector_type):
                if len(S) == 0:
                    flagNull = True
            for i in range(1,len(aux)):
                if isinstance(aux[i], Matrix_type):
                    if aux[i]._rows == 0 or aux[i]._cols == 0:
                        continue
                if isinstance(aux[i], Vector_type):
                    if len(aux[i]) == 0:
                        continue
                if flagNull:
                    if isinstance(aux[i], Vector_type):
                        S = aux[i]
                    elif isinstance(aux[i], Matrix_type):
                        S = aux[i]
                    else:
                        S = create_vector(D.subdomain(),[aux[i]])
                    flagNull = False
                    continue
                S = D._ext_vector_(S,aux[i])
            return S

def stack(*args):
    if len(args)< 1:
        return 'stack error: Parameters missing'
    if (isinstance(args[0], list)) and (len(args) == 1):
        return stack(create_vector(args[0]))
    else:
        if len(args) == 1:
            return args[0]
        else:
            LD = []
            aux = []
            for i in range(len(args)):
                if isinstance(args[i], int):
                    LD.append(ZZ())
                    aux.append(args[i])
                elif isinstance(args[i], float):
                    LD.append(QQ(ZZ()))
                    aux.append(args[i]>>QQ(ZZ()))
                elif isinstance(args[i], Symbol_type):
                    LD.append(Poly(ZZ()))
                    aux.append(args[i]._to_monomial_())
                elif isinstance(args[i], list):
                    LD.append(max_domain(args[i]))
                    aux.append(create_vector(LD[-1],args[i]))
                else:
                    LD.append(args[i].domain())
                    aux.append(args[i])
            D = max_domain(LD)
            if D._matrix_structure:
                D = D
            else:
                D = Matrix(D)
            if isinstance(aux[0], Vector_type):
                S = aux[0]
            elif isinstance(aux[0], Matrix_type):
                S = aux[0]
            else:
                S = create_vector(D.subdomain(),[aux[0]])
            flagNull = False
            if isinstance(S,Matrix_type):
                if S._rows == 0 or S._cols == 0:
                    flagNull = True
            if isinstance(S, Vector_type):
                if len(S) == 0:
                    flagNull = True
            for i in range(1,len(aux)):
                if isinstance(aux[i], Matrix_type):
                    if aux[i]._rows == 0 or aux[i]._cols == 0:
                        continue
                if isinstance(aux[i], Vector_type):
                    if len(aux[i]) == 0:
                        continue
                if flagNull:
                    if isinstance(aux[i], Vector_type):
                        S = aux[i]
                    elif isinstance(aux[i], Matrix_type):
                        S = aux[i]
                    else:
                        S = create_vector(D.subdomain(),[aux[i]])
                    flagNull = False
                    continue
                S = D._ext_with_rows_(S,aux[i])
            return S
    
def splice(*args):
    if len(args)< 1:
        return 'splice error: Parameters missing'
    if (isinstance(args[0], list)) and (len(args) == 1):
        return splice(create_vector(args[0]))
    else:
        if len(args) == 1:
            return args[0]
        else:
            LD = []
            aux = []
            for i in range(len(args)):
                if isinstance(args[i], int):
                    LD.append(ZZ())
                    aux.append(args[i])
                elif isinstance(args[i], float):
                    LD.append(QQ(ZZ()))
                    aux.append(args[i]>>QQ(ZZ()))
                elif isinstance(args[i], Symbol_type):
                    LD.append(Poly(ZZ()))
                    aux.append(args[i]._to_monomial_())
                elif isinstance(args[i], list):
                    LD.append(max_domain(args[i]))
                    aux.append(create_vector(LD[-1],args[i]))
                else:
                    LD.append(args[i].domain())
                    aux.append(args[i])
            D = max_domain(LD)
            if D._matrix_structure:
                D = D
            else:
                D = Matrix(D)
            if isinstance(aux[0], Vector_type):
                S = aux[0]
            elif isinstance(aux[0], Matrix_type):
                S = aux[0]
            else:
                S = create_vector(D.subdomain(),[aux[0]])
            flagNull = False
            if isinstance(S,Matrix_type):
                if S._rows == 0 or S._cols == 0:
                    flagNull = True
            if isinstance(S, Vector_type):
                if len(S) == 0:
                    flagNull = True
            for i in range(1,len(aux)):
                if isinstance(aux[i], Matrix_type):
                    if aux[i]._rows == 0 or aux[i]._cols == 0:
                        continue
                if isinstance(aux[i], Vector_type):
                    if len(aux[i]) == 0:
                        continue
                if flagNull:
                    if isinstance(aux[i], Vector_type):
                        S = aux[i]
                    elif isinstance(aux[i], Matrix_type):
                        S = aux[i]
                    else:
                        S = create_vector(D.subdomain(),[aux[i]])
                    flagNull = False
                    continue
                S = D._ext_with_columns_(S,aux[i])
            return S

def ncols(a):
    return a.ncols()

def nrows(a):
    return a.nrows()

def get_row(s,H):
    return H.find_index_row(s)

def get_col(s,H):
    return H.find_index_col(s)

row_index = get_row # H[s,:]
col_index = get_col
    
def extension(K,poly,symbol = 'x', name = ''):
    E = Fq(K,poly,symbol,name)
    x = E.variable();
    return [E,x]

def variable(f):
    if isinstance(f,Poly_type):
        return f.variable()
    if isinstance(f,Monomial_type):
        return f.variable()
    if isinstance(f,Symbol_type):
        return f.cp()
    if isinstance(f,int):
        return None
    if isinstance(f,float):
        return None
    return f._structure.variable()

def gcd(f0,f1): 
    if isinstance(f0,int) and isinstance(f1,int):
        return ZZ().gcd(f0,f1)
    if isinstance(f0,int):
        f00 = f0>>f1.domain()
    else:
        f00 = f0
    if isinstance(f1,int):
        f11 = f1>>f0.domain()
    else:
        f11 = f1
    if f11.domain() == f11.domain():
        return f00.domain().gcd(f00,f11)
    if f11.domain().is_sub_domain(f00.domain()):
        return f11.domain().gcd(f00>>f11.domain(),f11)
    if f00.domain().is_sub_domain(f11.domain()):
        return f11.domain().gcd(f00,f11>>f00.domain())
    D = max_domain(f00.domain(),f11.domain())
    return D.gcd(f00>>D,f11>>D)
    
def is_square_free(f):
    f1 = f.der()
    d = gcd(f,f1)
    if d == 1: return True
    return False


    

#def is_irreducible(f, K=''):
#    # error if K is not a field and f not univariate polynomial
#    # if coefficient_ring(f) != K: f = reduce(f,K)
#    # (in case f is given with integer coefficients)
#
#    if K == '': 
#        K = f.domain().subdomain()
#    if (f.domain().subdomain().is_sub_domain(K)) and f.domain().subdomain()!=K:
#        for i in range(len(f._value)):
#            if not K.is_element(f._value[i].transform_to_min_domain()):
#                raise TypeError("This polynomial is not defined in this domain")
#        aux = []
#        aind = []
#        for i in range(len(f._value)):
#            aux.append(pull(f._value[i],K))
#            aind.append(f._index[i])
#        Dom = Poly(K)
#        f = Poly_type(aux,Dom,aind)
#    if not K.is_sub_domain(f._structure.subdomain()):
#        raise TypeError("This polynomial cannot be computed in this domain")
#    if K.characteristic() == 0:
#        faux = f
#        while True:
#            lc = faux.constant_coeff()>>ZZ()
#            if isinstance(lc,int):
#                break
#            faux = lc
#        
#        D = divisors(lc)
#        for d in D:
#            if f.evaluate(d) == 0:
#                return False
#        return True
#    
#    n =  f.degree()#n = degree(f)
#    D = f.domain()
#    if n == 0: return False
#    if n == 1: return True
#    X = D.variable()  #X = variable(f)
#    q = K.cardinal()#q = cardinal(K)
#    #r = D.variable_d(q**n)#(X**(q**n))
#
#    #r = md_power(X,q**n,f)
#    r = f.domain()._md_double_pow_(X,q,n,f)
#    #r = md_power(xv,q**(n-1),f)#r%f#remainder(X**(q**n),f)  # later we will optimize this step
#
#    if r != X: return False
#    P = set(prime_factors(n))
#
#    for p in P:
#        #r = D.variable_d(q**(n//p))
#
#        r = md_power(X,q**(n//p),f)
##        r = (X**(q**(n//p))) % f #remainder(X**(q**(n//p))-X, f) # and also this one
#
#        a1 = D.gcd(r-X,f)
#
#        if a1 != 1: 
#            return False
#    return True

# cal definir is_finite_field(K)
def is_field(K): 
    if isinstance(K,Base.Base):
        return K.is_field()
    return False




def hohner(A,x):
    if list(A) == []: return 0
    a = A[0]
    for c in A[1:]:
        a = a*x+c
    return a
#
def polynomial(A,x): 
    if list(A) == []: return 0
    B=[A[j] for j in range(len(A)-1,-1,-1)]
    b = B[0]
    for c in B[1:]:
        b = b*x+c
    return b
    


def domain(a):
    if isinstance(a,int):
        return ZZ()
    elif isinstance(a,float):
        return QQ(ZZ())
    elif isinstance(a,Base.Base):
        return a.cp()
    elif isinstance(a,Symbol_type):
        return Poly(ZZ())
    else:
        return a.domain()


def generator(A): return A.element([1,0])
def degree(A, variable = None): 
    if isinstance(A, Poly_type):
        return A.degree(variable)
    if isinstance(A,Monomial_type):
        return A.degree(variable)
    if isinstance(A,Symbol_type):
        return A.degree(variable)
    return A.degree()
#def field(x): return x.min_domain()
def is_sub_domain(A,B):
    return B.is_sub_domain(A)

def pull(x,K = ''):
    if K == '':
        return x.transform_to_min_domain()
    if isinstance(x,Matrix_type):
        if x.ncols() != 1 or x.nrows() != 1:
            if K == x.domain():
                return x.cp()
            if x.domain().is_sub_domain(K):
                aux = []
                D = Matrix(K)
                for i in range(x.nrows()):
                    aux2 = []
                    for j in range(x.ncols()):
                        aux2.append(pull(x[i][j],K))
                    aux.append(aux2)
                return Matrix_type(D,aux)
            if K.is_sub_domain(x.domain().subdomain()):
                aux = []
                D = Matrix(K)
                for i in range(x.nrows()):
                    aux2 = []
                    for j in range(x.ncols()):
                        aux2.append(x[i][j]>>K)
                    aux.append(aux2)
                return Matrix_type(D,aux)
    if isinstance(x,Vector_type):
        if len(x)!= 1:
            if K == x.domain():
                return x.cp()
            if x.domain().is_sub_domain(K):
                aux = []
                D = Matrix(K)
                for i in range(len(x)):
                    aux.append(pull(x[i],K))
                if x._direction:
                    return Vector_type(D,aux)
                else:
                    return Vector_type(D,aux,'c')
            if K.is_sub_domain(x.domain().subdomain()):
                D = Matrix(K)
                aux = []
                for i in range(len(x)):
                    aux.append(x[i]>>K)
                if x._direction:
                    return Vector_type(D,aux)
                else:
                    return Vector_type(D,aux,'c')
    if isinstance(x,Symbol_type):
        return x.cp()
    if isinstance(x,Monomial_type):
        if x.mdegree() == 0:
            return pull(x._value,K)
        if K == x.domain():
            return x.cp()
        if x.domain().is_sub_domain(K):
            D = Poly(K)
            if x._varOrder== None:
                return Monomial_type(pull(x._value,K),x._variables[:],x._index[:],D)
            else:
                return Monomial_type(pull(x._value,K),x._variables[:],x._index[:],D,x._varOrder[:])
        if K.is_sub_domain(x.domain().subdomain()):
            D = Poly(K)
            if x._varOrder== None:
                return Monomial_type(x._value >> K,x._variables[:],x._index[:],D)
            else:
                return Monomial_type(x._value>>K,x._variables[:],x._index[:],D,x._varOrder[:])
        if x.domain().characteristic() == 0 and K.characteristic() != 0:
            y = x.change_characteristic(K.characteristic())
            return pull(y,K)
    if isinstance(x,Poly_type):
        if K == x.domain():
            return x.cp()
        if x.domain().is_sub_domain(K):
            aux = []
            D = Poly(K)
            md = 0
            if x._varOrder== None:
                for i in range(len(x._value)):
                    s = pull(x._value[i]._value,K)
                    if s != 0:
                        aux.append(Monomial_type(s,x._variables[:],x._index[:],D))
                        if md < aux[-1]._totIndex:
                            md = aux[-1]._totIndex
                if len(aux) == 0:
                    return Monomial_type(0>>K,x._variables[:],x._index[:],D)
                elif len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,x._variables[:],D,aux[0]._degree,md)
            else:
                for i in range(len(x._value)):
                    s = pull(x._value[i]._value,K)
                    if s != 0:
                        aux.append(Monomial_type(s,x._variables[:],x._index[:],D,x._varOrder[:]))
                        if md < aux[-1]._totIndex:
                            md = aux[-1]._totIndex
                if len(aux) == 0:
                    return Monomial_type(0>>K,x._variables[:],x._index[:],D)
                elif len(aux) == 1:
                    return aux[0]
                else:
                    vord = D._get_own_order_monomials_(aux)
                    return Poly_type(aux,x._variables[:],D,aux[0]._degree,md,vord,x._varOrder[:])
        if K.is_sub_domain(x.domain().subdomain()):
            aux = []
            D = Poly(K)
            md = 0
            if x._varOrder== None:
                for i in range(len(x._value)):
                    s = x._value[i]._value >> K
                    if s != 0:
                        aux.append(Monomial_type(s,x._variables[:],x._value[i]._index[:],D))
                        if md < aux[-1]._totIndex:
                            md = aux[-1]._totIndex
                if len(aux) == 0:
                    return Monomial_type(0>>K,x._variables[:],[0]*len(x._variables),D)
                elif len(aux) == 1:
                    return aux[0]
                else:
                    return Poly_type(aux,x._variables[:],D,aux[0]._degree,md)
            else:
                for i in range(len(x._value)):
                    s = x._value[i]._value >> K
                    if s != 0:
                        aux.append(Monomial_type(s,x._variables[:],x._value[i]._index[:],D,x._varOrder[:]))
                        if md < aux[-1]._totIndex:
                            md = aux[-1]._totIndex
                if len(aux) == 0:
                    return Monomial_type(0>>K,x._variables[:],[0]*len(x._variables),D)
                elif len(aux) == 1:
                    return aux[0]
                else:
                    vord = D._get_own_order_monomials_(aux)
                    return Poly_type(aux,x._variables[:],D,aux[0]._degree,md,vord,x._varOrder[:])
        if x.domain().characteristic() == 0 and K.characteristic() != 0:
            y = x.change_characteristic(K.characteristic())
            return pull(y,K)
            
    if isinstance(x,Fq_type):
        if len(x._index) > 1 or x._index[0] != 0:
            if K == x.domain():
                return Fq_type(x._value,x.domain(),x._index)
            if K.is_sub_domain(x.domain()):
                return x>>K
            if x.domain().is_sub_domain(K):
                KK = field(x)
                if K.is_sub_domain(KK):
                    return x.transform_to_min_domain()>>K
                if KK == K:
                    return x.transform_to_min_domain()
                return x
        else:
            return pull(x._value[0],K)
    if K.characteristic() == 0:
        return x>>K
    if x.domain().characteristic() == K.characteristic():
        if is_Zn(K):
            if field(x) == K:
                return x.transform_to_min_domain()
            return x
        Q = field(x)
        if K.is_sub_domain(Q):
            return (x.transform_to_min_domain())>>K #Falta mirar nomenclatura
        raise TypeError("This element cannot be projected to this subdomain")
    if is_Zn(K):
        return x.change_characteristic(K.characteristic())
    raise TypeError("This element cannot be projected to this space")
# s'ha de definir extension_base(A), generator(A), degree(A)


def get_irreducible_polynomial(K,r,variable='X'):
    [P,x] = univariate_polynomial_ring(K,variable)
    if r==1: return x
    f = x**r
    ST = cardinal(K)**r-1
    while True:
        k = rd_int(0,ST)
        g = f + pick_element(k,P,variable)
        if is_irreducible(g,K):
            return g
        #2k=k+1



'''
factors a monic f over K with n>0 into "distinct degree factors"
'''
def distinct_degree_factoring(f,K):
    #return (f,K)
    G = []
    x = variable(f)
    p = x
    #return x
    h = f
    q = cardinal(K)
    r = 0
    while True:
        r += 1
        p = md_power(p,q,f)
        g = gcd(p-x,h)
        if g != 1:        
            G += [(g,r)]
        h = h // g
        if degree(h) < 2*r+2:
            return G+[(h,r+1)]
#
DDF = distinct_degree_factoring

def equal_degree_splitting(f,r):
    #show('EDS')
    P = domain(f)
    K = K_(f)
    p = characteristic(K)
    q = cardinal(K)
    k = ifactor(q)[p]
    n = degree(f)
    car2 = False
    if p==2:
        car2 = True
    while True:
        while True:
            h = rd_pick(P,n-1, f.variables_name())
            if degree(h)>0: break 
        g = gcd(h,f)
        if degree(g)>0: return g
        if car2: h = md_trace(h,k*r,f)
        else: h = md_power(h,(q**r-1)//2,f)
        g = gcd(h-1,f)
        if 0<degree(g)<degree(f): 
            return g
EDS = equal_degree_splitting

def equal_degree_factoring(f,r):
    P = domain(f)
    K = K_(f)
    q = cardinal(K)
    if is_irreducible(f,K): return [f]
    g = EDS(f,r);    
    h = f//g
    return EDF(g,r)+EDF(h,r)
#
EDF = equal_degree_factoring 

def factor(f,K=Zn(2)):
    if K_(f)==ZZ(): 
        a = f.to_list()
        [_,x] = univariate_polynomial_ring(K,'x')
        ff = hohner(a,x)
        ff = ff >> K
    else:
        ff = f
    if K.is_sub_domain(K_(ff)):
        o = 0>>K
        ff = ff + o
    if K_(ff) != K: return 'Error: wrong relative field'
    x = variable(ff)
    q = cardinal(K_(ff))
    P = []
    h = x
    v = ff
    r = 0
    while True:
        h = md_power(h,q,ff) 
        g = gcd(h-x,v)
        r += 1
        if g!= 1:
            G = equal_degree_factoring(g,r)
            for u in G:
                e = 0
                while v % u == 0:
                   v = v//u
                   e += 1
                P += [(u,e)]
        if v==1: return P
        
def rd_pick(P,k,symb = ['x']): 
    if P._poly_structure:
        return P._rd(k,symb)
    return P._rd(k)  

def order(a,b=''):
    if isinstance(a,int):
        if isinstance(b,int):
            if gcd(a,b) > 1:
                return "a is not invertible mod b"
            ZD = Zn(b)
            a = a>>ZD
            return a.order()
        if isinstance(b,Base.Base):
            a = a>>b
            return a.order()
    if isinstance(a,Base_type.Base_type):
        if b == '':
            return a.order()
        if isinstance(b,Base.Base):
            a = a*(1>>b)
            return a.order()
    if isinstance(a,float):
        a = a>>QQ(ZZ())
        return order(a,b)
    if isinstance(a,Symbol_type):
        x = a._to_monomial_()
        return order(x,b)
    raise TypeError("Error parameters")  

def geometric_series(x,n,s0=1):
    if s0==1: return create_vector([x**k for k in range(n)])
    return create_vector([s0*x**k for k in range(n)])

def merge(A,B):
    L = A
    A1 = [a[0] for a in A]
    for b in B:
        if A1.count(b[0]) == 0: L += [b]
        else: 
            j = A1.index(b[0])
            L[j] = (L[j][0],L[j][1]+b[1])
    return L   

def find_roots(f,K=''):
    if isinstance(f,Symbol_type):
        return [0]
    if  K=='': K = f.domain().subdomain()
    f = pull(f>>ZZ(),K)
    x = variable(f)
    #h = variable_d(q,f.domain())
    h = md_power(x,cardinal(K),f)
    g = gcd(h-x,f)
    r = degree(g)
    if r == 0: return []
    U = EDF(g,1)
    R = [-u.constant_coeff() for u in U]
    return R
#
roots = find_roots
#def roots(f,F):
#    f = f+(0>>F)
#    X = variable(f)
#    R = EDF(f,1)
#    return [X-r for r in R] 

       
# Versió de pfactor. Es podria fer una versió
# que agafés el millor de les dues. Si hagués
# pensat en mirar-la, en lloc de fer-la de nou,
# potser ja estaria...
# n a positive integer
def cyclotomic_polynomial(n, K= ZZ()):
    P = set(prime_factors(n))
    (Z,X) = univariate_polynomial_ring(K,'X')
    #(Z,X) = univariate_polynomial_ring(Zn(3),'X')
    Q = X-1; e=n
    for p in P: 
        Q = (Q.evaluate(X**p))//Q 
        e = e//p 
    return Q.evaluate(X**e)

def frobenius(K): 
    if isinstance(K,int): return lambda x: x**K
    if K.is_field(): return frobenius(cardinal(K))
#
F_=frobenius

# Conjugates of an element x in L over subfield K
def conjugates(x,K=''):
    if K == '':
        K = prime_field(K_(x))
    q = cardinal(K)
    X = [x]
    y = x**q
    while y != x:
        X += [y]
        y = y**q
    return X


def conjugate_classes(E,K=Zn(2)):
    X = [k>>K for k in range(cardinal(K))]
    C = [[x] for x in X]
    for j in range(cardinal(K),cardinal(E)):
        x = pick(j,E)
        if X.count(x)==0:
            c = conjugates(x,K)
            X += c
            C += [c]
    return C
    
def minimal_polynomial(x,K=Zn(2),T='X'):
    C = conjugates(x,K)
    E = x.domain()
    [P,X]=univariate_polynomial_ring(E,T)
    m=1
    for y in C:
        m *= (X-y)
    return m

#def trace(x,K=Zn(2)): return sum(conjugates(x,K))
def nm(x,K=Zn(2)): 
    Y = conjugates(x,K)
    N = 1
    for y in Y:
        N *= y
    return N

def Set(K): return [pick(j,K) for j in range(cardinal(K))]
X_ = Set

def period(f):
    if isinstance(f,Poly_type):
        m = degree(f)
        x = variable(f)
        q = cardinal(K_(f))
        D = divisors(q**m-1)
        for d in D:
            if md_power(x,d,f) == 1: return d
            #if (x**d-1)%f == 0: break
        return 'Error'
    if isinstance(f,Monomial_type):
        m = degree(f)
        x = variable(f)
        q = cardinal(K_(f))
        D = divisors(q**m-1)
        for d in D:
            if md_power(x,d,f) == 1: return d
            #if (x**d-1)%f == 0: break
        return 'Error'
    if isinstance(f,Symbol_type):
        return period(f._to_monomial_())
    if isinstance(f,Fq_type):
        return f.order()
    if isinstance(f,Zn_type):
        return f.order()
exponent = period
#
def is_primitive(f):
    e = period(f)
    m = degree(f)
    q = cardinal(K_(f))
    if e == q**m-1: return True
    return False

def primitive_root(K):
    q = cardinal(K); n = q-1
    for j in range(2,q):
        if order(element(j,K))==n:
            return element(j,K)
    return "primitive_root error: root not found"

#20190201
def md_trace(h,m,f):
    T = h
    for j in range(1,m):
        h = h**2 % f
        T += h % f
    return T    

def product(X):
    P=1
    for x in X: P *= x
    return P

def cyclotomic_splitting(n,K=Zn(2)):
    q = cardinal(K)
    m = order(q,n)
    f=get_irreducible_polynomial(K,m,'x')
    x = variable(f)
    [F,t] = extension(K,f,'t','F')
    r = (q**m-1)//n
    a = primitive_root(F)
    w = a**r
    show('order(w) =',order(w))
    C = cyclotomic_classes(n,q)
    return [product([(x-w**j) for j in c]) for c in C]


# n-th cyclotomic polynomial Q_n
# It is a different version than the one in pfactor
# called cyclotomic_polynomial
def cyclotomic_polynomial_moebius(n):
    [_,X]=univariate_polynomial_ring(ZZ())
    if n==1: return X**0
    mu = mu_moebius
    D = divisors(n)
    Dp = [d for d in D if mu(d)==1]
    Dm = [d for d in D if mu(d)==-1]
    P = 1
    for d in Dp:
        P *= X**(n//d)-1
    for d in Dm:
        P = P/(X**(n//d)-1)
    return P

# Aquest funció fa el mateix que cyclotomic_polynomial
# de pfactor. Té l'avantage que li pots passar el
# nom de la variable. La seva forma és independent
# del cos --però la seva factorització no, és clar.
def fast_cyclotomic(n,x='x'):
    if n==0: return 1
    F = ifactor(n)
    P = list(F.keys())
    m = p = P[0]
    [_,x] = univariate_polynomial_ring(ZZ(),x)
    f = sum([x**j for j in range(p)])
    for p in P[1:]:
        m *= p
        f = evaluate(f,x**p)/f
    return evaluate(f,x**(n//m))




