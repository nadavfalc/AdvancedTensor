'''
Package CC: Computational codes
'''

import PyM.PyWIT.PyECC.PyECC_globals as PGl
from PyM.PyWIT.PyECC.Alternant_Code import *
from PyM.PyWIT.PyECC.GX import *
from PyM.PyWIT.PyECC.PyECC_sets import *
from PyM.PyWIT.PyECC.Poly import _order_symbols_
#


## Global symbols
Q_=QQ(ZZ(),'Q')
Z_ = ZZ()
vector = row_vector = create_vector
vector_append = append_vector
matrix = create_matrix

def show(*L,nl=1):
    print(*L,nl*'\n')


## Domain/Type functions

def belongs(x,X): return X.is_element(x)

EPS = PGl.EPS
def notnil(x): return abs(x)>=EPS
def nil(x): return abs(x)<EPS
def nonneg(x): return (notnil(x) and x>0) or nil(x)

def is_number(x):
    if isinstance(x, (int, float, complex,QQ_type)) and not isinstance(x, bool): 
        return True
    else: return False
#
def is_real(x):
    if is_number(x) and not isinstance(x,complex): return True
    else: return False

is_pair = ispair

def base(B): return domain(B).subdomain()

#def max_domain(X):
#    if isinstance(X[0],int):
#        d = ZZ()
#    else:
#        d = X[0].domain()
#    for x in X[1:]:
#        if isinstance(x,int): continue
#        d1 = x.domain()
#        if d.is_sub_domain(d1): continue
#        if d1.is_sub_domain(d): 
#            d = d1
#            continue
#        return 'Error: max domain does not exist'
#    return d  


def find_zero(f,a,b,tol):
    if f(a)*f(b) < 0:
        ini = a
        fin = b
        vant = f(a)
        if abs(vant) < tol:
            return a
        while True:
            mig = (ini+fin)/2
            if abs(f(mig)) < tol:
                pzero = mig
                break
            if f(mig)*f(ini) > 0:
                ini = mig
                vant = f(mig)
            else:
                fin = mig
        return pzero        
    else:
        delta = 1000
        vant = f(a)
        if abs(vant) < tol:
            return a
        L = [vant]
        V = [a]
        for i in range(delta-1):
            vact = f(a+(i+1)*(b-a)/delta)
            if vact*vant < 0:
                return find_zero(f,a+i*(b-a)/delta,a+(i+1)*(b-a)/delta,tol)
            elif abs(vact) < tol:
                return a+(i+1)*(b-a)/delta
            L.append(vact)
            V.append(a +(i+1)*(b-a)/delta )
        if vact < 0:
            for i in range(1,len(L)-1):
                if L[i] > L[i-1] and L[i]>L[i+1]:
                    fmig = L[i]
                    fesq = L[i-1]
                    fdre = L[i+1]
                    mig = V[i]
                    esq = V[i-1]
                    dre = V[i+1]
                    while(abs(fmig-fesq)>tol/10):
                        v1 = (mig+esq)/2
                        f1 = f(v1)
                        v2 = (dre+mig)/2
                        f2 = f(v2)
                        if (f1<fmig)and(fmig>f2):
                            fesq = f1
                            fdre = f2
                            esq = v1
                            dre = v2
                        elif (f1 < fmig)and(fmig < f2):
                            fesq = fmig
                            fmig = f2
                            esq = mig
                            mig = v2
                        else:
                            fdre = fmig
                            fmig = f1
                            dre = mig
                            mig = v1
                    if (abs(fmig)<10*tol):
                        return mig
        else:
            for i in range(1,len(L)-1):
                if L[i] < L[i-1] and L[i]<L[i+1]:
                    fmig = L[i]
                    fesq = L[i-1]
                    fdre = L[i+1]
                    mig = V[i]
                    esq = V[i-1]
                    dre = V[i+1]
                    while(abs(fmig-fesq)>tol/10):
                        v1 = (mig+esq)/2
                        f1 = f(v1)
                        v2 = (dre+mig)/2
                        f2 = f(v2)
                        if (f1>fmig)and(fmig<f2):
                            fesq = f1
                            fdre = f2
                            esq = v1
                            dre = v2
                        elif (f1 > fmig)and(fmig > f2):
                            fesq = fmig
                            fmig = f2
                            esq = mig
                            mig = v2
                        else:
                            fdre = fmig
                            fmig = f1
                            dre = mig
                            mig = v1
                    if (abs(fmig)<10*tol):
                        return mig
        return 'find_zero: Not zero found'


def is_prime_field(A):
    if isinstance(A,Base.Base):
        if is_field(A):
            if is_Zn(A):
                return True
    return False

def is_finite_field(A):
    if isinstance(A,Base.Base):
        if is_field(A):
            if characteristic(A)!=0:
                return True
    return False

def is_ring(A):
    if isinstance(A, Base.Base):
        return A._ring
    return False


def lift(a):  
    x = pull(a)
    
    A = x.domain()
    
    if A._poly_structure:
        if isinstance(a,Symbol_type):
            return a.cp()
        if isinstance(a,Poly_type):
            aux = []
            md = 0
            D = domain(lift(a._value[0]._value))
            PD = Poly(D)
            for val in x._value:
                M = lift(val)
                if M._value != 0:
                    aux.append(M)
                    if aux[-1]._totIndex > md:
                        md = aux[-1]._totIndex
            if len(aux) == 0:
                return Monomial_type(0 >> D, a._variables[:], [0]*len(a._variables),PD)
            elif len(aux) == 1:
                return aux[0]
            else:
                if a._varOrder == None:
                    return Poly_type(aux,a._variables[:],PD,aux[0]._degree,md)
                else:
                    vord = PD._get_own_order_monomials_(aux)
                    return Poly_type(aux,a._variables[:],PD,aux[0]._degree,md,vord,a._varOrder[:])
        elif isinstance(a,Monomial_type):
            if isinstance(a._value,int):
                return a.cp()
            else:
                D = domain(lift(a._value[0]._value))
                PD = Poly(D)
                if a._varOrder == None:
                    s = lift(a._value)
                    if s != 0:
                        return Monomial_type(s,a._variables[:],a._index[:],PD)
                    else:
                        return Monomial_type(0 >> D, a._variables[:], [0]*len(a._variables),PD)
                else:
                    s = lift(a._value)
                    if s != 0:
                        return Monomial_type(s,a._variables[:],a._index[:],PD,a._varOrder[:])
                    else:
                        return Monomial_type(0 >> D, a._variables[:], [0]*len(a._variables),PD)
    elif A._power_structure:
        aux = []
        for val in x._value:
            aux.append(lift(val))
        D = max_domain(aux)
        DD = Power_Series(D,A._symbol)
        return Power_Series_type(aux,DD,x._index)
    elif A._matrix_structure:
        if isinstance(x,Vector_type):
            aux = []
            for val in x._value:
                aux.append(lift(val))
            D = max_domain(aux)
            DD = Matrix(D)
            if x._direction:
                return DD.Vector_type(aux,'r')
            else:
                return DD.Vector_type(aux,'c')
        elif isinstance(x,Matrix_type):
            aux = []
            aux3 = []
            for val in x._value:
                aux2 = []
                for val2 in val:
                    ll = lift(val2)
                    aux2.append(ll)
                    aux3.append(ll)
                aux.append(aux2)
            D = max_domain(aux3)
            DD = Matrix(D)
            return DD.element(aux)
    elif A._fraction_structure:
        aux = []
        for val in x._value:
            aux.append(lift(val))
        D = max_domain(aux)
        DD = QQ(D)
        return QQ_type(DD,aux[0],aux[1])
    elif is_Zn(A):
        return x._value
    elif isinstance(A,Fq):
        D = A._poly.domain()
        v = A._poly._variables
        return Poly_type(x._value[:],D,v,[[inv] for inv in x._index],1,goodValues = 1)
    else:
        raise NotImplementedError
        
    return x

def push(a,D): #canviar el nom a push
    if isinstance(a,int):
        return a>>D
    if isinstance(a,float):
        aux = a>>QQ(ZZ())
        return aux._push_to_union(D)
    if isinstance(a,Base_type.Base_type):
        return a._push_to_union(D)
    if isinstance(a,Symbol_type):
        return a._push_to_union(D)
    if isinstance(a,list):
        aux = []
        for i in a:
            aux.append(push(i,D))
        return aux
#




## Numerical and set functions
def continuous_fraction(x,n):
    X = x
    F = []
    for i in range(n):
        f = floor(X)
        F += [f]
        if f==X:
            return F + [0]*(n-len(F));
        X = 1 /(X-f)
    return F

def continuous_fraction_value(F):
    r = F[len(F)-1]
    u = 1>>Q_
    for f in range(len(F)-2,-1,-1):
        r = F[f] + u/r
    return r>>Q_

# multiplication up to b bits
def mult(n,k,b):
    Z = Zn(2**b)
    n = n>>Z; k = k>>Z
    return lift(n*k)

# division up to b bits, k odd
def quot(n,k,b):    
    Z = Zn(2**b)
    n = n>>Z; k = k>>Z
    return lift(n/k)

# power up to b bits
def bpow(n,k,b):    
    return power(n,k,2**b)


def is_prime_power(n): 
    l = len(ifactor(n))
    return l==1

def next_q(n):
    # Assume x is a positive integer
    while not is_prime_power(n):
        n+=1
    return n

def is_square_free_n(n):
    return mu_moebius(n) != 0
    
# fast check whether n == x**k
def power_check(n,x,k):       
    f = blen(2*n)
    if x==1: 
        if n==1: return True
        else: return False
    P = x**k
    b = 1
    while True:
        m = 2**b
        r = power(x,k,m)
        if n % m != r: return False
        if b >= f:
            if n == P: return True
            else: return False
        b = min(2*b,f)


def nroot(n,k,b):
    if k==2: return nsqroot(n,b)
    if even(n): return "n has to be odd"
    if even(k): return "k has to be either odd or 2"
    if b==1: return 1
    B = []
    while b>1:
        B = [b]+B
        if b%2: b+=1
        b = b//2
    r = 1
    for b in B:
        r0 = mult(r,k+1,b)
        r1 = mult(n,bpow(r,k+1,b),b)
        r = quot(r0-r1,k,b)
    return r

'''
def tnroot(d):
    n = ri(d)
    if n%2==0: n+=1
    k = ri(3)
    if k%2==0: k+=1
    nu = rd_int(1,3)
    b = ri(nu)
    r = nroot(n,k,b)
    show(r,k,n,b)
    if ((r**k*n)%(2**b))!=1: return(r,n,b)
    return 1
'''    

def nsqroot(n,b):
    if even(n): return "n has to be odd"
    if b==1: 
        if (n%4)==1: return 1
        else: return 0
    if b==2: 
        if (n%8)==1: return 1
        else: return 0
    B = []
    while b>2:
        B = [b]+B
        if b%2: b = (b+1)//2
        else: b = 1+b//2
    if (n%8)!=1: return 0
    r = 1
    for b in B:
        r0 = mult(r,3,b+1)
        r1 = mult(n,bpow(r,3,b+1),b+1)
        r = ((r0-r1)//2)%(2**b)
        if r==0: return 0
    return r

'''    
def tsqroot(d):
    n = ri(d)
    if n%2==0: n+=1
    nu = rd_int(1,2)
    b = ri(nu)
    r = nsqroot(n,b)
    show(r,n,b)
    if r==0: return 0
    if ((r**2*n)%(2**(b+1)))!=1: return(r,n,b)
    return 1
'''

def inverse(n,m):
    if gcd(n,m)>1: return "n has no inverse mod m"
    return lift(1/(n>>Zn(m)))
   
'''    
def tinverse(d):
    n = ri(d)
    m = ri(d)
    if n%2==0 and m%2==0: m+=1
    f = gcd(n,m)
    if f>1: 
        show(n,m,f)
        return "Bad choice, try again"
    i = inverse(n,m)
    show(i,n,m)
    if i*n % m != 1: show((i,n,m))
    return 1
'''

def qceiling(f,k):
    (q,r)=(f//k,f%k)
    if r==0: b=q
    else: b=q+1
    return b

def is_power(n,k): 
    if n%2==0: return "n has to be odd"
    if k%2==0 and k>2: return "k has to be odd or 2"
    f = blen(2*n)
    b = qceiling(f,k)
    #show(f,b)
    y = inverse(n,2**(b+1))
    if k==2: 
        r = nsqroot(y,b)
        #show("r =",r)
        if r==0: return 0
    else: 
        r = nroot(y,k,b)
    if power_check(n,r,k):
        return r
    if k==2 and power_check(n,2**b-r,k):
        return 2**b-r
    return 0


def irr(q,t):
    N = 0
    D = divisors(t)
    for d in D:
        N += mu_moebius(d)*q**(t//d)
    return N//t
    
    
def hd(x,y):
    d = 0
    l = len(x)
    if len(y) != l: return 'hd Error: vectors do not have the same length'
    for j in range(l):
        if x[j] != y[j]:
            d += 1
    return d

def cyclic_shift(x): 
    if isinstance(x,Vector_type):
        x = list(x)
    return [x[-1]] + x[:-1]

def cyclic_shifts(x):
    X = [x]
    if isinstance(x,Vector_type):
        x = list(x)
    y = cyclic_shift(x)
    while y != x:
        X = X + [y]
        y = cyclic_shift(y)
    return list(X)

'''
def tis_power():
    k = 401
    x = ri(4)
    while x%2==0:
        x = ri(4)
    n = x**k+46
    show (x,k,n)
    return is_power(n,k)
'''    


def even(n):
    if isinstance(n,QQ_type):
        if n.denominator() == 1:
            return even(n.numerator())
        else:
            return "odd error:Not integer number"
    if n%2: return False
    return True
def odd(n):
    if isinstance(n,QQ_type):
        if n.denominator() == 1:
            return odd(n.numerator())
        else:
            return "odd error:Not integer number"
    if n%2: return True
    return False
def primes_less_than(f):
    P = []
    for p in range(3,f,2):
        if is_prime(p): P += [p]
    if f>2: P = [2]+P
    return P
    

def is_perfect_power(n):
    if n<2: return 'n has to be at least 2'
    if even(n): return 'n has to be odd'
    f = blen(2*n)
    b = qceiling(f,2)
    #show(f,b)
    y = nroot(n,1,b+1)
    P = primes_less_than(f)
    for p in P:
        x = is_power(n,p)
        if x>0: return (x,p)
    return (n,1)

# Computes K
def K_(f): 
    if isinstance(f,Base_type.Base_type):
        return f.K_()
    if isinstance(f,Base.Base):
        return f.K_()
    if isinstance(f,int):
        return ZZ()
    if isinstance(f,float):
        return QQ(ZZ())
    if isinstance(f,Symbol_type):
        return f.K_()
    if isinstance(f,Alternant_Code):
        return f.K_()
    if isinstance(f,list):
        return K_(vector(f))
    return "Type not found"

        
#

## Vector and list functions
# col_vector
def col_vector(K,X): return vector(K,X,'c')

    

# reverse a list
def reverse(x): 
    y = [x[j] for j in range(len(x)-1,-1,-1)]
    if isinstance(x,Vector_type): return push(vector(y),domain(x))
    return y

# support and histogram of a vector

def support(x): return [j for j in range(len(x)) if x[j]!=0]
def histogram(e):
    I = support(e)
    if I == []: return []
    v = [e[i] for i in I]
    return I, vector(v)
pattern = histogram

# weight of a vector
def wt(x):
    w = 0
    for t in x:
        if t != 0:
            w +=1
    return w

# k-th coefficient of convolution of a and b (lists or vectors)
def convolution(a,b,k=''):
    if k=='': 
        C = [convolution(a,b,k) for k in range(len(a)+len(b)-1)]
        if isinstance(a,Vector_type) or isinstance(b,Vector_type):
            return vector(C)
        return C
    r = len(a); s = len(b)
    if r==0 or s==0: return 0
    m = r+s-1
    if k>m: return 0
    ck = 0
    for j in range(max(0,k-s+1),min(r-1,k)+1):
        ck += a[j]*b[k-j]
    return ck

def eps(m,n,j='',k=''): 
    if j=='':
        j=n
        if j>=0 and j<m: return j*[0]+[1]+(m-j-1)*[0]
        return m*[0]
    else:
        A = create_matrix(0,m,n)
        if j<0 or j>=m:
            return A
        if k<0 or k>=n:
            return A
        A[j,k]=1
        return A


## Finite field functions

#def period(x):
#    q = cardinal(domain(x))
#    D = divisors(q-1)
#    for d in D:
#        if x**d == 1: return d
##
#ord = period   

## Polynomial functions (univariate)
def trailing_degree(f, var = None):
    return f.trailing_degree(var)

def trailing_coeff(f, var = None):
    return f.trailing_coeff(var)
    
def leading_coeff(f, var = None):
    return f.leading_coeff(var)

def leading_term(f, var = None):
    return f.leading_term(var)

def trailing_term(f, var = None):
    return f.trailing_term(var)

def constant_coeff(f, var = None):
    return f.constant_coeff(var)

def zero_degree_term(f):
    return f.zero_degree_term()

def derivative(f, symb = None): return f.der(symb)
der = derivative



def coeffs(f, var=None): 
    if isinstance(f, Poly_type):
        return f.coeffs(var)
    if isinstance(f, int):
        return f
    if isinstance(f, float):
        a = f>>QQ(ZZ())
        return a
    if isinstance(f,Base_type.Base_type):
        return f.cp()

def quo_rem(r0,r1):
    if isinstance(r0,int):
        Domr0 = ZZ()
    else:
        Domr0 = r0.domain()
    if isinstance(r1,int):
        Domr1 = ZZ()
    else:
        Domr1 = r1.domain()
    if Domr1 == ZZ():
        Domr1 = Domr0
        r1 = r1>>Domr1
    if Domr0 == ZZ():
        Domr0 = Domr1
        r0 = r0>>Domr0
    if Domr0.is_sub_domain(Domr1):
        Dom = Domr0
    elif Domr1.is_sub_domain(Domr0):
        Dom = Domr1
    elif Domr1._poly_structure:
        if Domr0._poly_structure:
            if K_(Domr1).is_sub_domain(K_(Domr0)):
                Dom = Domr1
                r0 = r0 >> Domr1
            elif K_(Domr0).is_sub_domain(K_(Domr1)):
                Dom = Domr0
                r1 = r1 >> Domr0
            else:
                raise TypeError("Types must be comparables")
        else:
            if Domr0.is_sub_domain(K_(Domr1)):
                Dom = Poly(Domr0)
                r0 = r0>>Dom
                r1 = r1>>Dom
            else:
                raise TypeError("Types must be comparables")
    elif Domr0._poly_structure:
        if Domr1.is_sub_domain(K_(Domr0)):
            Dom = Poly(Domr1)
            r0 = r0>>Dom
            r1 = r1>>Dom
        else:
            raise TypeError("Types must be comparables")
    else:
        raise TypeError("Types must be comparables")
    if not Dom._poly_structure:
        Dom = Poly(Dom)
        r0 = r0>>Dom
        r1 = r1>>Dom
    return Dom._div_res_(r0,r1)

def sugiyama(r0,r1,t):
    v0 = 0; v1 = 1
    while t <= degree(r1):
        [q,r] = quo_rem(r0,r1)
        v = v0-q*v1
        r0 = r1; r1 = r; v0 = v1; v1 = v
    return (v1,r1)

## Polynomial functions (bivariate)


#univariate_polynomial_ring(K,'x')
def bivariate_polynomial_ring(K,x='x',y='y',name=''):
    Kxy = Poly(K)
    vx = Kxy.variable([x])
    vy = Kxy.variable([y])
    return [Kxy,vx,vy]

def multivariate_polynomial_ring(K,*args, name = ''):
    aux = []
    PK = Poly(K, name)
    if isinstance(args[0],list):
        aux = PK.variables(args[0])
        return [PK] +  aux
    for i in range(len(args)):
        aux.append(PK.variable([args[i]]))
    return [PK] + aux
#

def multivariate_named_polynomial_ring(K,x = 'x', d = 1, name = ''):
    xaux = []
    for i in range(1,d+1):
        xaux += [x + i.__str__()]
    aux = multivariate_polynomial_ring(K,xaux,name = name)
    return aux

def polynomial_ring(K,*args, name=''):
    if len(args) == 1:
        if isinstance(args[0],list):
            if len(args[0]) == 1:
                S = univariate_polynomial_ring(K,args[0],name)
                return S
            else:
                if len(args[0]) == 2:
                    if isinstance(args[0][1], int):
                        return multivariate_named_polynomial_ring(K,args[0][0],args[0][1],name)
                    else:
                        return multivariate_polynomial_ring(K,args[0],name = name)
                else:
                    return multivariate_polynomial_ring(K,args[0],name = name)
        else:
            S = univariate_polynomial_ring(K,args[0], name)
            return S
    else:
        if len(args) == 2:
            if isinstance(args[1], int):
                return multivariate_named_polynomial_ring(K,args[0],args[1],name)
            else:
                return multivariate_polynomial_ring(K,*args,name = name)
        else:
            return multivariate_polynomial_ring(K,*args,name = name)
PR = polynomial_ring

def get_subpolynomial_degree(f,d):
    return f.get_subpolynomial_degree(d)

def variables(f):
    return f.variables()

def mdegree(f):
    return f.mdegree()
total_degree = mdegree

def skeleton(f):
    return f.skeleton()

def partial_derivative(f,x):
    if isinstance(f,int):
        return 0
    if not f.domain()._poly_structure:
        return 0
    return f.der(x)
#
pder = partial_derivative

def lexicographic_order():
    PyECC_set('POL_LEX',1)

def expand_polynomials_on():
    PyECC_set('POL_EXP',1)

def expand_polynomials_off():
    PyECC_set('POL_EXP',0)

def polynomial_variables_order(L):
    if not isinstance(L,list):
        A = [L]
    else:
        A = L[:]
    aux = []
    for i in A:
        if isinstance(i,str):
            aux.append(i)
        if isinstance(i,Poly_type):
            n = i._print_()
            aux.append(n)
    PyECC_set('POL_LEX',0)
    PyECC_set('POL_ORD', aux)

def add_variable(f,x):
    f.add_variable(x)

def evaluate (f,pos,value = ''):
    if value == '':
        flagF = False
        if isinstance(f,int):
            g = clone(f)
            flagF = True
        elif isinstance(f, float):
            g = f>>QQ(ZZ())
            flagF = True
        elif isinstance(f, Zn_type):
            g = clone(f)
            flagF = True
        elif isinstance(f, Fq_type):
            g = clone(f)
            flagF = True
        elif isinstance(f, Power_Series_type):
            raise NotImplementedError
        elif isinstance(f,Vector_type):
            aux = []
            for i in range(len(f)):
                aux.append(evaluate(f[i],pos))
            if f._direction:
                return vector(aux)
            else:
                return vector(aux,'c')
        elif isinstance(f,Matrix_type):
            aux = []
            for i in range(f._rows):
                aux2 = []
                for j in range(f._cols):
                    aux2.append(evaluate(f[i,j],pos))
                aux.append(aux2)
            return matrix(aux)
        elif isinstance(f,QQ_type):
            aux = f._value[:]
            aux = [evaluate(aux[0],pos), evaluate(aux[1],pos)]
            return (aux[0]/aux[1])>>f.domain()
        elif isinstance(f, Poly_type):
            flagF = False
    
        
        if flagF:        
            if isinstance(pos, Poly_type):
                return g
            if isinstance(pos, Monomial_type):
                return g
            if isinstance(pos, Symbol_type):
                return g
            elif isinstance(pos,Vector_type):
                if pos._direction:
                    return vector(g,len(pos))
                else:
                    return vector(f,len(pos),'c')
            elif isinstance(pos, list):
                return [g]*len(pos)
            elif isinstance(pos, int):
                return g
            elif isinstance(pos, float):
                return g
            elif isinstance(pos,Zn_type):
                return g
            elif isinstance(pos,Fq_type):
                return g
            elif isinstance(pos, QQ_type):
                return g
            elif isinstance(pos,Power_Series_type):
                return g
            elif isinstance(pos, Matrix_type):
                return matrix(g,pos._rows, pos._cols)
            else:
                return "evaluate error: Wrong variable parameters"
    
    
        if isinstance(f, Poly_type):
            if isinstance(pos,dict):
                return f.eval_poly(list(pos.values()), list(pos.keys()))
            if isinstance(pos, Poly_type):
                return f.eval_poly(pos)
            if isinstance(pos, Monomial_type):
                return f.eval_poly(pos)
            if isinstance(pos, Symbol_type):
                return f.eval_poly(pos)
            elif isinstance(pos,Vector_type):
                aux = f.eval_poly(pos._value)
                if pos._direction:
                    return vector(aux)
                else:
                    return vector(aux, 'c')

            elif isinstance(pos, list):
                return [f.eval_poly(posi) for posi in pos]
            elif isinstance(pos, int):
                return f.eval_poly(pos)
            elif isinstance(pos, float):
                return f.eval_poly(pos>>QQ(ZZ()))
            elif isinstance(pos,Zn_type):
                return f.eval_poly(pos)
            elif isinstance(pos,Fq_type):
                return f.eval_poly(pos)
            elif isinstance(pos, QQ_type):
                return f.eval_poly(pos)
            elif isinstance(pos,Power_Series_type):
                return f.eval_poly(pos)
            elif isinstance(pos, Matrix_type):
                aux = []
                for i in range(pos._rows):
                    aux2 = []
                    for j in range(pos._cols):
                        aux2.append(f.eval_poly(pos[i,j]))
                    aux.append(aux2)
                return matrix(aux)
            else:
                return "evaluate error: Wrong pos parameters"
        if isinstance(f, Monomial_type):
            if isinstance(pos,dict):
                return f.eval_monomial(list(pos.values()), list(pos.keys()))
            if isinstance(pos, Poly_type):
                return f.eval_monomial(pos)
            if isinstance(pos, Monomial_type):
                return f.eval_monomial(pos)
            if isinstance(pos, Symbol_type):
                return f.eval_monomial(pos)
            elif isinstance(pos,Vector_type):
                aux = f.eval_monomial(pos._value)
                if pos._direction:
                    return vector(aux)
                else:
                    return vector(aux, 'c')

            elif isinstance(pos, list):
                return [f.eval_monomial(posi) for posi in pos]
            elif isinstance(pos, int):
                return f.eval_monomial(pos)
            elif isinstance(pos, float):
                return f.eval_monomial(pos>>QQ(ZZ()))
            elif isinstance(pos,Zn_type):
                return f.eval_monomial(pos)
            elif isinstance(pos,Fq_type):
                return f.eval_monomial(pos)
            elif isinstance(pos, QQ_type):
                return f.eval_monomial(pos)
            elif isinstance(pos,Power_Series_type):
                return f.eval_monomial(pos)
            elif isinstance(pos, Matrix_type):
                aux = []
                for i in range(pos._rows):
                    aux2 = []
                    for j in range(pos._cols):
                        aux2.append(f.eval_monomial(pos[i,j]))
                    aux.append(aux2)
                return matrix(aux)
            else:
                return "evaluate error: Wrong pos parameters"
        
        if isinstance(f, Symbol_type):
            if isinstance(pos,dict):
                return f.eval_symbol(list(pos.values()), list(pos.keys()))
            if isinstance(pos, Poly_type):
                return f.eval_symbol(pos)
            if isinstance(pos, Monomial_type):
                return f.eval_symbol(pos)
            if isinstance(pos, Symbol_type):
                return f.eval_symbol(pos)
            elif isinstance(pos,Vector_type):
                aux = f.eval_symbol(pos._value)
                if pos._direction:
                    return vector(aux)
                else:
                    return vector(aux, 'c')

            elif isinstance(pos, list):
                return [f.eval_symbol(posi) for posi in pos]
            elif isinstance(pos, int):
                return f.eval_symbol(pos)
            elif isinstance(pos, float):
                return f.eval_symbol(pos>>QQ(ZZ()))
            elif isinstance(pos,Zn_type):
                return f.eval_symbol(pos)
            elif isinstance(pos,Fq_type):
                return f.eval_symbol(pos)
            elif isinstance(pos, QQ_type):
                return f.eval_symbol(pos)
            elif isinstance(pos,Power_Series_type):
                return f.eval_symbol(pos)
            elif isinstance(pos, Matrix_type):
                aux = []
                for i in range(pos._rows):
                    aux2 = []
                    for j in range(pos._cols):
                        aux2.append(f.eval_symbol(pos[i,j]))
                    aux.append(aux2)
                return matrix(aux)
            else:
                return "evaluate error: Wrong pos parameters"
            
    #3 parametres    
    if isinstance(pos, Poly_type):
        ppos = [pos]
    if isinstance(pos, Monomial_type):
        ppos = [pos]
    if isinstance(pos, Symbol_type):
        ppos = [pos]
    elif isinstance(pos,Vector_type):
        ppos = pos._value[:]
    elif isinstance(pos, list):
        ppos = pos
    elif isinstance(pos, str):
        ppos = [pos]
    else:
        return "evaluate error: Wrong pos parameters"
    if isinstance(value,int):
        pvalue = [value]
    elif isinstance(value, float):
        pvalue = [value >> QQ(ZZ())]
    elif isinstance(value, Zn_type):
        pvalue = [value]
    elif isinstance(value, Fq_type):
        pvalue = [value]
    elif isinstance(value, Power_Series_type):
        pvalue = [value]
    elif isinstance(value, Vector_type):
        pvalue = copy.deepcopy(value._value)
    elif isinstance(value, Poly_type):
        pvalue = [value]
    elif isinstance(value, Monomial_type):
        pvalue = [value]
    elif isinstance(value, Symbol_type):
        pvalue = [value]
    elif isinstance(value, list):
        pvalue = value
    else:
        return "evaluate error: Wrong value parameters"
    if len(ppos) == len(pvalue):
        if isinstance(f, int):
            return copy.deepcopy(f)
        if isinstance(f, Zn_type):
            return f.cp()
        if isinstance(f, Fq_type):
            return f.cp()
        if isinstance(f, Vector_type):
            aux = []
            for i in range(len(aux)):
                aux.append(evaluate(f[i],ppos,pvalue))
            if f._direction:
                return Vector_type(copy.deepcopy(f.domain()),aux)
            else:
                return Vector_type(copy.deepcopy(f.domain()),aux,'c')
        if isinstance(f, Matrix_type):
            aux = []
            for i in range(f._rows):
                aux2 = []
                for j in range(f._cols):
                    aux2.append(f[i,j],ppos,pvalue)
                aux.append(aux2)
            return Matrix_type(copy.deepcopy(f.domain()),aux)
        if isinstance(f,QQ_type):
            aux = [evaluate(f._value[0],ppos,pvalue), evaluate(f._value[1],ppos,pvalue)]
            return (aux[0]/aux[1])>>f.domain()
        if isinstance(f,Power_Series_type):
            raise NotImplementedError
        if isinstance(f, Poly_type):
            return f.eval_poly(pvalue,ppos)
        if isinstance(f, Monomial_type):
            return f.eval_monomial(pvalue,ppos)
        if isinstance(f, Symbol_type):
            return f.eval_symbol(pvalue,ppos)
    else:
        return "evaluate error: Different number of variables and values"

def symbol(*L):
    if isinstance(L[0],list):
        if len(L[0]) == 2:
            if isinstance(L[0][1], int):
                x = L[0][0]
                n = L[0][1]
                aux = []
                for i in range(1,n+1):
                    aux.append(Symbol_type(x + i.__str__()))
                return aux
        aux = []
        for i in L[0]:
            aux.append(Symbol_type(i))
        return aux
    else:
        if len(L) == 1:
            return Symbol_type(L[0])
        if len(L) == 2:
            if isinstance(L[1], int):
                x = L[0]
                n = L[1]
                aux = []
                for i in range(1,n+1):
                    aux.append(Symbol_type(x + i.__str__()))
                return aux
        aux = []
        for i in range(len(L)):
            aux.append(Symbol_type(L[i]))
        return aux
                        

def clone(A):
    return copy.deepcopy(A)                   
        

#def evaluate(f,*args): 
#    if len(args) == 1:
#        if isinstance(args[0],list):
#            args = args[0]
#        elif isinstance(args[0],Vector_type):
#            args = (args[0])[:]
#    if len(args) > f.domain()._multi_variable:
#        raise TypeError("Too many variable to evaluate")
#    aux = f
#    aux2 = list(args)
#    aux2.reverse()
#    for i in aux2:
#        aux = aux.eval_poly(i)
#    return aux

## List and Vector functions
def serie_const_func(x,a):
    if x == 0:
        return a
    return 0

def serie_from_list(x,a):
    if x>-1 and x<len(a):
        return a[x]
    return 0

def power_series(a,symbol = 'x',point = 0, n = 3):
    if isinstance(a,int):
        f = lambda x: serie_const_func(x,a)
        D = Power_Series(ZZ(),symbol)
        return Power_Series_type(D,f,point,n)
    if isinstance(a,float):
        aux = a>>QQ(ZZ())
        return power_series(aux,symbol,point,n)
    if isinstance(a,Poly_type):
        aux = a.to_list()
        aux.reverse()
        f = lambda x: serie_from_list(x,aux)
        D = Power_Series(a.domain().subdomain(),symbol)
        return Power_Series_type(D,f,point,n)
    
    if isinstance(a,Monomial_type):
        aux = a.to_list()
        aux.reverse()
        f = lambda x: serie_from_list(x,aux)
        D = Power_Series(a.domain().subdomain(),symbol)
        return Power_Series_type(D,f,point,n)
    
    if isinstance(a,Symbol_type):
        aux = a.to_list()
        aux.reverse()
        f = lambda x: serie_from_list(x,aux)
        D = Power_Series(a.domain().subdomain(),symbol)
        return Power_Series_type(D,f,point,n)
    if isinstance(a,Power_Series_type):
        f = a._value
        if point != a._point:
            raise TypeError("Point not match")
        D = Power_Series(a.domain().subdomain(),symbol)
        return Power_Series_type(D,f,point,n)
    if isinstance(a,Base_type.Base_type):
        f = lambda x: serie_const_func(x,a)
        D = Power_Series(a.domain(),symbol)
        return Power_Series_type(D,f,point,n)
    if isinstance(a,list):
        D = max_domain(a)
        f = lambda x: serie_from_list(x,a)
        D = Power_Series(D,symbol)
        return Power_Series_type(D,f,point,n)
        
def is_zero(s):
    if isinstance(s,int):
        if s==0:
            return True
        else:
            return False
    if isinstance(s,float):
        a = s>>Q_
        return is_zero(a)
    if isinstance(s,Symbol_type):
        return s.is_zero()
    if isinstance(s,Base_type.Base_type):
        return s.is_zero()
    for x in s:
        if not is_zero(x): return False
    return True

def prod(a,b):
    return vector(a)*vector(b).transpose()
    
def dot(a,b):
    s = 0
    J = range(min(len(a),len(b)))
    for j in J:
        s += a[j]*b[j]
    return s
    
def zero_positions(f,a):
    return [j for j in range(len(a)) if evaluate(f,a[j])==0]

def rd(K,s=1,r=''): 
    X=[]
    for k in range(s):
        if r == '': 
            X += [K._rd()]
        else: 
            X += [K._rd(r)]
    return X

def rd_nonzero(K,s=1,r=''): 
    X=[]
    for k in range(s):
        if r == '': 
            X += [K._rd_nonzero()]
        else: 
            X += [K._rd_nonzero(r)]
    return X

def rd_vector(K,n):
    x = []
    for j in range(n):
        x += rd(K)
    return vector(x)
    
#def rd_error_vector(K,n='',s=1):
#    if n=='': n = cardinal(K)-1
#    v = create_vector(K,n)
#    S = rd_choice(list(range(n)),s)
#    E = rd_nonzero(K,s)
#    if s == 1:
#        v[S] = E
#    else:
#        for i in range(s):
#            v[S[i]] = E[i]
#    return v
#
def error_vector(n,E):
    #K = (E[0][1]).domain()
    D = max_domain([E[k][1] for k in range(len(E))])
    #print(D)
    e = vector(D,n)
    for k in range(len(E)):
        e[E[k][0]] = E[k][1]>>D
    return e

def rd_linear_combination(G,K=''):
    if K=='': K = K_(G)
    k = len(G)
    u = vector(rd(K,k))
    return u*G

def null_vector(K = ''):
    if K == '':
        A = ZZ()
    else:
        A = K
    return vector(A,[])

def null_matrix(K = ''):
    if K == '':
        A = ZZ()
    else:
        A = K
    return matrix(A,[])

def invert_entries(a):
    return create_vector([1/t for t in a])

## Matrix functions
def I_(n, K=Zn(2)): 
    I = create_matrix(K,n,n)
    for k in range(n):
        I[k,k]=1
    return I

def subs_element(A,x,y):
    if isinstance(A,Vector_type):
        for i in range(len(A)):
            if A[i] == x:
                A[i] = y
    if isinstance(A,Matrix_type):
        for i in range(A.nrows()):
            for j in range(A.ncols()):
                if A[i,j] == x:
                    A[i,j] = y
                
    

def combination(n,k):
    N = list(range(n))
    return sorted(rd_choice(N,k))

def rd_permutation(n):
    return rd_choice(list(range(n)),n)
permutation = rd_permutation

def rd_permutation_matrix(n):
    N = list(range(n))
    p = rd_permutation(n)
    return permutation_matrix(p)
    
def permutation_matrix(p):
    if isinstance(p, int):
        return rd_permutation_matrix(p)
    n = len(p)
    P = create_matrix(ZZ(),n,n)
    for j in range(n):
        P[j,p[j]] = 1
    return P

def scramble_matrix(A,k):
    U = create_matrix(A,k,k)
    L = create_matrix(A,k,k)
    for i in range(k):
        U[i,i] = L[i,i] = 1
        for j in range(i+1,k):
            U[i,j]=rd(A)
            L[j,i]=rd(A)
    P = permutation_matrix(k)
    return P*U*L

# auxiliary function for rd_GL: Given a list L, returns (e,r,a) where e = e_r,  r is a random 
#index of the list L and a is vector [0...0,rdnonzero(K),rd_vector(K,n-r-1)] 
def rd_insert(A,r):
    n = ncols(A)
    if r<0 or r>n: return "r has to be in 0..n"
    K = K_(A)
    B = matrix(K,n+1,n+1)
    B[0,r] = 1
    for j in range(1,n+1):
        B[j,r] = rd(K)
    if r==0:
        B[1:,1:] = A[:,:]
    elif r==n:
        B[1:,0:n] = A[:,:] 
    else:
       B[1:,0:r] = A[:,0:r] 
       B[1:,r+1:n+1] = A[:,r:n] 
    return B
  
def rd_extend(A):
    n = ncols(A); K = K_(A)
    v = rd_nonzero_vector(K,n+1)
    r = 0
    for j in range(n+1):
        if v[j]!=0: 
            r = j
            break
    A = rd_insert(A,r)
    x = v[r]
    for j in range(r+1,n+1):
        A[:,j] = A[:,j]+ A[:,r]*v[j]
    A[:,r] = x*A[:,r]
    return A

def rd_GL(n,F=Zn(2)):
    a = rd_nonzero(F)[0]  
    A = matrix([a])
    for _ in range(2,n+1):
        A = rd_extend(A)
    return A
    
def rd_nonzero_vector(K,n):
    while True:
        v = rd_vector(K,n)
        if not is_zero(v): return v

def transpose(H): 
    if not hasattr(H,'transpose'):
        return H
    return H.transpose()

def kernel(A):
    return A.domain().kernel(A)
right_kernel = kernel

def shape(A):
    if isinstance(A,Matrix_type):
        return (A.nrows(), A.ncols())
    if isinstance(A,Vector_type):
        if A._direction:
            return (1,len(A))
        else:
            return (len(A),1)
    
def rank(A): return A.rank()

def det(A): return A.det()

#def nrows(M): return len(M)
#def ncols(M): return len(M.T) # (M.shape)[1]

def alternant_matrix(h,a,r):
    h = vector(h); a = vector(a)
    V = vandermonde_matrix(a,r)
    for j in range(r):
        V[j] = h * V[j]
    return V    


def blow(A,K):
    return A.domain().blow(A,K)

def prune(A):
    return A.domain().prune(A)

def hankel_matrix(s,r='', c=''):
    t = len(s)//2
    if r == '': r = t
    if c == '': c = t + 1
    if (r + c - 1 > len(s)):
        return "Wrong number of parameters"
    F = max_domain(s)
    S = matrix(F,r,c)
    for j in range(r):
        S[j] = s[j:j+c]
    return S

# Special Gauss-Jordan
def GJ(S):
    S.pivot_matrix()
    for i in range(min(ncols(S),nrows(S))):
        if S.Piv[i,i] == 0:
            break
        l = i+1
    A = S.Piv[:l,:l]
    b = S.Piv[:l,l]
    a = domain(A).solve_pivot_matrix_system(A,b)
    return a

def solve_linear_system(G,a):
    return domain(G).solve_linear_system(G,a)

def solve_pivot_matrix_system(G,a):
    return domain(G).solve_pivot_matrix_system(G,a)

## Decoders
def alternant_decoder(y,C):
    if not isinstance(y,list) and not isinstance(y,Vector_type):
        return 'AD error: wrong input type'
    K = K_(C)
    if K_(y) != K:
        if is_sub_domain(K,K_(y)):
            if isinstance(y,list): y = vector(K,y)
            if not isinstance(y,Vector_type):
                return "AD error: Argument is not a vector"
        else:
            return "AD error: Fields of argument and code incompatible"
    else: 
        if isinstance(y,list): y = vector(K,y)
        if not isinstance(y,Vector_type):
            return "AD error: Argument is not a vector"
    
    h = h_(C)
    if len(y) != len(h):
        return "AD error: Vector argument has wrong length"
    r = r_(C); a = a_(C); b = b_(C)
    H = H_(C); K1 = K_(H)
    y1 = pull(y.cp(),K)
    s = y1 * transpose(H)
    
    if is_zero(s): 
        print("AD: Input is a code vector")
        return y1
    [_,z] = univariate_polynomial_ring(K1,'z','K[z]')
    S = polynomial(s,z)
    
    if degree(S)< r//2 or trailing_degree(S)>= r//2:
        return "AD: Non-decodable vector: extraneous syndrome"
        
    (L,E)=sugiyama(z**r,S,r//2)
    
    M = zero_positions(L,b)
    
    if len(M)< degree(L): return 'AD error: undecodable vector: too few roots'
    
    if K == Zn(2):
        for j in range(len(M)):
            y1[M[j]] += 1 
        show('Corrected positions:',M)
        return y1
    
    if len(M)<degree(L):
        return "AD error: non-decodable vector: too few locator roots"
    
    L = der(L)

    err = []
    for m in M:
        x = (a[m]*evaluate(E,b[m]))/(h[m]*evaluate(L,b[m]))
        x = pull(x,K)
        err += [-x]
        y1[m] = y1[m]+x
    show("err =",err)
    if is_zero(y1*H.transpose()):
        print("AD: Successful decoding with error table:\n",M,'\n',vector(err))
        return y1
    else: return "AD error: output vector has non-zero syndrome"
# alias
BMS = AD = alternant_decoder

# C is a code, t>0 the number of errors introduced, a suitable finite field
def AD_checker(C,t,K=''):
    if K=='': K = K_(C)
    if G_(C)!= None:
        G = G_(C)
    else: 
        H = prune(blow(H_(C),K))
        G = left_kernel(transpose(H))
    n = ncols(G)
    x = rd_linear_combination(G,K)
    e = rd_error_pattern(n,t,K)
    x1 = alternant_decoder(x+e,C)
    if isinstance(x1,Vector_type):
        if x1==x:
            show('AD_checker: trial successful')
            return matrix([x,e,x+e])
        else:
            show('AD_checker: undectetable error')
            return matrix(x,e,x+e,x1)
    else:
        show('AD_checker: decoder error')
        return matrix([x,e])
decoder_trial = AD_checker

# Short hand decoder trial
def SHDT(C,s,K=''):
    if K == '': K = K_(C)
    if G_(C)!= None: G = G_(C)
    else: return "decoder_trial: G is not defined"
    n = ncols(G)
    x = rd_linear_combination(G,K)
    e = rd_error_vector(K,n,s)
    y = x+e
    x1 = AD(y,C)
    if x1==x: return [list(x),list(e),list(y)]
    if isinstance(x1,str): return [list(x),list(e)]
    return [list(x),list(e),list(y), list(x1)]


# The received vector y can be a list or a vector  
# C is an alternant code #
def PGZ(y,C):
    if isinstance(y,list): y = vector(K_(C),y)
    if not isinstance(y,Vector_type):
        return "PGZ> Argument is not a vector"
    h = h_(C)
    if len(y) != len(h):
        return "PGZ> Vector argument has wrong length"
    r = r_(C)
    alpha = a_(C); beta = b_(C)
    H = H_(C)
    K = K_(C)
    s = y*H.transpose()
    if is_zero(s): 
        print("PGZ> Input is a code vector")
        return pull(y.cp(),K)
    
    S = hankel_matrix(s)
    S.pivot_matrix()
    l = min(nrows(S),ncols(S))
    for i in range(min(ncols(S),nrows(S))):
        if S.Piv[i,i] == 0:
            l = i
            break
    A =  S.Piv[:l,:l]
    b = -S.Piv[:l,l]
    a = domain(A).solve_pivot_matrix_system(A,b)
      
    a = a.to_list()
    K1 = K_(H)
    [_,z] = univariate_polynomial_ring(K1,'z','K1[z]')
    L = hohner(a+[1],z)
    R = [s for s in beta if evaluate(L,s)==0]
    if len(R) < l:
        return "PGZ> Defective error location"
    M = [beta.to_list().index(r) for r in R]
    s1 = s.to_list(); #s1.reverse()
    s1 = polynomial(s1,z)
    E = md_mul(L,s1,z**r)  # L * S(z) mod z^r
    L = L.der()            # only de derivative is needed
    # Correct with Forney's formula
    w = [-alpha[m]*evaluate(E,beta[m])/(h[m]*evaluate(L,beta[m])) for m in M]                                
    show("PGZ> Error positions and values", M, vector(w))
    y1 = y.cp()
    for m in M:                           
        y1[m] += alpha[m]*evaluate(E,beta[m])/(h[m]*evaluate(L,beta[m]))
    return pull(y1,K)
    
def PGZm(y,C):
    if isinstance(y,list): y = vector(K_(C),y)
    if not isinstance(y,Vector_type):
        return "PGZm> Argument is not a vector"
    h = h_(C)
    if len(y) != len(h):
        return "PGZm> Vector argument has wrong length"
    r = r_(C)
    alpha = a_(C)
    H = H_(C)
    K = K_(C)
    #y1 = vector([t for t in y])
    s = y*H.transpose()   # y1 before
    #show("len(s)",len(s))
    if is_zero(s): 
        print("PGZm> Input is a code vector")
        return pull(y.cp(),K)
    S = hankel_matrix(s)
    c0 = S[:,0]   # keep the first column of S
    S.pivot_matrix()
    l = min(ncols(S),nrows(S))
    for i in range(min(ncols(S),nrows(S))):
        if S.Piv[i,i] == 0:
            l = i
            break
    A = S.Piv[:l,:l]
    b = -S.Piv[:l,l]
    a = solve_pivot_matrix_system(A,b)
    a = a.to_list();
        
    K1 = K_(H)
    [_,z] = univariate_polynomial_ring(K1,'z','K1[z]')
    L = polynomial(a+[1],z)
    R = [s for s in alpha if evaluate(L,s)==0]
    if len(R) < l:
        return "PGZm> Defective error location"
    M = [alpha.to_list().index(r) for r in R]
    h1 = [h[m] for m in M]
    V = alternant_matrix(h1,R,l)
    v = c0[:l]
    w = solve_linear_system(V,v)
    y1 = y.cp()
    show("PGZm> Error positions and values", M, w.transpose())
    for j in range(len(M)):
        y1[M[j]]-=w[j]
    return pull(y1,K)

'''
def RSID(a,k,y):
    n = len(a); t = (n-k)//2
    P = transpose(vandermonde(a,n-t)); Q = transpose(vandermonde(a,n-k-t))
    for j in range(n):
        Q[j] = y[j]*Q[j]
    PQ = kernel(splice(P,Q))[:,0]
    P = PQ[:(n-t)]
    Q = PQ[(n-t):]
    show('P:',P)
    show('Q:',Q)
    [_,T] = univariate_polynomial_ring(field(a),'T')
    P = polynomial(P,T)
    Q = polynomial(Q,T)
    show('P:',P)
    show('Q:',Q)
    if Q == 0: return "RSID: Error"
    if not P%Q: return "RSID: Error"
    f = -P//Q
    return vector([evaluate(f,t) for t in a])
'''

def left_kernel(A):
    return transpose(kernel(transpose(A)))

def trace(A):
    [m,n] = shape(A)
    r = min(m,n)
    return sum([A[j,j] for j in range(r)])

def Tr(x,K='',L=''):
    k = prime_field(K_(x))
    if K=='': K = k
    Cx = conjugates(x,K)
    if L=='': s = 1
    else: s = dim(L,K)//len(Cx)
    return s * sum(Cx)

##20181019
    

def parity_completion(G): 
    return splice(G,[-sum(g) for g in G])


def index(s,H):
    n = len(H)
    for j in range(n):
        if s==H[j]: return j
    return 0
    
def P_e(n,t,p):
    err = 0
    for j in range(t+1,n+1):
        err += binom(n,j) * p**j * (1-p)**(n-j)
    return round(err,10)

def erf(n,t,p):
    if p==0: return 1
    return round(P_e(n,t,p)/p,10)

def volume(n,r,q=2):
    v = 0
    b = q-1
    for i in range(r+1):
        v += binom(n,i) * b**i
    return v
    
def ub_sphere(n,d,q=2):
    return floor(q**n/volume(n,floor((d-1)/2),q))
   
def lb_gilbert(n,d,q=2):   
    return ceil(q**n/volume(n,d-1,q))

def isbn(x):
    r = (sum([(j+1)*x[j] for j in range(len(x))])) % 11
    if r==10: return 'X'
    return r    

def flip(x,j=''):   
    n = len(x)
    y = clone(x)
    if j == '':
        for j in range(len(x)):
            y[j] = 1 - y[j]
        return y
    if isinstance(j,int):
        if 0<=j and j<n:
            y[j] = 1-y[j]
            return y
        else: 
            return "flip error: index out of range"
    if isinstance(j,list):
        for k in j:
            if 0<=k and k<n:
                y[k] = 1-y[k]
            else:
                return "flip error: index out of range"
        return y
    if isinstance(j,Vector_type):
        return flip(x,list(j))
    return "flip error: wrong data"

def hadamard_matrix_recursive(n):  
    if n==0: return matrix([[1]])
    if n > 0: 
        H = hadamard_matrix_recursive(n-1)
        k = nrows(H)
        return stack(splice(H,H), splice(H,-H))

def paley_matrix(F):
    if characteristic(F)==2:
        return "paley_matrix error: the characteristic of F must not be 2"
    q = cardinal(F)
    def x(j): return element(j,F)
    def chi(a): return legendre(a,F)
    return matrix([[chi(x(i)-x(j)) for j in range(q)] for i in range(q)])

def hadamard_matrix_paley(K):
    q = cardinal(K)
    if q%4 != 3:
        return "hadamard_matrix_FF Error: {} is not 3 mod 4".format(q)
    A = I_(q,Z_) + paley_matrix(K)
    u = matrix(q*[1])
    x = splice(matrix([1]),u)
    X = splice(transpose(u),-A)
    return stack(x,X)

def conference_matrix(K):  
    q = cardinal(K)
    e = ((q-1)//2)%2
    e==(-1)**e
    S = paley_matrix(K)
    u = matrix(q*[1])
    x = splice(matrix([0]),u)
    X = splice(transpose(e*u),S)
    return stack(x,X)

def hadamard_matrix_finite_field(K): 
    q = cardinal(K)
    e = ((q-1)//2)%2
    #e = (-1)**e
    I = I_(q+1,Z_); 
    C = conference_matrix(K)
    if e==0:
        return stack(splice(I+C,-I+C),splice(-I+C,-I-C))
    return hadamard_matrix_paley(K)

def hadamard_code(X):  
    if isinstance(X,int):
        X = hadamard_matrix_recursive(X)
    elif isinstance(X,Matrix_type):
        pass
    elif is_finite_field(X):             
        X = hadamard_matrix_finite_field(X)
    
    C = stack(X,-X)
    subs_element(C,-1,0)
    return C

def hadamard_decoder(y,H):
    if len(y) != len(H):
        return "hadamard_decoder error: length mismatch"
    r = len(H)//2
    z = unbin_(y)
    for h in H:
        c = dot(z,h)
        if abs(c) > r:
            if c > 0:
                return bin_(h)
            else:
                return bin_(-h)
    return "hadamard_decoder error: undecodable"    

def hadamard_matrix(X):
    if isinstance(X,int):
        return hadamard_matrix_recursive(X)
    return hadamard_matrix_finite_field(X)
    

def RSID(a,k,y):   
    n = len(a); t = (n-k)//2
    P = vandermonde(a,n-t); Q = vandermonde(a,n-k-t+1)
    for j in range(n):
        Q[:,j] = y[j] * Q[:,j] 
    PQ = left_kernel(stack(P,Q))[0]
    P = PQ[:(n-t)]; Q = PQ[(n-t):]
    if is_zero(Q): return "RSID Error: Q = 0"
    [_,T] = univariate_polynomial_ring(K_(a),'T')
    P = polynomial(P,T); Q = polynomial(Q,T)
    if not P % Q: return "RSID Error: Q does not divide P"
    f = -P//Q
    return evaluate(f,a)


def weight_enumerator_mds(n,k,q=2):
    if min(k,n-k)<2 or max(k,n-k) >= q:
        return "weight_enumerator_mds error: wrong parameters"
    d = n-k+1
    W = [1] + (d-1)*[0]
    for i in range(d,n+1):
        w = 0
        for j in range(i-d+1):
            w += (-1)**j * binom(i-1,j) * q**(i-d-j)
        w *= binom(n,i)*(q-1)
        W += [w]
    return W
    
def min_weight(X):
    return min([wt(x) for x in X])
    
def min_weights(X):
    m = min_weight(X)
    return [x for x in X if wt(x)==m]

def submatrix(M,J, I=''):
    if (I == ''):
        return M[:,J]
    else:
        return M[I,J]
    
def lb_gilbert_varshamov(n,d,q=2):
    v = volume(n-1,d-2,q)
    k = 0
    while q**(n-k) > v:
        k += 1
    return q**(k-1)
    
def normalized_hamming_matrix(r,F=Zn(2)):
    q = cardinal(F); n = (q**r-1)//(q-1)
    K = [element(j,F) for j in range(q)]
    B = [(r-1)*[K[0]]+[K[1]]]
    H = B
    for j in range(1,r):
        S = []
        for b in B:
            s = []
            b = b[1:]
            for t in K:
                s = s+[b+[t]]
            S = S+s
        B=S
        H = H+S
    return transpose(matrix(H))

def hamming_weight_enumerator(r,q=2):
    e = q**(r-1); f = e*q
    n = (e-1)//(q-1)
    [_,T] = univariate_polynomial_ring(Q_,'T')
    a = 1 + (q-1)*T
    return a**n * (a**e + (f-1) * (1-T)**e) // f

HWE = hamming_weight_enumerator

def macwilliams(n,k,A,q=2):  ###
    [_,T] = univariate_polynomial_ring(Q_,'T')
    Q = QQ(_)
    a = 1 + (q-1)*T >> Q
    B = a**n * evaluate(A>>Q,(1-T)/a) // q**k  # ??
    return pull(B)
    
def tensor(A,B):   
    (ra,ca) = shape(A)
    (rb,cb) = shape(B)
    m = ra*rb; n = ca*cb
    T = matrix(Z_,m,n)
    for i in range(m):
        for j in range(n):
            (ia,ja) = (i//rb,j//cb)
            (ib,jb) = (i%rb,j%cb)
            T[i,j] = A[ia,ja]*B[ib,jb]
    return T
    
def legendre(x,F):   
    if x==0: return 0
    q = cardinal(F)
    p = characteristic(F)
    if p==2:
        return "legendre: {} is even, so legendre is not defined".format(q)
    u = x**((q-1)//2)
    if u==1: return 1
    else: return -1
    
    
def QR(F):   
    X = Set(F)[1:]
    q = len(X)+1
    if q%2==0: return X
    R = []
    for x in X:
        if x**((q-1)//2)==1: R += [x]
    return R
#
def QNR(F):
    X = Set(F)[1:]
    q = len(X)+1
    NR = []
    if q%2==0: return NR
    for x in X:
        if x**((q-1)//2)!= 1: NR += [x]
    return NR


##20181021
    

def bin_(x):
    v = vector(Z_,len(x))
    for j in range(len(v)):
        if x[j] == -1:
            v[j] = 0
        else:
            v[j] = x[j]
    return v

def unbin_(x):
    v = vector(Z_,len(x))
    for j in range(len(v)):
        if x[j]==0:
            v[j] = -1
        else:
            v[j] = x[j]
    return v

def paley_code(F):
    q = cardinal(F)
    I = I_(q,Z_)
    U = matrix(1,q)
    S = paley_matrix(F)
    z = matrix(0,1,q)
    C = stack(z,I+S+U)
    C = stack(C,-S+I+U)
    C = C//2
    u = matrix(1,1,q)
    C = stack(C,u)
    return C

## Bit-product of integers and Hadamard matrices

def bit_product(a,b,r=''):
    a = reverse(bin(a)[2:]); b = reverse(bin(b)[2:])
    l = min(len(a),len(b))
    if r=='': r = l
    else: r = min(l,r)
    a = [int(x) for x in a][:r]; b = [int(x) for x in b][:r]
    return  dot(a,b) % 2
#
bdot = bit_product

def hadamard(r):
    if not isinstance(r,int): return "Error hadamard: parameter is not an integer"
    if not r>0: return "Error hadamard: {} is not a positive integer".format(r)
    n = 2**r
    return matrix([[(-1)**bit_product(i,j,n) for j in range(n)] for i in range(n)])

def ub_griesmer(n,d,q=2):
    P = 1
    while n > 0:
        n = n - ceil(d/P)
        P = P*q
    return P

def ub_plotkin(n,d,q=2):
    b = (q-1)/q
    if d > b*n:
        if q==2:
            if d % 2:
                n += 1; d += 1
            return 2 * floor(d/(2*d-n))
        else:
            return floor(d/(d-b*n))
    m = (b*n-d)/b
    if m == int(m):
        m +=1
    else:
        m = ceil(m)
    return int(q**m * ub_plotkin(n-m,d,q))

def ub_elias(n,d,q=2):
    b = (q-1)/q
    def f(r):
        return r**2-2*b*n*r+b*n*d
    R = range( floor(b*n) +1 )
    M = q**n
    for r in R:
        if f(r)>0:
            A = floor((b*n*d/f(r))*q**n/volume(n,r,q))
            if A<M: M = A
    if q==2 and d%2:
        M = min(M,ub_elias(n+1,d+1))
    return M

def ub_johnson(n,d,w=''):
    if isinstance(w,int):
        h = floor((d+1)/2)
        b = 1
        if h <= w:
            for i in range(w-h,-1,-1):
                b = floor(b*(n-i)/(w-i))
            return b
        else:
            return 1
    if d%2==0:
        n = n-1; d=d-1
    t = floor(d/2)
    x = binom(n,t+1)-binom(d,t)*ub_johnson(n,d,d)
    x = x/floor(n/(t+1))
    for i in range(t+1):
        x += binom(n,i)
    return floor(2**n/x)

def newton(j):
    [_,x] = univariate_polynomial_ring(Q_,'x','Q[x]')
    N = 1
    for i in range(1,j+1):
        N *= x-i+1
    return N/factorial(j)
    
def krawtchouk(k,n,q=2):
    if k==0: return 1
    [_,x] = univariate_polynomial_ring(Q_,'x','Q[x]')
    def n_(j):
        N = 1
        for i in range(1,j+1):
            N *= (x-i+1)
        return N/factorial(j) 
    def m_(j):
        M = 1
        for i in range(1,j+1):
            M *= n-x-i+1
        return M/factorial(j)
    S = 0
    for j in range(k+1):
        sj = (-1)**j * n_(j)*m_(k-j)*(q-1)**(k-j)
        S += sj
    return S

def lies(x,a=0,b=1):
    if not isinstance(x,float): 
        return False
    if x<a or isinstance(b,float) and x>b:
        return False
    return True

def entropy(x,q=2):
    if not lies(x,0,1):
        return "Error entropy: wrong parameter"
    if not isinstance(q,int): 
        return "Error entropy: q is not an integer"
    if q<2: 
        return "Error entropy: q<2"
    def lg(x): return log(x,q)
    if x==0*1.0: return 0.0
    elif x==1.0: return log(q-1,q)
    else: return x*lg(q-1) - x*lg(x) - (1-x)*lg(1-x)
    
def singleton(x):
    if not lies(x,0,1):
        return "Error singleton: wrong parameter"
    return 1-x

def plotkin(x):
    if not lies(x,0,1/2):
        return "Error plotkin: wrong parameter"
    return 1-2*x

def gilbert(x):
    if not lies(x,0,1):
        return "Error gilbert: wrong parameter"
    if x<=1/2: return 1-entropy(x)
    else: return 0.0

def hamming(x):
    if not lies(x,0,1):
        return "Error hamming: wrong parameter"
    return 1-entropy(x/2)
    
def elias(x):
    if not lies(x,0,0.5):
        return "Error elias: wrong parameter"
    return 1-entropy((1-sqrt(1-2*x))/2)
    
def vanlint(x):
    if not lies(x,0,0.5):
        return "Error vanlint: wrong parameter"
    return entropy(1/2 - sqrt(x*(1-x)))
    
def mceliece(x):
    if not lies(x,0,0.5):
        return "Error mceliece: wrong parameter"
    def g(t): return entropy((1-sqrt(1-t))/2)
    h = lambda u: 1+g(u**2)-g(u**2+2*x*u+2*x)
    a = 0.005; b=1-2*x-a; du=0.04*(1-2*x)
    return min([h(a+j*du) for j in range(25)])


##20181024
def is_invertible(x):
    if isinstance(x, int):
        if x==1 or x == -1:
            return True
        else:
            return False
    if isinstance(x,float):
        if x == 0:
            return False
        else:
            return True
    return x.is_invertible()    


#20181129
def rd_pick(P,k): return P._rd(k)

#20190121
def collect(f,var):
    
    if isinstance(f,int):
        if isinstance(var,str):
            [K,x] = polynomial_ring(Z_,'x')
            g = f >> K
        else:
            x = var
            g = f >> x.domain()
    elif isinstance(f, float):
        if isinstance(var,str):
            [K,x] = polynomial_ring(Q_,'x')
            g = (f >> Q_) >> K
        else:
            x = var
            g = (f >> Q_) >> x.domain()
    elif isinstance(var,str):
        if f.domain()._poly_structure:
            x = f.domain().variable(var)
            g = f
        else:
            [K,x] = polynomial_ring(f.domain(),var)
            g = f>>K
    elif not var.domain().is_sub_element(f):
        x = var
        g = clone(f) >> max_domain(x.domain(),f.domain())
    else:
        x = var
        g = clone(f)
    #print(type(g))
    if (not isinstance(g, Poly_type)) and (not isinstance(g, Monomial_type)) and (not isinstance(g, Symbol_type)):
        return 'collect error: Needs a polynomial'
    
    return g.collect(var)
    
def expand_poly():
    PyECC_set('POL_EXP',1)
def not_expand_poly():
    PyECC_set('POL_EXP',0)
#20190131
def is_polynomial(f): return (isinstance(f,Poly_type)) or (isinstance(f,Monomial_type)) or (isinstance(f,Symbol_type))
def is_vector(f): return isinstance(f,Vector_type)
def types_off(): PyECC_set('TYPES',0)
def types_on(): PyECC_set('TYPES',1)
def toggle_types(): 
    if TYPES==0: PyECC_set('TYPES',1)
    else: PyECC_set('TYPES',0)
def decimals(n):
    PyECC_set('Q_DEC',n)
def nodecimals():
    PyECC_set('Q_DEC',-1)

def pad(L,r):
    if r <= len(L): return L[:r]
    return L+(r-len(L))*[0]

def to_list(f):
    if is_polynomial(f): return coeffs(f)
    if is_vector(f): return list(f)
    return f

coeffs_dec = coeffs

def coeffs_inc(f,var=None):
    cf = coeffs(f,var)
    return reverse(cf)

def parity_completion(G): 
    return splice(G,vector([-sum(g) for g in G]))

def left_parity_completion(G): 
    return splice(vector([-sum(g) for g in G]), G)

def cyclic_reduction(f,n):
    r = degree(f)
    a = coeffs_inc(f)
    for j in range(n,r+1):
        c = a[j]
        if c==0: continue
        f = f-c*x**j; f = f+c*x**(j%n)
    return f

def cyclic_product(f,g,n): return cyclic_reduction(f*g,n)

# cyclic matrix of k rows formed with the coeffs of f
def cyclic_matrix(f,k):
    if k<1: return 'Error: {} is not positive'.format(k)
    if is_polynomial(f): f = coeffs_inc(f)
    else: f = to_list(f)
    f = f+(k-1)*[0]
    C = [f]
    for _ in range(1,k):
        f = [0] + f[:-1]
        C += [f]
    return matrix(C)

def cyclic_generating_matrix(g,n):
    if is_polynomial(g): k = n-degree(g)
    else: k = n - (len(g) - 1)
    if k<1: return 'Error: {} less than {}'.format(k,degree(g))
    return cyclic_matrix(g,k)

def cyclic_control_matrix(h,n):
    if is_polynomial(h): h = coeffs(h)
    else: h = reverse(to_list(h))
    return cyclic_generating_matrix(h,n)

def cyclic_normalized_matrix(g,n):
    x = variable(g); k = n-degree(g)
    G = []
    for j in range(k):
        r = - x**(n-k+j) % g
        r = coeffs_inc(r)
        r = pad(r,n-k)
        G += [r]
    return splice(matrix(G),I_(k,Z_))

def cyclic_normalized_control_matrix(g,n):
    G = cyclic_normalized_matrix(g,n)
    R = G[:,list(range(degree(g)))]
    return splice(I_(degree(g),Z_),transpose(-R))

#20190201
# To serch the value associeted to a key k
# in a talbe T given as a list of pairs
def lookup(k,T):
    for (a,b) in T:
        if a==k: return b
    return 0

# Decoder for cyclic codes. Presupposes having
# computed a meggit table E
def meggitt_decoder(y,g):
    x = variable(g)
    s = remainder(y,g)
    if s==0: return y
    j = 0
    while lookup(s,E)==0:
        j += 1; s = remainder(x*s,g)
    e = lookup(s,E)/x**j
    return y-e    

##20190203
    
def bezout(m,n):
    (a,b,d) = extended_euclidean_algorithm(m,n)
    return [d,a,b]
    
def test_bezout(r=10):
    m = ri(r); n = ri(r)
    [d,a,b] = bezout(m,n)
    #show(m)
    if d == a*m+b*n:
        return 'Success'
    else: return ('Fail:',m,n,d)
    
def test_gcd(m=100, n=100):
    m = ri(m); n = ri(n)
    d = gcd(m,n)
    show(d)
    if gcd(m//d,n//d)==1:
        return 'Success'
    else: return 'Failure'

# 20190205

def sqroot(a,F):
    [_,Z] = univariate_polynomial_ring(F,'Z')
    return pull(equal_degree_splitting(Z**2-a,1)-Z, F)
    
# 190211
def index_table(x): return [(0*x,'_')] + [(x**j,j) for j in range(order(x))]

#def inverse(T): return [(b,a) for (a,b) in T]


# 16/2/2019 

def prime_field(F):
    while not is_prime_field(F):
        F = base(F)
    return F


#  17/2/2019

def dim(L,K=''):
    if hasattr(L,'_dim'):
        return L._dim
    p = characteristic(L)
    d = ifactor(cardinal(L))[p]
    if K=='': return d
    if not is_sub_domain(K,L):
        return 'dim error: subfield condition False'
    dK = ifactor(cardinal(K))[p]
    return d//dK  
    
def next_q(r,t=''):
    # Assume x is a positive real number and t a positive integer
    if t=='': n = ceiling(r)
    else: n = ceiling(1+2*t/(1-r))
    while not is_prime_power(n):
        n+=1
    return n

def PRS_G(w,k):
    n = order(w)
    return matrix([geometric_series(w**j,n) for j in range(k)])
    
def mattson_solomon_matrix(w):
    if isinstance(w,Zn) or isinstance(w,Fq):
        w = primitive_root(w)
    n = order(w)
    return matrix([[w**(i*j) for j in range(n)] for i in range(n)])    
    
#20190222
def invert_table(T): return [(b,a) for (a,b) in T]

#20190223
def rd_choice(X, m=1): 
    if isinstance(X,int):
        if X > 0: X = list(range(X))
        else: return 'rd_choice error: {} is not positive'.format(X)
    if isinstance(X,Vector_type): X = list(X)
    if not (isinstance(X,list) or isinstance(X,str)):
        return 'rd_choice error: wrong type of first parameter'
    if len(X)==0: return None
    if len(X)<m: 
        return "rd_choice error: not enough items to choose from"
    Y = clone(X)
    Z = []
    for j in range(m):
        t = rd_int(0,len(Y)-1)
        y = Y[t]
        Z += [y]
        Y = Y[:t] + Y[t+1:]
    if isinstance(X,str): Z = string(Z)
    return Z
rd_selection = rd_sample = rd_choice

def vec(x,K=Zn(2)): return(vector(K,x))

def GF(q,alpha='x'):
    if not is_prime_power(q): 
        return 'GF error: cardinal is not a prime power'
    P = prime_factors(q)
    p = P[0]; m = len(P)
    f = get_irreducible_polynomial(Zn(p),m)
    [F,x] = extension(Zn(p),f,alpha)
    return [F,x]

# 190223

def string(C):
    s = ''
    for c in C: s += c
    return s
    

def rd_error_pattern(n,w,K=Zn(2)):
    # n and w positive integers, w<=n, K a finite ring
    if not isinstance(n,int) or not isinstance(w,int):
        return 'rd_error_pattern error: #1 and #2 must be integers'
    if n<=0 or w<=0 or w>n:
        return 'rd_error_pattern error: #1 and #2 must positve and #2 <= #1'
    E = rd_nonzero(K,w)
    J = rd_choice(n,w)
    e = vec(n,K)
    for k in range(w):
        e[J[k]] = E[k]
    return e
# alias
rd_error_vector = rd_error_pattern


        
# Nota: index_table(t) s'ha modificat (v. supra): en lloc de representar
# 0 per '_', ho he deixat en '' (null string).

# Replaces the non-zero field elements by their logs with respect
# to a primitive root according to a precomputed table E. 
# Conventionally, 0 is represented by null string ''
        
#20190227
def shorthand(X,E):
    if isinstance(X,Vector_type):
        return [lookup(x) for x in X]
    if isinstance(X,Matrix_type):
        S = []
        for x in X:
            r =[]
            for a in x:
                r += [lookup(a,E)]
            S+=[r]
        return S
    return lookup(X,E)


#20191107
def take(v,n):
    if isinstance(n,list):
        if isinstance(v,list):
            c = []
            for l in n:
                c += [v[l]]
            return c
        return v[n]
    if isinstance(n,int):
        if n >= 0:
            if n > len(v):
                return v[:]
            return v[:n]
        if n < 0:
            if -n > len(v):
                return v[:]
            if isinstance(v,list):
                return v[n:]
            else:
                return v[len(v) + n:]

# def take(v,L):
#     if isinstance(L,slice):
#         if isinstance(v,list):
#             return v[L]
#         a0 = L.start
#         if a0 == None:
#             a0 = 0
#             fvec = True
#         if a0 < 0:
#             a0 = len(v) + a0
#         a1 = L.stop
#         if a1 == None:
#             a1 = len(v)
#             fvec = True
#         if a1 < 0:
#             a1 = len(v) + a1
#         a2 = L.step
#         if a2 == None:
#             a2 = 1
#         indexs = list(range(a0,a1,a2))
#         return v[indexs]
#     if isinstance(L,list):
#         if isinstance(v,list):
#             c = []
#             for l in L:
#                 c += [v[l]]
#             return c
#         return v[L]
#     if isinstance(L,int):
#         if L >= 0:
#             return v[:L]
#         if L < 0:
#             if isinstance(v,list):
#                 return v[L:]
#             else:
#                 return v[len(v) + L:]


# Returns the list of pairs (c,m) such that
# c*m is a non-zero monomial of the polynomial p.
# p can be a constant.
def polynomial_parts(p):
    if isinstance(p,int):
        return [(p,1)]
    if isinstance(p,float):
        f = p>>Q_
        return polynomial_parts(f)
    if not is_polynomial(p):
        return [(p,1>>domain(p))]
    if isinstance(p,Symbol_type):
        return [(1,p._to_monomial_())]
    if isinstance(p,Monomial_type):
        s = p.cp()
        s._value = 1 >> p.domain().subdomain()
        return [(p._value, s)]
    CM =[]
    for i in p._value:
        CM+= polynomial_parts(i)
    return CM

## Synonyms
hpart = get_subpolynomial_degree
tdegree = total_degree = mdegree

# Monomial on a list or vector of expressions with given exponents  
def monomial(X,E):
    r = len(X)
    if len(E)<r:
        return 'Error: not enough exponents'
    m = 1
    for j in range(r):
        m *= X[j]**E[j]
    return m

def monomial_parts(m):
    if not is_polynomial(m): return 'monomial_parts: not a monomial'
    if isinstance(m,Poly_type): return 'monomial_parts: not a monomial'
    if isinstance(m,Symbol_type): return monomial_parts(m._to_monomial_())
    if len(m._variables)==0: return 'monomial_parts: not a monomial'
    E = m._index[:]
    V = variables(m)
    return (E,V)

def monomial_substitution(m,D):
    if not is_polynomial(m): return m
    if len(m._variables)==0: return m
    if isinstance(m,Poly_type): return 'monomial_substitution: not a monomial'
    E,V = monomial_parts(m)
    VS=[]
    for v in V:
        if v in D:
            VS.append(D[v])
        else:
            VS.append(v)
    return monomial(VS,E)

def polynomial_substitution(y,D):
    if not is_polynomial(y): return y
    if len(y._variables)==0: return y
    CM = polynomial_parts(y)
    x = 0
    for c,m in CM:
        x += c*monomial_substitution(m,D)
    return x

def cyclic_matrix(v,k):
    K = K_(v)
    if k < 0:
        return 'cyclic_matrix: number of rows must be non negative'
    if k == 0:
        return null_matrix(K)
    r = list(v)+(k-1)*[0>>K]
    C = [r]
    for _ in range(1,k):
        r = [0>>K] + r[:-1]
        C += [r]
    return matrix(C)
    
def vector_resultant(v,w):
    n = len(v)
    m = len(w)
    Cv = cyclic_matrix(v,m-1)
    Cw = cyclic_matrix(w,n-1)
    R = stack(Cv,Cw)
    return det(R)

#20200603

def simplex(cost,R):
    V = variables(cost)
    for r in R:
        V += variables(r)
    V = list(set(V))
    i = 0
    while(i < len(V) - 1):
        o = _order_symbols_(V[i],V[i+1])
        if o > 0:
            aux = V[i]
            V[i] = V[i+1]
            V[i+1] = aux  
            if (i > 0):
                i -= 1
            else:
                i += 1
            continue
        i += 1
    A = []
    b = []
    nc = len(R)
    for i in range(nc):
        L = []
        for v in V:
            cf = coeffs(R[i],v)
            if len(cf) > 1:
                L.append(float(cf[0]._value))
            else:
                L.append(float(0))
        L += [float(0)]*i + [float(1)] + [float(0)]*(nc-i-1)
        A.append(L)
        b.append(float(-zero_degree_term(R[i])._value))
    
    c = []   
    for v in V: 
        cfc = coeffs(cost,v)
        if len(cfc) > 1:
            c.append(float(cfc[0]._value))
        else:
            c.append(float(0))
    c += [float(0)]*nc + [float(-zero_degree_term(R[i])._value)]
        
#    A = matrix(A)
#    b = vector(b)
#    c = vector(c)
    
    tableau = _initialTableau_simplex_(c, A, b)

    while _canImprove_simplex_(tableau):
        pivot = _findPivotIndex_simples_(tableau)
        _pivotAbout_simplex_(tableau, pivot)

    t,s = tableau, _primalSolution_simplex_(tableau)
    sd = dict(s)
    
    sol = {}
    for i in range(len(V)):
        if i in sd:
            sol[V[i]] = sd[i]>>QQ(ZZ())
        else:
            sol[V[i]] = 0>>QQ(ZZ())
    o = evaluate(cost,sol)
    return o._value>>QQ(ZZ()),sol
constrained_maximum = simplex

def _initialTableau_simplex_(c, A, b):
    tableau = [row[:] + [x] for row, x in zip(A, b)]
    tableau.append(c[:] + [0])
    return tableau

def _canImprove_simplex_(tableau):
    lastRow = tableau[-1]
    return any(x > 0 for x in lastRow[:-1])

def _variableValueForPivotColumn_simplex_(tableau, column):
    pivotRow = [i for (i, x) in enumerate(column) if x == 1][0]
    return tableau[pivotRow][-1]

def _primalSolution_simplex_(tableau):
    # the pivot columns denote which variables are used
    columns = _transpose_tableau_simplex_(tableau)
    indices = [j for j, col in enumerate(columns[:-1]) if _isPivotCol_simplex_(col)]
    return [(colIndex, _variableValueForPivotColumn_simplex_(tableau, columns[colIndex])) for colIndex in indices]
 
def _objectiveValue_simplex_(tableau):
    return -(tableau[-1][-1])

def _isPivotCol_simplex_(col):
    return (len([c for c in col if c == 0]) == len(col) - 1) and sum(col) == 1

def _column_tableau_simplex_(A, j):
   return [row[j] for row in A]

def _transpose_tableau_simplex_(A):
    return [_column_tableau_simplex_(A, j) for j in range(len(A[0]))]

def _moreThanOneMin_simplex_(L):
    import heapq
    if len(L) <= 1:
        return False

    x,y = heapq.nsmallest(2, L, key=lambda x: x[1])
    return x == y

def _findPivotIndex_simples_(tableau):
    # pick minimum positive index of the last row
    column_choices = [(i,x) for (i,x) in enumerate(tableau[-1][:-1]) if x > 0]
    column = min(column_choices, key=lambda a: a[1])[0]

    # check if unbounded
    if all(row[column] <= 0 for row in tableau):
        raise Exception('Linear program is unbounded.')

    # check for degeneracy: more than one minimizer of the quotient
    quotients = [(i, r[-1] / r[column])
        for i,r in enumerate(tableau[:-1]) if r[column] > 0]

    if _moreThanOneMin_simplex_(quotients):
        raise Exception('Linear program is degenerate.')

    # pick row index minimizing the quotient
    row = min(quotients, key=lambda x: x[1])[0]

    return row, column

def _pivotAbout_simplex_(tableau, pivot):
    i,j = pivot

    pivotDenom = tableau[i][j]
    tableau[i] = [x / pivotDenom for x in tableau[i]]

    for k,row in enumerate(tableau):
        if k != i:
            pivotRowMultiple = [y * tableau[k][j] for y in tableau[i]]
            tableau[k] = [x - y for x,y in zip(tableau[k], pivotRowMultiple)]


def LP(n,d,R = []):
    if odd(d):
        n = n + 1
        d = d + 1
    K = [krawtchouk(k,n) for k in range(1,n//2 + 1)]
    A = [[evaluate(k,0)] + [evaluate(k,j) for j in range(d,n+1,2)] for k in K]
    A = -matrix(A)
    alpha = symbol('a',n)
    alpha = [1] + alpha[d - 1:n:2]
    alpha = vector(alpha,'c')
    RT = A*alpha
    RT = [r for r in RT] + R
    cost = sum(alpha)
#    show('cost:',cost)
#    show('RT:',RT)
    m1,m2 = simplex(cost,RT)
    b = floor(m1)
    if odd(b):
        beta = stack([1-1/(b>>Q_)],alpha[1:])
        RT = A*beta
        RT = [r for r in RT] + R
#        show('cost:',cost)
#        show('RT:',RT)
        m1,m2 = simplex(cost,RT)
    return floor(m1),m2
    
def floor(r):
    if isinstance(r,int) or isinstance(r,float):
        if r>= 0:
            return int(r)
        if r < 0:
            if int(r) == r:
                return int(r)
            else:
                return int(r) - 1
    if isinstance(r,QQ_type):
        return r.numerator()//r.denominator()





























