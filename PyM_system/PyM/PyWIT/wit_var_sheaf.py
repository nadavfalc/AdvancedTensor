# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:35:11 2019

@author: narcis
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:31:35 2018

@author: UPC-ESAII
"""

from PyM.PyWIT.wit_power import *
#import PyECC_globals as PGl

#########################
### 1. Sheaf constructors
#########################

# Sheaf of rank r, Chern character c, on variety X of dim d
def SH(r, c, d=None, name = None):
    if c==[]:
        ch = [0]
        if d == None:
            nm = d
        elif isinstance(d,VAR):
            ld = dim(d)
        else:
            nm = name
    else:
        if d == None:
            ld = len(c)
            nm = d
        elif isinstance(d,VAR):
            ld = dim(d)
            nm = name
        else:
            ld = d
            nm = name
        ch = c2p(p2c(c),ld)
    return Sheaf_type(r,ch,ld,nm)


# Sheaf of rank r, Chern class c, on variety X with dim d
def sheaf(r,c,d='', name = ''):
    if isinstance(d,str):
        ld = len(c)
        nm = d
    elif isinstance(d,VAR):
        ld = dim(d)
        nm = name
    else:
        ld = d
        nm = name
    cc = c2p(c,ld)
    return SH(r,cc,ld,nm)
    
## Bundle of rank n with Chern classes [c1,c2,...] on variety X
def bundle(n,c,d='',name=''):
    
    if isinstance(n,int)and(n>0):
        if isinstance(d,str):
            _,*vc = polynomial_ring(Q_,c,n)
            return sheaf(n,vc,len(vc),d)
        elif isinstance(d,VAR):
            _,*vc = polynomial_ring(Q_,c,min(n,dim(d)))
            return sheaf(n,vc,dim(d),name)
        else:
            _,*vc = polynomial_ring(Q_,c,n)
            return sheaf(n,vc,d,name)
    else:
        if isinstance(d,str):
            _,*vc = polynomial_ring(Q_,c,PGl.DIM)
            return sheaf(n,vc,PGl.DIM,d)
        elif isinstance(d,VAR):
            _,*vc = polynomial_ring(Q_,c,dim(d))
            return sheaf(n,vc,dim(d),name)
        else:
            _,*vc = polynomial_ring(Q_,c,d)
            return sheaf(n,vc,d,name)


#Trivial bundle of rank n
def trivial_bundle(n,d='',name=''):
    if isinstance(d,str):
        return sheaf(n,[0],n,d)
    if isinstance(d,VAR):
        return sheaf(n,[0],dim(d), name)
    else:
        return sheaf(n,[0],d,name)
    


# The line bundle o(d)
def o_(d,k='',name=''):
    if isinstance(k,str):
        return sheaf(1,[d],PGl.DIM,k)
    elif isinstance(k,VAR):
        return sheaf(1,[d],dim(k))
    else:
        return sheaf(1,[d],k,name)

## Sheaf O_X
def O_(d='', name=''):
    if isinstance(d,str):
        return o_(0,PGl.DIM,d)
    else:
        return o_(0,d,name)
    
# Adams opperations
def adams(k, x):
    if isinstance(x,list):
        return [k**i*x[i] for i in range(x)]
    if isinstance(x, vector):
        return vector([k**i*x[i] for i in range(x)])
    if isinstance(x, Sheaf_type):
        return SH(x._rank, [k**i*x._ch[i] for i in range(len(x._ch))])

def ch(F):
    return F._ch
    
def rk(F):
    return F._rank

def chern_character(F,T):
    return F.chern_character(T)

def change_ch_2_size(F,k):
    return F.change_ch_2_size(k)
    
def chern_vector(F, k = ''):
    return F.chern_vector(k)
    
def chern(F, k, d = None):
    return F.chern(k,d)

def segre_vector(F, d = None):
    return F.segre_vector(d)

def segre(F, k, d = None):
    return F.segre(k,d)
        
def todd_vector(F, d = None):
    return F.todd_vector(d)
    
def todd(F, k):
    return F.todd(k)
    
def todd_character(F, T):
    return F.todd_character(T)  
Td=todd_character 
    
def Hom(F,G, k = ''):
    return F.Hom(G,k)

# End functor
def End(F):
    return F.Hom(F)

hom=Hom

# Exterior powers of sheaf
def wedge(F, p, d = ''):
    return F.wedge(p,d)

def symm(k,E): return E.symm(k)
    
def koszul(F):
    return F.koszul()

#Class sheaf
# rank: non-negative integer
# ch: list or vector
# name: string
class Sheaf_type:
    def __init__(self, rank, ch, d_variety = None, name = None):
        if d_variety == None:
            self._d_var = PGl.DIM
        else:
            self._d_var = d_variety
        self._rank = rank
        if isinstance(ch, list):
            if len(ch) == 0:
                self._ch = vector([0])
            else:
                self._ch = vector(ch)
        else:
            self._ch = ch
        if name == None or name == '':
            self._name = 'SH(' + self._rank.__str__() + ', ' + self._ch._print_() + ')'
        else:
            self._name = name
     
    def __str__(self):
        return self._name
     
    def __repr__(self):
        return self.__str__()
    
    ######################
    ### 2. Basic functions
    ######################
    
    def rank(self):
        return self._rank
    
    def ch(self):
        return self._ch
    
    def change_ch_2_size(k):
        self._ch = chpad(self._ch,k)
    
    def chern_character(self, T):
        if isinstance(T,Poly_type):
            return self._rank + T*polynomial(self._ch,T)
        if isinstance(T,Symbol_type):
            return self._rank + T*polynomial(self._ch,T)
        if isinstance(T,Monomial_type):
            return self._rank + T*polynomial(self._ch,T)
        if isinstance(T,str):
            [P,TT] = polynomial_ring(K_(self._ch),T)
            return self._rank + TT*polynomial(self._ch,TT)
        else:
            return 'Wrong Parameter T'
    
    # Total chern class of sheaf, in vector form
    def chern_vector(self, k = ''):
        if k == '':
            return p2c(self._ch,self._rank)
        else:
            return p2c(self._ch,min(k,self._rank))
    
    # Particular Chern classes
    def chern(self, k, d = None):
        if d == None:
            d = len(self._ch)
        if isinstance(k,int):
            if k<0 or k>d:
                return 0
            elif k == 0:
                return 1
            else:
                return self.chern_vector()[k-1]
        elif isinstance(k,list):
            return [chern(j,d) for j in k]
        elif isinstance(k,Vector_type):
            return vector([chern(j,d) for j in k])
        elif isinstance(k,range):
            return [chern(j,d) for j in k]

    # Total Segre class
    def segre_vector(self, d = None):
        if d == None:
            d = len(self._ch)
        return invert_vector(self.chern_vector(d),d)

    # Particular Segre classes
    def segre(self, k, d = None):
        if d == None:
            d = len(self._ch)
        if isinstance(k,int):
            if k<0 or k>d:
                return 0
            elif k == 0:
                return 1
            else:
                return self.segre_vector()[k-1]
        elif isinstance(k,list):
            return [segre(j,d) for j in k]
        elif isinstance(k,Vector_type):
            return vector([segre(j,d) for j in k])
        elif isinstance(k,range):
            return [segre(j,d) for j in k]
        


    # Todd classes
    def todd_vector(self, d = None):
        if  d== None:
            d = len(self._ch)
        if d <= 0:
            return []
        a = pad(self._ch,d)
        t = c2p(todd_numbers(d))
        T = p2c([factorial(i+1)*a[i]*t[i] for i in range(d)])
        return T  #[value(factor(u)) for u in T]
    
    def todd(self, k):
        if k<0:
            return 0
        elif k == 0:
            return 1
        else:
            todd_vector(k)[k]
    
    def todd_character(self, T):
        if isinstance(T,Poly_type):
            return 1 + T*polynomial(self.todd_vector(),T)
        if isinstance(T,Monomial_type):
            return 1 + T*polynomial(self.todd_vector(),T)
        if isinstance(T,Symbol_type):
            return 1 + T*polynomial(self.todd_vector(),T)
        elif isinstance(T,str):
            [P,TT] = polynomial_ring(K_(self._ch),T)
            return 1 + TT*polynomial(self.todd_vector(),TT)
    
    ##############################
    ### 3. Operations with sheaves
    ##############################
        
    def dual(self):
        return SH(self._rank,witdual(self._ch))
    
    # Determinant bundle
    def determinant(self, d = ''):
        if d == '':
            ld = len(self._ch)
        else:
            ld = d
        return o_(self._ch[0],ld)

    def osum(self,G,k = ''):
        if k == '':
            lk = max(len(self._ch),len(G._ch))
        else:
            lk = k
        m = max(len(self._ch),len(G._ch))
        chS = chpad(self._ch,m)
        chG = chpad(G._ch,m)
        return SH(self._rank + G._rank,pad(chS + chG,lk))
        #return SH(self._rank + G._rank,pad(self._ch,min(k, max(len(self._ch),len(G._ch)))) + pad(G._ch,min(k, max(len(self._ch),len(G._ch)))))

    def __add__(self, G):
        if isinstance(G, Sheaf_type):
            return self.osum(G, max(len(self._ch), len(G._ch)))
        else:
            try:
                return G.__radd__(self)
            except:
                raise TypeError("I do not know how to opperate these elements")

    def __radd_(self, G):
        if isinstance(G, Sheaf_type):
            return G.osum(self, max(len(G._ch), len(self._ch)))
        else:
            try:
                return G.__add__(self)
            except:
                raise TypeError("I do not know how to opperate these elements")
        
    def odiff(self, G, k = ''):
        if k == '':
            lk = max(len(self._ch),len(G._ch))
        else:
            lk = k
        m = max(len(self._ch),len(G._ch))
        chS = chpad(self._ch,m)
        chG = chpad(G._ch,m)
        return SH(self._rank - G._rank,pad(chS - chG,lk))
        # return SH(self._rank - G._rank,pad(self._ch,min(lk, max(len(self._ch),len(G._ch)))) - pad(G._ch,min(lk, max(len(self._ch),len(G._ch)))))

    def __sub__(self, G):
        if isinstance(G, Sheaf_type):
            return self.odiff(G, max(len(self._ch), len(G._ch)))
        else:
            try:
                return G.__rsub__(self)
            except:
                raise TypeError("I do not know how to opperate these elements")

    def __rsub_(self, G):
        if isinstance(G, Sheaf_type):
            return G.odiff(self, max(len(G._ch), len(self._ch)))
        else:
            try:
                return G.__sub__(self)
            except:
                raise TypeError("I do not know how to opperate these elements")
    
    def quotient(self,G,k = ''):
        return self.odiff(G, k)

    def __floordiv__(self,G):
        return self.odiff(G)
    
    def __truediv__(self, G):
        return self.odiff(G)

    def tensor(self, G, k = ''):
        if k == '':
            lk = max(len(self._ch),len(G._ch))
        else:
            lk = k
        m = len(self._ch) + len(G._ch) - 1    
        rp = self._rank * G._rank
        
        F = self
        chF = c2p(p2c(F._ch),m)
        chG = c2p(p2c(G._ch),m)
        cnv = vector([convolution(chF,chG,i-2) for i in range(1,m+1)])
        return SH(rp, pad(G._rank * chF + F._rank * chG + cnv,lk))

    def __mul__(self, G):
        if isinstance(G, Sheaf_type):
            return self.tensor(G)
        else:
            try:
                return G.__rmul__(self)
            except:
                raise TypeError("I do not know how to opperate these elements")
    
    def __rmul__(self, G):
        if isinstance(G, Sheaf_type):
            return G.tensor(self)
        else:
            try:
                return G.__mul__(self)
            except:
                raise TypeError("I do not know how to opperate these elements")

    #Hom functor
    def Hom(self,G, k = ''):
        F = self.dual()
        return F.tensor(G,k)

    # End functor
    def End(self):
        return self.Hom(F)

    hom=Hom

    # Exterior powers of sheaf
    def wedge(self, p, d = ''):
        r = self._rank
        if p < 0 or p > r:
            return sheaf(0,[0],d)
        elif p == 0:
            return sheaf(1,[0],d)
        elif p == 1:
            return SH(self._rank, self._ch)
        elif p == r:
            return self.determinant(d)
        elif p > r-p:
            E = self.wedge(r-p,d)
            return self.determinant(d)*E.dual()
        else:
            E = adams(p,self)
            W = (-1)**(p+1)*E._rank
            for m in range(1,p-1):
                E = self.wedge(m)*adams(p-m,self)
                W = W + (-1)**(p-m+1) * E._ch
            return SH(binomial(r,p),W/p)
    
    def koszul(self):
        K = self._ch
        n = len(K)
        for j in range(2, self._rank):
            E = self.wedge(j,n)
            K = K + (-1)**j * E._ch
        return SH(0,K)

    # Symmetric powers of a sheaf (I)
    def symm(self, p):
        if isinstance(p, int):
            if p < 0:
                return sheaf(0,[0])
            if p == 0:
                return sheaf(1,[0])
            if p == 1:
                return self
            r = self._rank
            E = self * self.symm(p-1)
            S = E._ch
            for k in range(2, min(r,p)+1):
                E = self.wedge(k) * self.symm(p-k) 
                S = S + (-1)**(k+1) * chpad(E._ch,len(S))
            return SH(binom(r+p-1,p),S)
#        if isinstance(p, Base_type):
#            O = sheaf(1,[0], self._rank)
#            if p == 0:
#                return O
#            if p == 1:
#                return self
#            r = self._rank
#            s = self.segre(range(1,r))
#            s = [1]+s
#            DIM = DIM + r - 1
#            E = self.dual() * o_()
    
#t=todd_vector( (dual(E) * o_(a)) / O )
#  t= SH(1,t) * o_(p*a)
#  t=ch\t
#  t=1+T*polynomial(t,T)
#  DIM=DIM-r+1
#  S=sigma (s.(i+1)*coefficient(t,a,r-1+i)) with i in 0..DIM
#  S=S/T^(r-1)
#  S=dense_coefficient_vector(S,T)
#  SH(S.1,tail(S))
#


##Variety

#####################################
### 1. Constructing general varieties 
#####################################

## VAR constructors
#  Variety denoted X, of dim n (n=DIM by default)
def variety(n, Table = {}, name = ''):
    return VAR(n,Table,name)

### Construction of Different type of Varieties


# To define a new algebraic_curve, of genus g and
# with P as the class of a point. Announce new current variety
# if y=true and new variety if y=false
def algebraic_curve(g, P, name = ''):
    [_,P] = polynomial_ring(Q_,P)
    n = 1
    d = {
        'kind': 'algebraic_curve',
        'gcs': vector([P]),
        'degs': vector([1]),
        'relations': [P*P],
        'pt': P,
        'monomial_values': {P: 1},
        'tan_bundle': SH(1,vector([collect((2-2*g)*P,P)])), ##FALTA FER BÉ AIXÒ [collect((2-2*g)*P,{P})]
        'todd_class': vector([collect((1-g)*P,P)])
    }
    return [VAR(1,d,name),P]




# # To define a new projective space of dim n
# # and hyperplane class named h. Proclaims the
# # construction, plainly if y=false, making the
# # space current if y=true
def projective_space(n, h, name = ''):
    [_,h] = polynomial_ring(Q_,h)
    S = sheaf(n,vector([binom(n+1,j)*h**j for j in range(1,n+1)]),n)
    if n == 1:
        kd = 'projective_line'
    elif n == 2:
        kd = 'projective_plane'
    else:
        kd = 'projective_space of dimension ' + n.__str__()    
    d = {
        'kind': kd,
        'gcs': vector([h]),
        'degs': vector([1]),
        'relations': [h**(n+1)],
        'monomials': [[h**j] for j in range(1,n+1)],
        'monomial_values': {h**n : 1},
        'pt': h**n,
        'tan_bundle': S,
        'todd_class': S.todd_vector(n)
    }
    return [VAR(n,d,name),h] 

# Defines the grassmanninan of k-planes
# in n dimensional projective space, with
# chern classes c.1,c2,... for the tautological bundle.
# Proclaims the construction, plainly if y=false,
# and making it current if y=true
    
def grassmannian(k, n, c = 'c', name = ''):
    dim = (k+1)*(n-k)
    S = bundle(k+1, c, dim)
    T = trivial_bundle(n+1,dim)
    Q = T/S.dual()
    #show("Q = ",Q)
    C = S.chern_vector()
    M = [monomial_list(C,range(1,k+2),j) for j in range(1,dim+1)]
    R = take(invert_vector(C,n+1),-k-1)#take(invert_vector(C,dim), -k) 
    
    pt = C[k]**(n-k)
    d = {
        'kind': 'grassmannian', 
        'tautological_bundle': S,
        'tautological_quotient': Q,
        'gcs': C,
        'degs': vector([l for l in range(1,k+2)]),
        'relations': R,
        'pt': pt,
        'monomials': M,
        'monomial_values': find_monomial_values((k+1)*(n-k),C,list(range(1,k+2)),R,pt),
        'tan_bundle': Q.hom(S.dual()),
        'todd_class': Q.hom(S.dual()).todd_vector()
        
    }
    return [VAR(dim,d,name)]+C._value[:]

def canonical_class(X):
    return -X._tan_bundle.chern(1)
    
def gcs(V):
    return V._gcs

def kind(V):
    return V._kind

def dim(V):
    return V._dim

def degs(V):
    return V._degs

def monomials(V):
    return V._monomials

def monomial_values(V):
    return V._monomial_values

def monomial_value(V,m): 
    M = monomial_values(V)
    if m in M:
        return M[m]
    else:
        return '{0} undefined in {1}'.format(m,'variety')

def pt(V):
    return V._pt

def relations(V):
    return V._relations
    
def tan_bundle(V):
    return V._tan_bundle

def todd_class(V):
    return V._todd_class

def basis(V):
    return V._basis

def dual_basis(V):
    return V._dual_basis

def tautological_bundle(V):
    return V._tautological_bundle

def tautological_quotient(V):
    return V._tautological_quocient
    
'''
The function find_monomial_values is the cornerstone for the
the proper functioning of the Integral function. It assigns numerical
values to all monomials of degree d in the variables h={h.1,...,h.m},
with weights w={w.1,...,w.m}, assuming the relations
R={R.1,...,R.n} (each relation is a homogeneous polynomial
in the variables h, taking into account the weights).
The numbers assigned are normalized so that the monomial P
representing a point has value 1.
'''

def find_monomial_values(d, h, w, R, P):

#  m  number of variables h
#  n  number of relations R
#  M  list of monomials up to degree d
# Md  list of monomials of degree d
#  N  number of monomials of degree d
# Rd  relations in degree d generated by R
#  r  the length of Rd
#  i  to hold a temporary degree
#  x  identifier stem for variables x1,...,xN
#  X  the list {x1,...,xN}
#  S  substitution of Md.j by xj
# eq  linearization via S of the monomials of Rd
#  A  matrix of coeffs of linearized rels
#  k  kernel of A
#  pt variety point 

    if isinstance(h,list):
        h = vector(h)

    m=len(h)
    n=len(R)

    #show("h =",h)
    #show("w =",w)

    # Monomials of degree up to d, grouped by degrees
    M=[monomial_list(h,w,j) for j in range(1,d+1)]
    
    Md=M[d-1]; N=len(Md)
    #show("Md=",Md)
    #show('M = ',M)

    # Find all relations in degree d generated by R
    #show("R =",R)
    Rd=[]
    for p in R:
        p.add_variables(h)
        i=wdeg(p,w)
        if i == 0:
            continue;
        if i==d:
            Rd = Rd + [p]
        if i<d:
            Rd += [p*q for q in M[d-i-1]]
    r = len(Rd)
    #show("Rd=",Rd)

    # Case with no relations in top degree
    if r==0:
        #if N==1 then return({{Md.1,1}})
        if N==1: 
            return {Md[0]:1}
        else:
             return 'No relations in top degree were found'

    # and so relations become a list of linear equations
    A = []
    E = partitions(d,w)
    #show("E =",E)
    for p in Rd:
        p.add_variables(h)
        D = p.to_dict()
        #show(p)
        #show(D)
        a = []
        for e in E:
            if tuple(e) in D:
                a += [D[tuple(e)]]
            else:
                a += [0]
        A += [a]
    A = matrix(A)
    #show("A=",A)
    
    #show('R = ',R)

    # Solve the system (find a basis of ker A)
    K=kernel(A)
    K = transpose(K)[0]
    #show("K=",K)
    # Get value of point
    DP = P.to_dict()
    j = E.index(list(list(DP.keys())[0]))
    pt=K[j]
    #show("pt=",pt)
    K = K / pt
    
    # Assign values to all monomials
    #{(Md.j=>(k.j/pt)) with j in range(Md)}
    MV = dict(zip(Md,K))
    #show('MV = ',MV)
    return MV

#Class VAR
class VAR:
    def __init__(self, n, T, name = ''):
        self._dim = n
        if 'gcs' in T:
            self._gcs = vector(T['gcs'])
        else:
            self._gcs = null_vector()
        if 'kind' in T:
            self._kind = T['kind']
        else:
            self._kind = None
        if 'degs' in T:
            self._degs = T['degs']
        else:
            self._degs = None
        if 'monomials' in T:
            self._monomials = T['monomials']
        else:
            self._monomials = []
        if 'monomial_values' in T:
            self._monomial_values = T['monomial_values']
        else:
            self._monomial_values = null_vector()
        if 'pt' in T:
            self._pt = T['pt']
        else:
            self._pt = None
        if 'relations' in T:
            self._relations = T['relations']
        else:
            self._relations = None
        if 'tan_bundle' in T:
            self._tan_bundle = T['tan_bundle']
        else:
            self._tan_bundle = None
        if 'todd_class' in T:
            self._todd_class  = T['todd_class']
        else:
            self._todd_class  = null_vector()
        if 'basis' in T:
            self._basis = T['basis']
        else:
            self._basis = None
        if 'dual_basis' in T:
            self._dual_basis = T['dual_basis']
        else:
            self._dual_basis = None
        if 'tautological_bundle' in T:
            self._tautological_bundle = T['tautological_bundle']
        else:
            self._tautological_bundle = None
        if 'tautological_quotient' in T:
            self._tautological_quocient = T['tautological_quotient']
        else:
            self._tautological_quocient = None
        if name == '' or name == None:
            self._name = None
        else:
            self._name = name
            
    def gcs(self):
        return self._gcs
    
    def kind(self):
        return self._kind
    
    def dim(self):
        return self._dim
    
    def degs(self):
        return self._degs
    
    def monomials(self):
        return self._monomials
    
    def monomial_values(self):
        return self._monomial_values
    
    def tautological_bundle(self):
        return self._tautological_bundle
    
    def tautological_quotient(self):
        return self._tautological_quocient
    
    def pt(self):
        return self._pt
    
    def relations(self):
        return self._relations
        
    def tan_bundle(self):
        return self._tan_bundle
    
    def todd_class(self):
        return self._todd_class
    
    def basis(self):
        return self._basis
    
    def dual_basis(self):
        return self._dual_basis
    
    def todd_character(self, T):
        if self._todd_class == None:
            return 0
        if isinstance(T,Poly_type):
            return 1 + T*polynomial(self._todd_class,T)
        if isinstance(T,Symbol_type):
            return 1 + T*polynomial(self._todd_class,T)
        if isinstance(T,Monomial_type):
            return 1 + T*polynomial(self._todd_class,T)
        elif isinstance(T,str):
            [P,TT] = polynomial_ring(K_(self._todd_class),T)
            return 1 + TT*polynomial(self._todd_class,TT)
    
    def __str__(self):
        if self._name == None:
            res = self._print_variety_()
            return res
        else:
            return self._name
     
    def __repr__(self):
        return self.__str__()
    
    def _print_variety_(self):
        res = 'Variety of dimension ' + self._dim.__str__()
        if self._gcs != null_vector():
            res += ', gcs = ' + self._gcs._print_()
        if self._kind != None:
            res += ', kind = '+ self._kind
            
        if self._degs != None:
            res += ', degs = ' +  self._degs.__str__()
        
        if self._monomials != []:
            res += ', monomials = ' + self._monomials.__str__()
        
        if self._monomial_values != null_vector():
            res += ', monomials values = ' + self._monomial_values.__str__()
        
        if self._pt != None:
            res += ', pt = ' + self._pt.__str__()
        
        if self._relations != None:
            res += ', relations = ' + self._relations.__str__()
        
        if self._tan_bundle != None:
            res += ', tan budle = ' + self._tan_bundle.__str__()
        
        if self._todd_class != null_vector():
            res += ', todd class = ' + self._todd_class.__str__()
        
        if self._basis != None:
            res += ', basis = ' + self._basis.__str__()
        
        if self._dual_basis != None:
            res += ', dual basis = ' + self._dual_basis.__str__()
        return res
        
  
## Constucts the cartesian product of two varieties
def cartesian_product(X,Y, name = ''):
    dx = X._dim
    dy = Y._dim
    d = dx + dy
    Dc = {}
    Dc['gcs'] = splice(X._gcs,Y._gcs)
    Dc['degs'] = splice(X._degs, Y._degs)
    if (X._monomial_values != null_vector()) and (Y._monomial_values != null_vector()):
        Dc['monomial_values'] = {m1*m2: X._monomial_values[m1]*Y._monomial_values[m2] for m1 in X._monomial_values[m1].keys() for m2 in Y._monomial_values[m2].keys()}
    if (X._pt != None)and(Y._pt != None):
        Dc['pt'] = X._pt*Y._pt
    if (X._tan_bundle != None) and (Y._tan_bundle != None):
        Dc['tan_bundle'] = osum(X._tan_bundle, Y._tan_bundle)
    if (X._basis != None) and (Y._basis != None):
        B1 = X._basis
        b1 = len(B1)
        B2 = Y._basis
        b2 = len(B2)
        T = []
        for k in range(d+1):
            BB = []
            for i in range(1,min(k+2,b1+1)):
                ii = k + 2 - i
                if ii > b1:
                    continue
            BB += [m1*m2 for m1 in B1[i] for m2 in B2[ii]]
            T += [BB]
        Dc['basis'] = T
    if (X._dual_basis != None) and (Y._dual_basis != None):
        B1 = X._dual_basis
        b1 = len(B1)
        B2 = Y._dual_basis
        b2 = len(B2)
        T = []
        for k in range(d+1):
            BB = []
            for i in range(1,min(k+2,b1+1)):
                ii = k + 2 - i
                if ii > b1:
                    continue
            BB += [m1*m2 for m1 in B1[i] for m2 in B2[ii]]
            T += [BB]
        Dc['dual_basis'] = T
        
    return VAR(d,Dc,name)

def betti(X,n = None):
    if n == None:
        [betti(X,i) for i in range(X._dim + 1)]
    else:
        return len(X._basis[i])

def additive_basis(X):
    d = X._dim
    m = d//2
    k = 0
    
    XB = []
    XDB = []
    #show(X._monomials)
    while k <= d-k:
        if k == 0:
            M1 = [1]
            M2 = X._monomials[d-1]
        else:
            M1 = X._monomials[k-1]
            if k < d:
                M2 = X._monomials[d-k-1]
            else:
                M2 = [1]
        
        A = []
        for m1 in M1:
            a = []
            for m2 in M2:
                if m1*m2 in X._monomial_values:
                    a.append(X._monomial_values[m1*m2])
                else:
                    a.append(0)
            A.append(a)
        A = matrix(A)
                    
        (B,I,J) = domain(A).rank_block_ext(A)
        B = transpose(1/B)
        M1 = [M1[i] for i in I]
        V = vector([M2[j] for j in J])
        M2 = B*V
        if isinstance(M2, Vector_type):
            M2 = M2.to_list()
        else:
            M2 = [M2]
        XB += [M1]
        XDB += [M2]
        k += 1
    if odd(d):
        m=m+1
    for k in range(d//2 + 1):
        XB += [XDB[m-k]]
        XDB += [XB[m-k]]
    X._basis = XB
    X._dual_basis = XDB         

def _integral(p,T):
    I = 0
    CM = polynomial_parts(p)
    for (c,m) in CM:
        if m in T:
            v = T[m]
            I += c*v
    return I

def integral(V,p):
    return _integral(p,monomial_values(V))      
                
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
