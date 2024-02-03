## wit_mor: defines the class of morphisms

from PyM.PyWIT.wit_var_sheaf import *

def pairing(p,T,v):
    C = coeffs_inc(p,v)
    I = C[0]
    CM = [(C[i],v**i) for i in range(1,len(C))]
    for (c,m) in CM:
        if m in T:
            v = T[m]
            I += c*v
    return I

def codimension(x,X):
    if hasattr(X,'_gcs'):
        if X._gcs != None:
            return wdeg(x,degs(X))
    return total_degree(x)

def vector_bundle(P):
    return P._vector_bundle
    
def fiber_dim(P):
    return P._fiber_dim
    
def base_dim(P):
    return P._base_dim

def lowerdata(P):
    return P._lowerdata

def lowerstar(f,x):
    return f.lowerstar(x)
    
def section(P):
    return P._section
    
def projective_bundle(X,E,h='h',y=False, name = ''):
    P = PROJ(X,E,h,y,name)
    h = P._gcs[0]
    return P,h

#Conics contained in a 5ic in P4 ---> 609250
#[G24,*c] = grassmannian(2,4,'c')
#S = tautological_bundle(G24)
#Q = symm(2,S.dual())
#P,z = projective_bundle(G24,Q,'z')
#S3 = symm(3,S.dual())
#S5 = symm(5,S.dual())
#F = S5 / (S3 * o_(-z,P.dim()))
#C = chern(F,rk(F))  
#c = lowerstar(P,C)
#N = integral(G24,c)
#show('N = ',N)

#Class PROJ (Projective bundle)
class PROJ:
    def __init__(self, X, E, h = 'h', y = False, name = ''):
        n = dim(X); r = rk(E); s = E.segre_vector()
        self._base_variety = X
        self._vector_bundle = E
        self._fiber_dim = r-1
        self._base_dim = n
        self._dim = n+r-1
        _,h = polynomial_ring(Q_,h)
        self._gcs = vector([h])
        self._degs = vector([1])
        LD = {h**(r-1):1}
        for i in range(1,n+1):
            LD[h**(i+r-1)]=s[i-1]
        self._lowerdata=LD
        self._kind = 'projective_bundle'
        self._section = h**(r-1)
        
        if pt(X)==None:
            self._pt = None
        else:
            self._pt = h**(r-1)*pt(X)
        if name == '' or name == None:
            self._name = None
        else:
            self._name = name
        
        if y:
            self._tan_bundle = E.dual()*o_(h)/O_()
            self._todd_class = self._tan_bundle.todd_vector()
        else:
            self._tan_bundle = None
            self._todd_class = None
            
    def gcs(self):
        return self._gcs
    
    def kind(self):
        return self._kind
    
    def dim(self):
        return self._dim
    vardim = dim
    
    def degs(self):
        return self._degs
    
    def pt(self):
        return self._pt
        
    def vector_bundle(self):
        return self._vector_bundle
        
    def fiber_dim(self):
        return self._fiber_dim
        
    def base_dim(self):
        return self._base_dim
    
    def lowerdata(self):
        return self._lowerdata
        
    def section(self):
        return self._section
    
    def tan_bundle(self):
        return self._tan_bundle
    
    def todd_class(self):
        return self._todd_class
    
    def lowerstar(self,p):
        return pairing(collect(p,self._gcs[0]),self._lowerdata,self._gcs[0])
    
    def __str__(self):
        if self._name == None:
            res = self._print_projective_bundle_()
            return res
        else:
            return self._name
     
    def __repr__(self):
        return self.__str__()
    
    def _print_projective_bundle_(self):
        res = 'Projective bundle of base dimension ' + self._base_dim.__str__() \
        + ', fiber dimension ' + self._fiber_dim.__str__() \
        + ', and hyperplane class ' + self._gcs[0]._print_()
        return res
   

#Class MOR (Morphism class)
def substitution(x,y):
    if len(x)>len(y): return 'substitution: too many entries x'
    return dict(zip(x,y))
    
def source(f): return f._source
def target(f): return f._target
def upperstardata(f): return f._upperstardata
def compose(g,f):
    V = upperstardata(g)
    D = [f.upperstar(v) for v in V.values()]
    return morphism(source(f),target(g),D)

def upperstar(f,y):
    return f.upperstar(y)

def pullback(f,F):
    X = F.ch()
    if (dim(f.source()) < len(X)):
        X = take(X,dim(f.source()))
    C = [upperstar(f,x) for x in X]
    return SH(rk(F),C)

def normal_bundle(f):
    F = pullback(f,tan_bundle(f.target()))
    G = tan_bundle(f.source())
    return F.quotient(G,dim(f.source()))

def double_points(f):
    N = normal_bundle(f)
    r = rk(N)
    return upperstar(f,lowerstar(f,1)) - chern(N,r)

def morphism(X,Y,e, name =''):
    return MOR(X,Y,e, name)

#upperstar(f,y): f.upperstardata(y)

class MOR:
    def __init__(self, X, Y, e, name = ''):
        self._source = X
        self._target = Y
        self._dim = dim(X)-dim(Y)
        self._upperstardata = substitution(gcs(Y),e)
        if name == '' or name == None:
            self._name = None
        else:
            self._name = name
    def source(self): return self._source
    def target(self): return self._target
    def dim(self): return self._dim
    def upperstardata(self): return self._upperstardata
    
    def upperstar(self,y):               ########################
        D = self._upperstardata
        return polynomial_substitution(y,D)                      ########################
    
    def lowerstar(self,x):
        X = self.source()
        Y = self.target()
        m = dim(Y); n = dim(X)
        k = codimension(x,X)
        if n-k > m: return 0
        if basis(Y) == None:
            additive_bases(Y)
        B1 = Y._basis[m-n+k]
        B2 = Y._dual_basis[m-n+k]
        return sum([integral(X,x*self.upperstar(B1[i]))*B2[i] for i in range(len(B1))])
    
    def __str__(self):
        if self._name == None:
            return self._print_morphism_()
        else:
            return self._name
     
    def __repr__(self):
        return self.__str__()
    
    def _print_morphism_(self):
        res = 'Morphism from variety of dimension ' + self._source._dim.__str__() \
        + ' to a space of dimension ' + self._target._dim.__str__() \
        + '; its dimension is ' + self._dim.__str__()
        return res

##Class bundle section and functions
def inclusion(B):
    return B._inclusion

def cl(B):
    return B._cl        

def bundle_section(X,F, name = ''):
    B = BSECT(X,F, name)
    i = B._inclusion
    return B,i
    
class BSECT(VAR):
    def __init__(self, X, F, name = ''):
        r = rk(F)
        dim = dim(X) - r
        T = {
            'kind': 'bundle_section',
            'gcs': gcs(X),
            'degs': degs(X)            
        }
        VAR.__init__(self,dim,T,name)
        i = morphism(self,X,gcs(X))
        i._kind = 'bundle_section_inclusion'
        self._inclusion = i
        self._cl = chern(F,r)
        self._tan_bundle = pullback(i,tan_bundle(X))/F
        Md = monomial_list(gcs(X),degs(X),dim)
        self._monomial_values = {x:integral(X,x*self._cl) for x in Md}
    
    def inclusion(self):
        return self._inclusion
    
    def cl(self):
        return self._cl

    

##Class Blowup and functions
def blowup_locus(W):
    return W._blowup_locus

def locus_inclusion(W):
    return W._locus_inclusion

def exceptional_class(W):
    return W._exceptional_class

def blowup_map(W):
    return W._blowup_map

def blowup(i,e='e',name=''):
    W = BLOWUP(i,e,name)
    return [W,exceptional_class(W)]

def blowup_points(k,X,e='e'):
    i = morphism(X,X,gcs(X))
    Y = BLOWUP(i,e,name)
    ce = exceptional_class(Y)
    [_,*ve] = PR(K_(ce),'e',k)
    E = vector(ve)
    L = o_(sum(E),X)
    T=trivial_bundle(dim(X),dim(X))
    Y._dim = dim(X)
    Y._blowup_locus = variety(0)
    Y._exceptional_class = sum(E)
    Y._gcs = splice(gcs(X),E)
    Y._degs = splice(gcs(X),vector(1,k))
    Y._pt = pt(X)
    Y._monomial_values = monomial_values(X)
    if (relations(X) != None) and (pt(X) != None):
        
        Y._relations = splice(relations(X),vector([x**dim(Y) + (-1)**dim(Y)*pt(Y) for x in E]))
        xe = []
        for x in gcs(X):
            for ee in E:
                xe.append(x*e)
        for i in range(k-1):
            for j in range(i+1,k):
                xe.append(E[i]*E[j])
        Y._relations = splice(Y._relations,vector(xe))
        Y._monomial_values = find_monomial_values(dim(Y),gcs(Y),degs(Y),relations(Y),pt(Y))
    else:
        Y._relations = None
    if tan_bundle(X) == None:
        Y._tan_bundle = None
    else:
        Y._tan_bundle = (tan_bundle(X) - T) + (T*L) + L - o_(0,X)
    f = morphism(Y,X,gcs(X))
    Y._blowup_map = f
    return [Y]+E.to_list()
    
class BLOWUP(VAR):
    def __init__(self, i, e = 'e', name = ''):
        X = target(i)
        Y = source(i)
        n = dim(X)
        m = dim(Y)
        
        [_,ve] = PR(K_(gcs(X)),e)
        gcsW = splice(ve,gcs(X))
        degsW = splice(1,degs(X))
        N = normal_bundle(i)
        V = monomial_values(X)
        for k in range(1,n):
            M = monomials(X)[n-k-1]
            for x in M:
                V[ve**k*x] = (-1)**(k-1)*integral(Y,segre(N,k-n+m,dim(Y))*upperstar(i,x))
        V[ve**n] = (-1)**(n-1)*integral(Y,segre(N,m,dim(Y)))
        
        T = {
            'kind': 'blowup',
            'gcs': gcsW,
            'degs': degsW,
            'pt': pt(X),
            'monomials': [monomial_list(gcsW,degsW,j) for j in range(1,n+1)],
            'monomial_values': V
            }
        
        VAR.__init__(self,n,T,name)
        f = morphism(self,X,gcs(X))
        self._blowup_locus = Y
        self._locus_inclusion = i
        self._exceptional_class = ve
        self._blowup_map = f
        
        def blowup_locus(self):
            return self._blowup_locus
        
        def locus_inclusion(self):
            return self._locus_inclusion
        
        def exceptional_class(self):
            return self._exceptional_class
        
        def blowup_map(self):
            return self._blowup_map


''' To do.
+ Ch = chern_character
+ Td = todd_character
+ compose(g,f) or compose(f,g)
+ pullback(f,F)
+ lowerstar(P,x)
+ lowerstar(f,x)
+ normal_bundle(f)

+ chi(X,F)
+ HRR(X,F)

- bundle_section(X, F, i)

- blowup(i,e,f)
- blowup_points(k,X,f,Y,e)

- double_points(f)

up(f:MOR):= make_current(source(f))
up(P:Projective_bundle):=make_current(P)
down(f:MOR):= make_current(target(f))
down(P:Projective_bundle):= make_current(base_variety(P))


'''




