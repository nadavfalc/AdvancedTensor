## Lie circle geometry

from PyM.PyWIT.PyECC.CC import *

# Utilities to cope with the symbolic use of floats

#EPS = 1/10**10
#def notnil(x): return abs(x)>=EPS
#def nil(x): return abs(x)<EPS
#def nonneg(x): return (notnil(x) and x>0) or nil(x)
#
## Q_type is the type of PyECC rational numbers
#def is_number(x):
#    if isinstance(x, (int, float, complex,QQ_type)) and not isinstance(x, bool): 
#        return True
#    else: return False
##
#def is_real(x):
#    if is_number(x) and not isinstance(x,complex): return True
#    else: return False

is_pair = ispair

# Lie vector of the circle with center M and radius r in the Euclidean plane.
# Points are encoded as circles of radius 0.
# For nonzero r, the orientation of the circle is encoded as the sign of r. 
# By default, r=0, so it gives the Lie vector of M, which by default is (0,0)
def Lie_vector(M=(0,0),R=0):
    a, b = M
    if is_pair(R):
        u,v = R
        d = a*u+b*v
        #from math import sqrt
        return [d,-d,a,b,sqrt(a**2+b**2)] 
    if is_real(R):
        u = (1+a**2+b**2-R**2)/2
        return [u,1-u,a,b,R]
    return 'Lie_vector: wrong parameters' 
#
LV = lie_vector=Lie_vector

def orientation(X):
    r = X[4]
    if r>0: return 1
    elif r<0: return -1
    else: return "orientation: object has no orientation"

# Lie metric for 5 and 4 vectors
def Lie_metric(X,Y):
    if len(X)==len(Y)==4: # make them 5-vectors
        X[4] = Y[4] = 0
    x1,x2,x3,x4,x5 = X
    y1,y2,y3,y4,y5 = Y
    lm = -x1*y1 + x2*y2 + x3*y3 + x4*y4 - x5*y5
    if nil(lm): return 0
    return lm
#
LM = lie_metric = Lie_metric


# Lie quadratic form for 5 or 4 vectors
def Lie_form(X): return Lie_metric(X,X)
#
LF = lie_form=Lie_form

def lie_gram_matrix(*S):
    M = [[Lie_metric(X,Y) for X in S] for Y in S]
    return M 


# Useful predicates
def is_lie_vector(X): 
    return len(X)==5 and nil(lie_form(X)) and notnil(dot(X,X))

def lie_type(X):
    if len(X)!=5 or nil(dot(X,X)): 
        return 'Lie_type: object must be a non-zero 5-vector'
    if notnil(Lie_form(X)): return -1
    elif nil(X[4]): return 0
    elif nil(X[0]+X[1]): return 1
    else: return 2
#
LT = Lie_type = lie_type
    
def is_point(X):
    if Lie_type(X)==0: return True
    else: return False

def is_line(X):
    if Lie_type(X)==1: return True
    else: return False

def is_circle(X): 
    if Lie_type(X)==2: return True
    else: return False

# Some basic functions 
def normalize(A):
    X = list(A[:])
    if not is_lie_vector(X): return 'normalize: object is not a Lie vector'
    if notnil(X[4]) and notnil(X[0]+X[1]):
        s = X[0]+X[1]
        return [x/s for x in X]
    else: return X
    
def reorient(A): 
    X = A[:]
    if not isinstance(X,list): return 'reorient: parameter is not a list'
    if len(X)!=5: return 'reorient: wrong length of object'
    X[4] = -X[4]
    return X
#
mirror = reorient


# Back from Lie objects to Euclidean objects

def circle(X):
    if is_circle(X): 
        Xn = normalize(X[:])
        return [(Xn[2],Xn[3]),Xn[4]]
    else: return 'circle: the object does not correspond to a circle'

def line(X):
    if is_line(X): 
        d = X[0]
        a,b = X[2:4]
        if notnil(a): return [(a,b),(d/a,0)]
        elif notnil(b): return [(a,b),(0,d/b)]
        else: return 'line: no non-zero orthogonal vector'
    else: return 'line: object is not a Lie line'

def point(X):
    if is_point(X):
        return (X[2],X[3])
    else: return 'point: this message should not occur'


# Lie angular metric
def Lie_angular_metric(X,Y):
    x1,x2,x3,x4,x5 = X
    if x5==0: return 'Lie_angular_metric: first object is not a circle'
    y1,y2,y3,y4,y5 = Y
    if y5==0: return 'Lie_angular_metric: second object is not a circle'
    x1=x1/x5; x2=x2/x5; x3=x3/x5; x4=x4/x5 
    y1=y1/y5; y2=y2/y5; y3=y3/y5; y4=y4/y5 
    lam = -x1*y1 + x2*y2 + x3*y3 + x4*y4
    if nil(lam): return 0
    else: return lam
# 
LAM = lie_angular_metric = Lie_angular_metric

def crossing_type(X,Y):
    x1,x2,x3,x4,x5 = X
    if x5==0: return 'crossing_type: first object is not a circle'
    y1,y2,y3,y4,y5 = Y
    if y5==0: return 'crossing_type: second object is not a circle'
    o = orientation(X)*orientation(Y)
    s = Lie_angular_metric(X,Y)
    if o==1:
        if s<-1: return 'external circles'
        elif s==-1: return 'external touching'
        elif s<1:
            if s==0: return 'orthogonal crossing'
            else: return 'crossing circles'
        elif s==1: return 'internal touching'
        elif s>1: return 'nested circles'
    else: 
        if s>1: return 'external circles'
        elif s==1: return 'external touching'
        elif s>-1:
            if s==0: return 'orthogonal crossing'
            else: return 'crossing circles'
        elif s==-1: return 'internal touching'
        elif s<-1: return 'nested circles'
    return 'crossing_type: could not decide'
        

def crossing_angle(X,Y):
    x1,x2,x3,x4,x5 = X
    if x5==0: return 'crossing_angle: first object is a point'
    y1,y2,y3,y4,y5 = Y
    if y5==0: return 'crossing_angle: second object is a point'
    from math import acos
    am = Lie_angular_metric(X,Y)
    if am<-1 or am>1:
        return 'crossing_angle: objects do not cross' 
    return acos(Lie_angular_metric(X,Y))
    

# Solving a quadratic homogeneous equation ax2+2bxy+cy2, given a,b,c
def solve_quadratic(a,b,c):
    if nil(a**2+b**2+c**2): 
        return 'solve_quadratic: at least one entry should be nonzero'
    if notnil(a):
        D = b**2-4*a*c
        from math import sqrt
        if nonneg(D):
            return [((-b+sqrt(D))/(2*a),1),((-b-sqrt(D))/(2*a),1)]
        else: 
            print('solve_quadratic: no real roots')
            return 0
    elif notnil(c): return [(1,0),(1,-2*b/c)]
    else: return [(0,1),(1,0)]


# Given three linearly independent 5-vectors, which may or many not
# belong to the Lie quadric, this function finds the intersections
# of the plane Lie-orthogonal to X, Y, Z
def Lie_section(A,B,C):
    X = A[:]; X[0]=-X[0]; X[4]=-X[4]
    Y = B[:]; Y[0]=-Y[0]; Y[4]=-Y[4]
    Z = C[:]; Z[0]=-Z[0]; Z[4]=-Z[4]
    import sympy
    M = sympy.Matrix([X,Y,Z])
    K = M.nullspace()
    #X = vector(X); Y = vector(Y); Z = vector(Z)
    #M = stack(X,Y,Z)
    #K = kernel(M)
    #if ncols(K)>2: return 'lie_section: Infinite solutions'
    #else: 
    #    v = K[:,0]
    #    w = K[:,1]
    
    if len(K)>2: return 'lie_section: Infinite solutions'
    else: [v,w] = K
    G = lie_gram_matrix(v,w) 
    st = solve_quadratic(G[0][0],2*G[0][1],G[1][1])
    if st==0: 
        print('lie_section: imaginary Lie objects')
        return 0
    s,t = st
    s1,s2 = s; t1,t2=t
    S = list(s1*v + s2*w); T = list(t1*v + t2*w)
    return normalize(S), normalize(T)
#
lie_section = Lie_section
   
### Examples
## Definition of the three circles as Lie vectors.
## Changing the sign of the second radius, or of the third, or of both
## yields the four pairs of solutions. 
#X = Lie_vector((6,4.25),3.5) 
#Y = Lie_vector((8,7),3.2)
#Z = Lie_vector((3.75,7.75),2.75)
#
## plotting
#close('all')
#
#fig1=plt.figure(figsize=(6,6))
#plt.axis('off')
#
#cfr(circle(X), lw=2)
#cfr(circle(Y), lw=2)
#cfr(circle(Z), lw=2)
#
#ST = lie_section(X,Y,Z)
#if ST == 0: print('no real circles of kind +++')
#else: 
#    S, T = ST
#    cfr(circle(S),lw=3,color='r')
#    cfr(circle(T),lw=3,color='b')
#plt.show()
#
## plotting
##close('all')
#
#fig2 = plt.figure(figsize=(6,6))
#plt.axis('off')
#
##show(X,Y,Z,Yb,Zb)
#
#Zb = mirror(Z)
#
#cfr(circle(X), lw=2)
#cfr(circle(Y), lw=2)
#cfr(circle(Zb), lw=2)
#
#ST = lie_section(X,Y,Zb)
#if ST == 0: print('no real circles of kind ++-')
#else: 
#    S, T = ST
#    cfr(circle(S),lw=3,color='r')
#    cfr(circle(T),lw=3,color='b')
#plt.show()
#
## plotting
##close('all')
#
#fig3 = plt.figure(figsize=(6,6))
#plt.axis('off')
#
#Yb = mirror(Y)
#cfr(circle(X), lw=2)
#cfr(circle(Yb), lw=2)
#cfr(circle(Z), lw=2)
#
#ST = lie_section(X,Yb,Z)
#if ST == 0: print('no real circles of kind +-+')
#else: 
#    S, T = ST
#    cfr(circle(S),lw=3,color='r')
#    cfr(circle(T),lw=3,color='b')
#plt.show()
#
#fig4 = plt.figure(figsize=(6,6))
#plt.axis('off')
#cfr(circle(X), lw=2)
#cfr(circle(Yb), lw=2)
#cfr(circle(Zb), lw=2)
#
#ST = lie_section(X,Yb,Zb)
#if ST == 0: print('no real circles of kind +--')
#else: 
#    S, T = ST
#    cfr(circle(S),lw=3,color='r')
#    cfr(circle(T),lw=3,color='b')
#plt.show()






