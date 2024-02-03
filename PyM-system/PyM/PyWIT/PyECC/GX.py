## GX : tools for geometric illustration

import numpy as np
#import scipy.optimize as opt
import matplotlib.pyplot as plt

# synonyms
cos = np.cos; sin = np.sin; 
acos = np.arccos; asin = np.arcsin; atan = np.arctan
pi = np.pi; sqrt = np.sqrt
L = np.log2
array=np.array; dot = np.dot
ls = np.linspace
plot = plt.plot
close = plt.close
display=plt.show
exp = np.exp

## Constructive Geometry


# a point is a pair
def ispair(X):
    return isinstance(X, tuple) and len(X)==2

# degs to rads and viceversa
def rad(x): return x*pi/180
def deg(t): return t*180/pi

# points on segments
def point(A,B,t): return ((1-t)*A[0]+t*B[0],(1-t)*A[1]+t*B[1]) 
def mid_point(A,B): return point(A,B,1/2)

# point from polar coordinates r, t
def ptGX(O,r,t): return (O[0]+r*cos(t*pi/180),O[1]+r*sin(t*pi/180))

# the argument of a point
def arg(X): 
    x = X[0]; y = X[1]
    if x**2+y**2 == 0: return 'arg: the argument of '+str(X)+ ' is not defined' 
    if x == 0:
        if y>0: return pi/2
        else: return 3*pi/2
    if y == 0:
        if x>0: return 0
        else: return pi
    if y>0:
        if x>0: return atan(y/x)
        else: return pi-atan(-y/x) 
    if y<0:
        if x<0: return pi+atan(y/x)
        else: return 2*pi-atan(-y/x)

# distance of two points
def distance(A,B): return sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)

# If K = [O,r] is the circle of radius r, centered at the origin O,
# and X is a point on K, this function returs the end point Y
# of the chord XY of length h, in the counterclockwise sense
def end_chord(K,X,h):
    r = K[1]
    a = arg(X)
    t = 2*asin(h/(2*r))
    return (r*cos(t+a), r*sin(t+a))
    
def normGX(v): return sqrt(dot(v,v))
def versor(v): 
    l = normGX(v)
    if l == 0: return 'versor: the versor of v = '+str(v)+' is not defined'
    return v/l
     
def vectorGX(A,B): return array([B[0]-A[0],B[1]-A[1]])

def proj(X,A,B):
    u = versor(vectorGX(A,B))
    x = vectorGX(A,X)
    p=dot(x,u)*u
    return (A[0]+p[0],A[1]+p[1])

def parallel_proj(X,v,A,B):
    w = vectorGX(A,B)
    M = np.transpose(np.vstack([v,w]))
    x = vectorGX(A,X)
    (_,t) = np.linalg.solve(M,x)
    w = t*w
    return (A[0]+w[0], A[1]+w[1])

def intersect(A,B,X,Y): return parallel_proj(A,vectorGX(A,B),X,Y)

def walk(A,B,d):
    u = versor(vectorGX(A,B))
    v = d*u
    return (A[0]+v[0],A[1]+v[1])

def area(A,B,C):
    f = lambda X,Y: X[0]*Y[1]-X[1]*Y[0]
    return (f(A,B)+f(B,C)+f(C,A))/2
def angle(A,B,C):
    (b1,b2) = (B[0]-A[0],B[1]-A[1])
    (c1,c2) = (C[0]-A[0],C[1]-A[1])
    nb = sqrt(b1*b1+b2*b2)
    nc = sqrt(c1*c1+c2*c2)
    return acos((b1*c1+b2*c2)/(nb*nc))*180/pi


    
 ## Constructive graphics   

def bullet(A, size=6, color='k', shape='o'): 
    return plot([A[0]],[A[1]], markersize=size, color=color, marker=shape)
    
def bullets(*P,size=6, color='k', shape='o'):
    for p in P:
        bullet(p,size=size, color=color, shape=shape)
    show()

def seg(A,B, lw=1, dashing='-', color='k'): 
    return plot([A[0],B[0]],[A[1],B[1]],linewidth=lw, linestyle=dashing,color=color)
def segs(*S,lw=1,dashing='-', color='k'):
    for s in S:
        seg(s[0],s[1],lw=lw, dashing=dashing, color=color)

def ruler(A,B,lw=1, dashing='-', color='k', left_inc=10, right_inc=10):
    u = versor(vectorGX(A,B))
    d = distance(A,B)
    incl = left_inc; incr = right_inc
    vl = -(d*left_inc/100)*u
    vr = d*(1+right_inc/100)*u
    return seg((A[0]+vl[0],A[1]+vl[1]),(A[0]+vr[0],A[1]+vr[1]), lw=lw, color=color, dashing=dashing)
    

def cfr(K, lw=1, dashing='-', color='k'):
    t = ls(0,2*pi,360)
    a = K[0][0]; b = K[0][1]; r = K[1]
    return plot(a+r*cos(t), b+r*sin(t), linestyle=dashing, linewidth=lw, color=color)
    
## labeling facilities

def lable(P,Text,fs=18, dx=0, dy=0, color='k'):
    return plt.text(P[0]+dx, P[1]+dy, Text, fontsize=fs, color=color)

def North(X,Text, fs=18, color='k'):
    return lable(X, Text, dx=-0.03, dy = 0.05, fs=fs, color=color)
def NNE(X,Text, fs=18, color='k'):
    return lable(X, Text, dy = 0.04, fs=fs, color=color)
def NE(X,Text, fs=18, color='k'):
    return lable(X, Text, dx=0.03, dy = 0.035, fs=fs, color=color)
def Est(X,Text, fs=18, color='k'):
    return lable(X, Text, dx = 0.03, dy = -0.02, fs=fs, color=color)
def SE(X,Text, fs=18, color='k'):
    return lable(X, Text, dy = -0.1, fs=fs, color=color)
def South(X,Text, fs=18, color='k'):
    return lable(X, Text, dx = -0.03, dy = -0.12, fs=fs, color=color)
def SW(X,Text, fs=18, color='k'):
    return lable(X, Text, dx = -0.12, dy = -0.08, fs=fs, color=color)
def West(X,Text, fs=18, color='k'):
    return lable(X, Text, dx = -0.12, dy = -0.02, fs=fs, color=color)
def NNW(X,Text, fs=18, color='k'):
    return lable(X, Text, dx=-0.1, dy = 0.03, fs=fs, color=color)
    