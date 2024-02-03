# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 18:43:14 2020

@author: narcis
"""

from PyM.PyWIT.wit_chow import *
from PyM.PyWIT.wit_lie import *

def resultant(f,g, vars = None):
    if vars == None:
        V = variables(product(variables(f)+variables(g)))
        if len(V) == 2:
            x,y = V
        elif len(V)==1:
            v = coeffs_inc(f)
            w = coeffs_inc(g)
            return vector_resultant(v,w)
        elif len(V) < 1:
            return 0
        else:
            x = V[-2]
            y = V[-1]
    else:
        if len(vars) == 1:
            x = vars[0]
            v = coeffs_inc(f,x)
            w = coeffs_inc(g,x)
            return vector_resultant(v,w)
        else:
            x = vars[-2]
            y = vars[-1]
    m = degree(f,y); n=degree(g,y)
    v = coeffs_inc(f,y); w = coeffs_inc(g,y)
    return vector_resultant(v,w)

def discriminant(f): return resultant(f,der(f))

def _vector_resultant_ext(v,w):
    n = len(v)
    m = len(w)
    Cv = cyclic_matrix(v,m-1)
    Cw = cyclic_matrix(w,n-1)
    R = stack(Cv,Cw)
    return R

def resultant_ext(f,g,vars = None):
    if vars == None:
        V = variables(product(variables(f)+variables(g)))
        if len(V) == 2:
            x,y = V
        elif len(V)==1:
            v = coeffs_inc(f)
            w = coeffs_inc(g)
            R = _vector_resultant_ext(v,w)
            Lf = [cofactor(R,i,0) for i in range(degree(g))]
            Lg = [cofactor(R,i + degree(g),0) for i in range(degree(f))]
            a = polynomial(Lf,V[0])
            b = polynomial(Lg,V[0])
            return det(R),a,b
        elif len(V) < 1:
            return 0
        else:
            x = V[-2]
            y = V[-1]
    else:
        if len(vars) == 1:
            x = vars[0]
            v = coeffs_inc(f,x)
            w = coeffs_inc(g,x)
            R = _vector_resultant_ext(v,w)
            Lf = [cofactor(R,i,0) for i in range(degree(g,x))]
            Lg = [cofactor(R,i + degree(g,x),0) for i in range(degree(f,x))]
            a = polynomial(Lf,x)
            b = polynomial(Lg,x)
            return det(R),a,b
        else:
            x = vars[-2]
            y = vars[-1]
    m = degree(f,y); n=degree(g,y)
    v = coeffs_inc(f,y); w = coeffs_inc(g,y)
    R = _vector_resultant_ext(v,w)
    Lf = [cofactor(R,i,0) for i in range(n)]
    Lg = [cofactor(R,i + n,0) for i in range(m)]
    a = polynomial(Lf,y)
    b = polynomial(Lg,y)
    return det(R),a,b

def cosubmatrix(A,i,j):
    m,n = shape(A)
    if i<0 or j<0 or i>= n or j>=n: return 'cosubmatrix: wrong indices'
    I = list(range(0,i))+list(range(i+1,m))
    J = list(range(0,j))+list(range(j+1,n))
    return A[I,J]
#    M = stack(A[:i,:],A[i+1:,:])
#    M = splice(M[:,:j],M[:,j+1:])
#    return M

    
def cofactor(A,i,j): 
    m,n = shape(A)
    if m!=n: return 'cofactor: matrix is not square'
    if i<0 or j<0 or i>= n or j>=n: return 'cofactor: wrong indices'
    return (-1)**(i+j)*det(cosubmatrix(A,i,j))
    
'''
The function imult(f,g,O) computes the intersection
multiplicity of the plane curves f(x,y)=0 and g(x,y)=0
at the point O =[ a,b], which by default is the origin O = [0,0].
The algoritm is extracted from Fulton's book Algebraic Curves 
(see http://www.math.lsa.umich.edu/~wfulton/CurveBook.pdf,
Section 3.3, including the Example on page 40).
The auxiliary functions used in the body of imult are
included below, and afterwords we work out Fulton's example.
'''

def imult(f,g,O=[0,0]):
    if f==0 or g==0: return 'Infinity'
    if len(O)==3:
        a,b,c = O
        if c!=0: O = [a/c,b/c]
        else: return imultinfinity(f,g,O)
    # if evaluate(f,O)!=0 or evaluate(g,O)!=0: return 0
    # V = variables(product(variables(f)+variables(g)))
    # print(V)
    # if len(V) == 2:
    #     x,y = V
    # elif len(V)==1:
    #     return 'Infinity'
    # elif len(V) < 1:
    #     return 0
    # elif len(V) > 2:
    #     return 'imult: too many variables'
    # else:
    #     x,y = V
    [_,x,y] = PR(K_(f),'x','y')
    if evaluate(f,[x,y],O)!=0 or evaluate(g,[x,y],O)!=0: return 0
    a,b = O
    if a!=0 or b!=0:
        f = evaluate(f,[x,y],[x+a,y+b])
        g = evaluate(g,[x,y],[x+a,y+b])
    f = f + x*y - x*y
    g = g + x*y - x*y
    f0 = constant_coeff(f,y); r = degree(f0)
    g0 = constant_coeff(g,y); s = degree(g0)
    if f0==0:
        if g0==0:
            return 'Infinity'
        else:
            return trailing_degree(g0) + imult(f/y,g)
    else: # f0!=0
        if g0==0:
            return trailing_degree(f0)+imult(f,g/y)
        else: # g0!=0
            c0 = leading_coeff(f0)
            d0 = leading_coeff(g0)
            if r<=s:
                return imult(f,c0*g-d0*x**(s-r)*f)
            else:
                return imult(d0*f-c0*x**(r-s)*g,g)
    return 'error'

# If f is a polynomial, the following function homogenizes it
# with respect to a specified new variable 
def homogenize(f,z=''):
    if z=='': _,z=polynomial_ring(Q_,'z')
    n = total_degree(f)
    F = 0
    for k in range(n+1):
        h = hpart(f,k)
        if h!=0: F += h*z**(n-k)
    return collect(F,z)
    
# _,x,y=polynomial_ring(Q_,'x','y')
# f = -x**3-x**2+y**2
# F = homogenize(f)
# show("F =",F)

# def dehomogenize(F,z):
#     # subs_element(F,1,z) o subs_element(F,z,1) no va
#     return evaluate(F,z,1)

def imultinfinity(f,g,O):
    VF = variables(f)
    VG = variables(g)
    if (VF == VG):
        x = VF[0]
        y = VF[1]
    else:
        if len(VF) == 2:
            x = VF[0]
            y = VF[1]  
        else:
            x = VG[0]
            y = VG[1]
    a,b,_ = O
    _,z = polynomial_ring(Q_,'z')
    F = homogenize(f,z); G= homogenize(g,z)
    #show(F)
    #show(G)
    if b!=0:
        a = a/b; b = 1
        f = evaluate(F,[x,y,z],[x,1,y])
        #show(f)
        g = evaluate(G,[x,y,z],[x,1,y])
        #show(g)
        return imult(f,g,[a,0])
    a=1; b=0
    f = evaluate(F,[x,y,z],[1,x,y])
    g = evaluate(G,[x,y,z],[1,x,y])
    return imult(f,g,[0,0])    


# # Test. The difference should be x**k
# [Qx,x] = polynomial_ring(Q_,'x')
# k=120
# bp=bernoulli_polynomial(k)
# cp=evaluate(bp,x-1)
# show(bp-cp)


# F = Q_   # in this case it yields 14
# F = Zn(5) # in this case it yields 18
# [Fxy,x,y] = bivariate_polynomial_ring(F,'x','y','Fxy')
# 
# f = (x**2+y**2)**3-4*x**2*y**2
# g = (x**2+y**2)**2 +3*x**2*y -y**3
# show(imult(f,g))
