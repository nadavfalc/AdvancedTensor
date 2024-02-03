## wit-power

from PyM.PyWIT.PyECC.CC import *

def inclast(p):
    p[-1] += 1
    return p

def pad(c,r):
    if r <= len(c):
        return c[:r]
    if isinstance(c,list): return c+(r-len(c))*[0]
    return vector(pad(list(c),r))
    
def witdual(c):
    if len(c)==0: return c[:]
    d = c[:]
    for i in range(0,len(c),2):
        d[i] = -c[i]
    return d

def vprod(x,y,d=''):
    if d=='': d = len(x)+len(y)
    if d==0: return []
    z = d*[0]
    #print(z)
    for j in range(min(d,len(y))): z[j]=y[j]
    #print(z)
    for j in range(min(d,len(x))):
        z[j] += x[j]
        for k in range(min(d-1-j,len(y))):
            z[j+k+1] += x[j]*y[k]
        #print(z)
    return vector(z)
    
# x = [1,2,3]
# y = [7,8,9,10]

def vpower(x,m,d=''):
    if d=='': d=len(x)
    if m==0: return []
    if m==1: return pad(x,d)
    z = m*[0]
    while m>0:
        if m%2: z = vprod(z,x,d)
        x = vprod(x,x,d)
        m = m//2
    return vector(z)

# Given list or vector [c0,c1,...,c(r-1)], computes list or vector  
# [s0,s1,..,s(r-1)s such that 1+s0*t+s1*t**2+...+s(r-1)t**r is inverse of
# 1+c0*t+c1*t**2+...+c(r-1)t**r mod t**(r+1)
def invert_vector(c, r=''): 
    if isinstance(c,Vector_type): return vector(invert_vector(list(c),r))
    if len(c)==0: return invert_vector([0],r)
    if r==0: return []
    if r=='': r=len(c)
    c = pad(c,r)
    s = [-c[0]]
    #if r==1: return s
    for k in range(1,r):
        s += [-c[k]-convolution(s,c,k-1)]
    return s 
 
# Taylor coeffs of x/(1-e^{-x}) at 0
# 1/2, 1/12, 0, -1/720, 0, 1/30240, 0, -1/1209600, 
# 0, 1/47900160, 0, -691/1307674368000, 0, 1/74724249600, 0, 
# -3617/10670622842880000, 0, 43867/5109094217170944000, 0, 
# -174611/802857662698291200000
def todd_numbers(n):
    u = 1>>Q_
    return invert_vector(witdual([u/factorial(j+1) for j in range(1,n+1)]))

def bernoulli_numbers(N):
    u = 1>>Q_
    if N==0: return [u]
    B = invert_vector([u/factorial(j+1) for j in range(1,N+1)])
    return [u]+[factorial(k+1)*B[k] for k in range(N)]
BV_ = bernoulli_numbers

def bernoulli_number(N):
    if N==0: return 1>>Q_
    elif N==1: return -1/2>>Q_
    else:
        if odd(N): return 0>>Q_
        else: return bernoulli_numbers(N)[-1]
B_=bernoulli_number

def bernoulli_polynomial(m,x='x'):
    u = 1>>Q_
    [Qx,x] = polynomial_ring(Q_,x)
    r = m//2
    B = bernoulli_numbers(r)
    P = (x**(m+1)+(u*(m+1)/2)*x**m)/(m+1)
    S = sum([(-1)**(k-1)*binom(m+1,2*k)*B[k-1]/(m+1)*x**(m+1-2*k) for k in range(1,r+1)])
    return P+S
    
def high_bernoulli_numbers(N,k):
    u = 1>>Q_
    if N==0: return [u]
    B = invert_vector([u/factorial(j+1) for j in range(1,N+1)])
    B = vpower(B,k,N)
    return [u]+[factorial(j+1)*B[j] for j in range(N)]
#    
HBs_=high_bernoulli_numbers

#B_5_4=-9 // ha funcionat!

# B(N,k) 
def high_bernoulli_number(N,k):
    if N==0: return 1>>Q_
    return high_bernoulli_numbers(N,k)[-1]
#
HB_=high_bernoulli_number

def zhe(n,j):
    k=0
    for _ in range(j):
        k=2**k
        print(k)
    #print(k)
    return HB_(n,k)/(2**(2*n))

def qfloat(x):
    return x.numerator()/x.denominator()


## Functions s2n, n2s, c2p, p2c

# For the elementary symmetric polynomials (s) and the Newton sums (n) 
# of the variables [x1,x2,...], s2n expresses the n's in
# terms of the s's and n2s goes the other way around
# See Fulton, p.56. Example: s2n([s1,s2])=[s1,s1^2-2s2]
# n2s([n1,n2])=[n1,(n1^2-n2)/2]

def s2n(s, r=''):
    #if isinstance(s,Vector_type): return vector(s2n(list(s),r))
    if r == '': r = len(s)
    if s==[] or r==0: return []
    n = [s[0]]
    if r==1: return n
    x = [t for t in s]
    x = pad(witdual(x),r)
    for k in range(1,r):
        nk = -(convolution(x,n,k-1) + (k+1)*x[k])
        n += [nk]
    return vector(n)

def n2s(n,r=''):
    #if isinstance(n,Vector_type): return vector(n2s(list(n),r))
    if r == '' or r>len(n): r = len(n)
    if n==[] or r==0: return []
    y = [t for t in n]
    s = [-y[0]]
    if r==1: return witdual(s)
    for k in range(1,r):
        sk = -(convolution(s,y,k-1) + y[k])/(k+1)
        s += [sk]
    return vector(witdual(s))


# If c=[c1,c2,...] and [n1,n2,...] are the elementary
# symmetric polynomials and the Newton sums of the
# variables [x1,x2,...], c2p maps c to p=[n1/1!,n2/2!,...]
# and p2c goes the other way around
#c2p(c:Vector):= c2p(c,length(c));

def c2p(c,r=''):
    #if isinstance(c,Vector_type): return vector(c2p(list(c),r))
    n = s2n(c,r)
    if isinstance(n,list):
        n = vector(n)
    if K_(n) == ZZ():
        return vector([QQ_type(Q_,n[k],factorial(k+1)) for k in range(len(n))])
    return vector([n[k]/factorial(k+1) for k in range(len(n))])

def p2c(p,r=''):
    #if isinstance(p,Vector_type): return vector(p2c(list(p),r))
    return n2s([factorial(k+1)*p[k] for k in range(len(p))],r)


# Examples


#[P,x1,x2,x3] = polynomial_ring(Z_,'x1','x2','x3')
#x = [x1,x2,x3]
#s= [s1,s2,s3]=[x1+x2+x3,x1*x2+x1*x3+x2*x3,x1*x2*x3]
#n = s2n(s); show('n =',n)

'''
s= [s1,s2,s3]=[x1+x2+x3,x1*x2+x1*x3+x2*x3,x1*x2*x3]
c=vector(s)

show('c =',c)
n = s2n(s); show('n =',n)
s1=n2s(n); show('s1 =',s1)  # s1=s
p = c2p(c); show('p =',p)
c1=p2c(p); show('c1 =',c1)  # s1=s
x = [x1,x2,x3]
y = s2n(x); show('y =',y)
show(n2s(y)); show('x =',x)
p = c2p(x); show('p =',p)
show(p2c(p))

'''
# Monomial on a list or vector of expressions with given exponents  
def monomial(X,E):
    r = len(X)
    if len(E)<r:
        return 'Error: not enough exponents'
    m = 1
    for j in range(r):
        m *= X[j]**E[j]
    return m
    
# Example
# show(monomial([x1,x2,x3],[5,7,11]))    

# Given a list or vector d of positive interger, and an integer n,
# find the list of tuples of non-negative integers [m1,...,mr]
# such that  m1*d1+...+mr*dr=n # (partitions of n with weights
# d1,...,dr)   
def partitions(n,d):
    L = []
    r = len(d)
    if n < min(d): return []
    if r==1:
        if n % d[0] == 0:
            return [[n//d[0]]]
        else: return []
    l = partitions(n,d[:r-1])
    if len(l) > 0:
        L = [p+[0] for p in l]
    l = partitions(n-d[-1],d)
    if len(l) > 0:
        l = [inclast(p) for p in l]
        L += l
    if d[-1]==n:
        L += [(r-1)*[0] + [1]]
    return L

# Example:
# partitions(5,[1,2,3])
# => [[5, 0, 0], [3, 1, 0], [1, 2, 0], [2, 0, 1], [0, 1, 1]]

    
    
# Given variables x={x1,x2,...,xr} of degrees
# d={d1,d2,...,dr}, get the list of monomials
# x1^m1 ... xr^mr such that m1*d1+...+mr*dr=n.
# The list of the lists {m1,...,mr} is returned by
# partitions(d,n).
def monomial_list(X,d,n):
    E = partitions(n,d)
    return [monomial(X,e) for e in E]

# Example
#  vector(monomial_list(vector(x),[1,2,3],5))
# => [x1**5, x1**3*x2, x1*x2**2, x1**2*x3, x2*x3] 
#     :: Vector[Q[x1][x2][x3]]


# Elementary symmetric polynomials of the expressions x
# of degree k
def symmetric_polynomial(x,k):
    n = len(x)
    if n==0 or k>n or k<0: return 0
    if k==0: return 1
    if k==1: return sum(x)
    y = x[1:]
    return x[0]*symmetric_polynomial(y,k-1)+symmetric_polynomial(y,k)
    
def symmetric_polynomials(x,K=''):
    if K=='': K = range(1+len(x))
    if isinstance(K,int): K = range(K+1)
    return vector([symmetric_polynomial(x,k) for k in K])

# Examples
# symmetric_polynomials(x)
# => [1, (x3 + x2 + x1), ((x2 + x1)*x3 + x1*x2), x1*x2*x3]
# symmetric_polynomials(x,[1,3])
# => [(x3 + x2 + x1), x1*x2*x3]


def newton_sum(x,k):
    if k<0: return 'Error: k must be non-negative'
    if k==0: return len(x)
    return sum([t**k for t in x])

def newton_sums(x,K):
    if isinstance(K,int): 
        if K<0: return 'Error: K cannot be a negative integer'
        K = range(K+1)
    return vector([newton_sum(x,k) for k in K])

def stirling_numbers(n):
    [Kx,x] = polynomial_ring(Z_,'x')
    f=1>>Kx
    for k in range(n):
        f*=x-k
    return reverse(coeffs(f))
    
def stirling_number_1st(n,k):
    if n<0 or k>n or k<0: return 'Erroneous parameters' 
    return stirling_numbers(n)[k]

def stirling_number_2nd(n,k):
    s = [(-1)**(k-j)*binom(k,j)*j**n for j in range(k+1)]
    return sum(s)//factorial(k)

def wdeg(p,r):
    if not is_polynomial(p):
        return 0
    x = variables(p)
    Tx = p.domain().variable(['Tx'])
    pp = evaluate(p,x,[x[i]**r[i]*Tx**r[i] for i in range(len(x))])
    if not is_polynomial(pp):
        return 0
    return degree(pp,'Tx')

def chpad(v,r):
    return c2p(p2c(v),r)
    

'''
'''    
