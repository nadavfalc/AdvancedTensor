## Integer routines
from functools import reduce
import random as rd

def a2b(a,b,d=1):
    if b < a: return list(range(a,b-1,d))
    return list(range(a,b+1,d))

def rd_int(a=0, b=1):
    return rd.randint(a,b)
#random integer of n decimal digits
def ri(n): return rd_int(10**(n-1),10**n)

def remainder(m,n): return m % n
#rem = remainder

def dlen(n):
    from math import log10
    return int(log10(n))+1

def blen(n):
    from math import log2
    return int(log2(n))+1

def floor(r): 
    if r>= 0:
        return int(r)
    if r < 0:
        if int(r) == r:
            return int(r)
        else:
            return int(r) - 1

def ceil(r):
    if r>= 0:
        if int(r)==r: return int(r)
        else: return int(r)+1
    else:
        return int(r)
ceiling = ceil

def frac(x): return x - floor(x)

# factorial n
def factorial(n):
    if n == 0 or n == 1: return 1
    f = 1
    for k in range(2,n+1):
        f *= k
    return f

# binomial (n,m)
def binom(n,m):
    if m<0 or m>n: return 0
    if m>n//2+1:
        m=n-m
    b =1
    for j in range(m):
        b*=(n-j)
    for j in range(2,m+1):
        b=b//j
    return b

# N-th Fibonacci number
def fib(N):
    f0, f1 = 1, 1
    if N==0 or N == 1: return 1
    for _ in range(2,N+1):
        f0, f1 = f1, f0+f1
    return f1
fibonacci = fib

# N-th Lucas numbers U_N(P,Q), V_N(P,Q)
def lucas_number(P,Q,x0,x1,N):
    if N==0: return x0
    if N==1: return x1
    for _ in range(2,N+1):
        x0, x1 = x1, P*x1-Q*x0
    return x1
def lucas_U(P,Q,N): return lucas_number(P,Q,0,1,N)
def lucas_V(P,Q,N): return lucas_number(P,Q,2,P,N)
def lucas(N): return lucas_V(1,-1,N)
def pell(N): return lucas_U(2,-1,N)
def pell_lucas(N): return lucas_V(2,-1,N)
def jacobsthal(N): return lucas_U(1,-2,N)
def jacobsthal_lucas(N): return lucas_V(1,-2,N)
def mersenne(N): return 2**N-1  # lucas_U(3,2,N)

#greatest common divisor of a sequence of integers N
def igcd(*N):
    def GCD(m,n):
        r0, r1 = m, n
        while r1 != 0:
            r0, r1 = r1, r0%r1
        return r0
    return abs(reduce(GCD,N))

#least common multiple of a sequence of integers N
def ilcm(*N):
    def LCM(m,n): return m*n//igcd(m,n)
    return abs(reduce(LCM,N))

# x**d mod n
def power(x,d,n):
    p = 1
    if x == 1: return p
    x = x%n
    while d >= 1:
        if d%2: 
            p = (p*x)%n
        x = (x**2)%n
        d = d//2
    return p


# Jacobi symbol (a/n), --Crandall-Pomerance, p. 98
def jacobi(a,n):
    # n has to be a positive odd integer, while a is any integer
    a = a%n
    s = 1
    while a!=0:
        while not a&1:
            a //= 2
            if n%8 in (3,5):
                s = -s
        (a,n) = (n,a)
        if a%4 == n%4 == 3:
            s = -s
        a = a%n
    if n==1: return s
    return 0
    
'''        
# Jacobi symbol (a/n)
def jacobi_old(a,n):
    # n has to be a positive odd integer, while a is any integer
    if n<=0 or not n&1 : return 'denominator is not a positive odd integer' 
    if abs(a)>n: a = a%n
    if igcd(a,n) != 1: return 0
    j = 1
    # reduction to the case a odd
    while a%2==0:
        j *= j2(n)
        a = a//2
    # case a = 1
    if a==1: return j
    # use reciprocity
    return j*rec(a,n)*jacobi(n,a)
# Jacobi symbol (2/n)
def j2(n):
    n = n%8
    if n == 1 or n == 7: return 1
    return -1
# reciprocity sign
def rec(m,n):
    if m%4==1 or n%4==1: return 1
    return -1    

def test(N,f):
    P = [2]
    for n in range(3,N,2):
        if f(n): P += [n]
    return P
'''

def is_square(n):
  x = n // 2
  R = {x}
  while x * x != n:
    x = (x + (n // x)) // 2
    if x in R: return False
    R.add(x)
  return True
  
def extended_euclidean_algorithm(a, b):
   if abs(b) > abs(a):
      (x,y,d) = extended_euclidean_algorithm(b, a)
      return (y,x,d)
 
   if abs(b) == 0:
      return (1, 0, a)
 
   x1, x2, y1, y2 = 0, 1, 1, 0
   while abs(b) > 0:
      q, r = divmod(a,b)
      x = x2 - q*x1
      y = y2 - q*y1
      a, b, x2, x1, y2, y1 = b, r, x1, x, y1, y
 
   return (x2, y2, a)
