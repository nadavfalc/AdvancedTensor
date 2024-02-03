## Factor by the Pollard rho algorithm

from PyM.PyWIT.PyECC.arithmetic.primes import *
from PyM.PyWIT.PyECC.arithmetic.cdi import *

def pollard(n):
        if(n%2==0):
            return 2;
        x=rd.randrange(2,1000000)
        c=rd.randrange(2,1000000)
        y=x
        d=1
        while(d==1):
            x=(x*x+c)%n
            y=(y*y+c)%n
            y=(y*y+c)%n
            d=igcd(x-y,n)
            if(d==n):
                break;
        return d

# Get the prime factors of a number
def prime_factors(x):
    if rabin_miller(x): return [x]
    P=[]
    while x>1:
        d = pollard(x)
        if d == x: continue
        P += prime_factors(d)
        x //= d
        if rabin_miller(x): 
            P+=[x]
            break
    return P

def list_2_dict(P):
    D = dict()
    S = set(P)
    for p in S:
        D.update({p:P.count(p)})
    return D

def ifactor(n): return list_2_dict(prime_factors(n))      
    
def divisors(n):
    F = ifactor(n)
    D=[1]
    
    for p in F:
        L = len(D)
        for j in range(L):
            dj = D[j]
            for i in range(F[p]):
                dj *= p
                D.append(dj)
    D.sort()
    return D

def phi_euler(n):
    F=ifactor(n)
    a = 1
    for p in F:
        a = a * p**(F[p]-1)*(p-1)
    return a

def order(a,n):
    if igcd(a,n) != 1: return 'order is not defined' 
    D = divisors(phi_euler(n)) 
    for d in D:
        if power(a,d,n) == 1: 
            return d

def tau(n):
    return len(divisors(n))

def sigma(n):
    return sum(divisors(n))

def sigmax(n):
    return sigma(n)-n

def mu_moebius(n):
    F=ifactor(n)
    if any(e>1 for e in F.values()):
        return 0
    r = len(F)%2
    if r: return -1
    else: return 1

def lambda_carmichael(n):
    if n == 1: return 1
    F=ifactor(n)
    P=sorted(list(F.keys()))
    L = []
    if 2 in P:
        e = F[2]
        if e<=2: L += [2**(e-1)]
        else: L += [2**(e-2)]
        P = P[1:]
    for p in P: 
        L += [(p-1)*p**(F[p]-1)]
    return ilcm(*L)    
        
# cyclotomic classes
def cyclotomic_class(j,n,q=2):
    if igcd(q,n) != 1: return 'cyclotomic class not defined'
    j = j%n
    C = [j]
    p = (j*q)%n
    while p != j:
        C += [p]
        p = (p*q)%n
    return C


def cyclotomic_classes(n,q=2):
    if igcd(q,n) != 1: return 'cyclotomic classes not defined'
    L = [[0]]
    J = list(range(1,n))
    while J != []:
        k = J[0]
        C = cyclotomic_class(k,n,q)
        L += [C]
        J = list(set(J)-set(C))
    return L



    
    