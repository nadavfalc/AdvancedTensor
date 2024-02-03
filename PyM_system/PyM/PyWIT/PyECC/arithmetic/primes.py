## Prime routines

from PyM.PyWIT.PyECC.arithmetic.ZZ import *

import random as rd

def rd_int(a=0, b=1):
    return rd.randint(a,b)

# list of small primes
small_primes =  [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
'''
small_primes =  [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,
223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,
349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,
479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,
619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,
769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,
929,937,941,947,953,967,971,977,983,991,997] 
'''

# v. Crandall-Pomerance-2005, Algorithm 3.5.6
# Rabin-Miller primality test
def rabin_miller(n, k=25):
    #if not isinstance(n,int):
    #    return 'parameter n is not an integer'
    if n<0: n = -n
    if 0 <= n <= 1: return False
    if n in small_primes: return True
    for p in small_primes:
        if n%p == 0: return False
    d = n-1
    s = 0
    while d & 1 == 0:
        s += 1; d = d//2 
    for _ in range(k):
        a = rd_int(2,n-1)
        if not sprp(n,a): return False
    return True

# Crandall-Pomerance-2005, algorithm 3.5.13
# n must be odd > 1
def miller(n):
    from math import log, e
    #gamma = 0.5772156649015329
    def ln(x): return log(x,e)
    #c = e**gamma/ln(2)
    w = min(ceil(2*(ln(n))**2),n-1)
    for a in a2b(2,w):
        if not sprp(n,a): return False
    return True

def is_prime(n,method='rabin_miller'):
    if method == 'rabin_miller': return rabin_miller(n)
    if method == 'BPSW': return BPSW(n)
    if method == 'miller': return miller(n)

def next_prime(x):
    while not rabin_miller(x):
        x += 1
    return x
    
# Select at random a prime number of n decimal digits
def rd_prime(n):
    while True:
        r = ri(n)
        if even(r): r += 1
        if is_prime(r): return r

# Picking a pair of distinct n-digit prime numbers    
def rsa_pair(n,m=3):
    while True:
        p = rd_prime(n); q = rd_prime(n)
        if igcd(m,p-1)>1 or igcd(m,q-1)>1 or p==q: 
            continue
        else: break
    return (p,q)

# def rd_prime(n=5):
#     return next_prime(ri(n))

# v. Crandall-Pomerance-2005, Algorithm 3.5.2
def strong_probable_prime(n,a):
    d = n-1; s = 0
    while d & 1 == 0:
        s += 1; d = d//2 
    x = power(a,d,n)
    if x == 1 or x == n-1: return True
    for r in range(1,s):
        x = x**2 % n
        if x == n-1: return True
    return False
#
sprp = strong_probable_prime

# v. Crandall-Pomerance-2005, Definition 3.5.3
def strong_pseudo_prime(n,a):
    return sprp(n,a) and not is_prime(n)
spsp = strong_pseudo_prime

## Lucas_based tests

'''
Selfridge specified the following parameters for the generation of the Lucas sequences: 
P = 1 and Q = (1 - D)/4, where D is the first integer in the sequence 
{5, -7, 9, -11, 13, -15, ...} for which GCD(D,N)=1 and the Jacobi symbol (D|N) = -1. 
Note, however, that if N is a perfect square, no such D will exist, and the search for D would continue all the way to ±sqrt(N), at which point GCD(D,N)=sqrt(N) would expose N as composite. Consequently, the algorithm also presumes the presence of a preliminary check for perfect squares (as well as even integers and integers < 3). 

From https://en.wikipedia.org/wiki/Baillie-PSW_primality_test and link https://oeis.org/A217255:
The first ten strong Lucas pseudoprimes (with Lucas parameters P = 1, Q = -1) are
5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519, 75077, 97439, 100127, 113573, 115639, 130139, 155819, 158399, 161027, 162133, 176399, 176471, 189419, 192509, 197801, 224369, 230691, 231703, 243629, 253259, 268349, 288919, 313499, 324899 (sequence A217255 in OEIS). 


def lucas_prp(n,a,b):
    if n<2: return 'lucas_prp: n has to be >1'
    D = a**2-4*b
    if is_square(D): return 'lucas_prp: D is a square'
    print(n, 2*a*b*D)
    if igcd(n,2*a*b*D)!=1: return 'lucas_prp: n and 2*a*b*D are not relatively prime' 
    return D
    
#P2 = [n for n in range(3,91000,2) if spsp(n,2)]
#P3 = [n for n in range(3,100000,2) if spsp(n,3)]
#P5 = [n for n in range(3,80000,2) if spsp(n,5)]
'''

# Test Baillie-Pomerance-Selfridge-Wagstaff
def BPSW(n):
    if n in small_primes: return True
    for p in small_primes:
        if n%p == 0: return False
    if not sprp(n,2): return False
    s = 1; k = 2
    while True:
        D = s*(2*k+1)
        if jacobi(D,n) == -1: break
        s = -s; k +=1
    P = 1; Q = (1-D)//4
    (v0,v1)=lucas_chain_V(P,Q,n+1,n)
    U = (2*v1-P*v0)%n
    if U: return False
    return True
baillie_pomerance_selfridge_wagstaff = BPSW

def lucas_chain_V(P,Q,m,n):
    mb = bin(m)[2:]
    v0 = 2; v1 = P
    j = 0
    for b in mb:
        Qj = power(Q,j,n)
        if b == '1':
            v0,v1 = (v0*v1)%n-(Qj*P)%n, (v1**2-2*Qj*Q)%n
        else:
            v0,v1 = (v0**2-2*Qj)%n,(v0*v1)%n-(Qj*P)%n 
        j = 2*j+int(b)
    return (v0,v1)

def lprp(n):
    if n%2==0 or n<2: return False
    if is_square(n): return False
    s = 1; k = 2
    while True:
        D = s*(2*k+1)
        if jacobi(D,n) == -1: break
        s = -s; k +=1
    P = 1; Q = (1-D)//4
    (v0,v1)=lucas_chain_V(P,Q,n+1,n)
    U = (2*v1-P*v0)%n
    if U: return False
    return True
    
def lpsp(n): return lprp(n) and not is_prime(n)
    
# Tests whether the b-th Mersenne number 2**n-1 is prime
def lucas_lehmer(n):
    if n==2: return True
    if not is_prime(n): return False
    N = 2**n-1
    S = 4
    for _ in range(2,n):
        S = (S*S-2) % N
    if S==0: return True
    return False
#
LL = lucas_lehmer


            