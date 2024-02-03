from PyM import *
from code import Code

class ReedSoloWrapper(Code):
    def __init__(self, n=16, r=8, base=2, f_extension=4, 
                 phi_extension = 5):
        super().__init__(base ** f_extension, n,n-r-1)
        self.n = n
        self.r = r
        self.F_base = Zn(base)
        
        F_poly = get_irreducible_polynomial(self.F_base, f_extension)
        [self.F,b]= extension (self.F_base, F_poly, 'a')

        Phi_poly = get_irreducible_polynomial(self.F, base * (phi_extension - f_extension))
        [self.Phi,b]= extension (self.F, Phi_poly, 'a')
        
        self.a = [element(i,self.Phi) for i in range(1, n + 1)]
        self.h = self.a
        self.C = alternant_code(self.h, self.a, self.r, self.F)
        self.H = prune(blow(H_(self.C),self.F))
        if G_(self.C)!= None:
            self.G = G_(self.C)
        else: 
            self.G = left_kernel(transpose(self.H))

    def encode(self, u: np.array) -> np.array:
        # Generate random data of length k
        u = vector(self.F, list(u))
        return u * self.G

    def decode(self, y: np.array, e) -> np.array:
        
        erasure_locations = [i for i,x in enumerate(e) if x == 1]
        print("erasure_locations", erasure_locations)
        bad_H = submatrix(self.H, erasure_locations)
        no_error_locations = [i for i,x in enumerate(e) if x == 0]
        print("no_error_locations", no_error_locations)
        good_H = submatrix(self.H, no_error_locations)

        y_good = [y[i] for i in no_error_locations]
        print(type(y_good))
        print(vector(self.F, y_good))
        y_good = vector(self.F, y_good)
        syn = good_H * y_good
        missing_bits = solve_linear_system(bad_H, -syn)
        c = y
        for i, loc in enumerate(erasure_locations):
            c[loc] = missing_bits[i]
        x_ = c[:-self.r]
        return x_, c
    def get_rate(self) -> np.float32:
        """
        returns the calculation of the current codes rate
        """
        return (self.n - self.r) / self.n
    
    
C = ReedSoloWrapper()
print(len(C.G))
u = np.ones(C.k)
c = C.encode(u)
print(c)
print(len(c))
e = rd_error_vector(C.n, 3, C.F_base)

u_, c_ = C.decode(c, e)
print(u_)
print(c_)
print(c == c_)