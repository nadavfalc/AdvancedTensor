from PyM import *
from code import Code


def create_dict_of_field(F: Zn, f_size: int) -> dict:
    F_map = {}
    for i in range(f_size):
        F_map[element(i, F)] = i
    return F_map

def map_binary_vector_to_extension_field(x: np.array):
    decimal_number = 0
    for i, bit in enumerate(x):
        decimal_number += bit << (len(x) - i - 1)
    return decimal_number    
    
def local_primitive_root(K):
    q = cardinal(K); n = q-1
    print("n is", n)
    for j in range(2,q):
        el = element(j,K)
        ord_el = order(el)
        print("in index ", j, "order", ord_el)
        if ord_el==n:
            return element(j,K)
    return "primitive_root error: root not found"
    
def local_BCH_Code(alpha,d,l=1,K=''):
    print("kaki-1")
    n = order(alpha)
    print("order of alpha is", n)
    h = geometric_series(alpha**l,n)
    print("kaki1")
    a = h
    print("kaki2")
    return AC(h,a,d-1,K)

def extend_field_by_m(q, m):
    Fq = Zn(q)
    F_poly = get_irreducible_polynomial(Fq, m)
    print("F poly= ", F_poly)
    [F,X]= extension (Fq, F_poly, 'a')
    return F, X

class ReedSoloWrapper(Code):
    def __init__(self, n=15, r=8, base=2, f_extension=4, 
                 phi_extension = 5):
        super().__init__(base ** f_extension, n,n-r-1)
        self.n = n
        self.r = r
        self.k = self.n - self.r
        self.F_base = Zn(base)
        if f_extension == 1:
            self.F = self.F_base
            self.F_size = base
        else:
            F_poly = get_irreducible_polynomial(self.F_base, f_extension)
            [self.F,b]= extension (self.F_base, F_poly, 'a')
            self.F_size = base ** f_extension
        print("order of prime in F is", order(local_primitive_root(self.F)))
        
        self.create_systematic_G()
        
        """Phi_poly = get_irreducible_polynomial(self.F, base * (phi_extension - f_extension))
        [self.Phi,b]= extension (self.F, Phi_poly, 'a')
        print("order of prime in Phi is", order(local_primitive_root(self.Phi)))

        self.a = [element(i,self.Phi) for i in range(1, n + 1)]
        self.h = self.a
        alpha = local_primitive_root(self.Phi)
        print("order of alpha", order(alpha))
        self.C = local_BCH_Code(alpha, self.r + 1,  1, self.F)
        print("crack1")
        self.H = H_(self.C)
        print(self.H)
        self.H = prune(blow(H_(self.C),self.F))
        print(self.H)
        print(G_(self.C))
        self.G = left_kernel(transpose(self.H))"""
        # self.C = alternant_code(self.h, self.a, self.r, self.F)
        #self.H = prune(blow(H_(self.C),self.F))
        #if G_(self.C)!= None:
        #    self.G = G_(self.C)
        #else: 
        #    self.G = left_kernel(transpose(self.H))
        """print("crack2")
        self.F_map = create_dict_of_field(self.F, base ** f_extension)
        self.Phi_map = create_dict_of_field(self.Phi, base ** phi_extension)
        print("crack3")"""
        
    def create_systematic_G(self):
        one_ = element(1,self.F)
        def _f_x(x, alpha, skip_i=-1):
            res = one_
            for i in range(self.k):
                if skip_i != i:
                    res *= (x - alpha ** i)
            return res

        def _f_x_i(x, alpha, i):
            mone = _f_x(x, alpha)
            mechane = (x - alpha ** i)
            return mone / mechane
        
        # find primitive element in the RS field
        for i in range(3, self.F_size):
            elem = element(i, self.F)
            for j in range(1, self.F_size - 1):
                if elem ** j == one_:
                    break
                if j == self.F_size - 2:
                    alpha = elem
                    break
        # create the RS field using the primitive element
        a = geometric_series(alpha, self.n) 
        C = RS(a,self.k)
        # convert G to be systematic
        for i in range(self.k):
            for j in range(self.n):
                mone = _f_x(alpha **j, alpha, i)
                mechane = _f_x(alpha ** i, alpha, i)
                val = mone / mechane
                G_(C)[i,j] = val
        self.G = G_(C)
        self.H = H_(C)
        return C
            
            
    def encode(self, u: np.array):
        # u is an array of integers of the syndromes
        # need to convert each number to its element
        u = [element(u_i, self.F) for u_i in u]
        # after converting each one 
        u = vector(u)
        return self.encode_vector(u)

    def encode_vector(self, u_vec) -> np.array:
        G = self.G
        word = u_vec
        
        return u_vec * self.G
    
    
    def decode(self, y: np.array, e: np.array) -> np.array:
        
        erasure_locations = [i for i,x in enumerate(e) if x != 0]
        bad_H = submatrix(self.H, erasure_locations)
        no_error_locations = [i for i,x in enumerate(e) if x == 0]
        good_H = submatrix(self.H, no_error_locations)

        y_good = [y[i] for i in no_error_locations]
 
        y_good = vector(self.F, y_good)
        syn = good_H * y_good
        missing_bits = solve_linear_system(bad_H, -syn)
        c = y
        for i, loc in enumerate(erasure_locations):
            c[loc] = missing_bits[i]
        x_ = c[:self.k]

        return x_, c
    
    
    
    def get_rate(self) -> np.float32:
        """
        returns the calculation of the current codes rate
        """
        return (self.n - self.r) / self.n
    
    
if __name__ == "__main__":
    C = ReedSoloWrapper(n=16, r=8, base=17, f_extension=1)
    num_of_erasures = 7
    for i in range(0, 300):
        u = rd_error_vector(C.k, C.k // 2, C.F_base)
        c = C.encode_vector(u)
        """print("pip", C.G._value)
        exit()"""
        e = rd_error_vector(C.n, num_of_erasures, C.F)
        print(c)
        y = c + e
        x_, c_ = C.decode(y, e)
        assert all(x_[i] == u[i] for i in range(len(x_))), (u, x_)
        
    print("Test is successful")
    
    



