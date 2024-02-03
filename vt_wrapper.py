from code_with_syndrome import CodeWithSyndrome
from vt import VTCode, decimal_to_binary_array
import numpy as np

class VTWrapper(CodeWithSyndrome):
    def __init__(self, q, n, k, systematic_locs=[]):
        if systematic_locs == []:
            self.systemaic_locs = [2 ** i for i in range(int(np.log2(n)))]
        super().__init__(q=q, n=n, k=k, self.systemaic_locs)

    def encode(self, x: np.array, syndrome: np.array) -> np.array:
        """
        encode a codeword of length k to a VT_0(self.n) codeword with systematic positions received as inputs
        """
        vt = VTCode(n=self.n, q=self.q, a=syndrome, systematic_locs=self.systematic_locs)
        # assert len(x) == self.k, ("kaki", len(x), self.k)
        c_ = np.zeros(self.n)
        # mod n and not mod n+1 because in cs we start from 0
        c_i = 0
        
        for i in range(self.n):
            if i not in self.vt.systematic_locs:
                c_[i] = x[c_i]
                c_i += 1        
                
        # i+1 and not mod i because in cs we start from 0
        syn_c_ = sum((i+1) * c_i for i, c_i in enumerate(c_)) % (self.n + 1)
        for i in range(2 ** len(self.vt.systematic_locs)):
            bin_a = decimal_to_binary_array(i, pad=len(self.vt.systematic_locs))
            syn = (np.sum(np.multiply(self.vt.systematic_locs+np.ones(len(self.vt.systematic_locs)), bin_a))) % (self.n + 1)
            if (syn + syn_c_) % (self.n + 1) == self.syndrome:
                c_[self.vt.systematic_locs] = bin_a
                return c_
        assert True, "no codeword found"
        
    def decode(self, x: np.array, syndrome: np.array) -> np.array:
        """
        Receives an errorneous sequence
        Returns the decoded sequence
        """
        vt = VTCode(n=self.n, q=self.q, a=syndrome, systematic_locs=self.systematic_locs)
        return vt.decode(x)

    def get_rate(self) -> np.float32:
        """
        returns the calculation of the current codes rate
        """
        return self.k / self
    
    def calculate_syndrome(self, x: np.array) -> np.array:
        """
        Receives an input sequence
        Returns the syndrome of the input sequence
        """
        return [x[i] * (i + 1) for i in range(len(x))]