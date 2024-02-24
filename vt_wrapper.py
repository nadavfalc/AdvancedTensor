from code_with_syndrome import CodeWithSyndrome
from vt import *
import numpy as np


def remove_elements_by_indexes(A, B):
    # Sort the indexes in descending order to avoid index shifting when removing elements
    new_A = []    
    # Remove elements from list A based on the indexes in list B
    for i, a in enumerate(A):
        if i + 1 not in B:
            new_A.append(a)
    
    return new_A


class VTWrapper(CodeWithSyndrome):
    def __init__(self, q, n, systematic_locs=[]):
        k = find_k(n, q)
        super().__init__(q, n, k)
        self.systematic_locs = []
        if systematic_locs == []:
            i = 0
            while q ** i < n:
                self.systematic_locs.append(q ** i)
                i += 1
        else:
            self.systematic_locs = systematic_locs
        
    def encode(self, x: np.array, syndrome: int) -> np.array:
        """
        encode a codeword of length k to a VT_0(self.n) codeword with systematic positions received as inputs
        """
        vt = VTCode(self.n, self.q, syndrome, systematic_locs=self.systematic_locs) 
        return vt.encode(x)
        
    def decode(self, c: np.array, syndrome: int) -> np.array:
        """
        Receives an errorneous sequence
        Returns the decoded sequence
        """
        vt = VTCode(self.n, self.q, syndrome, systematic_locs=self.systematic_locs) 
        x, c = vt.decode(c)
        return x, c

    def get_rate(self) -> np.float32:
        """
        returns the calculation of the current codes rate
        """
        return self.k / self.n
    
    def calculate_syndrome(self, x: np.array) -> int:
        """
        Receives an input sequence
        Returns the syndrome of the input sequence
        """
        return int(sum([x[i] * (i + 1) for i in range(len(x))])) % (self.n + 1)
    
    
    