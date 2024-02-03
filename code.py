import abc
import string
import numpy as np

class Code(abc.ABC):
    
    def __init__(self, q, n, k):
        """
        Args:
            q (int): q-ary code
            n (int): length of the code
            k (int): number of information bits
        """
        self.q = q
        self.n = n
        self.k = k
        
    @abc.abstractmethod
    def encode(self, x: np.array) -> np.array:
        """
        Receives an input sequence 
        Returns the encoded input sequence
        """
        pass
    
    @abc.abstractmethod
    def decode(self, x: np.array, erasure_locations=np.empty(0)) -> np.array:
        """
        Receives an errorneous sequence
        Returns the decoded sequence
        """
        pass

    @abc.abstractmethod
    def get_rate(self) -> np.float32:
        """
        returns the calculation of the current codes rate
        """
        pass
    