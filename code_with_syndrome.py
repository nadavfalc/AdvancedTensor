from code import *
import numpy as np

class CodeWithSyndrome(Code):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @abc.abstractmethod
    def encode(self, x: np.array, syndrome: np.array) -> np.array:
        """
        Receives an input sequence 
        Returns the encoded input sequence
        """
        pass
    
    @abc.abstractmethod
    def decode(self, x: np.array, syndrome: np.array) -> np.array:
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
    
    @abc.abstractmethod
    def calculate_syndrome(self, x: np.array) -> np.array:
        """
        Receives an input sequence
        Returns the syndrome of the input sequence
        """
        pass