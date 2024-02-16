from AdvancedTensor.code import Code
import numpy as np
from segment import Segment
from typing import Callable
from code_with_syndrome import CodeWithSyndrome
from PyM import *


class TensorProductCodes(Code):
    def __init__(self, q, inner_code: CodeWithSyndrome, outer_code: Code, number_of_constraintless_segments: int, k: int,
                 map_segment_to_symbol: Callable):
        super().__init__(q=q, n=inner_code.n * outer_code.n, k=k)
        self.inner_code = inner_code
        self.outer_code = outer_code
        self.row_length = self.inner_code.n
        self.column_length = self.outer_code.n
        self.number_of_constraintless_segments = number_of_constraintless_segments
        self.segments = []
        self.syndrome_array = np.empty((self.column_length, np.log2(self.row_length)))
        
    def encode(self, u: np.array) -> np.array:
        """
        Encode the received codeword into segments with the inner code.
        args:
            u (np.array): information bits
        returns:
            x (np.array): encoded codeword"""
        end_location_of_constraintless_segments = self.number_of_constraintless_segments * self.row_length
        # insert freely the first constraintless segments
        self.segments = [Segment(information_bits=u[i:i+self.row_length],
                                 coded_segment=u[i:i+self.row_length])
                                 for i in range(0, end_location_of_constraintless_segments, self.row_length)
                                 ]
        # calculate the constraintless segments syndromes
        for i in range(end_location_of_constraintless_segments):
            self.syndrome_array[i] = self.inner_code.calculate_syndrome(self.segments[i])
        # look at each syndrome as a symbol in the extended code and decode them
        symbol_array_from_syndroms = self.syndrome_array # maybe add some conversion
        decoded_symbols = self.outer_code.encode(symbol_array_from_syndroms)
        self.syndrome_array = decoded_symbols # add oposite side conversion
        
        start_location_of_constraint_segments = end_location_of_constraintless_segments
        self.segments += [Segment(information_bits = u[i:i+self.inner_code.k], 
                                  coded_segment = self.inner_code.encode(u[i:i+self.inner_code.k]))
                          for i in range(start_location_of_constraint_segments, len(u), self.inner_code.k)
                          ]
        assert sum([len(segment.information_bits) for segment in self.segments]) == len(u)
        assert sum([len(segment.coded_segment) for segment in self.segments]) == self.n
    
    
    def _create_tensor_matrix(self, y: np.array) -> np.array:
        """
        This function determines the segments out of the received codeword.
        args:
            y (np.array): received codeword
        returns:
            tensor_matrix (np.array): matrix of segments,each line is a different segment
        """
        pass
        
    def decode(self, y: np.array) -> np.array:
        """
        This function receives a noisy codeword and returns the decoded codeword
        """
    
        tensor_matrix = self._create_tensor_matrix(y)
        
