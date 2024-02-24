from code import Code
import random
import numpy as np
from segment import Segment
from typing import Callable
from code_with_syndrome import CodeWithSyndrome
from PyM import *
from vt_wrapper import *
from reedsolo_wrapper import *


def convert_field_element_to_index(el, F):
    i = 0
    while 1 == 1:
        if element(i, F) == el:
            return i
        
        i += 1


def decimal_to_binary_bits(decimal_number, final_len):
    # Convert decimal number to binary string
    binary_string = bin(decimal_number)[2:]
    
    # Pad the binary string with zeros to ensure it has 8 bits
    padded_binary_string = binary_string.zfill(final_len)
    binary_array = np.array([int(bit) for bit in padded_binary_string])
    return binary_array


def binary_array_to_decimal(binary_array):
    # Convert the binary array to a binary string
    decimal_number = 0
    for i, x in enumerate(binary_array):
        decimal_number += int(x * 2 ** (len(binary_array) - i - 1))

    return decimal_number



class TensorProductCodes(Code):
    def __init__(self, q, inner_code: CodeWithSyndrome, outer_code: Code):
        self.inner_code = inner_code
        self.outer_code = outer_code
        self.number_of_constraintless_segments = self.outer_code.k
        self.row_length = self.inner_code.n
        self.column_length = self.outer_code.n
        k = self.inner_code.n * self.outer_code.k 
        k += self.inner_code.k * (self.outer_code.n - self.outer_code.k)
        super().__init__(q=q, n=inner_code.n * outer_code.n, k=k)
        Z2 = Zn(2)
        # creating 2^8
        F_poly = get_irreducible_polynomial(Z2, 8)
        [self.F_info,b]= extension (Z2, F_poly, 'a')
        self.segments = []
        self.syndrome_array = np.empty(self.column_length)
        self.create_info_bits_rs()
        
        
    def encode(self, u: np.array) -> np.array:
        """
        Encode the received codeword into segments with the inner code.
        args:
            u (np.array): information bits
        returns:
            x (np.array): encoded codeword"""
        
        end_location_of_constraintless_segments = self.number_of_constraintless_segments * self.inner_code.n
        # u[i:i+self.inner_code.n], self.inner_code.systematic_locs
        # insert freely the first constraintless segments
        self.segments = [Segment(information_bits=u[i:i+self.inner_code.n],
                                 coded_segment=u[i:i+self.inner_code.n])
                                 for i in range(0, end_location_of_constraintless_segments, self.inner_code.n)
                                 ]
        # calculate the constraintless segments syndromes
        syndrome_array_known_symbols = [0] * self.outer_code.k
        for i in range(self.number_of_constraintless_segments):
            segment_i_syn = self.inner_code.calculate_syndrome(self.segments[i].coded_segment)

            syndrome_array_known_symbols[i] = segment_i_syn #element(segment_i_syn, self.outer_code.F)
        # look at each syndrome as a symbol in the extended code and decode them
        decoded_symbols = self.outer_code.encode(syndrome_array_known_symbols)
        int_syndromes = [convert_field_element_to_index(decoded_symbols[i], self.outer_code.F) for i in range(self.outer_code.n)]
        self.syndrome_array = decoded_symbols # add oposite side conversion
        start_location_of_constraint_segments = end_location_of_constraintless_segments
        
        self.segments += [Segment(information_bits = u[i:i+self.inner_code.k], 
                                  coded_segment = self.inner_code.encode(u[i:i+self.inner_code.k], 
                                                                         int_syndromes[self.outer_code.k + (i - start_location_of_constraint_segments) // self.inner_code.k ]
                                                                         )
                                  )
                          for i in range(start_location_of_constraint_segments, len(u), self.inner_code.k)
                          ]
    
    def decode_rs(self, y: np.array, e: np.array) -> np.array:
        erasure_locations = [i for i,x in enumerate(e) if x != 0]
        bad_H = submatrix(self.H_info, erasure_locations)
        no_error_locations = [i for i,x in enumerate(e) if x == 0]
        good_H = submatrix(self.H_info, no_error_locations)

        y_good = [y[i] for i in no_error_locations]

        y_good = vector(self.F_info, y_good)
        syn = good_H * y_good
        missing_bits = solve_linear_system(bad_H, -syn)
        c = y
        for i, loc in enumerate(erasure_locations):
            print(loc, i)
            c[loc] = missing_bits[i]
        x_ = c[:23]
        return x_, c
    
    
    def create_info_bits_rs(self):
        F_info = self.F_info
        n = self.k // 8
        k = n - 4
        a = [element(j,F_info) for j in range(2,256) if is_prime(j)][:n]
        h = [element(j,F_info) for j in range(1,1+len(a))]

        self.H_info = vandermonde(a, 4)
        self.G_info = left_kernel(transpose(self.H_info))
    
    
    def _create_tensor_matrix(self, y: np.array) -> np.array:
        """
        This function determines the segments out of the received codeword.
        args:
            y (np.array): received codeword
        returns:
            tensor_matrix (np.array): matrix of segments,each line is a different segment
        """
        pass
        
    def decode(self, y, segments_estimation=[]) -> np.array:
        """
        This function receives a noisy codeword and returns the decoded codeword
        """
        assert segments_estimation != [], "Needs to be implemented!"
        
        final_result = np.empty(0)
        phantom_syn_vector = [0] * self.outer_code.n
        
        erasure_vector = [0] * self.outer_code.n
        
        for i, s in enumerate(segments_estimation):
            if len(s) == self.inner_code.n:
                # no deletions in this segment
                seg_syn = self.inner_code.calculate_syndrome(s)
                phantom_syn_vector[i] = element(seg_syn, self.outer_code.F)
            else:
                erasure_vector[i] = 1
        if sum(erasure_vector) > self.outer_code.r:
            print("too many errors!!")
            return np.zeros(self.k)
        _, phantom_syn_vector  = self.outer_code.decode(phantom_syn_vector, erasure_vector)
        
        # from this point we have the syndromes of all the segment
        phantom_syn_vector_int = [convert_field_element_to_index(el, self.outer_code.F) for el in phantom_syn_vector]
        
        error_locations = np.empty(0)
        for i, s in enumerate(segments_estimation):
            cur_seg_syndrome = phantom_syn_vector_int[i]
            decoded_segment_bits = []
            decoded_vt_bits = []
            if len(s) < self.inner_code.n - 1:
                # the segment is corrupted. we want to remember it and save the erasued locations
                # based on wether that segment was encoded with vt or not
                print(f"segment {i} has >1 errors")
                cur_decoded_result_length = len(final_result)
                if i < self.number_of_constraintless_segments:
                    # all of this segments bits were information bits
                    num_of_erased_bits = self.inner_code.n
                else:
                    # this segment information bits were encoded to vt codeword
                    num_of_erased_bits = self.inner_code.k
                new_error_locations = [j for j in range(
                    cur_decoded_result_length, cur_decoded_result_length+num_of_erased_bits
                    )]
                # add the new error locations
                error_locations = np.concatenate([error_locations, new_error_locations])
                # add bits for structure of final result
                final_result = np.concatenate([final_result, np.array([-1] * num_of_erased_bits)])
            else:
                # there were 0 or 1 deletions and we know the vt syndrome, so 
                # these are the segments we can decode
                # there was one deletion, we need to use the vt decoder.
                if len(s) == self.inner_code.n:
                    decoded_vt_bits, _ = self.inner_code.decode(s, cur_seg_syndrome)
                    decoded_segment_bits = s
                else:
                    decoded_vt_bits, decoded_segment_bits = self.inner_code.decode(s, cur_seg_syndrome)
                if i < self.number_of_constraintless_segments:
                    information_bits_of_seg = decoded_segment_bits
                else:
                    information_bits_of_seg = decoded_vt_bits
                final_result = np.concatenate([final_result, information_bits_of_seg])
            
        # at this point we have the information bits one after the other, with -1 in locations of erasures.
        # no we prepare to use the inner RS code. we need to create a vector.
        final_result_in_symbols = []
        erasure_locations = []
        assert len(final_result) == self.k, (len(final_result), self.k)
        for i in range(0, len(final_result), 8):
            seg_bits = final_result[i: i+8]
            if -1 in seg_bits:
                erasure_symbol = 1
                new_symbol = element(0, self.F_info)
            else:
                bits_to_integer = binary_array_to_decimal(seg_bits)
                new_symbol = element(bits_to_integer, self.F_info)
                erasure_symbol = 0
            final_result_in_symbols += [new_symbol]
            erasure_locations += [erasure_symbol]
        if sum(erasure_locations) > 4:
            # cant be recovered
            print("erasures2:", sum(erasure_locations))
            return np.zeros(self.k)
        
        x_, c = self.decode_rs(vector(final_result_in_symbols), erasure_locations)
        x_binary = np.empty(0)
        for el in x_:
            symbol_in_bits = convert_field_element_to_index(el, self.F_info)
            seymbol_bits = decimal_to_binary_array(symbol_in_bits, 8)
            x_binary = np.concatenate([x_binary, seymbol_bits])
        return x_binary
    
    def get_rate(self):
        return (self.k - 32)  / self.n 
    
    
def make_errors_in_tensor(segments, p_d):
    
    for i in range(len(segments)):
        cur_segment = segments[i]
        num_dels = 0
        for j in range(len(cur_segment)):
            random_number = random.random()
            if random_number < p_d:
                print(f"error in {i, j}")
                segments[i] = np.delete(segments[i], j - num_dels)
                num_dels += 1
                
    return segments
    

def generate_information_bits(tc: TensorProductCodes):
    """
    Split the information bits to 8'th.
    there are 216 bits of information.
    lets split them to 8's. There are 27 segments like this ->
    GRS of 2^8 of length 
    """
    
    H = tc.H_info
    G = tc.G_info
    F_info = tc.F_info
    
    info_bits = np.random.randint(2, size=tc.k - 4 * 8, dtype=np.uint8)
    
    word = []
    for i in range(0, len(info_bits), 8):
        in_bits = info_bits[i:i+8]
        packed_bytes = np.packbits(in_bits)
        decimal_number = int.from_bytes(packed_bytes, byteorder='big')
        word.append(element(decimal_number,F_info))

    word = vector(word)
    c = word * G
    
    return c
    
    



def test_code(N, p_d, tc):
    
    
    F_info = tc.F_info
    total_success = 0
    for i in range(N):
        # create the information bits- currently are a rs code.
        info_bits = generate_information_bits(tc)
        binary_bits = np.empty(0)
        for el in info_bits:
            element_in_int = convert_field_element_to_index(el, F_info)
            element_in_binary_representation = decimal_to_binary_bits(element_in_int, 8)
            binary_bits = np.concatenate([binary_bits, element_in_binary_representation])
        assert len(binary_bits) == tc.k, (len(binary_bits), tc.k)
        print("Information bits to encode")
        tc.encode(binary_bits)

        info_segments = [s.information_bits for s in tc.segments]
        info_segments = [int(item) for sublist in info_segments for item in sublist]
        encoded_segments = [s.coded_segment for s in tc.segments]
        for s in encoded_segments:
            assert len(s) == tc.inner_code.n
        segments_with_errors = make_errors_in_tensor(encoded_segments, p_d)
        
        decoded_segments = tc.decode(segments_with_errors, segments_with_errors)
        
        if len(decoded_segments) != 184:
            print("len(decoded_segments)", len(decoded_segments))
            print("tc.outer_code.n=", tc.outer_code.n)
            print(f"decoding {i} unssucsesful")
            continue
        
        
        
        is_equal = True
        
        for j in range(184):
            if int(info_segments[j]) != int(decoded_segments[j]):
                print("diff in index", j)
                is_equal = False
        
        if is_equal:            
            print(f"test successful {i}")
            total_success += 1
        else:
            print(f"Mistake in decoding {i}")
            for i in range(0, tc.number_of_constraintless_segments * 16, 16):
                print("<<<<<<Segment ", i // 16)
                print(info_segments[i:i+16])
                print(decoded_segments[i:i+16])
                print(">>>>>>")
            for i in range(tc.number_of_constraintless_segments * 16, 184, 11):
                print("<<<<<<Segment ", tc.number_of_constraintless_segments + i // 16)
                print(info_segments[i:i+11])
                print(decoded_segments[i:i+11])
                print(">>>>>>")

        if i % 10 == 0:
            print("rate:", tc.get_rate())
            print(f"total length= {tc.n}")
            FER = (i + 1 - total_success) / (i + 1)
            print(f"FER= {FER}")
    return FER

if __name__ == "__main__":
    q = 2
    n_inner = 16
    inner_code = VTWrapper(q, n_inner) 
    k_inner = inner_code.k

    assert k_inner != 0, (k_inner, inner_code.n)
    n_outer = 16
    r_outer = 8
    k_outer = n_outer - r_outer
    outer_extension = 1
    outer_code = ReedSoloWrapper(n_outer, r_outer, n_outer + 1, outer_extension)
    
    tc = TensorProductCodes(q, inner_code, outer_code)
    FERS = {}
    FERS[0.01] = 0.16206666666666666
    for p_d in [0.02, 0.03, 0.015, 0.01]:
        FERS[p_d] = test_code(10000, p_d, tc)
        print(FERS[p_d])
    print(FERS)


