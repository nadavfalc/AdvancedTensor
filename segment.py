import numpy as np

class Segment:
    def __init__(self, information_bits, coded_segment=None):
        self.information_bits = np.array(information_bits)
        self.coded_segment = np.array(coded_segment)

