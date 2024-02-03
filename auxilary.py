from math import comb
import numpy as np
import pytest


import numpy as np

def convert_to_rank_sequence(x: np.array):
    """
    Convert an array to a rank sequence.

    Args:
        x (np.array): The input array.

    Returns:
        np.array: The rank sequence array.
    """
    x = np.insert(x, 0, 0)
    x = np.append(x, 1)
    r = np.zeros(1)
    for i in range(1, len(x)):
        r = np.append(r, (r[i - 1] if x[i] == x[i - 1] else r[i - 1] + 1))
    return r


syndromes = [
    lambda x: sum((i + 1) * x_i for i, x_i in enumerate(x)),
    lambda x: sum(comb(i + 1, 2) for i, x_i in enumerate(x)),
    lambda x: sum((i + 2) * r_i for i, r_i in enumerate(convert_to_rank_sequence(x))),
    lambda x: sum(comb(r_i, 2) for i, r_i in enumerate(1, convert_to_rank_sequence(x))),
]


def activate_syndrome(syndrome_index, x: np.array):
    """
    Activates a syndrome on the input array.

    Args:
        syndrome_index (int): The index of the syndrome to activate. Must be in the range [1, 4].
        x (np.array): The input array.

    Returns:
        np.array: The result of activating the specified syndrome on the input array.
    Raises:
        ValueError: If the syndrome index is not in the range [1, 4].
    """
    if syndrome_index not in range(1, 5):
        raise ValueError("Syndrome index must be in range [1, 4]")
    return syndromes[syndrome_index - 1](x)


def main():
    x = np.array([0,0,1,0,0,0,1,1,1,0,1,0])
    print(convert_to_rank_sequence(x))
    
    

if __name__ == "__main__":
    main()

