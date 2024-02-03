import math
from code import Code
import numpy as np

def levenshtein_distance(s, t):
    """Calculates the Levenshtein distance between two strings.

    Args:
        s (str): The first string.
        t (str): The second string.

    Returns:
        int: The Levenshtein distance between s and t.
    """
    m = len(s)
    n = len(t)
    d = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        d[i][0] = i
    for j in range(n + 1):
        d[0][j] = j
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            if s[i - 1] == t[j - 1]:
                d[i][j] = d[i - 1][j - 1]
            else:
                d[i][j] = min(d[i - 1][j], d[i][j - 1], d[i - 1][j - 1]) + 1
    return d[m][n]


def hamming_distance(a: np.array, b: np.array) -> int:
    """
    Calculates the Hamming distance between two strings.

    Args:
        s1 (np.array): The first np.array.
        s2 (np.array): The second np.array.

    Returns:
        int: The Hamming distance between s1 and s2.
    """
    if len(a) != len(b):
        print("len difference=", len(a)-len(b))
        raise ValueError("Strings must be of equal length")
    return np.count_nonzero(np.not_equal(a, b))



def success_rate(X: dict, Y: dict) -> float:
    """
    Returns the percentage of distinct strings with the same keys between the two dictionaries.

    Args:
        X (dict): The first dictionary.
        Y (dict): The second dictionary.

    Returns:
        float: The percentage of distinct strings with the same keys between X and Y.
    """
    num_matches = 0
    num_keys = len(set(X.keys()) | set(Y.keys()))
    for key in set(X.keys()) & set(Y.keys()):
        if X[key] == Y[key]:
            num_matches += 1
    return (num_matches / num_keys) * 100 if num_keys > 0 else 0.0

    
def calculate_NER(X: dict, Y: dict) -> float:
    """returns the nucleotide error rate between the input and the output sequences

    Args:
        X (dict): The first dictionary.
        Y (dict): The second dictionary.

    Returns:
        float: NER value
    """
    m = len(X[1])
    M = len(X)
    print(len(X), len(Y), len(X[1]), len(Y[1]))
    total_distance = 0
    for id in X.keys():
        total_distance += hamming_distance(X[id], Y[id])
    
    return total_distance / (M * m)


def calculate_average_edit_distance(X: dict, Y: dict) -> float:
    """returns the nucleotide error rate between the input and the output sequences

    Args:
        X (dict): The first dictionary.
        Y (dict): The second dictionary.

    Returns:
        float: NER value
    """
    m = len(X[1])
    M = len(set(X.keys()))
    total_distance = 0
    for id in X.keys():
        total_distance += levenshtein_distance(X[id], Y[id])
    
    return total_distance / (M * m)


def encode_batch(c: Code, X: dict) -> dict:
    """
    Encodes each string value in X using c.encode() method.

    Args:
        c (Code): An instance of the Code class.
        X (dict{int: str}): A dictionary that maps integer keys to string values.

    Returns:
        dict{int: str}: A dictionary that maps the same integer keys to the encoded string values.
    """
    Y = {}
    for k, v in X.items():
        Y[k] = c.inner_encoder(v)
    return Y


def decode_batch(c: Code, X: dict, error_locations=None) -> dict:
    """
    Encodes each string value in X using c.encode() method.

    Args:
        c (Code): An instance of the Code class.
        X (dict{int: str}): A dictionary that maps integer keys to string values.

    Returns:
        dict{int: str}: A dictionary that maps the same integer keys to the encoded string values.
    """
    Y = {}
    for k, v in X.items():
        Y[k] = c.inner_decoder(v, error_locations)
    return Y


if __name__ == '__main__':
    print(hamming_distance("ACGGG", "ACGCA"))
    print(levenshtein_distance("", "ACGCA"))
    print(success_rate({1:"kaki", 33:"shimon"}, {1:"kaki", 33:"shion"}))
    print(calculate_NER({1:"12345", 33:"shimon"}, {1:"12345", 33:"shifon"}))