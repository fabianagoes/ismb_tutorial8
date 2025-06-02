import numpy as np

def zcurve_features(sequences, max_length):
    """
    Compute Z-curve features for a list of sequences.

    Args:
        sequences (list of str): DNA/RNA sequences.
        max_length (int): Desired sequence length for padding.

    Returns:
        list of np.ndarray: Z-curve vectors (padded) per sequence.
    """
    features = []
    for seq in sequences:
        seq = seq.upper()
        R = Y = M = K = W = S = 0
        x, y, z = [], [], []

        for nucle in seq:
            if nucle in "AG":
                R += 1
            else:
                Y += 1
            x.append(R - Y)

            if nucle in "AC":
                M += 1
            else:
                K += 1
            y.append(M - K)

            if nucle in "ATU":
                W += 1
            else:
                S += 1
            z.append(W - S)

        concat = x + y + z
        padded = np.pad(concat, (0, 3 * (max_length - len(seq))), 'constant')
        features.append(padded)
    
    return features

def integer_features(sequences, max_length):
    """
    Convert sequences to integer-mapped vectors.

    A:2, C:1, G:3, T/U:0

    Args:
        sequences (list of str): DNA/RNA sequences.
        max_length (int): Desired length for padding.

    Returns:
        list of np.ndarray: Integer-encoded padded vectors.
    """
    mapping_dict = {"A": 2, "C": 1, "G": 3, "T": 0, "U": 0}
    features = []

    for seq in sequences:
        seq = seq.upper()
        mapping = [mapping_dict.get(nuc, 0) for nuc in seq]
        padded = np.pad(mapping, (0, max_length - len(seq)), 'constant')
        features.append(padded)
    
    return features
