import numpy as np

def classifical_chaos_features(sequences, max_length):
    """
    Computes Classifical Chaos Game Representation (CGR) features for each sequence.

    Args:
        sequences (list of str): List of nucleotide sequences.
        max_length (int): Padding target length for sequences.

    Returns:
        list of np.ndarray: Each array contains CGR features (padded) for a sequence.
    """
    features = []

    for seq in sequences:
        seq = seq.upper()
        Sx, Sy = [], []

        for nucle in seq:
            if nucle == "A":
                Sx.append(1)
                Sy.append(1)
            elif nucle == "C":
                Sx.append(-1)
                Sy.append(-1)
            el
