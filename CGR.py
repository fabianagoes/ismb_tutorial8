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
            elif nucle in ["T", "U"]:
                Sx.append(-1)
                Sy.append(1)
            else:  # G or unknown
                Sx.append(1)
                Sy.append(-1)

        CGR_x = []
        CGR_y = []

        for i in range(len(Sx)):
            if i == 0:
                CGR_x.append(0.5 * Sx[i])
                CGR_y.append(0.5 * Sy[i])
            else:
                CGR_x.append(0.5 * Sx[i] + 0.5 * CGR_x[i - 1])
                CGR_y.append(0.5 * Sy[i] + 0.5 * CGR_y[i - 1])

        concat = CGR_x + CGR_y
        padded = np.pad(concat, (0, 2 * (max_length - len(seq))), 'constant')
        features.append(padded)

    return features
