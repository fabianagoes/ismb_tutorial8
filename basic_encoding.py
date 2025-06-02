from itertools import product
from collections import OrderedDict

def generate_kmer_features(sequences, kmax, alphabet="ACGT"):
    """
    Generates k-mer frequency features from a list of sequences.

    Args:
        sequences (list of str): List of nucleotide sequences.
        kmax (int): Maximum size of k-mers to generate.
        alphabet (str): Set of characters allowed in the k-mers (default: "ACGT").

    Returns:
        list of dict: Each dict contains k-mer frequencies for one sequence.
    """
    def get_all_kmers(k, alphabet):
        return [''.join(p) for p in product(alphabet, repeat=k)]

    def count_kmers(seq, k, kmer_template):
        counts = kmer_template.copy()
        total_windows = len(seq) - k + 1
        if total_windows <= 0:
            return {kmer: 0.0 for kmer in kmer_template}  # Handle short sequences
        for i in range(total_windows):
            kmer = seq[i:i+k]
            if kmer in counts:
                counts[kmer] += 1
        return {kmer: count / total_windows for kmer, count in counts.items()}

    all_features = []
    for seq in sequences:
        seq = seq.upper()
        seq_features = OrderedDict()
        for k in range(1, kmax + 1):
            kmer_template = {kmer: 0 for kmer in get_all_kmers(k, alphabet)}
            freqs = count_kmers(seq, k, kmer_template)
            seq_features.update(freqs)
        all_features.append(seq_features)
    
    return all_features

import numpy as np

def one_hot_encode_sequences(sequences, alphabet="ACGT"):
    """
    Converts a list of sequences into one-hot encoded numpy arrays.

    Args:
        sequences (list of str): List of sequences (all same length recommended).
        alphabet (str): Alphabet used for encoding (default: "ACGT").

    Returns:
        np.ndarray: 3D array of shape (num_sequences, seq_length, alphabet_size)
    """
    char_to_index = {char: idx for idx, char in enumerate(alphabet)}
    seq_len = max(len(seq) for seq in sequences)
    num_seqs = len(sequences)
    one_hot = np.zeros((num_seqs, seq_len, len(alphabet)), dtype=int)

    for i, seq in enumerate(sequences):
        for j, char in enumerate(seq.upper()):
            if char in char_to_index:
                one_hot[i, j, char_to_index[char]] = 1
            # else it remains zero (unknown character or padding)

    return one_hot

