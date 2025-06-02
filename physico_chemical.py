def gc_content_features(sequences):
    """
    Computes GC content and related features for a list of sequences.

    Args:
        sequences (list of str): List of nucleotide sequences.

    Returns:
        list of dict: Each dict contains GC features for one sequence.
    """
    features = []
    for seq in sequences:
        seq = seq.upper()
        g = seq.count("G")
        c = seq.count("C")
        gc = g + c
        length = len(seq)
        gc_ratio = gc / length if length > 0 else 0
        features.append({
            "GC_Content": gc_ratio,
            "G_Count": g,
            "C_Count": c,
            "Length": length
        })
    return features

def global_descriptor_features(sequences, alphabet="ACGT"):
    """
    Computes global descriptor features: composition, transitions, and distribution.

    Args:
        sequences (list of str): List of RNA sequences.
        alphabet (str): Set of allowed characters (default: ACGT).

    Returns:
        list of dict: One feature dict per sequence.
    """
    features = []

    for seq in sequences:
        seq = seq.upper()
        seq_len = len(seq)
        feature_dict = {}

        # Composition: frequency of each nucleotide
        for nt in alphabet:
            feature_dict[f"Comp_{nt}"] = seq.count(nt) / seq_len if seq_len > 0 else 0

        # Transitions: frequency of NT1 -> NT2 (e.g., A->C, G->U, etc.)
        for i in range(len(seq) - 1):
            pair = seq[i:i+2]
            if all(nt in alphabet for nt in pair):
                feature_dict[f"Trans_{pair}"] = feature_dict.get(f"Trans_{pair}", 0) + 1
        # Normalize transitions
        for key in list(feature_dict):
            if key.startswith("Trans_"):
                feature_dict[key] /= (seq_len - 1) if seq_len > 1 else 1

        # Distribution: first, 25%, 50%, 75%, and last position index of each nt
        for nt in alphabet:
            positions = [i for i, base in enumerate(seq) if base == nt]
            if positions:
                percentiles = [positions[0],                             # First
                               positions[int(len(positions)*0.25)],     # 25%
                               positions[int(len(positions)*0.5)],      # 50%
                               positions[int(len(positions)*0.75)],     # 75%
                               positions[-1]]                           # Last
                for i, p in enumerate(percentiles):
                    feature_dict[f"Dist_{nt}_{i}"] = p / seq_len
            else:
                for i in range(5):
                    feature_dict[f"Dist_{nt}_{i}"] = 0.0

        features.append(feature_dict)
    
    return features
