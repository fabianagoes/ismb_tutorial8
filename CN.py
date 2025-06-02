import numpy as np
import pandas as pd
from igraph import Graph, ADJ_UNDIRECTED


def patterns(seq, win):
    """Generate k-mers with window size `win` from the input sequence."""
    for i in range(len(seq)):
        j = len(seq) if i + win > len(seq) else i + win
        yield seq[i:j]
        if j == len(seq):
            break


def extract_graph_metrics(graph):
    """Extract topological features from the given graph."""
    metrics = []
    metrics.append(np.mean(graph.betweenness(directed=False)))
    metrics.append(np.mean(graph.degree()))
    metrics.append(graph.assortativity_degree(directed=False))
    metrics.append(max(graph.degree()))
    metrics.append(min(graph.degree()))
    metrics.append(np.std(graph.degree()))
    metrics.append(graph.average_path_length(directed=False, unconn=False))
    metrics.append(graph.transitivity_avglocal_undirected())
    metrics.append(graph.transitivity_undirected())
    metrics.append(graph.ecount())
    metrics.append(graph.motifs_randesu_no(size=3))
    metrics.append(graph.motifs_randesu_no(size=4))
    return metrics


def complex_network_features(sequences, ksize=3, threshold=20):
    """
    Computes complex network features from a list of sequences using iGraph.

    Args:
        sequences (list of str): Input nucleotide sequences.
        ksize (int): K-mer size (default 3).
        threshold (int): Threshold limit for graph simplification.

    Returns:
        list of dict: Each dict contains complex network features for a sequence.
    """
    all_features = []

    for seq in sequences:
        seq = seq.upper()
        kmer_list = list(patterns(seq, ksize))
        seen = set()
        g = Graph()

        # Build original graph
        for i in range(len(kmer_list) - 1):
            a, b = kmer_list[i], kmer_list[i + 1]
            for node in (a, b):
                if node not in seen:
                    g.add_vertex(name=node)
                    seen.add(node)
            g.add_edge(a, b)

        # Extract features for thresholds
        seq_features = {}
        for t in range(1, threshold):
            adj_matrix = pd.DataFrame(g.get_adjacency().data)
            adj_matrix_thresh = np.where(adj_matrix < t, 0, adj_matrix)
            g_thresh = Graph.Adjacency(adj_matrix_thresh.astype(int).tolist(), mode=ADJ_UNDIRECTED)

            if g_thresh.ecount() < 1:
                for i in range(12):
                    seq_features[f"{i}_T{t}"] = 0
                continue

            metrics = extract_graph_metrics(g_thresh)
            for i, val in enumerate(metrics):
                seq_features[f"{i}_T{t}"] = val

        all_features.append(seq_features)

    return all_features
