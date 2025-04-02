#!/usr/bin/env python3

"""
Purpose: Clusters entropy profiles, assigns probabilities to SV calls.
Author: Gutierrez-Escobar-AJ
"""

import argparse
import os
from typing import List, Tuple
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cophenet
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE
from sklearn.neighbors import KernelDensity
from dtaidistance import dtw

# ========== HELPER FUNCTIONS ==========
def load_entropy_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required_cols = {'Sequence_ID', 'Window_Mid', 'Entropy'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in input file: {required_cols - set(df.columns)}")
    return df

def compute_distance_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, List[str]]:
    pivoted = df.pivot(index='Sequence_ID', columns='Window_Mid', values='Entropy').fillna(0)
    labels = pivoted.index.tolist()
    series = [pivoted.loc[label].values for label in labels]
    dist_matrix = np.zeros((len(series), len(series)))
    for i in range(len(series)):
        for j in range(i+1, len(series)):
            d = dtw.distance_fast(series[i], series[j])
            dist_matrix[i, j] = dist_matrix[j, i] = d
    return dist_matrix, labels

def optimal_num_clusters(dist_matrix: np.ndarray, max_clusters: int = 10) -> int:
    scores = {}
    for k in range(2, min(max_clusters + 1, len(dist_matrix))):
        linkage_matrix = linkage(squareform(dist_matrix), method='ward')
        labels = fcluster(linkage_matrix, k, criterion='maxclust')
        score = silhouette_score(dist_matrix, labels, metric='precomputed')
        scores[k] = score
    best_k = max(scores, key=scores.get)
    print(f"[✔] Optimal number of clusters: {best_k} (Silhouette = {scores[best_k]:.2f})")
    return best_k

def plot_entropy_heatmap(df: pd.DataFrame, output_path: str):
    pivot = df.pivot(index='Sequence_ID', columns='Window_Mid', values='Entropy').fillna(0)
    sns.clustermap(pivot, cmap='viridis', figsize=(12, 12), method='ward')
    plt.title("Entropy Heatmap with Hierarchical Clustering")
    plt.savefig(output_path)
    plt.close()
    print(f"[✔] Heatmap saved: {output_path}")

def save_newick(linkage_matrix, labels, output_path):
    from io import StringIO
    from scipy.cluster.hierarchy import to_tree

    def build_newick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return f"{leaf_names[node.id]}:{parentdist - node.dist}{newick}"
        else:
            left = build_newick(node.left, newick, node.dist, leaf_names)
            right = build_newick(node.right, newick, node.dist, leaf_names)
            return f"({left},{right}):{parentdist - node.dist}{newick}"

    tree = to_tree(linkage_matrix, rd=False)
    newick_str = build_newick(tree, ";", tree.dist, labels)
    with open(output_path, "w") as f:
        f.write(newick_str + "\n")
    print(f"[✔] Newick tree saved: {output_path}")

def estimate_sv_probabilities(df: pd.DataFrame) -> pd.Series:
    kde = KernelDensity(kernel='gaussian', bandwidth=0.1).fit(df[['Entropy']])
    log_probs = kde.score_samples(df[['Entropy']])
    probs = np.exp(log_probs)
    norm_probs = (probs - probs.min()) / (probs.max() - probs.min())
    return norm_probs

def cluster_and_visualize(df: pd.DataFrame, out_prefix: str):
    print("[+] Computing DTW distance matrix...")
    dist_matrix, labels = compute_distance_matrix(df)

    print("[+] Running Agglomerative Clustering...")
    linkage_matrix = linkage(squareform(dist_matrix), method='ward')
    coph_corr, _ = cophenet(linkage_matrix, squareform(dist_matrix))
    print(f"[✔] Cophenetic correlation: {coph_corr:.2f}")

    print("[+] Generating optimal cluster count...")
    n_clusters = optimal_num_clusters(dist_matrix)
    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

    cluster_map = dict(zip(labels, cluster_labels))
    df['Cluster'] = df['Sequence_ID'].map(cluster_map)

    print("[+] Estimating SV probabilities using KDE...")
    df['SV_Prob'] = estimate_sv_probabilities(df)

    print("[+] Generating entropy heatmap...")
    plot_entropy_heatmap(df, f"{out_prefix}_heatmap.png")

    print("[+] Generating linkage for Newick tree...")
    save_newick(linkage_matrix, labels, f"{out_prefix}.nwk")

    print("[+] Generating t-SNE visualization...")
    perplexity = min(30, len(dist_matrix) - 1)
    tsne = TSNE(n_components=2, metric='precomputed', init='random', random_state=42, perplexity=perplexity)
    coords = tsne.fit_transform(dist_matrix)
    tsne_df = pd.DataFrame(coords, columns=['x', 'y'])
    tsne_df['Cluster'] = cluster_labels
    tsne_df['Sequence_ID'] = labels
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=tsne_df, x='x', y='y', hue='Cluster', palette='tab20', legend='full')
    plt.title("t-SNE of DTW Clustered Entropy Profiles")
    plt.savefig(f"{out_prefix}_tsne.png")
    print(f"[✔] t-SNE plot saved: {out_prefix}_tsne.png")
    plt.close()

    # Save metadata summary
    print("[+] Saving cluster metadata summary...")
    summary = df.groupby('Cluster').agg(
        n_sequences=('Sequence_ID', 'nunique'),
        entropy_mean=('Entropy', 'mean'),
        entropy_std=('Entropy', 'std'),
        avg_sv_prob=('SV_Prob', 'mean')
    ).reset_index()
    summary.to_csv(f"{out_prefix}_cluster_summary.tsv", sep='\t', index=False)
    print(f"[✔] Cluster summary saved: {out_prefix}_cluster_summary.tsv")

    # Save full table with labels
    df.to_csv(f"{out_prefix}_labeled.csv", index=False)
    print(f"[✔] Full table saved: {out_prefix}_labeled.csv")

# ========== CLI ENTRY ==========
def main():
    parser = argparse.ArgumentParser(description="Entropy Profile Clustering with DTW + Agglomerative")
    parser.add_argument("--input", required=True, help="Input CSV file with entropy profiles")
    parser.add_argument("--output-prefix", required=True, help="Prefix for all output files")
    args = parser.parse_args()

    df = load_entropy_table(args.input)
    cluster_and_visualize(df, out_prefix=args.output_prefix) 

if __name__ == "__main__":
    main()
    
# python3 entropy_clustering.py --input simulated_fragments_results.csv --output-prefix simulated_fragments_clustered    

