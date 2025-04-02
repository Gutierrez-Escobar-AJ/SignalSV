#!/usr/bin/env python3

"""
Pilot version of SpacekmerSDI adapted for gene fragments.
Performs sliding window spaced k-mer entropy profiling.
Includes per-sequence or global Z-score, SV export, dynamic Z threshold, and plotting with SV annotations.
Author: Gutierrez-Escobar-AJ
"""

import argparse
import os
from typing import List
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

# ========== SPACED K-MER EXTRACTION ==========
def spaced_kmers(seq: str, mask: str) -> List[str]:
    k = len(mask)
    valid_kmers = []
    valid_nucs = {'A', 'T', 'C', 'G'}

    for i in range(len(seq) - k + 1):
        try:
            kmer = [seq[i + j] for j in range(k) if mask[j] == '1']
            if all(base.upper() in valid_nucs for base in kmer):
                valid_kmers.append(''.join(kmer))
        except IndexError:
            continue

    return valid_kmers

# ========== ENTROPY CALCULATION ==========
def compute_entropy(kmers: List[str]) -> float:
    if not kmers:
        return 0.0
    counts = Counter(kmers)
    freqs = np.array(list(counts.values())) / len(kmers)
    entropy = -np.sum(freqs * np.log2(freqs + 1e-10))
    return entropy

# ========== SLIDING WINDOW PROCESSING ==========
def process_sequence_windows(seq: str, mask: str, window_size: int, step_size: int,
                             seq_id: str) -> List[List]:
    results = []
    for start in range(0, len(seq) - window_size + 1, step_size):
        window = seq[start:start + window_size]
        kmers = spaced_kmers(window, mask)
        entropy = compute_entropy(kmers)
        skipped = (len(window) - len(mask) + 1) - len(kmers)
        results.append([seq_id, start, start + window_size, entropy, len(kmers), skipped])
    return results

# ========== MAIN PIPELINE ==========
def process_fasta(fasta_path: str, mask: str, window_size: int, step_size: int) -> pd.DataFrame:
    all_results = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        seq_id = record.id
        res = process_sequence_windows(seq, mask, window_size, step_size, seq_id)
        all_results.extend(res)

    df = pd.DataFrame(all_results, columns=[
        'Sequence_ID', 'Window_Start', 'Window_End',
        'Entropy', 'SpacedKmer_Count', 'Skipped_Kmers'])
    df["Window_Mid"] = (df["Window_Start"] + df["Window_End"]) // 2
    return df

# ========== Z-SCORE DETECTION ==========
def add_zscores(df: pd.DataFrame, z_thresh: float = 2.0, mode: str = 'global') -> pd.DataFrame:
    if mode == 'global':
        mean = df['Entropy'].mean()
        std = df['Entropy'].std(ddof=0)
        df['Zscore'] = (df['Entropy'] - mean) / std if std > 0 else 0.0
    elif mode == 'sequence':
        zscores = []
        for sid in df['Sequence_ID'].unique():
            subset = df[df['Sequence_ID'] == sid]
            mean = subset['Entropy'].mean()
            std = subset['Entropy'].std(ddof=0)
            subset_zscores = (subset['Entropy'] - mean) / std if std > 0 else 0.0
            zscores.extend(subset_zscores)
        df['Zscore'] = zscores
    
    df['SV_Flag'] = df['Zscore'].apply(
        lambda z: 'high_entropy_region' if z > z_thresh else 'low_entropy_region' if z < -z_thresh else 'normal'
    )
    return df

# ========== SV CALL EXPORT ==========
def export_sv_calls(df: pd.DataFrame, output_file: str, z_thresh: float = 2.0):
    sv_df = df[df['SV_Flag'] != 'normal'][[
        'Sequence_ID', 'Window_Start', 'Window_End', 'Entropy', 'Zscore', 'SV_Flag']]
    sv_df.to_csv(output_file, sep='\t', index=False)
    print(f"SV calls exported to: {output_file}")

# ========== PLOTTING ==========
def plot_entropy(df: pd.DataFrame, output_path: str = "entropy_plot.png", show_sv: bool = False) -> None:
    plt.figure(figsize=(15, 10))
    for seq_id, group in df.groupby("Sequence_ID"):
        # Convert pandas Series to numpy arrays
        x = group["Window_Mid"].values
        y = group["Entropy"].values
        plt.plot(x, y, label=seq_id, marker='o', markersize=3)
        
        if show_sv:
            flagged = group[group['SV_Flag'] != 'normal']
            for _, row in flagged.iterrows():
                alpha = 0.2 if row['SV_Flag'] == 'low_entropy_region' else 0.3
                color = 'blue' if row['SV_Flag'] == 'low_entropy_region' else 'red'
                plt.axvspan(row['Window_Start'], row['Window_End'], color=color, alpha=alpha)

    plt.xlabel("Window Midpoint (bp)")
    plt.ylabel("Entropy (bits)")
    plt.title("Spaced k-mer Entropy Profiles")
    plt.legend(title="Sequence", loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Entropy plot saved to: {output_path}")

# ========== CLI ENTRY ==========
def main():
    parser = argparse.ArgumentParser(description="Sliding window spaced k-mer entropy profiler")
    parser.add_argument("--fasta", required=True, help="Input FASTA file with gene fragments")
    parser.add_argument("--mask", required=True, help="Binary mask for spaced k-mer (e.g., 1101)")
    parser.add_argument("--window", type=int, default=50, help="Window size")
    parser.add_argument("--step", type=int, default=25, help="Step size between windows")
    parser.add_argument("--output", default="entropy_profile.csv", help="Output CSV file")
    parser.add_argument("--plot", action="store_true", help="Generate entropy line plot")
    parser.add_argument("--svtrack", action="store_true", help="Overlay SV regions on entropy plot")
    parser.add_argument("--zscore", action="store_true", help="Enable z-score calculation and SV export")
    parser.add_argument("--zthresh", type=float, default=2.0, 
                       help="Z-score threshold for SV detection (default: 2.0)")
    parser.add_argument("--zmode", choices=['global', 'sequence'], default='global',
                       help="Z-score mode. 'global': all windows; 'sequence': per-sequence")
    
    args = parser.parse_args()

    # Validations
    if not all(c in {'0', '1'} for c in args.mask):
        raise ValueError("Mask must contain only 0s and 1s")
    if args.window < len(args.mask):
        raise ValueError(f"Window size ({args.window}) must be ≥ mask length ({len(args.mask)})")

    df = process_fasta(args.fasta, args.mask, args.window, args.step)

    if args.zscore:
        df = add_zscores(df, z_thresh=args.zthresh, mode=args.zmode)
        sv_output = os.path.splitext(args.output)[0] + "_SVcalls.tsv"
        export_sv_calls(df, sv_output, z_thresh=args.zthresh)

    df.to_csv(args.output, index=False)
    print(f"Entropy profile saved to: {args.output}")

    if args.plot:
        plot_entropy(df, show_sv=args.svtrack)

if __name__ == "__main__":
    main()
    
# python3 SpacekmerEntropy.py --fasta fragments.fasta --mask 1101 --window 50 --step 25 --output results.csv --plot --svtrack --zscore --zmode global
