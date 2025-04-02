#!/usr/bin/env python3

"""
Purpose: Aligns ground truth SVs to entropy-based sliding windows.
Author: Gutierrez-Escobar-AJ
"""

import argparse
import pandas as pd

def map_sv_to_windows(truth_file, window_size, step_size, output_file):
    df = pd.read_csv(truth_file, sep='\t')
    records = []

    for _, row in df.iterrows():
        seq_id = row["Sequence_ID"]
        pos = int(row["Position"])
        svtype = row["SVTYPE"]
        svlen = row["SVLEN"]
        ref = row["REF"]
        alt = row["ALT"]

        # Loop over windows that might include the SV position
        for start in range(0, pos + window_size, step_size):
            end = start + window_size
            if start <= pos < end:
                records.append({
                    "Sequence_ID": seq_id,
                    "Position": pos,
                    "REF": ref,
                    "ALT": alt,
                    "SVTYPE": svtype,
                    "SVLEN": svlen,
                    "Window_Start": start,
                    "Window_End": end
                })
                break  # Stop after finding the first matching window

    out_df = pd.DataFrame(records)
    out_df.to_csv(output_file, sep='\t', index=False)
    print(f"[✔] Mapped SVs to windows: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map SVs to entropy windows")
    parser.add_argument("--truth", required=True, help="Input SV truth file (simulated_sv_truth.tsv)")
    parser.add_argument("--window", type=int, default=50, help="Window size")
    parser.add_argument("--step", type=int, default=25, help="Step size")
    parser.add_argument("--output", default="simulated_sv_truth_windows.tsv", help="Output file")
    args = parser.parse_args()

    map_sv_to_windows(args.truth, args.window, args.step, args.output)

# python3 map_truth_to_windows.py --truth simulated_sv_truth.tsv --window 50 --step 25 --output simulated_sv_truth_windows.tsv
