#!/usr/bin/env python3

"""
Purpose: Visualizes the distribution of predicted SV probabilities around gold standard SVs.
Author: Gutierrez-Escobar-AJ
"""

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# ========== FUNCTION ==========
def evaluate_probabilities(pred_path, truth_path, output_path, slack=50):
    # Load predicted and truth
    predicted = pd.read_csv(pred_path)
    truth = pd.read_csv(truth_path, sep='\t')

    # Build list of truth positions per Sequence_ID
    truth_positions = truth.groupby('Sequence_ID')['Position'].apply(list).to_dict()

    labels = []
    for _, row in predicted.iterrows():
        sid = row['Sequence_ID']
        mid = row['Window_Mid']
        found = False

        if sid in truth_positions:
            for pos in truth_positions[sid]:
                if abs(mid - pos) <= slack:
                    found = True
                    break
        labels.append("NearTruth" if found else "Other")

    predicted['Label'] = labels

    # Save mapped probabilities
    predicted.to_csv(output_path, sep='\t', index=False)
    print(f"[✔] Probability mapping saved: {output_path}")

    # Plot KDE with quantiles
    plt.figure(figsize=(10, 7))
    sns.kdeplot(data=predicted, x='SV_Prob', fill=True, alpha=0.5)

    # Draw quantiles
    for q in [0.25, 0.5, 0.75]:
        x = predicted['SV_Prob'].quantile(q)
        style = ':' if q != 0.5 else '--'
        plt.axvline(x=x, linestyle=style, color='k', label=f"{q:.2f} Quantile")

    plt.title("SV Probability Density: NearTruth vs Other")
    plt.xlabel("SV_Prob")
    plt.legend()
    plt.savefig(output_path.replace('.tsv', '.png'))
    plt.close()

# ========== MAIN ==========
def main():
    parser = argparse.ArgumentParser(description="Evaluate SV probability distributions near gold-standard SVs")
    parser.add_argument('--predicted', required=True, help='CSV of predicted SVs with SV_Prob and Window_Mid')
    parser.add_argument('--truth', required=True, help='TSV of gold standard SVs with Position')
    parser.add_argument('--output', required=True, help='Output TSV for merged probability data')
    parser.add_argument('--slack', type=int, default=50, help='Slack window in bp')
    args = parser.parse_args()

    evaluate_probabilities(
        pred_path=args.predicted,
        truth_path=args.truth,
        output_path=args.output,
        slack=args.slack
    )

if __name__ == '__main__':
    main()


# python3 evaluate_prob_distributions.py --predicted simulated_fragments_clustered_labeled.csv --truth simulated_sv_truth_windows.tsv --output sv_probability_distribution.tsv --slack 50

