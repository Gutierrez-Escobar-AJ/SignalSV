
# SignalSV

**SignalSV** is a minimalist, Python-powered framework for the detection and benchmarking of **small structural variants (SVs)** using entropy profiling, clustering, and probabilistic modeling. Designed for speed, modularity, and scientific rigor, SignalSV simulates gold-standard SVs, profiles spaced k-mer entropy, clusters profiles, and evaluates SV predictions — all in one streamlined pipeline.

<p align="center">
  <img src="logo/signal_sv_logo.png" alt="SignalSV Logo" width="220"/>
</p>

---

## 🚀 Features

- ⚙️ SV simulation with ground-truth annotations (insertions, deletions, substitutions)
- 🧬 Spaced k-mer entropy profiling of sequences
- 📊 Z-score–based detection of entropy anomalies
- 🧠 DTW + Ward-linkage clustering of profiles
- 🌡️ KDE-based SV probability estimation
- 🧪 Ground-truth benchmarking with soft-window matching
- 📈 Automated t-SNE plots, heatmaps, and tree visualizations
- 🌐 HTML dashboard (optional) for elegant summary reporting

---

## 📦 Installation

SignalSV requires Python ≥ 3.10. We recommend using Conda:

1. **Clone the repository**
   ```bash
   git clone https://github.com/your-username/SignalSV.git
   cd SignalSV
   ```

2. **Create the environment**
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate the environment**
   ```bash
   conda activate SignalSV
   ```

---

## ⚡ Quick Start

Once installed, run the full pipeline using:

```bash
./SignalSV.sh \
  --prefix tp53_run \
  --fasta example_data/tp53_ref.fasta \
  --mask 1101 \
  --window 50 \
  --step 25 \
  --zmode global \
  --num-svs 15 \
  --slack 50 \
  --prob-thresh 0.75
```

All results will be saved to a folder: `tp53_run_SignalSV/`.

---

## 📁 Output Overview

| File | Description |
|------|-------------|
| `simulated_with_svs.fasta` | FASTA with injected small SVs |
| `simulated_sv_truth.tsv` | Ground truth for the SV simulation |
| `entropy_results.csv` | Raw entropy profiles |
| `entropy_results_SVcalls.tsv` | Predicted SV windows (Z-score–based) |
| `clustered_labeled.csv` | Clustered profiles with assigned SV probabilities |
| `clustered.nwk` | Phylogenetic-style tree based on entropy profiles |
| `clustered_heatmap.png`, `clustered_tsne.png` | Visualizations of cluster structure |
| `sv_probability_distribution.tsv` | Probability mapping of predicted SVs vs. truth |
| `entropy_plot.png` | Whole-genome entropy plot (optional) |

---

## 🧪 Example Data

We provide a reference sequence in:
```
example_data/tp53_ref.fasta
```

Use this to quickly explore the pipeline and test all features.

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---

## 🙌 Acknowledgements

SignalSV is developed with ❤️ by researchers passionate about structural variant detection in microbial and genomic systems.
