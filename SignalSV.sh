#!/bin/bash

# ========== USAGE ========== #
usage() {
    echo ""
    echo "Usage: $0 \\"
    echo "  --prefix <run_prefix> \\"
    echo "  --fasta <reference.fasta> \\"
    echo "  --mask <spaced_mask> \\"
    echo "  --window <window_size> \\"
    echo "  --step <step_size> \\"
    echo "  --zmode <global|sequence> \\"
    echo "  --num-svs <int> \\"
    echo "  --slack <bp_slack> \\"
    echo "  --prob-thresh <float>"
    echo ""
    exit 1
}

# ========== DEFAULTS ========== #
ZFLAG="--zscore"

# ========== PARSE ARGS ========== #
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --prefix) PREFIX="$2"; shift ;;
        --fasta) FASTA="$2"; shift ;;
        --mask) MASK="$2"; shift ;;
        --window) WINDOW="$2"; shift ;;
        --step) STEP="$2"; shift ;;
        --zmode) ZMODE="$2"; shift ;;
        --num-svs) NUM_SVS="$2"; shift ;;
        --slack) SLACK="$2"; shift ;;
        --prob-thresh) PROB_THRESH="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check required
if [[ -z "$PREFIX" || -z "$FASTA" || -z "$MASK" || -z "$WINDOW" || -z "$STEP" || -z "$ZMODE" || -z "$NUM_SVS" || -z "$SLACK" || -z "$PROB_THRESH" ]]; then
    echo "[!] Missing required arguments!"
    usage
fi

# ========== OUTPUT DIR ========== #
OUTDIR="${PREFIX}_SignalSV"
mkdir -p "$OUTDIR"

echo "[✔] Running SignalSV pipeline with prefix: $PREFIX"
echo "[✔] Output folder: $OUTDIR"

# ========== STEP 1: Simulate SVs ========== #
echo "[1/5] Simulating SVs..."
python3 sv_gold_standard_simulator.py \
  --reference "$FASTA" \
  --output-fasta "$OUTDIR/simulated_with_svs.fasta" \
  --output-truth "$OUTDIR/simulated_sv_truth.tsv" \
  --header "$PREFIX" \
  --num-svs "$NUM_SVS"

# ========== STEP 2: Entropy Profiling ========== #
echo "[2/5] Entropy profiling..."
python3 SpacekmerEntropy.py \
  --fasta "$OUTDIR/simulated_with_svs.fasta" \
  --mask "$MASK" \
  --window "$WINDOW" \
  --step "$STEP" \
  --output "$OUTDIR/entropy_results.csv" \
  --plot \
  --svtrack \
  $ZFLAG \
  --zmode "$ZMODE"

# FIX: move the entropy plot to the output folder
if [[ -f "entropy_plot.png" ]]; then
    mv entropy_plot.png "$OUTDIR/entropy_plot.png"
fi

# ========== STEP 3: Clustering + KDE ========== #
echo "[3/5] Clustering + KDE probability modeling..."
python3 entropy_clustering.py \
  --input "$OUTDIR/entropy_results.csv" \
  --output-prefix "$OUTDIR/clustered"

# ========== STEP 4: Map Truth ========== #
echo "[4/5] Mapping truth SVs to entropy windows..."
python3 map_truth_to_windows.py \
  --truth "$OUTDIR/simulated_sv_truth.tsv" \
  --window "$WINDOW" \
  --step "$STEP" \
  --output "$OUTDIR/simulated_sv_truth_windows.tsv"

# ========== STEP 5: Evaluate Probabilities ==========
echo "[5/5] Evaluating SV probability distributions..."
python3 evaluate_prob_distributions.py \
  --predicted "$OUTDIR/clustered_labeled.csv" \
  --truth "$OUTDIR/simulated_sv_truth_windows.tsv" \
  --output "$OUTDIR/sv_probability_distribution.tsv" \
  --slack "$SLACK"

echo "[✔] All done! Check output at: $OUTDIR"

# chmod +x SignalSV.sh
# ./SignalSV.sh --prefix test_run --fasta tp53_ref.fasta --mask 1101 --window 50 --step 25 --zmode global --num-svs 15 --slack 50 --prob-thresh 0.75

