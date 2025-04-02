#!/usr/bin/env python3

"""
SignalSV | Gold Standard Simulator

Simulates small structural variants (SVs) including insertions, deletions, and substitutions
from a reference FASTA using BioPython.

Outputs:
1. A FASTA with reference and SV-mutated fragments
2. A TSV file with the ground truth SVs

Author: Gutierrez-Escobar-AJ
"""

import random
import argparse
from Bio import SeqIO

def simulate_svs(ref_seq, header_prefix, num_svs, sv_size_range):
    length = len(ref_seq)
    fragments = [(f">{header_prefix}_Ref\n", ref_seq)]
    truth_table = []

    for idx in range(num_svs):
        pos = random.randint(100, length - 100)
        svlen = random.randint(*sv_size_range)
        svtype = random.choice(["INS", "DEL", "SUB"])
        new_seq = list(ref_seq)
        label, alt_seq = None, ""

        if svtype == "INS":
            ins_seq = ''.join(random.choices("ACGT", k=svlen))
            new_seq[pos] = new_seq[pos] + ins_seq
            label = f"{header_prefix}_Ins{svlen}bp_{pos+1}"
            alt_seq = new_seq[pos]
            truth_table.append((label, pos+1, ref_seq[pos], alt_seq, svtype, svlen))

        elif svtype == "DEL":
            if pos + svlen < length:
                ref = ''.join(new_seq[pos:pos+svlen])
                new_seq[pos] = ref[0]
                for i in range(1, svlen):
                    new_seq[pos + i] = ""
                label = f"{header_prefix}_Del{svlen}bp_{pos+1}"
                truth_table.append((label, pos+1, ref, ref[0], svtype, svlen))

        elif svtype == "SUB":
            if pos + svlen < length:
                ref = ''.join(new_seq[pos:pos+svlen])
                alt = ''.join(random.choices("ACGT", k=svlen))
                for i in range(svlen):
                    new_seq[pos + i] = alt[i]
                label = f"{header_prefix}_Sub{svlen}bp_{pos+1}"
                truth_table.append((label, pos+1, ref, alt, svtype, svlen))

        if label:
            mutated = ''.join(new_seq)
            fragments.append((f">{label}\n", mutated))

    return fragments, truth_table

def write_fasta(fragments, output_path):
    with open(output_path, "w") as f:
        for header, seq in fragments:
            f.write(header)
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def write_truth_table(truth_table, output_path):
    with open(output_path, "w") as out:
        out.write("Sequence_ID\tPosition\tREF\tALT\tSVTYPE\tSVLEN\n")
        for row in truth_table:
            out.write("\t".join(map(str, row)) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Simulate small SVs from a reference FASTA.")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--output-fasta", default="simulated_sv_fragments.fasta", help="Output FASTA with mutated sequences")
    parser.add_argument("--output-truth", default="simulated_sv_truth.tsv", help="Ground truth TSV for SVs")
    parser.add_argument("--header", default="TP53", help="Prefix for sequence headers")
    parser.add_argument("--num-svs", type=int, default=100, help="Number of simulated SVs")
    parser.add_argument("--min-size", type=int, default=1, help="Minimum SV size")
    parser.add_argument("--max-size", type=int, default=50, help="Maximum SV size")
    args = parser.parse_args()

    record = next(SeqIO.parse(args.reference, "fasta"))
    ref_seq = str(record.seq).upper()

    fragments, truth_table = simulate_svs(
        ref_seq=ref_seq,
        header_prefix=args.header,
        num_svs=args.num_svs,
        sv_size_range=(args.min_size, args.max_size)
    )

    write_fasta(fragments, args.output_fasta)
    write_truth_table(truth_table, args.output_truth)

    print("[✔] Simulation complete:")
    print(f"     - FASTA: {args.output_fasta}")
    print(f"     - TSV  : {args.output_truth}")

if __name__ == "__main__":
    main()

# python3 sv_gold_standard_simulator.py --reference tp53_ref.fasta --output-fasta simulated_sv_fragments.fasta --output-truth simulated_sv_truth.tsv --header TP53 --num-svs 30 --min-size 2 --max-size 25
