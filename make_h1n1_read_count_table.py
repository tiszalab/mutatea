#!/usr/bin/env python3
"""
Build a combined table per subtype of:
  - Wastewater raw read pairs per month (from H*N*_sample_sheet.txt)
  - Clinical consensus sequence counts per month (from monthly FASTA files)

Produces one CSV per subtype: H1N1, H3N2, H5N1.
"""

import csv
import os
import sys
from collections import defaultdict

BASE_DIR    = "/data/tisza/analyses/crm/mutatea"
RUN_DIR     = f"{BASE_DIR}/test_runs/awm_new_flu/mutatea"
SUBTYPES    = ["H1N1", "H3N2", "H5N1"]


def build_table(subtype):
    ww_sample_sheet    = f"{BASE_DIR}/sample_sheet/{subtype}_sample_sheet.txt"
    ww_metadata_csv    = f"{RUN_DIR}/{subtype}_align/metadata_files/metadata_wastewater_combined.csv"
    clinical_fasta_dir = f"{RUN_DIR}/{subtype}_align/alignment_files/clinical/fastas_month"
    output_csv         = f"{BASE_DIR}/{subtype}_read_counts_by_month.csv"

    print(f"\n{'='*50}", flush=True)
    print(f"  {subtype}", flush=True)
    print(f"{'='*50}", flush=True)

    # ── 1. Parse sample sheet ─────────────────────────────────────────────
    print("Parsing wastewater sample sheet...", flush=True)
    samples = []
    with open(ww_sample_sheet) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or line.startswith(subtype):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                samples.append((parts[0], parts[2]))
    print(f"  {len(samples)} samples in sample sheet")

    # ── 2. Parse metadata: sample_id -> Month_Year ────────────────────────
    print("Parsing wastewater metadata...", flush=True)
    sample_to_month = {}
    with open(ww_metadata_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample_to_month[row["Sample_ID"]] = row["Month_Year"]
    print(f"  {len(sample_to_month)} samples in metadata")

    # ── 3. Count wastewater samples per month ─────────────────────────────
    print("Counting wastewater samples...", flush=True)
    ww_samples_by_month = defaultdict(int)
    missing_files = 0
    no_metadata   = 0

    for sid, r1_path in samples:
        if not os.path.isfile(r1_path):
            missing_files += 1
            continue
        month = sample_to_month.get(sid)
        if month is None:
            no_metadata += 1
            continue
        ww_samples_by_month[month] += 1

    if missing_files:
        print(f"  WARNING: {missing_files} R1 files not found on disk", file=sys.stderr)
    if no_metadata:
        print(f"  NOTE: {no_metadata} samples skipped (no metadata match)")
    print(f"  Done. {len(ww_samples_by_month)} months with wastewater data")

    # ── 4. Count clinical consensus sequences per month ───────────────────
    print("Counting clinical consensus sequences...", flush=True)
    clinical_seqs_by_month = {}
    for fname in sorted(os.listdir(clinical_fasta_dir)):
        if not fname.endswith(".fasta"):
            continue
        month = fname.replace(".fasta", "")
        fpath = os.path.join(clinical_fasta_dir, fname)
        n_seqs = sum(1 for line in open(fpath) if line.startswith(">"))
        clinical_seqs_by_month[month] = n_seqs
    print(f"  Done. {len(clinical_seqs_by_month)} months with clinical data")

    # ── 5. Build rows ─────────────────────────────────────────────────────
    all_months = sorted(set(ww_samples_by_month.keys()) | set(clinical_seqs_by_month.keys()))
    rows = []
    for month in all_months:
        parts = month.split("_")
        mm   = parts[0] if len(parts) == 2 else ""
        yyyy = parts[1] if len(parts) == 2 else ""
        rows.append({
            "Month_Year":              month,
            "Month":                   mm,
            "Year":                    yyyy,
            "WW_samples":              ww_samples_by_month.get(month, 0),
            "Clinical_consensus_seqs": clinical_seqs_by_month.get(month, ""),
        })

    # ── 6. Total row ──────────────────────────────────────────────────────
    total_ww   = sum(r["WW_samples"] for r in rows)
    total_clin = sum(r["Clinical_consensus_seqs"] for r in rows if r["Clinical_consensus_seqs"] != "")
    rows.append({
        "Month_Year":              "Total",
        "Month":                   "",
        "Year":                    "",
        "WW_samples":              total_ww,
        "Clinical_consensus_seqs": total_clin,
    })

    # ── 7. Write CSV ──────────────────────────────────────────────────────
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved: {output_csv}")
    print(f"Rows: {len(rows) - 1} months + 1 total row")

    # ── 8. Preview ────────────────────────────────────────────────────────
    print("\nPreview:")
    header = f"{'Month_Year':<12} {'WW_samples':>10} {'Clinical_consensus_seqs':>23}"
    print(header)
    print("-" * len(header))
    for row in rows:
        clin = str(row["Clinical_consensus_seqs"]) if row["Clinical_consensus_seqs"] != "" else "N/A"
        print(f"{row['Month_Year']:<12} {row['WW_samples']:>10} {clin:>23}")


for subtype in SUBTYPES:
    build_table(subtype)
