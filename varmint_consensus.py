#!/usr/bin/env python3
"""
Generate consensus FASTA files from yearly merged BAM files using varmint.
Iterates over all subtypes (H1N1, H3N2, H5N1), sources (clinical, wastewater),
and years (derived from BAM filenames in each bams_year directory).

Outputs:
  - Consensus FASTA: test_runs/consensus/{subtype}.{source}.{year}.fasta
  - Variants TSV:    test_runs/consensus/{subtype}.{source}.{year}.tsv
"""

# executable: python3 /data/tisza/analyses/crm/mutatea/varmint_consensus.py

import subprocess
from pathlib import Path

# --- adjust from here ---

VARMINT = "/mmfs1/home/u255582/.local/bin/varmint"
BASE_DIR = Path("/data/tisza/analyses/crm/mutatea/test_runs/yearly_time/mutatea")
OUT_DIR = Path("/data/tisza/analyses/crm/mutatea/test_runs/consensus")

SUBTYPES = ["H1N1", "H3N2", "H5N1"]

# source name -> bams_year subpath within each {subtype}_align directory
SOURCES = {
    "clinical":    Path("alignment_files/clinical/bams_year"),
    "wastewater":  Path("alignment_files/wastewater/bams_merged/bams_year"),
}

CONSENSUS_AF = 0.5   # allele frequency threshold for IUPAC ambiguity codes
THREADS = 4          # parallel threads per varmint run

# --- end of adjustable section ---

for subtype in SUBTYPES:
    subtype_dir = BASE_DIR / f"{subtype}_align"
    ref_dir = subtype_dir / "reference_files"

    fna_files = sorted(ref_dir.glob("*.fna"))
    gff_files = sorted(ref_dir.glob("*.gff"))

    if not fna_files or not gff_files:
        print(f"[WARN] Missing reference files for {subtype} in {ref_dir}, skipping.")
        continue

    ref_fasta = fna_files[0]
    ref_gff = gff_files[0]

    for source, bams_subpath in SOURCES.items():
        bams_dir = subtype_dir / bams_subpath

        if not bams_dir.exists():
            print(f"[WARN] {bams_dir} does not exist, skipping {subtype}/{source}.")
            continue

        OUT_DIR.mkdir(parents=True, exist_ok=True)

        bam_files = sorted(bams_dir.glob("*.sort.bam"))
        if not bam_files:
            print(f"[WARN] No BAM files found in {bams_dir}, skipping {subtype}/{source}.")
            continue

        for bam in bam_files:
            year = bam.name.replace(".sort.bam", "")
            stem = f"{subtype}.{source}.{year}"
            consensus_out = OUT_DIR / f"{stem}.fasta"
            variants_out = OUT_DIR / f"{stem}.tsv"

            cmd = [
                VARMINT,
                "--bam", str(bam),
                "--ref", str(ref_fasta),
                "--gff", str(ref_gff),
                "--out", str(variants_out),
                "--consensus", str(consensus_out),
                "--consensus-af", str(CONSENSUS_AF),
                "--threads", str(THREADS),
            ]

            print(f"[{subtype}/{source}/{year}] running varmint...")
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"[ERROR] {subtype}/{source}/{year}:\n{result.stderr.strip()}")
            else:
                print(f"[OK]    consensus -> {consensus_out}")
                print(f"[OK]    variants  -> {variants_out}")
