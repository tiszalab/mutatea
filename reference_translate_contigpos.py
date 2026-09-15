#!/usr/bin/env python3
"""
Translate reference FASTA nucleotide sequences to amino acid sequences
using CDS annotations from per-subtype GFF files

Input:  /data/tisza/analyses/crm/mutatea/test_input_data/clinical_input_data_{subtype}/*.fna
Output: /data/tisza/analyses/crm/mutatea/test_input_data/clinical_input_data_{subtype}_aa.fna

Each output file contains one AA record per annotated CDS, with headers:
  >{subtype}.{gene_name}
"""

# executable: python3 /data/tisza/analyses/crm/mutatea/reference_translate.py

import re
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq

# --- adjust from here ---
REFERENCE_DIR = Path(f"/data/tisza/analyses/crm/mutatea/test_input_data/clinical_input_data_H3N2")
OUT_DIR = REFERENCE_DIR / "amino_acids"
BASE_DIR = Path("/data/tisza/analyses/crm/mutatea/test_runs/yearly_time/mutatea")
SUBTYPES = ["H1N1", "H3N2", "H5N1"]

# --- end of adjustable section ---


def parse_cds_from_gff(gff_path):
    """
    Parse CDS features from a GFF3 file.

    Groups multi-exon CDS by the GFF 'ID' attribute (e.g. 'cds-UUB85550.1'),
    which is shared across all exons of the same protein.

    Returns:
        cds_id_to_parts : dict  cds_id -> sorted list of (contig, start0, end, phase)
        cds_id_to_gene  : dict  cds_id -> gene name string
    """
    cds_id_to_parts = defaultdict(list)
    cds_id_to_gene = {}

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue

            contig = fields[0]
            start0 = int(fields[3]) - 1   # GFF3 is 1-based → 0-based
            end    = int(fields[4])        # GFF3 end is inclusive; use as Python slice end
            phase  = int(fields[7]) if fields[7] != "." else 0
            attrs  = fields[8]

            cds_id_m = re.search(r'ID=([^;]+)', attrs)
            gene_m   = re.search(r'gene=([^;]+)', attrs)

            if not cds_id_m:
                continue

            cds_id    = cds_id_m.group(1)
            gene_name = gene_m.group(1) if gene_m else cds_id

            cds_id_to_parts[cds_id].append((contig, start0, end, phase))
            if cds_id not in cds_id_to_gene:
                cds_id_to_gene[cds_id] = gene_name

    # Sort each CDS's exons by start position
    for cds_id in cds_id_to_parts:
        cds_id_to_parts[cds_id].sort(key=lambda x: x[1])

    return cds_id_to_parts, cds_id_to_gene


def translate_cds(parts, seq_dict):
    """
    Concatenate CDS exon sequences from seq_dict and translate to amino acids.

    Applies the phase of the first exon as a 5' trim before translation.
    Returns the amino acid string, or None if a contig is missing.
    """
    first_phase = parts[0][3]
    nt_parts = []
    for contig, start0, end, _phase in parts:
        if contig not in seq_dict:
            return None
        nt_parts.append(str(seq_dict[contig].seq[start0:end]))

    nt_seq = "".join(nt_parts)[first_phase:]
    return str(Seq(nt_seq).translate(to_stop=True))


# ── main ──────────────────────────────────────────────────────────────────────

OUT_DIR.mkdir(parents=True, exist_ok=True)

for subtype in SUBTYPES:
    subtype_dir = BASE_DIR / f"{subtype}_align"
    gff_files = sorted((subtype_dir / "reference_files").glob("*.gff"))
    if not gff_files:
        print(f"[WARN] No GFF for {subtype}, skipping.")
        continue

    cds_parts, cds_genes = parse_cds_from_gff(gff_files[0])

    for fasta_path in sorted(REFERENCE_DIR.glob(f"*.fna")):
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

        records = []
        for cds_id, exons in cds_parts.items():
            gene_name = cds_genes[cds_id]
            aa_seq = translate_cds(exons, seq_dict)
            if aa_seq is None:
                print(f"[WARN] contig missing for {cds_id} in {fasta_path.name}, skipping.")
                continue
            records.append(f">{subtype}.{gene_name}\n{aa_seq}")

        out_path = OUT_DIR / f"{subtype}.aa.fasta"
        with open(out_path, "w") as fh:
            fh.write("\n".join(records) + "\n")

        print(f"[OK] {out_path.name}  ({len(records)} proteins)")
