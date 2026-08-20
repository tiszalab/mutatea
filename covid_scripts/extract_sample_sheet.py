#!/usr/bin/env python3

# crm: need to add bam file extraction

# parse given file directories to identify file paths of associated reads
# save file paths of the reads as a txt file in specified sample sheet format

"""
Executable code:
python extract_sample_sheet.py \
    --pathogen sars_cov2 \
    --input_dir /data/contract/TEPHI
"""

import argparse
import os
from pathlib import Path
from typing import List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Finding file paths of pathogen reads")
    p.add_argument("--pathogen", required=True, help="Pathogen name")
    p.add_argument("--input_dir", required=True, help="Directory containing pathogen reads")
    p.add_argument("--output", "-o", default=None, help="Output sample sheet filename (default: ./<pathogen>_sample_sheet.txt)")
    return p.parse_args()


# recursive search for all fasta/fastq files matching the pathogen name
def find_wastewater_reads(input_dir: str, pathogen: str) -> List[str]:
    extensions = ('.fasta', '.fastq', '.fastq.gz', '.sort.bam')

    # excluding folders I know contain junk reads
    excluded_dirs = {'smk_v016_oldrun', 'smk_v016_undetermined', 'smk_v016_undetermined_samples'}

    all_files = []
    for root, dirs, files in os.walk(input_dir):
        # want to exclude the random test folder from the pools
        if any(part.startswith('test_p') for part in Path(root).parts):
            dirs.clear()
            continue
        if any(part in excluded_dirs for part in Path(root).parts):
            dirs.clear()
            continue
        for fname in files:
            if pathogen in fname and any(fname.endswith(ext) for ext in extensions):
                all_files.append(os.path.join(root, fname))
    if not all_files:
        raise FileNotFoundError(f"No FASTA, FASTQ, FASTQ.GZ, or SORT.BAM files found for {pathogen} in {input_dir}")
    return sorted(all_files)


def main():
    args = parse_args()

    input_dir = args.input_dir
    pathogen = args.pathogen
    # set output file
    outputfile = Path(args.output) if args.output else Path(f"./{pathogen}_sample_sheet.txt")

    # find all matching read files
    all_files = find_wastewater_reads(input_dir, pathogen)

    # verify files exist
    valid_files = [f for f in all_files if os.path.exists(f)]
    missing = [f for f in all_files if not os.path.exists(f)]
    if missing:
        print(f"{len(missing)} file(s) could not be verified:")
        for f in missing:
            print(f"  {f}")

    # identify sample type, organize file paths into rows (different by sample type)
    pairs: dict = {}
    for file_path in sorted(valid_files):
        fname = os.path.basename(file_path)
        sample_id = fname[:fname.index(pathogen)].rstrip('._-').split('.')[0]
        if file_path.endswith('.fasta'):
            sample_type = 'fasta'
        elif file_path.endswith('.sort.bam'):
            sample_type = 'bam'
        elif file_path.endswith('.fastq'):
            sample_type = 'fastq'
        elif file_path.endswith('.fastq.gz'):
            sample_type = 'fastq.gz'
        
        # organize the samples into key; include source directory so a duplicate sample ID from a different dir can be kept as a separate row
        key = (sample_id, sample_type, os.path.dirname(file_path))
        if key not in pairs:
            pairs[key] = [None, None]
        if sample_type == 'bam':
            pairs[key][0] = file_path
        elif '.R1.' in os.path.basename(file_path) or '_R1.' in os.path.basename(file_path) or '_1.' in os.path.basename(file_path):
            pairs[key][0] = file_path
        elif '.R2.' in os.path.basename(file_path) or '_R2.' in os.path.basename(file_path) or '_2.' in os.path.basename(file_path):
            pairs[key][1] = file_path
        else:
            pairs[key][0] = file_path

    # check for duplicate sample_ids (same id from multiple directories) and warn
    sample_ids = [sid for (sid, _st, _src) in pairs]
    seen = set()
    duplicates = [sid for sid in sample_ids if sid in seen or seen.add(sid)]
    if duplicates:
        print(f"WARNING: Duplicate sample_id(s) found: {', '.join(sorted(set(duplicates)))}")

    # filter to only allow one type of sample per sample sheet
    n_bam    = sum(1 for (_, st, _src), _        in pairs.items() if st == 'bam')
    n_paired = sum(1 for (_, st, _src), (r1, r2) in pairs.items() if st != 'bam' and r1 and r2)
    n_single = sum(1 for (_, st, _src), (r1, r2) in pairs.items() if st != 'bam' and not (r1 and r2))

    read_types_present = []
    if n_bam > 0:
        read_types_present.append('bam')
    if n_paired > 0:
        read_types_present.append('paired')
    if n_single > 0:
        read_types_present.append('single')
    if len(read_types_present) > 1:
        print(f"WARNING: Mixed sample types found: {n_bam} BAM files, {n_paired} paired reads, {n_single} single reads. The sample sheet should only contain one sample type.")

    if n_bam:
        print(f"Found {n_bam} BAM files for {pathogen} from input directory: {input_dir}")
    if n_paired:
        print(f"Found {n_paired} paired reads for {pathogen} from input directory: {input_dir}")
    if n_single:
        print(f"Found {n_single} single reads for {pathogen} from input directory: {input_dir}")

    # write sample sheet
    with open(outputfile, 'w') as f:
        f.write(f"{pathogen} samples identified from input directory: {input_dir}\n")
        f.write("# sample_id\tsample_type\tfile_path_1\tfile_path_2\n")
        # crm: need to test to make sure if there is one paired read it isn't read as as single read
        for (sample_id, stype, _src_dir), (r1, r2) in sorted(pairs.items()):
            if stype == 'bam':
                sample_type = 'bam'
            else:
                sample_type = 'paired' if r1 and r2 else 'single'
            f.write(f"{sample_id}\t{sample_type}\t{r1 or ''}\t{r2 or ''}\n")

    print("Please check the outputted file to confirm you are only including samples you want processed")

if __name__ == "__main__":
    main()
