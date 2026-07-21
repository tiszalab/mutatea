#!/usr/bin/env python3

# crm: need to add bam file extraction
# crm: need to add a filter to only allow one type of sample per sample sheet
# crm: change all inputs for pathogen to not be case-sensitive

# parse given file directories to identify file paths of associated reads
# save file paths of the reads as a txt file in specified sample sheet format

"""
Executable code:
python extract_sample_sheet.py \
    --pathogen H1N1 \
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
    return p.parse_args()


# recursive search for all fasta/fastq files matching the pathogen name
def find_wastewater_reads(input_dir: str, pathogen: str) -> List[str]:
    extensions = ('.fasta', '.fastq', '.fastq.gz')

    all_files = []
    for root, dirs, files in os.walk(input_dir):
        # crm: want to exclude that random test folder from the pools
        if any(part.startswith('test_p') for part in Path(root).parts):
            dirs.clear()
            continue
        # crm: making sure the pathogen input is not case sensitive
        for fname in files:
            if pathogen.lower() in fname.lower() and any(fname.endswith(ext) for ext in extensions):
                all_files.append(os.path.join(root, fname))
    if not all_files:
        raise FileNotFoundError(f"No FASTA, FASTQ, or FASTQ.GZ files found for {pathogen} in {input_dir}")
    return sorted(all_files)


def main():
    args = parse_args()

    input_dir  = args.input_dir
    pathogen   = args.pathogen
    # set default output file 
    outputfile = Path(f"./{pathogen}_sample_sheet.txt")

    # find all matching read files
    all_files = find_wastewater_reads(input_dir, pathogen)

    # verify files exist
    valid_files = [f for f in all_files if os.path.exists(f)]
    missing = [f for f in all_files if not os.path.exists(f)]
    if missing:
        print(f"{len(missing)} file(s) could not be verified:")
        for f in missing:
            print(f"  {f}")

    # pair R1 the R2 reads into same row
    pairs: dict = {}
    for file_path in sorted(valid_files):
        # crm: making sure the pathogen input is not case sensitive
        fname_lower = os.path.basename(file_path).lower()
        sample_id = os.path.basename(file_path)[:fname_lower.index(pathogen.lower())].rstrip('._-').split('.')[0]
        sample_type = 'fasta' if file_path.endswith('.fasta') else 'fastq'
        key = (sample_id, sample_type)
        if key not in pairs:
            pairs[key] = [None, None]
        if '.R1.' in os.path.basename(file_path) or '_R1.' in os.path.basename(file_path):
            pairs[key][0] = file_path
        elif '.R2.' in os.path.basename(file_path) or '_R2.' in os.path.basename(file_path):
            pairs[key][1] = file_path
        else:
            pairs[key][0] = file_path

    n_paired = sum(1 for r1, r2 in pairs.values() if r1 and r2)
    n_single = sum(1 for r1, r2 in pairs.values() if not (r1 and r2))
    if n_paired:
        print(f"Found {n_paired} paired reads for {pathogen} from input directory: {input_dir}")
    if n_single:
        print(f"Found {n_single} single reads for {pathogen} from input directory: {input_dir}")

    # write sample sheet
    with open(outputfile, 'w') as f:
        f.write(f"{pathogen} samples identified from input directory: {input_dir}\n")
        f.write("# sample_id\tsample_type\tfile_path_1\tfile_path_2\n")
        # crm: need to test to make sure if there is one paired read it isn't read as as single read
        for (sample_id, _), (r1, r2) in sorted(pairs.items()):
            sample_type = 'paired' if r1 and r2 else 'single'
            f.write(f"{sample_id}\t{sample_type}\t{r1 or ''}\t{r2 or ''}\n")

    print("Please check the outputted file to confirm you are only including samples you want processed")

if __name__ == "__main__":
    main()
