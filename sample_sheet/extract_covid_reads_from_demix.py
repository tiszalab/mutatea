#!/usr/bin/env python3

# crm: outdated, we are not running the files from Esviritu, the covid reads are already pulled out as paired reads in another directory

# parse Freyja demix files to identify COVID-positive samples
# find BAM files from those specific samples 
# save file paths of the BAM files in a txt file

"""
Executable code:
python extract_covid_reads_from_demix.py \
    --pools_dir /data/service/Pools/EsViritu \
    --outputfile ./covid_bam_files.txt
"""

import argparse
import re
import shutil
import tempfile
from pathlib import Path
from typing import List, Tuple, Set, Dict
import pysam


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Finding file paths of pre-aligned COVID bam files from demix samples")
    p.add_argument("--pools_dir", required=True, help="Directory containing pool folders")
    p.add_argument("--outputfile", required=True, help="Output file path for COVID-positive samples list")    
    p.add_argument("--covid_ref", default="NC_045512.2", help="COVID reference accession (default: NC_045512.2)")
    return p.parse_args()


# determine if each demix.out file found covid-positive samples
def parse_demix_file(demix_path: Path) -> bool:
    try:
        with open(demix_path, 'r') as f:
            content = f.read()
        
        summarized_match = re.search(r'summarized\s+\[(.*?)\]', content)
        if not summarized_match:
            return False
        
        summarized_str = summarized_match.group(1)
        is_positive = bool(summarized_str.strip() and summarized_str != '[]')
        return is_positive
        
    except Exception as e:
        print(f"Error parsing {demix_path}: {e}")
        return False


# extract the sampleID from the file path of the demix.out file names
def extract_sample_id_from_demix(demix_path: Path) -> str:
    filename = demix_path.name
    match = re.match(r'([A-Z0-9]+)\.p\d+\.demix\.out', filename)
    return match.group(1) if match else filename.replace('.demix.out', '')


# find all the demix.out files in the given pools directory
def find_all_demix_files(pools_dir: Path) -> List[Path]:
    demix_files = []
    pool_count = 0
    
    for pool_dir in pools_dir.iterdir():
        if not pool_dir.is_dir() or not pool_dir.name.startswith('p'):
            continue
        
        pool_count += 1
        if pool_count % 10 == 0:
            print(f"  Scanned {pool_count} pools")
            
        for freyja_dir in pool_dir.rglob('*freyja*'):
            if freyja_dir.is_dir():
                demix_files.extend(freyja_dir.glob('*.demix.out'))
    return demix_files


def main():
    args = parse_args()
    
    pools_dir = Path(args.pools_dir)
    outputfile = Path(args.outputfile)
    
    # find all demix files
    demix_files = find_all_demix_files(pools_dir)
    print(f"Found {len(demix_files)} demix files")
    
    # identify COVID-positive samples
    covid_positive_samples = []
    total_samples = 0
    
    print("Analyzing demix files for COVID-positive samples")
    # goes through each demix_file
    for i, demix_path in enumerate(demix_files, 1):
        # print line every 500 demix files to show progress
        if i % 500 == 0:
            print(f"  Processed {i}/{len(demix_files)} demix files")
        
        # add to the counter with every loop
        total_samples += 1
        # calls function in to analyze each demix file, returns if the demix file contains covid reads (positive)
        is_positive = parse_demix_file(demix_path)
        
        # if the demix file reported covid
        if is_positive:
            sample_id = extract_sample_id_from_demix(demix_path)
            pool_id = demix_path.parent.parent.name
            
            # swap demix.out with sort.bam to get corresponding BAM file path
            bam_file = str(demix_path).replace('.demix.out', '.sort.bam')
            
            # return the following values for each demix file
            covid_positive_samples.append({
                'sample_id': sample_id,
                'pool_id': pool_id,
                'bam_file': bam_file
            })
    
    print(f"Found COVID reads in {len(covid_positive_samples)} samples out of total {total_samples} samples")
    
    # check that the bam actually exists
    valid_samples = []
    missing_bam_files = []
    
    for sample in covid_positive_samples:
        bam_path = Path(sample['bam_file'])
        if bam_path.exists():
            valid_samples.append(sample)
        else:
            missing_bam_files.append(sample['bam_file'])
    
    print(f"Verified BAM files: {len(valid_samples)} valid, {len(missing_bam_files)} missing")
    
    # print all BAM files that were given but did not exist for downstream troubleshooting
    if missing_bam_files:
        print("Missing BAM files:")
        for bam_file in missing_bam_files:
            print(f"  {bam_file}")

    # prompt the user to check that their file before plugging it into mutatea 
    print("Please check the outputted file to confirm you are only including samples you want processed")
    
    # save list of covid samples as txt file
    samples_file = outputfile
    with open(samples_file, 'w') as f:
        f.write("# BAM files identified from Freyja demix results\n")
        f.write("# sample_id\tsample_type\tfile_path_1\tfile_path_2\n")
        
        for sample in sorted(valid_samples, key=lambda x: (x['pool_id'], x['sample_id'])):
            sample_id = sample['sample_id']
            sample_type = "bam"
            file_path_1 = sample['bam_file']
            file_path_2 = ""  # empty on purpose for bam files
            
            f.write(f"{sample_id}\t{sample_type}\t{file_path_1}\t{file_path_2}\n")
        return

if __name__ == "__main__":
    main()
