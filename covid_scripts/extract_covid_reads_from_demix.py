#!/usr/bin/env python3

# crm: outdated, we are not running the files from Esviritu, the covid reads are already pulled out as paired reads in another directory

# parse Freyja demix files to identify COVID-positive samples
# find COVID-positive BAM files from those specific samples 
# save file paths of the BAM files in a txt file

"""
Executable code:
python extract_covid_reads_from_demix.py \
    --pools_dir /data/service/Pools/EsViritu \
    --outputfile ./covid_bam_files.txt \
    --min_coverage 0.01
"""

### NTF
# crm: need to prompt user to check and exclude folders without valid poolID (e.g. "old_run")
# crm: add check for positive demix.out having associated sort.bam files
##### crm: need to check that each bam file you're listing is real 

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
    p.add_argument("--min_coverage", type=float, default=0.01, help="Minimum coverage threshold (default: 0.01)")
    p.add_argument("--covid_ref", default="NC_045512.2", help="COVID reference accession (default: NC_045512.2)")
    return p.parse_args()


def parse_demix_file(demix_path: Path, min_coverage: float) -> Tuple[bool, float]:
    """parse a single demix.out file and determine if sample is COVID-positive."""
    try:
        with open(demix_path, 'r') as f:
            content = f.read()
        
        coverage_match = re.search(r'coverage\s+([\d.]+)', content)
        if not coverage_match:
            return False, 0.0
        
        coverage = float(coverage_match.group(1))
        
        summarized_match = re.search(r'summarized\s+\[(.*?)\]', content)
        if not summarized_match:
            return False, coverage
        
        summarized_str = summarized_match.group(1)
        is_positive = bool(summarized_str.strip() and summarized_str != '[]')
        
        return is_positive and coverage >= min_coverage, coverage
        
    except Exception as e:
        print(f"Error parsing {demix_path}: {e}")
        return False, 0.0


def extract_sample_id_from_demix(demix_path: Path) -> str:
    """extract sample ID from demix file path"""
    filename = demix_path.name
    match = re.match(r'([A-Z0-9]+)\.p\d+\.demix\.out', filename)
    return match.group(1) if match else filename.replace('.demix.out', '')


def find_all_demix_files(pools_dir: Path) -> List[Path]:
    """find all demix.out files in the pools directory structure"""
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
        # calls function in to analyze each demix file, returns if the demix file contains covid reads (positive) and the coverage of the sample
        is_positive, coverage = parse_demix_file(demix_path, args.min_coverage)
        
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
                'coverage': coverage,
                'bam_file': bam_file
            })
    
    print(f"Found COVID reads in {len(covid_positive_samples)} samples out of total {total_samples} samples")
    
    # crm: add a check to make sure that you can get the bam files from the given demix.out files
    # save list of covid samples as txt file
    samples_file = outputfile
    with open(samples_file, 'w') as f:
        f.write("# COVID-positive samples identified from Freyja demix results\n")
        f.write("# Format: sample_id pool_id coverage bam_file\n")
        
        for sample in sorted(covid_positive_samples, key=lambda x: (x['pool_id'], x['sample_id'])):
            f.write(f"{sample['sample_id']} {sample['pool_id']} {sample['coverage']:.4f} {sample['bam_file']}\n")
        return


if __name__ == "__main__":
    main()
