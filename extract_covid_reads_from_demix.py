#!/usr/bin/env python3

"""
Integrated script to:
1. Parse Freyja demix files to identify COVID-positive samples
2. Extract COVID reads from those specific samples 
3. Save reads to tempdir in mutatea folder

Usage:
python extract_covid_reads_from_demix.py \
    --pools_dir /data/service/Pools/EsViritu \
    --tempdir ./covid_reads_temp \
    --min_coverage 0.01
"""

### NTF
# crm: need to adjust to be able to exclude folders without valid poolID (e.g. "old_run")
# crm: could remove min_len

import argparse
import re
import shutil
import tempfile
from pathlib import Path
from typing import List, Tuple, Set, Dict
import pysam


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract COVID reads from demix-positive samples")
    p.add_argument("--pools_dir", required=True, help="Directory containing pool folders")
    p.add_argument("--tempdir", required=True, help="Temporary directory to store COVID reads")
    p.add_argument("--min_coverage", type=float, default=0.01, help="Minimum coverage threshold (default: 0.01)")
    p.add_argument("--min_len", type=int, default=100, help="Minimum read length (default: 100)")
    p.add_argument("--covid_ref", default="NC_045512.2", help="COVID reference accession (default: NC_045512.2)")
    p.add_argument("--dry_run", action="store_true", help="Only identify samples, don't extract reads")
    return p.parse_args()


def parse_demix_file(demix_path: Path, min_coverage: float) -> Tuple[bool, float]:
    """Parse a single demix.out file and determine if sample is COVID-positive."""
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
    """Extract sample ID from demix file path"""
    filename = demix_path.name
    match = re.match(r'([A-Z0-9]+)\.p\d+\.demix\.out', filename)
    return match.group(1) if match else filename.replace('.demix.out', '')


def find_all_demix_files(pools_dir: Path) -> List[Path]:
    """Find all demix.out files in the pools directory structure"""
    demix_files = []
    pool_count = 0
    
    for pool_dir in pools_dir.iterdir():
        if not pool_dir.is_dir() or not pool_dir.name.startswith('p'):
            continue
        
        pool_count += 1
        if pool_count % 10 == 0:
            print(f"  Scanned {pool_count} pools...")
            
        for freyja_dir in pool_dir.rglob('*freyja*'):
            if freyja_dir.is_dir():
                demix_files.extend(freyja_dir.glob('*.demix.out'))
    
    return demix_files


def find_sample_fastq_files(pools_dir: Path, pool_id: str, sample_id: str) -> Tuple[Path, Path]:
    """Find FASTQ files for a specific sample using SampleListRaw"""
    samplelist_file = pools_dir / pool_id / "SampleListRaw"
    
    if not samplelist_file.exists():
        print(f"SampleListRaw file not found for pool {pool_id}")
        return None, None
    
    try:
        with open(samplelist_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 3 and parts[0] == sample_id:
                    r1_path = Path(parts[1])
                    r2_path = Path(parts[2])
                    
                    # Check if files exist
                    if r1_path.exists() and r2_path.exists():
                        return r1_path, r2_path
                    else:
                        print(f"FASTQ files not found for {sample_id}: {r1_path}, {r2_path}")
                        return None, None
        
        print(f"Sample {sample_id} not found in SampleListRaw for pool {pool_id}")
        return None, None
        
    except Exception as e:
        print(f"Error reading SampleListRaw for pool {pool_id}: {e}")
        return None, None


def extract_covid_reads_from_bam(bam_path: Path, covid_ref: str, min_len: int, output_file: Path) -> int:
    """Extract COVID reads from a BAM file and save as FASTQ"""
    try:
        # Check if BAM is indexed
        bai_file = bam_path.with_name(bam_path.name + ".bai")
        if not bai_file.exists():
            print(f"Indexing {bam_path.name}...")
            pysam.index(str(bam_path))
        
        bam = pysam.AlignmentFile(bam_path, "rb")
        
        # Check if COVID reference exists in BAM
        if covid_ref not in bam.references:
            print(f"Warning: COVID reference {covid_ref} not found in {bam_path.name}")
            bam.close()
            return 0
        
        read_count = 0
        with open(output_file, 'w') as out_handle:
            for read in bam.fetch(covid_ref):
                if read.is_unmapped or read.query_alignment_length < min_len:
                    continue
                
                if read.query_sequence is None or read.qual is None:
                    continue
                
                # Write FASTQ format
                sample_id = re.search(r"virus_pathogen_database\.mmi\.(.*?)\.p\d+", bam_path.name)
                sample_name = sample_id.group(1) if sample_id else "UNKNOWN"
                
                out_handle.write(f"@{sample_name}:{read.query_name}\n")
                out_handle.write(f"{read.query_sequence}\n+\n{read.qual}\n")
                read_count += 1
        
        bam.close()
        return read_count
        
    except Exception as e:
        print(f"Error extracting reads from {bam_path.name}: {e}")
        return 0


def main():
    args = parse_args()
    
    pools_dir = Path(args.pools_dir)
    tempdir = Path(args.tempdir)
    
    # Create tempdir structure immediately
    tempdir.mkdir(parents=True, exist_ok=True)
    (tempdir / "logs").mkdir(exist_ok=True)
    (tempdir / "reads").mkdir(exist_ok=True)
    
    print(f"Using temporary directory: {tempdir}")
    
    # Find all demix files
    demix_files = find_all_demix_files(pools_dir)
    print(f"Found {len(demix_files)} demix files")
    
    # Identify COVID-positive samples
    covid_positive_samples = []
    total_samples = 0
    
    print("Analyzing demix files for COVID-positive samples...")
    for i, demix_path in enumerate(demix_files, 1):
        if i % 500 == 0:
            print(f"  Processed {i}/{len(demix_files)} demix files...")
        
        total_samples += 1
        is_positive, coverage = parse_demix_file(demix_path, args.min_coverage)
        
        if is_positive:
            sample_id = extract_sample_id_from_demix(demix_path)
            pool_id = demix_path.parent.parent.name
            
            covid_positive_samples.append({
                'sample_id': sample_id,
                'pool_id': pool_id,
                'coverage': coverage,
                'demix_file': str(demix_path)
            })
    
    print(f"Found {len(covid_positive_samples)} COVID-positive samples out of {total_samples} total")
    
    # Save list of positive samples
    samples_file = tempdir / "covid_positive_samples.txt"
    with open(samples_file, 'w') as f:
        f.write("# COVID-positive samples identified from Freyja demix results\n")
        f.write("# Format: sample_id pool_id coverage demix_file\n")
        
        for sample in sorted(covid_positive_samples, key=lambda x: (x['pool_id'], x['sample_id'])):
            f.write(f"{sample['sample_id']} {sample['pool_id']} {sample['coverage']:.4f} {sample['demix_file']}\n")
    
    if args.dry_run:
        print("Dry run mode - not extracting reads")
        print(f"\nCOVID-positive samples identified:")
        print("=" * 60)
        for sample in sorted(covid_positive_samples, key=lambda x: (x['pool_id'], x['sample_id'])):
            print(f"Pool {sample['pool_id']}: Sample {sample['sample_id']} (coverage: {sample['coverage']:.4f})")
        print(f"\nCOVID-positive samples saved to: {samples_file}")
        return
    
    # Extract reads from COVID-positive samples
    reads_dir = tempdir / "reads"
    extraction_log = tempdir / "logs" / "extraction_log.txt"
    
    total_reads_extracted = 0
    successful_extractions = 0
    
    with open(extraction_log, 'w') as log:
        log.write("COVID Read Extraction Log\n")
        log.write("=" * 50 + "\n\n")
        
        for i, sample in enumerate(covid_positive_samples, 1):
            sample_id = sample['sample_id']
            pool_id = sample['pool_id']
            
            print(f"Processing sample {i}/{len(covid_positive_samples)}: {sample_id} (pool {pool_id})")
            
            # Find BAM files for this sample
            bam_files = find_bam_files_for_sample(pools_dir, pool_id, sample_id)
            
            if not bam_files:
                log.write(f"No BAM files found for sample {sample_id} in pool {pool_id}\n")
                print(f"  Warning: No BAM files found for {sample_id}")
                continue
            
            for bam_file in bam_files:
                output_file = reads_dir / f"{pool_id}.{sample_id}.covid.fastq"
                
                print(f"  Extracting reads from {bam_file.name}")
                read_count = extract_covid_reads_from_bam(bam_file, args.covid_ref, args.min_len, output_file)
                
                if read_count > 0:
                    total_reads_extracted += read_count
                    successful_extractions += 1
                    log.write(f"SUCCESS: {sample_id} - {read_count} reads extracted from {bam_file.name}\n")
                    print(f"    Extracted {read_count} COVID reads")
                else:
                    log.write(f"NO READS: {sample_id} - no COVID reads found in {bam_file.name}\n")
                    print(f"    No COVID reads found")
    
    # Summary
    print(f"\nExtraction Summary:")
    print(f"- COVID-positive samples: {len(covid_positive_samples)}")
    print(f"- Successful extractions: {successful_extractions}")
    print(f"- Total COVID reads extracted: {total_reads_extracted:,}")
    print(f"- Reads saved to: {reads_dir}")
    print(f"- Logs saved to: {extraction_log}")
    
    # Create summary file
    summary_file = tempdir / "extraction_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"COVID Read Extraction Summary\n")
        f.write(f"=" * 40 + "\n\n")
        f.write(f"Pools directory: {pools_dir}\n")
        f.write(f"Output directory: {tempdir}\n")
        f.write(f"Minimum coverage: {args.min_coverage}\n")
        f.write(f"Minimum read length: {args.min_len}\n")
        f.write(f"COVID reference: {args.covid_ref}\n\n")
        f.write(f"Results:\n")
        f.write(f"- Total samples analyzed: {total_samples}\n")
        f.write(f"- COVID-positive samples: {len(covid_positive_samples)}\n")
        f.write(f"- Successful extractions: {successful_extractions}\n")
        f.write(f"- Total COVID reads extracted: {total_reads_extracted:,}\n")


if __name__ == "__main__":
    main()
