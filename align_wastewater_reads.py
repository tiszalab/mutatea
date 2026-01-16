from multiprocessing import Pool

# defining helper function
def _align_wastewater_reads(args):

    # set the arguments it will be expecting
    read_file, pool_id, reference_fasta, pool_output_dir, minimap2_path, samtools_path, threads = args
    
# creating subprocess lines for both single and paired reads
    # paired reads
    if isinstance(read_file, tuple):
        r1_file, r2_file = read_file
        sample_name = os.path.basename(r1_file).split(".")[0]
        output_bam = os.path.join(pool_output_dir, f"{sample_name}.{pool_id}.sort.bam")
        cmd = f"{minimap2_path} -t {threads} -ax sr {reference_fasta} {r1_file} {r2_file} | {samtools_path} view -@ {threads} -bS | {samtools_path} sort -@ {threads} -o {output_bam}"
    # single reads
    else:
        filename = os.path.basename(read_file)
        parts = filename.split(".")
        sample_name = parts[1] if len(parts) >= 3 else parts[0]
        output_bam = os.path.join(pool_output_dir, f"{sample_name}.{pool_id}.sort.bam")
        cmd = f"{minimap2_path} -t {threads} -ax map-ont {reference_fasta} {read_file} | {samtools_path} view -@ {threads} -bS | {samtools_path} sort -@ {threads} -o {output_bam}"
    
# will run the subprocess lines for minimap2 and samtools view/sort/index
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        subprocess.run([samtools_path, "index", output_bam], check=True, capture_output=True)
        return f"Success: {sample_name}"
    except subprocess.CalledProcessError as e:
        return f"Error processing {sample_name}: {e}"

def align_wastewater_reads(reads_by_pool: dict, reference_dir: str, pools: str, threads: int = 4, max_workers: int = 2) -> dict:
    # find paths of imported functions
    minimap2_path = shutil.which("minimap2")
    samtools_path = shutil.which("samtools")
    
    if not minimap2_path or not samtools_path:
        raise RuntimeError("minimap2 or samtools not found in PATH. Please ensure they are installed and accessible.")
    
    # load in reference fasta
    reference_fasta = glob.glob(os.path.join(reference_dir, "*.fna"))[0]
    
    # Prepare all tasks
    # create empty list
    tasks = []
    for pool_id, read_files in reads_by_pool.items():
        pool_output_dir = os.path.join(pools, pool_id)
        os.makedirs(pool_output_dir, exist_ok=True)
        
        # removing enumerate line
        for read_file in read_files:
            # for all reads, append the arguments to the tasks
            tasks.append((read_file, pool_id, reference_fasta, pool_output_dir, 
                         minimap2_path, samtools_path, threads))
    
    # Process in parallel
    # print line is now saying number of tasks run with number of max_workers, not number of reads/total per pool
    print(f"\nAligning {len(tasks)} samples using {max_workers} parallel workers\n")
    
    # run 
    with Pool(processes=max_workers) as pool:
        results = pool.map(_align_single_sample, tasks)
    
    # Print results
    for result in results:
        if result.startswith("Error"):
            print(result)
    
    print()