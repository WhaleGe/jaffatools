'''
Module for filtering FASTQ files based on read names in BAM files
'''
import gzip
import os
import pysam
import logging
from multiprocessing import Pool, cpu_count

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extract_read_names(bam_file):
    """
    Extract read names from BAM file, removing 'lib1' suffix
    
    Args:
        bam_file (str): Path to the BAM file
        
    Returns:
        set: Set of read names from the BAM file
    """
    logging.info(f"Extracting read names from BAM file {bam_file}")
    read_names = set()
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                # Remove 'lib1' suffix from read name
                name = read.query_name
                if name.endswith('lib1'):
                    name = name[:-4]  # remove 'lib1'
                read_names.add(name)
        logging.info(f"Extracted {len(read_names)} unique read names from BAM file")
        return read_names
    except Exception as e:
        logging.error(f"Error extracting read names: {str(e)}")
        raise

def count_fastq_records(fastq_file):
    """
    Count the number of records in FASTQ file
    
    Args:
        fastq_file (str): Path to the FASTQ.gz file
        
    Returns:
        int: Number of records in the FASTQ file
    """
    logging.info(f"Counting records in FASTQ file {fastq_file}")
    count = 0
    with gzip.open(fastq_file, 'rt') as f:
        while True:
            # Each record has four lines
            header = f.readline()
            if not header:
                break
            f.readline()  # sequence line
            f.readline()  # '+' line
            f.readline()  # quality score line
            count += 1
    logging.info(f"FASTQ file contains {count} records")
    return count

def process_fastq_chunk(args):
    """
    Process a chunk of FASTQ file
    
    Args:
        args (tuple): (chunk_id, fastq_file, start_pos, chunk_size, read_names)
        
    Returns:
        list: List of matching FASTQ records as (header, seq, plus, qual) tuples
    """
    chunk, fastq_file, start_pos, chunk_size, read_names = args
    output_records = []
    
    with gzip.open(fastq_file, 'rt') as f:
        # Jump to starting position
        f.seek(start_pos)
        
        for _ in range(chunk_size):
            # Read four lines (one FASTQ record)
            header = f.readline().strip()
            if not header:
                break  # End of file
                
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            if not header or not seq or not plus or not qual:
                break  # Incomplete record
                
            # Extract read name (remove '@' prefix and content after space)
            read_name = header[1:].split()[0]
            
            # Check if read name is in BAM file
            if read_name in read_names:
                output_records.append((header, seq, plus, qual))
    
    return output_records

def filter_fastq_parallel(fastq_file, read_names, output_file, chunk_size=1000000, processes=None):
    """
    Filter FASTQ file in parallel, keeping only reads that are in the provided read_names set
    
    Args:
        fastq_file (str): Path to the input FASTQ.gz file
        read_names (set): Set of read names to keep
        output_file (str): Path to the output FASTQ.gz file
        chunk_size (int): Number of FASTQ records to process in each chunk
        processes (int): Number of parallel processes to use
        
    Returns:
        int: Number of matching records written to output file
    """
    if processes is None:
        processes = min(4, cpu_count())
        
    logging.info(f"Filtering FASTQ file using {processes} processes")
    
    # Find file size
    fastq_size = os.path.getsize(fastq_file)
    
    # Create list of task parameters
    tasks = []
    with gzip.open(fastq_file, 'rt') as f:
        pos = 0
        chunk = 0
        while True:
            start_pos = f.tell()
            
            # Read four lines to check if end of file
            header = f.readline()
            if not header:
                break
            f.readline()  # sequence
            f.readline()  # '+'
            f.readline()  # quality
            
            # Skip chunk_size-1 records
            for _ in range(chunk_size - 1):
                for _ in range(4):  # 4 lines per record
                    if not f.readline():
                        break
            
            tasks.append((chunk, fastq_file, start_pos, chunk_size, read_names))
            chunk += 1
            
            # Record progress
            if chunk % 10 == 0:
                current_pos = f.tell()
                progress = (current_pos / fastq_size) * 100
                logging.info(f"Preparation progress: {progress:.2f}%")
    
    # Process in parallel using process pool
    logging.info(f"Starting parallel processing of {len(tasks)} chunks")
    with Pool(processes=processes) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(process_fastq_chunk, tasks)):
            results.extend(result)
            if (i + 1) % 10 == 0 or (i + 1) == len(tasks):
                progress = ((i + 1) / len(tasks)) * 100
                logging.info(f"Processing progress: {progress:.2f}%, found {len(results)} matching records")
    
    # Write results to output file
    logging.info(f"Writing {len(results)} matching records to {output_file}")
    with gzip.open(output_file, 'wt') as out:
        for header, seq, plus, qual in results:
            out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
            
    return len(results)