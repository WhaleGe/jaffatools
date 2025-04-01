'''
Command-line interface for fastq_filter module
'''
import argparse
import logging
from multiprocessing import cpu_count
from jaffatools.fastq_filter import extract_read_names, filter_fastq_parallel

def main():
    """Command-line entry point for fastq_filter"""
    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Parse arguments
    parser = argparse.ArgumentParser(description='Filter FASTQ file based on BAM file')
    parser.add_argument('-b', '--bam', required=True, help='input BAM file')
    parser.add_argument('-f', '--fastq', required=True, help='input FASTQ.gz file')
    parser.add_argument('-o', '--output', required=True, help='output FASTQ.gz file')
    parser.add_argument('-c', '--chunk_size', type=int, default=1000000, 
                        help='number of FASTQ records to process in each chunk, default 1,000,000')
    parser.add_argument('-p', '--processes', type=int, default=min(4, cpu_count()),
                        help=f'number of processes to use, default min(4, {cpu_count()})')
    
    args = parser.parse_args()
    
    try:
        # Extract read names from BAM
        read_names = extract_read_names(args.bam)
        
        # Filter FASTQ file and write to output
        filter_fastq_parallel(args.fastq, read_names, args.output, args.chunk_size, args.processes)
        
        logging.info("Processing complete")
        return 0
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    exit(main())