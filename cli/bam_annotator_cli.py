'''
Command-line interface for bam_annotator module
'''
import argparse
import sys
from jaffatools.bam_annotator import extract_nb_tags, annotate_jaffa_results

def main():
    """Command-line entry point for bam_annotator"""
    # Parse arguments
    parser = argparse.ArgumentParser(description='Annotate JAFFA results with NB:Z tags from BAM')
    parser.add_argument('-b', '--bam', required=True, help='input BAM file')
    parser.add_argument('-j', '--jaffa', required=True, help='input JAFFA results CSV file')
    parser.add_argument('-o', '--output', default="jaffa_results_nb.csv", 
                        help='output annotated CSV file (default: jaffa_results_nb.csv)')
    
    # If no arguments provided, use sys.argv
    if len(sys.argv) == 1:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args()
    
    try:
        print("Extracting read names and NB:Z: tags from BAM file...")
        read_to_nb = extract_nb_tags(args.bam)
        print(f"Successfully extracted NB:Z: tags for {len(read_to_nb)} reads from BAM file")
        
        print("Processing JAFFA results CSV file...")
        df, match_count = annotate_jaffa_results(args.jaffa, read_to_nb, args.output)
        print(f"Successfully matched NB:Z: tags for {match_count} contigs")
        print(f"Results saved to {args.output}")
        
        return 0
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    exit(main())