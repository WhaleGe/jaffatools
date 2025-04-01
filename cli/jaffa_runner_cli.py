'''
Command-line interface for jaffa_runner module
'''
import argparse
import sys
from jaffatools.jaffa_runner import run_jaffa

def main():
    """Command-line entry point for jaffa_runner"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run JAFFA fusion gene detection")
    
    # Original parameters
    parser.add_argument("-m", "--memory", default="72GB", 
                        help="Memory allocation (default: 72GB)")
    parser.add_argument("-n", "--threads", type=int, default=12, 
                        help="Number of threads (default: 12)")
    parser.add_argument("-f", "--sample-name", default="sample", 
                        help="Output file name of report")
    parser.add_argument("-j", "--jaffa-path", 
                        default="/path/to/your/JAFFA-version-2.3", 
                        help="Path to JAFFA installation")
    parser.add_argument("-i", "--input", 
                        default="/path/to/your/sample.fastq.gz", 
                        help="Input FASTQ file")
    
    # Additional bpipe parameters
    parser.add_argument("-b", "--branch", 
                        help="Comma separated list of branches to limit execution to")
    parser.add_argument("-d", "--dir", 
                        help="Output directory")
    parser.add_argument("-l", "--resource", action="append", 
                        help="Place limit on named resource (format: resource=value)")
    parser.add_argument("-L", "--interval", 
                        help="The default genomic interval to execute pipeline for (samtools format)")
    parser.add_argument("-p", "--param", action="append", 
                        help="Defines a pipeline parameter, or file of parameters via @<file> (format: param=value)")
    parser.add_argument("-r", "--report", action="store_true", 
                        help="Generate an HTML report / documentation for pipeline")
    parser.add_argument("-R", "--report-template", 
                        help="Generate report using named template")
    parser.add_argument("-t", "--test", action="store_true", 
                        help="Test mode")
    parser.add_argument("-u", "--until", 
                        help="Run until stage given")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Print internal logging to standard error")
    parser.add_argument("-y", "--yes", action="store_true", 
                        help="Answer yes to any prompts or questions")
    
    args = parser.parse_args()
    
    # Run JAFFA with specified parameters
    success = run_jaffa(
        memory=args.memory,
        threads=args.threads,
        sample_name=args.sample_name,
        jaffa_path=args.jaffa_path,
        fastq_file=args.input,
        branch=args.branch,
        output_dir=args.dir,
        report=args.report,
        report_template=args.report_template,
        resources=args.resource,
        interval=args.interval,
        parameters=args.param,
        test=args.test,
        until=args.until,
        verbose=args.verbose,
        yes=args.yes
    )
    
    # Exit with appropriate status code
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())