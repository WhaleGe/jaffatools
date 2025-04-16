'''
Command-line interface for JAFFA visualization using R
'''
import argparse
import sys
import os
import subprocess

def r_plot_fusion_circos(csv_file, genome_version='hg38', output_file=None):
    """
    Generate Circos plots to visualize fusion genes using R
    
    Args:
        csv_file (str): Path to JAFFA result CSV file
        genome_version (str): Genome version ('hg19', 'hg38', 'mm10')
        output_file (str, optional): Path to output image file
    
    Returns:
        bool: Returns True if successful, False if failed
    """
    try:
        import subprocess
        import os
        
        # Get R script path
        r_script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'jaffa_plots.R')
        
        # Get data directory path
        data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
        
        # Build command
        cmd = ['Rscript', r_script_path, 'circos', csv_file, genome_version]
        if output_file:
            cmd.append(output_file)
        # Add data directory path
        cmd.append(data_dir)
        
        # Execute R script
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode != 0:
            print(f"R script execution error: {process.stderr}")
            return False
        
        print(process.stdout)
        return True
    
    except Exception as e:
        print(f"Error when executing R script: {str(e)}")
        return False

def r_plot_fusion_per_cell(csv_file, output_file=None):
    """
    Generate statistical plots of fusion gene counts per cell using R
    
    Args:
        csv_file (str): Path to JAFFA result CSV file
        output_file (str, optional): Path to output image file
    
    Returns:
        bool: Returns True if successful, False if failed
    """
    try:
        # Get R script path
        r_script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'jaffa_plots.R')
        
        # Build command
        cmd = ['Rscript', r_script_path, 'cell_stats', csv_file]
        if output_file:
            cmd.append(output_file)
        
        # Execute R script
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode != 0:
            print(f"R script execution error: {process.stderr}", file=sys.stderr)
            return False
        
        print(process.stdout)
        return True
    
    except Exception as e:
        print(f"Error when executing R script: {str(e)}", file=sys.stderr)
        return False

def r_plot_fusion_classification(csv_file, output_file=None):
    """
    Generate statistical plots of fusion gene classifications using R
    
    Args:
        csv_file (str): Path to JAFFA result CSV file
        output_file (str, optional): Path to output image file
    
    Returns:
        bool: Returns True if successful, False if failed
    """
    try:
        # Get R script path
        r_script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'jaffa_plots.R')
        
        # Build command
        cmd = ['Rscript', r_script_path, 'class_stats', csv_file]
        if output_file:
            cmd.append(output_file)
        
        # Execute R script
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode != 0:
            print(f"R script execution error: {process.stderr}", file=sys.stderr)
            return False
        
        print(process.stdout)
        return True
    
    except Exception as e:
        print(f"Error when executing R script: {str(e)}", file=sys.stderr)
        return False

def check_r_installation():
    """
    Check if R and required R packages are installed
    
    Returns:
        bool: Returns True if R is installed and can execute simple commands, otherwise False
    """
    try:
        # Try to execute R version check
        process = subprocess.run(['Rscript', '--version'], capture_output=True, text=True)
        if process.returncode != 0:
            print("R installation not detected. Please install R and ensure 'Rscript' is executable from the command line.", file=sys.stderr)
            return False
        
        # Check for required R packages
        r_code = 'if(!all(c("data.table", "ggplot2", "RCircos", "gtools", "dplyr") %in% installed.packages()[,"Package"])) { cat("Missing required R packages"); quit(status=1) }'
        process = subprocess.run(['Rscript', '-e', r_code], capture_output=True, text=True)
        
        if process.returncode != 0:
            print("Missing required R packages. Please run the following command in R to install necessary packages:", file=sys.stderr)
            print("install.packages(c('data.table', 'ggplot2', 'RCircos', 'gtools', 'dplyr'))", file=sys.stderr)
            return False
        
        return True
    
    except FileNotFoundError:
        print("R installation not detected. Please install R and ensure 'Rscript' is executable from the command line.", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error checking R installation: {str(e)}", file=sys.stderr)
        return False

def main():
    """Command-line entry point for jaffa_viz"""
    # Check if R is installed
    if not check_r_installation():
        print("Please install R and required R packages before running this tool.", file=sys.stderr)
        return 1
    
    parser = argparse.ArgumentParser(description='Visualize JAFFA fusion gene detection results using R')
    
    # Subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Circos plot command
    circos_parser = subparsers.add_parser('circos', help='Draw Circos plots of fusion genes')
    circos_parser.add_argument('-i', '--input', required=True, help='JAFFA result CSV file')
    circos_parser.add_argument('-g', '--genome', choices=['hg19', 'hg38', 'mm10'], default='hg38',help='Genome version (default: hg38)')
    circos_parser.add_argument('-o', '--output', required=True, help='Output image file (.pdf)')
    
    # Cell statistics command
    cell_stats_parser = subparsers.add_parser('cell-stats', help='Draw statistical plots of fusion gene counts per cell')
    cell_stats_parser.add_argument('-i', '--input', required=True, help='JAFFA result CSV file')
    cell_stats_parser.add_argument('-o', '--output', required=True, help='Output image file (.pdf/.png)')
    
    # Classification statistics command
    class_stats_parser = subparsers.add_parser('class-stats', help='Draw statistical plots of fusion gene classifications')
    class_stats_parser.add_argument('-i', '--input', required=True, help='JAFFA result CSV file')
    class_stats_parser.add_argument('-o', '--output', required=True, help='Output image file (.pdf/.png)')
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    try:
        if args.command == 'circos':
            # Verify file exists
            if not os.path.exists(args.input):
                print(f"Error: Input file '{args.input}' does not exist", file=sys.stderr)
                return 1
            
            # Draw Circos plot
            success = r_plot_fusion_circos(args.input, args.genome, args.output)
            if success:
                print(f"Circos plot saved to {args.output}")
            else:
                print("Failed to generate Circos plot", file=sys.stderr)
                return 1
        
        elif args.command == 'cell-stats':
            # Verify file exists
            if not os.path.exists(args.input):
                print(f"Error: Input file '{args.input}' does not exist", file=sys.stderr)
                return 1
            
            # Draw cell fusion gene count statistics
            success = r_plot_fusion_per_cell(args.input, args.output)
            if success:
                print(f"Cell fusion gene count statistics plot saved to {args.output}")
            else:
                print("Failed to generate cell fusion gene count statistics plot", file=sys.stderr)
                return 1
        
        elif args.command == 'class-stats':
            # Verify file exists
            if not os.path.exists(args.input):
                print(f"Error: Input file '{args.input}' does not exist", file=sys.stderr)
                return 1
            
            # Draw fusion gene classification statistics
            success = r_plot_fusion_classification(args.input, args.output)
            if success:
                print(f"Fusion gene classification statistics plot saved to {args.output}")
            else:
                print("Failed to generate fusion gene classification statistics plot", file=sys.stderr)
                return 1
        
        return 0
    
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
