'''
Module for running JAFFA fusion gene detection pipeline
'''
import subprocess
import os
import time

def run_jaffa(memory="72GB", threads=12, sample_name="sample", 
              jaffa_path="/path/to/your/JAFFA-version-2.3", 
              fastq_file="/path/to/your/sample.fastq.gz",
              branch=None, output_dir=None, report=False, report_template=None, 
              resources=None, interval=None, parameters=None, test=False, 
              until=None, verbose=False, yes=False):
    """
    Run JAFFA fusion gene detection pipeline
    
    Args:
        memory (str): Memory allocation
        threads (int): Number of threads
        sample_name (str): Sample name for output files
        jaffa_path (str): Path to JAFFA installation
        fastq_file (str): Path to input FASTQ file
        branch (str, optional): Comma separated list of branches to limit execution to
        output_dir (str, optional): Output directory
        report (bool, optional): Generate an HTML report
        report_template (str, optional): Generate report using named template
        resources (list, optional): Place limit on named resources
        interval (str, optional): Default genomic interval to execute pipeline for
        parameters (list, optional): Pipeline parameters
        test (bool, optional): Test mode
        until (str, optional): Run until stage given
        verbose (bool, optional): Print internal logging to standard error
        yes (bool, optional): Answer yes to any prompts or questions
        
    Returns:
        bool: True if JAFFA completed successfully, False otherwise
    """
    # Check if input file exists
    if not os.path.exists(fastq_file):
        raise FileNotFoundError(f"Input FASTQ file '{fastq_file}' does not exist")
    
    # Construct the command
    bpipe_cmd = os.path.join(jaffa_path, "tools/bin/bpipe")
    jaffal_groovy = os.path.join(jaffa_path, "JAFFAL.groovy")
    
    # Start with base command
    command = [bpipe_cmd, "run"]
    
    # Add options
    if memory:
        command.extend(["-m", memory])
    
    if threads:
        command.extend(["-n", str(threads)])
    
    if sample_name:
        command.extend(["-f", sample_name])
    
    if branch:
        command.extend(["-b", branch])
    
    if output_dir:
        command.extend(["-d", output_dir])
    
    if resources:
        for resource in resources:
            command.extend(["-l", resource])
    
    if interval:
        command.extend(["-L", interval])
    
    if parameters:
        for param in parameters:
            command.extend(["-p", param])
    
    if report:
        command.append("-r")
    
    if report_template:
        command.extend(["-R", report_template])
    
    if test:
        command.append("-t")
    
    if until:
        command.extend(["-u", until])
    
    if verbose:
        command.append("-v")
    
    if yes:
        command.append("-y")
    
    # Add script and input file
    command.append(jaffal_groovy)
    command.append(fastq_file)
    
    # Record start time
    start_time = time.time()
    
    try:
        # Execute the command
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Collect output
        stdout_lines = []
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                stdout_lines.append(output.strip())
                print(output.strip())  # Print real-time output
        
        # Get return code
        return_code = process.poll()
        
        # Print any errors
        if return_code != 0:
            stderr = process.stderr.read()
            print(f"Error occurred (return code {return_code}):")
            print(stderr)
            return False
        
        # Calculate elapsed time
        elapsed_time = time.time() - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        print(f"JAFFA analysis completed successfully in {int(hours)}h {int(minutes)}m {int(seconds)}s")
        return True
        
    except Exception as e:
        print(f"Error running JAFFA: {str(e)}")
        return False