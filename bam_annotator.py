'''
Module for annotating JAFFA results with NB:Z tags from BAM files
'''
import pysam
import pandas as pd

def extract_nb_tags(bam_file):
    """
    Extract read names and NB:Z tags from BAM file
    
    Args:
        bam_file (str): Path to the BAM file
        
    Returns:
        dict: Dictionary mapping read names to NB:Z tag values
    """
    read_to_nb = {}
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                read_name = read.query_name
                # Remove "lib1" from the end of read name
                if read_name.endswith("lib1"):
                    read_name = read_name[:-4]
                
                # Extract NB:Z tag value
                nb_tag = None
                for tag, value in read.get_tags():
                    if tag == "NB":
                        nb_tag = value
                        break
                
                if nb_tag:
                    read_to_nb[read_name] = nb_tag
        
        return read_to_nb
    except Exception as e:
        raise RuntimeError(f"Unable to process BAM file: {e}")

def annotate_jaffa_results(jaffa_file, read_to_nb, output_file=None):
    """
    Annotate JAFFA results CSV with NB:Z tag values
    
    Args:
        jaffa_file (str): Path to the JAFFA results CSV file
        read_to_nb (dict): Dictionary mapping read names to NB:Z tag values
        output_file (str, optional): Path to the output CSV file. If None, returns DataFrame only.
        
    Returns:
        pd.DataFrame: Annotated JAFFA results DataFrame
    """
    try:
        # Read CSV file (with header)
        df = pd.read_csv(jaffa_file)
        
        # Add new column for NB:Z tag values
        df['nb_tag'] = None
        
        # Find matches and fill in NB:Z tag values
        match_count = 0
        for idx, row in df.iterrows():
            contig = str(row['contig'])
            if contig in read_to_nb:
                df.at[idx, 'nb_tag'] = read_to_nb[contig]
                match_count += 1
        
        # Save results if output file provided
        if output_file:
            df.to_csv(output_file, index=False)
        
        return df, match_count
    
    except Exception as e:
        raise RuntimeError(f"Error processing JAFFA file: {e}")