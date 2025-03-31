'''
JaffaTools - A Python package for working with JAFFA fusion gene detection
'''

__version__ = '0.1.0'
__author__ = 'Huang Waidong'
__email__ = 'wdhuang927@gmail.com'

# Import main functions for easy access
from .fastq_filter import extract_read_names, filter_fastq_parallel
from .bam_annotator import extract_nb_tags, annotate_jaffa_results
from .jaffa_runner import run_jaffa

# Define what's available when using `from jaffatools import *`
__all__ = [
    'extract_read_names',
    'filter_fastq_parallel',
    'extract_nb_tags',
    'annotate_jaffa_results',
    'run_jaffa'
]

