
import os
from os.path import splitext

from helpers.fastq import (
    read_fastq,
    write_fastq,
    check_gc,
    check_length,
    check_quality
)

from helpers.bio_files import (
    read_fasta,
    get_best_blast_hits
)



def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | int | float = (0, 100),
    length_bounds: tuple | int | float = (0, 2**32),
    quality_threshold: float = 0
) -> None:
    
    '''
    Filters fastq file - checks GC, read length and avg quality
    
    
    params:

    input_fastq - path to the init file
    
    output_fastq - name of the output file
    
    gc_bounds - allowed interval of GC (%)
    
    length_bounds - allowed interval od read leangth
    
    quality_threshold - min average quality (Q)

    The resut is being written into the folder "filtered"
    '''
    
    # create directory if doesn't exist
    os.makedirs("filtered", exist_ok = True)

    # create path for the output file in filtered folder
    output_file_path = os.path.join("filtered", output_fastq)

    # checking gc, len, Q -> if passes, writes file
    with open(output_file_path, "w") as output_file:
        for read in read_fastq(input_fastq):
            _, seq, _, quality = read

            if (
                check_gc(seq, gc_bounds)
                and check_length(seq, length_bounds)
                and check_quality(quality, quality_threshold)
            ):
                write_fastq(output_file, read)



def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str | None = None
) -> None:
    
    '''
    Arranges many line FASTA into one line
    
    params:

    input_fasta - path to the intial multi line FASTA-file
    
    output_fasta - path to the result file (should be created automatically, 
    if not explicitly given)

    New FASTA-file is created.
    '''
    # changing the file name for returning (taking from initial)
    if output_fasta is None:
        output_fasta = splitext(input_fasta)[0] + "_oneline.fasta"

    with open(output_fasta, "w") as output:
        
        # multiline -> one line
        for header, sequence in read_fasta(input_fasta):
            output.write(header + "\n")
            output.write(sequence + "\n")




def parse_blast_output(input_file: str, 
                       output_file: str) -> None:
    
    '''
    Takes best alignments from BLAST
    
    params:

    input_file_path - txt_file path
    output_file_path - resting file path

    Saves proteins' names - one per line
    '''
    
    best_hits = get_best_blast_hits(input_file)

    with open(output_file, "w") as file:
        for protein in sorted(best_hits):
            file.write(protein + "\n")


