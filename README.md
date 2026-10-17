Scripts for handling bioinf data in such formats as fastq, fasta, blast table. 
Scripts use only the Python standard libraries.

The main script is located in `main.py`.

It contains the following functions:

1) `filter_fastq` - filters FASTQ reads by GC content, length of sequence, and average quality score (Phred+33);

2) `convert_multiline_fasta_to_oneline` - converts multi-line FASTA sequences into single-line sequences;

3) `parse_blast_output` - extracts the best protein hit for each query from a BLAST report and saves the descriptions in alphabetical order.


Folder "helpers" has supporting scripts for the main one.

