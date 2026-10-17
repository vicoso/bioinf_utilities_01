
def read_fasta(input_fasta: str):
    
    '''
    Reads multi-line FASTA file and returns one line FASTA files
    
    params:
    
    input_fasta - initial fasta multi-lined seq
    '''
    
    with open(input_fasta) as file:
        header_str = None

        for line in file:
            
            # get rid of \n
            line = line.strip()
            
            # checking if header
            if line.startswith(">"):
                
                # for checking if we meet a new header
                if header_str is not None:
                    
                    # glueing into one line
                    yield header_str, "".join(seq_list)

                header_str = line
                seq_list = []

            elif header_str is not None:
                seq_list.append(line)
                

        yield header_str, "".join(seq_list)


def get_best_blast_hits(input_file: str) -> list[str]:
    
    '''
    Takes out best descriptions from BLAST-report
    
    params:
    input_file - path to the file with BLAST results
    '''
    
    best_hits_list = []
    looking_for_hit = False
    description_end = None

    with open(input_file) as file:
        for line in file:
            
            # looking for table start
            if line.startswith("Sequences producing significant alignments:"):
                looking_for_hit = True
                description_end = None
                continue
            
            # skip irrelevant lines
            if not looking_for_hit:
                continue

            # calculate description end
            if line.lstrip().startswith("Scientific"):
                description_end = line.index("Scientific")
                continue
            
            # skip empty lines
            if description_end is None or not line.strip():
                continue

            # skip column name
            if line.startswith("Description"):
                continue

            # extract wanted line
            best_hits_list.append(line[:description_end].strip())
            looking_for_hit = False

    return best_hits_list
