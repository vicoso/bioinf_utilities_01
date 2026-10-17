
def read_fastq(path: str):
    
    '''
    Reads fastq file per read
    
    params:
    path - path for file
    '''
    with open(path) as file:
        while True:
            name = file.readline()
            
            if len(name) == 0:
                break
            
            seq = file.readline()
            separator = file.readline()
            quality = file.readline()

            yield name, seq, separator, quality
            
def check_bounds(calculated_val: float, 
                 bounds: tuple | int | float) -> bool:
    
    '''
    Checks if tthe given value fits the bounds
    
    params:
    
    calculated_val - can be percentage or actual numbers of something of interest
    
    bounds - allowed thresholds for the values, can be percentage or actual numbers
    '''
    if isinstance(bounds, (int, float)):
        return 0 <= calculated_val <= bounds
    
    else:
        return bounds[0] <= calculated_val <= bounds[1]




def check_gc(seq: str, 
             bounds: tuple | int | float) -> bool:
    
    '''
    Checks GC of the read
    
    params:
    
    seq - a string of letters (big or small)
    bounds - allowed percentage of GC
    
    '''
    # make all letters big
    seq = seq.strip().upper()
    
    
    # no read -> don't use
    if len(seq) == 0:
        return False

    gc_count_per = (seq.count("G") + seq.count("C")) / len(seq) * 100
    
    gc_check_result = check_bounds(gc_count_per, bounds)

    return gc_check_result



def check_length(seq: str, 
                 bounds: tuple | int | float) -> bool:
    
    '''
    Checks length of the sequence (seq)
    
    params:
    
    seq - 2nd line of the read
    
    bounds - allowed thresholds
    '''
    length = len(seq.strip())

    return check_bounds(length, bounds)


def check_quality(quality: str, 
                  threshold: float) -> bool:
    
    '''
    checks avg quality of the read according to Phred+33
    
    params:
    
    quality - what we get rom fastq file (4th line per read)
    
    threshold - allowed score thrshold
    '''
    quality = quality.strip()

    if len(quality) == 0:
        return False

    # recalculating from fastq into ASCII, and thenn - 33 into Phred+33
    scores = [ord(char) - 33 for char in quality]
    mean_quality = sum(scores) / len(scores)

    return mean_quality >= threshold
 
 
def write_fastq(file, 
                read: tuple) -> None:
    
    '''
    Writes down one read at a time in FASTQ file
    
    params:
    
    file - fastq file we want to write down
    
    read - read (4 lines) fron the initial fastq file
    '''
    file.writelines(read)
