import numpy as np
import re

def alphanum_sort(l): 
    """
    Taken directly from here: https://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
    
    Input:
        - l: an array or list for sorting
        
    Output:
        A sorted version of 'l' accounting for the usual integer and alphabetic ordering
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key = alphanum_key)

def chrom_to_num(name):
    """
    Convert Hi-C scaffold (chromosome) name into integer index. All DNAZoo Hi-C scaffold names have the 
    form 'HiC_scaffold_{i}' where {i} is an integer. The human and mouse genomes have chromosome names
    'chr{i}' where {i} is an integer or is X, for X chromosomes. We map X chromosomes to 999.
    
    Input:
        - name: a string with the name of the Hi-C scaffold / chromosome, either in the form 'HiC_scaffold_{i}', 
        'chr{i}', or 'chrX', where {i} is an integer
    
    Output:
        An integer index of the Hi-C scaffold / chromosome, either the integer in the name or 999 for 
        explicitly named X chromosomes
    """
    if 'HiC' in name:
        return int(name.split('_')[-1])
    else:
        ID = name[3:]
        if ID == 'X':
            return 999
        else:
            return int(ID)

def num_to_chrom(num, sp):
    """
    Convert Hi-C scaffold (chromosome) integer index into name. The 'sp' variable should be in the 
    form 'Genus_species'.
    
    Input:
        - num: integer Hi-C scaffold / chromosome index, where 999 is an X chromosome
        - sp: species name, in the form 'Genus_species'
        
    Output:
        Chromosome name, in the form 'HiC_scaffold_{num}', 'chr{num}', or 'chrX'
    """
    if (sp == 'Homo_sapiens') or (sp == 'Mus_musculus'):
        if num == 999:
            return 'chrX'
        else:
            return 'chr%i'%num
    else:
        return 'HiC_scaffold_%i'%num

def convert_to_absolute_indices(blocks, dot_plot):
    """
    Convert a list of synteny blocks with separate gene indexing for each chromosome into a list of synteny blocks
    with absolute indexing (which simply indexes all genes in an organism from 1 to the total number of genes).
    
    Input:
        - blocks: list of N X 4 arrays, with
            - column 1: chromosome index for species A 
            - column 2: gene index for species A within specified chromosome
            - column 3: chromosome index for species B
            - column 4: gene index for species B within specified chromosome
        - dot_plot: a dictionary containing dot plot results for species A and B. 
                    Note: need to add where to find a complete description of dot plot format.
    
    Output:
        List of N X 2 arrays, where each array contains the absolute coordinates of a block, containing N genes.
    """
    chroms = np.vstack(list(dot_plot['data'].keys()))
    chrom_names1 = alphanum_sort(np.unique(chroms[:,0]))
    chrom_names2 = alphanum_sort(np.unique(chroms[:,1]))
    
    lens1 = [dot_plot['data'][(name1,chrom_names2[0])]['homology_matrix'].toarray().shape[1] for name1 in chrom_names1]
    lens2 = [dot_plot['data'][(chrom_names1[0],name2)]['homology_matrix'].toarray().shape[0] for name2 in chrom_names2]
    
    cumulative_lens1 = np.cumsum([0] + lens1)
    cumulative_lens2 = np.cumsum([0] + lens2)
    
    converted = []
    for block in blocks:
        if block[0,0] == 999:
            chrom1 = len(lens1)
        else:
            chrom1 = block[0,0]
        if block[0,2] == 999:
            chrom2 = len(lens2)
        else:
            chrom2 = block[0,2]
        
        add1 = cumulative_lens1[chrom1-1]
        add2 = cumulative_lens2[chrom2-1]
      
        converted_block = np.hstack([block[:,1:2]+add1,block[:,3:4]+add2])
        converted.append(converted_block)
        
    return converted

def interval_overlap(I1, I2):
    """
    Take two "intervals" of integers and return the size of the overlap. Inputs are lists [a,b] where both a and b are integers
    and a <= b, representing all integers a, a + 1, ..., b - 1, b. The output is the cardinality of the intersection of the two 
    input sets.

    Input:
        - I1: A list with two elements, [a,b] where a and b are integers and a <= b.
        - I2: A list with two elements, [a,b] where a and b are integers and a <= b.
    
    Output:
        The number of integers in the overlap between the two input "intervals".
    """
    if (I1[1] < I1[0]) or (I2[1] < I2[0]):
        raise ValueError("These aren't invervals written as [a,b] with a <= b")
    if (I1[1] < I2[0]) or (I2[1] < I1[0]):
        return 0
    elif I1[0] == I2[0]:
        return np.min([I1[1], I2[1]]) - I1[0]
    elif I1[1] == I2[1]:
        return I1[1] - np.max([I1[0], I2[0]]) + 1
    elif (I1[0] < I2[0]) and (I1[1] < I2[1]):
        return I1[1] - I2[0] + 1
    elif (I1[0] < I2[0]) and (I1[1] > I2[1]):
        return I2[1] - I2[0] + 1
    elif (I1[0] > I2[0]) and (I1[1] < I2[1]):
        return I1[1] - I1[0] + 1
    elif (I1[0] > I2[0]) and (I1[1] > I2[1]):
        return I2[1] - I1[0] + 1

def load_dagchainer_results(dagchainer_file):
    """
    Load the output of DAGchainer into a list of arrays, each array describing a synteny block.
    Note: This needs better documentation. 
    """
    with open(dagchainer_file,'r') as f:
            lines = f.readlines()
    breaks = np.where(np.array(['#' in l for l in lines]) == True)[0]

    blocks = []
    scores = []
    for n in range(len(breaks)):
            scores.append(int(lines[breaks[n]].split('=')[1].split(' ')[1].split('.')[0]))
            try:
                    block_list = lines[breaks[n]+1:breaks[n+1]]
            except IndexError:
                    block_list = lines[breaks[n]+1:]
            fixed_block = []
            for block in block_list:
                    block_cols = np.array(block.split('\t'))[[0,2,4,6]]
                    if 'B' in block_cols[0]:
                            block_row = np.array([chrom_to_num(block_cols[2].strip('AB')),block_cols[3],chrom_to_num(block_cols[0].strip('AB')),block_cols[1]]).astype(int)
                    else:
                            block_row = np.array([chrom_to_num(block_cols[0].strip('AB')),block_cols[1],chrom_to_num(block_cols[2].strip('AB')),block_cols[3]]).astype(int)
                    fixed_block.append(block_row)
            fixed_block = np.vstack(fixed_block)
            blocks.append(fixed_block)
    return blocks, scores
