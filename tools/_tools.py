
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
        return np.min([I1[1], I2[1]]) - I1[0] + 1
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

def block_slope(block):
    x = block[:,1].astype(int)
    y = block[:,3].astype(int)
    slope = (y[-1] - y[0]) / (x[-1] - x[0])
    return slope

def get_chrom_info(species_data):
    chrom_info = {}
    chroms = alphanum_sort(np.unique(species_data[:,0]))
    for n, chrom in enumerate(chroms):
        chrom_info[chrom] = {}
        chrom_info[chrom]['number'] = n + 1
        chrom_info[chrom]['size'] = species_data[species_data[:,0] == chrom].shape[0]
    return chrom_info

def parse_tree(tree_string):
    tree = Phylo.read(StringIO(tree_string),'newick')
    return tree
