
import numpy as np
from GenomicTools.tools import *

def block_to_string_relative(block):
    """
    Convert a synteny block in relative indices (with chromosome names and indices relative to chromosomes) into strings.

    Input:
        - block: N X 4 array, synteny block

    Output:
        - bs: string, synteny block with pairs of genes are separated by '_'
        - bsT: string, synteny block with species swapped and with pairs of genes are separated by '_'
        - bsTf: string, synteny block with species swapped, in opposite order, and with pairs of genes are separated by '_'
    """
    b1 = np.char.add(block[:,0].astype(str),'-')
    b2 = np.char.add(b1,block[:,1].astype(str))
    b3 = np.char.add(block[:,2].astype(str),'-')
    b4 = np.char.add(b3,block[:,3].astype(str))
    bb1 = np.char.add(b2,'-')
    bb2 = np.char.add(bb1,b4)
    bs = '_'.join(bb2)

    bb1T = np.char.add(b4,'-')
    bb2T = np.char.add(bb1T,b2)
    bsT = '_'.join(bb2T)

    bb1Tf = np.char.add(b4[::-1].astype(str),'-')
    bb2Tf = np.char.add(bb1Tf,b2[::-1].astype(str))
    bsTf = '_'.join(bb2Tf)

    return bs, bsT, bsTf

def block_to_string_absolute(block):
    """
    Convert a synteny block in absolute indices (unique indices for all genes) into strings.

    Input:
        - block: N X 2 array, synteny block

    Output:
        - bs: string, synteny block with pairs of genes are separated by '_'
        - bsT: string, synteny block with species swapped and with pairs of genes are separated by '_'
        - bsTf: string, synteny block with species swapped, in opposite order, and with pairs of genes are separated by '_'
    """
    b1 = np.char.add(block[:,0].astype(str),'-')
    b2 = np.char.add(b1,block[:,1].astype(str))
    bs = '_'.join(b2)

    b1T = np.char.add(block[:,1].astype(str),'-')
    b2T = np.char.add(b1T,block[:,0].astype(str))
    bsT = '_'.join(b2T)

    b1Tf = np.char.add(block[::-1,1].astype(str),'-')
    b2Tf = np.char.add(b1Tf,block[::-1,0].astype(str))
    bsTf = '_'.join(b2Tf)

    return bs, bsT, bsTf

def block_palindrome(block):
    """
    Is the synteny block a palindrome when comparing a species with itself?

    Input:
        - block: N X 4 array, synteny block in relative indices

    Output:
        - palindrome: boolean
    """
    if block[0,0] == block[0,2]:
        indicesA = '-'.join(block[:,1])
        indicesB = '-'.join(block[::-1,3])
        if indicesA == indicesB:
            return True
        else:
            return False
    else:
        return False

def block_palindromoid(block):
    """
    Is the synteny block a palindromoid (is it matching a region partially with itself) when comparing a species with itself?

    Input:
        - block: N X 4 array, synteny block in relative indices

    Output:
        - palindromoid: boolean
    """
    if block[0,0] == block[0,2]:
        gene_intersection = set(block[:,1]).intersection(block[:,3])
        if len(gene_intersection) > 0:
            return True
        else:
            return False
    else:
        return False

def blocks_symmetric(blocks):
    """
    For a set of synteny blocks, which blocks have their species-swapped or order-reversed counterpart in the set?

    Input:
        - blocks: list of N X 4 array, synteny blocks in relative indices

    Output:
        - symm: list of boolean, is swapped or flipped counterpart in the set of blocks?
    """
    n_block_columns = blocks[0].shape[1]
    string_blocks = []
    transpose_string_blocks = []
    transpose_flip_string_blocks = []

    if n_block_columns == 4:
        for b in blocks:
            bs, bsT, bsTf = block_to_string_relative(b)
            string_blocks.append(bs)
            transpose_string_blocks.append(bsT)
            transpose_flip_string_blocks.append(bsTf)

    elif n_block_columns == 2:
        for b in blocks:
            bs, bsT, bsTf = block_to_string_absolute(b)
            string_blocks.append(bs)
            transpose_string_blocks.append(bsT)
            transpose_flip_string_blocks.append(bsTf)

    symm = []
    for n in range(len(string_blocks)):
        transpose = (transpose_string_blocks[n] in string_blocks)
        transpose_flip = (transpose_flip_string_blocks[n] in string_blocks)
        symm.append(transpose or transpose_flip)

    return symm

def self_diagonal_block(block):
    """
    Is a synteny block found by comparing a genome to itself just mathching genes with themselves?
    
    Input:
        - block: N X 4 array, synteny block in relative indices

    Output:
        - is_self_diagonal: boolean
    """
    n_block_columns = block.shape[1]
    if n_block_columns == 4:
        A0 = np.char.add(block[:,0],'-')
        A = np.char.add(A0,block[:,1])
        B0 = np.char.add(block[:,2],'-')
        B = np.char.add(B0,block[:,3])
    elif n_block_columns == 2:
        A = block[:,0]
        B = block[:,1]
    if np.all(A == B):
        return True
    else:
        return False

def synteny_overlap(blockA, blockB):
    """
    In a comparison of two genomes, how much do two synteny blocks overlap? The input can be in either relative
    or absolute indices, but must be consistent between the two blocks.

    Input:
        - blockA: N X 4 array or N X 2 array, synteny block in relative indices or absolute indices
        - blockB: N X 4 array or N X 2 array, synteny block in relative indices or absolute indices

    Output:
        - overlap: integer, number of genes overlapping between the two blocks
    """
    if (blockA.shape[1] == 4) and (blockB.shape[1] == 4):
        if (blockA[0,0] != blockB[0,0]) or (blockA[0,2] != blockB[0,2]):
            return 0
        else:
            x_overlap = len(set(blockA[:,1]).intersection(set(blockB[:,1])))
            y_overlap = len(set(blockA[:,3]).intersection(set(blockB[:,3])))
            if (x_overlap > 0) and (y_overlap > 0):
                return np.min([x_overlap, y_overlap])
            else:
                return 0
    elif (blockA.shape[1] == 2) and (blockB.shape[1] == 2):
        x_overlap = len(set(blockA[:,0]).intersection(set(blockB[:,0])))
        y_overlap = len(set(blockA[:,1]).intersection(set(blockB[:,1])))
        if (x_overlap > 0) and (y_overlap > 0):
            return np.min([x_overlap, y_overlap])
        else:
            return 0
    else:
        raise ValueError("blockA and blockB need to be both in absolute or relative coordinates.")
        
def return_overlapping_block_pairs(blocks, overlap_threshold):
    """
    For a set of synteny blocks, return all overlapping pairs.

    Input:
        - blocks: list of N X 4 array, synteny blocks in relative indices
        - overlap_threshold: integer, how many genes must two blocks overlap by to be "overlapping"?

    Output:
        - overlaps: list of arrays with three elements, the indices of the two overlapping blocks and how many genes they overlap by
    """
    overlaps = []
    for i in range(len(blocks)):
        for j in range(i+1,len(blocks)):
            overlap = synteny_overlap(blocks[i],blocks[j])
            if (overlap >= overlap_threshold):
                overlaps.append([i,j,overlap])
    return overlaps
    
def block_complexity(absolute_block, species_data_A, species_data_B):
    """
    Return complexity information from a synteny block in absolute indices.

    Input:
        - absolute block: N X 2 array, synteny block in absolute indices
        - species_data_A: N X 12 array, species data array with complexity information added for species A
        - species_data_B: N X 12 array, species data array with complexity information added for species B

    Output:
        - complexity_A: N X 3 array, complexity information for synteny block genes in species A
        - complexity_B: N X 3 array, complexity information for synteny block genes in species B
    """
    complexity_A = species_data_A[absolute_block[:,0]-1,7:10]
    complexity_B = species_data_B[absolute_block[:,1]-1,7:10]
    return complexity_A, complexity_B

def return_low_complexity(absolute_blocks, species_data_A, species_data_B, entropy_fraction_threshold = .75, entropy_window_size = 10):
    """
    For a set of synteny blocks in absolute indices, return low complexity blocks.

    Input:
        - absolute_blocks: list of N X 2 arrays, where N varies from block to block
        - species_data_A: N X 12 array, species data array with complexity information added for species A
        - species_data_B: N X 12 array, species data array with complexity information added for species B
        - entropy_fraction_threshold: float, threshold for Shannon entropy below which a block is "low complexity" (Default = .75)
        - entropy_window_size: integer, window size with which to calculate average sliding window entropy (Default = 10)

    Output:
        - low_complexity: list of integers, block indices with low complexity
    """
    low_complexity = []
    threshold = entropy_fraction_threshold * np.log2(entropy_window_size)
    for n, block in enumerate(absolute_blocks):
        complexity_A, complexity_B = block_complexity(block, species_data_A, species_data_B)
        belowA = np.any(complexity_A[:,1].astype(float) < threshold)
        belowB = np.any(complexity_B[:,1].astype(float) < threshold)
        if belowA or belowB:
            low_complexity.append(n)
    return low_complexity

def convert_synteny_relative_to_absolute_indices(synteny_blocks, chrom_info_A, chrom_info_B):
    """
    Convert syntney blocks in relative indices to absolute indices.

    Input:
        - synteny_blocks: list of N X 4 arrays, synteny blocks in relative indices
        - chrom_info_A: dictionary, chromosome information for species A
        - chrom_info_B: dictionary, chromosome information for species B

    Output:
        - absolute blocks: list of N X 2 arrays, synteny blocks in absolute indices
    """
    chrom_locs_A = np.cumsum([0] + [chrom_info_A[key]['size'] for key in chrom_info_A.keys()])
    chrom_locs_B = np.cumsum([0] + [chrom_info_B[key]['size'] for key in chrom_info_B.keys()])
    if type(synteny_blocks[0][0,0]) in [int,np.int_]:
        abs_A = {n+1:s for n, s in enumerate(chrom_locs_A[:-1])}
        abs_B = {n+1:s for n, s in enumerate(chrom_locs_B[:-1])}
    elif type(synteny_blocks[0][0,0]) in [str,np.str_]:
        abs_A = {alphanum_sort(chrom_info_A.keys())[n]:s for n, s in enumerate(chrom_locs_A[:-1])}
        abs_B = {alphanum_sort(chrom_info_B.keys())[n]:s for n, s in enumerate(chrom_locs_B[:-1])}
    else:
        raise ValueError("Something is wrong with the format of synteny_blocks.")
    absolute_blocks = []
    for block in synteny_blocks:
        absolute_block = []
        for b in block:
            absolute_block.append([abs_A[b[0]]+int(b[1]),abs_B[b[2]]+int(b[3])])
        absolute_blocks.append(np.vstack(absolute_block))
    return absolute_blocks

def convert_synteny_absolute_to_relative_indices(synteny_blocks, chrom_info_A, chrom_info_B):
    """
    Convert synteny blocks in absolute indices to relative indices.

    Input:
        - synteny_blocks: list of N X 2 arrays, synteny blocks in absolute indices
        - chrom_info_A: dictionary, chromosome information for species A
        - chrom_info_B: dictionary, chromosome information for species B

    Output:
        - relative blocks: list of N X 4 arrays, synteny blocks in relative indices
    """
    chrom_locs_A = np.cumsum([0] + [chrom_info_A[key]['size'] for key in chrom_info_A.keys()])
    chrom_locs_B = np.cumsum([0] + [chrom_info_B[key]['size'] for key in chrom_info_B.keys()])
    abs_A = {alphanum_sort(chrom_info_A.keys())[n]:s for n, s in enumerate(chrom_locs_A[:-1])}
    abs_B = {alphanum_sort(chrom_info_B.keys())[n]:s for n, s in enumerate(chrom_locs_B[:-1])}
    abs_int_A = {n+1:s for n, s in enumerate(chrom_locs_A[:-1])}
    abs_int_B = {n+1:s for n, s in enumerate(chrom_locs_B[:-1])}
    relative_blocks = []
    for block in synteny_blocks:
        relative_block = []
        chromA = np.digitize(block[0,0]-1,chrom_locs_A)
        chromB = np.digitize(block[0,1]-1,chrom_locs_B)
        for b in block:
            relative_block.append([chromA,b[0]-abs_int_A[chromA],chromB,b[1]-abs_int_B[chromB]])
        relative_blocks.append(np.vstack(relative_block))
    return relative_blocks
