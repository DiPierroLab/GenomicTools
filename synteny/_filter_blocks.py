
import numpy as np
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *
from scipy.signal import fftconvolve

def supported_by_nanosynteny(block, species_data_A, species_data_B, maps_A, maps_B, nanosynteny_minsize, filter_by_unique_orthogroups = True):
    """
    Input:
    - block: numpy array (k X 4), where k is the number of genes in the block, in relative coordinates.
    - species_data_A: numpy array (N X 12), species data for species A.
    - species_data_B: numpy array (M X 12), species data for species B.
    - maps_A: [cc_map, inv_cc_map, shift_map, unshift_map], output of create_shift_map(species_data_A, windowsize).
        - cc_map: gene index -> cc index
        - inv_cc_map: cc index -> set of gene indices
        - shift_map: cc index -> compressed gene index
        - unshift_map: compressed gene index -> cc index
    - maps_B: [cc_map, inv_cc_map, shift_map, unshift_map], output of create_shift_map(species_data_B, windowsize), see above.
    - nanosynteny_minsize: int, minimum number of genes in a nanosynteny block.
    - filter_by_unique_orthogroups: bool, if True the supporting nanosynteny block must have nanosynteny_minsize unique orthogroups within it for the synteny block in question to be "supported".
    
    Output:
    - supported: bool, is the synteny block supported by nanosynteny?
    """
    chromA = block[0,0]
    chromB = block[0,2]
    x = block[:,1].astype(int)
    y = block[:,3].astype(int)
    xdiff = np.abs(np.diff(x))
    ydiff = np.abs(np.diff(y))
    kernel = np.ones(nanosynteny_minsize - 1)
    xconv = np.round(fftconvolve(xdiff,kernel)[1:])
    yconv = np.round(fftconvolve(ydiff,kernel)[1:])
    xconv[-1] = xconv[-2]
    yconv[-1] = yconv[-2] 
    
    if filter_by_unique_orthogroups == False:
        supported = np.any((xconv == nanosynteny_minsize - 1)*(yconv == nanosynteny_minsize - 1))
        return supported
    
    elif filter_by_unique_orthogroups == True:
        cc_maps, inv_cc_maps, shift_maps, unshift_maps = maps_A
        inv_cc_map = inv_cc_maps[chromA]
        unshift_map = unshift_maps[chromA]
        supported_spots = (xconv == nanosynteny_minsize - 1)*(yconv == nanosynteny_minsize - 1).astype(int)
        edges = ((xdiff == 1) * (ydiff == 1)).astype(int)
        A = np.diag(np.zeros(block.shape[0])) + np.diag(edges,1) + np.diag(edges,-1)
        G = ntx.Graph(A)
        nano_ccs = list(ntx.connected_components(G))
        supported = False
        for nano_cc in nano_ccs:
            nano_cc_genes = block[list(nano_cc)][:,1]
            species_data_chrom = species_data_A[species_data_A[:,0] == chromA]
            nano_cc_ogs = np.unique([species_data_chrom[inv_cc_map[unshift_map[int(g)]]-1][0,4] for g in nano_cc_genes])            
            if len(nano_cc_ogs) >= nanosynteny_minsize:
                supported = True
                break
        return supported    
   
def synteny_overlap(blockA, blockB, spA, spB):
    """
    This function takes two synteny blocks from the same homology matrix (comparing the same two species) and returns 
    the overlap in number of genes. 
    
    Each synteny blocks must be in the format:
        - first column: chromosome from species A
        - second column: relative gene indices from species A
        - third column: chromosome from species B
        - fourth column: relative gene indices from species B
    
    Input:
        - blockA: N X 4 numpy array, synteny block in relative indices
        - blockB: M X 4 numpy array, synteny block in relative indices
        - spA: str, name of species A
        - spB: str, name of species B
    Output:
        - overlap: int, number of genes by which the two blocks overlap
    """
    if (blockA.shape[1] == 4) and (blockB.shape[1] == 4):
        species_same = (spA == spB)
        same_chroms_case_1 = ((blockA[0,0] == blockB[0,0]) and (blockA[0,2] == blockB[0,2]) and (blockA[0,0] == blockA[0,2]))
        same_chroms_case_2 = ((blockA[0,0] == blockB[0,0]) and (blockA[0,2] == blockB[0,2]) and (blockA[0,0] != blockA[0,2]))
        same_chroms_case_3 = ((blockA[0,0] == blockB[0,2]) and (blockA[0,2] == blockB[0,0]) and (blockA[0,0] != blockA[0,2]))
        if species_same:
            # In this case, we're looking at self-synteny, comparing a genome with itself (under the assumption that we have one genome per species)
            if same_chroms_case_1:
                # Case 1a: With self-synteny, both blocks result from comparison of a chromosome with itself
                x_overlap1 = len(set(blockA[:,1]).intersection(set(blockB[:,1])))
                y_overlap1 = len(set(blockA[:,3]).intersection(set(blockB[:,3])))
                x_overlap2 = len(set(blockA[:,1]).intersection(set(blockB[:,3])))
                y_overlap2 = len(set(blockA[:,3]).intersection(set(blockB[:,1])))
                overlap1 = np.min([x_overlap1, y_overlap1])
                overlap2 = np.min([x_overlap1, y_overlap1])
                overlap = np.max([overlap1, overlap2])
                return overlap
            elif same_chroms_case_2: 
                # Case 2a: With self-synteny, the first chromosomes (first and second columns) and second chromosomes (third and fourth columns) match across blocks
                x_overlap = len(set(blockA[:,1]).intersection(set(blockB[:,1])))
                y_overlap = len(set(blockA[:,3]).intersection(set(blockB[:,3])))
                overlap = np.min([x_overlap, y_overlap])
                return overlap
            elif same_chroms_case_3:
                # Case 3a: With self-synteny, the first chromosome (first and second columns) of one block and the second chromosome (third and fourth columns) of the other match, and same for the remaining two chromosomes
                x_overlap = len(set(blockA[:,1]).intersection(set(blockB[:,3])))
                y_overlap = len(set(blockA[:,3]).intersection(set(blockB[:,1])))
                overlap = np.min([x_overlap, y_overlap])
                return overlap
            else:
                # Case 4a: The two blocks cannot overlap since they don't share both chromosomes
                return 0
        else:
            if same_chroms_case_2:
                # Case 2b: Here we have cross-species synteny, and we have matching chromosomes in matching columns
                x_overlap = len(set(blockA[:,1]).intersection(set(blockB[:,1])))
                y_overlap = len(set(blockA[:,3]).intersection(set(blockB[:,3])))
                overlap = np.min([x_overlap, y_overlap])
                return overlap
            else:
                # Case 4b: The two blocks cannot overlap since they don't share both chromosomes
                return 0
    else:
        raise ValueError("blockA and blockB must both be in relative coordinates.")
        
def return_overlapping_duplication_pairs(blocks, spA, spB, overlap_threshold):
    """
    This function takes a list of synteny blocks from the same homology matrix (comparing the same two species) and returns 
    indices of pairs of blocks that overlap by at least overlap_threshold genes. 
    
    Input:
        - blocks: list of N X 4 numpy arrays, synteny blocks in relative indices
        - spA: str, name of species A
        - spB: str, name of species B
        - overlap_threshold: int, number of genes by which two blocks must overlap for a conflict to be found
    Output:
        - overlap: numpy array, first two columns have indices of overlapping blocks, third column has the overlap in number of genes
    """
    overlaps = []
    for i in range(len(blocks)):
        for j in range(i+1,len(blocks)):
            overlap = synteny_overlap(blocks[i],blocks[j],spA,spB)
            if (overlap >= overlap_threshold):
                overlaps.append([i,j,overlap])
    if len(overlaps) > 0:
        return np.vstack(overlaps)
    else:
        return np.array(overlaps)

def remove_duplication_overlaps(condensed_duplications, overlap_threshold = 1):
    """
    This function takes as input a list of duplications (the duplications output of identify_duplications), which are a subset of the 
    self-synteny blocks for a species, in condensed indices. It produces a list of condensed duplications with conflicts generated
    by overlapping blocks removed - the largest block in a connected component of overlapping blocks is chosen to be retained.
    
    Input:
        - condensed_duplications: list of N X 4 numpy arrays, synteny blocks in relative, condensed indices
        - overlap_threshold: int, number of genes by which two blocks must overlap for a conflict to be found
    Output:
        - overlap: numpy array, first two columns have indices of overlapping blocks, third column has the overlap in number of genes
    """
    spA = 'species'
    spB = 'species'
    remove_list = []
    N = len(condensed_duplications)
    overlaps = return_overlapping_duplication_pairs(condensed_duplications,spA,spB,overlap_threshold)
    for o in overlaps:
        b1 = condensed_duplications[o[0]]
        b2 = condensed_duplications[o[1]]
        if b1.shape[0] < b2.shape[0]:
            remove_list.append(o[0])
        else:
            remove_list.append(o[1])
    duplications_no_overlaps = [condensed_duplications[n] for n in range(N) if n not in remove_list]
    return duplications_no_overlaps 
