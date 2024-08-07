import numpy as np
import copy
from scipy.stats import pearsonr
from scipy import sparse
from GenomicTools.tools import *
from GenomicTools.dot_plots import *

def permute_homology_matrix(homology_matrix, species_data_A, species_data_B, perm_A, perm_B):
    chrom_info_A = get_chrom_info(species_data_A)
    chrom_info_B = get_chrom_info(species_data_B)

    chromsA = alphanum_sort(list(chrom_info_A.keys()))
    chromsB = alphanum_sort(list(chrom_info_B.keys()))
    
    chrom_locs_A = np.cumsum([0] + [chrom_info_A[c]['size'] for c in chromsA])
    chrom_locs_B = np.cumsum([0] + [chrom_info_B[c]['size'] for c in chromsB])
    chrom_indices_A = {c:chrom_locs_A[i:(i+2)] for i,c in enumerate(chromsA)}
    chrom_indices_B = {c:chrom_locs_B[i:(i+2)] for i,c in enumerate(chromsB)}    

    dot_matrix = []
    for n, chromA in enumerate(chromsA):
        dot_matrix_col = []
        for chromB in chromsB:
            if (chromA,chromB) in homology_matrix.keys():
                dot_matrix_col.append(homology_matrix[(chromA,chromB)])
        if len(dot_matrix_col) > 0:
            dot_matrix.append(sparse.hstack(dot_matrix_col))
    dot_matrix = sparse.csr_matrix(sparse.vstack(dot_matrix))
    permuted_dot_matrix = dot_matrix[perm_A,:][:,perm_B]

    rx = pearsonr(np.array(dot_matrix.sum(0))[0,1:],np.array(dot_matrix.sum(0))[0,:-1])[0]
    ry = pearsonr(np.array(dot_matrix.sum(1))[1:,0],np.array(dot_matrix.sum(1))[:-1,0])[0]
    rx_perm = pearsonr(np.array(permuted_dot_matrix.sum(0))[0,1:],np.array(permuted_dot_matrix.sum(0))[0,:-1])[0]
    ry_perm = pearsonr(np.array(permuted_dot_matrix.sum(1))[1:,0],np.array(permuted_dot_matrix.sum(1))[:-1,0])[0]

    permuted_homology_matrix = copy.deepcopy(homology_matrix)
    for chromA in chromsA:
        for chromB in chromsB:
            a1, a2 = chrom_indices_A[chromA]
            b1, b2 = chrom_indices_B[chromB]
            permuted_homology_matrix[(chromA,chromB)] = permuted_dot_matrix[a1:a2,:][:,b1:b2]

    permuted_species_data_A = np.copy(species_data_A)    
    permuted_species_data_B = np.copy(species_data_B)
    
    permuted_species_data_A = np.hstack([species_data_A[:,:4],species_data_A[perm_A,4:]])
    permuted_species_data_B = np.hstack([species_data_B[:,:4],species_data_B[perm_B,4:]])

    return permuted_homology_matrix, permuted_species_data_A, permuted_species_data_B, np.array([rx, ry, rx_perm, ry_perm])
