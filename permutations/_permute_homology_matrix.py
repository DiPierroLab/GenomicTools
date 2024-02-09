import numpy as np
import copy
from scipy import sparse
from GenomicTools.tools import *
from GenomicTools.dotplots import *

def permute_homology_matrix(homology_matrix, species_data_A, species_data_B, perm_A, perm_B):
    chromsAB = np.vstack(list(homology_matrix.keys()))
    chromsA = alphanum_sort(np.unique(chromsAB[:,0]))
    chromsB = alphanum_sort(np.unique(chromsAB[:,1]))
    dummy_chrom_A = chromsA[0]
    dummy_chrom_B = chromsB[0]

    chrom_info_A = get_chrom_info(species_data_A)
    chrom_info_B = get_chrom_info(species_data_B)
    chrom_locs_A = np.cumsum([0] + [chrom_info_A[c]['size'] for c in chromsA])
    chrom_locs_B = np.cumsum([0] + [chrom_info_B[c]['size'] for c in chromsB])
    chrom_indices_A = {c:chrom_locs_A[i:(i+2)] for i,c in enumerate(chromsA)}
    chrom_indices_B = {c:chrom_locs_B[i:(i+2)] for i,c in enumerate(chromsB)}    

    dot_matrix = []
    for n, chromA in enumerate(chromsA):
        dot_matrix_col = []
        for chromB in chromsB:
            dot_matrix_col.append(homology_matrix[(chromA,chromB)])
        dot_matrix.append(sparse.vstack(dot_matrix_col)) 
    dot_matrix = sparse.csr_matrix(sparse.hstack(dot_matrix)) 
    permuted_dot_matrix = dot_matrix[perm_A,:][:,perm_B]

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

    return permuted_homology_matrix, permuted_species_data_A, permuted_species_data_B
