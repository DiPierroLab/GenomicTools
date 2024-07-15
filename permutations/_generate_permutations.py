import numpy as np
from scipy import sparse
from GenomicTools.tools import *
from GenomicTools.dot_plots import *

def generate_permutations(chrom_info_A, chrom_info_B, swap_batches, batch_size = 1000):
    """
    It's easy to generate and execute disjoint swaps. For small numbers of random swaps, the probability the swaps will not be disjoint is small. So, for large numbers of swaps, we can generate a series of disjoint swaps, which should be faster than doing one swap at a time.
    """
    if (swap_batches != np.inf) and ((type(swap_batches) == int) or (type(swap_batches) == float)):
        swap_batches = int(swap_batches)
        if swap_batches < 0:
            raise ValueError("If 'swap_batches' is an integer, it must be >= 0.")
    else:
        swap_batches = np.inf
    
    chromsA = alphanum_sort(list(chrom_info_A.keys()))
    chromsB = alphanum_sort(list(chrom_info_B.keys()))

    N_genes_A = np.sum([chrom_info_A[c]['size'] for c in chromsA])
    N_genes_B = np.sum([chrom_info_B[c]['size'] for c in chromsB])

    if swap_batches != np.inf:
        base_perm_A = np.arange(N_genes_A)
        base_perm_B = np.arange(N_genes_B)
        perm_A = np.copy(base_perm_A)
        perm_B = np.copy(base_perm_B)
        for batch in range(swap_batches):
            swaps_A = np.random.choice(N_genes_A,size=(batch_size,2),replace=False)
            swaps_B = np.random.choice(N_genes_B,size=(batch_size,2),replace=False)
            perm_A[swaps_A[:,0]] = base_perm_A[swaps_A[:,1]]
            perm_A[swaps_A[:,1]] = base_perm_A[swaps_A[:,0]]
            perm_B[swaps_B[:,0]] = base_perm_B[swaps_B[:,1]]
            perm_B[swaps_B[:,1]] = base_perm_B[swaps_B[:,0]]
            base_perm_A = np.copy(perm_A)
            base_perm_B = np.copy(perm_B)
    else:
        perm_A = np.random.permutation(N_genes_A)
        perm_B = np.random.permutation(N_genes_B)
    
    return perm_A, perm_B
