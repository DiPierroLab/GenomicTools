import numpy as np
from scipy import sparse
from GenomicTools.tools import *
from GenomicTools.dotplots import *

def generate_permutations(dot_plot_result, swap_batches, batch_size = 1000):
    if (swap_batches != np.inf) and ((type(swap_batches) == int) or (type(swap_batches) == float)):
        swap_batches = int(swap_batches)
        if swap_batches < 0:
            raise ValueError("If 'swap_batches' is an integer, it must be >= 0.")
    else:
        swap_batches = np.inf
    
    chroms = np.vstack(list(dot_plot_result['data'].keys()))
    chromsA = alphanum_sort(np.unique(chroms[:,0]))
    chromsB = alphanum_sort(np.unique(chroms[:,1]))

    if dot_plot_result['species2'] in ['Homo_sapiens', 'Mus_musculus']:
        dummy_chrom2 = 'chr1'
    else:
        dummy_chrom2 = 'HiC_scaffold_1'
    if dot_plot_result['species1'] in ['Homo_sapiens', 'Mus_musculus']:
        dummy_chrom1 = 'chr1'
    else:
        dummy_chrom1 = 'HiC_scaffold_1'

    N_genes_A = np.sum([dot_plot_result['data'][(c,dummy_chrom2)]['homology_matrix'].shape[1] for c in chromsA])
    N_genes_B = np.sum([dot_plot_result['data'][(dummy_chrom1,c)]['homology_matrix'].shape[0] for c in chromsB])

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

# If state is an integer, it's the number of swaps the permutation is from the native permutation (operationally, not the minimum swaps). It's easy to generate and execute disjoint swaps. For small numbers of random swaps, the probability the swaps will not be disjoint is small. So, for large numbers of swaps, we can generate a series of disjoint swaps, which should be faster than doing one swap at a time.
