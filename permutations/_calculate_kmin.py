
import numpy as np
from scipy import stats

def calculate_density(dot_plot, chrom_info_A, chrom_info_B):
    """
    Calculate the density of dots in a dot plot. This is total number of dots (number of homology relations between genes)
    divided by the total number of possible dots (the number of genes in species A times the number of genes in species B).

    Input:
        - dot_plot: N X 9 array, dot plot array
        - chrom_info_A: dictionary, chromosome information dictionary for species A
        - chrom_info_B: dictionary, chromosome information dictionary for species B

    Output:
        - rho: float, dot density
        - N_A: total number of genes in species A
        - N_B: total number of genes in species B
    """
    N_A = np.sum([chrom_info_A[c]['size'] for c in chrom_info_A.keys()])
    N_B = np.sum([chrom_info_B[c]['size'] for c in chrom_info_B.keys()])
    rho = dot_plot.shape[0] / (N_A * N_B)
    return rho, N_A, N_B

def k_min(alpha, rho, N_A, N_B, conservative = False):
    """
    Calculate the minimum nanosynteny block size, from a percolation theory model.

    Input:
        - alpha: float between 0 and 1, number controlling statistical significance. We use 0.05.
        - rho: float, dot density
        - N_A: total number of genes in species A
        - N_B: total number of genes in species B
        - conservative: boolean, be conservatve in rounding k_min?

    Output:
        - k_min: integer, the minimum significant nanosynteny block size
    """
    N = 2 * N_A * N_B
    k = np.log(1-(1-alpha)**(1/N)) / np.log(rho)
    if conservative:
        return int(np.ceil(k))
    else:
        return int(np.round(k))

def synteny_block_probability_permuted(n, k, rho, N_A, N_B):
    """
    Calculate the probability of finding n nansynteny blocks consisting of k genes.

    Input:
        - n: integer, number of synteny blocks
        - k: integer, synteny block size
        - rho: float, dot density
        - N_A: total number of genes in species A
        - N_B: total number of genes in species B

    Output:
        - P_n: probability of finding n nansynteny blocks consisting of k genes in a dot plot with density rho
          and with N_A and N_B genes in the two genomes.
    """
    N = 2 * N_A * N_B
    return stats.binom.pmf(n, N, rho**k)
