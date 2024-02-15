
import numpy as np
from scipy import stats

def calculate_density(dot_plot, chrom_info_A, chrom_info_B):
    N_A = np.sum([chrom_info_A[c]['size'] for c in chrom_info_A.keys()])
    N_B = np.sum([chrom_info_B[c]['size'] for c in chrom_info_B.keys()])
    rho = dot_plot.shape[0] / (N_A * N_B)
    return rho, N_A, N_B

def k_min(alpha, rho, N_A, N_B, conservative = False):
    N = 2 * N_A * N_B
    k = np.log(1-(1-alpha)**(1/N)) / np.log(rho)
    if conservative:
        return int(np.ceil(k))
    else:
        return int(np.round(k))

def synteny_block_probability_permuted(n, k, rho, N_A, N_B):
    N = 2 * N_A * N_B
    return stats.binom.pmf(n, N, rho**k)
