
import numpy as np
from scipy.sparse import csr_matrix
from ._tools import *

def convert_dot_plot_to_homology_matrix(dot_plot, chrom_info_A, chrom_info_B): 
    chromsA = alphanum_sort(np.unique(dot_plot[:,0]))
    chromsB = alphanum_sort(np.unique(dot_plot[:,2]))
    homology_matrix = {}
    for chromA in chromsA:
        for chromB in chromsB:
            dot_plot_chromAB = dot_plot[(dot_plot[:,0] == chromA)*(dot_plot[:,2] == chromB)]
            M = np.zeros((chrom_info_A[chromA]['size'],chrom_info_B[chromB]['size']))
            M[dot_plot_chromAB[:,1].astype(int),dot_plot_chromAB[:,3].astype(int)] = 1
            homology_matrix[(chromA,chromB)] = csr_matrix(M)
    return homology_matrix

def convert_homology_matrix_to_dot_plot(homology_matrix):
    dot_plot = []
    for chromAB in homology_matrix.keys():
        chromA, chromB = chromAB
        dotsA, dotsB = np.where(homology_matrix[chromAB].toarray() == 1)
        N = dotsA.shape[0]
        dot_plot_chromAB = np.vstack([np.array(N*[chromA]),dotsA+1,np.array(N*[chromB]),dotsB+1]).T
        dot_plot.append(dot_plot_chromAB)
    dot_plot = np.vstack(dot_plot)
    return dot_plot
