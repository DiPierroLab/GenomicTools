import numpy as np
from scipy.sparse import csr_matrix
from GenomicTools.tools import *

def assemble_homology_matrix(spA, spB, species_data_A, species_data_B, chrom_info_A, chrom_info_B):
    """
    This function creates a homology matrix, forming a comparison between two distinct species, or one species and itself. Genes with an unassigned orthogroup should have an orthogroup index of 0. The function explicitly constructs a homology matrix where all genes with orthogroup 0 have no homology relations except when this is a self-comparison and we are comparing a gene with itself (in the strictest sense).

    Input: 
    - spA: string, name of species A ('Genus_species')
    - spB: string, name of species B ('Genus_species')
    - species_data_A: numpy array (N X 12), species data for species A
    - species_data_B: numpy array (N X 12), species data for species B
    - chrom_info_A: dictionary, chromosome information for species A
    - chrom_info_B: dictionary, chromosome information for species B

    Output:
    - homology_matrix: dictionary of scipy csr sparse matrices, the keys of the dictionary are tuples of chromosome names (strings), with chromosomes from the first species first. For example, if we are comparing human and mouse, the keyword ('chr1','chr3') returns the homology matrix for human chromosome 1 and mouse chromosome 3 (using Gencode chromosome names for these two species).
    """
    chromsA = alphanum_sort(list(chrom_info_A.keys()))
    chromsB = alphanum_sort(list(chrom_info_B.keys()))
    homology_matrix = {}
    for chromA in chromsA:
        for chromB in chromsB:
            genesA = species_data_A[species_data_A[:,0] == chromA]
            genesB = species_data_B[species_data_B[:,0] == chromB]
            ngenesA = genesA.shape[0]
            ngenesB = genesB.shape[0]
            MA = np.vstack(ngenesB * [genesA[:,4:5].T])
            MB = np.hstack(ngenesA * [genesB[:,4:5]])
            zerosA = 1 - (MA == 0)
            zerosB = 1 - (MB == 0)
            M = (MA == MB).astype(int) * zerosA * zerosB
            if (spA == spB) and (chromA == chromB):
                np.fill_diagonal(M,1)
            homology_matrix[(chromA,chromB)] = csr_matrix(M)
    return homology_matrix
