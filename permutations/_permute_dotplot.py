import numpy as np
from scipy import sparse
from GenomicTools.tools import *
from GenomicTools.dotplots import *

def permute_dotplot(dot_plot_result, og_data_A, og_data_B, perm_A, perm_B):
    chroms = np.vstack(list(dot_plot_result['data'].keys()))
    chromsA = tools.alphanum_sort(np.unique(chroms[:,0]))
    chromsB = tools.alphanum_sort(np.unique(chroms[:,1]))
    
    chrom_locsA = np.cumsum([0] + [dot_plot_result['data'][(c,'chr1')]['homology_matrix'].shape[1] for c in chromsA])
    chrom_locsB = np.cumsum([0] + [dot_plot_result['data'][('chr1',c)]['homology_matrix'].shape[0] for c in chromsB])
    chromA_indices = {c:chrom_locsA[i:(i+2)] for i,c in enumerate(chromsA)}
    chromB_indices = {c:chrom_locsB[i:(i+2)] for i,c in enumerate(chromsB)}    

    dot_matrix = []
    ogs_A_all = []
    gene_names_A_all = []
    gene_labels_A_all = []
    ogs_B_all = []
    gene_names_B_all = []
    gene_labels_B_all = []
    for n, chromA in enumerate(chromsA):
        dot_matrix_col = []
        ogs_A_all.append(og_data_A[chromA]['ogs'])
        gene_names_A_all.append(og_data_A[chromA]['gene_names'])
        gene_labels_A_all.append(og_data_A[chromA]['gene_labels'])
        for chromB in chromsB:
            dot_matrix_col.append(dot_plot_result['data'][(chromA,chromB)]['homology_matrix'])
            if n == 0:
                ogs_B_all.append(og_data_B[chromB]['ogs'])
                gene_names_B_all.append(og_data_B[chromB]['gene_names'])
                gene_labels_B_all.append(og_data_B[chromB]['gene_labels'])
        dot_matrix.append(sparse.vstack(dot_matrix_col)) 
    dot_matrix = sparse.csr_matrix(sparse.hstack(dot_matrix).T) 
    ogs_A_all = np.hstack(ogs_A_all)
    gene_names_A_all = np.hstack(gene_names_A_all)
    gene_labels_A_all = np.hstack(gene_labels_A_all)
    ogs_B_all = np.hstack(ogs_B_all)
    gene_names_B_all = np.hstack(gene_names_B_all)
    gene_labels_B_all = np.hstack(gene_labels_B_all)

    permuted_dot_matrix = dot_matrix[perm_A,:][:,perm_B]
    permuted_ogs_A_all = ogs_A_all[:,perm_A]
    permuted_gene_names_A_all = gene_names_A_all[perm_A]
    permuted_gene_labels_A_all = gene_labels_A_all[perm_A]
    permuted_ogs_B_all = ogs_B_all[:,perm_B]
    permuted_gene_names_B_all = gene_names_B_all[perm_B]
    permuted_gene_labels_B_all = gene_labels_B_all[perm_B]

    permuted_dot_plot_result = {}
    for key in dot_plot_result:
        permuted_dot_plot_result[key] = dot_plot_result[key]
    permuted_dot_plot_result['perm1'] = perm_A
    permuted_dot_plot_result['perm2'] = perm_B

    for chromA in chromsA:
        for chromB in chromsB:
            a1, a2 = chromA_indices[chromA]
            b1, b2 = chromB_indices[chromB]
            permuted_dot_plot_result['data'][(chromA,chromB)]['homology_matrix'] = permuted_dot_matrix[a1:a2,:][:,b1:b2].T

    permuted_og_data_A = {}
    for key in og_data_A.keys():
        permuted_og_data_A[key] = og_data_A[key]
    permuted_og_data_B = {}
    for key in og_data_B.keys():
        permuted_og_data_B[key] = og_data_B[key]
    
    for chromA in chromsA:
        a1, a2 = chromA_indices[chromA]
        permuted_og_data_A[chromA]['ogs'] = permuted_ogs_A_all[:,a1:a2]
        permuted_og_data_A[chromA]['gene_names'] = permuted_gene_names_A_all[a1:a2]
        permuted_og_data_A[chromA]['gene_labels'] = permuted_gene_labels_A_all[a1:a2]
    for chromB in chromsB:
        b1, b2 = chromB_indices[chromB]
        permuted_og_data_B[chromB]['ogs'] = permuted_ogs_B_all[:,b1:b2]
        permuted_og_data_B[chromB]['gene_names'] = permuted_gene_names_B_all[b1:b2]
        permuted_og_data_B[chromB]['gene_labels'] = permuted_gene_labels_B_all[b1:b2]

    return permuted_dot_plot_result, permuted_og_data_A, permuted_og_data_B
