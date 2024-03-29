
import numpy as np
from scipy import stats
from GenomicTools.tools import *

def contingency_dot_plot(chromA, chromB, dots_input, chrom_info_A, chrom_info_B, microsynteny_input = True):
    if microsynteny_input == True:
        chromA = chrom_info_A[chromA]['number']
        chromB = chrom_info_B[chromB]['number']
        dot_plot = np.vstack(dots_input)
    else:
        dot_plot = np.copy(dots_input)
    a = dot_plot[(dot_plot[:,0]==chromA)*(dot_plot[:,2]==chromB)].shape[0]
    b = dot_plot[(dot_plot[:,0]==chromA)].shape[0]
    c = dot_plot[(dot_plot[:,2]==chromB)].shape[0]
    d = dot_plot[(dot_plot[:,0]!=chromA)*(dot_plot[:,2]!=chromB)].shape[0]
    C = np.array([[a,b],[c,d]])
    return C

def calculate_macrosynteny(alpha, dots_input, chrom_info_A, chrom_info_B, microsynteny_input = True, calculate_metric = True):
    chrom_names_A = alphanum_sort(chrom_info_A.keys())
    chrom_names_B = alphanum_sort(chrom_info_B.keys())
    chrom_locs_A = np.cumsum([0] + [chrom_info_A[key]['size'] for key in chrom_names_A])
    chrom_locs_B = np.cumsum([0] + [chrom_info_B[key]['size'] for key in chrom_names_B])
    chrom_edges_A = {c:chrom_locs_A[n:(n+2)] for n, c in enumerate(chrom_names_A)}
    chrom_edges_B = {c:chrom_locs_B[n:(n+2)] for n, c in enumerate(chrom_names_B)}
    macrosynteny_chrom_pairs = []
    H_mat_macro = np.zeros([len(chrom_info_A.keys()),len(chrom_info_B.keys())])
    H_mat_macro_max = np.zeros([len(chrom_info_A.keys()),len(chrom_info_B.keys())])
    for nA, chromA in enumerate(chrom_names_A):
        for nB, chromB in enumerate(chrom_names_B):
            C_mat = contingency_dot_plot(chromA,chromB,dots_input,chrom_info_A,chrom_info_B,microsynteny_input=microsynteny_input)
            FET = stats.fisher_exact(C_mat,alternative='greater')
            if FET.pvalue < alpha:
                macrosynteny_chrom_pairs.append([chromA,chromB])
                H_mat_macro[nA,nB] = C_mat[0,0]
            H_mat_macro_max[nA,nB] = np.diff(chrom_edges_A[chromA])[0] * np.diff(chrom_edges_B[chromB])[0]
    if calculate_metric == True:
        H_mat_macro_null = np.zeros([len(chrom_info_A.keys()),len(chrom_info_B.keys())])
        np.fill_diagonal(H_mat_macro_null,H_mat_macro.sum(np.argmax(H_mat_macro.shape)))
        p_null = H_mat_macro_null / H_mat_macro_null.sum()   
        p_actual = H_mat_macro / H_mat_macro.sum()
        p_max = H_mat_macro_max / np.sum(H_mat_macro_max)
        H_null = - np.sum(np.nan_to_num(p_null * np.log(p_null)))
        H_actual = - np.sum(np.nan_to_num(p_actual * np.log(p_actual)))
        H_max = - np.sum(np.nan_to_num(p_max * np.log(p_max)))
        metric = (H_actual - H_null) / (H_max - H_null)
        return macrosynteny_chrom_pairs, H_mat_macro, metric
    else:
        return macrosynteny_chrom_pairs, H_mat_macro
