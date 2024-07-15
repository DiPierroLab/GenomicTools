
import numpy as np
from scipy.sparse import csr_matrix
from ._tools import *
from ._io import *
import gzip
import pickle as pkl

def convert_dot_plot_to_homology_matrix(dot_plot, chrom_info_A, chrom_info_B):
    chromsA = alphanum_sort(list(chrom_info_A.keys()))
    chromsB = alphanum_sort(list(chrom_info_B.keys()))
    homology_matrix = {}
    for chromA in chromsA:
        for chromB in chromsB:
            dot_plot_chromAB = dot_plot[(dot_plot[:,0] == chromA)*(dot_plot[:,2] == chromB)]
            M = np.zeros((chrom_info_A[chromA]['size'],chrom_info_B[chromB]['size']))
            M[dot_plot_chromAB[:,1].astype(int)-1,dot_plot_chromAB[:,3].astype(int)-1] = 1
            homology_matrix[(chromA,chromB)] = csr_matrix(M)
    return homology_matrix

def convert_homology_matrix_to_dot_plot(homology_matrix, labels = None, additional_data = None):
    dot_plot = []
    for chromAB in homology_matrix.keys():
        chromA, chromB = chromAB
        dotsA, dotsB = np.where(homology_matrix[chromAB].toarray() == 1)
        N = dotsA.shape[0]
        dot_plot_chromAB = np.vstack([np.array(N*[chromA]),dotsA+1,np.array(N*[chromB]),dotsB+1]).T
        dot_plot.append(dot_plot_chromAB)
    dot_plot = np.vstack(dot_plot)
    empty = np.array(dot_plot.shape[0]*[['']])
    if additional_data is None:
        dot_plot = np.hstack([dot_plot]+5*[empty])
    else:
        n_additional_data = len(additional_data)
        to_add = []
        for i in range(n_additional_data):
            to_add.append(additional_data[i].reshape(dot_plot.shape[0],1))
        to_add += (5-i-1)*[empty]
        dot_plot = np.hstack([dot_plot]+to_add)
    if labels is None:
        if additional_data is not None:
            raise ValueError("Label your additional data!")
        else:
            labels = np.array(['chromosome name A','relative index A','chromosome name B','relative index B','empty 1','empty 2','empty 3','empty 4','empty 5'])
    return dot_plot, labels
    
def convert_save_old_species_data_to_new(species_data_file):
    sp = species_data_file.split('/')[-1].split('-')[0]
    with open(species_data_file,"rb") as f:
        spdat = pkl.load(f)
    spdat_list = []
    n_tot = 0
    for nc, chrom in enumerate(alphanum_sort(list(spdat.keys()))):
        N = spdat[chrom]['gene_names'].shape[0]
        chrom_name = np.array(N*[chrom])
        chrom_num = np.array(N*[nc+1])
        relative_index = np.arange(1,N+1)
        absolute_index = n_tot + relative_index
        orthogroup = spdat[chrom]['ogs']
        gene_name = spdat[chrom]['gene_names']
        gene_label = spdat[chrom]['gene_labels']
        empty = np.array(N*[''])
        sppdat_chrom = np.vstack([chrom_name,chrom_num,relative_index,absolute_index,orthogroup,gene_name,gene_label,empty,empty,empty,empty,empty])
        spdat_list.append(sppdat_chrom.T)
        n_tot += N    
    spdat_arr = np.vstack(spdat_list)
    save_species_data(sp+'_species_data.csv',sp,spdat_arr,gzip_file=True)
