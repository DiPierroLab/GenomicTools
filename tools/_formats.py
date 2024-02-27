
import numpy as np
from scipy.sparse import csr_matrix
from ._tools import *
from ._io import *
import gzip
import pickle as pkl

def convert_dot_plot_to_homology_matrix(dot_plot, chrom_info_A, chrom_info_B): 
    chromsA = alphanum_sort(np.unique(dot_plot[:,0]))
    chromsB = alphanum_sort(np.unique(dot_plot[:,2]))
    homology_matrix = {}
    for chromA in chromsA:
        for chromB in chromsB:
            dot_plot_chromAB = dot_plot[(dot_plot[:,0] == chromA)*(dot_plot[:,2] == chromB)]
            M = np.zeros((chrom_info_A[chromA]['size'],chrom_info_B[chromB]['size']))
            M[dot_plot_chromAB[:,1].astype(int)-1,dot_plot_chromAB[:,3].astype(int)-1] = 1
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

def convert_save_old_dot_plot_to_new(dot_plot_file):
    spA, spB = dot_plot_file.split('/')[-1].split('-')[0:2]
    with open(dot_plot_file,"rb") as f:
        hmat = pkl.load(f)
    hmat_lines = []
    for chrom_pair in hmat['data'].keys():
        x,y = np.where(hmat['data'][chrom_pair]['homology_matrix'].A.T == 1)
        chromA,chromB = chrom_pair
        out = np.vstack([np.array(x.shape[0]*[chromA]),x+1,np.array(x.shape[0]*[chromB]),y+1]).T
        hmat_lines.append(out)
    hmat_lines = np.vstack(hmat_lines)
    empty = np.vstack((9-hmat_lines.shape[1])*[np.array(hmat_lines.shape[0]*[''])]).T
    hmat_lines = np.hstack([hmat_lines,empty])
    save_dot_plot(spA+'-'+spB+'-dotplot.csv',spA,spB,hmat_lines,gzip_file=True)
    
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
