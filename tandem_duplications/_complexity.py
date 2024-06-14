
import numpy as np
import networkx as ntx
from GenomicTools.tools import *
from ._shift_maps import *
from ._shift_unshift import *

def entropy(p):
    H = - np.nan_to_num(p * np.log2(p)).sum()
    return H

def window_entropy_avg(x, w):
    L = x.shape[0]
    N = L - w + 1
    H_avg = np.zeros(L)
    coverage = np.zeros(L)
    for n in range(N):
        og, c = np.unique(x[n:(n+w)],return_counts=True)
        H_window = entropy(c / c.sum())
        H_avg[n:(n+w)] += H_window
        coverage[n:(n+w)] += 1
    return H_avg / coverage

def label_tandem(species_data, cc_maps, inv_cc_maps):
    tandem = []
    for gene in species_data:
        chrom = gene[0]
        index = int(gene[2])
        tandem_size = len(inv_cc_maps[chrom][cc_maps[chrom][index]])
        if tandem_size > 1:
            tandem.append(1)
        else:
            tandem.append(0)
    tandem = np.array(tandem)
    return tandem

def find_entropy(species_data, chrom_info, w):
    H = []
    for chrom in alphanum_sort(chrom_info.keys()):
        H += list(window_entropy_avg(species_data[species_data[:,0] == chrom][:,4], w))
    H = np.array(H)
    return H

def find_tandem_entropy(species_data, chrom_info, cc_maps, w):
    tandem_H = []
    for chrom in alphanum_sort(chrom_info.keys()):
        tandem_cc = []
        for gene in species_data[species_data[:,0] == chrom]:
            index = int(gene[2])
            tandem_cc.append(cc_maps[chrom][index])
        tandem_H += list(window_entropy_avg(np.array(tandem_cc),w))  
    tandem_H = np.array(tandem_H)
    return tandem_H

def add_complexity_info_to_species_data(sp, species_data, chrom_info, entropy_window = 10, tandem_window = 1):
    cc_maps, inv_cc_maps, shift_maps, unshift_maps = create_shift_map(species_data, tandem_window)
    tandem = label_tandem(species_data,cc_maps,inv_cc_maps)
    H = find_entropy(species_data,chrom_info,entropy_window)
    tandem_H = find_tandem_entropy(species_data,chrom_info,cc_maps,entropy_window)
    species_data[:,7] = tandem
    species_data[:,8] = H
    species_data[:,9] = tandem_H
    labels = 'chromosome name,chromosome number,relative index,absolute index,orthogroup,gene name,MAKER gene label,tandem dup,H-w%i,tandem H-w%i,empty 3,empty 4,empty 5\n'%(entropy_window,entropy_window)
    return sp, species_data, labels, chrom_info
