
import numpy as np
import networkx as ntx
from GenomicTools.tools import *

def create_shift_map(species_data, windowsize):    
    """
    Input:
    - species data
    - windowsize: window for joining genes (e.g. 1 means two genes in the same OG will be joined 
                  only if they have indices i and i+1).
    
    Output:
    - cc_map: gene index -> cc index
    - inv_cc_map: cc index -> set of gene indices
    - shift_map: cc index -> compressed gene index
    - unshift_map: compressed gene index -> cc index
    """
    cc_maps = {}
    inv_cc_maps = {}
    shift_maps = {}
    unshift_maps = {}
    chroms = alphanum_sort(np.unique(species_data[:,0]))
    for chrom in chroms:
        data_chrom = species_data[species_data[:,0] == chrom]
        n_dots = data_chrom.shape[0]
        ogs = data_chrom[:,4:5].astype(int)
        indices = data_chrom[:,2:3].astype(int)
        M_ogs = ogs * np.ones(n_dots) - ogs.T * np.ones(n_dots)
        M_ind = indices * np.ones(n_dots) - indices.T * np.ones(n_dots)
        A = ((M_ogs == 0) * (np.abs(M_ind) <= windowsize)).astype(int)
        G = ntx.Graph(A)
        CC = ntx.connected_components(G)
        cc_map = {}
        shift_map = {}
        inv_cc_map = {}
        unshift_map = {}
        multiplicity_map = {}
        for n, cc in enumerate(CC):
            inv_cc_map[n] = np.array(list(cc))
            for rep in cc:
                cc_map[rep+1] = n
                shift_map[n] = indices.flatten()[n]
                unshift_map[indices.flatten()[n]] = n
        cc_maps[chrom] = cc_map
        inv_cc_maps[chrom] = inv_cc_map
        shift_maps[chrom] = shift_map
        unshift_maps[chrom] = unshift_map
    return [cc_maps, inv_cc_maps, shift_maps, unshift_maps]
