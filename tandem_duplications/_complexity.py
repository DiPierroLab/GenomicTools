
import numpy as np
import networkx as ntx
from GenomicTools.tools import *
from ._shift_maps import *
from ._shift_unshift import *

def entropy(p):
    """
    Calculate the Shannon entropy (in bits) of a PMF.

    Input:
        - p: array, a probability mass function with all positive entries, summing to 1

    Output:
        - H: float, Shannon entropy in bits
    """
    H = - np.nan_to_num(p * np.log2(p)).sum()
    return H

def window_entropy_avg(x, w):
    """
    Sliding window entropy of an array of integers. The input x is a 1-dimensional array (orthogroup indices or tandem 
    duplication indices in our case) of integers, and we examine each (counting overlapping) window of w consecutive elements 
    of x. For each window, we find an empirical PMF, counting the frequency of each present integer in the w window. We then 
    calculate the Shannon entropy, and count it towards each slot in x covered by the window. The output is then the average 
    at each slot of x over all windows that cover that slot.

    Input:
        - x: 1-dimensional array, an array of integers
        - w: integer, window size

    Output:
        - H_average: 1-dimensional array, average sliding window Shannon entropy of x
    """
    L = x.shape[0]
    N = L - w + 1
    H_avg = np.zeros(L)
    coverage = np.zeros(L)
    for n in range(N):
        og, c = np.unique(x[n:(n+w)],return_counts=True)
        H_window = entropy(c / c.sum())
        H_avg[n:(n+w)] += H_window
        coverage[n:(n+w)] += 1
    H_average = H_avg / coverage
    return H_average

def label_tandem(species_data, cc_maps, inv_cc_maps):
    """
    Label genes as being in (or not being in) a tandem duplication block.

    Input:
        - species_data: N X 12 array, species data array
        - cc_map: gene index -> cc index
        - inv_cc_map: cc index -> set of gene indices

    Output:
        - tandem: 1-dimensional array of N elements, binary entries indicating if each gene is in
          a tandem duplication block.
    """
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
    """
    Find the average sliding window entropy for the orthogroups of a genome.

    Input:
        - species_data: N X 12 array, species data array
        - chrom_info: dictionary, chromosome information
        - w: integer, sliding window size

    Output:
        - H: 1-dimensional array of N elements, the average sliding window entropy  
    """
    H = []
    for chrom in alphanum_sort(chrom_info.keys()):
        H += list(window_entropy_avg(species_data[species_data[:,0] == chrom][:,4], w))
    H = np.array(H)
    return H

def find_tandem_entropy(species_data, chrom_info, cc_maps, w):
    """
    Find the average sliding window entropy for the tandem duplication indices of a genome.

    Input:
        - species_data: N X 12 array, species data array
        - chrom_info: dictionary, chromosome information
        - w: integer, sliding window size

    Output:
        - H: 1-dimensional array of N elements, the average sliding window entropy  
    """
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
    """
    Calculate and add complexity information (orthogroup and tandem duplication entropies) to species data.

    Input:
        - sp: string, binomial species name
        - species_data: N X 12 array, species data array
        - chrom_info: dictionary, chromosome information
        - entopy_window: integer, sliding window size (Default = 10)
        - tandem_window: integer, how many genes away can a tandem duplicate be (Default = 1)

    Output:
        - sp: string, binomial species name
        - species_data: N X 12 array, species data array with complexity information added
        - labels: string of 12 comma separated labels for new species data array columns
        - chrom_info: dictionary, chromosome information
    """
    cc_maps, inv_cc_maps, shift_maps, unshift_maps = create_shift_map(species_data, tandem_window)
    tandem = label_tandem(species_data,cc_maps,inv_cc_maps)
    H = find_entropy(species_data,chrom_info,entropy_window)
    tandem_H = find_tandem_entropy(species_data,chrom_info,cc_maps,entropy_window)
    species_data[:,7] = tandem
    species_data[:,8] = H
    species_data[:,9] = tandem_H
    labels = 'chromosome name,chromosome number,relative index,absolute index,orthogroup,gene name,MAKER gene label,tandem dup,H-w%i,tandem H-w%i,empty 4,empty 5\n'%(entropy_window,entropy_window)
    return sp, species_data, labels, chrom_info
