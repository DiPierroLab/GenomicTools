
import numpy as np
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def create_species_data(sp, annotation_file, orthogroups_file, calculate_duplications = True, entropy_window = 10, tandem_window = 1):
    """
    Create a species data array from a .gff format annotation. 

    Input:
    - sp: string, species name
    - annotation_file: string, name of .gff annotation file to load
    - orthogroups_file: string, name of .tsv orthogroups file from OrthoFinder
    - calculate_duplications: boolean, add tandem duplication information and calculate entropies?
    - entropy_window: int, size of sliding window (in number of genes) for calculating entropy
    - tandem_window: int, size of window for identifying tandem duplications
    
    Output:
    - sp: string, species name
    - sp_data: N X 12 numpy array, where N is the number of genes
    - labels: string, labels for species data columns separated by commas
    - chrom_info: dictionary, chromosome data
    """
    with open(orthogroups_file,"r") as f:
        og = f.readlines()
    bins = np.cumsum([0] + [len(i) for i in og[1:]])[:-1]
    og_str = ''.join(og[1:])
    annotation = np.loadtxt(annotation_file,dtype=str,delimiter='\t')
    chroms = alphanum_sort(np.unique(annotation[:,1]))
    sp_data_list = []
    n_tot = 0
    for nc, chrom in enumerate(chroms):
        annotation_chrom = annotation[annotation[:,1] == chrom]
        order = np.argsort(annotation_chrom[:,2].astype(int))
        annotation_chrom = annotation_chrom[order]        
        N = annotation_chrom.shape[0]
        chrom_name = np.array(N*[chrom])
        chrom_num = np.array(N*[nc+1])
        relative_index = np.arange(1,N+1)
        absolute_index = n_tot + relative_index
        gene_name = annotation_chrom[:,0]
        gene_label = annotation_chrom[:,-1]
        orthogroup = []
        for g in gene_name:
            g_og = np.digitize(og_str.find(g),bins)
            orthogroup.append(g_og)
        orthogroup = np.array(orthogroup)
        empty = np.array(N*[''])
        sp_data_chrom = np.vstack([chrom_name,chrom_num,relative_index,absolute_index,orthogroup,gene_name,gene_label,empty,empty,empty,empty,empty])
        sp_data_list.append(sp_data_chrom.T)
        n_tot += N
    sp_data = np.vstack(sp_data_list)
    chrom_info = get_chrom_info(sp_data)

    if calculate_duplications == False:
        labels = 'chromosome name,chromosome number,relative index,absolute index,orthogroup,gene name,MAKER gene label,empty 1,empty 2,empty 3,empty 4,empty 5\n'
    elif calculate_duplications == True:
        cc_maps, inv_cc_maps, shift_maps, unshift_maps = create_shift_map(sp_data, tandem_window)
        tandem = label_tandem(sp_data,cc_maps,inv_cc_maps)
        H = find_entropy(sp_data,chrom_info,entropy_window)
        tandem_H = find_tandem_entropy(sp_data,chrom_info,cc_maps,entropy_window)
        sp_data[:,7] = tandem
        sp_data[:,8] = H
        sp_data[:,9] = tandem_H
        labels = 'chromosome name,chromosome number,relative index,absolute index,orthogroup,gene name,MAKER gene label,tandem dup,H-w%i,tandem H-w%i,empty 3,empty 4,empty 5\n'%(entropy_window,entropy_window)    

    return sp, sp_data, labels, chrom_info
