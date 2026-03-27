
import numpy as np
import gzip
from ._tools import * 

def load_dot_plot(dot_plot_file):
    """
    Load a dot plot file.

    Input:
        - dot_plot_file: string, path to a dot_plot file

    Output:
        - spA: string, binomial name of species A in dot plot
        - spB: string, binomial name of species B in dot plot
        - dot_plot: N X 9 array, dot plot array
        - labels: 1-dimensional array with 9 entries, labels for dot plot array columns
    """
    if dot_plot_file.split('.')[-1] == 'gz':
        with gzip.open(dot_plot_file,'rb') as f:
            lines = f.readlines()
            lines = [line.decode() for line in lines]
    else:
        with open(dot_plot_file,'r') as f:
            lines = f.readlines()
    sp = lines[0].rstrip('\n')
    spA, spB = sp.split(',')
    labels = np.array(lines[1].rstrip('\n').split(','))
    dot_plot = []
    for line in lines[2:]:
        dot_plot.append(np.array(line.rstrip('\n').split(',')))
    dot_plot = np.vstack(dot_plot)
    return spA, spB, dot_plot, labels

def load_species_data(species_data_file):
    """
    Load species data file:

    Input:
        - species_data_file: string, path to a species data file

    Output:
        - sp: string, binomial name of species
        - species_data: N X 12, species data array
        - labels: 1-dimensional array with 12 entries, labels for species data array columns
        - chrom_info: dictionary, a chromosome information dictionary
    """
    if species_data_file.split('.')[-1] == 'gz':
        with gzip.open(species_data_file,'rb') as f:
            lines = f.readlines()
            lines = [line.decode() for line in lines]
    else:
        with open(species_data_file,'r') as f:
            lines = f.readlines()
    sp = lines[0].rstrip('\n')
    labels = np.array(lines[1].rstrip('\n').split(','))
    species_data = []
    for line in lines[2:]:
        species_data.append(np.array(line.rstrip('\n').split(',')))
    species_data = np.vstack(species_data)    
    chrom_info = get_chrom_info(species_data)
    return sp, species_data, labels, chrom_info

def load_synteny_blocks(synteny_blocks_file):
    """
    Load a synteny_blocks list.

    Input:
        - synteny_blocks_file: string, path to a synteny blocks file

    Output:
        - spA: string, binomial name of species A in dot plot
        - spB: string, binomial name of species B in dot plot
        - synteny_blocks: list of N X 9 arrays, where N can vary by entry. This contains pairs of genes in synteny.
        - labels: 1-dimensional array with 9 entries, labeling the columns of synteny_blocks entries
    """
    if synteny_blocks_file.split('.')[-1] == 'gz':
        with gzip.open(synteny_blocks_file,'rb') as f:
            lines = f.readlines()
            lines = [line.decode() for line in lines]
    else:
        with open(synteny_blocks_file,'r') as f:
            lines = f.readlines()
    sp = lines[0].rstrip('\n')
    spA, spB = sp.split(',')
    labels = np.array(lines[1].rstrip('\n').split(','))
    synteny_blocks = []
    for line in lines[2:]:
        if line == '#\n':
            try:
                synteny_blocks.append(np.vstack(block))
            except NameError:
                pass
            block = []
        else:
            try:
                block.append(np.array(line.rstrip('\n').split(',')))
            except NameError:
                raise ValueError("Your synteny block file format seems to be wrong.")
    if lines[-1] != '#\n':
        synteny_blocks.append(np.vstack(block))
    return spA, spB, synteny_blocks, labels

def save_dot_plot(dot_plot_file, spA, spB, dot_plot, labels = None, gzip_file = False):
    """
    Save a dot plot array to a .csv file.

    Input:
        - dot_plot_file: string, path to a dot plot file name to save to
        - spA: string, binomial name of species A in dot plot
        - spB: string, binomial name of species B in dot plot
        - dot_plot: N X 9 array, dot plot array
        - labels: 1-dimensional array with 9 entries, labels for dot plot array columns
          Default: None
        - gzip_file: boolean, save compressed dotplot file with gzip?

    Output:
        - Saved dot plot file.
    """
    if labels is None:
        labels = 'chromosome name A,relative index A,chromosome name B,relative index B,empty 1,empty 2,empty 3,empty 4,empty 5\n'
    else:
        if len(labels) == 9:
            labels = ','.join(labels)+'\n'
        else:
            raise ValueError("The labels list/array for dot plots must be length 9.")
    lines = []
    lines.append(','.join([spA,spB])+'\n')
    lines.append(labels)
    for line in dot_plot:
        lines.append(','.join(line)+'\n')
    if gzip_file == True:
        with gzip.open(dot_plot_file+'.gz','wb') as f:
            gzip_lines = [line.encode() for line in lines]
            f.writelines(gzip_lines)
    else:
        with open(dot_plot_file,'w') as f:
            f.writelines(lines)

def save_species_data(species_data_file, sp, species_data, labels = None, gzip_file = False):
    """
    Save a species data array to a .csv file.

    Input:
        - species_data_file: string, path to a species data file name to save to
        - sp: string, binomial name of species
        - species_data: N X 12 array, species data array
        - labels: 1-dimensional array with 12 entries, labels for species data array columns
          Default: None
        - gzip_file: boolean, save compressed dotplot file with gzip?

    Output:
        - Saved species data file.
    """
    if labels is None:
        labels = 'chromosome name,chromosome number,relative index,absolute index,orthogroup,gene name,MAKER gene label,empty 1,empty 2,empty 3,empty 4,empty 5\n'
    else:
        if len(labels) == 12:
            labels = ','.join(labels)+'\n'
        else:
            raise ValueError("The labels list/array for species data must be length 12.")    
    lines = []
    lines.append(sp+'\n')
    lines.append(labels)
    for line in species_data:
        lines.append(','.join(line)+'\n')
    if gzip_file == True:
        with gzip.open(species_data_file+'.gz','wb') as f:
            gzip_lines = [line.encode() for line in lines]
            f.writelines(gzip_lines)
    else:
        with open(species_data_file,'w') as f:
            f.writelines(lines)

def save_synteny_blocks(synteny_blocks_file, spA, spB, synteny_blocks, chrom_info_A, chrom_info_B, labels = None, gzip_file = False):
    """
    Save a synteny blocks list to a .csv file.

    Input:
        - synteny_blocks_file: string, path to a synteny blocks file name to save to
        - spA: string, binomial name of species A in dot plot
        - spB: string, binomial name of species B in dot plot
        - synteny_blocks: list of N X 9 arrays, synteny blocks list
        - chrom_info_A: dictionary, chromosome information dictionary for species A
        - chrom_info_B: dictionary, chromosome information dictionary for species B
        - labels: 1-dimensional array with 9 entries, labels for dot plot array columns
          Default: None
        - gzip_file: boolean, save compressed dotplot file with gzip?

    Output:
        - Saved dot plot file.
    """
    if labels is None:
        labels = 'chromosome name A,relative index A,chromosome name B,relative index B,empty 1,empty 2,empty 3,empty 4,empty 5\n'
    else:
        if len(labels) == 9:
            labels = ','.join(labels)+'\n'
        else:
            raise ValueError("The labels list/array for synteny blocks must be length 9.")    
    try:
        dummy_var = int(synteny_blocks[0][0,0])
        chrom_int = True
    except ValueError:
        chrom_int = False    
    if chrom_int == True:
        namesA = list(chrom_info_A.keys())
        numbersA = [chrom_info_A[chrom]['number'] for chrom in namesA]
        chrom_num_to_nameA = {numbersA[i]:namesA[i] for i in range(len(namesA))}
        namesB = list(chrom_info_B.keys())
        numbersB = [chrom_info_B[chrom]['number'] for chrom in namesB]
        chrom_num_to_nameB = {numbersB[i]:namesB[i] for i in range(len(namesB))}
    lines = []
    lines.append(','.join([spA,spB])+'\n')
    lines.append(labels)
    for block in synteny_blocks:
        lines.append('#\n')
        str_block = block.astype(str)
        if chrom_int == True:
            str_block[:,0] = chrom_num_to_nameA[block[0,0]]
            str_block[:,2] = chrom_num_to_nameB[block[0,2]]
        for line in str_block:
            lines.append(','.join(line)+'\n')
    if gzip_file == True:
        with gzip.open(synteny_blocks_file+'.gz','wb') as f:
            gzip_lines = [line.encode() for line in lines]
            f.writelines(gzip_lines)
    else:
        with open(synteny_blocks_file,'w') as f:
            f.writelines(lines)
