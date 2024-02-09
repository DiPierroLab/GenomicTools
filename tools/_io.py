
import numpy as np
import gzip
from ._tools import * 

def get_chrom_info(species_data):
    chrom_info = {}
    chroms = alphanum_sort(np.unique(data[:,0]))
    for n, chrom in enumerate(chroms):
        chrom_info[chrom] = {}
        chrom_info[chrom]['number'] = n + 1
        chrom_info[chrom]['size'] = species_data[species_data[:,0] == chrom].shape[0]
    return chrom_info

def load_dot_plot(dot_plot_file):
    d = np.loadtxt(dot_plot_file,skiprows=1,dtype=str,delimiter=',')
    labels = d[0]
    dot_plot = d[1:]
    with open(dot_plot_file,'r') as f:
        sp = f.readlines()[0].rstrip('\n')
    spA, spB = sp.split(',')
    return spA, spB, dot_plot, labels

def load_species_data(species_data_file):
    d = np.loadtxt(species_data_file,skiprows=1,dtype=str,delimiter=',')
    labels = d[0]
    species_data = d[1:]
    with open(dot_plot_file,'r') as f:
        sp = f.readlines()[0].rstrip('\n')
    chrom_info = get_chrom_info(species_data)
    return sp, species_data, labels, chrom_info

def load_synteny_blocks(synteny_blocks_file):
    with open(synteny_blocks_file,'r') as f:
        lines = f.readlines()
    spA, spB = lines[0].rstrip('\n').split(',')
    labels = np.array(lines[1].strip('\n').split(','))
    synteny_blocks = []
    for line in lines[2:]:
        if line == '#':
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
    return spA, spB, synteny_blocks, labels

def save_dot_plot(dot_plot_file, spA, spB, dot_plot, labels = None, gzip_file = False):
    if labels == None:
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
            f.writelines(lines)
    else:
        with open(dot_plot_file,'w') as f:
            f.writelines(lines)

def save_species_data(species_data_file, sp, species_data, labels = None, gzip_file = False):
    if labels == None:
        labels = 'chromosome name,chromosome number,relative index,absolute index,orthogroup,gene name,MAKER gene label,empty 1,empty 2,empty 3,empty 4,empty 5' 
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
            f.writelines(lines)
    else:
        with open(species_data_file,'w') as f:
            f.writelines(lines)

def save_synteny_blocks(synteny_blocks_file, spA, spB, synteny_blocks, labels = None, gzip_file = False):
    if labels == None:
        labels = 'chromosome name A,relative index A,chromosome name B,relative index B,empty 1,empty 2,empty 3,empty 4,empty 5\n'
    else:
        if len(labels) == 9:
            labels = ','.join(labels)+'\n'
        else:
            raise ValueError("The labels list/array for synteny blocks must be length 9.")    
    lines = []
    lines.append(','.join([spA,spB])+'\n')
    lines.append(labels)
    for block in synteny_blocks:
        lines.append('#')
        for line in block:
            lines.append(','.join(line)+'\n')
    if gzip_file == True:
        with gzip.open(synteny_blocks_file+'.gz','wb') as f:
            f.writelines(lines)
    else:
        with open(synteny_blocks_file,'w') as f:
            f.writelines(lines)
