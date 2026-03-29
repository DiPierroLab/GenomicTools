
import numpy as np
from ._identify_blocks import *
from ._filter_blocks import *
from ._block_tools import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *
from GenomicTools.permutations import *
from multiprocessing import Pool
from functools import partial

def find_nanosynteny(condensed_dots, species_data_A, species_data_B, chrom_info_A, chrom_info_B, maps_A, maps_B, nanosynteny_minsize, check_for_nanosynteny_support = True, verbose = False):
    """
    Primary function for finding nanosynteny in a comparison of two genomes. The main input is the condensed (i.e. TDECs instead of
    genes) dot plot.

    Input:
        - condensed_dots: N X 9 array, array of condensed dot plot
        - species_data_A: N X 12 array, species data for species A
        - species_data_B: N X 12 array, species data for species B
        - chrom_info_A: dictionary, chromosome information dictionary for species A
        - chrom_info_B: dictionary, chromosome information dictionary for species B
        - maps_A: list of 4 dictionaries, shift/unshift maps for species A
        - maps_B: list of 4 dictionaries, shift/unshift maps for species B
        - nanosynteny_minsize: integer, minimum number of genes in a nanosynteny block (usually 3)
        - check_for_nanosynteny_support: boolean, double check that each block contains nanosynteny (Default = True)
        - verbose: boolean, be noisy, printing out progress? (Default = False)

    Output:
        - blocks: list of N X 4 arrays, nanosynteny blocks, N varies from block to block
    """
    chromsA = np.unique(condensed_dots[:,0])
    chromsB = np.unique(condensed_dots[:,2])
    blocks = []
    N = chromsA.shape[0] * chromsB.shape[0]
    n = 0
    for chromA in chromsA:
        for chromB in chromsB:
            condensed_dots_AB = condensed_dots[(condensed_dots[:,0] == chromA) * (condensed_dots[:,2] == chromB)]            
            if condensed_dots_AB.shape[0] > 0:
                blocks += find_nanosynteny_chromosome_pair(condensed_dots_AB, species_data_A, species_data_B, chrom_info_A, chrom_info_B, maps_A, maps_B, nanosynteny_minsize, check_for_nanosynteny_support)
            n += 1
            if verbose:
                print('%i / %i'%(n,N),end='\r',flush=True)
    return blocks

def find_nanosynteny_parallel(condensed_dots, species_data_A, species_data_B, chrom_info_A, chrom_info_B, maps_A, maps_B, nanosynteny_minsize, check_for_nanosynteny_support = True, n_proc = 1):
    """
    Use multiprocessing to find nanosynteny in a comparison of two genomes. The main input is the condensed (i.e. TDECs instead of
    genes) dot plot.

    Input:
        - condensed_dots: N X 9 array, array of condensed dot plot
        - species_data_A: N X 12 array, species data for species A
        - species_data_B: N X 12 array, species data for species B
        - chrom_info_A: dictionary, chromosome information dictionary for species A
        - chrom_info_B: dictionary, chromosome information dictionary for species B
        - maps_A: list of 4 dictionaries, shift/unshift maps for species A
        - maps_B: list of 4 dictionaries, shift/unshift maps for species B
        - nanosynteny_minsize: integer, minimum number of genes in a nanosynteny block (usually 3)
        - check_for_nanosynteny_support: boolean, double check that each block contains nanosynteny (Default = True)
        - n_proc = integer, number of processes (Default = 1)

    Output:
        - blocks: list of N X 4 arrays, nanosynteny blocks, N varies from block to block
    """
    chromsA = np.unique(condensed_dots[:,0])
    chromsB = np.unique(condensed_dots[:,2])
    condensed_dots_all_chrom_pairs = []
    for chromA in chromsA:
        for chromB in chromsB:
            condensed_dots_AB = condensed_dots[(condensed_dots[:,0] == chromA) * (condensed_dots[:,2] == chromB)]            
            if condensed_dots_AB.shape[0] > 0:
                condensed_dots_all_chrom_pairs.append(condensed_dots_AB)
    n_pairs = len(condensed_dots_all_chrom_pairs)

    with Pool(processes=n_proc) as pool:
        find_blocks = partial(find_nanosynteny_chromosome_pair,species_data_A=species_data_A,species_data_B=species_data_B,chrom_info_A=chrom_info_A,chrom_info_B=chrom_info_B,maps_A=maps_A,maps_B=maps_B,nanosynteny_minsize=nanosynteny_minsize,check_for_nanosynteny_support=check_for_nanosynteny_support)
        blocks = pool.imap(find_blocks,condensed_dots_all_chrom_pairs,chunksize=round(n_pairs/n_proc)+1)
        pool.close()
        pool.join()

    return blocks

def find_microsynteny(condensed_dots, nanosynteny_blocks, chrom_info_A, chrom_info_B, max_distance, distance_cutoff, nanosynteny_minsize, verbose = False):
    """
    Primary function for finding microsynteny in a comparison of two genomes. The main input is the condensed (i.e. TDECs instead of
    genes) dot plot.

    Input:
        - condensed_dots: N X 9 array, array of condensed dot plot
        - nanosynteny_blocks: list of N X 4 arrays, nanosynteny blocks, N varies from block to block
        - chrom_info_A: dictionary, chromosome information dictionary for species A
        - chrom_info_B: dictionary, chromosome information dictionary for species B
        - max_distance: integer, the maximum distance two dots can be from each other and still be incorporated into microsynteny
          as part of a path between nanosynteny blocks.
        - distance_cutoff: integer, the maximum distance (in genes) that two nanosyteny blocks can be from each other and be 
          connected into microsynteny.
        - nanosynteny_minsize: integer, minimum number of genes in a nanosynteny block (usually 3)
        - verbose: boolean, be noisy, printing out progress? (Default = False)

    Output:
        - blocks: list of N X 4 arrays, microsynteny blocks, N varies from block to block
    """
    chromsA = np.unique(condensed_dots[:,0])
    chromsB = np.unique(condensed_dots[:,2])
    blocks = []
    N = chromsA.shape[0] * chromsB.shape[0]
    n = 0
    for chromA in chromsA:
        for chromB in chromsB:
            condensed_dots_AB = condensed_dots[(condensed_dots[:,0] == chromA) * (condensed_dots[:,2] == chromB)]            
            if condensed_dots_AB.shape[0] > 0:
                nano_AB = [block for block in nanosynteny_blocks if (block[0,0] == chromA) * (block[0,2] == chromB)]
                if len(nano_AB) > 0:
                    blocks += find_microsynteny_chromosome_pair(condensed_dots_AB, nano_AB, max_distance, distance_cutoff, nanosynteny_minsize, chrom_info_A, chrom_info_B)
            n += 1
            if verbose:
                print('%i / %i'%(n,N),end='\r',flush=True)
    return blocks

def find_microsynteny_parallel(condensed_dots, nanosynteny_blocks, chrom_info_A, chrom_info_B, max_distance, distance_cutoff, nanosynteny_minsize, n_proc = 1):
    """
    Primary function for finding microsynteny in a comparison of two genomes. The main input is the condensed (i.e. TDECs instead of
    genes) dot plot.

    Input:
        - condensed_dots: N X 9 array, array of condensed dot plot
        - nanosynteny_blocks: list of N X 4 arrays, nanosynteny blocks, N varies from block to block
        - chrom_info_A: dictionary, chromosome information dictionary for species A
        - chrom_info_B: dictionary, chromosome information dictionary for species B
        - max_distance: integer, the maximum distance two dots can be from each other and still be incorporated into microsynteny
          as part of a path between nanosynteny blocks.
        - distance_cutoff: integer, the maximum distance (in genes) that two nanosyteny blocks can be from each other and be 
          connected into microsynteny.
        - nanosynteny_minsize: integer, minimum number of genes in a nanosynteny block (usually 3)
        - n_proc = integer, number of processes (Default = 1)

    Output:
        - blocks: list of N X 4 arrays, microsynteny blocks, N varies from block to block
    """
    chromsA = np.unique(condensed_dots[:,0])
    chromsB = np.unique(condensed_dots[:,2])
    condensed_dots_nanosynteny_all_chrom_pairs = []
    for chromA in chromsA:
        for chromB in chromsB:
            condensed_dots_AB = condensed_dots[(condensed_dots[:,0] == chromA) * (condensed_dots[:,2] == chromB)]
            if condensed_dots_AB.shape[0] > 0:
                nano_AB = [block for block in nanosynteny_blocks if (block[0,0] == chromA) * (block[0,2] == chromB)]
                if len(nano_AB) > 0:
                    condensed_dots_nanosynteny_all_chrom_pairs.append([condensed_dots_AB, nano_AB])                
    n_pairs = len(condensed_dots_nanosynteny_all_chrom_pairs)
    microsynteny_chrom_pair = lambda x: find_microsynteny_chromosome_pair(*x, max_distance, distance_cutoff, nanosynteny_minsize, chrom_info_A, chrom_info_B)    

    with Pool(processes=n_proc) as pool:
        find_blocks = partial(find_microsynteny_chromosome_pair)
        blocks = pool.imap(find_blocks,condensed_dots_nanosynteny_all_chrom_pairs,chunksize=round(n_pairs/n_proc)+1)
        pool.close()
        pool.join()    

    return blocks
