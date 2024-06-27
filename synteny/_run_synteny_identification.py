
import numpy as np
from ._identify_blocks import *
from ._filter_blocks import *
from ._convolution_filters import *
from ._block_tools import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *
from GenomicTools.permutations import *

def find_nanosynteny(condensed_dots, species_data_A, species_data_B, chrom_info_A, chrom_info_B, maps_A, maps_B, nanosynteny_minsize, verbose = False):
    chromsA = np.unique(condensed_dots[:,0])
    chromsB = np.unique(condensed_dots[:,2])
    blocks = []
    N = chromsA.shape[0] * chromsB.shape[0]
    n = 0
    for chromA in chromsA:
        for chromB in chromsB:
            condensed_dots_AB = condensed_dots[(condensed_dots[:,0] == chromA) * (condensed_dots[:,2] == chromB)]            
            if condensed_dots_AB.shape[0] > 0:
                blocks += find_nanosynteny_chromosome_pair(condensed_dots_AB, species_data_A, species_data_B, chrom_info_A, chrom_info_B, maps_A, maps_B, nanosynteny_minsize)
            n += 1
            if verbose:
                print('%i / %i'%(n,N),end='\r',flush=True)
    return blocks

def find_microsynteny(condensed_dots, nanosynteny_blocks, chrom_info_A, chrom_info_B, max_distance, distance_cutoff, nanosynteny_minsize, verbose = False):
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
