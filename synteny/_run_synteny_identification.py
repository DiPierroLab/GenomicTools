
import numpy as np
from ._identify_blocks import *
from ._filter_blocks import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def extend_dots(dots, synteny_blocks, maxdist):
    dots_xy = dots[:,np.array([1,3])].astype(int)
    xd = dots_xy[:,0:1]
    yd = dots_xy[:,1:2]
    block_dots_xy = np.vstack(synteny_blocks)[:,np.array([1,3])].astype(int)
    xb = block_dots_xy[:,0:1]
    yb = block_dots_xy[:,1:2]
    xdists = xb * np.ones((xb.shape[0],xd.shape[0])) - xd.T * np.ones((xb.shape[0],xd.shape[0]))
    xclose = ((np.abs(xdists) <= maxdist).sum(0) > 0)
    ydists = yb * np.ones((yb.shape[0],yd.shape[0])) - yd.T * np.ones((yb.shape[0],yd.shape[0]))
    yclose = ((np.abs(ydists) <= maxdist).sum(0) > 0)
    extended_dots = dots[xclose * yclose]
    return extended_dots

def run_synteny_identification(dot_plot, species_data_A, species_data_B, chrom_info_A, chrom_info_B, params):
    chromsA = alphanum_sort(np.unique(dot_plot[:,0]))
    chromsB = alphanum_sort(np.unique(dot_plot[:,2]))
    maps_A = create_shift_map(species_data_A, params['tandem_windowsize'])
    maps_B = create_shift_map(species_data_B, params['tandem_windowsize'])
    shifted_dots = shift_dots(dot_plot, species_data_A, species_data_B, maps_A, maps_B)
    synteny_blocks = []
    for chromA in chromsA:
        for chromB in chromsB:
            shifted_dots_AB = shifted_dots[(shifted_dots[:,0] == chromA)*(shifted_dots[:,2] == chromB)]
            convolved_dots_AB = identify_blocks_chrom_pair(shifted_dots_AB,params['condensed_block_minsize'],1,maps_A,maps_B)
            if len(convolved_dots_AB) == 0:
                continue
            else:
                extended_dots_AB = extend_dots(shifted_dots_AB, convolved_dots_AB, params['dot_maxdist'])
                blocks_AB = identify_blocks_chrom_pair(extended_dots_AB,params['condensed_block_minsize'],params['dot_maxdist'],maps_A,maps_B)
                # print(len(convolved_dots_AB),extended_dots_AB.shape[0],len(blocks_AB))
                i = 0
                while True:
                    try:
                        dots_blocks_AB = np.vstack(blocks_AB)
                    except ValueError:
                        break
                    n_dots_in_blocks = dots_blocks_AB.shape[0]
                    extended_dots_AB = extend_dots(shifted_dots_AB, blocks_AB, params['dot_maxdist'])
                    blocks_AB = identify_blocks_chrom_pair(extended_dots_AB,params['condensed_block_minsize'],params['dot_maxdist'],maps_A,maps_B)
                    dots_blocks_AB = np.vstack(blocks_AB)
                    i += 1
                    if dots_blocks_AB.shape[0] == n_dots_in_blocks:
                        break
                    if i > params['max_iterations']:
                        raise ValueError('Surpassed maximum allowed iterations.')
                if len(blocks_AB) > 0:
                    synteny_blocks += blocks_AB
    unshifted_blocks = unshift_synteny_blocks(synteny_blocks, maps_A, maps_B, params['nanosynteny_minsize'])
    unshifted_blocks_int = []
    for b in unshifted_blocks:
        block_int = []
        for i in b:
            block_int.append(np.array([chrom_info_A[i[0]]['number'],i[1],chrom_info_B[i[2]]['number'],i[3]]).astype(int))
        unshifted_blocks_int.append(np.vstack(block_int))
    return unshifted_blocks
    #fixed_blocks = fix_blocks(unshifted_blocks_int, params['block_minsize'], params['block_overlap_threshold'], params['large_block_dot_threshold'])
    #return fixed_blocks
