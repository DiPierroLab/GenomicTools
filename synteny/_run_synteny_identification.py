
import numpy as np
from ._identify_blocks import *
from ._filter_blocks import *
from ._convolution_filters import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def run_synteny_identification(dot_plot, species_data_A, species_data_B, chrom_info_A, chrom_info_B, params):
    if 'tandem_windowsize' not in params.keys():
        params['tandem_windowsize'] = 1
    if 'block_overlap_threshold' not in params.keys():
        params['block_overlap_threshold'] = 1
    if 'dist_cutoff' not in params.keys():
        params['dist_cutoff'] = 10 * params['dot_maxdist']
    if 'return_absolute_blocks' not in params.keys():
        params['return_absolute_blocks'] = False 
    if 'verbose' not in params.keys():
        params['verbose'] = False
    chromsA = alphanum_sort(np.unique(dot_plot[:,0]))
    chromsB = alphanum_sort(np.unique(dot_plot[:,2]))
    maps_A = create_shift_map(species_data_A, params['tandem_windowsize'])
    maps_B = create_shift_map(species_data_B, params['tandem_windowsize'])
    shifted_dots = shift_dots(dot_plot, species_data_A, species_data_B, maps_A, maps_B)
    synteny_blocks = []
    I = 0
    for chromA in chromsA:
        for chromB in chromsB:
            if params['verbose'] == True:
                print(chromA,chromB,str(I)+"/%i"%(len(chromsA)*len(chromsB)),end='\r',flush=True)
            I += 1
            shifted_dots_AB = shifted_dots[(shifted_dots[:,0] == chromA)*(shifted_dots[:,2] == chromB)]
            H = np.zeros((chrom_info_A[chromA]['size'],chrom_info_B[chromB]['size']))
            H[shifted_dots_AB[:,1].astype(int)-1,shifted_dots_AB[:,3].astype(int)-1] = 1
            Mr = nanosynteny_convolve_dot_plot(H,params['block_minsize'])
            if Mr.sum() == 0:
                continue
            else:
                R = convolve_deconvolve_maxdist_dot_plot(Mr, H, params['dist_cutoff'])
                n_active_dots = R.sum()
                while True:
                    R = convolve_deconvolve_maxdist_dot_plot(R, H, params['dist_cutoff'])
                    if R.sum() == n_active_dots:
                        break
                    else:
                        n_active_dots = R.sum()
                xr, yr = np.where(R == 1)
                shifted_dots_AB = np.vstack([np.array(xr.shape[0]*[chromA]),xr+1,np.array(yr.shape[0]*[chromB]),yr+1]).T
                blocks_AB = identify_blocks_chrom_pair(shifted_dots_AB,params['block_minsize'],params['dot_maxdist'],maps_A,maps_B)
                if len(blocks_AB) > 0:
                    synteny_blocks += blocks_AB
    if len(synteny_blocks) == 0:
        return synteny_blocks
    unshifted_blocks = unshift_synteny_blocks(synteny_blocks, maps_A, maps_B, params['block_minsize'])
    unshifted_blocks_int = []
    for b in unshifted_blocks:
        block_int = []
        for i in b:
            block_int.append(np.array([chrom_info_A[i[0]]['number'],i[1],chrom_info_B[i[2]]['number'],i[3]]).astype(int))
        unshifted_blocks_int.append(np.vstack(block_int))
    fixed_blocks = fix_blocks(unshifted_blocks_int, params['block_minsize'],params['block_overlap_threshold'])
    
    if params['return_absolute_blocks'] == True:
        chrom_locs_A = np.cumsum([0] + [chrom_info_A[key]['size'] for key in chrom_info_A.keys()])
        chrom_locs_B = np.cumsum([0] + [chrom_info_B[key]['size'] for key in chrom_info_B.keys()])
        abs_A = {alphanum_sort(chrom_info_A.keys())[n]:s for n, s in enumerate(chrom_locs_A[:-1])}
        abs_B = {alphanum_sort(chrom_info_B.keys())[n]:s for n, s in enumerate(chrom_locs_B[:-1])}
        abs_int_A = {n+1:s for n, s in enumerate(chrom_locs_A[:-1])}
        abs_int_B = {n+1:s for n, s in enumerate(chrom_locs_B[:-1])}    
        absolute_blocks = []
        for block in fixed_blocks:
            absolute_block = []
            for b in block:
                absolute_block.append([abs_int_A[b[0]]+int(b[1]),abs_int_B[b[2]]+int(b[3])])
            absolute_blocks.append(np.vstack(absolute_block))
        return absolute_blocks
    else:
        return fixed_blocks
