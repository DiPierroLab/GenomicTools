
import numpy as np
from GenomicTools.tools import *

def shift_dots(dot_plot, species_data_A, species_data_B, maps_A, maps_B):
    cc_maps_A, inv_cc_maps_A, shift_maps_A, unshift_maps_A = maps_A
    cc_maps_B, inv_cc_maps_B, shift_maps_B, unshift_maps_B = maps_B
    chromsA = alphanum_sort(np.unique(dot_plot[:,0]))
    chromsB = alphanum_sort(np.unique(dot_plot[:,2]))    
    shifted_dots = [] 
    for chromA in chromsA:
        shiftA = shift_maps_A[chromA.rstrip('A')]
        ccA = cc_maps_A[chromA.rstrip('A')]
        for chromB in chromsB:
            shiftB = shift_maps_B[chromB.rstrip('B')]
            ccB = cc_maps_B[chromB.rstrip('B')]
            dotsAB = dot_plot[(dot_plot[:,0] == chromA)*(dot_plot[:,2] == chromB)][:,np.array([1,3])].astype(int)
            if dotsAB.shape[0] > 0:
                preshifted_dots_A = np.array([shiftA[ccA[key]] for key in dotsAB[:,0]])
                preshifted_dots_B = np.array([shiftB[ccB[key]] for key in dotsAB[:,1]])
                pair, indices = np.unique([str(a)+'-'+str(b) for a,b in zip(preshifted_dots_A,preshifted_dots_B)],return_index=True)
                shifted_dots_A = preshifted_dots_A[np.sort(indices)]
                shifted_dots_B = preshifted_dots_B[np.sort(indices)]
                N = shifted_dots_A.shape[0]
                shifted_dots.append(np.vstack([np.array(N*[chromA]),shifted_dots_A,np.array(N*[chromB]),shifted_dots_B]).T)
    shifted_dots = np.vstack(shifted_dots)
    return shifted_dots

def unshift_synteny_blocks(synteny_blocks, maps_A, maps_B):
    if type(synteny_blocks) != list:
        raise ValueError("The input synteny_blocks must be a list of synteny block arrays.")    
    cc_maps_A, inv_cc_maps_A, shift_maps_A, unshift_maps_A = maps_A
    cc_maps_B, inv_cc_maps_B, shift_maps_B, unshift_maps_B = maps_B
    unshifted_synteny_blocks = []
    for block in synteny_blocks:
        block = block[np.argsort(block[:,1]).astype(int)]
        slope = block_slope(block)
        chromA, chromB = block[0,np.array([0,2])]
        unshifted_dots_A = []
        unshifted_dots_B = []
        half_block_size = (block.shape[0]/2.)
        for nb, b in enumerate(block):
            cc_block_A = inv_cc_maps_A[chromA][unshift_maps_A[chromA][int(b[1])]]
            cc_block_B = inv_cc_maps_B[chromB][unshift_maps_B[chromB][int(b[3])]]
            dimA = cc_block_A.shape[0]
            dimB = cc_block_B.shape[0]
            if (dimA == 1) and (dimB == 1):
                unshifted_dots_A.append(cc_block_A[0])
                unshifted_dots_B.append(cc_block_B[0])
            else:
                dim = np.min([dimA, dimB])
                if (slope < 0) and (nb > half_block_size):
                    unshifted_dots_A += list(np.sort(cc_block_A)[:dim])
                    unshifted_dots_B += list(np.sort(cc_block_B)[:dim][::-1])
                elif (slope < 0) and (nb <= half_block_size):
                    unshifted_dots_A += list(np.sort(cc_block_A)[-dim:])
                    unshifted_dots_B += list(np.sort(cc_block_B)[-dim:][::-1])
                elif (slope > 0) and (nb > half_block_size):
                    unshifted_dots_A += list(np.sort(cc_block_A)[:dim])
                    unshifted_dots_B += list(np.sort(cc_block_B)[:dim])
                elif (slope > 0) and (nb <= half_block_size):
                    unshifted_dots_A += list(np.sort(cc_block_A)[-dim:])
                    unshifted_dots_B += list(np.sort(cc_block_B)[-dim:])
                else:
                    raise ValueError("Something is up with this synteny block.")
        unshifted_dots_A = np.array(unshifted_dots_A)
        unshifted_dots_B = np.array(unshifted_dots_B) 
        unshifted_length = unshifted_dots_A.shape[0]
        unshifted_block = np.vstack([np.array(unshifted_length*[chromA]),unshifted_dots_A,np.array(unshifted_length*[chromB]),unshifted_dots_B]).T
        unshifted_synteny_blocks.append(unshifted_block)
    return unshifted_synteny_blocks 
