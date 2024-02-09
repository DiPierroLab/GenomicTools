import numpy as np
import networkx as ntx
from GenomicTools.tools import *
from ._shift_maps import *
from ._filter_blocks import *

def damnchainer(dots, max_dist, minsize = 2):

    chromsA = np.unique(dots[:,0])
    chromsB = np.unique(dots[:,2])
    
    DAGs = []
    
    damn_blocks = []
    for chromA in chromsA:
        for chromB in chromsB:
            dotsAB = dots[(dots[:,0] == chromA) * (dots[:,2] == chromB)]
            x = dotsAB[:,1].astype(int)
            y = dotsAB[:,3].astype(int)
            nAB = x.shape[0]
    
            Dx = np.abs((x.reshape([nAB,1]) * np.ones([nAB,nAB]) - x.reshape([1,nAB]) * np.ones([nAB,nAB])).astype(int))
            Dy = np.abs((y.reshape([nAB,1]) * np.ones([nAB,nAB]) - y.reshape([1,nAB]) * np.ones([nAB,nAB])).astype(int))
            D = ((Dx <= max_dist) * (Dy <= max_dist)).astype(int)

            Ax = (x.reshape([nAB,1]) * np.ones([nAB,nAB]) < x.reshape([1,nAB]) * np.ones([nAB,nAB])).astype(int)
            Ay = (y.reshape([nAB,1]) * np.ones([nAB,nAB]) < y.reshape([1,nAB]) * np.ones([nAB,nAB])).astype(int)
            Ap = Ax * Ay
            Gp = ntx.DiGraph(D * Ap)
            
            Ax = (x.reshape([nAB,1]) * np.ones([nAB,nAB]) < x.reshape([1,nAB]) * np.ones([nAB,nAB])).astype(int)
            Ay = (y.reshape([nAB,1]) * np.ones([nAB,nAB]) > y.reshape([1,nAB]) * np.ones([nAB,nAB])).astype(int)
            Am = Ax * Ay
            Gm = ntx.DiGraph(D * Am)
            
            CCp = list(ntx.weakly_connected_components(Gp))
            CCm = list(ntx.weakly_connected_components(Gm))
  
            for ccp in CCp:
                G_block_p = ntx.subgraph(Gp,ccp)
                block_ind = ntx.dag_longest_path(G_block_p)
                block = dotsAB[np.array(block_ind)]                
                if block.shape[0] >= minsize:
                    col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
                    col2 = block[:,1:2].astype(int)
                    col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
                    col4 = block[:,3:4].astype(int)
                    damn_blocks.append(np.hstack([col1,col2,col3,col4]))
                
                    DAGs.append(G_block_p)
                
            for ccm in CCm:
                G_block_m = ntx.subgraph(Gm,ccm)
                block_ind = ntx.dag_longest_path(G_block_m)
                block = dotsAB[np.array(block_ind)]
                if block.shape[0] >= minsize:
                    col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
                    col2 = block[:,1:2].astype(int)
                    col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
                    col4 = block[:,3:4].astype(int)
                    damn_blocks.append(np.hstack([col1,col2,col3,col4]))

                    DAGs.append(G_block_m)
    
    return damn_blocks#, DAGs

def get_the_damn_synteny_blocks(dot_plot_result, og_data_A, og_data_B, params, return_absolute = True):
    spA = dot_plot_result['species1']
    spB = dot_plot_result['species2']
    
    cc_map_A, inv_cc_map_A, shift_map_A, unshift_map_A = create_shift_map(og_data_A,params['shift_windowsize'])
    cc_map_B, inv_cc_map_B, shift_map_B, unshift_map_B = create_shift_map(og_data_B,params['shift_windowsize'])
    
    shifted_dots = shift_input_data(dot_plot_result, og_data_A, og_data_B, cc_map_A, cc_map_B, shift_map_A, shift_map_B)

    convolved_dots = nanosynteny_convolution_step(shifted_dots,params['convolution_kernel_size'],params['damn_maxdist'])

    extended_dots = extension_colvolution_step(convolved_dots, shifted_dots, maxdist)
    blocks = damnchainer(shifted_dots, params['damn_maxdist'], minsize = params['damn_minsize']) 
    shifted_dots = np.vstack(blocks)
    

    unshifted_blocks = unshift_blocks(blocks, dot_plot_result, inv_cc_map_A, inv_cc_map_B, unshift_map_A, unshift_map_B)    

    cc_map_A, inv_cc_map_A, shift_map_A, unshift_map_A = create_shift_map(og_data_A,params['shift_windowsize_post'])
    cc_map_B, inv_cc_map_B, shift_map_B, unshift_map_B = create_shift_map(og_data_B,params['shift_windowsize_post'])
    unshifted_blocks_filtered_for_orphan_dup_blocks = []
    for b in unshifted_blocks:
        b_condensed = [str(cc_map_A[num_to_chrom(i[0],spA)][i[1]])+'-'+str(cc_map_B[num_to_chrom(i[2],spB)][i[3]]) for i in b]
        if np.unique(b_condensed).shape[0] >= params['min_allowed_condensed_dots']:
            unshifted_blocks_filtered_for_orphan_dup_blocks.append(b)
    unshifted_blocks = unshifted_blocks_filtered_for_orphan_dup_blocks    
    
    fixed_blocks = fix_blocks(unshifted_blocks, params['damn_minsize'], params['overlap_threshold'], params['large_block_dot_threshold'])
    
    if return_absolute == False:
        return fixed_blocks
    elif return_absolute == True:
        return convert_to_absolute_indices(fixed_blocks, dot_plot_result)
    else:
        raise ValueError("Not a valid option for 'return_absolute'...")


# x=params['convolution_kernel_size'], maxdist=params['damn_maxdist']
# shifted_dots.append(np.vstack([chrom_numA,shifted_dots_A,chrom_numB,shifted_dots_B]).T.astype(int))

def nanosynteny_convolution_step(shifted_dots, x): 
    chromsA = alphanum_sort(np.unique(shifted_dots[:,0]))
    chromsB = alphanum_sort(np.unique(shifted_dots[:,2]))
    convolved_dots = []
    for chromA in chromsA:
        for chromB in chromsB:
             dataAB = shifted_dots[(shifted_dots[:,0] == chromA) * (shifted_dots[:,2] == chromB)] 
             shifted_dots_A = dataAB[:,1]
             shifted_dots_B = dataAB[:,3]
             M = np.zeros([shifted_dots_A.max()+1,shifted_dots_B.max()+1])
             M[shifted_dots_A, shifted_dots_B] = 1
             Cp, Cm = nanosynteny_convolve_dot_plot(M, x)
             Mr = nanosynteny_deconvolve_dotplot(M, Cp, Cm, x, maxdist)
             shifted_dots_A, shifted_dots_B = np.where(Mr == 1)
            
             convolved_dots.append(np.vstack([dataAB[:,0],shifted_dots_A,dataAB[:,2],shifted_dots_B]).T.astype(int))

    return np.vstack(convolved_dots)

def extension_colvolution_step(convolved_dots, shifted_dots, maxdist):
    chromsA = alphanum_sort(np.unique(shifted_dots[:,0]))
    chromsB = alphanum_sort(np.unique(shifted_dots[:,2]))
    convolved_dots = []
    for chromA in chromsA:
        for chromB in chromsB:
             dataAB = shifted_dots[(shifted_dots[:,0] == chromA) * (shifted_dots[:,2] == chromB)]
             shifted_dots_A = dataAB[:,1]
             shifted_dots_B = dataAB[:,3]
             M = np.zeros([shifted_dots_A.max()+1,shifted_dots_B.max()+1])
             M[shifted_dots_A, shifted_dots_B] = 1
             Cp, Cm = nanosynteny_convolve_dot_plot(M, x)
             Mr = nanosynteny_deconvolve_dotplot(M, Cp, Cm, x, maxdist)
             shifted_dots_A, shifted_dots_B = np.where(Mr == 1)
            
             convolved_dots.append(np.vstack([dataAB[:,0],shifted_dots_A,dataAB[:,2],shifted_dots_B]).T.astype(int))

    return np.vstack(convolved_dots)
