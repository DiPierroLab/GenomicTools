import numpy as np
import networkx as ntx
from scipy.signal import fftconvolve
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *
from ._filter_blocks import *

def identify_blocks_chrom_pair(dots, minsize, maxdist, maps_A, maps_B):
    chromsA = np.unique(dots[:,0])
    chromsB = np.unique(dots[:,2])
    if (chromsA.shape[0] != 1) or (chromsB.shape[0] != 1):
        raise ValueError("The input 'dots' should be all the dots between two chromosomes.")
    chromA = chromsA[0]
    chromB = chromsB[0]
    cc_maps_A, inv_cc_maps_A, shift_maps_A, unshift_maps_A = maps_A
    cc_maps_B, inv_cc_maps_B, shift_maps_B, unshift_maps_B = maps_B
    unshift_map_A = unshift_maps_A[chromA]
    inv_cc_map_A = inv_cc_maps_A[chromA]
    unshift_map_B = unshift_maps_B[chromB]    
    inv_cc_map_B = inv_cc_maps_B[chromB]

    blocks = []
    x = dots[:,1].astype(int)
    y = dots[:,3].astype(int)
    n = x.shape[0]
    
    Ax = (x.reshape([n,1]) * np.ones([n,n]) < x.reshape([1,n]) * np.ones([n,n])).astype(int)
    Ay = (y.reshape([n,1]) * np.ones([n,n]) < y.reshape([1,n]) * np.ones([n,n])).astype(int)
    Ap = Ax * Ay
            
    Ax = (x.reshape([n,1]) * np.ones([n,n]) < x.reshape([1,n]) * np.ones([n,n])).astype(int)
    Ay = (y.reshape([n,1]) * np.ones([n,n]) > y.reshape([1,n]) * np.ones([n,n])).astype(int)
    Am = Ax * Ay
    
    Dx = np.abs((x.reshape([n,1]) * np.ones([n,n]) - x.reshape([1,n]) * np.ones([n,n])).astype(int))
    Dy = np.abs((y.reshape([n,1]) * np.ones([n,n]) - y.reshape([1,n]) * np.ones([n,n])).astype(int))
    D = ((Dx <= maxdist) * (Dy <= maxdist)).astype(int)

    cc_sizes_x = np.array([inv_cc_map_A[unshift_map_A[i]].shape[0] for i in x])
    cc_sizes_y = np.array([inv_cc_map_B[unshift_map_B[i]].shape[0] for i in y])
    cc_sizes = np.vstack([cc_sizes_x,cc_sizes_y]).min(0)
    S = np.ones([n,n]) * cc_sizes

    Gp = ntx.DiGraph(Ap * D * S)   
    Gm = ntx.DiGraph(Am * D * S)

    CCp = list(ntx.weakly_connected_components(Gp))
    CCm = list(ntx.weakly_connected_components(Gm))
  
    for ccp in CCp:
        G_block_p = ntx.subgraph(Gp,ccp)
        block_ind = ntx.dag_longest_path(G_block_p)
        block = dots[np.array(block_ind)]
        supported = supported_by_nanosynteny(block, minsize)               
        if supported: #block.shape[0] >= minsize:
            col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
            col2 = block[:,1:2].astype(int)
            col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
            col4 = block[:,3:4].astype(int)
            blocks.append(np.hstack([col1,col2,col3,col4]))
                
    for ccm in CCm:
        G_block_m = ntx.subgraph(Gm,ccm)
        block_ind = ntx.dag_longest_path(G_block_m)
        block = dots[np.array(block_ind)]
        supported = supported_by_nanosynteny(block, minsize) 
        if supported: #block.shape[0] >= minsize:        
            col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
            col2 = block[:,1:2].astype(int)
            col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
            col4 = block[:,3:4].astype(int)
            blocks.append(np.hstack([col1,col2,col3,col4]))

    return blocks
