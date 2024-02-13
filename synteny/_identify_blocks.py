import numpy as np
import networkx as ntx
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *
from ._filter_blocks import *

def identify_blocks_chrom_pair(dots, minsize = 2, filter_maxdist = False, maxdist = 1):
    chromsA = np.unique(dots[:,0])
    chromsB = np.unique(dots[:,2])
    if (chromsA.shape[0] != 1) or (chromsB.shape[0] != 1):
        raise ValueError("The input 'dots' should be all the dots between two chromosomes.")
    chromA = chromsA[0]
    chromB = chromsB[0]
    
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
    
    if filter_maxdist == True:
        Dx = np.abs((x.reshape([n,1]) * np.ones([n,n]) - x.reshape([1,n]) * np.ones([n,n])).astype(int))
        Dy = np.abs((y.reshape([n,1]) * np.ones([n,n]) - y.reshape([1,n]) * np.ones([n,n])).astype(int))
        D = ((Dx <= maxdist) * (Dy <= maxdist)).astype(int)
        Ap *= D
        Am *= D
         
    Gp = ntx.DiGraph(Ap)   
    Gm = ntx.DiGraph(Am)

    CCp = list(ntx.weakly_connected_components(Gp))
    CCm = list(ntx.weakly_connected_components(Gm))
  
    for ccp in CCp:
        G_block_p = ntx.subgraph(Gp,ccp)
        block_ind = ntx.dag_longest_path(G_block_p)
        block = dots[np.array(block_ind)]                
        if block.shape[0] >= minsize:
            col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
            col2 = block[:,1:2].astype(int)
            col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
            col4 = block[:,3:4].astype(int)
            blocks.append(np.hstack([col1,col2,col3,col4]))
                
    for ccm in CCm:
        G_block_m = ntx.subgraph(Gm,ccm)
        block_ind = ntx.dag_longest_path(G_block_m)
        block = dots[np.array(block_ind)]
        if block.shape[0] >= minsize:
            col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
            col2 = block[:,1:2].astype(int)
            col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
            col4 = block[:,3:4].astype(int)
            blocks.append(np.hstack([col1,col2,col3,col4]))

    return blocks
