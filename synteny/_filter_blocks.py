import numpy as np
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *
from ._identify_blocks import *

def split_plus_minus_chrom_vs_chrom(blocksAB, minsize=2):
    b_plus = []
    plus_order = []
    b_minus = []
    minus_order = []
    for b in blocksAB:
        b_sorted = b[np.argsort(b[:,1].astype(int))]
        s = block_slope(b_sorted)
        if (s > 0) and (b_sorted.shape[0] >= minsize):
            b_plus.append(b_sorted)
            plus_order.append(b_sorted[0,1])
        elif (s < 0) and (b_sorted.shape[0] >= minsize):
            b_minus.append(b_sorted)
            minus_order.append(b_sorted[0,1])
    b_plus_sorted = [b_plus[i] for i in np.argsort(plus_order)]
    b_minus_sorted = [b_minus[i] for i in np.argsort(minus_order)]

    return b_plus_sorted, b_minus_sorted

def define_forbidden_zones(b_plus_sorted, b_minus_sorted):
    plus_zones_x = []
    plus_zones_y = []
    for b in b_plus_sorted:
        bx = b[:,1].astype(int)
        by = b[:,3].astype(int)
        xmin = bx.min()
        xmax = bx.max()
        ymin = by.min()
        ymax = by.max()
        plus_zones_x.append([[xmin,xmax]])
        plus_zones_y.append([[ymin,ymax]])

    if len(plus_zones_x) > 0:
        plus_zones_x = np.vstack(plus_zones_x)
        plus_zones_y = np.vstack(plus_zones_y)
    else:
        plus_zones_x = np.array([])
        plus_zones_y = np.array([])

    minus_zones_x = []
    minus_zones_y = []
    for b in b_minus_sorted:
        bx = b[:,1].astype(int)
        by = b[:,3].astype(int)
        xmin = bx.min()
        xmax = bx.max()
        ymin = by.min()
        ymax = by.max()        
        minus_zones_x.append([[xmin,xmax]])
        minus_zones_y.append([[ymin,ymax]])

    if len(minus_zones_x) > 0:
        minus_zones_x = np.vstack(minus_zones_x)
        minus_zones_y = np.vstack(minus_zones_y)
    else:
        minus_zones_x = np.array([])
        minus_zones_y = np.array([])

    return plus_zones_x, plus_zones_y, minus_zones_x, minus_zones_y

def same_slope_conflict(b_sorted, zones, minsize, overlap_threshold):
    zones_x, zones_y = zones
    N = zones_x.shape[0]
    remove = []
    for ni in range(N):
        si = b_sorted[ni].shape[0]
        if si < minsize:
            remove.append(ni)
            continue
        for nj in range(ni+1,N):
            if (ni in remove) or (nj in remove):
                continue
            olx = interval_overlap(zones_x[ni],zones_x[nj])
            oly = interval_overlap(zones_y[ni],zones_y[nj])
            if (olx >= overlap_threshold) and (oly >= overlap_threshold):
                sj = b_sorted[nj].shape[0]
                if sj < minsize:
                    remove.append(nj)
                    continue
                if si > sj:
                    remove.append(nj)
                elif sj > si:
                    remove.append(ni)
                elif si == sj:
                    remove.append(ni)
                    remove.append(nj)    
    return remove

def different_slope_conflict(b_plus_sorted, b_minus_sorted, zones, pp_remove, mm_remove, minsize=3, overlap_threshold=1):
    plus_zones_x, plus_zones_y, minus_zones_x, minus_zones_y = zones
    nplus = plus_zones_x.shape[0]
    nminus = minus_zones_x.shape[0] 
    pp_splits = []
    mm_splits = []
    pp_split_remove = []
    mm_split_remove = []
    for ni in range(nplus):
        for nj in range(nminus):
            if (ni in pp_remove) or (nj in mm_remove):
                continue
            olx = interval_overlap(plus_zones_x[ni],minus_zones_x[nj])
            oly = interval_overlap(plus_zones_y[ni],minus_zones_y[nj])
            if (olx >= overlap_threshold) or (oly >= overlap_threshold):
                si = b_plus_sorted[ni].shape[0]
                sj = b_minus_sorted[nj].shape[0]

                xi1, xi2 = plus_zones_x[ni]
                yi1, yi2 = plus_zones_y[ni]
                xcondition = (b_minus_sorted[nj][:,1] >= xi1) * (b_minus_sorted[nj][:,1] <= xi2)
                ycondition = (b_minus_sorted[nj][:,3] >= yi1) * (b_minus_sorted[nj][:,3] <= yi2)
                ol_genesj = b_minus_sorted[nj][xcondition * ycondition]
                ol_genesjx = b_minus_sorted[nj][xcondition]
                ol_genesjy = b_minus_sorted[nj][ycondition]
                supportedj = supported_by_nanosynteny(ol_genesj,minsize)               

                xj1, xj2 = minus_zones_x[nj]
                yj1, yj2 = minus_zones_y[nj]
                xcondition = (b_plus_sorted[ni][:,1] >= xj1) * (b_plus_sorted[ni][:,1] <= xj2)
                ycondition = (b_plus_sorted[ni][:,3] >= yj1) * (b_plus_sorted[ni][:,3] <= yj2)
                ol_genesi = b_plus_sorted[ni][xcondition * ycondition]
                ol_genesix = b_plus_sorted[ni][xcondition]
                ol_genesiy = b_plus_sorted[ni][ycondition]
                supportedi = supported_by_nanosynteny(ol_genesi,minsize)
 
                if si < sj:
                    if (ol_genesjx.shape[0] < si) and (ol_genesjy.shape[0] < si) and (not supportedj):
                        mm_split_remove.append(nj)
                        if ol_genesjx.shape[0] == 0:
                            x_split1 = (xi1 + xi2) / 2.
                            x_split2 = (xi1 + xi2) / 2.
                        else:
                            x_split1 = ol_genesjx[:,1].min()
                            x_split2 = ol_genesjx[:,1].max()
                        mm_splits.append([nj,x_split1,x_split2])
                    else:
                        pp_remove.append(ni)
                elif si > sj:
                    if (ol_genesix.shape[0] < sj) and (ol_genesiy.shape[0] < sj) and (not supportedi):
                        pp_split_remove.append(ni)
                        if ol_genesix.shape[0] == 0:
                            x_split1 = (xj1 + xj2) / 2.
                            x_split2 = (xj1 + xj2) / 2.
                        else:
                            x_split1 = ol_genesix[:,1].min()
                            x_split2 = ol_genesix[:,1].max()
                        pp_splits.append([ni,x_split1,x_split2])
                    else:
                        mm_remove.append(nj)
                else:
                    pp_remove.append(ni)
                    mm_remove.append(nj)
    pp_remove = pp_remove + list(np.unique(pp_split_remove))
    mm_remove = mm_remove + list(np.unique(mm_split_remove))
    return pp_remove, mm_remove, pp_splits, mm_splits

def find_resolve_conflicts(blocksAB, minsize=3, overlap_threshold=1, large_block_dot_threshold=3):
    b_plus_sorted, b_minus_sorted = split_plus_minus_chrom_vs_chrom(blocksAB, minsize=minsize)
    plus_zones_x, plus_zones_y, minus_zones_x, minus_zones_y = define_forbidden_zones(b_plus_sorted, b_minus_sorted)

    pp_remove = same_slope_conflict(b_plus_sorted, [plus_zones_x, plus_zones_y], minsize, overlap_threshold)
    mm_remove = same_slope_conflict(b_minus_sorted, [minus_zones_x, minus_zones_y], minsize, overlap_threshold)                   
                        
    zones = [plus_zones_x, plus_zones_y, minus_zones_x, minus_zones_y]
    pp_remove, mm_remove, pp_splits, mm_splits = different_slope_conflict(b_plus_sorted, b_minus_sorted, zones, pp_remove, mm_remove, minsize=minsize, overlap_threshold=overlap_threshold)

    pp_split_blocks = []
    mm_split_blocks = []
    
    if len(pp_splits) > 0:
        pp_splits_stacked = np.vstack(pp_splits)
        pp_splits_n = np.unique(pp_splits_stacked[:,0])
        for pp_n in pp_splits_n.astype(int):
            pp_block = b_plus_sorted[pp_n]
            pp_splits_locs = np.sort(pp_splits_stacked[pp_splits_stacked[:,0] == pp_n][:,1:])
            for split in pp_splits_locs:
                pp_piece = pp_block[pp_block[:,1] < split[0]]
                pp_block = pp_block[pp_block[:,1] > split[1]]
                pp_split_blocks.append(pp_piece)
            pp_split_blocks.append(pp_block)

    if len(mm_splits) > 0:
        mm_splits_stacked = np.vstack(mm_splits)
        mm_splits_n = np.unique(mm_splits_stacked[:,0])
        for mm_n in mm_splits_n.astype(int):
            mm_block = b_minus_sorted[mm_n]
            mm_splits_locs = np.sort(mm_splits_stacked[mm_splits_stacked[:,0] == mm_n][:,1:])
            for split in mm_splits_locs:
                mm_piece = mm_block[mm_block[:,1] < split[0]]
                mm_block = mm_block[mm_block[:,1] > split[1]]
                mm_split_blocks.append(mm_piece)
            mm_split_blocks.append(mm_block)

    nplus = plus_zones_x.shape[0]
    nminus = minus_zones_x.shape[0]
    p_keep = set(np.arange(nplus)) - set(pp_remove)
    m_keep = set(np.arange(nminus)) - set(mm_remove)
    p_final = []
    m_final = []
    if len(p_keep) > 0:
        p_final = [b_plus_sorted[i] for i in p_keep]
        p_final += pp_split_blocks
    else:
        p_final += pp_split_blocks

    if len(m_keep) > 0:
        m_final = [b_minus_sorted[i] for i in m_keep]
        m_final += mm_split_blocks
    else:
        m_final += mm_split_blocks
    both_final_init = p_final + m_final              
    both_final = [b[np.argsort(b[:,1])] for b in both_final_init if b.shape[0] >= minsize]
    return both_final

def fix_blocks(blocks, minsize, overlap_threshold, large_block_dot_threshold):
    stacked_blocks = np.vstack(blocks)
    #chromsA = alphanum_sort(np.unique(stacked_blocks[:,0]))
    #chromsB = alphanum_sort(np.unique(stacked_blocks[:,2]))
    chromsA = np.sort(np.unique(stacked_blocks[:,0]))
    chromsB = np.sort(np.unique(stacked_blocks[:,2]))

    fixed_blocks = []
    for chromA in chromsA:
        for chromB in chromsB:
            blocksAB = [b for b in blocks if (b[0,0] == chromA) * (b[0,2] == chromB)]
            fixed_blocksAB = find_resolve_conflicts(blocksAB, minsize, overlap_threshold, large_block_dot_threshold)
            fixed_blocks += fixed_blocksAB
    return fixed_blocks    
