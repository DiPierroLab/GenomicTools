
def split_plus_minus_chrom_vs_chrom(blocksAB, minsize=2):
    b_plus = []
    plus_order = []
    b_minus = []
    minus_order = []
    for b in blocksAB:
        s = (b[-1,3] - b[0,3]) / (b[-1,1] - b[0,1])
        if (s > 0) and (b.shape[0] >= minsize):
            b_plus.append(b)
            plus_order.append(b[0,1])
        elif (s < 0) and (b.shape[0] >= minsize):
            b_minus.append(b)
            minus_order.append(b[0,1])
    b_plus_sorted = [b_plus[i] for i in np.argsort(plus_order)]
    b_minus_sorted = [b_minus[i][::-1] for i in np.argsort(minus_order)]

    return b_plus_sorted, b_minus_sorted

def define_forbidden_zones(b_plus_sorted, b_minus_sorted):
    plus_zones_x = []
    plus_zones_y = []
    for b in b_plus_sorted:
        xmin = b[:,1].min()
        xmax = b[:,1].max()
        ymin = b[:,3].min()
        ymax = b[:,3].max()
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
        xmin = b[:,1].min()
        xmax = b[:,1].max()
        ymin = b[:,3].min()
        ymax = b[:,3].max()
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
        ri = np.max([np.max(np.diff(b_sorted[ni][:,1])), np.max(np.diff(b_sorted[ni][:,3]))])
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
                    rj = np.max([np.max(np.diff(b_sorted[nj][:,1])), np.max(np.diff(b_sorted[nj][:,3]))])
                    if ri > rj:
                        remove.append(ni)
                    elif rj > ri:
                        remove.append(nj)
                    else:
                        remove.append(ni)
                        remove.append(nj)
    
    return remove

def find_resolve_conflicts(blocksAB, minsize=2, overlap_threshold=1, large_block_dot_threshold=3):
    b_plus_sorted, b_minus_sorted = split_plus_minus_chrom_vs_chrom(blocksAB, minsize=minsize)
    plus_zones_x, plus_zones_y, minus_zones_x, minus_zones_y = define_forbidden_zones(b_plus_sorted, b_minus_sorted)

    nplus = plus_zones_x.shape[0]
    nminus = minus_zones_x.shape[0]

    pp_remove = same_slope_conflict(b_plus_sorted, [plus_zones_x, plus_zones_y], minsize, overlap_threshold)
    mm_remove = same_slope_conflict(b_minus_sorted, [minus_zones_x, minus_zones_y], minsize, overlap_threshold)                   
                        
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
            if (olx >= overlap_threshold) and (oly >= overlap_threshold):
                si = b_plus_sorted[ni].shape[0]
                ri = np.max([np.max(np.diff(b_plus_sorted[ni][:,1])), np.max(np.diff(b_plus_sorted[ni][:,3]))])
                sj = b_minus_sorted[nj].shape[0]
                rj = np.max([np.max(np.diff(b_minus_sorted[nj][:,1])), np.max(np.diff(b_minus_sorted[nj][:,3]))])
                if si < sj:
                    xi1, xi2 = plus_zones_x[ni]
                    yi1, yi2 = plus_zones_y[ni]
                    xcondition = (b_minus_sorted[nj][:,1] >= xi1) * (b_minus_sorted[nj][:,1] <= xi2)
                    ycondition = (b_minus_sorted[nj][:,3] >= yi1) * (b_minus_sorted[nj][:,3] <= yi2)
                    ol_genesj = b_minus_sorted[nj][xcondition*ycondition]
                    ol_genesjx = b_minus_sorted[nj][xcondition*(ycondition == False)]
                    ol_genesjy = b_minus_sorted[nj][(xcondition == False)*ycondition]
                    if (ol_genesj.shape[0] < si) and (ol_genesjx.shape[0] < large_block_dot_threshold) and (ol_genesjy.shape[0] < large_block_dot_threshold):
                        mm_split_remove.append(nj)
                        if ol_genesj.shape[0] == 0:
                            x_split1 = (xi1 + xi2) / 2.
                            x_split2 = (xi1 + xi2) / 2.
                        else:
                            x_split1 = ol_genesj[:,1].min()
                            x_split2 = ol_genesj[:,1].max()
                        mm_splits.append([nj,x_split1,x_split2])
                    else:
                        pp_remove.append(ni)
                elif si > sj:
                    xj1, xj2 = minus_zones_x[nj]
                    yj1, yj2 = minus_zones_y[nj]
                    xcondition = (b_plus_sorted[ni][:,1] >= xj1) * (b_plus_sorted[ni][:,1] <= xj2)
                    ycondition = (b_plus_sorted[ni][:,3] >= yj1) * (b_plus_sorted[ni][:,3] <= yj2)
                    ol_genesi = b_plus_sorted[ni][xcondition * ycondition]
                    ol_genesix = b_plus_sorted[ni][xcondition*(ycondition == False)]
                    ol_genesiy = b_plus_sorted[ni][(xcondition == False)*ycondition]
                    if (ol_genesi.shape[0] < sj) and (ol_genesix.shape[0] < large_block_dot_threshold) and (ol_genesiy.shape[0] < large_block_dot_threshold):
                        pp_split_remove.append(ni)
                        if ol_genesi.shape[0] == 0:
                            x_split1 = (xj1 + xj2) / 2.
                            x_split2 = (xj1 + xj2) / 2.                        
                        else:
                            x_split1 = ol_genesi[:,1].min()
                            x_split2 = ol_genesi[:,1].max()
                        pp_splits.append([ni,x_split1,x_split2])
                    else:
                        mm_remove.append(nj)
                else:
                    if ri > rj:
                        pp_remove.append(ni)
                    elif ri < rj:
                        mm_remove.append(nj)
                    else:
                        pp_remove.append(ni)
                        mm_remove.append(nj)
                        
    pp_remove = pp_remove + list(np.unique(pp_split_remove))
    mm_remove = mm_remove + list(np.unique(mm_split_remove))

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

    ###
    both_final_init = p_final + m_final              
    both_final = [b[np.argsort(b[:,1])] for b in both_final_init if b.shape[0] >= minsize]
    ###
    
    return both_final

def fix_blocks(blocks, minsize, overlap_threshold, large_block_dot_threshold, filterblocks=False):
    stacked_blocks = np.vstack(blocks)
    chromsA = np.unique(stacked_blocks[:,0])
    chromsB = np.unique(stacked_blocks[:,2])

    fixed_blocks = []
    for chromA in chromsA:
        for chromB in chromsB:
            blocksAB = [b for b in blocks if (b[0,0] == chromA) * (b[0,2] == chromB)]
            fixed_blocksAB = find_resolve_conflicts(blocksAB, minsize, overlap_threshold, large_block_dot_threshold)
            fixed_blocks += fixed_blocksAB
    
    if filterblocks == True:
        final_blocks = []
        for b in fixed_blocks:
            if b.shape[0] >= minsize:
                hx = np.unique(b[:,1],return_counts=True)
                hy = np.unique(b[:,3],return_counts=True)
                repx = hx[0][np.where(hx[1] > 1)[0]]
                repy = hy[0][np.where(hy[1] > 1)[0]]
                indices = list(range(b.shape[0]))
                remove = []
                if (repx.shape[0] > 0):
                    for r in repx:
                        reps = np.where(b[:,1] == r)[0]
                        Rsq = []
                        for s in reps:
                            x0, y0 = b[s-1][np.array([1,3])]
                            x1, y1 = b[s][np.array([1,3])]
                            x2, y2 = b[s+1][np.array([1,3])]
                            x_exp = np.mean([x0,x2])
                            y_exp = np.mean([y0,y2])
                            Rsq.append((x1-x_exp)**2+(y1-y_exp)**2)
                        remove += list(set(indices) - set(list(np.argmin(Rsq))))
                if (repy.shape[0] > 0):
                    for r in repy:
                        reps = np.where(b[:,3] == r)[0]
                        Rsq = []
                        for s in reps:
                            x0, y0 = b[s-1][np.array([1,3])]
                            x1, y1 = b[s][np.array([1,3])]
                            x2, y2 = b[s+1][np.array([1,3])]
                            x_exp = np.mean([x0,x2])
                            y_exp = np.mean([y0,y2])
                            Rsq.append((x1-x_exp)**2+(y1-y_exp)**2)
                        remove += list(set(indices) - set(list(np.argmin(Rsq))))
                final_blocks.append(b[np.array(list(set(indices) - set(remove)))])
        return final_blocks
    else:
        return fixed_blocks
