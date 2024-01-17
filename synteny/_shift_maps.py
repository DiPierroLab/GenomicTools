
def create_shift_map(data, windowsize):    
    """
    - data: OG data
    - windowsize: window for joining genes (e.g. 1 means two genes in the same OG will be joined 
                  only if they have indices i and i+1).
    """
    
    cc_maps = {}
    inv_cc_maps = {}
    shift_maps = {}
    unshift_maps = {}
    
    for chrom in data.keys():
        ogs = data[chrom]['ogs']
        locs = np.arange(1,ogs.shape[1]+1).reshape(ogs.shape)
        
        M_ogs = ogs * np.ones(ogs.shape[1]) - ogs.T * np.ones(ogs.shape[1])
        M_loc = locs * np.ones(ogs.shape[1]) - locs.T * np.ones(ogs.shape[1])

        A = ((M_ogs == 0) * (np.abs(M_loc) <= windowsize)).astype(int)
        G = ntx.Graph(A)
        CC = ntx.connected_components(G)

        cc_map = {}
        shift_map = {}
        inv_cc_map = {}
        unshift_map = {}
        multiplicity_map = {}
        for n, cc in enumerate(CC):
            inv_cc_map[n] = np.array(list(cc))
            for rep in cc:
                cc_map[rep] = n
                shift_map[n] = locs.flatten()[n]
                unshift_map[locs.flatten()[n]] = n

        cc_maps[chrom] = cc_map
        inv_cc_maps[chrom] = inv_cc_map
        shift_maps[chrom] = shift_map
        unshift_maps[chrom] = unshift_map
    
    # cc_map: gene index -> cc index
    # inv_cc_map: cc index -> set of gene indices
    # shift_map: cc index -> compressed gene index
    # unshift_map: compressed gene index -> cc index
    
    return cc_maps, inv_cc_maps, shift_maps, unshift_maps

def shift_input_data(dot_plot_result, og_data_A, og_data_B, cc_map_A, cc_map_B, shift_map_A, shift_map_B):
    data = dot_plot_result['data']
    chromsAB_all = np.vstack(list(data.keys()))
    chromsA = alphanum_sort(np.unique(chromsAB_all[:,0]))
    chromsB = alphanum_sort(np.unique(chromsAB_all[:,1]))
    
    shifted_dots = [] 
    for chromA in chromsA:
        shiftA = shift_map_A[chromA.rstrip('A')]
        ccA = cc_map_A[chromA.rstrip('A')]
        for chromB in chromsB:
            shiftB = shift_map_B[chromB.rstrip('B')]
            ccB = cc_map_B[chromB.rstrip('B')]
            dataAB = data[(chromA,chromB)]['homology_matrix'].toarray()
            dotsAB = np.vstack(np.where(dataAB == 1)).T
            if dotsAB.shape[0] > 0:
                preshifted_dots_A = np.array([shiftA[ccA[key]] for key in dotsAB[:,1].astype(int)])
                preshifted_dots_B = np.array([shiftB[ccB[key]] for key in dotsAB[:,0].astype(int)])
                pair, indices = np.unique([str(a)+'-'+str(b) for a,b in zip(preshifted_dots_A,preshifted_dots_B)],return_index=True)
                shifted_dots_A = preshifted_dots_A[np.sort(indices)]
                shifted_dots_B = preshifted_dots_B[np.sort(indices)]
                
                chrom_numA = np.array(shifted_dots_A.shape[0] * [chrom_to_num(chromA.rstrip('A'))])
                chrom_numB = np.array(shifted_dots_B.shape[0] * [chrom_to_num(chromB.rstrip('B'))])
                shifted_dots.append(np.vstack([chrom_numA,shifted_dots_A,chrom_numB,shifted_dots_B]).T)
                
    return np.vstack(shifted_dots)

def block_slope(block):
    lr = stats.linregress(block[:,1],block[:,3])
    return lr.slope
    
def unshift_blocks(blocks, dot_plot_result, inv_cc_map_A, inv_cc_map_B, unshift_map_A, unshift_map_B):
    spA = dot_plot_result['species1']
    spB = dot_plot_result['species2']
    unshifted_blocks = []
    for block in blocks:
        block = block[np.argsort(block[:,1])]
        slope = block_slope(block)
        chromA = num_to_chrom(block[0,0],spA)
        chromB = num_to_chrom(block[0,2],spB)
        unshifted_dots_A = []
        unshifted_dots_B = []
        half_block_size = (block.shape[0]/2.)
        for nb, b in enumerate(block):
            cc_block_A = inv_cc_map_A[chromA][unshift_map_A[chromA][b[1]]]
            cc_block_B = inv_cc_map_B[chromB][unshift_map_B[chromB][b[3]]]
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
                    print(block)
        unshifted_dots_A = np.array(unshifted_dots_A)
        unshifted_dots_B = np.array(unshifted_dots_B)
        
        unshifted_length = unshifted_dots_A.shape[0]
        unshifted_block = np.vstack([np.array(unshifted_length*[block[0,0]]),unshifted_dots_A,np.array(unshifted_length*[block[0,2]]),unshifted_dots_B]).T
        unshifted_blocks.append(unshifted_block)
    
    return unshifted_blocks


