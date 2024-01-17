def assemble_dotplot_chrom_permutations(block_matrix, perm1 = None, perm2 = None):
    """
    Permutations should be formatted as a dictionary, with the names of chromosomes as keys and a list
    with the index of the chromosome first and the sign second as the value.
    """
    chrom_pairs = np.vstack(list(block_matrix.keys()))
    
    chroms1 = alphanum_sort(np.unique(chrom_pairs[:,0]))
    nchroms1 = len(chroms1)
    base_perm1 = {i:[chroms1[i],1] for i in range(nchroms1)}
    
    chroms2 = alphanum_sort(np.unique(chrom_pairs[:,1]))
    nchroms2 = len(chroms2)
    base_perm2 = {i:[chroms2[i],1] for i in range(nchroms2)}
    
    if perm1 == None:
        perm1 = {}
        for key in base_perm1.keys():
            perm1[key] = base_perm1[key]
    elif type(perm1) == list:
        perm1_list = list(np.copy(perm1))
        perm1 = {np.abs(perm1_list[i]):[chroms1[i],np.sign(perm1_list[i])] for i in range(nchroms1)}    

    if perm2 == None:
        perm2 = {}
        for key in base_perm2.keys():
            perm2[key] = base_perm2[key]
    elif type(perm2) == list:
        perm2_list = list(np.copy(perm2))
        perm2 = {np.abs(perm2_list[i]):[chroms2[i],np.sign(perm2_list[i])] for i in range(nchroms2)}
        
    blocks = []
    for i in range(nchroms1):
        block_i = []
        c1 = perm1[i][0]
        sign1 = perm1[i][1]
        for j in range(nchroms2):
            c2 = perm2[j][0]
            sign2 = perm2[j][1]
            M = block_matrix[(c1,c2)]['homology_matrix'][::sign1,::sign2].T
            block_i.append(M.toarray())
        blocks.append(block_i)
    
    chrom_locs1 = []
    for block1 in blocks:
        chrom_locs1.append(block1[0].shape[0])
    chrom_locs1 = np.cumsum([0] + chrom_locs1)
    
    chrom_locs2 = []
    for block2 in blocks[0]:
        chrom_locs2.append(block2.shape[1])
    chrom_locs2 = np.cumsum([0] + chrom_locs2)
    
    B = np.block(blocks)
    dots = np.vstack(np.where(B == 1)).T
    del B
    return dots, chrom_locs1, chrom_locs2, perm1, perm2
