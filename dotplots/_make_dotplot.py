def make_dotplot(dots, chrom_locs1, chrom_locs2, perm1, perm2, synteny_blocks = None, xlim = None, ylim = None, xfrac_lim = None, yfrac_lim = None, zoom_specific = None, highlight_zoom = False, label1 = None, label2 = None, block_color = 'r', min_block_size = 2, fontsize = 24, line_width = 2, dot_size = .5):
    size_x = 20
    
    if zoom_specific != None:
        chrA_zoom, indA_zoom, sizeA_zoom, chrB_zoom, indB_zoom, sizeB_zoom = zoom_specific
        aspect1 = sizeA_zoom
        aspect2 = sizeB_zoom
    elif zoom_specific == None:
        if (type(xlim) == int) and (xfrac_lim == None):
            aspect1 = chrom_locs1[xlim] - chrom_locs1[xlim-1]
        elif (type(xlim) == int) and (xfrac_lim != None):
            aspect1 = (chrom_locs1[xlim] - chrom_locs1[xlim-1]) * (np.max(xfrac_lim)-np.min(xfrac_lim))        
        else:
            aspect1 = chrom_locs1[-1]

        if (type(ylim) == int) and (yfrac_lim == None):
            aspect2 = chrom_locs2[ylim] - chrom_locs2[ylim-1]
        elif (type(ylim) == int) and (yfrac_lim != None):
            aspect2 = (chrom_locs2[ylim] - chrom_locs2[ylim-1]) * (np.max(yfrac_lim)-np.min(yfrac_lim))
        else:
            aspect2 = chrom_locs2[-1]  
    size_y = size_x * (aspect2 / aspect1)
    if size_y > 25:
        size_y = 25
        size_x = size_y * (aspect1 / aspect2)
    
    f, ax = plt.subplots(figsize=(size_x,size_y))
    ax.scatter(dots[:,0],dots[:,1],s=dot_size,zorder=0)

    nc1 = 0
    x_label_locs = []
    x_labels = []
    for chrom1 in chrom_locs1:
        plt.plot([chrom1,chrom1],[chrom_locs2[0],chrom_locs2[-1]],c='k',lw=.5,zorder=10)
        if nc1 > 0:
            x_label_locs.append((chrom_locs1[nc1] + chrom_locs1[nc1 - 1])/2)
            if perm1[nc1-1][1] < 0:
                x_labels.append("-"+perm1[nc1-1][0])
            else:
                x_labels.append(perm1[nc1-1][0])
        nc1 += 1

    nc2 = 0
    y_label_locs = []
    y_labels = []
    for chrom2 in chrom_locs2:
        plt.plot([chrom_locs1[0],chrom_locs1[-1]],[chrom2,chrom2],c='k',lw=.5,zorder=10)
        if nc2 > 0:
            y_label_locs.append((chrom_locs2[nc2] + chrom_locs2[nc2 - 1])/2)
            if perm2[nc2-1][1] < 0:
                y_labels.append("-"+perm2[nc2-1][0])
            else:
                y_labels.append(perm2[nc2-1][0])
        nc2 += 1  
    
    ax.set_xticks(x_label_locs)
    ax.set_xticklabels(x_labels,rotation=90, fontsize=round(.5*fontsize));
    ax.xaxis.set_ticks_position('none')
    ax.set_yticks(y_label_locs)
    ax.set_yticklabels(y_labels, fontsize=round(.5*fontsize));
    ax.yaxis.set_ticks_position('none')
    
    if label1 != None:
        ax.set_xlabel(r"%s"%label1, fontsize=fontsize, labelpad=10)
        ax.xaxis.set_label_position('top') 
    
    if label2 != None:    
        ax.set_ylabel(r"%s"%label2, fontsize=fontsize, rotation=270, labelpad=35)
        ax.yaxis.set_label_position('right') 
    
    if zoom_specific != None:
        chrA_zoom0 = chrom_locs1[chrA_zoom-1] + indA_zoom
        ax.set_xlim([chrA_zoom0-sizeA_zoom,chrA_zoom0+sizeA_zoom])
    elif xlim == None:
        ax.set_xlim([chrom_locs1[0],chrom_locs1[-1]])
    elif (type(xlim) == int) and (xfrac_lim == None):
        ax.set_xlim(chrom_locs1[xlim-1:xlim+1])
    elif (type(xlim) == int) and (xfrac_lim != None):
        xlims = chrom_locs1[xlim-1:xlim+1]
        xchrom_min = np.min(xlims)
        xchrom_max = np.max(xlims)
        xlen = xchrom_max - xchrom_min
        xmin = xchrom_min + xfrac_lim[0] * xlen
        xmax = xchrom_min + xfrac_lim[1] * xlen
        ax.set_xlim([xmin,xmax])
    else:
        ax.set_xlim(xlim)
    
    if zoom_specific != None:
        chrB_zoom0 = chrom_locs2[chrB_zoom-1] + indB_zoom
        ax.set_ylim([chrB_zoom0-sizeB_zoom,chrB_zoom0+sizeB_zoom])
        if highlight_zoom == True:
            ax.scatter([chrA_zoom0],[chrB_zoom0],marker="o",facecolors='none',edgecolors='blue',s=10*dot_size)
    elif ylim == None:
        ax.set_ylim([chrom_locs2[0],chrom_locs2[-1]])
    elif (type(ylim) == int) and (yfrac_lim == None):
        ax.set_ylim(chrom_locs2[ylim-1:ylim+1])
    elif (type(ylim) == int) and (yfrac_lim != None):
        ylims = chrom_locs2[ylim-1:ylim+1]
        ychrom_min = np.min(ylims)
        ychrom_max = np.max(ylims)
        ylen = ychrom_max - ychrom_min
        ymin = ychrom_min + yfrac_lim[0] * ylen
        ymax = ychrom_min + yfrac_lim[1] * ylen
        ax.set_ylim([ymin,ymax])
    else:
        ax.set_ylim(ylim)   
    
    if synteny_blocks != None:
        for synteny_block in synteny_blocks:
            if synteny_block.shape[0] >= min_block_size:
                if block_color == 'r':
                    ax.plot(synteny_block[:,0],synteny_block[:,1],lw=line_width,c='r',zorder=5)
                elif block_color == 'many':
                    ax.plot(synteny_block[:,0],synteny_block[:,1],lw=line_width,zorder=5)
                else:
                    raise ValueError("Specify a valid block color")
