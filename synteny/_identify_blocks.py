
import numpy as np
import networkx as ntx
from ._convolution_filters import *
from ._filter_blocks import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def find_nanosynteny_chromosome_pair(condensed_dots, species_data_A, species_data_B, chrom_info_A, chrom_info_B, maps_A, maps_B, nanosynteny_minsize, check_for_nanosynteny_support = True):
    nano_dots = get_nano_dots(condensed_dots, nanosynteny_minsize, chrom_info_A, chrom_info_B)
    if nano_dots.shape[0] > 0:
        Gm, Gp = create_chromosome_pair_dot_dags(nano_dots, 1)
        blocks = identify_blocks_chrom_pair(nano_dots, Gm, Gp, species_data_A, species_data_B, maps_A, maps_B, nanosynteny_minsize, check_for_nanosynteny_support)
    else:
        blocks = []
    return blocks

def find_microsynteny_chromosome_pair(condensed_dots, nanosynteny_blocks, max_distance, distance_cutoff, nanosynteny_minsize, chrom_info_A, chrom_info_B):
    nano_dots, nano_neighbor_dots = get_nano_neighborhood_dots(condensed_dots, nanosynteny_minsize, distance_cutoff, chrom_info_A, chrom_info_B)
    if (len(nanosynteny_blocks) > 0) and (nano_neighbor_dots.shape[0] > 0):
        blocks = connect_blocks(nanosynteny_blocks, nano_neighbor_dots, max_distance, distance_cutoff)
    else:
        blocks = []
    return blocks
    
def identify_blocks_chrom_pair(condensed_dots, Gm, Gp, species_data_A, species_data_B, maps_A, maps_B, nanosynteny_minsize, check_for_nanosynteny_support = True):
    chromsA = np.unique(condensed_dots[:,0])
    chromsB = np.unique(condensed_dots[:,2])
    if (chromsA.shape[0] != 1) or (chromsB.shape[0] != 1):
        print(chromsA.shape[0], chromsB.shape[0])
        raise ValueError("The input 'condensed_dots' should only be the dots between two chromosomes.")
    chromA = chromsA[0]
    chromB = chromsB[0]

    blocks = []
    CCp = list(ntx.weakly_connected_components(Gp))
    CCm = list(ntx.weakly_connected_components(Gm))
    
    for ccp in CCp:
        if len(ccp) >= nanosynteny_minsize:
            G_block_p = ntx.subgraph(Gp,ccp)
            block_ind = ntx.dag_longest_path(G_block_p)
            if len(block_ind) >= nanosynteny_minsize:
                block = condensed_dots[np.array(block_ind)]
                if check_for_nanosynteny_support:
                    supported = supported_by_nanosynteny(block, species_data_A, species_data_B, maps_A, maps_B, nanosynteny_minsize)
                else:
                    supported = True
                if supported:
                    col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
                    col2 = block[:,1:2].astype(int)
                    col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
                    col4 = block[:,3:4].astype(int)
                    blocks.append(np.hstack([col1,col2,col3,col4]))

    for ccm in CCm:
        if len(ccm) >= nanosynteny_minsize:
            G_block_m = ntx.subgraph(Gm,ccm)
            block_ind = ntx.dag_longest_path(G_block_m)
            if len(block_ind) >= nanosynteny_minsize:
                block = condensed_dots[np.array(block_ind)]
                if check_for_nanosynteny_support:
                    supported = supported_by_nanosynteny(block, species_data_A, species_data_B, maps_A, maps_B, nanosynteny_minsize)
                else:
                    supported = True
                if supported:
                    col1 = np.array(block.shape[0] * [chromA]).reshape([block.shape[0],1])
                    col2 = block[:,1:2].astype(int)
                    col3 = np.array(block.shape[0] * [chromB]).reshape([block.shape[0],1])
                    col4 = block[:,3:4].astype(int)
                    blocks.append(np.hstack([col1,col2,col3,col4]))

    return blocks

def get_block_extremes(block):
    x = block[:,1].astype(int)
    y = block[:,3].astype(int)
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(y)
    y_max = np.max(y)
    x_span = np.arange(x_min,x_max+1)
    y_span = np.arange(y_min,y_max+1)
    return x, y, x_min, x_max, y_min, y_max, x_span, y_span
    
def block_distance(blockA, blockB):
    if (blockA[0,0] != blockB[0,0]) or (blockA[0,2] != blockB[0,2]):
        # Different chromosomes
        return np.inf
    
    slopeA = np.sign(block_slope(blockA))
    slopeB = np.sign(block_slope(blockB))
    if slopeA != slopeB:
        # Incompatible slopes
        return np.inf
    else:
        slope = slopeA
    
    Ax, Ay, Ax_min, Ax_max, Ay_min, Ay_max, Ax_span, Ay_span = get_block_extremes(blockA)
    Bx, By, Bx_min, Bx_max, By_min, By_max, Bx_span, By_span = get_block_extremes(blockB)  
    
    x_overlap = len(set(Ax_span).intersection(set(Bx_span)))
    y_overlap = len(set(Ay_span).intersection(set(By_span)))
    if (x_overlap > 0) or (y_overlap > 0):
        # Overlap
        return np.inf
    
    if Ax_max < Bx_min:
        xdist = Bx_min - Ax_max
        interblock_x = [Ax_max, Bx_min]
        if Ay_max < By_min:    
            if slope > 0:
                ydist = By_min - Ay_max
                interblock_y = [Ay_max, By_min]
            elif slope < 0:
                # Slope and positioning mismatch
                return np.inf 
        elif Ay_max > By_min:
            if slope > 0:
                # Slope and positioning mismatch
                return np.inf
            elif slope < 0:
                ydist = Ay_min - By_max
                interblock_y = [By_max, Ay_min]
    else:
        # Wrong order
        return np.inf
    distance = np.max([xdist, ydist])  
    
    return distance, slope, interblock_x, interblock_y

def create_chromosome_pair_dot_dags(dots, max_distance):
    x = dots[:,1:2].astype(int)
    y = dots[:,3:4].astype(int)
    n = x.shape[0]
    o = np.ones([n,n])
    x_mat = x * o
    x_mat_T = x_mat.T
    y_mat = y * o
    y_mat_T = y_mat.T
    
    Ax = (x_mat < x_mat_T).astype(int)
    Aym = (y_mat > y_mat_T).astype(int)
    Am = Ax * Aym
    
    Ayp = (y_mat < y_mat_T).astype(int)
    Ap = Ax * Ayp

    Dx = np.abs((x_mat - x_mat_T).astype(int))
    Dy = np.abs((y_mat - y_mat_T).astype(int))
    D = ((Dx <= max_distance) * (Dy <= max_distance)).astype(int)

    Gm = ntx.DiGraph(Am * D)
    Gp = ntx.DiGraph(Ap * D)
    
    return Gm, Gp

def create_inter_nanosynteny_dot_dag(dots, G, slope, interblock_bounds, dots_not_in_blocks):
    interblock_x, interblock_y = interblock_bounds
    interx_min, interx_max = interblock_x
    intery_min, intery_max = interblock_y
    interx_condition = (dots[:,1].astype(int) > interx_min) * (dots[:,1].astype(int) < interx_max)
    intery_condition = (dots[:,3].astype(int) > intery_min) * (dots[:,3].astype(int) < intery_max)
    if slope < 0:
        block_i_condition = (dots[:,1].astype(int) == interx_min) * (dots[:,3].astype(int) == intery_max)
        block_j_condition = (dots[:,1].astype(int) == interx_max) * (dots[:,3].astype(int) == intery_min)
    else:
        block_i_condition = (dots[:,1].astype(int) == interx_min) * (dots[:,3].astype(int) == intery_min)
        block_j_condition = (dots[:,1].astype(int) == interx_max) * (dots[:,3].astype(int) == intery_max)       
    nodes = set(np.where(dots_not_in_blocks * interx_condition * intery_condition + block_i_condition + block_j_condition)[0])
    end_nodes = np.where(block_i_condition + block_j_condition)[0]
    G_inter = ntx.subgraph(G,nodes)
    
    return G_inter, end_nodes

def longest_dag_path(G, source, target):
    longest = ntx.dag_longest_path(G)
    if (source in longest) and (target in longest):
        return longest
    else:
        paths = [path for path in ntx.all_simple_paths(G,source=source,target=target) if (source in path) and (target in path)]
        if len(paths) > 0:
            longest = np.argmax([len(path) for path in paths])
            return paths[longest]
        else:
            return paths

def block_connection(G, end_nodes, connection_path = 'longest'):
    if connection_path not in ['longest', 'shortest']:
        raise ValueError("connection_path must be 'longest' or 'shortest'.")
    inter_cc_ij = list(ntx.weakly_connected_components(G))
    lower = end_nodes[0]
    upper = end_nodes[1]
    connected = False
    path = []
    for cc in inter_cc_ij:
        if (lower in cc) and (upper in cc):
            if connection_path == 'longest':
                path = longest_dag_path(G,lower,upper)
                if len(path) > 0:
                    connected = True
            elif connection_path == 'shortest':
                try:
                    path = ntx.shortest_path(G,source=lower,target=upper)
                    connected = True
                except ntx.NetworkXNoPath:
                    path = []
                    connected = False
            else:
                raise ValueError("connection_path must be 'longest' or 'shortest'.")
            break
    
    return connected, path

def find_block_distance_matrix(blocks, distance_cutoff):
    D_minus = np.ones([len(blocks),len(blocks)]) * np.inf
    D_plus = np.ones([len(blocks),len(blocks)]) * np.inf
    interblock_bounds = {}
    for i in range(len(blocks)):
        for j in range(len(blocks)):
            distance_info = block_distance(blocks[i], blocks[j])
            if type(distance_info) == float:
                if distance_info != np.inf:
                    raise ValueError("Something is wrong with the distance calculation.")
            else:
                d, slope, interblock_x, interblock_y = distance_info
                if slope < 0:
                    D_minus[i,j] = d
                elif slope > 0:
                    D_plus[i,j] = d
                interblock_bounds[(i,j)] = [interblock_x, interblock_y]
    DA_block_minus = np.nan_to_num(D_minus) * (D_minus <= distance_cutoff).astype(int)
    DA_block_plus = np.nan_to_num(D_plus) * (D_plus <= distance_cutoff).astype(int)
    return DA_block_minus, DA_block_plus, interblock_bounds

def find_dots_not_in_blocks(dots, blocks):
    blocks_stacked = np.vstack(blocks)[:,[1,3]].astype(int)
    dots_blocks = []
    for dot in dots:
        dot_in_block_x = (blocks_stacked[:,0] == dot[1].astype(int))
        dot_in_block_y = (blocks_stacked[:,1] == dot[3].astype(int))
        if np.sum(dot_in_block_x * dot_in_block_y) > 0:
            dots_blocks.append(True)
        else:
            dots_blocks.append(False)        
    dots_not_in_blocks = (1 - np.array(dots_blocks)).astype(bool)
    return dots_not_in_blocks

def connect_blocks(blocks, dots, max_distance, distance_cutoff):
    dots_not_in_blocks = find_dots_not_in_blocks(dots, blocks)
    DA_block_minus, DA_block_plus, interblock_bounds = find_block_distance_matrix(blocks, distance_cutoff)

    G_dots_minus, G_dots_plus = create_chromosome_pair_dot_dags(dots, max_distance)
    
    A_block_minus = np.zeros(DA_block_minus.shape)
    slope_minus = -1
    paths_minus = {}
    for i, node in enumerate(DA_block_minus):
        edges = np.where(node > 0)[0]
        if edges.shape[0] > 0:
            distances = DA_block_minus[i,edges]
            distances_argsort = np.argsort(distances)
            for j in edges[distances_argsort]: 
                G_ij, end_nodes = create_inter_nanosynteny_dot_dag(dots, G_dots_minus, slope_minus, interblock_bounds[(i,j)], dots_not_in_blocks)
                connected, path = block_connection(G_ij, end_nodes)
                if connected == True:
                    A_block_minus[i,j] = 1
                    paths_minus[(i,j)] = path
                    break
            
    A_block_plus = np.zeros(DA_block_plus.shape)
    slope_plus = 1
    paths_plus = {}
    for i, node in enumerate(DA_block_plus):
        edges = np.where(node > 0)[0]
        if edges.shape[0] > 0:
            distances = DA_block_plus[i,edges]
            distances_argsort = np.argsort(distances)
            for j in edges[distances_argsort]: 
                G_ij, end_nodes = create_inter_nanosynteny_dot_dag(dots, G_dots_plus, slope_plus, interblock_bounds[(i,j)], dots_not_in_blocks)
                connected, path = block_connection(G_ij, end_nodes)
                if connected == True:
                    A_block_plus[i,j] = 1
                    paths_plus[(i,j)] = path
                    break
    
    for j in np.where(A_block_minus.sum(0) > 1)[0]:
        edges = np.where(A_block_minus[:,j] == 1)[0]
        distances = DA_block_minus[edges,j]
        distances_argsort = np.argsort(distances)
        for n in distances_argsort[1:]:
            i = edges[n]
            A_block_minus[i,j] = 0
            del paths_minus[(i,j)]
            
    for j in np.where(A_block_plus.sum(0) > 1)[0]:
        edges = np.where(A_block_plus[:,j] == 1)[0]
        distances = DA_block_plus[edges,j]
        distances_argsort = np.argsort(distances)
        for n in distances_argsort[1:]:
            i = edges[n]
            A_block_plus[i,j] = 0
            del paths_plus[(i,j)]
    
    micro_blocks = {i:blocks[i] for i in range(len(blocks))}
    for i in range(A_block_minus.shape[0]):
        minus_edges = np.where(A_block_minus[i] == 1)[0]
        plus_edges = np.where(A_block_plus[i] == 1)[0]
        if (minus_edges.shape[0] == 0) and (plus_edges.shape[0] == 0):
            pass
        elif (minus_edges.shape[0] == 1) and (plus_edges.shape[0] == 0):
            j = minus_edges[0]
            path_ij = np.vstack([dots[k,:4] for k in paths_minus[(i,j)]])
            if path_ij.shape[0] < 2:
                raise ValueError("The path is too short, something is wrong.")
            if path_ij.shape[0] == 2:
                micro_blocks[j] = np.vstack([micro_blocks[i],micro_blocks[j]])
            else:
                micro_blocks[j] = np.vstack([micro_blocks[i],path_ij[1:-1],micro_blocks[j]])
            del micro_blocks[i] 
        elif (minus_edges.shape[0] == 0) and (plus_edges.shape[0] == 1):
            j = plus_edges[0]
            path_ij = np.vstack([dots[k,:4] for k in paths_plus[(i,j)]])
            if path_ij.shape[0] < 2:
                raise ValueError("The path is too short, something is wrong.")
            if path_ij.shape[0] == 2:
                micro_blocks[j] = np.vstack([micro_blocks[i],micro_blocks[j]])
            else:
                micro_blocks[j] = np.vstack([micro_blocks[i],path_ij[1:-1],micro_blocks[j]])
            del micro_blocks[i]            
        else:
            raise ValueError("A nanosynteny block shouldn't be connected to two others with different slopes!")
    
    return [micro_blocks[key] for key in micro_blocks.keys()]
