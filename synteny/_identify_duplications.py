
import numpy as np
import re
from Bio import Phylo
from io import StringIO
from ._identify_blocks import *
from ._run_synteny_identification import *
from ._block_tools import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def identify_duplications(condensed_blocks, species_data, chrom_info, nanosynteny_minsize, unshift_blocks = False):
    maps = create_shift_map(species_data)
    string_blocks = []
    duplications = []
    palindromes = []
    palindromoids = []
    for block in condensed_blocks:
        diagonal_block = self_diagonal_block(block)
        if not diagonal_block:
            bs, bsT, bsTf = block_to_string_relative(block)
            if (bsT not in string_blocks) and (bsTf not in string_blocks):
                palindrome = block_palindrome(block)
                palindromoid = block_palindromoid(block)
                if not palindrome:
                    if not palindromoid:
                        duplications.append(block)
                    else:
                        palindromoids.append(block)
                else:
                    palindromes.append(block)
            string_blocks.append(bs)
    if unshift_blocks == True:
        if len(duplications) > 0:
            duplications = unshift_synteny_blocks(duplications,maps,maps,nanosynteny_minsize)
        if len(palindromes) > 0:
            palindromes = unshift_synteny_blocks(palindromes,maps,maps,nanosynteny_minsize)
        if len(palindromoids) > 0:
            palindromoids = unshift_synteny_blocks(palindromoids,maps,maps,nanosynteny_minsize)
    return duplications, palindromes, palindromoids

def span(indices):
    return np.max(indices) - np.min(indices) + 1

def condense_og_sequence(og_sequence):
    og_sequence_condensed = [og_sequence[0]]
    for og in og_sequence[1:]:
        if og != og_sequence_condensed[-1]:
            og_sequence_condensed.append(og)
    return og_sequence_condensed

def merge_duplicate_regions(absolute_duplications, species_data, filter_for_all_unique_ogs = True):
    duplication_set = set(range(len(absolute_duplications)))
    merge_map = {}

    if filter_for_all_unique_ogs == True:
        for i, dupA in enumerate(absolute_duplications):
            dupA = species_data[dupA[:,0]-1]
            dupA_ogs = dupA[:,4]
            dupA_str_forward = '-'.join(['x'+str(i)+'x' for i in dupA_ogs])
            dupA_str_backward = '-'.join(['x'+str(i)+'x' for i in dupA_ogs[::-1]])
            condensed_dupA_ogs = condense_og_sequence(dupA_ogs)
            if np.unique(condensed_dupA_ogs).shape[0] < len(condensed_dupA_ogs):
                duplication_set.remove(i)
                merge_map[i] = -1
                
    for i, dupA in enumerate(absolute_duplications):
        dupA = species_data[dupA[:,0]-1]
        dupA_ogs = dupA[:,4]
        dupA_str_forward = '-'.join(['x'+str(i)+'x' for i in dupA_ogs])
        dupA_str_backward = '-'.join(['x'+str(i)+'x' for i in dupA_ogs[::-1]])
        for j, dupB in enumerate(absolute_duplications):
            if (i == j) or (j not in duplication_set):
                continue
            else:
                dupB = species_data[dupB[:,0]-1]
                dupB_ogs = dupB[:,4]
                dupB_str = '-'.join(['x'+str(i)+'x' for i in dupB_ogs])
                match_forward = re.search(dupA_str_forward,dupB_str)
                match_backward = re.search(dupA_str_backward,dupB_str)
                
                if (match_forward is not None) or (match_backward is not None):
                    duplication_set.remove(i)
                    merge_map[i] = j
                    break
        if i in duplication_set:
            merge_map[i] = i
            
    merged = [absolute_duplications[d] for d in duplication_set]
    fixed_point_map = {}
    for i in merge_map.keys():
        j = i
        path = [j]
        while True:
            k = merge_map[j]
            if (k == j) or (k == -1):
                fixed_point_map[i] = k
                break
            elif k in path:
                fixed_point_map[i] = np.min(path)
                break
            else:
                j = k
                path.append(j)
    return merged, fixed_point_map

def select_best_regions(index_list):
    N = len(index_list)
    selected = []
    for i in range(N):
        found_intersection = False
        for s in range(len(selected)):
            I = list(index_list[i])
            S = list(index_list[s])
            if len(set(I).intersection(set(S))) > 0:
                if span(I) < span(S):
                    null = selected.pop(s)
                    selected.append(i)
                found_intersection = True
        if found_intersection == False:
            selected.append(i)
    return [index_list[s] for s in selected]

def select_best_regions_all(og_indices):
    og_indices_selected = {}
    for k in og_indices.keys():
        og_indices_selected[k] = select_best_regions(og_indices[k])
    return og_indices_selected
    
def merge_overlapping_regions(indicesA, indicesB): 
    lenA = len(indicesA)
    lenB = len(indicesB)
    if lenA <= lenB:
        indexA = 1
        indexB = lenA - lenB - 1
    else:
        indexA = lenA - lenB + 1
        indexB = -1
    if np.all(np.array(indicesA[indexA:]) == np.array(indicesB[:indexB])):
        return list(indicesA) + list(indicesB[indexB:])
    else:
        return []
            
def merge_all_overlapping(og_indices):
    og_indices = {k:og_indices[k] for k in og_indices.keys()}
    shifts = np.sort(list(og_indices.keys()))
    for s in shifts[:-1]:
        if s+1 in og_indices.keys():
            ns1 = len(og_indices[s])
            ns2 = len(og_indices[s+1])          
            remove_s1 = []
            remove_s2 = []
            add_s2 = []
            for i in range(ns1):
                for j in range(ns2):
                    overlap = merge_overlapping_regions(og_indices[s][i],og_indices[s+1][j])
                    if len(overlap) > 0:
                        remove_s1.append(i)
                        remove_s2.append(j)
                        add_s2.append(overlap)
                        break
            og_indices[s] = [list(og_indices[s][i]) for i in range(ns1) if i not in remove_s1]
            og_indices[s+1] = [list(og_indices[s+1][i]) for i in range(ns2) if i not in remove_s2]
            og_indices[s+1] += add_s2
                
    return og_indices

def match_to_gene_indices(og_sequence, match, ogs_chrom_str, species_ogs_chrom):
    ngenes = len(match.group(1).split('-'))
    start = match.start()
    end = match.start()+len(match.group(1))
    indices = []
    start = ogs_chrom_str[:start].count('-')
    end = ogs_chrom_str[:end].count('-')
    i = 0
    for og_chrom_index in range(start,end+1):
        og_chrom = species_ogs_chrom[og_chrom_index]
        if (i < len(og_sequence)) and (og_sequence[i] == og_chrom):
            indices.append(og_chrom_index)
            i += 1
    indices = np.array(indices)
    return indices, ngenes

def merge_all_overlapping_forward_backward(merged_forward, merged_backward, species_data_chrom, min_match_size):
    forward = []
    for k in merged_forward.keys():
         forward += merged_forward[k]
    backward = []
    for k in merged_backward.keys():
         backward += merged_backward[k]
    
    forward_backward = []
    for a in forward:
        a = np.array(a)
        if np.unique(species_data_chrom[a][:,4]).shape[0] >= min_match_size:
            forward_backward.append(a)
    for a in backward:
        a = np.array(a)
        if np.unique(species_data_chrom[::-1][a][:,4]).shape[0] >= min_match_size:
            forward_backward.append(species_data_chrom.shape[0] - 1 - a[::-1])

    N = len(forward_backward)
    edges = []
    for i in range(N):
        for j in range(N):
            len_i = len(forward_backward[i])
            len_j = len(forward_backward[j])
            for k in range(1,np.min([len_i,len_j])):
                overlap = np.all(np.array(forward_backward[i][::-1][:k][::-1]) == np.array(forward_backward[j][:k]))
                if overlap == True:
                    edges.append([i,j])
                    break
    
    G = ntx.DiGraph()
    G.add_nodes_from(np.arange(N))
    G.add_edges_from(edges)
    final_regions = []
    for cc in ntx.weakly_connected_components(G):
        Gcc = ntx.subgraph(G,cc)
        path = ntx.dag_longest_path(Gcc)
        final_regions.append(list(path))

    copies = []
    for r in final_regions:
        copies_r = []
        for i in r:
            copies_r.append(species_data_chrom[forward_backward[i]])
        copies.append(copies_r)

    return copies

def find_og_seq(og_sequence, species_data, window_size, min_match_size = 3):
    N = len(og_sequence) 
    og_sequence = [str(i) for i in og_sequence]
    og_sequence_str = ['x'+str(i)+'x' for i in og_sequence]
    
    search_strings = {}
    og_subsequences = {}
    for i in range(0,N-min_match_size+1):
        og_subsequences[i] = og_sequence[i:i+min_match_size]
        og_subsequence_str = og_sequence_str[i:i+min_match_size]
        search_str_init = '.+?'.join(og_subsequence_str)
        search_strings[i] = '(?=('+search_str_init+'))'
    
    chroms = np.unique(species_data[:,0])
    copies = []
    for chrom in chroms:    
        species_data_chrom = species_data[species_data[:,0] == chrom]
        N_genes_chrom = species_data_chrom.shape[0]
        species_ogs_chrom = species_data_chrom[:,4]
        ogs_chrom_forward = '-'.join(['x'+str(i)+'x' for i in species_ogs_chrom])
        ogs_chrom_backward = '-'.join(['x'+str(i)+'x' for i in species_ogs_chrom[::-1]])
        
        og_indices_forward = {}
        og_indices_backward = {}
        for i in range(0,N-min_match_size+1):
            search_str = search_strings[i]
            match_forward_all = list(re.finditer(search_str,ogs_chrom_forward))
            match_backward_all = list(re.finditer(search_str,ogs_chrom_backward))
            
            og_indices_forward[i] = []
            for match in match_forward_all:
                indices, ngenes = match_to_gene_indices(og_subsequences[i], match, ogs_chrom_forward, species_ogs_chrom)
                if ngenes <= window_size:
                    og_indices_forward[i].append(indices)
            if len(og_indices_forward[i]) == 0:
                del og_indices_forward[i]
            else:
                og_indices_forward[i] = select_best_regions(og_indices_forward[i])
            
            og_indices_backward[i] = []
            for match in match_backward_all:
                indices, ngenes = match_to_gene_indices(og_subsequences[i], match, ogs_chrom_backward, species_ogs_chrom[::-1])
                if ngenes <= window_size:
                    og_indices_backward[i].append(indices)
            if len(og_indices_backward[i]) == 0:
                del og_indices_backward[i]
            else:
                og_indices_backward[i] = select_best_regions(og_indices_backward[i])
        
        selected_forward = select_best_regions_all(og_indices_forward)
        merged_forward = merge_all_overlapping(selected_forward) 
        selected_backward = select_best_regions_all(og_indices_backward)
        merged_backward = merge_all_overlapping(selected_backward) 
        
        copies += merge_all_overlapping_forward_backward(merged_forward, merged_backward, species_data_chrom, min_match_size)

    return copies

def find_duplications_across_species_pair(absolute_duplications, species_data_source, species_data_target):
    duplication_matches = {}
    for i, abs_dup in enumerate(absolute_duplications):
        dup_og = condense_og_sequence(species_data_source[abs_dup[:,0] - 1][:,4])
        dup_matches = find_og_seq(dup_og, species_data_target, 2*len(dup_og))
        n_matches = len(dup_matches)
        duplication_matches[i] = {}
        duplication_matches[i]['n_matches'] = n_matches
        duplication_matches[i]['genes'] = dup_matches
    return duplication_matches

def find_duplications_across_all_species(sp_source, all_duplications, all_species_data):
    sp_data_source = all_species_data[sp_source]['species_data']
    chrom_info_source = all_species_data[sp_source]['chrom_info']
    absolute_duplications = convert_synteny_relative_to_absolute_indices(all_duplications[sp_source],chrom_info_source,chrom_info_source)
    merged_absolute_duplications, merge_map = merge_duplicate_regions(absolute_duplications,all_species_data[sp_source]['species_data'], filter_for_all_unique_ogs=True)
    duplications_across_all_species = {}
    sp_list = list(all_species_data.keys())
    i = 0
    for sp_target in sp_list:
        sp_data_target = all_species_data[sp_target]['species_data']
        data = find_duplications_across_species_pair(merged_absolute_duplications, sp_data_source, sp_data_target)
        duplications_across_all_species[sp_target] = data
        i += 1
        #print(i,end='\r',flush=True)
        
    return duplications_across_all_species

def get_species_time_to_mrca(spA, spB, tree, time_tree):
    time_tree_sp = set([term.name for term in time_tree.get_terminals()])
    C = tree.common_ancestor(spA,spB)
    a = C.name
    if len(C) == 0:
        if spA == spB:
            t = 0
        else:
            raise ValueError("Odd species tree behavior...")
    elif len(C) == 2:
        terms1 = set([term.name for term in C[0].get_terminals()])
        terms2 = set([term.name for term in C[1].get_terminals()])
        set1 = time_tree_sp.intersection(terms1)
        set2 = time_tree_sp.intersection(terms2)
        if (len(set1) == 0) or (len(set2) == 0):
            t = np.nan
        else:
            s1 = list(set1)[0]
            s2 = list(set2)[0]
            t = time_tree.distance(s1,s2)/2
    else:
        raise ValueError("Species tree is not binary...")
    return a, t

def get_time_to_ancestor(ancestor, tree, time_tree):
    time_tree_sp = set([term.name for term in time_tree.get_terminals()])
    C = list(tree.find_clades(ancestor))[0]
    a = C.name
    if len(C) == 0:
        t = 0
    elif len(C) == 2:
        terms1 = set([term.name for term in C[0].get_terminals()])
        terms2 = set([term.name for term in C[1].get_terminals()])
        set1 = time_tree_sp.intersection(terms1)
        set2 = time_tree_sp.intersection(terms2)
        if (len(set1) == 0) or (len(set2) == 0):
            t = np.nan
        else:
            s1 = list(set1)[0]
            s2 = list(set2)[0]
            t = time_tree.distance(s1,s2)/2
    else:
        raise ValueError("Species tree is not binary...")
    return a, t

def duplication_lower_bound(sp_source, duplications_across_all_species, all_species_data, all_duplications, gene_tree_data, map_dups_to_species, tree, time_tree):
    N = len(duplications_across_all_species[sp_source])
    absolute_duplications_sp = convert_synteny_relative_to_absolute_indices(all_duplications[sp_source],
                                                                            all_species_data[sp_source]['chrom_info'],
                                                                            all_species_data[sp_source]['chrom_info'])
    merged_duplications_sp, merge_map = merge_duplicate_regions(absolute_duplications_sp,all_species_data[sp_source]['species_data'])
    timings_sp = synteny_duplication_timing_with_gene_tree(merged_duplications_sp, sp_source,
                                                           all_species_data[sp_source]['species_data'],
                                                           all_species_data[sp_source]['chrom_info'],
                                                           gene_tree_data, tree, map_dups_to_species,
                                                           filter_earliest_supported=True)

    timings = {}
    for i in range(N):
        n_copies_i = {sp:duplications_across_all_species[sp][i]['n_matches'] for sp in duplications_across_all_species.keys()}
        dup_sp_dlcpar, dup_time_dlcpar = get_time_to_ancestor(timings_sp[i], tree, time_tree)
        dup_sp_synteny = sp_source
        dup_time_synteny = 0
        for sp in n_copies_i.keys():
            if n_copies_i[sp] > 1:
                ancestor, t = get_species_time_to_mrca(sp_source, sp, tree, time_tree)
                if t > dup_time_synteny:
                    dup_sp_synteny = ancestor
                    dup_time_synteny = t

        timings[i] = {}
        timings[i]['time_synteny'] = dup_time_synteny
        timings[i]['node_synteny'] = dup_sp_synteny
        timings[i]['time_DLCpar'] = dup_time_dlcpar
        timings[i]['node_DLCpar'] = dup_sp_dlcpar

    return timings

def match_duplicate_regions(duplicate_regions_i_source, duplicate_regions_i_target, sp_source, sp_target):
    matches = []
    for na, a in enumerate(duplicate_regions_i_source['genes']):
        if sp_source == sp_target:
            target_regions = duplicate_regions_i_target['genes'][na+1:]
        else:
            target_regions = duplicate_regions_i_target['genes']
        if len(target_regions) > 0:
            for b in target_regions:
                A = a[0]
                B = b[0]
                if (sp_source == sp_target) and np.all(A[0,:3] == B[0,:3]):
                    pass
                else:
                    matchA = np.vstack([A[A[:,4] == B[i,4]] for i in range(B.shape[0])])
                    matchB = np.vstack([B[B[:,4] == A[i,4]] for i in range(A.shape[0])])
                    matches.append(np.hstack([matchA[:,[0,2,5]],matchB[:,[0,2,5]]]))
    return matches

def match_duplicate_regions_across_species(duplications_all_species, timings, sp_source):
    region_indices = np.sort(list(duplications_all_species[sp_source].keys()))
    sp_list = list(duplications_all_species.keys())
    matched_duplications = {}
    for r in region_indices:
        matched_duplications[r] = {}
        matched_duplications[r]['timing'] = timings[r]
        matched_duplications[r]['matches'] = {}
        source_regions = duplications_all_species[sp_source][r]
        for sp_target in sp_list:
            target_regions = duplications_all_species[sp_target][r]
            if target_regions['n_matches'] == 0:
                continue
            else:
                matches = match_duplicate_regions(source_regions, target_regions, sp_source, sp_target)
                matched_duplications[r]['matches'][sp_target] = matches
    return matched_duplications

def parse_tree(tree_string):
    tree = Phylo.read(StringIO(tree_string),'newick')
    return tree

def gene_coalescence_species(geneA, geneB, og, gene_tree, map_dups_to_species, species_tree):
    gene_mrca = gene_tree.common_ancestor(geneA,geneB) 
    key = (og, gene_mrca.name)
    if key in map_dups_to_species.keys():
        gene_mrca_sp, support = map_dups_to_species[key]
        support = float(support)
    else:
        gene_mrca_sp = 'Undetermined'
        support = np.nan
    return gene_mrca_sp, support

def earliest_supported(syntenic_duplication, support_threshold = .5):
    names = syntenic_duplication[:,0]
    supports = np.nan_to_num(syntenic_duplication[:,1].astype(float))
    well_supported = (supports >= support_threshold)
    if np.sum(well_supported) == 0:
        support_threshold = np.max(supports)
        well_supported = (supports >= support_threshold)
    well_supported_names = np.unique(names[well_supported])
    if 'Undetermined' in well_supported_names:
        well_supported_names.remove('Undetermined')
    terminals = [s for s in well_supported_names if '_' in s]
    non_terminals = list(set(well_supported_names) - set(terminals))
    if len(non_terminals) == 0:
        if len(terminals) > 1:
            raise ValueError("Something weird is going on...")
        else:
            name = terminals[0]
    else:
        oldest_num = np.min([int(s.lstrip('N')) for s in non_terminals])
        name = 'N' + str(int(oldest_num))
    return name

def synteny_duplication_timing_with_gene_tree(duplications, sp, species_data, chrom_info, gene_trees, species_tree, map_dups_to_species, filter_earliest_supported = True):
    if duplications[0].shape[1] == 4:
        absolute_duplications = convert_synteny_relative_to_absolute_indices(duplications,chrom_info,chrom_info)
    elif duplications[0].shape[1] == 2:
        absolute_duplications = duplications
    else:
        raise ValueError("The format of the duplications input is not supported.")
    gene_tree_duplications_timings = []
    for duplication in absolute_duplications:
        genesA = species_data[duplication[:,0]-1]
        genesB = species_data[duplication[:,1]-1]
        gene_tree_duplication_timings = []
        for n in range(duplication.shape[0]):
            og_genes = 'OG'+str(int(genesA[n,4])-1).zfill(7)
            if og_genes in gene_trees.keys():
                gene_tree = parse_tree(gene_trees[og_genes])
                geneA = genesA[n,5]
                for gene in gene_tree.get_terminals():
                    if geneA in gene.name:
                        geneA = gene.name
                        break
                geneB = genesB[n,5]
                for gene in gene_tree.get_terminals():
                    if geneB in gene.name:
                        geneB = gene.name
                        break
                timing, support = gene_coalescence_species(geneA,geneB,og_genes,gene_tree,map_dups_to_species,species_tree)
                gene_tree_duplication_timings.append([timing.split('.')[0],support])
            else:
                gene_tree_duplication_timings.append([sp, 1.0])
        timings_array = np.vstack(gene_tree_duplication_timings)
        if filter_earliest_supported == True:
            gene_tree_duplications_timings.append(earliest_supported(timings_array))
        else:
            gene_tree_duplications_timings.append(timings_array)
    return gene_tree_duplications_timings
