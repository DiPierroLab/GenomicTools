
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

def get_species_time_to_ancestor(ancestor, tree, time_tree):
    time_tree_sp = set([term.name for term in time_tree.get_terminals()])
    C = list(tree.find_clades(ancestor))[0]
    name = C.name
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
    return name, t

def get_single_node_time_to_mrca(sp, tree, time_tree):
    time_tree_sp = set([term.name for term in time_tree.get_terminals()])
    C = list(tree.find_clades(sp))[0]
    a = C.name

    if len(C) == 0:
        t = 0
    else:
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
    return t

def condense_og_sequence(og_sequence):
    og_sequence_condensed = [og_sequence[0]]
    for og in og_sequence[1:]:
        if og != og_sequence_condensed[-1]:
            og_sequence_condensed.append(og)
    return og_sequence_condensed

def generate_node_to_time_map(tree, time_tree):
    node_to_time_map = {}
    for node in tree.get_nonterminals() + tree.get_terminals():
        node = node.name
        node_to_time_map[node] = {}
        if node == 'N0':
            time = get_species_time_to_ancestor(node,tree,time_tree)[1]
            node_to_time_map[node]['parent'] = 'N0'
            node_to_time_map[node]['avg_time'] = time
            node_to_time_map[node]['avg_sps_from_root'] = 0
        else:
            parent_node = ([tree.root]+tree.root.get_path(node))[-2].name
            time = get_species_time_to_ancestor(node,tree,time_tree)[1]
            time_parent = get_species_time_to_ancestor(parent_node,tree,time_tree)[1]
            time = (time + time_parent) / 2.
            sps_from_root = tree.distance('N0',node)
            sps_from_root_parent = tree.distance('N0',parent_node)
            sps_from_root = (sps_from_root + sps_from_root_parent) / 2.
            node_to_time_map[node]['parent'] = parent_node
            node_to_time_map[node]['avg_time'] = time
            node_to_time_map[node]['avg_sps_from_root'] = sps_from_root
    return node_to_time_map

def duplication_timing_gene_pair_DLCpar(geneA, geneB, og, gene_tree, map_dups_to_species):
    gene_mrca = gene_tree.common_ancestor(geneA,geneB)
    key = (og, gene_mrca.name)
    if key in map_dups_to_species.keys():
        gene_mrca_sp, support = map_dups_to_species[key]
        support = float(support)
    else:
        gene_mrca_sp = 'Undetermined'
        support = np.nan
    return gene_mrca_sp, support

def duplication_timing_DLCpar(condensed_duplication, sp, all_species_data, gene_tree_data, map_dups_to_species, node_to_time_map):
    maps_sp = all_species_data[sp]['maps']
    chrA = condensed_duplication[0,0]
    chrB = condensed_duplication[0,2]
    sp_data_chrA = all_species_data[sp]['species_data'][all_species_data[sp]['species_data'][:,0] == chrA]
    sp_data_chrB = all_species_data[sp]['species_data'][all_species_data[sp]['species_data'][:,0] == chrB]
    condensed_gene_to_geneA = lambda x: maps_sp[1][chrA][maps_sp[3][chrA][x]]
    condensed_gene_to_geneB = lambda x: maps_sp[1][chrB][maps_sp[3][chrB][x]]
    duplication_node_candidates = []
    for gene_pair in condensed_duplication:
        condensed_geneA = gene_pair[1]
        condensed_geneB = gene_pair[3]
        genesA = condensed_gene_to_geneA(int(condensed_geneA))
        genesB = condensed_gene_to_geneB(int(condensed_geneB))
        og_str = 'OG'+str(int(sp_data_chrA[sp_data_chrA[:,2] == str(genesA[0])][0,4])-1).zfill(7)
        if og_str in gene_tree_data.keys():
            gene_tree = parse_tree(gene_tree_data[og_str])
            for geneA in genesA:
                gene_nameA = sp_data_chrA[sp_data_chrA[:,2] == str(geneA)][0,5]
                for geneB in genesB:
                    gene_nameB = sp_data_chrB[sp_data_chrB[:,2] == str(geneB)][0,5]
                    for gene in gene_tree.get_terminals():
                        if gene_nameA in gene.name:
                            gene_nameA = gene.name
                            break
                    for gene in gene_tree.get_terminals():
                        if gene_nameB in gene.name:
                            gene_nameB = gene.name
                            break
                    gene_mrca_sp, support = duplication_timing_gene_pair_DLCpar(gene_nameA, gene_nameB, og_str, gene_tree, map_dups_to_species)
                    gene_mrca_sp = gene_mrca_sp.split('.')[0]
                    if not np.isnan(support):
                        sps_to_root = node_to_time_map[gene_mrca_sp]['avg_sps_from_root']
                        time = node_to_time_map[gene_mrca_sp]['avg_time']
                        duplication_node_candidates.append([gene_mrca_sp,sps_to_root,time,support])
    duplication_node_candidates = np.vstack(duplication_node_candidates)
    supported_candidates = duplication_node_candidates[duplication_node_candidates[:,3].astype(float) >= .5]
    if supported_candidates.shape[0] > 0:
        earliest_supported = supported_candidates[supported_candidates[:,1].astype(float) == supported_candidates[:,1].astype(float).min()]
        DLCpar_node, DLCpar_sps_from_root, DLCpar_time, DLCpar_support = earliest_supported[earliest_supported[:,3].astype(float).argmax()]
    else:
        earliest_nonsupported = duplication_node_candidates[duplication_node_candidates[:,1].astype(float) == duplication_node_candidates[:,1].astype(float).min()]
        DLCpar_node, DLCpar_sps_from_root, DLCpar_time, DLCpar_support = earliest_nonsupported[earliest_nonsupported[:,3].astype(float).argmax()]
    return DLCpar_node, DLCpar_sps_from_root, DLCpar_time, DLCpar_support

def duplication_timing_parsimony(condensed_duplication, sp, all_species_data, species_tree, node_to_time_map):
    sp_list = list(all_species_data.keys())
    species_data = all_species_data[sp]['species_data']
    chrom = condensed_duplication[0,0]
    species_data_chrom = species_data[species_data[:,0] == chrom]
    dup_ogs = [species_data_chrom[species_data_chrom[:,-2] == tdec_index][0,4] for tdec_index in condensed_duplication[:,1]]
    dup_og_str_forward = '-'.join(dup_ogs)
    dup_og_str_backward = '-'.join(dup_ogs[::-1])
    
    multiple_copy_species = []
    for spB in sp_list:
        sp_data = all_species_data[spB]['species_data']
        chr_info = all_species_data[spB]['chrom_info']
        chroms = list(chr_info.keys())
        count = 0
        for chrom in chroms:
            ogs_chrom = condense_og_sequence(sp_data[sp_data[:,0] == chrom][:,4])
            ogs_chrom_str = '-'.join(ogs_chrom)
            forward = re.findall(dup_og_str_forward,ogs_chrom_str)
            backward = re.findall(dup_og_str_backward,ogs_chrom_str)
            count += len(forward) + len(backward)
            if count > 1:
                multiple_copy_species.append(spB)
                break
     
    ancestor = species_tree.common_ancestor(multiple_copy_species)
    return ancestor.name

def nanosynteny_duplication_timing(condensed_duplications, sp, all_species_data, gene_tree_data, species_tree, map_dups_to_species, node_to_time_map):
    nanosynteny_duplication_timing = {}    
    for dn, dup_block in enumerate(condensed_duplications):
        DLCpar_node, DLCpar_sps_from_root, DLCpar_time, DLCpar_support = duplication_timing_DLCpar(dup_block, sp, all_species_data, gene_tree_data, map_dups_to_species, node_to_time_map)
        parsimony_node = duplication_timing_parsimony(dup_block, sp, all_species_data, species_tree, node_to_time_map)
        parsimony_time = node_to_time_map[parsimony_node]['avg_time']
        parsimony_sps_from_root = node_to_time_map[parsimony_node]['avg_sps_from_root']
        
        if (float(DLCpar_sps_from_root) <= float(parsimony_sps_from_root)):
            consensus_node = DLCpar_node
            consensus_time = float(DLCpar_time)
            consensus_sps_from_root = float(DLCpar_sps_from_root)
        else:
            consensus_node = parsimony_node
            consensus_time = float(parsimony_time)
            consensus_sps_from_root = float(parsimony_sps_from_root)
           
        nanosynteny_duplication_timing[dn] = {}
        nanosynteny_duplication_timing[dn]['node'] = consensus_node
        nanosynteny_duplication_timing[dn]['time'] = consensus_time
        nanosynteny_duplication_timing[dn]['sps_from_root'] = consensus_sps_from_root
        nanosynteny_duplication_timing[dn]['condensed_genes'] = dup_block
            
    return nanosynteny_duplication_timing

def microsynteny_duplication_timing(nanosynteny_duplication_timing, condensed_self_microsynteny):
    self_micro = [a for a in condensed_self_microsynteny if (a[0,0] != a[0,2]) or (a[0,1] != a[0,3])]
    nano_dups_indices = np.sort(list(nanosynteny_duplication_timing.keys()))
    nano_dups_set = set(nano_dups_indices)
    links = []
    where_nano = {md:[] for md in range(len(self_micro))}
    for nd1 in nano_dups_indices:
        nano_dup1 = nanosynteny_duplication_timing[nd1]['condensed_genes']
        nano_genes1A = np.char.add(np.char.add(nano_dup1[:,0],'-'),nano_dup1[:,1])
        nano_genes1B = np.char.add(np.char.add(nano_dup1[:,2],'-'),nano_dup1[:,3])
        nano_genes1 = np.char.add(np.char.add(nano_genes1A,'-'),nano_genes1B)
        for nd2 in nano_dups_indices:
            nano_dup2 = nanosynteny_duplication_timing[nd2]['condensed_genes']
            nano_genes2A = np.char.add(np.char.add(nano_dup2[:,0],'-'),nano_dup2[:,1])
            nano_genes2B = np.char.add(np.char.add(nano_dup2[:,2],'-'),nano_dup2[:,3])
            nano_genes2 = np.char.add(np.char.add(nano_genes2A,'-'),nano_genes2B)
            if (nd1 in nano_dups_set) and (nd2 in nano_dups_set) and (nd1 < nd2):
                for md, micro_dup in enumerate(self_micro):
                    micro_genesA = np.char.add(np.char.add(micro_dup[:,0],'-'),micro_dup[:,1])
                    micro_genesB = np.char.add(np.char.add(micro_dup[:,2],'-'),micro_dup[:,3])
                    micro_genes = np.char.add(np.char.add(micro_genesA,'-'),micro_genesB)
                    nano1_in_micro = set(nano_genes1).issubset(set(micro_genes))
                    nano2_in_micro = set(nano_genes2).issubset(set(micro_genes))
                    if nano1_in_micro and nano2_in_micro:
                        links.append([nd1,nd2,md])
                        where_nano_link = np.zeros(micro_genes.shape[0])
                        for g1 in nano_genes1:
                            where_nano_link[np.where(micro_genes == g1)[0]] = 1
                        for g2 in nano_genes2:
                            where_nano_link[np.where(micro_genes == g2)[0]] = 1
                        where_nano[md].append(where_nano_link)
                        break
    
    nano_indices = set(nano_dups_indices)
    if len(links) > 0:
        links = np.vstack(links)    

    i = 0
    microsynteny_duplication_timing = {}
    micro_dup = [i]
    while len(nano_indices) > 0:
        if (len(links) > 0) and (i in links[:,0]) and (i in nano_indices):
            j = links[links[:,0] == i][:,1].min()
            micro_dup.append(j)
            nano_indices = nano_indices - {i}
            i = j
        else:
            microsynteny_duplication_timing[tuple(micro_dup)] = {}
            time = 0
            sps_from_root = .4
            node = 'NA'
            for d in micro_dup:
                dtime = float(nanosynteny_duplication_timing[d]['time'])
                dsps = float(nanosynteny_duplication_timing[d]['sps_from_root'])
                if  dtime >= time:
                    time = dtime
                    sps_from_root = dsps
                    node = nanosynteny_duplication_timing[d]['node']
            if len(micro_dup) > 1:
                md = links[links[:,0] == micro_dup[0]][0,2]
                genes = self_micro[md]
                is_in_nanosynteny = np.any(np.vstack(where_nano[md]) == 1, axis=0)
            else:
                genes = nanosynteny_duplication_timing[micro_dup[0]]['condensed_genes']
                is_in_nanosynteny = np.ones(genes.shape[0]).astype(bool)

            microsynteny_duplication_timing[tuple(micro_dup)]['node'] = node
            microsynteny_duplication_timing[tuple(micro_dup)]['time'] = time
            microsynteny_duplication_timing[tuple(micro_dup)]['sps_from_root'] = sps_from_root
            microsynteny_duplication_timing[tuple(micro_dup)]['condensed_genes'] = genes
            microsynteny_duplication_timing[tuple(micro_dup)]['is_in_nanosynteny'] = is_in_nanosynteny
            microsynteny_duplication_timing[tuple(micro_dup)]['is_palindromoid'] = block_palindromoid(genes)
            
            nano_indices = nano_indices - {i}
            if len(nano_indices) > 0:
                i = np.min(list(nano_indices))
            else:
                break
            micro_dup = [i]

    return microsynteny_duplication_timing
