
import numpy as np
import re
from ._identify_blocks import *
from ._run_synteny_identification import *
from ._block_tools import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def identify_duplications(blocks, species_data, chrom_info, nanosynteny_minsize, unshift_blocks = False):
    maps = create_shift_map(species_data)
    blocks = shift_synteny_blocks(blocks, species_data, species_data, maps, maps)
    string_blocks = []
    duplications = []
    palindromes = []
    palindromoids = []
    for block in blocks:
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

def refine_matches(search_str, matches, ogs_chrom_str):
    search_str_refine = '(?=('+search_str+'))'
    matches_refined = []
    for match in matches:
        refined_matches = re.finditer(search_str_refine,match.group())
        m_min = ''
        m_len_min = np.inf
        for m in refined_matches:
            m_start = m.start()
            m_group = m.group(1)
            m_len = len(m_group.split('-'))
            if m_len < m_len_min:
                m_len_min = m_len
                m_min = m_group
        matches_refined += list(re.finditer(m_min,ogs_chrom_str))
    return matches_refined

def find_og_sequence(og_sequence, species_data, window_size):
    og_sequence = [str(i) for i in og_sequence]
    og_sequence_str = ['x'+str(i)+'x' for i in og_sequence]
    search_str = '.+?'.join(og_sequence_str)
    chroms = np.unique(species_data[:,0])
    matches_forward = []
    matches_backward = []
    for chrom in chroms:
        species_data_chrom = species_data[species_data[:,0] == chrom]
        species_ogs_chrom = species_data_chrom[:,4]
        ogs_chrom_str_forward = '-'.join(['x'+str(i)+'x' for i in species_ogs_chrom])
        match_forward_all = list(re.finditer(search_str,ogs_chrom_str_forward))
        ogs_chrom_str_backward = '-'.join(['x'+str(i)+'x' for i in species_ogs_chrom[::-1]])
        match_backward_all = list(re.finditer(search_str,ogs_chrom_str_backward))
        
        match_forward_all = refine_matches(search_str, match_forward_all, ogs_chrom_str_forward)
        match_backward_all = refine_matches(search_str, match_backward_all, ogs_chrom_str_backward)
        
        for match_forward in match_forward_all:
            ngenes = len(match_forward.group().split('-'))
            if ngenes <= window_size:
                og_indices = []                   
                start = ogs_chrom_str_forward[:match_forward.start()].count('-')
                end = ogs_chrom_str_forward[:match_forward.end()].count('-')
                i = 0
                for og_chrom_index in range(start,end+1):
                    og_chrom = species_ogs_chrom[og_chrom_index]
                    if (i < len(og_sequence)) and (og_sequence[i] == og_chrom):
                        og_indices.append(og_chrom_index)
                        i += 1
                og_indices = np.array(og_indices)
                og_match = species_ogs_chrom[og_indices]
                matches_forward.append([chrom,og_indices,og_match])

        for match_backward in match_backward_all:
            ngenes = len(match_backward.group().split('-'))
            if ngenes <= window_size:
                og_indices = []                   
                start = ogs_chrom_str_backward[:match_backward.start()].count('-')
                end = ogs_chrom_str_backward[:match_backward.end()].count('-')
                i = 0
                for og_chrom_index in range(start,end+1):
                    og_chrom = species_ogs_chrom[::-1][og_chrom_index]
                    if (i < len(og_sequence)) and (og_sequence[i] == og_chrom):
                        og_indices.append(og_chrom_index)
                        i += 1
                og_indices = np.array(og_indices)
                og_match = species_ogs_chrom[::-1][og_indices]
                matches_backward.append([chrom,og_indices,og_match])
            
    return matches_forward, matches_backward

def og_subsequence_scan(og_sequence, species_data, window_size, min_match_size = 3):
    all_matches_forward = {}
    all_matches_backward = {}
    sequence_matches_forward, sequence_matches_backward = find_og_sequence(og_sequence, species_data, window_size)
    
    for n, match_forward in enumerate(sequence_matches_forward):
        chrom, match, og_match = match_forward
        if chrom not in all_matches_forward.keys():
            all_matches_forward[chrom] = {}
        all_matches_forward[chrom][n] = match
        
    for n, match_backward in enumerate(sequence_matches_backward):
        chrom, match, og_match = match_backward
        if chrom not in all_matches_backward.keys():
            all_matches_backward[chrom] = {}
        all_matches_backward[chrom][n] = match
        
    for i in range(len(og_sequence)-min_match_size+1):
        subsequence_matches_forward, subsequence_matches_backward = find_og_sequence(og_sequence[i:i+min_match_size], species_data, window_size)
        
        if len(subsequence_matches_forward) > 0:
            for match_forward in subsequence_matches_forward:
                chrom, match, og_match = match_forward

                if chrom in all_matches_forward.keys():
                    new_match = True
                    for key in all_matches_forward[chrom].keys():
                        key_indices = set(all_matches_forward[chrom][key])
                        combine_candidate = key_indices.intersection(set(match))
                        if len(combine_candidate) >=  min_match_size - 1:
                            all_matches_forward[chrom][key] = key_indices.union(set(match))
                            new_match = False
                    if new_match == True:
                        new_key = np.max(list(all_matches_forward[chrom].keys())) + 1
                        all_matches_forward[chrom][new_key] = set(match)
                else:
                    all_matches_forward[chrom] = {}
                    all_matches_forward[chrom][0] = set(match)

        if len(subsequence_matches_backward) > 0:
            for match_backward in subsequence_matches_backward:
                chrom, match, og_match = match_backward

                if chrom in all_matches_backward.keys():
                    new_match = True
                    for key in all_matches_backward[chrom].keys():
                        key_indices = set(all_matches_backward[chrom][key])
                        combine_candidate = key_indices.intersection(set(match))
                        if len(combine_candidate) >= min_match_size - 1:
                            all_matches_backward[chrom][key] = key_indices.union(set(match))
                            new_match = False
                    if new_match == True:
                        new_key = np.max(list(all_matches_backward[chrom].keys())) + 1
                        all_matches_backward[chrom][new_key] = set(match)
                else:
                    all_matches_backward[chrom] = {}
                    all_matches_backward[chrom][0] = set(match)
    
    chroms = list(set(all_matches_forward.keys()).union(all_matches_backward.keys()))
    matches = []
    for chrom in chroms:
        species_data_chrom = species_data[species_data[:,0] == chrom]
        
        if chrom in all_matches_forward.keys():
            for key in all_matches_forward[chrom].keys():
                og_indices = np.sort(list(all_matches_forward[chrom][key]))
                match_genes = species_data_chrom[og_indices]
                match_gene_ogs = np.unique(match_genes[:,4])
                if match_gene_ogs.shape[0] >= min_match_size:
                    matches.append(match_genes)
                
        if chrom in all_matches_backward.keys():
            for key in all_matches_backward[chrom].keys():
                og_indices = np.sort(list(all_matches_backward[chrom][key]))
                match_genes = species_data_chrom[::-1][og_indices]
                match_gene_ogs = np.unique(match_genes[:,4])
                if match_gene_ogs.shape[0] >= min_match_size:
                    matches.append(match_genes)
                
    return matches
