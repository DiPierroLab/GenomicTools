
import numpy as np
from ._identify_blocks import *
from ._run_synteny_identification import *
from ._block_tools import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def identify_duplications(condensed_dots, species_data, chrom_info, nanosynteny_minsize, unshift_blocks = False):
    maps = create_shift_map(species_data)
    blocks = find_nanosynteny(condensed_dots, species_data, species_data, chrom_info, chrom_info, maps, maps, nanosynteny_minsize)
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
