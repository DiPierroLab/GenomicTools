
import numpy as np
from ._identify_blocks import *
from ._run_synteny_identification import *
from ._block_tools import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def identify_duplications(dot_plot, species_data, chrom_info, nanosynteny_minsize):
    sp = 'Species'
    blocks, maps_A, maps_B = run_basic_nanosynteny_identification(dot_plot, species_data, species_data, chrom_info, chrom_info, sp, sp, nanosynteny_minsize, unshift_blocks = False)
    string_blocks = []
    duplications = []
    palindromes = []
    for block in blocks:
        diagonal_block = self_diagonal_block(block)
        if not diagonal_block:
            bs, bsT, bsTf = block_to_string_relative(block)
            if (bsT not in string_blocks) and (bsTf not in string_blocks):
                palindrome = block_palindrome(block)
                if not palindrome:
                    duplications.append(block)
                else:
                    palindromes.append(block)
            string_blocks.append(bs)
    if len(duplications) > 0:
        duplications = unshift_synteny_blocks(duplications,maps_A,maps_B,nanosynteny_minsize)
    if len(palindromes) > 0:
        palindromes = unshift_synteny_blocks(palindromes,maps_A,maps_B,nanosynteny_minsize)
    return duplications, palindromes
