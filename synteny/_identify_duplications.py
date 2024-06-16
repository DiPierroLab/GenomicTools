
import numpy as np
from ._identify_blocks import *
from ._run_synteny_identification import *
from GenomicTools.tools import *
from GenomicTools.tandem_duplications import *

def block_to_string_relative(block):
    b1 = np.char.add(block[:,0].astype(str),'-')
    b2 = np.char.add(b1,block[:,1].astype(str))
    b3 = np.char.add(block[:,2].astype(str),'-')
    b4 = np.char.add(b3,block[:,3].astype(str))            
    bb1 = np.char.add(b2,'-')
    bb2 = np.char.add(bb1,b4)
    bs = '_'.join(bb2)
    
    bb1T = np.char.add(b4,'-')
    bb2T = np.char.add(bb1T,b2)
    bsT = '_'.join(bb2T)   
    
    bb1Tf = np.char.add(b4[::-1].astype(str),'-')
    bb2Tf = np.char.add(bb1Tf,b2[::-1].astype(str))
    bsTf = '_'.join(bb2Tf)
    
    return bs, bsT, bsTf
    
def block_to_string_absolute(block):
    b1 = np.char.add(block[:,0].astype(str),'-')
    b2 = np.char.add(b1,block[:,1].astype(str))
    bs = '_'.join(b2)

    b1T = np.char.add(block[:,1].astype(str),'-')
    b2T = np.char.add(b1T,block[:,0].astype(str))
    bsT = '_'.join(b2T)

    b1Tf = np.char.add(block[::-1,1].astype(str),'-')
    b2Tf = np.char.add(b1Tf,block[::-1,0].astype(str))
    bsTf = '_'.join(b2Tf)
    
    return bs, bsT, bsTf

def blocks_symmetric(blocks):
    n_block_columns = blocks[0].shape[1]
    string_blocks = []
    transpose_string_blocks = []
    transpose_flip_string_blocks = []
    
    if n_block_columns == 4:
        for b in blocks:
            bs, bsT, bsTf = block_to_string_relative(b)
            string_blocks.append(bs)
            transpose_string_blocks.append(bsT)
            transpose_flip_string_blocks.append(bsTf)  
            
    elif n_block_columns == 2:    
        for b in blocks:
            bs, bsT, bsTf = block_to_string_absolute(b)
            string_blocks.append(bs)
            transpose_string_blocks.append(bsT)
            transpose_flip_string_blocks.append(bsTf)   
            
    symm = []
    for n in range(len(string_blocks)):
        transpose = (transpose_string_blocks[n] in string_blocks)
        transpose_flip = (transpose_flip_string_blocks[n] in string_blocks)
        symm.append(transpose or transpose_flip)
        
    return symm

def self_diagonal_block(block):
    n_block_columns = block.shape[1]
    if n_block_columns == 4:
        A0 = np.char.add(block[:,0],'-')
        A = np.char.add(A0,block[:,1])
        B0 = np.char.add(block[:,2],'-')
        B = np.char.add(B0,block[:,3])
    elif n_block_columns == 2:
        A = block[:,0]
        B = block[:,1]
    if np.all(A == B):
        return True
    else:
        return False

def block_palindrome(block):
    if block[0,0] == block[0,2]:
        indicesA = '-'.join(block[:,1]) 
        indicesB = '-'.join(block[::-1,3])
        if indicesA == indicesB:
            return True
        else:
            return False
    else:
        return False

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
