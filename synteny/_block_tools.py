
import numpy as np

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
