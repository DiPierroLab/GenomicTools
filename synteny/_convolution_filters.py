import numpy as np
from scipy.ndimage import convolve

def nanosynteny_convolve_dot_plot(M, x):
    x = int(x)
    k_plus = np.identity(x)[::-1]
    k_minus = np.identity(x)
    C_plus = convolve(M,k_plus,mode='constant',cval=0)[:-(x-1),:-(x-1)]
    Cp = (C_plus > x - 1).astype(int)
    C_minus = convolve(M,k_minus,mode='constant',cval=0)[:-(x-1),:-(x-1)]
    Cm = (C_minus > x - 1).astype(int)
    return Cp, Cm 

def nanosynteny_deconvolve_dot_plot(Cp, Cm, x):
    x = int(x)
    k_plus = np.identity(x)[::-1]
    k_minus = np.identity(x)
    kblocks_plus = np.vstack(np.where(Cp == 1)).T
    kblocks_minus = np.vstack(np.where(Cm == 1)).T
    M_shape = [Cm.shape[0]+x-1, Cm.shape[1]+x-1]
    Mp = np.zeros(M_shape)
    Mm = np.zeros(M_shape)
    for kbp in kblocks_plus:
        spot = kbp + np.vstack([np.arange(-1,x-1),np.arange(-1,x-1)]).T
        Mp[spot[::-1,0], spot[:,1]] += 1
    for kbm in kblocks_minus:
        spot = kbm + np.vstack([np.arange(-1,x-1),np.arange(-1,x-1)]).T
        Mm[spot[:,0], spot[:,1]] += 1
    Mr = (Mp + Mm > 0).astype(int)    
    return Mr

def maxdist_convolve_dot_plot(Mr, M, maxdist):
    k_ones = np.ones((maxdist+1,maxdist+1))
    I = convolve(Mr,k_ones)
    H = M * (I > 0).astype(int)
    return H

def nanosynteny_convolution_step_chrom_pair(shifted_dots, x):
    shifted_dots_A = shifted_dots[:,1].astype(int)
    shifted_dots_B = shifted_dots[:,3].astype(int)
    M = np.zeros([shifted_dots_A.max()+1,shifted_dots_B.max()+1])
    M[shifted_dots_A, shifted_dots_B] = 1
    Cp, Cm = nanosynteny_convolve_dot_plot(M, x)
    Mr = nanosynteny_deconvolve_dot_plot(Cp, Cm, x)
    shifted_dots_A, shifted_dots_B = np.where(Mr == 1)
    chromsA = np.array(shifted_dots_A.shape[0] * [shifted_dots[0,0]])
    chromsB = np.array(shifted_dots_B.shape[0] * [shifted_dots[0,2]])
    nanosynteny_convolved_dots = np.vstack([chromsA,shifted_dots_A+1,chromsB,shifted_dots_B+1]).T
    return nanosynteny_convolved_dots

def extension_convolution_step_chrom_pair(convolved_dots, shifted_dots, maxdist):
    shifted_dots_A = shifted_dots[:,1].astype(int)
    shifted_dots_B = shifted_dots[:,3].astype(int)
    convolved_dots_A = convolved_dots[:,1].astype(int)
    convolved_dots_B = convolved_dots[:,3].astype(int)
    M = np.zeros([shifted_dots_A.max()+1,shifted_dots_B.max()+1])
    M[shifted_dots_A, shifted_dots_B] = 1
    Mr = np.zeros([shifted_dots_A.max()+1,shifted_dots_B.max()+1])
    Mr[convolved_dots_A, convolved_dots_B] = 1
    H = maxdist_convolve_dot_plot(Mr, M, maxdist)
    extended_dots_A, extended_dots_B = np.where(H == 1)
    chromsA = np.array(extended_dots_A.shape[0] * [shifted_dots[0,0]])
    chromsB = np.array(extended_dots_B.shape[0] * [shifted_dots[0,2]])
    extended_dots = np.vstack([chromsA,extended_dots_A+1,chromsB,extended_dots_B+1]).T
    return extended_dots
