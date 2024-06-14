
import numpy as np
from scipy.signal import fftconvolve

def nanosynteny_convolve_dot_plot(M, minsize):
    """
    Finds all dots in a homology matrix M that are in nanosynteny blocks with at least minsize genes in them.

    Input:
    - M: numpy array (N X N), input homology matrix
    - minsize: int, the minimum number of genes allowed in a nanosynteny block

    Output:
    - Mr: numpy array (N X N), homology matrix encoding dots in nanosynteny
    """
    minsize = int(minsize)
    k_plus = np.identity(minsize)
    k_minus = k_plus[::-1]
    C_plus = np.round(fftconvolve(M,k_plus,mode='full'))
    C_minus = np.round(fftconvolve(M,k_minus,mode='full'))
    xp, yp = np.where(C_plus[(minsize-1):,(minsize-1):] == minsize)
    xm, ym = np.where(C_minus[(minsize-1):,(minsize-1):] == minsize)
    mr = np.zeros(M.shape)
    for n_plus in range(xp.shape[0]):
        mr[xp[n_plus]:(xp[n_plus]+minsize),yp[n_plus]:(yp[n_plus]+minsize)] += k_plus
    for n_minus in range(xm.shape[0]):
        mr[xm[n_minus]:(xm[n_minus]+minsize),ym[n_minus]:(ym[n_minus]+minsize)] += k_minus
    Mr = (mr > 0).astype(int)
    return Mr

def convolve_deconvolve_maxdist_dot_plot(Mr, M, maxdist):
    """
    Identify all dots in the homology matrix M that are within maxdist of any dot in the homology matrix Mr.
 
    Input:
    - Mr: numpy array (N X N), a homology matrix encoding active dots (usually those dots in synteny)
    - M: numpy array (N X N), a homology matrix 
    - maxdist: int, the maximum distance for adding dots. Distance is defined here as max(|x1-x2|,|y1-y2|) for dots (x1,y1) and (x2,y2).

    Output:
    - H: a homology matrix containing all dots in the homology matrix M that are within maxdist of any dot in the homology matrix Mr
    """
    maxdist = 2 * int(maxdist) + 1
    k_ones = np.ones((maxdist,maxdist))
    I = np.round(fftconvolve(Mr,k_ones,mode='same'))
    H = M * (np.abs(I) > 0).astype(int)
    return H
