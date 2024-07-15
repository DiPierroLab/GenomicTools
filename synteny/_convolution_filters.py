
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

def get_nano_dots(dots, nanosynteny_minsize, chrom_info_A, chrom_info_B):
    chromsA = np.unique(dots[:,0])
    chromsB = np.unique(dots[:,2])
    if (chromsA.shape[0] != 1) or (chromsB.shape[0] != 1):
        raise ValueError("The input 'dots' should be all the dots between two chromosomes.")
    chromA = chromsA[0]
    chromB = chromsB[0]

    H = np.zeros((chrom_info_A[chromA]['size'],chrom_info_B[chromB]['size']))
    H[dots[:,1].astype(int)-1,dots[:,3].astype(int)-1] = 1
    
    H_nano = nanosynteny_convolve_dot_plot(H,nanosynteny_minsize)
    xr, yr = np.where(H_nano)
    nano_dots = np.vstack([np.array(xr.shape[0]*[chromA]),xr+1,np.array(yr.shape[0]*[chromB]),yr+1]).T
    
    return nano_dots

def get_nano_neighborhood_dots(dots, nanosynteny_minsize, distance_cutoff, chrom_info_A, chrom_info_B):
    chromsA = np.unique(dots[:,0])
    chromsB = np.unique(dots[:,2])
    if (chromsA.shape[0] != 1) or (chromsB.shape[0] != 1):
        raise ValueError("The input 'dots' should be all the dots between two chromosomes.")
    chromA = chromsA[0]
    chromB = chromsB[0]

    H = np.zeros((chrom_info_A[chromA]['size'],chrom_info_B[chromB]['size']))
    H[dots[:,1].astype(int)-1,dots[:,3].astype(int)-1] = 1
    
    H_nano = nanosynteny_convolve_dot_plot(H,nanosynteny_minsize)
    xr, yr = np.where(H_nano)
    nano_dots = np.vstack([np.array(xr.shape[0]*[chromA]),xr+1,np.array(yr.shape[0]*[chromB]),yr+1]).T
    
    H_nano_neighborhood = convolve_deconvolve_maxdist_dot_plot(H_nano, H, distance_cutoff)
    xn, yn = np.where(H_nano_neighborhood)
    nano_neighborhood_dots = np.vstack([np.array(xn.shape[0]*[chromA]),xn+1,np.array(yn.shape[0]*[chromB]),yn+1]).T
    
    return nano_dots, nano_neighborhood_dots
