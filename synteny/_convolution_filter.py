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

def convolve_deconvolve_maxdist_dot_plot(Mr, M, maxdist):

    k_ones = np.ones((maxdist+1,maxdist+1))
    I = convolve(Mr,k_ones)
    H = M * (I > 0).astype(int)

    return H
