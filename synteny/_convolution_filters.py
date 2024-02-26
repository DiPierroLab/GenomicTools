
import numpy as np
from scipy.signal import fftconvolve

#def nanosynteny_convolve_dot_plot(M, x):
#    x = int(x)
#    k_plus = np.identity(x)[::-1]
#    k_minus = np.identity(x)
#    C_plus = fftconvolve(M,k_plus,mode='same')[:-(x-1),:-(x-1)]
#    Cp = (C_plus > x - 1).astype(int)
#    C_minus = fftconvolve(M,k_minus,mode='same')[:-(x-1),:-(x-1)]
#    Cm = (C_minus > x - 1).astype(int)
#    return Cp, Cm

#def nanosynteny_deconvolve_dot_plot(Cp, Cm, x):
#    x = int(x)
#    kblocks_plus = np.vstack(np.where(Cp == 1)).T
#    kblocks_minus = np.vstack(np.where(Cm == 1)).T
#    M_shape = [Cm.shape[0]+x-1, Cm.shape[1]+x-1]
#    Mp = np.zeros(M_shape)
#    Mm = np.zeros(M_shape)
#    for kbp in kblocks_plus:
#        spot = kbp + np.vstack([np.arange(-1,x-1),np.arange(-1,x-1)]).T
#        Mp[spot[::-1,0], spot[:,1]] += 1
#    for kbm in kblocks_minus:
#        spot = kbm + np.vstack([np.arange(-1,x-1),np.arange(-1,x-1)]).T
#        Mm[spot[:,0], spot[:,1]] += 1
#    Mr = (Mp + Mm > 0).astype(int)
#    return Mr

def nanosynteny_convolve_dot_plot(M, x):
    x = int(x)
    k_plus = np.identity(x)
    k_minus = np.identity(x)[::-1]
    C_plus = fftconvolve(M,k_plus,mode='same')
    Cp = (C_plus > x - .5).astype(int)
    C_minus = fftconvolve(M,k_minus,mode='same')
    Cm = (C_minus > x - .5).astype(int)
    if x % 2 != 0:
        a = int((x-1)/2)
        b = int((x-1)/2)
    else:
        a = int(x/2)
        b = int(x/2)-1
    xp,yp = np.where(Cp==1)
    xp = np.hstack([xp+i for i in range(-a,b+1)])
    yp = np.hstack([yp+i for i in range(-a,b+1)])
    xm,ym = np.where(Cm==1)
    xm = np.hstack([xm+i for i in range(-a,b+1)])
    ym = np.hstack([ym-i-1 for i in range(-a,b+1)])
    Mr = np.zeros(M.shape)
    Mr[xp,yp] = 1
    Mr[xm,ym] = 1
    return Mr

def convolve_deconvolve_maxdist_dot_plot(Mr, M, maxdist):
    maxdist = 2 * int(maxdist) + 1
    k_ones = np.ones((maxdist,maxdist))
    I = fftconvolve(Mr,k_ones,mode='same')
    H = M * (np.abs(I) > .5).astype(int)
    return H
