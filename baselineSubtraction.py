#!/usr/bin/env python
#
#
#
########################
# import ###############
########################
import math
#import pylab as pl
#from pylab import *
import sys 
import numpy as np
from operator import itemgetter, attrgetter
from astropy.io import fits
from scipy import asarray as ar,exp
#import pyfits
from scipy import ndimage, sparse
from scipy.sparse.linalg import spsolve
from scipy.signal import detrend

def baseline_als(y, lam, p, niter):
        '''
        Asymmetric Least Square algorithm.
        y       : signal
        lam     : smoothness, often takes values between 10**2 and 10**9
        p       : asymmetry value, often takes values between 0.001 and 0.1
        niter   : number of iterations 
        z       : simulated (smoothed) signal
        L       : length of y 
        D       : difference matrix
        w       : weights
        W       : diag(w) 

        np.eye  : Return a 2-D array with ones on the diagonal and zeros elsewhere
        '''
        L = len(y)
        D = sparse.csc_matrix(np.diff(np.eye(L), 2)) 
        w = np.ones(L)
        for i in xrange(niter):
                W = sparse.spdiags(w, 0, L, L)
                Z = W + lam * D.dot(D.transpose())
                z = spsolve(Z, w*y)
                w = p * (y > z) + (1-p) * (y < z)      # asymmetric
                #w = p * (y > z) + (p) * (y < z)       # symmetric
        return z


