#!/usr/bin/env python
#
import filchap
import numpy as np
from astropy.io import fits 

filamentNo = 1
filamentData = np.loadtxt('/home/suri/scripts/python/skeletonProcessing/filaments_north/averaged/filament' + str(filamentNo)+'.txt')

params = filchap.readParameters("./userDefined.param")
results = filchap.calculateWidth(filamentNo,filamentData,params,plotIndividualProfiles=True,printResults=True)
