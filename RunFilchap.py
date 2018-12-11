#!/usr/bin/env python
#
import filchap
import numpy as np
from astropy.io import fits 
import os
import math

# read parameter file
params		= filchap.readParameters("./userDefined.param")
beam_size       = params["beam_size"]


plotsDir	= './tutorial/plots/'
resultsDir	= './tutorial/results/'
filamentData	= np.loadtxt('./tutorial/filament.txt')
filamentNo	= 1		# in case you write a loop, filament no can be iterated over 	

filamentLength	= filchap.calculateLength(filamentData, params)
widthResults = filchap.calculateWidth(filamentNo,filamentData,filamentLength,params,plotIndividualProfiles=True,printResults=True, plotSaveLoc=plotsDir)
# save properties
np.savetxt(str(resultsDir) + '/filament'+ str(filamentNo) + '.txt', widthResults)
	




