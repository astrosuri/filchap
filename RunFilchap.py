#!/usr/bin/env python
#
import filchap
import numpy as np
from astropy.io import fits 
import os
params = filchap.readParameters("./userDefined.param")

numberOfFilaments = 500

for filamentNo in range(numberOfFilaments):
	
	try:
		filamentData = np.loadtxt('/path/to/filament/coordinates/filament' + str(filamentNo)+'.txt')
	
		if len(filamentData) >= 12 and filamentData.size > 12: 
			results = filchap.calculateWidth(filamentNo,filamentData,params,plotIndividualProfiles=True,printResults=True)
			np.savetxt('/save/directory/results' + str(filamentNo)+'.txt',results)
	
	except (IOError,UnboundLocalError) as e:	
		pass
