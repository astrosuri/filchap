#!/usr/bin/env python
#
#this is an example script on how filchap can be run
#
import filchap2D as filchap
import numpy as np
from astropy.io import fits 
import os
params = filchap.readParameters("./userDefined.param")

for filamentNo in range(1,30):
  try:
	filamentData = np.loadtxt('/path/to/filament/coordinates/filament' + str(filamentNo) + '.txt')
	if len(filamentData) >= 12 and filamentData.size > 12: 
		results = filchap.calculateWidth(filamentNo,filamentData,params,plotIndividualProfiles=False,printResults=True)
		#print results
		np.savetxt('./results/results' + str(filamentNo) + '.txt',results)
	
  except (IOError,UnboundLocalError) as e:

	pass
