#!/usr/bin/env python
#
import filchap
import numpy as np
from astropy.io import fits 
import os
params = filchap.readParameters("./userDefined.param")

for filamentNo in range(495):
	
	try:
		filamentData = np.loadtxt('/home/suri/scripts/python/skeletonProcessing/filaments_north/averaged/filament' + str(filamentNo)+'.txt')
		if not os.path.exists('/home/suri/development/filchap_1.0/c18o_north/width/c18owidth_north' + str(filamentNo) + '.txt'):
			if len(filamentData) >= 12 and filamentData.size > 12: 
				results = filchap.calculateWidth(filamentNo,filamentData,params,plotIndividualProfiles=True,printResults=True)
				print results
				np.savetxt('/home/suri/development/filchap_1.0/c18o_north/width/c18owidth_north' + str(filamentNo)+'.txt',results)
	
	except (IOError,UnboundLocalError) as e:
		
		pass
