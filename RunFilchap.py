#!/usr/bin/env python
#
import filchap3D_c18o_north as filchap
import numpy as np
from astropy.io import fits 
import os
import math

# read parameter file
params		= filchap.readParameters("./userDefined_3d_c18o_north.param")
beam_size       = params["beam_size"]
d		= params["distance"]
ra_dec_inc	= params["pixel_size"]

numberOfFilaments	= 600

for filamentNo in range(1,numberOfFilaments):
	try:
		plotsDir	= './plots/'
		resultsDir	= './results/'
		filamentData	= np.loadtxt('/home/suri/scripts/python/skeletonProcessing/filaments_north/filaments_north_old/averaged/filament' + str(filamentNo)+'.txt')
		print 'number of skeleton points: ', len(filamentData)
		# in case filament spines with certain lengths should be excluded 
		# here we calculate the length of the filament
		# we do not use the filaments that are shorter than 3 beamsizes
		ii = 0
		filamentLength = 0
		for ii in range(len(filamentData)-1):
			x_1	= filamentData[ii,0]
                	y_1	= filamentData[ii,1]
                	x_2 	= filamentData[ii+1,0]
                	y_2 	= filamentData[ii+1,1]
                	delta_x = (x_2 - x_1)**2
               		delta_y = (y_2 - y_1)**2
                	delta 	= math.sqrt(delta_x + delta_y)
                	theta 	= (delta*ra_dec_inc)/206265               # covert arcsec to radians
                	R 	= theta*d
                	filamentLength += R	
			
		# if the filament is longer than 3 beamsizes
		# we calculate it's properties 
		if filamentLength >= 3*beam_size: 
			results = filchap.calculateWidth(filamentNo,filamentData,filamentLength,params,plotIndividualProfiles=True,printResults=True, plotSaveLoc=plotsDir)
			# save properties
			np.savetxt('./results/filament'+ str(filamentNo) + '.txt', results)
	
	except (IOError,UnboundLocalError,IndexError,ZeroDivisionError) as e:
		print 'skipping filament', str(filamentNo)		
		pass



