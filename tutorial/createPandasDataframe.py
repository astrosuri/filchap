#!/usr/bin/env python
#this script converts the width results from filchap
# txt to pandas dataframe in hdf5 format
#
import numpy as np
from astropy.io import fits
from astropy.table import Table
import pandas as pd
#
#

t		= Table()
listIX  	= []
i		= 1
numberOfFilaments = 1

while i <= numberOfFilaments:
	print i
	try:
		fitDir = np.loadtxt('./results/filament' + str(i) + '.txt')
		
	
		if (len(fitDir) == 32 and fitDir.size == 32):
			
			gausAmp = fitDir[0]
			gausAmpErr = fitDir[1]
			gausCen = fitDir[2]
			gausCenErr = fitDir[3]
			gausWidth = abs(fitDir[4])
                	gausWidthErr = fitDir[5]
			plum2Amp = fitDir[6]
			plum2AmpErr = fitDir[7]
			plum2Cen = fitDir[8]
			plum2CenErr = fitDir[9]
                	plum2Width = abs(fitDir[10])
        	       	plum2WidthErr = fitDir[11]
			plum4Amp = fitDir[12]
			plum4AmpErr = fitDir[13]
			plum4Cen = fitDir[14]
			plum4CenErr =fitDir[15]
               		plum4Width = abs(fitDir[16])
               		plum4WidthErr = fitDir[17]
                 	momentWid = fitDir[18]
                  	skewness = fitDir[19]
                  	kurtosis = fitDir[20]
			chi_sq_gaus = fitDir[21]
			red_chi_sq_gaus = fitDir[22]
			chi_sq_plum2 = fitDir[23]
			red_chi_sq_plum2= fitDir[24]
			chi_sq_plum4 = fitDir[25]
			red_chi_sq_plum4 = fitDir[26]
			rangePix = fitDir[27]	
			xPix = fitDir[28]
			yPix = fitDir[29]
			zPix = fitDir[30]
			nPeaks = fitDir[31]
			listIX.extend((i,gausAmp,gausAmpErr,gausCen,gausCenErr,gausWidth,gausWidthErr,plum2Amp,plum2AmpErr,plum2Cen,plum2CenErr,plum2Width,plum2WidthErr,plum4Amp,plum4AmpErr,plum4Cen,plum4CenErr,plum4Width,plum4WidthErr,momentWid,skewness,kurtosis,chi_sq_gaus,red_chi_sq_gaus,chi_sq_plum2,red_chi_sq_plum2,chi_sq_plum4,red_chi_sq_plum4,rangePix,xPix,yPix,zPix,nPeaks))
				
			
				
		elif  fitDir.size > 32:
			for jj in enumerate(np.arange(0,len(fitDir))):
				ii = jj[0]
				gausAmp = fitDir[ii,0]
                                gausAmpErr = fitDir[ii,1]
                                gausCen = fitDir[ii,2]
                                gausCenErr = fitDir[ii,3]
                                gausWidth = abs(fitDir[ii,4])
                                gausWidthErr = fitDir[ii,5]
                                plum2Amp = fitDir[ii,6]
                                plum2AmpErr = fitDir[ii,7]
                                plum2Cen = fitDir[ii,8]
                                plum2CenErr = fitDir[ii,9]
                                plum2Width = abs(fitDir[ii,10])
                                plum2WidthErr = fitDir[ii,11]
                                plum4Amp = fitDir[ii,12]
                                plum4AmpErr = fitDir[ii,13]
                                plum4Cen = fitDir[ii,14]
                                plum4CenErr =fitDir[ii,15]
                                plum4Width = abs(fitDir[ii,16])
                                plum4WidthErr = fitDir[ii,17]
                                momentWid = fitDir[ii,18]
                                skewness = fitDir[ii,19]
                                kurtosis = fitDir[ii,20]
			 	chi_sq_gaus = fitDir[ii,21]
                                red_chi_sq_gaus = fitDir[ii,22]
                                chi_sq_plum2 = fitDir[ii,23]
                                red_chi_sq_plum2= fitDir[ii,24]
                                chi_sq_plum4 = fitDir[ii,25]
                                red_chi_sq_plum4 = fitDir[ii,26]
                                rangePix = fitDir[ii,27]
				xPix = fitDir[ii,28]
                        	yPix = fitDir[ii,29]
                        	zPix = fitDir[ii,30]
				nPeaks = fitDir[ii,31]
				listIX.extend((i,gausAmp,gausAmpErr,gausCen,gausCenErr,gausWidth,gausWidthErr,plum2Amp,plum2AmpErr,plum2Cen,plum2CenErr,plum2Width,plum2WidthErr,plum4Amp,plum4AmpErr,plum4Cen,plum4CenErr,plum4Width,plum4WidthErr,momentWid,skewness,kurtosis,chi_sq_gaus,red_chi_sq_gaus,chi_sq_plum2,red_chi_sq_plum2,chi_sq_plum4,red_chi_sq_plum4,rangePix,xPix,yPix,zPix,nPeaks))	
			
			
		i += 1
	except IOError,IndexError:
		i += 1
		pass



IX = np.array(listIX,dtype=float).reshape(-1,33)


t['filament'] = IX[:,0]
t['gausAmp'] = IX[:,1]
t['gausAmpErr'] = IX[:,2]
t['gausCen'] = IX[:,3]
t['gausCenErr'] = IX[:,4]
t['gausWidth'] = IX[:,5]
t['gausWidthErr'] = IX[:,6]
t['plum2Amp'] = IX[:,7]
t['plum2AmpErr'] = IX[:,8]
t['plum2Cen'] = IX[:,9]
t['plum2CenErr'] = IX[:,10]
t['plum2Width'] = IX[:,11]
t['plum2WidthErr'] = IX[:,12]
t['plum4Amp'] = IX[:,13]
t['plum4AmpErr'] = IX[:,14]
t['plum4Cen'] = IX[:,15]
t['plum4CenErr'] = IX[:,16]
t['plum4Width'] = IX[:,17]
t['plum4WidthErr'] = IX[:,18]
t['momentWid'] = IX[:,19]
t['skewness'] = IX[:,20]
t['kurtosis'] = IX[:,21]
t['chi_sq_gaus'] = IX[:,22]
t['red_chi_sq_gaus'] = IX[:,23]
t['chi_sq_plum2'] = IX[:,24]
t['red_chi_sq_plum2'] = IX[:,25]
t['chi_sq_plum4'] = IX[:,26]
t['red_chi_sq_plum4'] = IX[:,27]
t['rangePix'] = IX[:,28]
t['xPix'] = IX[:,29]
t['yPix'] = IX[:,30]
t['zPix'] = IX[:,31]
t['nPeaks'] = IX[:,32]


df = t.to_pandas()

store = pd.HDFStore('./widthResults.h5')
store['df'] = df
