#!/usr/bin/env python
#
#import packages
#

from __future__ import division
import math
import matplotlib
matplotlib.use('Agg')		# use this when running on screen
#matplotlib.use('Qt4Agg')	
import matplotlib.pyplot as pl
#pl.ion()
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from baselineSubtraction import baseline_als
from scipy.ndimage import gaussian_filter
from scipy.signal import argrelextrema
from scipy import stats 
import aplpy
import matplotlib.pyplot as plt
from math import sqrt
import warnings
warnings.filterwarnings("ignore")


matplotlib.rcParams.update({'font.size': 14.})
matplotlib.rcParams.update({'font.family':'serif'})
#

def gaus(arr,a,mu,sigma):
	return a*np.exp(-(arr-mu)**2/(2*sigma**2))

def plum2(arr,a0,x0,rflat):
	return (a0*rflat)/(1.+((arr-x0)/rflat)**2)**(1/2.) #1/((1+r**2)**p/2)

def plum4(arr,a0,x0,rflat):
	return (a0*rflat)/(1.+((arr-x0)/rflat)**2)**(3/2.)	

def calculateAveragingLength(f, d, ra_dec_inc, average_pc, ii_avg):
        '''
        function to calculate how many skeleton points are in
	one averaging length
	f               : numpy array of filament coordinates
        d               : distance to the object in parsec
        ra_dec_inc      : pixel increment in arcsec
        average_pc      : desired averaging length in parcsec
        '''
        length  = 0 
        while length <= average_pc:
		#while ii_avg < len(f)-2:
		
		x_1 = f[ii_avg,0]
		y_1 = f[ii_avg,1]
		x_2 = f[ii_avg+1,0]
		y_2 = f[ii_avg+1,1]
		delta_x = (x_2 - x_1)**2
		delta_y = (y_2 - y_1)**2
		delta = math.sqrt(delta_x + delta_y)
		theta = (delta*ra_dec_inc)/206265               # covert arcsec to radians
		R = theta*d
		length += R
		if ii_avg < len(f)-2:    
			ii_avg += 1
		else:
			break
			print 'done'
        return ii_avg 


def calculateWidth(filamentNo,filamentData, filamentLength, params, plotIndividualProfiles, printResults, plotSaveLoc):
        
	# reading user defined parameters 	
	npix    	= params["npix"] 
	avg_len 	= params["avg_len"] 
	avg_sep		= params["avg_sep"]
	
	niter   	= params["niter"]
	lam 		= params["lam"]
	pix 		= params["pixel_size"]
	dist 		= params["distance"] 
	fits_file	= params["fits_file"]
	add_fits_file	= params["peak_intensity"]
	dv 		= params["dv"]
	smooth		= params["smooth"]
	noise_level 	= params["noise_level"]
	int_thresh 	= params["int_thresh"]
	intensity	= fits.getdata(fits_file)
	peak_intensity	= fits.getdata(add_fits_file)
	beam_size	= params["beam_size"]
	
	
	resultsList	= []
	range_perp	= (np.arange(-npix,npix)*pix*dist/206265)
  	
 	n = 1

	# find how many skeleton points to iterate over for 1 averaging length
	jj_avg = 0
	nn_avg = 0
	nnList = [] 

	
	#check if the filament length is divisible by 3 beamsizes
	leftover	= filamentLength%(avg_len)

	
	# number of skeleton points that are divisible by 3 beamsizes
	# the lefover chunk is added then with the calculateAveraging funciton above. 
	numMainChunk = calculateAveragingLength(filamentData, dist, pix, (filamentLength-leftover), 0)	
	

	# here we create a list of skeleton points 
	# these points indicate up to which skeleton point the profiles should be integrated over
	# for example if a list is [15,26]
	# the profiles that are extracted perpendicular to the skeleton points from 0 to 15 would be 1 chunk
	# from 16 to 26 would be another chunk. Within each chunk, the profiles are averaged. 
	while jj_avg < numMainChunk-2:
        	nn_avg = calculateAveragingLength(filamentData, d=dist, ra_dec_inc=pix, average_pc=avg_len, ii_avg=nn_avg) 
		nnList.append((nn_avg))
        	jj_avg = nn_avg

       	averagingArray = np.array(nnList)

	
	# we divide the filament into smaller chunks so we can average 
	# the intensity profiles over these smaller chunks
	# the chunks are set to be 3 beam_size pieces and can be set with avg_len
	# in case the data is not divisble by avg_len there will be a left over chunk
	# these left over profiles will be averaged together		
	 
	print ' '
	print '######################################################################'
	print 'Calculating filament width using FilChaP'
	print 'Filament number		= ', filamentNo 	
	


	line_perpList			= []
	line_perpList2			= []	
	average_profile			= np.zeros((npix))
	average_profile2		= np.zeros((npix))	
	average_profileFull		= np.zeros((npix*2))
	stacked_profile			= np.zeros((len(filamentData),npix*2))
	stacked_profile_baseSubt 	= np.zeros((len(filamentData),npix*2))
	count_avg 			= 0

	# for each chunk we extract the intensity profiles
	for idx in averagingArray:
  		count_n		= 0	
		count_avg	+= 1
		
		
		print 'Number of chunks	= ', len(averagingArray)
		print 'Processed chunk		= ', count_avg, 'of', len(averagingArray) 	
		print '######################################################################'
		

		# this is a plotting bit	
		fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))	
		#'/home/suri/projects/carmaorion/analysis/disperse/c18o/0.2kms_resolution/han1_mask_imfit_c18o_pix_2_Tmb_noEdges_north_v4.7-13.1.peak.fits')
		fig.subplots_adjust(left=0.11,wspace=0.2, bottom=0.17)	
		ax1.imshow(peak_intensity, origin=0,interpolation="gaussian", vmax = 14, vmin = 0, extent=[0,peak_intensity.shape[1],0, peak_intensity.shape[0]])
	

		while n <= idx:
		
			count_n += 1	
			x = filamentData[:,0]
			y = filamentData[:,1]
			z = filamentData[:,2]

			ax1.plot(x,y,color='orange', marker='.', markersize=7., linestyle='None')
			v = [np.min(x)-120,np.max(x)+120,np.min(y)-120,np.max(y)+120]
			ax1.axis(v)
		
			x0 = int(x[n])
			y0 = int(y[n])
			z0 = int(z[n])
			r0 = np.array([x0,y0], dtype=float) 		# point on the filament
		
			#to save the data
			x0save = x0
			y0save = y0
			z0save = z0


			profileLSum = np.zeros(npix)      
			profileRSum = np.zeros(npix)
	
			
			###################################################################################
			# this is where we calculate the slices perpendicular to the spine of the filament
			# below loops are  part is implemented from Duarte-Cabral & Dobbs 2016.
			# we calculate distances in 2 separate for loops below
			# ################################################################################ 
			
			# find the tangent curve
			a = x[n+1] - x[n-1]
			b = y[n+1] - y[n-1]
			normal=np.array([a,b], dtype=float)
			#print a, b	
			#pl.plot(normal)
			#equation of normal plane: ax+by=const

			const = np.sum(normal*r0)
					
			# defining lists an array to be used below
			distance = np.zeros_like(intensity[0])
			distance2 = np.zeros_like(intensity[0])

			line_perp = np.zeros_like(intensity[0])
			line_perp2 = np.zeros_like(intensity[0])		
			
		
			# Loop 1: if the slope is negative
			if -float(b)/a > 0:
				npix            = params["npix"]
				print distance.shape
				if y0+npix < distance.shape[0] and x0+npix < distance.shape[1]:
					npix = npix
				#	print 'npix unmodified: ', npix
				elif y0+npix < distance.shape[0] and x0+npix >= distance.shape[1]:
					npix = distance.shape[1]-1- x0
				#	print 'npix modified: ', npix
				elif y0+npix >= distance.shape[0] and x0+npix < distance.shape[1]:
					npix = distance.shape[0]-1- y0
				#	print 'npix modified: ', npix
				elif y0+npix >= distance.shape[0] and x0+npix >= distance.shape[1]:
					npix_x = distance.shape[1]-1- x0
					npix_y = distance.shape[0]-1- y0
					if npix_x < npix_y:
						npix = npix_x
					else: 
						npix = npix_y
				
				for ii in range(y0-npix,y0+1):
					for jj in range(x0-npix,x0+1):
			
						distance[ii,jj]=((jj-x0)**2.+(ii-y0)**2.)**0.5 #distance between point (i,j) and filament
						if (distance[ii,jj] <  npix-1):
							dist_normal=(np.fabs(a*jj+b*ii-const))/(a**2+b**2)**0.5 #distance between point (i,j) and the normal 
							# take the point if it is in the vicinity of the normal (distance < 2 pix)
							if (dist_normal < 1):
								line_perp[ii,jj] = distance[ii,jj] #storing the nearby points
								line_perpList.extend((ii,jj,distance[ii,jj]))



				for ii in range(y0,y0+npix):
					for jj in range(x0,x0+npix):
	     
						distance2[ii,jj]=((jj-x0)**2.+(ii-y0)**2.)**0.5 
						if (distance2[ii,jj] <  npix-1):
							dist_normal2=(np.fabs(a*jj+b*ii-const))/(np.sum(normal*normal))**0.5 
							
							if (dist_normal2 < 1): 
								line_perp2[ii,jj] = distance2[ii,jj] 
								line_perpList2.extend((ii,jj,distance2[ii,jj]))			
							
			
			# Loop 2_ if the slope is positive
			elif -float(b)/a < 0:
				npix            = params["npix"]
				
				
				if y0+npix < distance.shape[0] and x0+npix < distance.shape[1]:
					npix = npix
					#print 'npix unmodified: ', npix
				elif y0+npix < distance.shape[0] and x0+npix >= distance.shape[1]:
					npix = distance.shape[1]-1- x0
					#print 'npix modified: ', npix
				elif y0+npix >= distance.shape[0] and x0+npix < distance.shape[1]:
					npix = distance.shape[0]-1- y0
					#print 'npix modified: ', npix
				elif y0+npix >= distance.shape[0] and x0+npix >= distance.shape[1]:
					npix_x = distance.shape[1]-1- x0
					npix_y = distance.shape[0]-1- y0
					if npix_x < npix_y:
						npix = npix_x
					else: 
						npix = npix_y
					#print 'npix modified: ', npix	
			
	
				for ii in range(y0,y0+npix):
					for jj in range(x0-npix,x0+1):
						distance[ii,jj]=((jj-x0)**2.+(ii-y0)**2.)**0.5
						if (distance[ii,jj] <  npix-1):
							dist_normal=(np.fabs(a*jj+b*ii-const))/(np.sum(normal*normal))**0.5
							if (dist_normal < 1):
								line_perp[ii,jj] = distance[ii,jj]
								line_perpList.extend((ii,jj,distance[ii,jj]))


				for ii in range(y0-npix,y0+1):
					for jj in range(x0, x0+npix):
						distance2[ii,jj]=((jj-x0)**2.+(ii-y0)**2.)**0.5
						if (distance2[ii,jj] <  npix-1):
							dist_normal2=(np.fabs(a*jj+b*ii-const))/(np.sum(normal*normal))**0.5
							if (dist_normal2 < 1):
								line_perp2[ii,jj] = distance2[ii,jj]
								line_perpList2.extend((ii,jj,distance2[ii,jj]))

			
			####################################################
			# now that we have the perpendicular slices ########
			# we can get the intensities along these slices ####
			####################################################
			npix            = params["npix"]

			perpendicularLine = np.array(line_perpList).reshape(-1,3)
			perpendicularLine2 = np.array(line_perpList2).reshape(-1,3)
				
			ax1.plot(perpendicularLine[:,1],perpendicularLine[:,0],'rx', markersize=4)
			ax1.plot(perpendicularLine2[:,1],perpendicularLine2[:,0],'rx', markersize=4)
			#print perpendicularLine
			#pl.show()
			
			for dd in range(0,npix):

				if (dd == 0):	
					# this is where the skeleton point is x0,y0,z0
					# sum the intensities of the velocity channel before and after
					profileLSum[dd] = np.sum([intensity[z0-1-1,y0-1,x0-1], intensity[z0-1,y0-1,x0-1], intensity[z0-1+1,y0-1,x0-1]]) 

				if (dd > 0):
					# this is where we have to get the list of the perpendicular points
					# it could be that close to the caculated perpendicular line
					# there are several points that have to same distance to the line
					# we take the mean intensity and some over 3 channels
					index_d = np.where((line_perp>dd-1) * (line_perp<=dd))	
		
					profileLSum[dd] = np.sum( [np.mean(intensity[int(z0)-1-1,index_d[0]-1,index_d[1]-1]), \
							np.mean(intensity[int(z0)-1,index_d[0]-1,index_d[1]-1]), \
							np.mean(intensity[int(z0)-1+1,index_d[0]-1,index_d[1]-1])] ) 
	
				
				# it could also be that what the perpendicular got was NaNs
				# in that case, ignore them
				# if not, average them 
				# the average profile is what we will use for the fitting 				
				if np.isnan(profileLSum[dd]) != True:
					average_profile[dd] += profileLSum[dd]
	

			for ddd in range(0,npix):        
				if (ddd == 0):
					# this is where the skeleton point is x0,y0,z0
					# sum the intensities of the velocity channel before and after	                       
					profileRSum[ddd] = np.sum([intensity[z0-1-1,y0-1,x0-1], intensity[z0-1,y0-1,x0-1], intensity[z0-1+1,y0-1,x0-1]]) 
		
				if (ddd > 0):
					# this is where we have to get the list of the perpendicular points
					# it could be that close to the caculated perpendicular line
					# there are several points that have to same distance to the line
					# we take the mean intensity and some over 3 channels
					index_d2 = np.where((line_perp2>ddd-1) * (line_perp2<=ddd))		
	
					profileRSum[ddd] = np.sum( [np.mean(intensity[int(z0)-1-1,index_d2[0]-1,index_d2[1]-1]), \
							np.mean(intensity[int(z0)-1,index_d2[0]-1,index_d2[1]-1]), \
							np.mean(intensity[int(z0)-1+1,index_d2[0]-1,index_d2[1]-1])] ) 
			
				
				if np.isnan(profileRSum[ddd]) != True:
					 average_profile2[ddd] += profileRSum[ddd]
	
			##############################################################
			# stack both sides of the intensity profiles #################
			##############################################################
		
			stacked_profile[n] = np.hstack((profileLSum[::-1],profileRSum))
			stacked_profile[n] = stacked_profile[n]*dv

		
			# subtract baselines from each of these profiles
		
			z = baseline_als(stacked_profile[n], lam, 0.01,niter)
			stacked_profile_baseSubt[n] = stacked_profile[n] -z	
			
			
			if plotIndividualProfiles == True:
				ax2.step(range_perp,stacked_profile_baseSubt[n],ls='-',color='k', lw=0.5, alpha=0.5)

			line_perpList, line_perpList2		= [], []	

			len_sep = 0
			print n, idx
			kk = n
	

			while len_sep <= avg_sep:
				sep_x	= (x[kk+1] - x[kk])**2
				sep_y	= (y[kk+1] - y[kk])**2
				del_xy 	= math.sqrt(sep_x + sep_y)
				dist_xy	= del_xy*pix*dist/206265
				len_sep += dist_xy
				print kk, len_sep,'//', avg_sep
				kk += 1
				n = kk
				
					
				
			print 'new n:', n
		#####################################################################################
		# exiting the first loop that allowed us to average a number of intensity profiles ##
		# this number is taken as three times the beam_size of the CARMA-NRO data ############
		# avg_length:12 , can be changed according to the used dataset. #####################
		# below, we fit this averaged profile to calculate the width ########################
		##################################################################################### 

		# this is the average radial intensity profile we need for the width calculation
		# so stack together both sides of the profile (- and +)
		# and multiply with the velocity channel width because it is an integrated intensity profile
		#print count_n	
		average_profile		= average_profile/count_n
		average_profile2	= average_profile2/count_n
		average_profileFull	= np.hstack((average_profile[::-1],average_profile2))
		average_profileFull	= average_profileFull*dv
		average_profile, average_profile2	= np.zeros((npix)), np.zeros((npix)) 	
			
		# subtract baseline from the averaged profile
		# and also smooth it 
		# the smoothed profile will be used to find dips and peaks

		z2 = baseline_als(average_profileFull, lam, 0.01,niter)
		y_base_subt = average_profileFull -z2
		y_base_subt_smooth = gaussian_filter(y_base_subt, sigma=smooth) #3 beam=12
			
				
		###################################################################################
		####### calculating minima ######################################################## 
		###################################################################################	
		# we calculate minima by looking at the minus side of the peak 
		# and to the plus side: minima left and right
		# in order to make sure the minima are global,
		# we put an integrated intensity threshold (at the moment 5*sigma) 
		# only minima that have values below this threshold will be taken into account
			
		minimaLeft = argrelextrema(y_base_subt_smooth[0:npix], np.less, order=6)
		minimaRight = argrelextrema(y_base_subt_smooth[npix:npix*2], np.less, order=6) 
		
		
		# following loops are where we decide which minima to use for the fit boundaries
		# in case there are multiple minima, the one close to the peak is selected
		# if there is no minima found, the entire range is used (from 0 to 2*npix).
			
		# left of the peak FILTER THE MINIMA HERE

		if len(minimaLeft[0]) > 1:
			
			b1 = minimaLeft[0][-1]
			ax2.axvline(x=range_perp[b1], ymin=0,ls='--',color='black', alpha=0.5)

		elif len(minimaLeft[0]) == 1:
			
			ax2.axvline(x=range_perp[minimaLeft[0][0]], ymin=0,ls='--',color='black', alpha=0.5)
			b1 = minimaLeft[0][0]
		else:

			b1 = 0
			ax2.axvline(x=range_perp[b1], ymin=0,ls='--',color='black', alpha=0.5)

			
		#right of the peak
		if len(minimaRight[0]) > 1:

			b2 = minimaRight[0][0]+npix
			ax2.axvline(x=range_perp[b2], ymin=0,ls='--',color='black', alpha=0.5)
		
		
		elif len(minimaRight[0]) == 1:
			ax2.axvline(x=range_perp[minimaRight[0][0]+npix], ymin=0,ls='--',color='black', alpha=0.5)
			b2 = minimaRight[0][0]+npix

		else:	
			b2 = 2*npix
			ax2.axvline(x=range_perp[b2-1], ymin=0,ls='--',color='black', alpha=0.5)

	
		# plot the averaged profile
		ax2.step(range_perp,y_base_subt,'k-', lw=2.0, alpha=1.0)
	 		
		# uncomment if you want to plot the smoothed average profile
		#pl.step(range_perp,y_base_subt_smooth,'g',lw=1.0, alpha=0.4)
		
		###################################################################################
		# here we calculate the number of peaks ###########################################
		# within our boundaries 		###########################################
		# this will help compare the number of peaks & shoulders to the width #############
		###################################################################################
		
		#Adopted from Seamus' peak finding. 
	
		print 'Finding Peaks'	
		ydata_og 	= y_base_subt[b1:b2]
		ydata		= y_base_subt_smooth[b1:b2]
		r		= range_perp[b1:b2]
		ny 		= len(ydata)
		minr 		= np.min(r)
		maxr 		= np.max(r)
		dr		= (maxr-minr)/ny

		noise			= noise_level/sqrt(3)/sqrt(count_n)
	
		limit		= 5*noise		# this is to check peak's significance	
		
		#derivatives
		dy		= np.zeros_like(ydata)
		for ii in range(0,ny-1):
			dy[ii] = (ydata[ii+1]-ydata[ii])/dr
		
		ddy 		= np.zeros_like(ydata)
		for ii in range(0,ny-2):
			ddy[ii] = (dy[ii+1]-dy[ii])/dr
		
		# work out the number of peaks and shoulders
		switch 		= np.zeros_like(ydata)
		decrease	= 0
		shoulder	= np.zeros_like(ydata)
		
		for ii in range(2,ny-2):
			# find a shoulder
			if(ddy[ii+1] > ddy[ii] and ddy[ii+2] > ddy[ii] and ddy[ii-1] > ddy[ii] and ddy[ii-2] > ddy[ii] and (ydata[ii]>limit or ydata[ii-1]>limit or ydata[ii+1]>limit)):
				
				shoulder[ii] = 1
			# find a peak
			if((dy[ii] < 0.0 and dy[ii-1]>0.0) and (ydata[ii]>limit or ydata[ii-1]>limit or ydata[ii+1]>limit)):

				switch[ii] = 1
		
		# check if there are any peaks detected	
		if( np.sum(switch) < 1 ):
			print "No peak was detected in this slice"
			print "Did I go wrong?"
			#return [[0,0,0],0]
		
		n_peaks = np.sum(switch)
		n_peaks = int(n_peaks)
		
		index = np.linspace(0,ny-1,ny)
		index = np.array(index,dtype=int)

		id_g = index[switch==1]
		cent_g = r[id_g]
		amp_g = ydata[id_g]
		
		is_shoulder = int(np.sum(shoulder)) - n_peaks
		newShoulderList = []

		if(is_shoulder > 0):
			
			# if there exists a shoulder we plot them with vertical dashed lines 
			shoulder_pos	= r[index[shoulder==1]]
			shoulder_amp	= ydata_og[index[shoulder==1]]
			print "Here are the shoulder positions", shoulder_pos
			
			for kk in range(len(shoulder_pos)):
				ax2.axvline(x=shoulder_pos[kk], ymin=0, ls='--', lw=1., color='g', alpha=0.5)
				newShoulderList.append(shoulder_pos[kk])
				'''
				if shoulder_amp[kk] > limit:
					if shoulder_pos[kk] > 0. and shoulder_pos[kk] > beam_size:
						pl.axvline(x=shoulder_pos[kk], ymin=0, ls='--', lw=1., color='g', alpha=0.5)	
						newShoulderList.append(shoulder_pos[kk])
					if shoulder_pos[kk] < 0. and shoulder_pos[kk] < -1*beam_size:
						pl.axvline(x=shoulder_pos[kk], ymin=0, ls='--', lw=1., color='g', alpha=0.5)	
						newShoulderList.append(shoulder_pos[kk])
				'''
		else:
			shoulder_pos    = []
			print 'I found no shoulders.'
		
		##################################################################################
		# finally calculating the width ################################################## 
		##################################################################################
		
		# initial guesses for the fits
		a		= np.amax(ydata_og)
		mu		= r[ np.argmax(ydata_og) ]
		pos_half	= np.argmin( np.abs( ydata_og-a/2 ) )
		sig 		= np.abs( mu - r[ pos_half] )
		p01 		= (a,mu,sig)
		p02 		= (a,mu,sig)


		#try:
		# 1st method: calculate moments	
		tot_2, tot_3, tot_4 = 0, 0, 0	
		try:
			for ii in range(len(r)): 
				tot_2 += ydata_og[ii]*(r[ii] - np.mean(r))**2
				tot_3 += ydata_og[ii]*(r[ii] - np.mean(r))**3
				tot_4 += ydata_og[ii]*(r[ii] - np.mean(r))**4
	       
			var = math.sqrt(tot_2/np.sum(ydata_og))
			mom3 = tot_3/np.sum(ydata_og)
			mom4 = tot_4/np.sum(ydata_og)
				
			FWHM_moments = var*2.35
			skewness = mom3/(var**3)
			kurtosis = mom4/(var**4) - 3

			#print 'moment:', FWHM_moments	
			# 2nd method: fit Gaussian and Plummer functions
	
			co_eff,var_matrix	= curve_fit(gaus,r, ydata_og,p0=p01,absolute_sigma=True)	
			#print co_eff
			co_eff3,var_matrix3	= curve_fit(plum2,r, ydata_og,p0=p02)
			#print co_eff3
			co_eff4,var_matrix4	= curve_fit(plum4,r, ydata_og,p0=p02)
			#print co_eff4
		except (RuntimeError,UnboundLocalError,ValueError,TypeError) as e:
			print 'fit cannot be calculated'
			pass
			
		#Calculate Chi-squared 
		num_freeParams		= 3  	        
		# gaussian fits
		chi_sq_gaus		= np.sum((ydata_og-gaus(r,*co_eff))**2) / noise**2
		red_chi_sq_gaus		= chi_sq_gaus / (len(ydata_og) - num_freeParams)			
		
		# plummer 2 fits
		chi_sq_plum2		= np.sum((ydata_og-plum2(r,*co_eff3))**2) / noise**2	
		red_chi_sq_plum2	= chi_sq_plum2 / (len(ydata_og) - num_freeParams)	
		
		# plummer 4 fits
		chi_sq_plum4		= np.sum((ydata_og-plum4(r,*co_eff4))**2) / noise**2
		red_chi_sq_plum4	= chi_sq_plum4 / (len(ydata_og) - num_freeParams)
		
		#fits
		fit			= gaus(range_perp,*co_eff)
		fit3			= plum2(range_perp,*co_eff3)
		fit4			= plum4(range_perp,*co_eff4)
	
		# fit standard deviation
		perr			= np.sqrt(np.diag(var_matrix))	
		perr3			= np.sqrt(np.diag(var_matrix3))
		perr4			= np.sqrt(np.diag(var_matrix4))	

		#pl.plot(r,y_base_subt)
		ax2.plot(range_perp,fit,ls='-.', color='#0000CD', lw=1.) 	 
		ax2.plot(range_perp,fit3,ls='-',color='#DAA520', lw=1.) 
		ax2.plot(range_perp,fit4,ls='--',color='red', lw=1.)	
		
		ax2.set_xlabel('Distance from the ridge [pc]')
		ax2.set_ylabel('Integrated Intensity [K.km/s]')
		#pl.axis('equal')	
		ax2.grid(True,alpha=0.3, ls='--')
		#ax2.gcf().subplots_adjust(bottom=0.15)		
		
		
		pl.savefig(str(plotSaveLoc) + 'filament'+ str(filamentNo)+'_slice' + str(idx) +'.png', dpi=300)
		#pl.show()
		pl.clf()
		rangePix = b2-b1
	 
		FWHM_plummer2 = 3.464*co_eff3[2]
		FWHM_plummer4 = 1.533*co_eff4[2]
		
		resultsList.extend((co_eff[0],perr[0],co_eff[1],perr[1],co_eff[2]*2.35,perr[2],co_eff3[0],perr3[0],co_eff3[1],perr3[1],FWHM_plummer2,perr3[2],co_eff4[0],perr4[0],co_eff4[1],perr4[1],FWHM_plummer4,perr4[2],FWHM_moments,skewness,kurtosis,chi_sq_gaus,red_chi_sq_gaus,chi_sq_plum2,red_chi_sq_plum2,chi_sq_plum4,red_chi_sq_plum4,rangePix,x0save,y0save,z0save,len(newShoulderList)))
	
		if printResults == True:
			print '###########################################################'
			print '############## Width Results ##############################'
			print ' '
			print 'Filament				=', filamentNo
			print 'FWHM (Second Moment)		=', FWHM_moments
			print 'FWHM (Gaussian Fit)		=', co_eff[2]*2.35
			print 'FWHM (Plummer 2)			=', FWHM_plummer2
			print 'FWHM (Plummer 4)			=', FWHM_plummer4 	
			print ' '
			print 'Skewness				=', skewness
			print 'Kurtosis				=', kurtosis
			print '###########################################################'	
	 

		
		resultsArray = np.array(resultsList).reshape(-1,32)
		#print resultsArray		
		
	return resultsArray


def readParameters(param_file):
	'''
	read parameters from the .param file
	this routine is taken from BTS (Clarke et al. 2018)
	https://github.com/SeamusClarke/BTS
	'''
	### The dictionaries for the type of variable and the variable itself
	type_of_var = {"npix"                        	: "int",
			"avg_len" 		     	: "float",
			"avg_sep"                       : "float",
			"fits_file"			: "str",
			"peak_intensity"		: "str",
			"distance"			: "float",
			"pixel_size"			: "float",
			"dv"				: "float",
			"noise_level"			: "float",
			"beam_size"			: "float",
			"int_thresh"			: "float",
			"niter"				: "int",
			"lam"				: "float",
			"smooth"			: "int"}
	param = {"npix"                 : 120,
		"avg_len" 		: 0.045,
		"avg_sep"               : 0.025,
		"fits_file"		: "c18o.fits",
		"peak_intensity"	: "/home/suri/projects/carmaorion/analysis/disperse/c18o/0.2kms_resolution/han1_mask_imfit_c18o_pix_2_Tmb_noEdges_north_v4.7-13.1.peak.fits",
		"distance"		: 388.0,
		"pixel_size"		: 2.0,
		"dv"			: 0.22,
		"noise_level"		: 0.47,
		"int_thresh"		: 0.4,
		"beam_size"		: 0.015,
		"niter"			: 100,
		"lam"			: 10e5,
		"smooth"		: 4}


        ### Open the file and read through, ignoring comments.

        with open(param_file) as f:

                for line in f:

                        if(line=="\n"):
                                continue

                        if(line[0]=="#"):
                                continue

                        words = line.split()	

                        try:
				var = type_of_var[words[0]]

                                if(var=="str"):
                                        param[words[0]]=words[2]
                                elif(var=="int"):
                                        param[words[0]]=np.int(words[2])
                                elif(var=="float"):
                                        param[words[0]]=np.float(words[2])
                                else:

                                        print "The variable is neither a string, float or integer. I don't know how to deal with this"

                        except KeyError:

                                print "There is no such parameter. Add it to the type_of_var and param dictionaries"


                f.close()

        ### Print the parameters to screen
	'''

        print " "
        print " "
        print " "
        print "######################################################################################"
        print "################ Parameters ##########################################################"
        print "######################################################################################"
        print " "
        print "############# Important two  #########################################################"
        print "Number of pixels in a perpendicular slice		= ", param["npix"]
	print "Number of points along the filament to be averaged	= ", param["avg_len"]
        print " "
	print "############# Data Specifics #########################################################"
	print "Fits file					= ", param["fits_file"]
        print "Distance to the cloud (pc)               	= ", param["distance"]
        print "Pixel size (arcsec)				= ", param["pixel_size"]	
	#print "Velocity resolution				=", param["dv"]
	print "Noise level of the averaged intensity profile	=", param["noise_level"]
	print "Intensity threshold for finding boundaries	=", param["int_thresh"]
	print " "
	print "############# For baseline subtraction  ##############################################"
	print "Number of iterations		= ", param["niter"]
        print "Lambda for ALS			= ", param["lam"]
	print "############# For finding minima #####################################################"
	print "Smooth over (pixels)		= ", param["smooth"]
	print " "
	'''
        return param
