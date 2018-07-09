#!/usr/bin/env python
##################################################################
# this is an algorithm that runs the width calculation script
# for different filaments in paralell
#
##################################################################
from mpi4py import MPI
import numpy as np
import filchap
from astropy.io import fits
comm = MPI.COMM_WORLD



#set up MPI parameters
pos= comm.Get_rank()
num_procs = comm.Get_size()	
num_filaments = 483.		# number of filaments
fil_per_proc = int(num_filaments//num_procs)
print ' '
print '##############################'
print '##############################'
print 'fil per proc is %i' %(fil_per_proc)
print ' '

offset = pos*fil_per_proc
leftover = num_filaments%num_procs

if leftover != 0. and pos==num_procs-1:
    fil_per_proc += leftover

print ''
print '###############################################################'
print "FilChaP the magical filament algorithm"
print 'running in %i cores'%num_procs
print "I'm processor %i and i have %i filaments "%(pos, fil_per_proc)
print '################################################################'
print ''

params = filchap.readParameters("./userDefined.param")

for jj in range(fil_per_proc):
	try:

        	filamentNo = jj+offset 
        	filamentData      = np.loadtxt('/home/suri/scripts/python/skeletonProcessing/filaments_north/averaged/filament' + str(filamentNo)+'.txt')
       	
		if len(filamentData) >= 12:
            		try:

                		results = filchap.calculateWidth(filamentNo,filamentData,params,plotIndividualProfiles=True,printResults=True)
                		np.savetxt('/home/suri/development/filchap_1.0/c18o_north/width/c18owidth_north' + str(filamentNo)+'.txt',results)

            		except (ValueError,UnboundLocalError):
                		print 'There is something wrong with the fitting.'
                		pass

	except (IOError,IndexError):
        	print 'Passing the filament.'
       		pass

