#!/usr/bin/env python

import matplotlib.pyplot as pl
import numpy as np
import matplotlib
import pandas as pd

matplotlib.rcParams.update({'font.size': 14.})
matplotlib.rcParams.update({'font.family':'serif'})

# read hdf file as pandas dataframe
df = pd.read_hdf('./widthResults.h5')

#create figure
fig1 = pl.figure()
ax1 = fig1.add_subplot(111)

# arrays to plot
# x-axis, number of slices is the length of the dataframe
xArray = np.arange(len(df))

# plot gausWidth, plum2Width, plum4Width and momentWid
# values from the dataframe
ax1.plot(xArray,df.gausWidth, '-', alpha=0.8, color='#0000CD', lw=1., label='Gaussian')
ax1.plot(xArray,df.plum2Width, '-', alpha=0.8, color='#DAA520', lw=1., label='Plummer (p=2)')
ax1.plot(xArray,df.plum4Width, '-', alpha=0.8, color='r', lw=1., label='Plummer (p=4)')
ax1.plot(xArray,df.momentWid, '-', alpha=0.9, color='k', lw=1., label='Second Moment')


pl.ylabel('Width [pc]')

pl.xlabel('Slice Number')
pl.grid(True, alpha=0.5)

pl.legend(loc='best', fontsize='small')
pl.show()

