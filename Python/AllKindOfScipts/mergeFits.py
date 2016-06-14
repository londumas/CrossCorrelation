# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >



import subprocess
import sys
import os

import numpy
import scipy
import astropy.io.fits as pyfits
#from iminuit import Minuit
import matplotlib.pyplot as plt
import re
#import cosmolopy.distance as cosmology
#import array
#import matplotlib.patches as mpatches
#import decimal
#import profile
import warnings
warnings.filterwarnings("error")



import myTools
from const_delta import *

raw_path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_Margala/'
all_t = [raw_path+'DR12_primery/DR12_primery.fits',raw_path+'DR12_reObs/DR12_reObs.fits']


### Get the FitsFile
all_t_file  = []
for el in all_t:
	all_t_file.append( pyfits.open(el,  memmap=True) )
	print el

size_1 = all_t_file[0][1].data.size
print size_1
size_2 = all_t_file[1][1].data.size
print size_2

### Get arrays of size of FitsFile
all_t_nrows = []
for el in all_t_file:
	all_t_nrows.append( el[1].data.shape[0] )
print all_t_nrows
nrowsTot = numpy.sum(all_t_nrows)
all_t_nrows = numpy.array(all_t_nrows)

print 
print '  ', all_t_nrows.size
print '  ', nrowsTot

### Set the Fits_File which will contain all
hdu = pyfits.BinTableHDU.from_columns(all_t_file[0][1].columns, nrows=nrowsTot)

### Set the values of each rows
first = 0
last  = all_t_nrows[0]
for i in range(1, all_t_nrows.size ):
	first += all_t_nrows[i-1]
	last  += all_t_nrows[i]
	print first, last
	for colname in all_t_file[0][1].columns.names:
		hdu.data[colname][first:last] = all_t_file[i][1].data[colname]

cat_tbhdu = hdu.data
cat_tbhdu['BIT'][size_1:] = 10000


## Map
#plt.plot(cat_tbhdu["RA"], cat_tbhdu["DEC"], linestyle="", marker="o")
#plt.show()
print cat_tbhdu['BIT']

tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
tbhdu.update()
tbhdu.writeto(raw_path+'DR12_coAdd/DR12_primery_reObs.fits', clobber=True)







