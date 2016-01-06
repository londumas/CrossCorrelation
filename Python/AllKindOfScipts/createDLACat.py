# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

### My tools
import myTools


import astropy.io.fits as pyfits
import math
import numpy
import os
import decimal ## To set the precision of the double values
import cosmolopy.distance as cosmology
import matplotlib.pyplot as plt

import myTools
from myTools import Get_TProfile
from const_delta import *



'''
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/List/DLA_DR12_Pasquier.txt')

size = data[:,0].size
print size
col_ra              = pyfits.Column(name='RA',  format='D', array=numpy.zeros(size), unit='deg')
col_de              = pyfits.Column(name='DEC', format='D', array=numpy.zeros(size), unit='deg')
col_zz              = pyfits.Column(name='Z',   format='D', array=data[:,3] )
col_plate           = pyfits.Column(name='PLATE',     format='J', array=data[:,1].astype(int) )
col_mjd             = pyfits.Column(name='MJD',       format='J', array=data[:,0].astype(int) )
col_fiber           = pyfits.Column(name='FIBERID',   format='J', array=data[:,2].astype(int) )
col_nh1             = pyfits.Column(name='NHI',       format='D', array=data[:,4] )
	
tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz,col_plate,col_mjd,col_fiber,col_nh1])
tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DLA_all.fits', clobber=True)
'''




cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DLA_all.fits')[1].data
print cat.size

print cat['RA']

## Map
plt.plot(cat['RA'],  cat['DEC'], linestyle="", marker="o")
plt.xlabel(r'$R.A. (\degree)$')
plt.ylabel(r'$Dec. (\degree)$')
myTools.deal_with_plot(False,False,True)
plt.show()
## Distribution redshift
plt.hist(cat['Z'], bins=50)
plt.xlabel(r'$z$')
plt.ylabel(r'$\#$')
myTools.deal_with_plot(False,False,True)
plt.show()
## Distribution 
plt.hist(cat['PLATE'], bins=50)
plt.xlabel(r'$z$')
plt.ylabel(r'$\#$')
myTools.deal_with_plot(False,False,True)
plt.show()
## Distribution 
plt.hist(cat['NHI'], bins=50)
plt.xlabel(r'$z$')
plt.ylabel(r'$\#$')
myTools.deal_with_plot(False,False,True)
plt.show()






