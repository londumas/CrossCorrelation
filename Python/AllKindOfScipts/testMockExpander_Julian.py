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




### File Jean-Marc
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation/TESTS/00/mocks-0.fits'
key = ['loglam','flux','mock_F','mock_contpca','ivar','model','mock_miscalib','and_mask','mock_ivar']

cat = pyfits.open(path)

print cat[2]



print cat[0].header
print cat[0].data
print cat[1].header

for i in numpy.arange(1,47):
	z  = cat[i].header['ZQSO']
	el = cat[i].data

	el['loglam'][ (el['loglam']!=0.) ] = numpy.power(10., el['loglam'][ (el['loglam']!=0.) ] )
	#el['loglam'][ (el['loglam']!=0.) ] /= 1.+z

	print el[ numpy.logical_and( (el['loglam']>=lambdaRFMin__) , (el['loglam']<lambdaRFMax__) ) ].size

	
	for j in key[1:2]:
		plt.plot(el['loglam'],el[j],label=j)
		myTools.deal_with_plot(False,False,True)
plt.show()
	

