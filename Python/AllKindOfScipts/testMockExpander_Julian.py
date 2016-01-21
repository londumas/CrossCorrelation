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
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation/Box_000/Simu_000/Data/mocks-0_10000.fits'
key = ['loglam','flux','mock_F','mock_contpca','ivar','model','mock_miscalib','and_mask','mock_ivar']

cat = pyfits.open(path, memmap=True)

print cat[0].header
print cat[0].data
print cat[1].header
print cat[1].data


nb = []
zA = []
difLOGLOBS = []
difLOBS = []
difLRF = []

for i in numpy.arange(1,100):
	z  = cat[i].header['ZQSO']
	el = cat[i].data

	el = el[ (el['LOGLAM']>0.) ]
	el = el[ (el['IVAR']>0.) ]
	el['LOGLAM'] = numpy.power(10., el['LOGLAM'] )
	el['LOGLAM'] /= (1.+z)

	'''
	difLOGLOBS += [ el['loglam'][i+1]-el['loglam'][i] for i in numpy.arange(el['loglam'].size-1) if el['loglam'][i]!=0. ]
	el['loglam'][ (el['loglam']!=0.) ] = numpy.power(10., el['loglam'][ (el['loglam']!=0.) ] )
	difLOBS    += [ el['loglam'][i+1]-el['loglam'][i] for i in numpy.arange(el['loglam'].size-1) if el['loglam'][i]!=0. ]
	el['loglam'][ (el['loglam']!=0.) ] /= 1.+z
	difLRF     += [ el['loglam'][i+1]-el['loglam'][i] for i in numpy.arange(el['loglam'].size-1) if el['loglam'][i]!=0. ]
	'''

	plt.plot(el['LOGLAM'],el['flux'], marker='o')
	#plt.plot(el['LOGLAM'],el['ivar'], marker='o')
	plt.plot(el['LOGLAM'],el['mock_contpca'])
	myTools.deal_with_plot(False,False,False)
	plt.show()

	'''
	nb += [ el['loglam'][ numpy.logical_and( (el['loglam']>=lambdaRFTemplateMin__) , (el['loglam']<lambdaRFTemplateMax__) ) ].size ]
	zA += [z]
	'''

nb = numpy.array(nb)
difLOGLOBS = numpy.array(difLOGLOBS)
difLOBS    = numpy.array(difLOBS)
difLRF     = numpy.array(difLRF)
print numpy.amax(nb)

myTools.deal_with_plot(False,False,False)
plt.show()

print difLOGLOBS
plt.hist(difLOGLOBS, bins=100)
plt.show()
plt.hist(difLOBS, bins=100)
plt.show()
plt.hist(difLRF, bins=100)
plt.show()

plt.hist(zA, bins=100)
plt.show()

















	

