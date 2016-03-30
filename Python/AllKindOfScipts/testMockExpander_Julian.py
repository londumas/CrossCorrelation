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



path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Program/LyAMockExpander/usefulFiles/PCA/pca_suzuki.fits'
cat1 = pyfits.open(path, memmap=True)[1].data


lambdaRF = cat1['LAMBDA'][0]
flux     = cat1['FLUX'][0]
normFactor = numpy.mean( flux[ numpy.logical_and( (lambdaRF>1275.), (lambdaRF<1295.) ) ] )
flux /= normFactor



'''
### File Jean-Marc
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Raw/mocks-0_1000.fits'
key = ['loglam','flux','mock_F','mock_contpca','ivar','model','mock_miscalib','and_mask','mock_ivar']

cat = pyfits.open(path, memmap=True)


print cat[1].header

lam  = numpy.array([])
lamO = numpy.array([])
fl   = numpy.array([])
pca  = numpy.array([])

isNotFinished = True
i = 1
while isNotFinished:

	try:
		z  = cat[i].header['ZQSO']
		el = cat[i].data
	except Exception,error:
		isNotFinished = False
		continue

	el = el[ (el['LOGLAM']>0.) ]
	el = el[ (el['IVAR']>0.) ]
	if (el.size == 0):
		i += 1
		print i, ' none ' 
		continue

	### Normalize
	tmp_logZ = numpy.log10(1.+z)
	if (el["FLUX"][ numpy.logical_and( (el["LOGLAM"]>log10lambdaRFNormaMin__+tmp_logZ), (el["LOGLAM"]<log10lambdaRFNormaMax__+tmp_logZ) ) ].size<=0.):
		i += 1
		print i, ' none ' 
		continue

	normFactor = numpy.mean( el["FLUX"][ numpy.logical_and( (el["LOGLAM"]>log10lambdaRFNormaMin__+tmp_logZ), (el["LOGLAM"]<log10lambdaRFNormaMax__+tmp_logZ) ) ] )
	if (normFactor<=0.):
		i += 1
		print i, ' none ' 
		continue

	el['FLUX'] /= normFactor

	el['LOGLAM'] = numpy.power(10., el['LOGLAM'] )
	el['LOGLAM'] /= (1.+z)
	#plt.plot(el['LOGLAM'],el['flux'], marker='o')
	#plt.plot(el['LOGLAM'],el['mock_contpca']/normFactor, color='cyan',label='one spectrum')
	print i, z

	#isNotFinished = False
	i += 1

	el = el[ el['LOGLAM']*(1.+z)>3600. ]
	el = el[ numpy.logical_and( el['LOGLAM']>1040.,el['LOGLAM']<1200. ) ]
	if (el["FLUX"].size<=0.):
		i += 1
		print i, ' none ' 
		continue

	lam  = numpy.append( el['LOGLAM'] ,lam )
	lamO = numpy.append( el['LOGLAM']*(1.+z) ,lamO )
	fl   = numpy.append( el['flux'] ,fl )
	pca  = numpy.append( el['mock_contpca']/normFactor ,pca )

	#if (i>100): isNotFinished = False


wei = numpy.ones(lam.size)


plt.plot(lambdaRF, flux, label='PCA input')
#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Run/hDeltaVsLambdaRF_LYA_0_0.txt')
#plt.errorbar(data[:,0]+1037.5, data[:,1], markersize=8,linewidth=2, label='template')
#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Run/template_LYA_0_0.txt')
#plt.errorbar(data[:,0], data[:,1], markersize=8,linewidth=2, label='mean flux')

xxx, yyy, eyyy, nyyy = Get_TProfile(lam,pca, 160,wei)	
plt.errorbar(xxx, yyy, label='mean PCA')
xxx, yyy, eyyy, nyyy = Get_TProfile(lam,fl, 160,wei)
plt.errorbar(xxx, yyy, label='v1575')
'''



cat = pyfits.open('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits', memmap=True)[1].data[0:10000]
flux      = cat["NORM_FLUX"][ cat['NORM_FLUX_IVAR'] >0.]/alphaStart__
lambdaObs = cat["LAMBDA_OBS"][ cat['NORM_FLUX_IVAR'] >0.]
lambdaRF  = numpy.dot(numpy.diag(1./(1.+cat['Z'])),cat['LAMBDA_OBS'])[ cat['NORM_FLUX_IVAR'] >0.]
weight    = numpy.ones(flux.size)
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,flux, 160,weight)	
plt.errorbar(xxx, yyy, label='data')

cat2       = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Data/delta.fits', memmap=True)[1].data[0:10000]
flux2      = cat2["NORM_FLUX"][ cat2['NORM_FLUX_IVAR'] >0.]/alphaStart__
lambdaObs2 = cat2["LAMBDA_OBS"][ cat2['NORM_FLUX_IVAR'] >0.]
lambdaRF2  = numpy.dot(numpy.diag(1./(1.+cat2['Z_VI'])),cat2['LAMBDA_OBS'])[ cat2['NORM_FLUX_IVAR'] >0.]
weight2    = numpy.ones(flux2.size)
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF2,flux2, 160,weight2)	
plt.errorbar(xxx, yyy, label='v1547')


cat3       = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Data/delta.fits', memmap=True)[1].data[0:10000]
flux3      = cat3["NORM_FLUX"][ cat3['NORM_FLUX_IVAR'] >0.]/alphaStart__
lambdaObs3 = cat3["LAMBDA_OBS"][ cat3['NORM_FLUX_IVAR'] >0.]
lambdaRF3  = numpy.dot(numpy.diag(1./(1.+cat3['Z'])),cat3['LAMBDA_OBS'])[ cat3['NORM_FLUX_IVAR'] >0.]
weight3    = numpy.ones(flux3.size)
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF3,flux3, 160,weight3)	
plt.errorbar(xxx, yyy, label='v1575')

myTools.deal_with_plot(False,False,True)
plt.show()



xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaObs,flux, 600,weight)	
plt.errorbar(xxx, yyy, label='data')
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaObs2,flux2, 600,weight2)	
plt.errorbar(xxx, yyy, label='v1547')
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaObs3,flux3, 600,weight3)	
plt.errorbar(xxx, yyy, label='v1575')

myTools.deal_with_plot(False,False,True)
plt.show()












	

