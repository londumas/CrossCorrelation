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


### Cosmology
cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}




### File Jean-Marc
cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Raw/mocks-0_1000.fits')
print cat[0].header
print cat[0].data
print cat[1].header
lambdaTST = []
fluxTST   = []
zTST = []
nbPixel = []
distPixel = []
distPixel2 = []
deltaDist = []
for i in range(1,10): ##47
	print i

	try:
		el = cat[i]
	except Exception,error:
		continue

	#print el.header
	z = el.header['ZQSO']
	el = el.data
	#if (el['flux'][ el['flux']<0.05 ].size==0): continue
	plt.errorbar( el['lambda']/(1.+z), el['flux'],marker='o') #/(1.+z)
	lambdaTST = numpy.append(el['lambda']/(1.+z), lambdaTST) #/(1.+z)
	fluxTST   = numpy.append(el['flux'], fluxTST)
	zTST = numpy.append(z,zTST)
	nbPixel = numpy.append(el['lambda'].size,nbPixel)
	print i, el['lambda'].size, z
	myTools.deal_with_plot(False,False,True)
	plt.ylim([ -0.1, 1.1])
	
	distPixel  = numpy.append(distPixel, [ el['lambda'][i+1]-el['lambda'][i] for i in range(0,el['lambda'].size-1) ])
	distPixel2 = numpy.append(distPixel2, [ (el['lambda'][i+1]+el['lambda'][i])*0.5 for i in range(0,el['lambda'].size-1) ])
	dist = cosmology.comoving_distance(el['lambda']/1215.67-1., **cosmo)*h__
	deltaDist = numpy.append( deltaDist, [ dist[i+1]-dist[i] for i in range(0,el['lambda'].size-1) ] )
	
	
	'''
	print el['lambda'][0],dist[0]
	for i in range(0,el['lambda'].size-1):
		print el['lambda'][i+1], dist[i+1], distPixel[i], deltaDist[i]
	'''
	
plt.show()

plt.hist(deltaDist,bins=100)
plt.show()

lambdaTST = numpy.asarray(lambdaTST)
print numpy.min(lambdaTST), numpy.max(lambdaTST)
fluxTST   = numpy.asarray(fluxTST)
weight    = numpy.ones(fluxTST.size)
print lambdaTST
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaTST,fluxTST, 150,weight)	
plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
myTools.deal_with_plot(False, False, False)
plt.show()

plt.hist(zTST)
plt.show()
plt.hist(nbPixel)
plt.show()
plt.hist(distPixel)
plt.show()


plt.errorbar(distPixel2,deltaDist, fmt="o")
myTools.deal_with_plot(False,False,False)
plt.xlabel(r'$(d1+d2)/2.$')
plt.ylabel(r'$d1-d2$')
plt.show()
print numpy.amin(deltaDist)
print numpy.amax(deltaDist)

cat = pyfits.open('/home/usr201/mnt/hdumasde/spectra-780-0.fits')
print cat[0].header
print cat[1].header
cat = cat[1]
print cat.header['ZQSO']
z = cat.header['ZQSO']
cat = cat.data
print cat.size

cat['loglam'] = numpy.power(10.,cat['loglam'])
cat['loglam'] /= (1.+z)

plt.plot( cat['loglam'], cat['flux'] ,label='flux')
plt.plot( cat['loglam'], cat['mock_F']  ,label='mock_F')
plt.plot( cat['loglam'], cat['mock_contpca']  ,label='mock_contpca')
plt.plot( cat['loglam'], cat['ivar']  ,label='ivar')
plt.plot( cat['loglam'], cat['model']  ,label='model')
plt.plot( cat['loglam'], cat['mock_miscalib']  ,label='mock_miscalib')
#plt.plot( cat['loglam'], cat['and_mask']  ,label='and_mask')
plt.plot( cat['loglam'], cat['mock_ivar']  ,label='mock_ivar')
myTools.deal_with_plot(False,False,True)
plt.show()

'''
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
'''





