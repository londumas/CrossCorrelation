# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

### My tools

import astropy.io.fits as pyfits
import math
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import subprocess
import sys
import os


import myTools
from const_delta import *

#path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/vXXXX/Box_000/Simu_000/Data/'
#name__ = 'delta_noNoise_noCont.fits'
name__ = 'delta.fits'
nbPixel = 645
sizeMax = 291271


def haveAlookForest():

	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_001/Data/delta.fits', memmap=True)[1].data[0:1000]
	#cat2 = pyfits.open('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits', memmap=True)[1].data[0:1000]

	"""
	print cat.size
	cat = cat[ (cat['Z'] != 0.) ]
	print cat.size

	print numpy.amin(cat['X']), numpy.amax(cat['X'])
	print numpy.amin(cat['Y']), numpy.amax(cat['Y'])
	print numpy.amin(cat['Z']), numpy.amax(cat['Z'])

	
	from myTools import Get_TProfile
	flux      = cat["NORM_FLUX"][ cat['NORM_FLUX_IVAR'] >0.]/alphaStart__
	template  = cat["TEMPLATE"][ cat['NORM_FLUX_IVAR'] >0.]
	ivar      = cat["NORM_FLUX_IVAR"][ cat['NORM_FLUX_IVAR'] >0.]
	delta     = cat["DELTA"][ cat['NORM_FLUX_IVAR'] >0.]
	lambdaRF  = numpy.dot(numpy.diag(1./(1.+cat['Z'])),cat['LAMBDA_OBS'])[ cat['NORM_FLUX_IVAR'] >0.]
	lambdaObs = cat["LAMBDA_OBS"][ cat['NORM_FLUX_IVAR'] >0.]
	weight    = numpy.ones(flux.size)
	###
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,flux, 160,weight)	
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()
	###
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,flux/template, 160,weight)	
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()
	###
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,delta, 160,weight)	
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()
	###
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,template, 160,weight)	
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()
	###
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,ivar, 160,weight)	
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()
	###
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaObs,flux, 3000,weight)	
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()

	ar_cut          = (cat['NORM_FLUX_IVAR']>0.)
	ar_lambdaRF     = numpy.dot(numpy.diag(1./(1.+cat['Z'])),cat['LAMBDA_OBS'])[ cat['NORM_FLUX_IVAR'] >0.]
	ar_lambdaObs    = cat['LAMBDA_OBS'][ ar_cut ]
	ar_flux         = cat['NORM_FLUX'][ ar_cut ]
	ar_flux_ivar    = cat['NORM_FLUX_IVAR'][ ar_cut ]
	ar_zi           = cat['LAMBDA_OBS'][ ar_cut ]/(lambdaRFLine__) -1.
	ar_weight       = cat['DELTA_WEIGHT'][ ar_cut ]

	print numpy.amin(ar_zi), numpy.amax(ar_zi)
	
	"""
	for i in range(0,cat.size):
		el = cat[i]
		lamndaRF = el['LAMBDA_OBS']/(1.+el['Z'])
		if (lamndaRF[0]>1041.): continue
		cut = el['NORM_FLUX_IVAR']>0.
		lambdaAAA = el['LAMBDA_OBS'][cut]

		template = (el['ALPHA']+el['BETA']*(lamndaRF[cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]
		print el['ALPHA'], el['BETA'], el['MEAN_FOREST_LAMBDA_RF'], el['Z']

		plt.errorbar(lambdaAAA,el['NORM_FLUX'][cut],yerr=numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='$Data$')

		plt.plot(lambdaAAA,numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='flux err')
		plt.plot(lambdaAAA,el['DELTA'][cut],label='delta')
		plt.plot(lambdaAAA,el['DELTA_WEIGHT'][cut],label='delta weight')		
		plt.plot(lambdaAAA,template,label=r'$QSO \, template$',color='red')
		
		myTools.deal_with_plot(False,False,True)
		plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
		plt.ylabel(r'$flux$', fontsize=40)
		plt.show()
	

	
	### delta vs. lambda_RF
	xxx, yyy, eyyy, nyyy = myTools.Get_TProfile(ar_lambdaRF,ar_flux, lambdaRFTemplateBinEdges__, ar_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$flux$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### delta vs. lambda_Obs
	xxx, yyy, eyyy, nyyy = myTools.Get_TProfile(ar_lambdaObs,ar_flux, lambdaObsBinEdges__, ar_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, fmt="o")
	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$flux$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.show()
	
	### Distribution redshift
	plt.hist(ar_zi, bins=50)
	plt.xlabel(r'$z_{pixel}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution norm_flux_ivar
	plt.hist(numpy.log(ar_flux_ivar), bins=50, histtype='step')
	plt.hist(numpy.log(cat2['NORM_FLUX_IVAR'][ cat2['NORM_FLUX_IVAR']>0. ]), bins=50, histtype='step')
	plt.xlabel(r'$norm \, flux \, ivar$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,True,False)
	plt.show()
	### Map
	plt.plot(cat['RA'],cat['DEC'], linestyle="", marker="o")
	plt.xlabel(r'$X [Mpc]$')
	plt.ylabel(r'$Y [Mpc]$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	
	### Distribution redshift
	plt.hist(cat['Z_VI'], bins=50)
	plt.xlabel(r'$z_{forest}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	
	return
def haveAlookQSO():	

	#	cat = pyfits.open(path__+'QSO.fits', memmap=True)[1].data
	cat   = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_test_new_file_composition/Box_000/Simu_000/Data/QSO_withRSD.fits', memmap=True)[1].data
	#cat2  = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Data/QSO_noRSD.fits', memmap=True)[1].data

	print cat.size
	
	print numpy.amin(cat['X']), numpy.amax(cat['X'])
	print numpy.amin(cat['Y']), numpy.amax(cat['Y'])
	print
	print (numpy.amax(cat['X'])+numpy.amin(cat['X']))/(4.5*0.7)
	print (numpy.amax(cat['Y'])+numpy.amin(cat['Y']))/(4.5*0.7)
	print
	print numpy.amin(cat['Z']), numpy.amax(cat['Z'])
	print 1280*1280
	print cat.size
	print 1.*cat.size/(1280.*1280.)

	print
	
	### Map
	plt.plot(cat['X'],cat['Y'], linestyle="", marker="o")
	plt.xlabel(r'$X \, [h^{-1}.Mpc]$')
	plt.ylabel(r'$Y \, [h^{-1}.Mpc]$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution X
	plt.hist(cat['X'], bins=50)
	plt.xlabel(r'$X$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution Y
	plt.hist(cat['Y'], bins=50)
	plt.xlabel(r'$Y$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	
	### Distribution redshift
	plt.hist(cat['Z'], bins=50)
	plt.xlabel(r'$z_{QSO}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	'''
	### Distribution redshift
	plt.hist(cat2['Z'], bins=50)
	plt.xlabel(r'$z_{QSO}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution redshift
	plt.hist(cat['Z']-cat2['Z'], bins=50)
	plt.xlabel(r'$z_{QSO}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	'''

	return


#cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation/Box_000/Simu_000/Data/delta.fits', memmap=True)[1].data
#print cat.size

#createEmptyFitsFile()
#copyToNormalFits()
#allNeededFileds()
#saveCatalogueQSO()
#saveCatalogueQSOFromAscii()
#haveAlookQSO()
#removeUselessLines()
haveAlookForest()

'''
col_xx              = pyfits.Column(name='X',  format='D', array=numpy.zeros(291271), unit='Mpc/h')
col_yy              = pyfits.Column(name='Y',  format='D', array=numpy.zeros(291271), unit='Mpc/h')
col_zz              = pyfits.Column(name='Z',  format='D', array=numpy.zeros(291271), unit='0 (redshift)')
	
tbhdu = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz])
tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/v1524/QSO_noRSD.fits', clobber=True)

'''

'''
nb = 20
nbQSO = 238929
a = numpy.zeros(nb)
for i in range(0,2):
	for j in range(0,10):
		print i, ' ', j
		cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_withRSD.fits', memmap=True)[1].data
		a[i*10+j] = cat[ (cat['Z']!=0.) ].size

plt.plot(numpy.arange(nb),(a-nbQSO)*100./nbQSO)
#plt.plot(numpy.arange(10),numpy.ones(10)*,color='red')
plt.show()
'''












