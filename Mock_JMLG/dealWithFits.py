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

path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/vXXXX/Box_000/Simu_000/Data/'
#name__ = 'delta_noNoise_noCont.fits'
name__ = 'delta.fits'
nbPixel = 645
sizeMax = 291271

def createEmptyFitsFile():

	tmp_nbBinForest     = str(nbPixel)+'D'

	plate                = pyfits.Column(name='PLATE',           format='J', array=numpy.zeros(sizeMax) )
	mjd                  = pyfits.Column(name='MJD',             format='J', array=numpy.zeros(sizeMax) )
	fiber                = pyfits.Column(name='FIBERID',         format='J', array=numpy.zeros(sizeMax) )
	
	ra                   = pyfits.Column(name='RA',              format='D', array=numpy.zeros(sizeMax))
	de                   = pyfits.Column(name='DEC',             format='D', array=numpy.zeros(sizeMax))
	zz                   = pyfits.Column(name='Z_VI',            format='D', array=numpy.zeros(sizeMax))
	nb                   = pyfits.Column(name='NB_PIXEL',        format='I', array=numpy.zeros(sizeMax))
	meanLambdaRF         = pyfits.Column(name='MEAN_FOREST_LAMBDA_RF', format='D', array=numpy.zeros(sizeMax), unit='angstrom')
	alpha1               = pyfits.Column(name='ALPHA_1',         format='D', array=numpy.ones(sizeMax) )
	beta1                = pyfits.Column(name='BETA_1',          format='D', array=numpy.zeros(sizeMax) )
	alpha2               = pyfits.Column(name='ALPHA_2',         format='D', array=numpy.ones(sizeMax) )
	beta2                = pyfits.Column(name='BETA_2',          format='D', array=numpy.zeros(sizeMax) )

	tmp_nbBinForest     = str(nbPixel)+'D'
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',       format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)), unit='angstrom' )
	lambdaRFForest      = pyfits.Column(name='LAMBDA_RF',        format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)),  unit='angstrom' )
	normFluxForest      = pyfits.Column(name='NORM_FLUX',        format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	normFluxIvarForest  = pyfits.Column(name='NORM_FLUX_IVAR',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	fluxDLA             = pyfits.Column(name='FLUX_DLA',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)))
	deltaForest         = pyfits.Column(name='DELTA',            format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	deltaIvarForest     = pyfits.Column(name='DELTA_IVAR',       format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	deltaWeight         = pyfits.Column(name='DELTA_WEIGHT',     format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)))
	template            = pyfits.Column(name='TEMPLATE',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)))	
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, ra, de, zz, nb, meanLambdaRF, alpha1, beta1, alpha2, beta2, lambdaForest, lambdaRFForest, normFluxForest, normFluxIvarForest, fluxDLA, deltaForest, deltaIvarForest, deltaWeight, template])
	cat_tbhdu = tbhdu.data

	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
	tbhdu.update()
	tbhdu.writeto(path__+name__, clobber=True)

	return
def removeUselessLines():
	'''
	'''

	print path__+name__
	file_cat = pyfits.open(path__+name__) #,mode='update')
	cat = file_cat[1].data[:885]
	print cat.size
	cat = cat[ (cat['Z_VI'] != 0.) ]
	print cat.size

	#pyfits.writeto(path__+name__, cat, clobber=True)

	return
def haveAlookForest():	

	print path__+name__
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Data/delta.fits', memmap=True)[1].data[:100]


	print cat.size
	cat = cat[ (cat['Z_VI'] != 0.) ]
	print cat.size

	print numpy.amin(cat['RA']), numpy.amax(cat['RA'])
	print numpy.amin(cat['DEC']), numpy.amax(cat['DEC'])
	print numpy.amin(cat['Z_VI']), numpy.amax(cat['Z_VI'])


	ar_cut          = (cat['NORM_FLUX_IVAR']>0.)
	ar_lambdaObs    = cat['LAMBDA_OBS'][ ar_cut ]
	ar_lambdaRF     = cat['LAMBDA_RF'][ ar_cut ]
	ar_flux         = cat['NORM_FLUX'][ ar_cut ]
	ar_zi           = cat['LAMBDA_OBS'][ ar_cut ]/(lambdaRFLine__) -1.
	ar_weight       = cat['DELTA_WEIGHT'][ ar_cut ]

	print numpy.amin(ar_zi), numpy.amax(ar_zi)


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
def saveCatalogueQSO():
	'''
		Save the list of quasar to a FITS file with one header but many lines
	'''

	### Constants defined by user
	path = path__+name__

	cat = pyfits.open(path, memmap=True)[1].data

	col_xx              = pyfits.Column(name='X',  format='D', array=cat['X'], unit='Mpc/h')
	col_yy              = pyfits.Column(name='Y',  format='D', array=cat['Y'], unit='Mpc/h')
	col_zz              = pyfits.Column(name='Z',  format='D', array=cat['Z'], unit='0 (redshift)')
	
	tbhdu = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz])
	tbhdu.writeto(path__+'QSO.fits', clobber=True)
	
	return
def haveAlookQSO():	

#	cat = pyfits.open(path__+'QSO.fits', memmap=True)[1].data
	cat   = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Data/QSO_withRSD.fits', memmap=True)[1].data
	cat2  = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Data/QSO_noRSD.fits', memmap=True)[1].data

	print cat.size
	
	print numpy.amin(cat['X']), numpy.amax(cat['X'])
	print numpy.amin(cat['Y']), numpy.amax(cat['Y'])
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
	return

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












