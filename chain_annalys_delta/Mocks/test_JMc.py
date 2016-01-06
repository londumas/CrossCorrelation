# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

### My tools

import sys
import astropy.io.fits as pyfits
import math
import numpy
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



import myTools
from const_delta import *

location__     = 'ICLUST'   ### "HOME" or "ICLUST"
pipeline__     = 'MOCK'     ### "DR12" or 'Guy' or 'Margala' or "MOCK" or 'eBOSS'
reObs__        = False      ### False or True
dataFitsType__ = 'spec'     ### What type of FITS file will we use: 'spPlate' or 'spec'

def copyToNormalFits():

	import warnings
	warnings.filterwarnings("error")

	### h=0.71, Omega_M = 0.27, Omega_L = 0.73.

	path = '/home/nfs/manip/mnt/bao/legoff/spectra-1280.fits'
	cat = pyfits.open(path, memmap=True)

	sizeMax             = len(cat)
	print '  size = ', len(cat)-1, "  sizeMax = ", sizeMax
	del cat

	nbPixel             = 645
	tmp_nbBinForest     = str(nbPixel)+'D'
	col_xx              = pyfits.Column(name='X',                      format='D',             array=numpy.zeros(sizeMax), unit='Mpc/h')
	col_yy              = pyfits.Column(name='Y',                      format='D',             array=numpy.zeros(sizeMax), unit='Mpc/h')
	col_zz              = pyfits.Column(name='Z',                      format='D',             array=numpy.zeros(sizeMax), unit='0 (redshift)' )
	col_normFlux        = pyfits.Column(name='NORM_FACTOR',            format='D',             array=numpy.zeros(sizeMax), unit='Flux')
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',             format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)), unit='Angstrom')
	normFluxForest      = pyfits.Column(name='FLUX',                   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	normFluxIvarForest  = pyfits.Column(name='FLUX_ERR',               format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	tbhdu = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz,col_normFlux, lambdaForest, normFluxForest, normFluxIvarForest])
	cat_tbhdu = tbhdu.data


	idx = 0
	for i in range(0,1000):
		if (idx >= sizeMax): break
		cat = pyfits.open(path, memmap=True)
		for j in range(0,1000):
			if (idx >= sizeMax): break

			el = cat[idx+1]

			### Get x, y, z
			cat_tbhdu[idx]['X'] = el.header['X']
			cat_tbhdu[idx]['Y'] = el.header['Y']
			cat_tbhdu[idx]['Z'] = el.header['ZQSO']
			zz = cat_tbhdu[idx]['Z']

			el = el.data
			### CCD and too many sky lines
			el = el[ (numpy.logical_and( (el["LAMBDA"]>=lambdaObsMin__) , (el["LAMBDA"]<lambdaObsMax__) )) ]

			### Get the normalisation factor
			try:
				normFactor = numpy.mean( el['FLUX'][ numpy.logical_and( (el["LAMBDA"]>=(1.+zz)*lambdaRFNormaMin__) , (el["LAMBDA"]<(1.+zz)*lambdaRFNormaMax__) ) ] )
			except Exception,error:
				tmp_command = "echo  \"" + "  Error in: 'normFactor = numpy.mean',  z = " + str(zz) + ", idx = " + str(idx) + "\""
				subprocess.call(tmp_command, shell=True)
				tmp_command = "echo  \"  " + str(error)  + "\n\""
				subprocess.call(tmp_command, shell=True)
				continue
			
			cat_tbhdu[idx]['NORM_FACTOR'] = normFactor

			### Only the forest
			el = el[ numpy.logical_and( (el["LAMBDA"]>=(1.+zz)*lambdaRFTemplateMin__) , (el["LAMBDA"]<(1.+zz)*lambdaRFTemplateMax__) ) ]
			nb = el.size

			cat_tbhdu[idx]['LAMBDA_OBS'][:nb] = el['LAMBDA']
			cat_tbhdu[idx]['FLUX'][      :nb] = el['FLUX']
			cat_tbhdu[idx]['FLUX_ERR'][  :nb] = el['DFLUX']
			idx += 1
			
		cat.close()
	
	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
	tbhdu.update()
	tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/mock_1280_0_0.fits', clobber=True)
	

def allNeededFileds():	

	#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/mock_0_0.fits'
	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/mock_1280_0_0.fits'
	cat = pyfits.open(path)[1].data
	sizeMax = cat.size
	nbPixel = cat['FLUX'][0].size
	print sizeMax, nbPixel

	ar_lambdaRF = numpy.array(cat['LAMBDA_OBS'])
	### Set values of columns
        step = 500
        iStart = 0
        iEnd   = step
        while (iStart<=sizeMax):
                print iStart, iEnd
		cat['FLUX'][iStart:iEnd]      = numpy.dot(numpy.diag(1./cat['NORM_FACTOR'][iStart:iEnd]),cat['FLUX'][iStart:iEnd])
		cat['FLUX_ERR'][iStart:iEnd]  = numpy.dot(numpy.diag(1./cat['NORM_FACTOR'][iStart:iEnd]),cat['FLUX_ERR'][iStart:iEnd])
                ar_lambdaRF[iStart:iEnd]      = numpy.dot(numpy.diag(1./(1.+cat['Z'][iStart:iEnd])),cat['LAMBDA_OBS'][iStart:iEnd])
		iStart = iEnd+1
                iEnd   += step
                if (iEnd>sizeMax):
                        iEnd = sizeMax

	ar_fluxIvar = 1./(cat['FLUX_ERR']**2.)
	ar_fluxIvar[ (numpy.isnan(ar_fluxIvar))  ] = 0.
	ar_fluxIvar[ (numpy.isinf(ar_fluxIvar))  ] = 0.

	### Remove forest with less than nbBinRFMin__ pixels in the forest
	cut_noForestPixel = numpy.logical_and( (ar_lambdaRF>=lambdaRFMin__) , (ar_lambdaRF<lambdaRFMax__) ).astype(int)
	len_forest = numpy.sum(cut_noForestPixel,axis=1)
	mean_lambdaRF_in_forest = numpy.sum(cut_noForestPixel*ar_lambdaRF,axis=1)/len_forest

	print len_forest

	### define fits file
	plate                = pyfits.Column(name='PLATE',           format='J', array=numpy.zeros(sizeMax) )
	mjd                  = pyfits.Column(name='MJD',             format='J', array=numpy.zeros(sizeMax) )
	fiber                = pyfits.Column(name='FIBERID',         format='J', array=numpy.zeros(sizeMax) )
	
	ra                   = pyfits.Column(name='RA',              format='D', array=cat['X'])
	de                   = pyfits.Column(name='DEC',             format='D', array=cat['Y'])
	zz                   = pyfits.Column(name='Z_VI',            format='D', array=cat['Z'] )
	nb                   = pyfits.Column(name='NB_PIXEL',        format='I', array=len_forest)
	meanLambdaRF         = pyfits.Column(name='MEAN_FOREST_LAMBDA_RF', format='D', array=mean_lambdaRF_in_forest, unit='angstrom')
	alpha1               = pyfits.Column(name='ALPHA_1',         format='D', array=numpy.ones(sizeMax) )
	beta1                = pyfits.Column(name='BETA_1',          format='D', array=numpy.zeros(sizeMax) )
	alpha2               = pyfits.Column(name='ALPHA_2',         format='D', array=numpy.ones(sizeMax) )
	beta2                = pyfits.Column(name='BETA_2',          format='D', array=numpy.zeros(sizeMax) )

	tmp_nbBinForest     = str(nbPixel)+'D'
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',       format=tmp_nbBinForest, array=cat['LAMBDA_OBS'], unit='angstrom' )
	lambdaRFForest      = pyfits.Column(name='LAMBDA_RF',        format=tmp_nbBinForest, array=ar_lambdaRF,  unit='angstrom' )
	normFluxForest      = pyfits.Column(name='NORM_FLUX',        format=tmp_nbBinForest, array=cat['FLUX'] )
	normFluxIvarForest  = pyfits.Column(name='NORM_FLUX_IVAR',   format=tmp_nbBinForest, array=ar_fluxIvar )
	fluxDLA             = pyfits.Column(name='FLUX_DLA',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)) )
	deltaForest         = pyfits.Column(name='DELTA',            format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)) )
	deltaIvarForest     = pyfits.Column(name='DELTA_IVAR',       format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)) )
	deltaWeight         = pyfits.Column(name='DELTA_WEIGHT',     format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)) )
	template            = pyfits.Column(name='TEMPLATE',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)) )	
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, ra, de, zz, nb, meanLambdaRF, alpha1, beta1, alpha2, beta2, lambdaForest, lambdaRFForest, normFluxForest, normFluxIvarForest, fluxDLA, deltaForest, deltaIvarForest, deltaWeight, template])
	cat_tbhdu = tbhdu.data

	print cat_tbhdu['NB_PIXEL']
	print cat_tbhdu.size
	cat_tbhdu = cat_tbhdu[ (cat_tbhdu['NB_PIXEL']>=nbBinRFMin__) ]
	print cat_tbhdu.size

	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
	tbhdu.update()

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/new_mock_1280_0_0.fits'
	tbhdu.writeto(path, clobber=True)

def haveAlookForest():	

	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/new_mock_1280_0_0.fits', memmap=True)[1].data
	print cat.size
	print numpy.amin(cat['RA']), numpy.amax(cat['RA'])
	print numpy.amin(cat['DEC']), numpy.amax(cat['DEC'])
	print numpy.amin(cat['Z_VI']), numpy.amax(cat['Z_VI'])

	ar_cut          = (cat['NORM_FLUX_IVAR']>0.)
	ar_lambdaObs    = cat['LAMBDA_OBS'][ ar_cut ]
	ar_lambdaRF     = cat['LAMBDA_RF'][ ar_cut ]
	ar_flux         = cat['NORM_FLUX'][ ar_cut ]
	ar_zi           = cat['LAMBDA_OBS'][ ar_cut ]/(lambdaRFLine__) -1.
	ar_weight       = cat['DELTA_WEIGHT'][ ar_cut ] #numpy.ones(shape=(cat.size,cat['NORM_FLUX'][0].size))[ ar_cut ]

	print numpy.amin(ar_zi), numpy.amax(ar_zi)

	### deltaLoglambda (mock_JMc)    = 0.00023025850929947467
	### deltaLoglambda (data)        = 0.00023002222870971423
	### deltaLoglambda (mock Julian) = 0.00023002222870971423,

	### delta vs. lambda_RF
	xxx, yyy, eyyy, nyyy = myTools.Get_TProfile(ar_lambdaRF,ar_flux, lambdaRFTemplateBinEdges__, ar_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$flux$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### delta vs. lambda_RF
	xxx, yyy, eyyy, nyyy = myTools.Get_TProfile(ar_lambdaObs,ar_flux, lambdaObsBinEdges__, ar_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
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
	plt.xlabel(r'$z_{qso}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()

	return
def saveCatalogueQSO():
	'''
	'''

	path = '/home/nfs/manip/mnt/bao/legoff/spectra-1280.fits'
	cat = pyfits.open(path, memmap=True)

	xx = []
	yy = []
	zz = []

	for el in cat[1:]:
		xx += [ el.header['X']    ]
		yy += [ el.header['Y']    ]
		zz += [ el.header['ZQSO'] ]

	col_xx              = pyfits.Column(name='X',  format='D', array=xx, unit='Mpc/h')
	col_yy              = pyfits.Column(name='Y',  format='D', array=yy, unit='Mpc/h')
	col_zz              = pyfits.Column(name='Z',  format='D', array=zz, unit='0 (redshift)')
	
	tbhdu = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz])
	tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/mock_QSO_catalogue_1280_0_0.fits', clobber=True)
	
	return
def haveAlookQSO():	

	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMc/mock_QSO_catalogue_1280_0_0.fits', memmap=True)[1].data
	print cat.size
	print numpy.amin(cat['X']), numpy.amax(cat['X'])
	print numpy.amin(cat['Y']), numpy.amax(cat['Y'])
	print numpy.amin(cat['Z']), numpy.amax(cat['Z'])

	### Map
	plt.plot(cat['X'],cat['Y'], linestyle="", marker="o")
	plt.xlabel(r'$X [Mpc]$')
	plt.ylabel(r'$Y [Mpc]$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution redshift
	plt.hist(cat['Z'], bins=50)
	plt.xlabel(r'$z_{qso}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()

	return

#copyToNormalFits()
#allNeededFileds()
#saveCatalogueQSO()
#haveAlookQSO()
haveAlookForest()
