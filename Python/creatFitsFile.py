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
from iminuit import Minuit
import matplotlib.pyplot as plt
import re
import cosmolopy.distance as cosmology
#import array
#import matplotlib.patches as mpatches
#import decimal
#import profile
#import warnings
#warnings.filterwarnings("error")

### My tools
import myTools
from myTools import Get_TProfile
### The Constants
from delta_const import *

### "HOME" or "ICLUST"
location__ = 'ICLUST'
### "DR12" or 'Guy' or 'Margala' or "MOCK" or 'eBOSS'
pipeline__ = 'DR12'
### False or True
reObs__ = False
### What type of FITS file will we use: 'spPlate' or 'spec'
dataFitsType__ = 'spec'


def make_all_Fits(iStart=0,iEnd=-1):

	tmp_string = "\n   ---- Starting ---- \n\n  iStart = " + str(iStart) + "  iEnd = " + str(iEnd) + "\n"

	tmp_command = "echo \"" + tmp_string + "\""
	subprocess.call(tmp_command, shell=True)

	### Path to the folder where all spectra are
	pathToSpec = Get_Path_To_Fits()
	
	### Get the catalogue
	data, plate_list = Get_Catalogue()

	sizeMax = data[:,0].size

	print '  sizeMax = ', sizeMax

	if (iEnd!=-1):
		if (iEnd>sizeMax):
			iEnd=sizeMax
		if (iStart>sizeMax):
			print '  iStart>sizeMax'
			return
		iEnd += 1
		data = data[iStart:iEnd]
	
	sizeMax = data[:,0].size

	tmp_command = "echo \"" + "  Final length = "+ str(sizeMax) + "\""
	subprocess.call(tmp_command, shell=True)

	tmp_command = "echo  \"" + "\n  Creating the FITS file" + "\""
	subprocess.call(tmp_command, shell=True)

	### Create a FITS file with only the usefull data
	plate                = pyfits.Column(name='PLATE',                 format='J', array=data[:,0] )
	mjd                  = pyfits.Column(name='MJD',                   format='J', array=data[:,1] )
	fiber                = pyfits.Column(name='FIBERID',               format='J', array=data[:,2] )
	
	ra                   = pyfits.Column(name='RA',                    format='D', array=data[:,3], unit='deg')
	de                   = pyfits.Column(name='DEC',                   format='D', array=data[:,4], unit='deg')
	zz                   = pyfits.Column(name='Z_VI',                  format='D', array=data[:,5] )
	col_normFlux         = pyfits.Column(name='NORM_FACTOR',           format='D', array=numpy.zeros(sizeMax), unit='Flux')
	
	tmp_nbBinForest     = str(nbBinRFMax__)+'D'
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',  format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Angstrom')
	lambdaRFForest      = pyfits.Column(name='LAMBDA_RF',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Angstrom')
	normFluxForest      = pyfits.Column(name='FLUX',        format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Flux' )
	normFluxIvarForest  = pyfits.Column(name='FLUX_IVAR',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Flux' )
		
	del data
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, ra, de, zz, col_normFlux, lambdaForest, lambdaRFForest, normFluxForest, normFluxIvarForest])
	cat_tbhdu = tbhdu.data
	
	tmp_command = "echo  \"" + "\n  Starting the loop\n\n" + "\""
	subprocess.call(tmp_command, shell=True)
	
	if (dataFitsType__=='spec'):

		### Read and write in the FITS file
		for el in cat_tbhdu:
	
			try:
				#cat = pyfits.open(pathToSpec + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
				cat = pyfits.open(pathToSpec + str(el['PLATE']) + "/spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
				### For DR12 Margala
				#cat = pyfits.open(pathToSpec + str(el['PLATE']) + "/corrected-spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
				### For mocks
				#cat = pyfits.open(pathToSpec + str(el['PLATE']) + "/mock-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
			except Exception,error:
				#tmp_string = pathToSpec + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
				tmp_string = pathToSpec + str(el['PLATE']) + "/spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
				#tmp_string = pathToSpec + str(el['PLATE']) + "/corrected-spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
				#tmp_string = pathToSpec + str(el['PLATE']) + "/mock-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
	
				tmp_command = "echo  \"" + "  File not found \n " + tmp_string + "\""
				subprocess.call(tmp_command, shell=True)
				tmp_command = "echo  \"  " + str(error)  + "\n\""
				subprocess.call(tmp_command, shell=True)
				continue
				
			### Apply cuts for bad pixels
			############################## 
			
			### CCD and too many sky lines
			cat = cat[ (numpy.logical_and( (cat["LOGLAM"]>=log10lambdaObsMin__) , (cat["LOGLAM"]<log10lambdaObsMax__) )) ]
			### Sky Lines
			for lines in skyLines__:
				cat = cat[ (numpy.logical_or( (cat["LOGLAM"]<=lines[0]) , (cat["LOGLAM"]>=lines[1]) )) ]
			cat = cat[ numpy.logical_and( numpy.logical_and( (cat["IVAR"]>0.), (cat["AND_MASK"]<bit16__)), (numpy.isfinite(cat["FLUX"])) ) ]
		
			### Get where to cut for lambda_RF_Norma
			tmp_logZ = numpy.log10(1.+el['Z_VI'])
			
			### Get the normalisation factor
			try:
				normFactor = numpy.mean( cat["FLUX"][ numpy.logical_and( (cat["LOGLAM"]>log10lambdaRFNormaMin__+tmp_logZ), (cat["LOGLAM"]<log10lambdaRFNormaMax__+tmp_logZ) ) ] )
			except Exception,error:
				tmp_command = "echo  \"" + "  Error in: 'normFactor = numpy.mean',  z = " + str(el['Z_VI']) + ", idx = " + str(cat_tbhdu[ (cat_tbhdu['NORM_FACTOR']!=0.) ].size+iStart) + "\""
				subprocess.call(tmp_command, shell=True)
				tmp_command = "echo  \"  " + str(error)  + "\n\""
				subprocess.call(tmp_command, shell=True)
				continue
	
			el['NORM_FACTOR'] = normFactor
		
			### Apply cuts to keep only the forest
			cat = cat[ numpy.logical_and( (cat["LOGLAM"]>=log10lambdaRFTemplateMin__+tmp_logZ) , (cat["LOGLAM"]<log10lambdaRFTemplateMax__+tmp_logZ) ) ]
			
			### Find the number of pixels
			tmp_lenCat = cat.size
			
			### Store the data
			el['FLUX'][:tmp_lenCat]       = cat["FLUX"]
			el['FLUX_IVAR'][:tmp_lenCat]  = cat["IVAR"]
			el['LAMBDA_OBS'][:tmp_lenCat] = cat["LOGLAM"]
	
	print "  Number of forest (init)                                : " + str(cat_tbhdu.size)
	
	### Set values of columns
	cat_tbhdu['LAMBDA_OBS'][ (cat_tbhdu['LAMBDA_OBS'] != 0.) ]   = numpy.power(10., cat_tbhdu['LAMBDA_OBS'][ (cat_tbhdu['LAMBDA_OBS'] != 0.) ])
	cat_tbhdu['LAMBDA_RF']                                       = numpy.dot(numpy.diag(1./(1.+cat_tbhdu['Z_VI'])),cat_tbhdu['LAMBDA_OBS'])

	### Remove forest with less than nbBinRFMin__ pixels in the forest
	cut_noForestPixel = numpy.logical_and( (cat_tbhdu['LAMBDA_RF']>=lambdaRFMin__) , (cat_tbhdu['LAMBDA_RF']<lambdaRFMax__) ).astype(int)
	len_forest = numpy.sum(cut_noForestPixel,axis=1)
	cat_tbhdu = cat_tbhdu[ (len_forest>=nbBinRFMin__) ]

	print "  Number of forest (len_forest>=nbBinRFMin__)            : " + str(cat_tbhdu.size)
	print

	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
	tbhdu.update()

	path = Get_Path_To_Save_New_Fits(iStart,iEnd)
	tbhdu.writeto(path, clobber=True)

	print "\n\n\n"
	
def Merge_Files():
	
	where = Get_Data_Fits_File()
	folder = where[0]
	path   = where[1]
	#path   = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_reObs/DR12_reObs.fits'
	#folder = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_reObs/'
	path   = "/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_primery/DR12_primery.fits"
	folder = "/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_primery/"

	if (pipeline__=="DR12" or pipeline__=='Guy' or pipeline__=='Margala'):
		if (reObs__):
			scheme = 'DR12_reObs_'
		else:
			scheme = 'allDR12_test_'
	elif (pipeline__=='MOCK'):
		scheme = "Mock__DR11_"
	lenScheme = len(scheme)
	
	### Get the list of files
	all_t = os.listdir(folder)	
	tmp_all_t = []
	for el in all_t:
                if (el[:lenScheme]==scheme):
                        tmp_all_t.append(el)
	all_t = tmp_all_t

	### Sort the list of files
	convert      = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	all_t = sorted(all_t, key = alphanum_key)

	### Get the FitsFile
	all_t_file  = []
	for el in all_t:
		all_t_file.append( pyfits.open(folder+el,  memmap=True) )
		print folder+el
	
	### Get arrays of size of FitsFile
	all_t_nrows = []
	for el in all_t_file:
		all_t_nrows.append( el[1].data.shape[0] )
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

	### Add some column to the data
	cat     = hdu.data
	sizeMax = cat.size

	### Remove forest with normFactor <= 0
	cat = cat[ (cat['NORM_FACTOR']>0.) ]
	sizeMax = cat.size
	print "  Number of forest (normFactor > 0)                      : " + str(sizeMax)
	
	### Set values of columns
	step = 500
	iStart = 0
	iEnd   = step
	while (iStart<=sizeMax):
		print iStart, iEnd
		cat['FLUX'][iStart:iEnd]      = numpy.dot(numpy.diag(1./cat['NORM_FACTOR'][iStart:iEnd]),cat['FLUX'][iStart:iEnd])
		cat['FLUX_IVAR'][iStart:iEnd] = numpy.dot(numpy.diag(cat['NORM_FACTOR'][iStart:iEnd]**2.),cat['FLUX_IVAR'][iStart:iEnd])
		iStart = iEnd+1
		iEnd   += step
		if (iEnd>sizeMax):
			iEnd = sizeMax
	
	cut_noForestPixel = numpy.logical_and( (cat['LAMBDA_RF']>=lambdaRFMin__) , (cat['LAMBDA_RF']<lambdaRFMax__) ).astype(int)
	len_forest = numpy.sum(cut_noForestPixel,axis=1)
	mean_lambdaRF_in_forest = numpy.sum(cut_noForestPixel*cat['LAMBDA_RF'],axis=1)/len_forest

	### Create a FITS file with only the usefull data
	plate                = pyfits.Column(name='PLATE',           format='J', array=cat['PLATE'] )
	mjd                  = pyfits.Column(name='MJD',             format='J', array=cat['MJD'] )
	fiber                = pyfits.Column(name='FIBERID',         format='J', array=cat['FIBERID'] )
	
	ra                   = pyfits.Column(name='RA',              format='D', array=cat['RA'],     unit='deg' )
	de                   = pyfits.Column(name='DEC',             format='D', array=cat['DEC'],    unit='deg' )
	zz                   = pyfits.Column(name='Z_VI',            format='D', array=cat['Z_VI'] )
	nb                   = pyfits.Column(name='NB_PIXEL',        format='I', array=len_forest )
	meanLambdaRF         = pyfits.Column(name='MEAN_FOREST_LAMBDA_RF', format='F', array=mean_lambdaRF_in_forest, unit='angstrom')

	tmp_nbBinForest     = str(nbBinRFMax__)+'D'
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',       format=tmp_nbBinForest, array=cat['LAMBDA_OBS'], unit='angstrom' )
	lambdaRFForest      = pyfits.Column(name='LAMBDA_RF',        format=tmp_nbBinForest, array=cat['LAMBDA_RF'],  unit='angstrom' )
	normFluxForest      = pyfits.Column(name='NORM_FLUX',        format=tmp_nbBinForest, array=cat['FLUX'] )
	normFluxIvarForest  = pyfits.Column(name='NORM_FLUX_IVAR',   format=tmp_nbBinForest, array=cat['FLUX_IVAR'] )
	deltaForest         = pyfits.Column(name='DELTA',            format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)) )
	deltaIvarForest     = pyfits.Column(name='DELTA_IVAR',       format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)) )
	deltaWeight         = pyfits.Column(name='DELTA_WEIGHT',     format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)) )
	template            = pyfits.Column(name='TEMPLATE',         format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)) )	
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, ra, de, zz, nb, meanLambdaRF, lambdaForest, lambdaRFForest, normFluxForest, normFluxIvarForest, deltaForest, deltaIvarForest, deltaWeight, template])
	cat_tbhdu = tbhdu.data
	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)

	## Map
	#plt.plot(cat_tbhdu["RA"], cat_tbhdu["DEC"], linestyle="", marker="o")
	#plt.show()
	
	tbhdu.update()
	tbhdu.writeto(path, clobber=True)
	#tbhdu.writeto('/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_primery/DR12_primery_forAnalys.fits', clobber=True)

	return
def Get_Only_Needed_Columns():
	'''
	'''
	
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_primery/DR12_primery_forAnalys.fits'
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_primery/DR12_primery.fits'
	print path
	file_cat = pyfits.open(path)
	cat = file_cat[1].data
	sizeMax = cat.size
	print sizeMax


	### Remove useless pixels
	cat['DELTA'][ numpy.logical_or( (cat['LAMBDA_RF']<lambdaRFMin__) , (cat['LAMBDA_RF']>=lambdaRFMax__) ) ] = 0.
	cat['DELTA_WEIGHT'][ numpy.logical_or( (cat['LAMBDA_RF']<lambdaRFMin__) , (cat['LAMBDA_RF']>=lambdaRFMax__) ) ] = 0.
	cat['DELTA'][ (cat['DELTA_IVAR']<=minDeltaIvar__) ] = 0.
	cat['DELTA_WEIGHT'][ (cat['DELTA_IVAR']<=minDeltaIvar__) ] = 0.

	### Get the Spherical and cartesian coordinate

	### Get an array to convert redshift into distance
	cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}
	distfunc, redfunc = cosmology.quick_distance_function(cosmology.comoving_distance, zmin=1.5, zmax=10., zstep=0.001, return_inverse=True, **cosmo)

	ar_ra = numpy.ones((sizeMax,nbBinRFMax__))
	ar_de = numpy.ones((sizeMax,nbBinRFMax__))

	### Create arays of R.A. and Dec.
	step = 500
	iStart = 0
	iEnd   = step
	while (iStart<=sizeMax):
		print iStart, iEnd
		ar_ra[iStart:iEnd] = numpy.dot(numpy.diag(cat['RA'][iStart:iEnd]),ar_ra[iStart:iEnd] )
		ar_de[iStart:iEnd] = numpy.dot(numpy.diag(cat['DEC'][iStart:iEnd]),ar_de[iStart:iEnd] )
		iStart = iEnd+1
		iEnd   += step
		if (iEnd>sizeMax):
			iEnd = sizeMax
	
	### Spherical coordinate
	ar_ra *= numpy.pi/180.
	ar_de  = (90.-ar_de)*numpy.pi/180.
	ar_rr  = distfunc(cat['LAMBDA_OBS']/lambdaRFLya__-1.)*h__

	### Cartesian coordinate
	ar_xxx = ar_rr*numpy.sin(ar_de)*numpy.cos(ar_ra)
	ar_yyy = ar_rr*numpy.sin(ar_de)*numpy.sin(ar_ra)
	ar_zzz = ar_rr*numpy.cos(ar_de)



	### Create a FITS file with only the usefull data
	plate                = pyfits.Column(name='PLATE',           format='J', array=cat['PLATE'] )
	mjd                  = pyfits.Column(name='MJD',             format='J', array=cat['MJD'] )
	fiber                = pyfits.Column(name='FIBERID',         format='J', array=cat['FIBERID'] )
	
	ra                   = pyfits.Column(name='RA',              format='D', array=cat['RA'],     unit='deg' )
	de                   = pyfits.Column(name='DEC',             format='D', array=cat['DEC'],    unit='deg' )
	zz                   = pyfits.Column(name='Z_VI',            format='D', array=cat['Z_VI'] )
	nb                   = pyfits.Column(name='NB_PIXEL',        format='I', array=cat['NB_PIXEL'] )

	tmp_nbBinForest     = str(nbBinRFMax__)+'D'
	col_rrr             = pyfits.Column(name='R_SPHERICAL',      format=tmp_nbBinForest, array=ar_rr,  unit='h^-1 Mpc')
	col_xxx             = pyfits.Column(name='X_CARTESIAN',      format=tmp_nbBinForest, array=ar_xxx, unit='h^-1 Mpc')
	col_yyy             = pyfits.Column(name='Y_CARTESIAN',      format=tmp_nbBinForest, array=ar_yyy, unit='h^-1 Mpc')
	col_zzz             = pyfits.Column(name='Z_CARTESIAN',      format=tmp_nbBinForest, array=ar_zzz, unit='h^-1 Mpc')
	col_delta           = pyfits.Column(name='DELTA',            format=tmp_nbBinForest, array=cat['DELTA'] )
	col_deltaWeight     = pyfits.Column(name='DELTA_WEIGHT',     format=tmp_nbBinForest, array=cat['DELTA_WEIGHT'] )
	
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, ra, de, zz, nb, col_rrr, col_xxx, col_yyy, col_zzz, col_delta, col_deltaWeight])
	cat_tbhdu = tbhdu.data
	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
	tbhdu.update()
	tbhdu.writeto('test2.fits', clobber=True)
	file_cat.close()
	#return


	### delta vs. lambda_RF
	tmp_a = []
	for i in range(0,cat_tbhdu.size):
		tmp_a += [ numpy.average(cat['DELTA'][i][ (cat['DELTA_WEIGHT'][i]>0.) ], weights=cat['DELTA_WEIGHT'][i][ (cat['DELTA_WEIGHT'][i]>0.) ]) ]

	plt.hist(tmp_a)
	plt.show()

	### Test if every thing is O.K.
	
	ar_cut              = (cat['DELTA_WEIGHT']>0.)
	ar_lambdaObs        = cat['LAMBDA_OBS'][ ar_cut ]
	ar_lambdaRF         = cat['LAMBDA_RF'][ ar_cut ]
	ar_delta            = cat['DELTA'][ ar_cut ]
	ar_delta_weight     = cat['DELTA_WEIGHT'][ ar_cut ]
	ar_z_pixel          = cat['LAMBDA_OBS'][ ar_cut ]/lambdaRFLya__-1.
	ar_rara             = ar_ra[ ar_cut ]*180./numpy.pi
	ar_dede             = 90.-ar_de[ ar_cut ]*180./numpy.pi
	ar_rrrr             = ar_rr[ ar_cut ]
	ar_xxxx             = ar_xxx[ ar_cut ]
	ar_yyyy             = ar_yyy[ ar_cut ]
	ar_zzzz             = ar_zzz[ ar_cut ]

	print numpy.average(ar_delta, weights=ar_delta_weight)


	### delta vs. lambda_RF
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_lambdaRF, ar_delta, lambdaRFTemplateBinEdges__, ar_delta_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o", label=" - ")
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### delta vs. lambda_Obs
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_lambdaObs, ar_delta, lambdaObsBinEdges__, ar_delta_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o", label=" - ")
	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### delta vs. z_pixel
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_z_pixel, ar_delta, 300, ar_delta_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o", label=" - ")
	plt.xlabel(r'$z_{pixel}$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### Distrivution of forest in the sky
	plt.plot(cat['RA'], cat['DEC'], linestyle="", marker="o", label=r'$ map of forest $')
	plt.xlabel(r'$R. A. \, [^\circ]$', fontsize=40)
	plt.ylabel(r'$Dec. \, [^\circ]$', fontsize=40)
	#plt.xlim([0.,360.])
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### delta vs. dist_rrrr
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_rrrr, ar_delta, 300, ar_delta_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o", label=" - ")
	plt.xlabel(r'$r \, [h^{-1}. Mpc]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### delta vs. dist_xxxx
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_xxxx, ar_delta, 300, ar_delta_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o", label=" - ")
	plt.xlabel(r'$x \, [h^{-1}. Mpc]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### delta vs. dist_yyyy
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_yyyy, ar_delta, 300, ar_delta_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o", label=" - ")
	plt.xlabel(r'$y \, [h^{-1}. Mpc]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### delta vs. dist_zzzz
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_zzzz, ar_delta, 300, ar_delta_weight)
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o", label=" - ")
	plt.xlabel(r'$z \, [h^{-1}. Mpc]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### delta vs. (ra,dec)
	xxx, yyy, zzz = myTools.Get_2DTProfile(ar_rara, ar_dede, ar_delta, 100, 100,ar_delta_weight)
	plt.imshow(xxx,origin='lower',extent=[0., 10., 0., 10.],interpolation='None')
	plt.xlabel(r'$R.A. \, [^\circ]$', fontsize=40)
	plt.ylabel(r'$Dec. \, [^\circ]$', fontsize=40)
	plt.grid(True)
	cbar = plt.colorbar()
	cbar.set_label(r'$<\delta >$',size=40)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
	
	file_cat.close()
	
	return
def Get_Template():
	'''
	'''
	
	### Number of loops
	nbLoop = 5
	
	### Init values for the weight
	init_eta = 1.
	init_sigma2LSS = 0.1

	### Print that it is starting
	tmp_command = "echo  \"  " +  "  ---- Begin ---- " + "\n\""
	subprocess.call(tmp_command, shell=True)
	
	### Get Data
	### ,mode='update' ##, memmap=True
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_primery/DR12_primery.fits'
	print path
	file_cat = pyfits.open(path,mode='update')
	print '  Nb total of forest    = ', file_cat[1].data.size
	cat = file_cat[1].data
	print '  Nb selected of forest = ', cat.size
	
	### Init the template and corrections
	eta               =  [ [init_eta] ]
	sigma2LSS         =  [ [init_sigma2LSS] ]
	centerRedshift    =  [ [3.] ]
	
	### Init the different needed plots (ex: template)
	coef_weight         = numpy.power(cat['LAMBDA_OBS']/(lambdaRFLya__*onePlusZ0__), halfGama__)
	cat['DELTA_WEIGHT'] = coef_weight/(sigma2LSS[0][0]+1./(eta[0][0]*cat['NORM_FLUX_IVAR']))
	ar_cut              = (cat['NORM_FLUX_IVAR']>minDeltaIvar__)
	ar_lambdaObs        = cat['LAMBDA_OBS'][ ar_cut ]
	ar_lambdaRF         = cat['LAMBDA_RF'][ ar_cut ]
	ar_flux             = cat['NORM_FLUX'][ ar_cut ]
	ar_weight           = cat['DELTA_WEIGHT'][ ar_cut ]
	
	### Template in 1D
	xxx, yyy, eyyy, nyyy = Get_TProfile(ar_lambdaRF,ar_flux, lambdaRFTemplateBinEdges__,ar_weight)

	sizeObs    = lambdaObsBinCenter__.size

	template          = [ (xxx,  yyy,   eyyy) ]
	residualRF        = [ (xxx,  numpy.zeros(xxx.size),  numpy.zeros(xxx.size)) ]
	meanTransRF       = [ (xxx,  numpy.ones(xxx.size),   numpy.zeros(xxx.size)) ]
	residualObs       = [ (lambdaObsBinCenter__, numpy.zeros(sizeObs), numpy.zeros(sizeObs)) ]
	meanTrandObs      = [ (lambdaObsBinCenter__, numpy.ones(sizeObs),  numpy.zeros(sizeObs)) ]
	varPipeline       = [ (deltaIvarBinCenter__, numpy.power(sigma2LSS[0]+1./(eta[0]*numpy.power(10,deltaIvarBinCenter__)),-1.), numpy.zeros(deltaIvarBinCenter__.size)) ]
	meanDelta = numpy.zeros(nbLoop+1)
	stdDelta  = numpy.zeros(nbLoop+1)
	
	
	### Template in 2D
	#mean, err, number = myTools.Get_2DTProfile(ar_lambdaRF,ar_lambdaObs/ar_lambdaRF-1.,ar_flux,lambdaRFTemplateBinEdges__,redshiftBinEdges__,ar_weight)
	#template2D = [ (mean, err, number) ]
	template2D    = [ (numpy.ones(1),numpy.ones(1), numpy.ones(1)) ]

	### Weight in 2D
	varPipeline2D     = [ (numpy.ones(1),numpy.ones(1), numpy.ones(1)) ]
	
	for loopIdx in numpy.arange(0,nbLoop):

		tmp_command = "echo  \"  ---- " + str(loopIdx) + " ---- \n\" "
		subprocess.call(tmp_command, shell=True)

		cat['TEMPLATE']     = numpy.interp(cat['LAMBDA_RF'], template[loopIdx][0], template[loopIdx][1])

		### Find delta of catalogue
		for el in cat:

			ar_cut        = numpy.logical_and( numpy.logical_and( (el['LAMBDA_RF']>=lambdaRFMin__) , (el['LAMBDA_RF']<lambdaRFMax__) ), (el['DELTA_IVAR']>minDeltaIvar__))
			tmp_lambda_rf = el['LAMBDA_RF'][ ar_cut ]/el['MEAN_FOREST_LAMBDA_RF']-1.
			tmp_norm_flux = el['NORM_FLUX'][ ar_cut ]
			tmp_template  = el['TEMPLATE'][ ar_cut ]
			tmp_weight    = el['DELTA_WEIGHT'][ ar_cut ]
	
			def chi2(alpha,beta):
				return numpy.sum(numpy.power((tmp_norm_flux-(alpha+beta*tmp_lambda_rf)*tmp_template),2.)*tmp_weight)

			m = Minuit(chi2, alpha=1.,error_alpha=1., beta=0.,error_beta=1.,print_level=-1, pedantic=False)
			m.migrad()

			el['TEMPLATE'] *= (m.values['alpha']+m.values['beta']*(el['LAMBDA_RF']/el['MEAN_FOREST_LAMBDA_RF']-1.) )

		### Set arrays
		################################################################

		cat['TEMPLATE']    *= numpy.interp(cat['LAMBDA_OBS'], meanTrandObs[loopIdx][0], meanTrandObs[loopIdx][1])*numpy.interp(cat['LAMBDA_RF'], meanTransRF[loopIdx][0], meanTransRF[loopIdx][1])
		cat['DELTA']        = cat['NORM_FLUX']/cat['TEMPLATE']-1.
		cat['DELTA_IVAR']   = cat['NORM_FLUX_IVAR']*cat['TEMPLATE']*cat['TEMPLATE']
		

		ar_cut          = numpy.logical_and((cat['NORM_FLUX_IVAR']>0.),(cat['DELTA_IVAR']>minDeltaIvar__))
		ar_lambdaObs    = cat['LAMBDA_OBS'][ ar_cut ]
		ar_lambdaRF     = cat['LAMBDA_RF'][ ar_cut ]
		ar_flux         = cat['NORM_FLUX'][ ar_cut ]
		ar_delta        = cat['DELTA'][ ar_cut ]
		ar_weight       = cat['DELTA_WEIGHT'][ ar_cut ]
		log_delta_ivar  = numpy.log10(cat['DELTA_IVAR'][ar_cut])
	
		### 1D: Template
		xxx, yyy, eyyy, nyyy = Get_TProfile(ar_lambdaRF,ar_flux, lambdaRFTemplateBinEdges__,ar_weight)
		template.append( (xxx,yyy,eyyy) )
	
		### 2D: Template
		#mean, err, number = myTools.Get_2DTProfile(ar_lambdaRF,ar_lambdaObs/ar_lambdaRF-1.,ar_flux,lambdaRFTemplateBinEdges__,redshiftBinEdges__,ar_weight)
		#template2D.append( (lambdaRFTemplateBinCenter__,redshiftBinCenter__,mean)  )
		
		### Delta vs. lambda_RF
		xxx, yyy, eyyy, nyyy = Get_TProfile(ar_lambdaRF,ar_delta, lambdaRFTemplateBinEdges__,ar_weight)
		residualRF.append( (xxx,yyy,eyyy) )
		
		### Delta+1 vs. lambda_RF
		if (loopIdx%2==1):
			meanTransRF.append( (xxx, (yyy+1.)*numpy.interp(xxx,meanTransRF[loopIdx][0], meanTransRF[loopIdx][1]), eyyy) )
		else:
			meanTransRF.append( meanTransRF[loopIdx] )

		### Remove delta at the edge of the forest
		ar_cut          = numpy.logical_and( (ar_lambdaRF>=lambdaRFMin__) , (ar_lambdaRF<lambdaRFMax__) )
		ar_lambdaObs    = ar_lambdaObs[ ar_cut ]
		ar_delta        = ar_delta[ ar_cut ]
		ar_weight       = ar_weight[ ar_cut ]
		log_delta_ivar  = log_delta_ivar[ ar_cut ]

		### Delta vs. lambda_Obs
		xxx, yyy, eyyy, nyyy = Get_TProfile(ar_lambdaObs,ar_delta, lambdaObsBinEdges__,ar_weight)
		residualObs.append( (xxx,yyy,eyyy) )
		
		
		### Delta+1 vs. lambda_Obs
		if (loopIdx%2==0):
			meanTrandObs.append( (xxx, (yyy+1.)*numpy.interp(xxx,meanTrandObs[loopIdx][0], meanTrandObs[loopIdx][1]), eyyy) )
		else:
			meanTrandObs.append( meanTrandObs[loopIdx] )

		
		### Variance in 2D
		mean, err, number = myTools.Get_2DTProfile(log_delta_ivar,ar_lambdaObs/lambdaRFLya__-1.,ar_delta,deltaIvarBinEdges__,redshiftWeightBinEdges__,ar_weight)

		mean  = 1./( (err**2.)*number)
		err   = mean/numpy.sqrt(number)
		varPipeline2D.append( (mean, err, number) )

		TMPcenterRedshift = []
		TMPsigma2LSS      = []
		TMPeta            = []
		
		for i in range(0,redshiftWeightBinCenter__.size):

			tmp_cut  = numpy.logical_and( (number[:,i]>2), (numpy.logical_and( (deltaIvarBinCenter__>-1.5), (deltaIvarBinCenter__<2.3) )) )
			tmp_xxx  = deltaIvarBinCenter__[ tmp_cut ]
			tmp_yyy  = mean[:,i][ tmp_cut ]
			tmp_eyyy = err[:,i][ tmp_cut ]
			
			if ( tmp_xxx.size < 15.):
				continue
			else:
				### Set minuit
				def chi2(tmp_TMPeta,tmp_TMPsigma2LSS):
					return numpy.sum( numpy.power((tmp_yyy-1./(tmp_TMPsigma2LSS+1./(tmp_TMPeta*numpy.power(10,tmp_xxx))))/tmp_eyyy ,2.) )
				m = Minuit(chi2, tmp_TMPeta=eta[0][0], error_tmp_TMPeta=0.1,  tmp_TMPsigma2LSS=sigma2LSS[0][0],limit_tmp_TMPsigma2LSS=(0.001,1000.), error_tmp_TMPsigma2LSS=0.1, print_level=-1, errordef=0.01) 	
				m.migrad()

				TMPcenterRedshift += [ redshiftWeightBinCenter__[i] ]
				TMPsigma2LSS      += [ m.values['tmp_TMPsigma2LSS'] ]
				TMPeta            += [ m.values['tmp_TMPeta'] ]
		
		centerRedshift += [ TMPcenterRedshift ]
		sigma2LSS      += [ TMPsigma2LSS ]
		eta            += [ TMPeta ]

		### Mean and variance
		meanDelta[loopIdx+1] =  numpy.average(ar_delta,weights=ar_weight )
		stdDelta[loopIdx+1]  =  numpy.average( (ar_delta-meanDelta[loopIdx+1])*(ar_delta-meanDelta[loopIdx+1]), weights=ar_weight )

		cat['DELTA_WEIGHT'] = coef_weight/(numpy.interp(cat['LAMBDA_OBS']/lambdaRFLya__-1., centerRedshift[loopIdx],sigma2LSS[loopIdx])+1./(numpy.interp(cat['LAMBDA_OBS']/lambdaRFLya__-1.,centerRedshift[loopIdx],eta[loopIdx])*cat['DELTA_IVAR']))



	data  = cat['DELTA_WEIGHT'][ (cat['DELTA_WEIGHT']>0.) ]
	data2 = cat['DELTA_IVAR'][ (cat['DELTA_IVAR']>0.) ]
	print '  Nb of pixels            = ', data.size
	print '  Nb of Forest            = ', cat.size
	print '  Nb of pixels per forest = ', 1.*data.size/cat.size
	print data2.size
	
	file_cat.close()
	


	### Save results
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Code/Python/Results/'
	aloopIdx = numpy.arange(0,nbLoop+1)
	numpy.savetxt(path + 'meanDelta.txt', zip(aloopIdx, meanDelta))
	numpy.savetxt(path + 'stdDelta.txt',  zip(aloopIdx, stdDelta))

	numpy.save(path + 'eta',       zip(centerRedshift,eta))
	numpy.save(path + 'sigma2LSS', zip(centerRedshift,sigma2LSS))
	numpy.save(path + 'template',          template)
	numpy.save(path + 'residualRF',        residualRF)
	numpy.save(path + 'meanTransRF',       meanTransRF)
	numpy.save(path + 'residualObs',       residualObs)
	numpy.save(path + 'meanTrandObs',      meanTrandObs)
	numpy.save(path + 'varPipeline',       varPipeline )

	numpy.save(path + 'template2D',        template2D)
	numpy.save(path + 'varPipeline2D',     varPipeline2D)
	
	print '  < delta >    = ', meanDelta[-1]
	print '  var delta    = ', stdDelta[-1]
	print '  eta          = ', eta[-1]
	print '  sigma^2_LSS  = ', sigma2LSS[-1]
	
	return
def Plot_Results():
	'''
	'''
	
	print
	print "  ---- Plots ----"
	print
	
	
	### Show only the last step or not
	lastStep = True
	minStep  = 0

	path = '/home/gpfs/manip/mnt/bao/hdumasde/Code/Python/Results/'
	#path = '/home/helion/Documents/Thèse/Results/RootFile/Results_backup_allDR12/'

	### Cross-correlation
	
	def plot_data():
		loopIdx = 0
		for el in data:
			if (loopIdx>=minStep):
				plt.errorbar(el[0], el[1] ,yerr=el[2], marker="o", label="step = " + str(loopIdx) )
			loopIdx += 1

	
	### Evolution of < delta > vs. loopIdx
	data = numpy.loadtxt(path+'meanDelta.txt')
	plt.errorbar(data[:,0], numpy.fabs(data[:,1]), marker="o", label="|< delta >|")
	plt.errorbar(data[:,0], data[:,1], marker="o", label="< delta >")
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$step$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	print '  < delta >    = ', data[-1,1]
	
	### Evolution of std_delta vs. loopIdx
	data = numpy.loadtxt(path+'stdDelta.txt')
	plt.errorbar(data[:,0], numpy.fabs(data[:,1]), marker="o", label="-")
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$step$', fontsize=40)
	plt.ylabel(r'$\sigma_{\delta}$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	print '  var delta    = ', data[-1,1]
	
	### eta vs. loopIdx
	data = numpy.load(path+'eta.npy')
	loopIdx = 0
	for el in data:
		if (loopIdx>=minStep):
			plt.plot(el[0], el[1], marker="o", label="step = " + str(loopIdx) )
		loopIdx += 1
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$z$', fontsize=40)
	plt.ylabel(r'$\eta$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	print '  eta          = ', data[-1,1]
	
	### sigma^2_LSS vs. loopIdx
	data = numpy.load(path+'sigma2LSS.npy')
	loopIdx = 0
	for el in data:
		if (loopIdx>=minStep):
			plt.plot(el[0], el[1], marker="o", label="step = " + str(loopIdx) )
		loopIdx += 1
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$z$', fontsize=40)
	plt.ylabel(r'$\sigma_{L.S.S.}^{2}$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	print '  sigma^2_LSS  = ', data[-1,1]

	### Flux vs. lambda_RF
	data = numpy.load(path+'template.npy')
	plot_data()
	plt.title(r'$Template$', fontsize=40)
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$Norm \, flux$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### Delta vs. lambda_RF
	data = numpy.load(path+'residualRF.npy')
	plot_data()
	plt.title(r'$Residual$', fontsize=40)
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### Delta+1 vs. lambda_RF
	data = numpy.load(path+'meanTransRF.npy')
	plot_data()
	plt.title(r'$Correction$', fontsize=40)
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$< \delta + 1 >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### Delta vs. lambda_Obs
	data = numpy.load(path+'residualObs.npy')
	plot_data()
	plt.title(r'$Residual$', fontsize=40)
	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$< \delta >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### Delta+1 vs. lambda_Obs
	data = numpy.load(path+'meanTrandObs.npy')
	plot_data()
	plt.title(r'$Correction$', fontsize=40)
	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$< \delta +1>$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	### varPipeline
	data       = numpy.load(path+'varPipeline2D.npy')
	data = data[1:]
	loopIdx = 0
	for el in data:
		if (loopIdx>=minStep):
			for i in range(0,redshiftWeightBinCenter__.size):
				cut = (el[2][:,i]>1)
				xxx = deltaIvarBinCenter__[ cut ]
				yyy = el[0][:,i][ cut ]
				eyy = el[1][:,i][ cut ]
				if (xxx.size>0):
					plt.errorbar(numpy.power(10,xxx), yyy, yerr=eyy, marker="o")
			loopIdx += 1
		plt.title(r'', fontsize=40)
		plt.xlabel(r'$\sigma_{pip.}^{-2}$', fontsize=40)
		plt.ylabel(r'$\sigma_{\delta}^{-2}$', fontsize=40)
		myTools.deal_with_plot(True, True,True)
		plt.show()

	
	'''
	data       = numpy.load(path+'template2D.npy')
	loopIdx = 0	
	for el in data:
		if (loopIdx>=minStep):
			#numpy.min(el[1]), numpy.max(el[1]), numpy.min(el[0]), numpy.max(el[0])
			plt.imshow(el[2],origin='lower',extent=[0., 10., 0., 10.],interpolation='None')
			plt.xlabel(r'$z$', fontsize=40)
			plt.ylabel(r'$\lambda_{R.F.}$', fontsize=40)
        		plt.show()
		
			for i in range(0,el[2][0].size):
				plt.plot(lambdaRFTemplateBinCenter__, el[2][:,i])
			plt.show()
		loopIdx += 1
	'''
	return



########################################################################
###  To get path
########################################################################
	
	
def Get_Path_To_Fits():
	'''
	Return the path to the catalogue
	'''

	if (pipeline__=="DR12"):
		if (location__ == "HOME"):
			path = "/home/helion/Documents/Thèse/Data/Fits/Spectra/spec-"
		elif (location__ == "ICLUST"):
			path = "/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_7_0Bailey/"
	elif (pipeline__=='Guy' or pipeline__=='Margala'):
		path = '/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_8_guy/spectra/'
	elif (pipeline__=="MOCK"):
		if (location__ == "ICLUST"):
			path = "/home/gpfs/manip/mnt0607/bao/MockV4/M3_0_" + str(chunckNb__) + '/' + str(simulNb__).zfill(3) + '/'
	elif (pipeline__ == 'eBOSS'):
			path = "/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_8_0/"
	
	print '  Path to specFiles = ', path	
	return path			
def Get_Path_To_Save_New_Fits(iStart=0,iEnd=-1):
	'''
	Return the path to the catalogue
	'''
	
	if (iEnd==-1):
		endString = '.fits'
	else:
		endString = '_'+str(iStart)+'_'+str(iEnd)+'.fits'

	if (location__ == "HOME"):
		return
	if (pipeline__=='DR12'):
		if (not reObs__):
			path = "/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/allDR12_test" + endString
		else:
			path = "/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_reObs" + endString
	elif (pipeline__=='Guy'):
		path = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_Guy' + endString
	elif (pipeline__=='Margala'):
		path = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_Margala' + endString
	elif (pipeline__=='eBOSS'):
		path = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/eBOSS' + endString
	elif (pipeline__=='MOCK'):
		path = "/home/gpfs/manip/mnt/bao/MockV4_Production_withOwnTemplate/M3_0_" + str(chunckNb__) + '_python/' + str(simulNb__).zfill(3) +'/Mock__DR11' + endString
		
	print '  Path where to save production fits file = ', path
	return path
def Get_Data_Fits_File():
	'''
	Return the folder and the path to the catalogue
	'''
	
	if (pipeline__=="DR12"):
		if (location__ == "HOME"):
			folder = '/home/helion/Documents/Thèse/Results/RootFile/'
			path   = folder + 'allDR12_test'
		elif (location__ == "ICLUST"):
			folder = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/'
			path   = folder + 'allDR12.fits'
	elif (pipeline__=="MOCK"):
		if (location__ == "ICLUST"):
			folder = '/home/gpfs/manip/mnt/bao/MockV4_Production_withOwnTemplate/M3_0_' + str(chunckNb__) + '_python/' + str(simulNb__).zfill(3) + '/'
			path   = folder + 'allMock.fits'
		
	print '  Path where to get data fits folder = ', folder
	print '  Path where to get data fits file   = ', path
	return [folder,path]
def Get_Catalogue():
	'''
	'''

	plate_mjd_fiber_list = []
	
	if (location__ == "HOME"):
		return
	
	if (pipeline__=="DR12" or pipeline__=='Guy' or pipeline__=='Margala'):
		path = "/home/gpfs/manip/mnt/bao/hdumasde/Lists/DR12Q_v2_10.fits"
	elif (pipeline__ == 'eBOSS'):
		path = "/home/gpfs/manip/mnt/bao/Spectra/spAll-v5_8_0.fits"
	elif (pipeline__=="MOCK"):
		path = "/home/gpfs/manip/mnt/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/List/DR11MultipleOfficial.txt"
	print path		
			
	if (pipeline__=='DR12' or pipeline__=='Guy' or pipeline__=='Margala'):
		cat = pyfits.open(path, memmap=True )[1].data
		print "  The size of the catalogue is           : " + str(cat.size)
		cat = cat[ cat["BAL_FLAG_VI"]==0]
		print "  We keep BAL_FLAG_VI == 0 , the size is : " + str(cat.size)
		cat = cat[ cat["Z_VI"]>minRedshift__]
		print "  We keep Z_VI > " + str(minRedshift__) + "  , the size is      : " + str(cat.size)
		cat = cat[ cat["Z_VI"]<=maxRedshift__]
		print "  We keep Z_VI <= " + str(maxRedshift__) + "  , the size is      : " + str(cat.size)
		

		if (dataFitsType__=='spPlate'):
			### Get the sorted list of plates
			plate_list = []
			for i in cat['PLATE']:
				if i not in plate_list:
					plate_list.append(i)
			plate_list = numpy.sort(plate_list)

			### Get the list of couple plate, mjd
			plate_mjd_list = []
			for i in plate_list:
				MJD_list = []
				for j in cat['MJD'][ cat['PLATE']==i ]:
					if j not in MJD_list:
						MJD_list.append(j)
				for j in MJD_list:
					plate_mjd_list.append( (i,j) )
			plate_mjd_list = numpy.asarray(plate_mjd_list)

			### Get a list where each line is a couple plate, mjd than all the fibers
			plate_mjd_fiber_list = []
		
			for i in plate_mjd_list:
				tmp_list2 = numpy.sort( cat['FIBERID'][ numpy.logical_and(cat['PLATE']==i[0],cat['MJD']==i[1]) ] )
				plate_mjd_fiber_list.append(numpy.concatenate([i,tmp_list2]))


		plate = cat['PLATE']
		mjd   = cat['MJD']
		fiber = cat['FIBERID']
		ra    = cat['RA']
		dec   = cat['DEC']
		z     = cat['Z_VI']
		totalObs = plate.size
		
		if (reObs__):
			plate_reObs = []
			mjd_reObs   = []
			fiber_reObs = []
			ra_reObs    = []
			dec_reObs   = []
			z_reObs     = []

			### Get re-observations
			cat   = cat[ cat['NSPEC_BOSS']>0 ]
			print '   Number of reobserved quasares = ', cat.size
			for el in cat:
				for i in range(0,el['NSPEC_BOSS']):
					plate_reObs.append(el['PLATE_DUPLICATE'][i])
					mjd_reObs.append(el['MJD_DUPLICATE'][i])
					fiber_reObs.append(el['FIBERID_DUPLICATE'][i])
					ra_reObs.append(el['RA'])
					dec_reObs.append(el['DEC'])
					z_reObs.append(el['Z_VI'])

			plate = plate_reObs
			mjd   = mjd_reObs
			fiber = fiber_reObs
			ra    = ra_reObs
			dec   = dec_reObs
			z     = z_reObs

		print '  Number of additionnals observations = ',plate.size-totalObs
		
		
	elif (pipeline__ == 'eBOSS'):
		redShiftKey = 'Z'
		cat = pyfits.open(path)[1].data
		print "  The size of the catalogue is           : ", cat.size
		cat = cat[ cat[redShiftKey]>minRedshift__]
		print "  We keep Z > " + str(minRedshift__) + "  , the size is      : ", cat.size
		cat = cat[ cat[redShiftKey]<=maxRedshift__]
		print "  We keep Z <= " + str(maxRedshift__) + "  , the size is      : ", cat.size
		cat = cat[ cat['CLASS'] == 'QSO' ]
		print '  We keep CLASS == QSO, the size is      : ', cat.size
		cat = cat[ cat['ZWARNING'] == 0 ]
		print '  We keep ZWARNING == 0, the size is      : ', cat.size
		cat = cat[ cat['Z_ERR'] > 0 ]
		print '  We keep Z_ERR > 0, the size is      : ', cat.size

		if (dataFitsType__=='spPlate'):
			### Get the sorted list of plates
			plate_list = []
			for i in cat['PLATE']:
				if i not in plate_list:
					plate_list.append(i)
			plate_list = numpy.sort(plate_list)

			### Get the list of couple plate, mjd
			plate_mjd_list = []
			for i in plate_list:
				MJD_list = []
				for j in cat['MJD'][ cat['PLATE']==i ]:
					if j not in MJD_list:
						MJD_list.append(j)
				for j in MJD_list:
					plate_mjd_list.append( (i,j) )
			plate_mjd_list = numpy.asarray(plate_mjd_list)

			### Get a list where each line is a couple plate, mjd than all the fibers
			plate_mjd_fiber_list = []
		
			for i in plate_mjd_list:
				tmp_list2 = numpy.sort( cat['FIBERID'][ numpy.logical_and(cat['PLATE']==i[0],cat['MJD']==i[1]) ] )
				plate_mjd_fiber_list.append(numpy.concatenate([i,tmp_list2]))

		plate = cat['PLATE']
		mjd   = cat['MJD']
		fiber = cat['FIBERID']
		ra    = cat['RA']
		dec   = cat['DEC']
		z     = cat['Z']

	elif (pipeline__=="MOCK"):
		
		tmp_data = numpy.loadtxt(path)
		plate = tmp_data[:,0].astype('int')
		mjd   = tmp_data[:,1].astype('int')
		fiber = tmp_data[:,2].astype('int')
		ra    = tmp_data[:,3]
		dec   = tmp_data[:,4]
		z     = tmp_data[:,5]

		print '  ', plate.size
		
		tmp_bool = numpy.logical_end( (z>minRedshift__), (z<=maxRedshift__) )
		plate = plate[tmp_bool]
		mjd   = mjd[tmp_bool]
		fiber = fiber[tmp_bool]
		ra    = ra[tmp_bool]
		dec   = dec[tmp_bool]
		z     = z[tmp_bool]

		print '  ', plate.size	
		
	data = numpy.asarray( zip(plate, mjd, fiber, ra, dec, z) )
	return data, plate_mjd_fiber_list

iStart = 0
iEnd   = 0

if (len(sys.argv)>=5):

	iStart     = int(sys.argv[1])
	iEnd       = int(sys.argv[2])
	chunckNb__ = int(sys.argv[3])
	simulNb__  = int(sys.argv[4])

#make_all_Fits(iStart,iEnd)
#Merge_Files()
Get_Only_Needed_Columns()
#Get_Template()
#Plot_Results()
#Convert_Fits_To_Ntuples()
#read_data_from_root()
