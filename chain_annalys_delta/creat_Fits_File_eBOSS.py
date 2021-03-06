# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >



import subprocess
import sys
import os


from scipy import interpolate
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
from iminuit import Minuit


import myTools
from const_delta import *


location__     = 'ICLUST'   ### "HOME" or "ICLUST"
pipeline__     = 'eBOSS'     ### "DR12" or 'Guy' or 'Margala' or "MOCK" or 'eBOSS'
reObs__        = False      ### False or True
dataFitsType__ = 'spPlate'     ### What type of FITS file will we use: 'spPlate' or 'spec'


def make_all_Fits(iStart=0,iEnd=-1):

	verbose = True

	print forest__

	tmp_string = "\n   ---- Starting ---- \n\n  iStart = " + str(iStart) + "  iEnd = " + str(iEnd) + "\n"

	tmp_command = "echo \"" + tmp_string + "\""
	subprocess.call(tmp_command, shell=True)

	### Path to the folder where all spectra are
	pathToSpec = Get_Path_To_Fits()
	
	### Get the catalogue
	data, plate_list = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Run/list_'+forest__+'.npy') #Get_Catalogue()

	### Get the flux vs. lambda_Obs
	removeSkyLines = scipy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/Calibration/calibration_flux_using_CIV_forest.txt')
	removeSkyLines = interpolate.interp1d(numpy.log10(3447.5+removeSkyLines[:,0]),removeSkyLines[:,1],bounds_error=False,fill_value=1)

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
	
	### Read and write in the FITS file
	for el in cat_tbhdu:
	
		try:
			### eBOSS
			tmpCat = pyfits.open(pathToSpec + 'spPlate-' + str(el['PLATE']) + '-' + str(el['MJD']) + '.fits', memmap=True)
		except Exception,error:
			#if (not verbose): continue
			tmp_string = pathToSpec + 'spPlate-' + str(el['PLATE']) + '-' + str(el['MJD']) + '.fits'
	
			tmp_command = "echo  \"" + "  File not found \n " + tmp_string + "\""
			subprocess.call(tmp_command, shell=True)
			tmp_command = "echo  \"  " + str(error)  + "\n\""
			subprocess.call(tmp_command, shell=True)

		nbPixels         = tmpCat[0].header['NAXIS1']
		stepPixels       = tmpCat[0].header['CD1_1']
		logLamFirstPixel = tmpCat[0].header['CRVAL1']

		### Create an objects with data
		dtype_out=[('LOGLAM','<f8'),('FLUX','<f8'),('IVAR','<f8'),('AND_MASK',numpy.int16)]
  		cat = numpy.ndarray(shape=nbPixels,dtype=dtype_out)
		cat["LOGLAM"]   = numpy.arange(0.,nbPixels)*stepPixels+logLamFirstPixel
		cat["FLUX"]     = tmpCat[0].data[el['FIBERID']-1]
		cat["IVAR"]     = tmpCat[1].data[el['FIBERID']-1]
		cat["AND_MASK"] = tmpCat[2].data[el['FIBERID']-1]
		
		### Get where to cut for lambda_RF_Norma
		tmp_logZ = numpy.log10(1.+el['Z_VI'])

		### See if there are pixels
		cutForest = numpy.logical_and( (cat["LOGLAM"]>=log10lambdaRFTemplateMin__+tmp_logZ) , (cat["LOGLAM"]<log10lambdaRFTemplateMax__+tmp_logZ) )
		if (cat[cutForest].size<=nbBinRFMin__): continue
				
		### Apply cuts for bad pixels
		############################## 

		### CCD and too many sky lines
		#cat = cat[ (numpy.logical_and( (cat["LOGLAM"]>=log10lambdaObsMin__) , (cat["LOGLAM"]<log10lambdaObsMax__) )) ]
		### Sky Lines
		for lines in skyLines__:
			cat = cat[ (numpy.logical_or( (cat["LOGLAM"]<=lines[0]) , (cat["LOGLAM"]>=lines[1]) )) ]
		cat = cat[ numpy.logical_and( numpy.logical_and( (cat["IVAR"]>0.), (cat["AND_MASK"]<bit16__)), (numpy.isfinite(cat["FLUX"])) ) ]

		### Get the devident coef to correct for sky residuals
		coef = removeSkyLines(cat["LOGLAM"])
		'''
		plt.errorbar(cat["LOGLAM"],cat['FLUX'])
		plt.errorbar(cat["LOGLAM"],cat['FLUX']/coef)
		plt.errorbar(cat["LOGLAM"],coef)
		myTools.deal_with_plot(False,False,False)
		plt.show()
		'''

		cat['FLUX'] /= coef
		cat['IVAR'] *= coef*coef
		
			
		### Get the normalisation factor
		try:
			normFactor = numpy.mean( cat["FLUX"][ numpy.logical_and( (cat["LOGLAM"]>log10lambdaRFNormaMin__+tmp_logZ), (cat["LOGLAM"]<log10lambdaRFNormaMax__+tmp_logZ) ) ] )
		except Exception,error:
			if (verbose):
				tmp_command = "echo  \"" + "  Error in: 'normFactor = numpy.mean',  z = " + str(el['Z_VI']) + ", idx = " + str(cat_tbhdu[ (cat_tbhdu['NORM_FACTOR']!=0.) ].size+iStart) + "\""
				subprocess.call(tmp_command, shell=True)
				tmp_command = "echo  \"  " + str(error)  + "\n\""
				subprocess.call(tmp_command, shell=True)
			if (not takeNormaZone__):
				continue
			else:
				if (verbose):
					tmp_command = "echo  \"  Keep it still  \n\""
					subprocess.call(tmp_command, shell=True)
				normFactor = 1.
			
		el['NORM_FACTOR'] = normFactor
		
		### Apply cuts to keep only the forest
		cat = cat[ numpy.logical_and( (cat["LOGLAM"]>=log10lambdaRFTemplateMin__+tmp_logZ) , (cat["LOGLAM"]<log10lambdaRFTemplateMax__+tmp_logZ) ) ]
			
		### Find the number of pixels
		tmp_lenCat = cat.size
			
		### Store the data
		el['LAMBDA_OBS'][:tmp_lenCat] = cat["LOGLAM"]
		el['FLUX'][:tmp_lenCat]       = cat["FLUX"]
		el['FLUX_IVAR'][:tmp_lenCat]  = cat["IVAR"]

		'''
		print tmp_lenCat
		plt.errorbar(el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.],el['FLUX'][el["FLUX_IVAR"]>0.])
		plt.errorbar(el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.],el["FLUX_IVAR"][el["FLUX_IVAR"]>0.])
		myTools.deal_with_plot(False,False,False)
		plt.show()
			
		print tmp_lenCat
		plt.errorbar(numpy.power(10., el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.])/(1.+el['Z_VI']),el['FLUX'][el["FLUX_IVAR"]>0.])
		plt.errorbar(numpy.power(10., el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.])/(1.+el['Z_VI']),el["FLUX_IVAR"][el["FLUX_IVAR"]>0.])
		print numpy.min( el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.]), numpy.amax(el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.])
		print numpy.min( numpy.power(10., el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.]))/(1.+el['Z_VI']), numpy.amax( numpy.power(10., el["LAMBDA_OBS"][el["FLUX_IVAR"]>0.]))/(1.+el['Z_VI'])
		myTools.deal_with_plot(False,False,False)
		plt.show()
		'''

	print "  Number of forest (init)                                : " + str(cat_tbhdu.size)
	
	### Set values of columns
	cat_tbhdu['LAMBDA_OBS'][ (cat_tbhdu['LAMBDA_OBS'] != 0.) ]   = numpy.power(10., cat_tbhdu['LAMBDA_OBS'][ (cat_tbhdu['LAMBDA_OBS'] != 0.) ])
	cat_tbhdu['LAMBDA_RF']                                       = numpy.dot(numpy.diag(1./(1.+cat_tbhdu['Z_VI'])),cat_tbhdu['LAMBDA_OBS'])

	### Remove forest with less than nbBinRFMin__ pixels in the forest
	cut_noForestPixel = numpy.logical_and( numpy.logical_and( (cat_tbhdu['LAMBDA_RF']>=lambdaRFMin__) , (cat_tbhdu['LAMBDA_RF']<lambdaRFMax__) ), cat_tbhdu['FLUX_IVAR']>0.).astype(int)
	len_forest = numpy.sum(cut_noForestPixel,axis=1)
	cat_tbhdu = cat_tbhdu[ (len_forest>=nbBinRFMin__) ]

	print numpy.amax(len_forest)
	print numpy.amin(len_forest)

	print "  Number of forest (len_forest>=nbBinRFMin__)            : " + str(cat_tbhdu.size)
	print

	print numpy.amax(cat_tbhdu['NORM_FACTOR'])
	print numpy.amin(cat_tbhdu['NORM_FACTOR'])
	print cat_tbhdu[ cat_tbhdu['NORM_FACTOR']==1.].size
	print 
	print cat_tbhdu.size
	cat_tbhdu = cat_tbhdu[ (cat_tbhdu['NORM_FACTOR']>=0.) ]
	print cat_tbhdu.size

	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
	tbhdu.update()

	path = Get_Path_To_Save_New_Fits(iStart,iEnd)
	tbhdu.writeto(path, clobber=True)

	print "\n\n\n"







def Merge_Files():
	

	path   = "/home/gpfs/manip/mnt/bao/hdumasde/Data/" + forest__ + "/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits"
	folder = "/home/gpfs/manip/mnt/bao/hdumasde/Data/" + forest__ + "/FitsFile_eBOSS_Guy/all_eBOSS_primery/"

	print path
	print folder

	if (pipeline__=="DR12" or pipeline__=='Guy' or pipeline__=='Margala'):
		if (reObs__):
			scheme = 'DR12_reObs_'
		if (pipeline__=='Guy'):
			scheme = 'DR12_Guy_'
		if (pipeline__=='Margala'):
			scheme = 'DR12_Guy_Margala_'
		else:
			scheme = 'allDR12_'
	elif (pipeline__ == 'eBOSS'):
		scheme = 'eBOSS_'
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
		print folder+el
		all_t_file.append( pyfits.open(folder+el,  memmap=True) )
	
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


	print numpy.amin( cat['LAMBDA_OBS'][cat['FLUX_IVAR']>0.])
	print numpy.amax( cat['LAMBDA_OBS'][cat['FLUX_IVAR']>0.])

	if (takeNormaZone__):
		cut_noForestPixel  = numpy.logical_and( numpy.logical_and( (cat['LAMBDA_RF']>=lambdaRFMin__) , (cat['LAMBDA_RF']<lambdaRFMax__) ), cat['FLUX_IVAR']>0.).astype(int)
        	len_forest         = numpy.sum(cut_noForestPixel,axis=1)
		meanFluxInForest   = numpy.sum(cut_noForestPixel*cat['FLUX'],axis=1)/len_forest

		### Find the coef between <z> vs. norma_zone
		cut_badFlux = numpy.logical_and( numpy.logical_and( cat['NORM_FACTOR']!=1.,meanFluxInForest>0.), cat['NORM_FACTOR']>0.)
		#coef = numpy.mean( cat['NORM_FACTOR'][cut_badFlux]/meanFluxInForest[cut_badFlux] )

		xxx = meanFluxInForest[cut_badFlux]
		xMax = numpy.amax(xxx)
		yyy = cat['NORM_FACTOR'][cut_badFlux]
		def chi2(a0,a1,a2):
			return numpy.sum(numpy.power( yyy-(a0*xxx+a1*xxx**2+a2*xxx**3) ,2.))
		m = Minuit(chi2,a0=1.,error_a0=1.,a1=0.,error_a1=1.,a2=0.,error_a2=1., print_level=-1, errordef=0.01) 	
		m.migrad()
		a0 = m.values['a0']
		a1 = m.values['a1']
		a2 = m.values['a2']

		cut_noNormaFlux = cat['NORM_FACTOR']==1.
		print a0,a1,a2, cat[cut_noNormaFlux].size, numpy.mean( cat['NORM_FACTOR'][cut_badFlux]/meanFluxInForest[cut_badFlux] )
		cat['NORM_FACTOR'][cut_noNormaFlux] = a0*meanFluxInForest[cut_noNormaFlux]+a1*meanFluxInForest[cut_noNormaFlux]**2+a2*meanFluxInForest[cut_noNormaFlux]**3

		plt.errorbar(meanFluxInForest[cut_badFlux],cat['NORM_FACTOR'][cut_badFlux],fmt='o')
		plt.errorbar(meanFluxInForest[cut_noNormaFlux],cat['NORM_FACTOR'][cut_noNormaFlux],fmt='o')
		plt.plot(numpy.arange(0.,xMax,0.1),a0*numpy.arange(0.,xMax,0.1)+a1*numpy.arange(0.,xMax,0.1)**2+a2*numpy.arange(0.,xMax,0.1)**3)
		plt.xlabel(r'$< meanFlux >$', fontsize=40)
		plt.ylabel(r'$NORM \, FACTOR$', fontsize=40)
		myTools.deal_with_plot(False,False,False)
		plt.show()

		'''
		cut = numpy.logical_and( cut_noNormaFlux,meanFluxInForest>10.)
		for el in cat[cut]:
			cut = el['FLUX_IVAR']>0.
			xxx = el['LAMBDA_OBS'][cut]
			plt.plot(xxx,el['FLUX'][cut])
			plt.plot(xxx,el['FLUX_IVAR'][cut])
			plt.plot( [ xxx[0],xxx[-1] ], [ el['NORM_FACTOR'],el['NORM_FACTOR'] ] )
			myTools.deal_with_plot(False,False,False)
			plt.show()
		'''

	print cat['NORM_FACTOR']
	plt.hist(cat['NORM_FACTOR'],bins=100)
	plt.show()
	print '  min z = ', numpy.min( cat['Z_VI'] )
	print '  max z = ', numpy.max( cat['Z_VI'] )


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
	plate               = pyfits.Column(name='PLATE',                 format='J', array=cat['PLATE'] )
	mjd                 = pyfits.Column(name='MJD',                   format='J', array=cat['MJD'] )
	fiber               = pyfits.Column(name='FIBERID',               format='J', array=cat['FIBERID'] )
	
	ra                  = pyfits.Column(name='RA',                    format='D', array=cat['RA'], unit='degree')
	de                  = pyfits.Column(name='DEC',                   format='D', array=cat['DEC'], unit='degree')
	zz                  = pyfits.Column(name='Z',                     format='D', array=cat['Z_VI'], unit='0 (redshift)')
	boolA               = pyfits.Column(name='BIT',                   format='I', array=numpy.zeros(sizeMax))
	meanLambdaRF        = pyfits.Column(name='MEAN_FOREST_LAMBDA_RF', format='D', array=mean_lambdaRF_in_forest,             unit='angstrom')
	alpha               = pyfits.Column(name='ALPHA',                 format='D', array=alphaStart__*numpy.ones(sizeMax), unit='0' )
	beta                = pyfits.Column(name='BETA',                  format='D', array=numpy.zeros(sizeMax),             unit='1/angstrom' )

	tmp_nbBinForest     = str(nbBinRFMax__)+'D'
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',       format=tmp_nbBinForest, array=cat['LAMBDA_OBS'], unit='angstrom' )
	normFluxForest      = pyfits.Column(name='NORM_FLUX',        format=tmp_nbBinForest, array=cat['FLUX'], unit='0')
	normFluxIvarForest  = pyfits.Column(name='NORM_FLUX_IVAR',   format=tmp_nbBinForest, array=cat['FLUX_IVAR'], unit='0')
	fluxDLA             = pyfits.Column(name='FLUX_DLA',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbBinRFMax__)), unit='0')
	template            = pyfits.Column(name='TEMPLATE',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbBinRFMax__)), unit='0')
	deltaForest         = pyfits.Column(name='DELTA',            format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='0')
	deltaWeight         = pyfits.Column(name='DELTA_WEIGHT',     format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='0')
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, ra, de, zz, boolA, meanLambdaRF, alpha, beta, lambdaForest, normFluxForest, normFluxIvarForest, fluxDLA, template, deltaForest, deltaWeight])
	cat_tbhdu = tbhdu.data
	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)


	
	## Map
	plt.plot(cat_tbhdu["RA"], cat_tbhdu["DEC"], linestyle="", marker="o")
	plt.show()
	### z distrib
	plt.hist(cat_tbhdu['Z'],bins=100)
	plt.show()
	### 
	#plt.plot(cat_tbhdu["LAMBDA_RF"][ cat_tbhdu['NORM_FLUX_IVAR'] >0.], cat_tbhdu["NORM_FLUX"][ cat_tbhdu['NORM_FLUX_IVAR'] >0.], linestyle="", marker="o")
	#plt.show()
	### 
	#plt.plot(cat_tbhdu["LAMBDA_OBS"][ cat_tbhdu['NORM_FLUX_IVAR'] >0.], cat_tbhdu["NORM_FLUX"][ cat_tbhdu['NORM_FLUX_IVAR'] >0.], linestyle="", marker="o")
	#plt.show()
	
	tbhdu.update()
	tbhdu.writeto(path, clobber=True)


	print cat_tbhdu["NORM_FLUX"][0].size





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
			path = "/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_9_1/"
	
	print '  Path to specFiles = ', path	
	return path			
def Get_Path_To_Save_New_Fits(iStart=0,iEnd=-1):
	'''
	Return the path to the catalogue
	'''

	base = '/home/gpfs/manip/mnt/bao/hdumasde/Data/' + forest__ + '/FitsFile_DR12_BAL'
	
	if (iEnd==-1):
		endString = '.fits'
	else:
		endString = '_'+str(iStart)+'_'+str(iEnd)+'.fits'



	if (location__ == "HOME"):
		return
	if (pipeline__=='DR12'):
		if (not reObs__):
			path = base + "/allDR12" + endString
		else:
			path = base + "/DR12_reObs" + endString
	elif (pipeline__=='Guy'):
		path = base + '_Guy/DR12_Guy' + endString
	elif (pipeline__=='Margala'):
		path = base + 'DR12_Margala' + endString
	elif (pipeline__=='eBOSS'):
		path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/' + forest__ + '/FitsFile_eBOSS_Guy/eBOSS' + endString
	elif (pipeline__=='MOCK'):
		path = "/home/gpfs/manip/mnt/bao/hdumasde/MockV4/M3_0_" + str(chunckNb__) + '/' + str(simulNb__).zfill(3) +'/Mock__DR11' + endString
		
	print '  Path where to save production fits file = ', path
	return path
def Get_Catalogue():
	'''
	'''

	
	if (pipeline__ == 'eBOSS'):
		redShiftKey = 'Z'
		cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_9_1.fits')[1].data

		#print cat.size
		#cat = cat[ (cat['MJD']>56837) ]
		#print cat.size

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

		cat = cat[ cat['THING_ID'] > 0 ]
		print '  We keep THING_ID > 0, the size is      : ', cat.size

		cat = cat[ cat['SPECPRIMARY'] ==1 ]
		print '  We keep SPECPRIMARY == 1, the size is      : ', cat.size


		### BOSS_TARGET1
		selection_BOSS_TARGET1 = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
		bits_BOSS_TARGET1 = 0
		for el in selection_BOSS_TARGET1:
			bits_BOSS_TARGET1 += 2**el
		### ANCILLARY_TARGET1
		selection_ANCILLARY_TARGET1 = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,50,51,52,53,54,55,58,59]
		bits_ANCILLARY_TARGET1 = 0
		for el in selection_ANCILLARY_TARGET1:
			bits_ANCILLARY_TARGET1 += 2**el
		### ANCILLARY_TARGET2
		selection_ANCILLARY_TARGET2 = [0,1,2,3,4,5,7,8,9,10,13,14,15,24,25,26,27,31,32,33,53,54,55,56]
		bits_ANCILLARY_TARGET2 = 0
		for el in selection_ANCILLARY_TARGET2:
			bits_ANCILLARY_TARGET2 += 2**el
		### EBOSS_TARGET0
		selection_EBOSS_TARGET0 = [10,11,12,13,14,15,16,17,18,20,22,30,31,33,34,35,40]
		bits_EBOSS_TARGET0 = 0
		for el in selection_EBOSS_TARGET0:
			bits_EBOSS_TARGET0 += 2**el
		### EBOSS_TARGET1
		selection_EBOSS_TARGET1 = [9,10,11,12,13,14,15,16,17,18,30,31]
		bits_EBOSS_TARGET1 = 0
		for el in selection_EBOSS_TARGET1:
			bits_EBOSS_TARGET1 += 2**el
		### EBOSS_TARGET2
		selection_EBOSS_TARGET2 = [0,2,4,20,21,23,24,25,26,27,31]
		bits_EBOSS_TARGET2 = 0
		for el in selection_EBOSS_TARGET2:
			bits_EBOSS_TARGET2 += 2**el

		print cat.size
		select = ((cat['BOSS_TARGET1'] & bits_BOSS_TARGET1 )>0 )            \
			| ((cat['ANCILLARY_TARGET1'] & bits_ANCILLARY_TARGET1 )>0)  \
			| ((cat['ANCILLARY_TARGET2'] & bits_ANCILLARY_TARGET2 )>0 ) \
			| ((cat['EBOSS_TARGET0'] & bits_EBOSS_TARGET0 )>0 )         \
			| ((cat['EBOSS_TARGET1'] & bits_EBOSS_TARGET1 )>0 )         \
			| ((cat['EBOSS_TARGET2'] & bits_EBOSS_TARGET2 )>0 )
		cat = cat[select]
		print cat.size

		

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

	## Map
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.plot(cat["RA"], cat["DEC"], linestyle="", marker="o")
	plt.xlabel("Right Ascension (degree)")
	plt.ylabel("Declination (degree)")
	plt.title("BOSS DR12")
	plt.show()
	## Distribution redshift
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.hist(cat[redShiftKey], bins=50)
	plt.xlabel("Z")
	plt.ylabel("#")
	plt.title("BOSS DR12")
	plt.show()

	data = numpy.asarray( zip(plate, mjd, fiber, ra, dec, z) )
	return data, plate_mjd_fiber_list

iStart = 0
iEnd   = 0

if (len(sys.argv)>=5):

	iStart     = int(sys.argv[1])
	iEnd       = int(sys.argv[2])
	chunckNb__ = int(sys.argv[3])
	simulNb__  = int(sys.argv[4])

#cat = Get_Catalogue()
#numpy.save('list_'+forest__,cat)

#make_all_Fits(iStart,iEnd)
Merge_Files()
