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
from const_delta import *

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

	### Correct the sky
	removeSkyLines = scipy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/hDeltaVsLambdaObs_CIV.txt_backup')
	removeSkyLines = interpolate.interp1d(numpy.log10(3547.5+removeSkyLines[:,0]),removeSkyLines[:,1],bounds_error=False,fill_value=0)

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
	
	zz                   = pyfits.Column(name='Z_VI',                  format='D', array=data[:,5] )
	
	tmp_nbBinForest     = str(5000)+'D'
	lambdaRFForest      = pyfits.Column(name='LAMBDA_RF',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,5000)), unit='Angstrom')
        lambdOBS            = pyfits.Column(name='LAMBDA_OBS',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,5000)), unit='Angstrom')
	normFluxForest      = pyfits.Column(name='FLUX',        format=tmp_nbBinForest, array=numpy.zeros((sizeMax,5000)), unit='Flux' )
		
	del data
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, zz, normFluxForest, lambdaRFForest, lambdOBS])
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

			### Sky Lines
			#cat = cat[ (numpy.logical_and( (cat["LOGLAM"]>=log10lambdaObsMin__) , (cat["LOGLAM"]<log10lambdaObsMax__) )) ]
			cat = cat[ cat["LOGLAM"]<log10lambdaObsMax__ ]
			for lines in skyLines__:
				cat = cat[ (numpy.logical_or( (cat["LOGLAM"]<=lines[0]) , (cat["LOGLAM"]>=lines[1]) )) ]
			cat = cat[ numpy.logical_and( numpy.logical_and( (cat["IVAR"]>0.), (cat["AND_MASK"]<bit16__)), (numpy.isfinite(cat["FLUX"])) ) ]

			### Get the devident coef to correct for sky residuals
			coef = removeSkyLines(cat["LOGLAM"])
			cat["FLUX"] /= coef

			### Find the number of pixels
			tmp_lenCat = cat.size
			
			### Store the data
			el['FLUX'][:tmp_lenCat]        = cat["FLUX"]
			el['LAMBDA_RF'][:tmp_lenCat]   = numpy.power(10.,cat["LOGLAM"])/(1.+el['Z_VI'])
			el['LAMBDA_OBS'][:tmp_lenCat]  = numpy.power(10.,cat["LOGLAM"])
	
	print "  Number of forest (init)                                : " + str(cat_tbhdu.size)

	cut       =(cat_tbhdu['FLUX']!=0.)
	lambdaRF  = cat_tbhdu['LAMBDA_RF'][ cut ]
	lambdaOBS = cat_tbhdu['LAMBDA_OBS'][ cut ]
	flux      = cat_tbhdu['FLUX'][ cut ]
	weight    = numpy.ones(flux.size)

	
	### f vs. lambda_RF
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,flux, 5000,weight)	
	plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()

	numpy.savetxt('spectrum3.txt',zip(xxx, yyy, eyyy))
	
	'''
	### f vs. lambda_OBS
        xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaOBS,flux, 5000,weight)
        plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
        myTools.deal_with_plot(False, False, False)
        plt.show()
	
	numpy.savetxt('f_vs_lambdaObs.txt',zip(xxx, yyy, eyyy))
	'''
	print "\n\n\n"


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
		path = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits"
	elif (pipeline__ == 'eBOSS'):
		path = "/home/gpfs/manip/mnt/bao/Spectra/spAll-v5_8_0.fits"
	elif (pipeline__=="MOCK"):
		path = "/home/gpfs/manip/mnt/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/List/DR11MultipleOfficial.txt"
	print path		
			
	if (pipeline__=='DR12' or pipeline__=='Guy' or pipeline__=='Margala'):
		cat = pyfits.open(path)[1].data
		print "  The size of the catalogue is           : " + str(cat.size)
		cat = cat[ cat["BAL_FLAG_VI"]==0]
		print "  We keep BAL_FLAG_VI == 0 , the size is : " + str(cat.size)
		

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

def plot():

	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/All_Spectrum/spectrum.txt')

	xxx = data[:,0]
	yyy = data[:,1]
	yer = data[:,2]

	plt.errorbar(xxx, yyy, yerr=yer, marker="o")
	myTools.deal_with_plot(False, False, False)
	plt.show()

iStart = 0
iEnd   = 0

if (len(sys.argv)>=5):

	print sys.argv[2]

	iStart     = int(sys.argv[1])
	iEnd       = int(sys.argv[2])
	chunckNb__ = int(sys.argv[3])
	simulNb__  = int(sys.argv[4])

	make_all_Fits(iStart,iEnd)

plot()
