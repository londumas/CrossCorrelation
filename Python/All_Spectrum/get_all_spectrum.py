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



import myTools
from const_delta import *


location__     = 'ICLUST'   ### "HOME" or "ICLUST"
pipeline__     = 'Guy'     ### "DR12" or 'Guy' or 'Margala' or "MOCK" or 'eBOSS'
reObs__        = False      ### False or True
dataFitsType__ = 'spec'     ### What type of FITS file will we use: 'spPlate' or 'spec'


def make_all_Fits(iStart=0,iEnd=-1):

	print forest__

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

	nbBinRFMax__ = 4844	
	tmp_nbBinForest     = str(nbBinRFMax__)+'D'
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',  format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Angstrom')
	lambdaRFForest      = pyfits.Column(name='LAMBDA_RF',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Angstrom')
	normFluxForest      = pyfits.Column(name='FLUX',        format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Flux' )
	normFluxIvarForest  = pyfits.Column(name='FLUX_IVAR',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbBinRFMax__)), unit='Flux' )
		
	del data
	
	tbhdu = pyfits.BinTableHDU.from_columns([plate,mjd,fiber,ra, de, zz, col_normFlux, lambdaForest, lambdaRFForest, normFluxForest, normFluxIvarForest])
	cat_tbhdu = tbhdu.data
	
	tmp_command = "echo  \"" + "\n  Starting the loop\n\n" + "\""
	subprocess.call(tmp_command, shell=True)
	
	if (dataFitsType__=='spec'):

		### Read and write in the FITS file
		for el in cat_tbhdu:
	
			try:
				### DR12
				#cat = pyfits.open(pathToSpec + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
				### DR12 Guy
				cat = pyfits.open(pathToSpec + str(el['PLATE']) + "/spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
				### DR12 Margala
				#cat = pyfits.open(pathToSpec + str(el['PLATE']) + "/corrected-spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
				### eBOSS
				#cat = pyfits.open(pathToSpec + 'spPlate-' + str(el['PLATE']) + '-' + str(el['MJD']) + '.fits', memmap=True)
				### mocks
				#cat = pyfits.open(pathToSpec + str(el['PLATE']) + "/mock-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
			except Exception,error:
				#tmp_string = pathToSpec + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
				tmp_string = pathToSpec + str(el['PLATE']) + "/spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
				#tmp_string = pathToSpec + str(el['PLATE']) + "/corrected-spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
				#tmp_string = pathToSpec + 'spPlate-' + str(el['PLATE']) + '-' + str(el['MJD']) + '.fits'
				#tmp_string = pathToSpec + str(el['PLATE']) + "/mock-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
	
				tmp_command = "echo  \"" + "  File not found \n " + tmp_string + "\""
				subprocess.call(tmp_command, shell=True)
				tmp_command = "echo  \"  " + str(error)  + "\n\""
				subprocess.call(tmp_command, shell=True)
				continue
				
			### Apply cuts for bad pixels
			############################## 
			
			### CCD and too many sky lines
			#cat = cat[ (numpy.logical_and( (cat["LOGLAM"]>=log10lambdaObsMin__) , (cat["LOGLAM"]<log10lambdaObsMax__) )) ]
			### Sky Lines
			#for lines in skyLines__:
			#	cat = cat[ (numpy.logical_or( (cat["LOGLAM"]<=lines[0]) , (cat["LOGLAM"]>=lines[1]) )) ]
			cat = cat[ numpy.logical_and( numpy.logical_and( (cat["IVAR"]>0.), (cat["AND_MASK"]<bit16__)), (numpy.isfinite(cat["FLUX"])) ) ]
		
			### Get where to cut for lambda_RF_Norma
			tmp_logZ = numpy.log10(1.+el['Z_VI'])
			
			### Get the normalisation factor
			try:
				normFactor = numpy.mean( cat["FLUX"][ numpy.logical_and( (cat["LOGLAM"]>log10lambdaRFNormaMin__+tmp_logZ), (cat["LOGLAM"]<log10lambdaRFNormaMax__+tmp_logZ) ) ] )
				el['NORM_FACTOR'] = normFactor
			except Exception,error:
				print '  mean of empty'
			

			### Apply cuts to keep only the forest
			#cat = cat[ numpy.logical_and( (cat["LOGLAM"]>=log10lambdaRFTemplateMin__+tmp_logZ) , (cat["LOGLAM"]<log10lambdaRFTemplateMax__+tmp_logZ) ) ]
			
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

	tbhdu = pyfits.BinTableHDU(data=cat_tbhdu)
	tbhdu.update()

	path = Get_Path_To_Save_New_Fits(iStart,iEnd)
	tbhdu.writeto(path, clobber=True)

	print "\n\n\n"







def Merge_Files():
	
	path   = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_al/DR12_primery.fits"
	folder = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_all/"
	print path
	print folder

	scheme = "DR12_Guy_"
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
	tbhdu = pyfits.BinTableHDU(data=cat)
	tbhdu.update()
	tbhdu.writeto(path, clobber=True)

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

	base = '/home/gpfs/manip/mnt/bao/hdumasde/Data/' + forest__ + '/FitsFile_DR12'
	
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
		path = base + '_Guy_all/DR12_Guy' + endString
	elif (pipeline__=='Margala'):
		path = base + 'DR12_Margala' + endString
	elif (pipeline__=='eBOSS'):
		path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/' + forest__ + '/FitsFile_eBOSS/eBOSS_' + endString
	elif (pipeline__=='MOCK'):
		path = "/home/gpfs/manip/mnt/bao/hdumasde/MockV4/M3_0_" + str(chunckNb__) + '/' + str(simulNb__).zfill(3) +'/Mock__DR11' + endString
		
	print '  Path where to save production fits file = ', path
	return path
def Get_Catalogue():
	'''
	'''

	plate_mjd_fiber_list = []
	
	if (location__ == "HOME"):
		return
	
	if (pipeline__=="DR12" or pipeline__=='Guy' or pipeline__=='Margala' or pipeline__=='MOCK'):
		path = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits"
	elif (pipeline__ == 'eBOSS'):
		path = "/home/gpfs/manip/mnt/bao/Spectra/spAll-v5_8_0.fits"
	#elif (pipeline__=="MOCK"):
	#	path = "/home/gpfs/manip/mnt/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/List/DR11MultipleOfficial.txt"
	print path		
			
	if (pipeline__=='DR12' or pipeline__=='Guy' or pipeline__=='Margala' or pipeline__=='MOCK'):
		cat = pyfits.open(path, memmap=True )[1].data
		print "  The size of the catalogue is           : " + str(cat.size)
		cat = cat[ cat["BAL_FLAG_VI"]==0]
		print "  We keep BAL_FLAG_VI == 0 , the size is : " + str(cat.size)
		cat = cat[ cat["Z_VI"]>0.]
		print "  We keep Z_VI > " + str(0.) + "  , the size is      : " + str(cat.size)
		

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
			plate_reObs = numpy.asarray([])
			mjd_reObs   = numpy.asarray([])
			fiber_reObs = numpy.asarray([])
			ra_reObs    = numpy.asarray([])
			dec_reObs   = numpy.asarray([])
			z_reObs     = numpy.asarray([])

			### Get re-observations
			cat   = cat[ cat['NSPEC_BOSS']>0 ]
			print '   Number of reobserved quasars = ', cat.size

			for el in cat:
				plate_reObs = numpy.append( plate_reObs, el['PLATE_DUPLICATE'][:el['NSPEC_BOSS']])
				mjd_reObs   = numpy.append( mjd_reObs,   el['MJD_DUPLICATE'][:el['NSPEC_BOSS']])
				fiber_reObs = numpy.append( fiber_reObs, el['FIBERID_DUPLICATE'][:el['NSPEC_BOSS']])
				ra_reObs    = numpy.append( ra_reObs,    [ el['RA']   ]*el['NSPEC_BOSS'])
				dec_reObs   = numpy.append( dec_reObs,   [ el['DEC']  ]*el['NSPEC_BOSS'])
				z_reObs     = numpy.append( z_reObs,     [ el['Z_VI'] ]*el['NSPEC_BOSS'])

				'''
				for i in range(0,el['NSPEC_BOSS']):
					plate_reObs.append(el['PLATE_DUPLICATE'][i])
					mjd_reObs.append(el['MJD_DUPLICATE'][i])
					fiber_reObs.append(el['FIBERID_DUPLICATE'][i])
					ra_reObs.append(el['RA'])
					dec_reObs.append(el['DEC'])
					z_reObs.append(el['Z_VI'])
				'''
			plate = plate_reObs.astype(int)
			mjd   = mjd_reObs.astype(int)
			fiber = fiber_reObs.astype(int)
			ra    = ra_reObs
			dec   = dec_reObs
			z     = z_reObs
			print plate
			print z

		print '  Number of additionnals observations = ',len(plate)
		
		
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
	'''
	elif (pipeline__=="MOCK"):
		
		tmp_data = numpy.loadtxt(path)
		plate = tmp_data[:,0].astype('int')
		mjd   = tmp_data[:,1].astype('int')
		fiber = tmp_data[:,2].astype('int')
		ra    = tmp_data[:,3]
		dec   = tmp_data[:,4]
		z     = tmp_data[:,5]

		print '  ', plate.size
		
		tmp_bool = numpy.logical_and( (z>minRedshift__), (z<=maxRedshift__) )
		plate = plate[tmp_bool]
		mjd   = mjd[tmp_bool]
		fiber = fiber[tmp_bool]
		ra    = ra[tmp_bool]
		dec   = dec[tmp_bool]
		z     = z[tmp_bool]

		print '  ', plate.size	
	'''
	data = numpy.asarray( zip(plate, mjd, fiber, ra, dec, z) )
	return data, plate_mjd_fiber_list

iStart = 0
iEnd   = 0

if (len(sys.argv)>=5):

	iStart     = int(sys.argv[1])
	iEnd       = int(sys.argv[2])
	chunckNb__ = int(sys.argv[3])
	simulNb__  = int(sys.argv[4])

make_all_Fits(iStart,iEnd)
#Merge_Files()
