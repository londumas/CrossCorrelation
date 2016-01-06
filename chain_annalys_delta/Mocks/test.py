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
import warnings
warnings.filterwarnings("error")


import myTools
from const_delta import *

location__     = 'ICLUST'   ### "HOME" or "ICLUST"
pipeline__     = 'MOCK'     ### "DR12" or 'Guy' or 'Margala' or "MOCK" or 'eBOSS'
reObs__        = False      ### False or True
dataFitsType__ = 'spec'     ### What type of FITS file will we use: 'spPlate' or 'spec'

def main():

	RA = []
	DE = []
	ZZ = []

	print chunck, simul

	### Get the catalogue
	data, plate_list = Get_Catalogue()
	data = data.astype(int)

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_Raw_temporary/M3_0_'+str(chunck)+'/rawlite/'

	for el in data:
		try:
			cat = pyfits.open(path + str(el[0]) + '/mockrawShort-' + str(el[0]) + '-' + str(el[1]) + '-' + str(el[2]).zfill(4) + '.fits.gz', memmap=True)
		except Exception,error:
			continue

		zz   = cat[11+simul].data['Z_DLA']
		size = zz.size
		
		RA += [cat[0].header['M_RA']]*size
		DE += [cat[0].header['M_DEC']]*size
		ZZ += zz.tolist()

	RA = numpy.array(RA)
	DE = numpy.array(DE)
	ZZ = numpy.array(ZZ)
	print RA.size
	print DE.size
	print ZZ.size

	col_ra              = pyfits.Column(name='RA',  format='D', array=RA, unit='deg')
	col_de              = pyfits.Column(name='DEC', format='D', array=DE, unit='deg')
	col_zz              = pyfits.Column(name='Z',   format='D', array=ZZ )
	tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/DLA_Catalogue/new_M3_0_'+str(chunck)+'__'+str(simul)+'.fits'
	print path
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

	base = "/home/gpfs/manip/mnt/bao/hdumasde/Data/" + forest__ + "/FitsFile_DR12/"
	
	if (iEnd==-1):
		endString = '.fits'
	else:
		endString = '_'+str(iStart)+'_'+str(iEnd)+'.fits'

	if (location__ == "HOME"):
		return
	if (pipeline__=='DR12'):
		if (not reObs__):
			path = base + "allDR12" + endString
		else:
			path = base + "DR12_reObs" + endString
	elif (pipeline__=='Guy'):
		path = base + 'DR12_Guy' + endString
	elif (pipeline__=='Margala'):
		path = base + 'DR12_Margala' + endString
	elif (pipeline__=='eBOSS'):
		path = base + 'eBOSS' + endString
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

chunck = int(sys.argv[1])
simul  = int(sys.argv[2])

main()
