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
import warnings
warnings.filterwarnings("error")

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



folder__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/All_Spectrum/Iterative_spectrum/Run1/"


def make_all_Fits(iStart=0,iEnd=-1,step=0):

	tmp_string = "\n   ---- Starting ---- \n\n  iStart = " + str(iStart) + "  iEnd = " + str(iEnd) + "\n"

	tmp_command = "echo \"" + tmp_string + "\""
	subprocess.call(tmp_command, shell=True)

	### Correct the sky
	removeSkyLines = scipy.loadtxt('/home/gpfs/manip/mnt/bao/hdumasde/Data/CIV/FitsFile_DR14/DR14_primery_for_calib/histos/hDeltaVsLambdaObs_CIV.txt')
	removeSkyLines = interpolate.interp1d(numpy.log10(3447.5+removeSkyLines[:,0]),removeSkyLines[:,1],bounds_error=False,fill_value=1)

	### Path to the folder where all spectra are
	pathToSpec = Get_Path_To_Fits()
	
	data, plate_list = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/All_Spectrum/list_spect.npy')
	#data = data[(data[:,5]<0.4)]

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
	del data
	tbhdu = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, zz])
	cat_tbhdu = tbhdu.data
	

	### Mean spectrum
	binEdge = numpy.arange(0.,15001.,1.)
	xxx = numpy.arange(0.,15000.,1.)
	yyy = numpy.zeros( xxx.size )
	eyy = numpy.zeros( xxx.size )
	www = numpy.zeros( xxx.size )
	nbb = numpy.zeros( xxx.size )


	### Template
	template = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/All_Spectrum/Iterative_spectrum/Run1/template_step_'+str(step)+'.txt')
	template = interpolate.interp1d(template[:,0],template[:,1],bounds_error=False,fill_value=1)


	tmp_command = "echo  \"" + "\n  Starting the loop\n\n" + "\""
	subprocess.call(tmp_command, shell=True)
	
	### Read and write in the FITS file
	i = -1
	for el in cat_tbhdu:
		if (i%100==0): print i
		i += 1

		try:
			cat = pyfits.open(pathToSpec + str(el['PLATE']) + "/spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits", memmap=True)[1].data
		except Exception,error:
			tmp_string = pathToSpec + str(el['PLATE']) + "/spec-" + str(el['PLATE']) + "-" + str(el['MJD']) + "-" + str(el['FIBERID']).zfill(4) + ".fits"
			tmp_command = "echo  \"" + "  File not found \n " + tmp_string + "\""
			subprocess.call(tmp_command, shell=True)
			continue
			
		### Sky Lines
		cat = cat[ (numpy.logical_and( (cat["LOGLAM"]>=log10lambdaObsMin__) , (cat["LOGLAM"]<log10lambdaObsMax__) )) ]
		cat = cat[ numpy.logical_and( numpy.logical_and( (cat["IVAR"]>0.), (cat["AND_MASK"]<bit16__)), (numpy.isfinite(cat["FLUX"])) ) ]
		for lines in skyLines__:
			cat = cat[ (numpy.logical_or( (cat["LOGLAM"]<=lines[0]) , (cat["LOGLAM"]>=lines[1]) )) ]
		if (cat.size<=50):
			tmp_command = "echo  \"" + " cat.size<=50 \""
			subprocess.call(tmp_command, shell=True)
			continue

		### Get the devident coef to correct for sky residuals
		coef = removeSkyLines(cat["LOGLAM"])
		cat["FLUX"] /= coef
		cat['IVAR'] *= coef*coef

		lamObs = numpy.power(10.,cat["LOGLAM"])
		lamRF  = lamObs/(1.+el['Z_VI'])
		flux   = cat["FLUX"]
		#weigh = numpy.ones(flux.size) #cat["IVAR"]
		weigh = numpy.power(1.+lamObs/1215.67-1.,1.9)/(0.1+1./cat["IVAR"])
		model = template(lamRF)

		a0 = fit_spectrum(flux,model,weigh)
		if (a0<=0.01):
			tmp_command = "echo  \"" + " a0<=0.01 \""
			subprocess.call(tmp_command, shell=True)
			continue
		'''
		print a0
		plt.plot(lamRF,flux)
		plt.plot(lamRF,a0*model)
		plt.show()
		'''

		delta = flux/(a0*model)-1.
		number, axisX, axisY = numpy.histogram2d(lamRF, flux, (binEdge,1), weights=numpy.ones(delta.size))
		weight, axisX, axisY = numpy.histogram2d(lamRF, flux, (binEdge,1), weights=weigh)
		mean,   axisX, axisY = numpy.histogram2d(lamRF, flux, (binEdge,1), weights=weigh*delta)
		err,    axisX, axisY = numpy.histogram2d(lamRF, flux, (binEdge,1), weights=weigh*(delta**2.))

		yyy += mean[:,0]
		eyy += err[:,0]
		www += weight[:,0]
		nbb += number[:,0]

	numpy.savetxt(folder__+'/spectrum_new_technic_step_'+str(step)+'_'+str(iStart)+'_'+str(iEnd)+'.txt',zip(xxx, yyy, eyy, www, nbb))
	
	print "\n\n\n"

	return


def fit_spectrum(f,m,w):

	def chi2(a0):
		return numpy.sum( numpy.power(f-a0*m,2)*w )

	minuit = Minuit(chi2,a0=1.,error_a0=1.,print_level=-1, errordef=0.01) 	
	minuit.migrad()

	return minuit.values['a0']
def get_delta(step):

	print folder__

	raw_scheme = 'spectrum_new_technic_step_'+str(step)+'_'

	### Mean spectrum
	xxx = numpy.arange(0.,15000.,1.)
	yyy = numpy.zeros( xxx.size )
	eyy = numpy.zeros( xxx.size )
	www = numpy.zeros( xxx.size )
	nbb = numpy.zeros( xxx.size )
	ddd = numpy.zeros( xxx.size )
	
	scheme = raw_scheme
	lenScheme = len(scheme)
		
	### Get the list of files
	all_t = os.listdir(folder__)	
	tmp_all_t = []
	for el in all_t:
		if (el[:lenScheme]==scheme):
			tmp_all_t.append(el)
	all_t = tmp_all_t
	
	### Sort the list of files
	convert      = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	all_t = sorted(all_t, key = alphanum_key)
	for i in range(len(all_t)):
		print all_t[i]
		data = numpy.loadtxt(folder__+all_t[i])
		yyy += data[:,1]
		eyy += data[:,2]
		www += data[:,3]
		nbb += data[:,4]

	cut = nbb>0.
	ddd[cut] = yyy[cut]/www[cut]
	numpy.savetxt(folder__+'/residual_step_'+str(step)+'.txt',zip(xxx, ddd))

	return
def get_new_spectrum(step):

	lambdaRFNormaMin = 1275.
	lambdaRFNormaMax = 1295.

	xxx = numpy.arange(0.,15000.,1.)
	if (step==0):
		template = numpy.ones( (xxx.size,2) )
		template[:,0] = xxx
		ddd = numpy.zeros(xxx.size)
	else:
		ddd = numpy.loadtxt(folder__+'/residual_step_'+str(step-1)+'.txt')[:,1]
		
		template = numpy.loadtxt(folder__+'/template_step_'+str(step-1)+'.txt')
		template[:,1] *= ddd+1.


		### get index of first and last not empty pixels
		idx   = numpy.arange(0.,15000.,1.).astype('int')
		idx2  = idx[ (ddd!=0.) ]
		first = idx2[0]
		last  = idx2[-1]
		template[:,1][ (idx<first) ] = template[first,1]
		template[:,1][ (idx>last) ]  = template[last,1]
		### Set values inside the spectrum without pixels
		tmp_template = interpolate.interp1d(template[:,0][(template[:,1]!=1.)],template[:,1][(template[:,1]!=1.)],bounds_error=False,fill_value=1)
		template[:,1][(template[:,1]==1.)] = tmp_template(template[:,0][(template[:,1]==1.)])
		### normalize the spectrum
		norma = numpy.mean( template[:,1][ numpy.logical_and(xxx>lambdaRFNormaMin,xxx<lambdaRFNormaMax) ] )
		template[:,1] /= norma
		### No pixels less or equal to zero
		template[:,1][ (template[:,1]<=0.01) ] = 0.01
		template[:,1][ (template[:,0]<=800.) ] = 0.01

	plt.plot(xxx, template[:,1] )
	plt.show()
	plt.plot(xxx, ddd )
	plt.show()

	numpy.savetxt(folder__+'/template_step_'+str(step)+'.txt',zip(template[:,0], template[:,1]))

	return
def plot_all_spectrum():

	xxx = numpy.arange(0.,15000.,1.)
	step_nb=3


	### template
	for i in range(1,step_nb):
		template = numpy.loadtxt(folder__+'/template_step_'+str(i)+'.txt')
		plt.plot(xxx, template[:,1] )
	plt.show()

	### residuals
	for i in range(step_nb-1):
		template = numpy.loadtxt(folder__+'/residual_step_'+str(i)+'.txt')
		plt.plot(xxx, template[:,1] )
	plt.show()

	return



########################################################################
###  To get path
########################################################################
	
	
def Get_Path_To_Fits():
	return "/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_8_guy/spectra/"	
	#return "/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_10_0/"
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

iStart = 0
iEnd   = 0
step   = 0
if (len(sys.argv)>=4):

	print sys.argv[2]

	iStart     = int(sys.argv[1])
	iEnd       = int(sys.argv[2])
	step       = int(sys.argv[3])

	make_all_Fits(iStart,iEnd,step)
elif (len(sys.argv)>=2):
	step       = int(sys.argv[1])
	if (step!=0): get_delta(step-1)
	get_new_spectrum(step)
else:
	plot_all_spectrum()




