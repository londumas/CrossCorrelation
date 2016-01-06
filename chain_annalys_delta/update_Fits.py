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

### My tools
import myTools
from myTools import Get_TProfile
### The Constants
from delta_const import *




def Get_Template(loopIdx):
	'''
	'''

	### Init values for the weight
	init_eta = 1.
	init_sigma2LSS = 0.1
	
	### Print that it is starting
	tmp_command = "echo  \"  " +  "  ---- Begin ---- " + "\n\""
	subprocess.call(tmp_command, shell=True)
	
	### Get Data
	### ,mode='update' ##, memmap=True
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/DR12_primery/DR12_primery_test.fits'
	print path
	file_cat = pyfits.open(path,mode='update')
	print '  Nb total of forest    = ', file_cat[1].data.size
	cat = file_cat[1].data
	print '  Nb selected of forest = ', cat.size


	### Load the histos
	path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/"
	template         = numpy.loadtxt(path + "template.txt")
	deltaVSLambdaRF  = numpy.loadtxt(path + 'deltaVSLambdaRF.txt')
	deltaVSLambdaObs = numpy.loadtxt(path + 'deltaVSLambdaObs.txt')
	eta	         = numpy.loadtxt(path + "eta.txt")
	sigma2LSS        = numpy.loadtxt(path + "sigma2LSS.txt")

	########################## TMP #############
	cat['TEMPLATE']     = numpy.interp(cat['LAMBDA_RF'], template[:,0], template[:,1])
	##cat['TEMPLATE']     = numpy.interp(cat['LAMBDA_OBS'], meanTrandObs[loopIdx][0], meanTrandObs[loopIdx][1])*numpy.interp(cat['LAMBDA_RF'], meanTransRF[loopIdx][0], meanTransRF[loopIdx][1])
	cat['DELTA']        = (cat['NORM_FLUX']/cat['TEMPLATE']-1.)/cat['FLUX_DLA']
	cat['DELTA_IVAR']   = cat['NORM_FLUX_IVAR']*cat['TEMPLATE']*cat['TEMPLATE']*cat['FLUX_DLA']*cat['FLUX_DLA']
	
	if (loopIdx==0):
		coef_weight         = numpy.power(cat['LAMBDA_OBS']/(lambdaRFLya__*onePlusZ0__), halfGama__)
		cat['DELTA_WEIGHT'] = coef_weight/(init_sigma2LSS+1./(init_eta*cat['DELTA_IVAR']))
	else:
		coef_weight         = numpy.power(cat['LAMBDA_OBS']/(lambdaRFLya__*onePlusZ0__), halfGama__)
		cat['DELTA_WEIGHT'] = coef_weight/(numpy.interp(cat['LAMBDA_OBS']/lambdaRFLya__-1., sigma2LSS[:,0],sigma2LSS[:,1])+1./(numpy.interp(cat['LAMBDA_OBS']/lambdaRFLya__-1.,eta[:,0],eta[:,1])*cat['DELTA_IVAR']))

	

	file_cat.close()

	return

Get_Template(int(sys.argv[1]))





