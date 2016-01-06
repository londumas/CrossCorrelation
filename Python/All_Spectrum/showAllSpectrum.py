# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#

import subprocess
import sys
import os
import time

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

### My tools
import myTools
from myTools import Get_TProfile
### The Constants
from const_delta import *



def main():

	print
	print "------ Start ------"
	print


	### f vs. lambda_OBS
	data = numpy.loadtxt('spectrum.txt')
	xxx  = data[:,0]
	yyy  = data[:,1]
	eyyy = data[:,2]
        plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
        myTools.deal_with_plot(False, False, False)
        plt.show()
	
	


	
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_all/DR12_Guy_0_5000.fits'
	cat = pyfits.open(path)[1].data
	print cat.size
	cat = cat



	### f vs. lambda_RF
	cut       = (cat['FLUX_IVAR']>0.)
	lambdaRF  = cat['LAMBDA_RF'][ cut ]
	flux      = cat['FLUX'][ cut ]
	weight    = numpy.ones(flux.size)
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRF,flux, 500,weight)
	plt.plot(xxx, yyy)
	myTools.deal_with_plot(False, False, False)
	plt.show()
	
	
	### f vs. lambda_OBS
	cut       = numpy.logical_and( (cat['FLUX_IVAR']>0.), (cat['LAMBDA_RF']>=lambdaRFNormaMax__) )
	lambdaOBS = cat['LAMBDA_OBS'][ cut ]
	flux      = cat['FLUX'][ cut ]
	weight    = numpy.ones(flux.size)
	min = numpy.amin(lambdaOBS)
	max = numpy.amax(lambdaOBS)
	lambdaObsBinEdges__ = numpy.arange(min, max+1., 1.)
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaOBS,flux, lambdaObsBinEdges__,weight)
	plt.plot(xxx, yyy)
	myTools.deal_with_plot(False, False, False)
	plt.show()
	


	### Get the previus histo to a function
	correction = interpolate.interp1d(xxx,yyy,bounds_error=False,fill_value=0)
	
	
	
	### f vs. lambda_OBS (in forest)
	cut       = numpy.logical_and( (cat['FLUX_IVAR']>0.), numpy.logical_and( (cat['LAMBDA_RF']>=lambdaRFMin__),(cat['LAMBDA_RF']<lambdaRFMax__))  )
	lambdaOBS = cat['LAMBDA_OBS'][ cut ]
	flux      = cat['FLUX'][ cut ]
	weight    = numpy.ones(flux.size)
	
	
	xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaOBS,flux, lambdaObsBinEdges__,weight)
	plt.plot(xxx, yyy)
	
	xxx2, yyy2, eyyy2, nyyy2 = Get_TProfile(lambdaOBS,correction(lambdaOBS), lambdaObsBinEdges__,weight)
	plt.plot(xxx2, yyy2, color='red')
	
	plt.plot(xxx, yyy/yyy2, color='green')
	
	myTools.deal_with_plot(False, False, False)
	plt.show()

	plt.plot(xxx, yyy)
	plt.show()
	plt.plot(xxx2, yyy2)
	plt.show()
	plt.plot(xxx, yyy/yyy2)
	plt.show()
	
	plt.plot(xxx, yyy)
	plt.plot(xxx, yyy/yyy2)
	plt.show()

	print " ------ End ------"
	print

main()








