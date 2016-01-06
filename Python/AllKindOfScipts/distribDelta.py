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


commandProd     = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Get_delta/bin/main.exe"
step      = 5000
nbSpectra = 200000

def main():

	print
	print "------ Start ------"
	print

	name = ['/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits',
		'/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Data/delta.fits',
		'/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/mock.fits',
		'/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Data/delta.fits'
		]
	label = ['data','simu','pipeline','corrected simu']

	saveHist = []
	saveVar  = []

	for i in numpy.arange(len(name)):
		print '\n\n'
		print label[i]
		### Data
		cat = pyfits.open(name[i])[1].data
		print cat.size
	
		cut = numpy.logical_and(numpy.logical_and( numpy.logical_and((cat['DELTA_IVAR'] > 0.), (cat['DELTA_WEIGHT']>0.)), numpy.logical_and((cat['NORM_FLUX_IVAR'] > 0.), (cat['FLUX_DLA']>=0.8)) ), numpy.logical_and( (cat['LAMBDA_RF']>=1040.), (cat['LAMBDA_RF']<1200.)))
		delta   = cat['DELTA'][cut].flatten()
		weights = cat['DELTA_WEIGHT'][cut].flatten()


		m = numpy.average(delta,weights=weights)
		v = numpy.average((delta-m)**2, weights=weights)

		print delta.size
		print m, v

		
		### Plot
		cut = numpy.abs(delta)<10.
		delta   = delta[cut]
		weights = weights[cut]

		print '  with cut'
		m = numpy.average(delta,weights=weights)
		v = numpy.average((delta-m)**2, weights=weights)
		print delta.size
		print m, v

		delta -= m
		saveVar += [v]
		

		hist, axisX = numpy.histogram(delta,bins=numpy.arange(-10.,10.,0.1),weights=weights)
	        xxx  = numpy.array([ axisX[j]+(axisX[j+1]-axisX[j])/2. for j in range(0,axisX.size-1) ])
		hist = numpy.asarray(zip(xxx,hist))
		hist[:,1] /= numpy.sum(hist[:,1])
		#hist[:,1] /= numpy.amax(hist[:,1])
		saveHist += [hist]

	print '\n\n'


        for i in numpy.arange(len(name)):
		plt.errorbar(saveHist[i][:,0],saveHist[i][:,1],fmt='o',label=label[i])
		print numpy.sum(saveHist[i][:,1][ numpy.abs(saveHist[i][:,0])<2.*saveVar[i] ]  )

	plt.grid(True, which='both')
	plt.xlabel(r'$\delta-<\delta>$', fontsize=40)
	plt.ylabel(r'$Nb/integral$', fontsize=40)
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	plt.yscale('log')
	plt.show()
	for i in numpy.arange(len(name)):
		plt.errorbar(saveHist[i][:,0],saveHist[i][:,1],label=label[i])
	plt.grid(True, which='both')
	plt.xlabel(r'$\delta-<\delta>$', fontsize=40)
	plt.ylabel(r'$Nb/integral$', fontsize=40)
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	plt.show()

	for i in numpy.arange(len(name)):
                plt.errorbar(saveHist[i][:,0],(saveHist[i][:,1]-saveHist[0][:,1])/saveHist[0][:,1],fmt='o',label=label[i])
	plt.grid(True, which='both')
        plt.xlabel(r'$\delta-<\delta>$', fontsize=40)
        plt.ylabel(r'$(nb_{i}-nb_{data})/nb_{data}$', fontsize=40)
        plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
        plt.show()

	return

main()





























