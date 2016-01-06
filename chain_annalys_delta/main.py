# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#

import subprocess
import sys
import os
import time

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
commandProd     = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Correlation/bin/main.exe"
step      = 2000
nbSpectra = 300000

def main():

	stepIdx = 1

	print
	print "------ Start ------"
	print


	'''
	### Get the templates
	tmp_command = commandProd
	subprocess.call(tmp_command, shell=True)
	'''

	if (stepIdx==1):
		### Do the fits
		nbSpectra = 300000
		step  = 2000
		first = 0
		last  = step
		while (first <= nbSpectra):
		
			tmp_command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Get_delta/bin/main.exe 0 0 " + str(first) + " " + str(last) +"\""
			subprocess.call(tmp_command, shell=True)
		
			tmp_command = "echo " + tmp_command
			subprocess.call(tmp_command, shell=True)
			tmp_command = "echo " + time.ctime()
			subprocess.call(tmp_command, shell=True)
		
			first = first + step
			last  = last  + step
		
			#myTools.isReadyForNewJobs(150, 430)
			time.sleep(0.5)


	
	if (stepIdx==2):
	
		forest = 'LYA'

		### Put files of alpha and beta together

		folder = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/'
		scheme = 'alphaAndBeta_'+forest+'_'
		lenScheme = len(scheme)

		tmp_command = "rm " + folder + scheme + 'all.txt'
		subprocess.call(tmp_command, shell=True)

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

		for el in all_t:
			print el

		data  = numpy.loadtxt(folder+all_t[0])
		idx   = data[:,0].astype(int)
		alpha = data[:,1]
		beta  = data[:,2]
		chi   = data[:,3]
	
		for el in all_t[1:]:	
			data = numpy.loadtxt(folder+el)
			idx    = numpy.append( idx,   data[:,0].astype(int) )
			alpha  = numpy.append( alpha, data[:,1]  )
			beta   = numpy.append( beta,  data[:,2]  )
			chi    = numpy.append( chi,   data[:,3]  )

		print '  nb spectra         : ', alpha.size
		print '  alpha              : ', alpha
		print '  beta               : ', beta
		print '  nb of chi2=inf     : ', chi[ (numpy.isinf(chi)) ].size
		print '  idx of chi2=inf    : ', idx[ (numpy.isinf(chi)) ]
		print '  nb of alpha==1     : ', chi[ (alpha==1.) ].size
		print '  idx of alpha==1    : ', idx[ (alpha==1.) ]
		print '  nb of beta==0      : ', chi[ (beta==0.) ].size
		print '  idx of beta==0     : ', idx[ (beta==0.) ]


		numpy.savetxt(folder+scheme+'all.txt', zip(idx,alpha,beta,chi))
		
		plt.hist(alpha[ (numpy.isfinite(chi)) ],bins=1000,label='alpha')
		plt.show()
		plt.hist(beta[ (numpy.isfinite(chi)) ],bins=1000,label='beta')
		plt.show()
		plt.hist(chi[ (numpy.isfinite(chi)) ],bins=1000,label='chi2')
		plt.show()
		
		
		### Put the new alpha and beta in fits file
		path     = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_DR12_Guy/DR12_primery/DR12_primery_test_PDFMocksJMC_meanLambda.fits'
		#path     = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits'
		file_cat = pyfits.open(path,mode='update')
		cat      = file_cat[1].data

		print '  nb spectra         : ', cat.size
		print '  alpha diff         : ', cat['ALPHA_2']-alpha
		print '  beta diff          : ', cat['BETA_2']-beta
		
		plt.hist(cat['ALPHA_2']-alpha,bins=1000,label='alpha')
		plt.show()
		plt.hist(cat['BETA_2']-beta,bins=1000,label='beta')
		plt.show()
		plt.hist(chi[ (numpy.isfinite(chi)) ]/cat['NB_PIXEL'][ (numpy.isfinite(chi)) ],bins=1000,label='chi2')
		plt.show()
		
		cat['ALPHA_2'] = alpha
		cat['BETA_2']  = beta

		file_cat.close()
		




	print
	print " ------ End ------"
	print

main()
