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


commandProd     = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Get_delta/bin/main.exe"
commandProd     = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Correlation/bin/main.exe"
step      = 2000
nbSpectra = 300000
stepIdx = int(sys.argv[1])

def main():

	forest = 'LYB'
	reObs = False
	eBOSS = False
	method = ''#'_method1'

	if (forest=='LYA'): alphaStart__ = 1.3
	elif (forest=='CIV' or forest=='SIIV' or forest=='LYB'): alphaStart__ = 1.

	print
	print "------ Start ------"
	print


	'''
	### Get the templates
	tmp_command = commandProd
	subprocess.call(tmp_command, shell=True)
	'''

	if (stepIdx==0):
		tmp_command = "time /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Get_delta/bin/main.exe 0 0 0 0 2 0 "
		subprocess.call(tmp_command, shell=True)
		
		tmp_command = "echo " + tmp_command
		subprocess.call(tmp_command, shell=True)
		tmp_command = "echo " + time.ctime()
		subprocess.call(tmp_command, shell=True)

	if (stepIdx==1):
		### Do the fits
		nbSpectra = 300000
		step  = 2000
		first = 0
		last  = step
		while (first <= nbSpectra):
		
			tmp_command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Get_delta/bin/main.exe 0 0 " + str(first) + " " + str(last) +" 2 1 \""
			subprocess.call(tmp_command, shell=True)
		
			tmp_command = "echo " + tmp_command
			subprocess.call(tmp_command, shell=True)
			tmp_command = "echo " + time.ctime()
			subprocess.call(tmp_command, shell=True)
		
			first = first + step
			last  = last  + step
		
			#myTools.isReadyForNewJobs(150, 430)
			time.sleep(0.2)


	
	if (stepIdx==2):

		### Put files of alpha and beta together

		if (reObs):
			folder = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_DR12_Guy/DR12_reObs/histos/'
		elif (eBOSS):
			folder = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_eBOSS_Guy/all_eBOSS_primery/histos/'
		else:
			folder = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_DR12_Guy/DR12_primery/histos'+method+'/'
		scheme = 'alphaAndBeta_'+forest+'_'
		lenScheme = len(scheme)
		print folder

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
		alphaErr = data[:,4]
		betaErr  = data[:,5]
		flag     = data[:,6].astype(int)
		nbPixel  = data[:,7].astype(int)

		for el in all_t[1:]:	
			data = numpy.loadtxt(folder+el)
			idx    = numpy.append( idx,   data[:,0].astype(int) )
			alpha  = numpy.append( alpha, data[:,1]  )
			beta   = numpy.append( beta,  data[:,2]  )
			chi    = numpy.append( chi,   data[:,3]  )
			alphaErr  = numpy.append( alphaErr, data[:,4]  )
			betaErr   = numpy.append( betaErr,  data[:,5]  )
			flag      = numpy.append( flag,   data[:,6].astype(int)  )
			nbPixel   = numpy.append( nbPixel,   data[:,7].astype(int)  )

		### change 'nan' values into negative ones for the errors
		alphaErr[ numpy.isinf(alphaErr) ] = -1.
		betaErr[ numpy.isinf(betaErr) ] = -1.
		alphaErr[ numpy.isnan(alphaErr) ] = -1.
		betaErr[ numpy.isnan(betaErr) ] = -1.


		print
		print
		print '  nb spectra         : ', alpha.size
		print '  alpha              : ', alpha
		print '  beta               : ', beta
		print '  < alpha >          : ', numpy.mean(alpha)
		print '  min(alpha)         : ', numpy.min(alpha)
		print '  max(alpha)         : ', numpy.max(alpha)
		print '  < beta >           : ', numpy.mean(beta)
		print '  min(beta)          : ', numpy.min(beta)
		print '  max(beta)          : ', numpy.max(beta)
		print
		print
		##
		print '  nb pixel <=0              : ', chi[ (nbPixel<=0.) ].size
		print '  idx nbPixel <=0           : ', idx[ (nbPixel<=0.) ]
		print '  nbPixel of (nbPixel <=0)  : ', nbPixel[ (nbPixel<=0.) ]
		print
		##
		print '  nb of flag!=0      : ', chi[ (flag!=0) ].size
		print '  idx of flag!=0     : ', idx[ (flag!=0) ]
		##
		print '  nb of chi2=inf     : ', chi[ (numpy.isinf(chi)) ].size
		print '  idx of chi2=inf    : ', idx[ (numpy.isinf(chi)) ]
		##
		print '  nb of chi2==0      : ', chi[ (chi==0.) ].size
		print '  idx of chi2==0     : ', idx[ (chi==0.) ]
		##
		print '  nb of alpha==alphaStart__     : ', chi[ (alpha==alphaStart__) ].size
		print '  idx of alpha==alphaStart__    : ', idx[ (alpha==alphaStart__) ]
		##
		print '  nb of beta==0      : ', chi[ (beta==0.) ].size
		print '  idx of beta==0     : ', idx[ (beta==0.) ]
		##
		print '  nb of alphaErr<=0.      : ', chi[ (alphaErr<=0.) ].size
		print '  idx of alphaErr<=0.     : ', idx[ (alphaErr<=0.) ]
		##
		print '  nb of betaErr<=0.      : ', chi[ (betaErr<=0.) ].size
		print '  idx of betaErr<=0.     : ', idx[ (betaErr<=0.) ]
		##
		print '  nb of alphaErr>=alpha      : ', chi[ (alphaErr>=numpy.abs(alpha)) ].size
		print '  idx of alphaErr>=alpha     : ', idx[ (alphaErr>=numpy.abs(alpha)) ]
		##
		print '  nb of betaErr>=beta      : ', chi[ (betaErr>=numpy.abs(beta)) ].size
		print '  idx of betaErr>=beta     : ', idx[ (betaErr>=numpy.abs(beta)) ]

		saveAlpha = numpy.array(alpha)

		alpha[ (alpha==1.) ] = alphaStart__
		beta[ (alpha==1.) ]  = 0.
		alpha[ (nbPixel<=0.) ] = -600.
		beta[ (nbPixel<=0.) ]  = -600.
		alpha[ (flag!=0) ] = -300.
		beta[ (flag!=0) ]  = -300.
		alpha[ (alphaErr<=0.) ] = -400.
		beta[ (alphaErr<=0.) ]  = -400.
		alpha[ (betaErr<=0.) ] = -500.
		beta[ (betaErr<=0.) ]  = -500.
		alpha[ (numpy.isinf(chi)) ] = -100.
		beta[ (numpy.isinf(chi)) ]  = -100.
		alpha[ (numpy.isnan(chi)) ] = -100.
		beta[ (numpy.isnan(chi)) ]  = -100.
		alpha[ (chi==0.) ] = -200.
		beta[ (chi==0.) ]  = -200.


		numpy.savetxt(folder+scheme+'all.txt', zip(idx,alpha,beta,chi))

		cut = numpy.logical_and( chi>0., numpy.logical_and( numpy.logical_and(  numpy.logical_and( alpha!=alphaStart__, beta!=0.),  numpy.abs(alpha)<=39.5 ), numpy.abs(beta)<=0.25 ) )
		plt.hist(alpha[cut],bins=1000,label='alpha', log=True)
		plt.show()
		#plt.hist(numpy.abs(alpha[ (chi>0.) ])/alphaErr[ (chi>0.) ],bins=1000,label='alpha')
		#plt.show()
		#plt.errorbar(chi[ (chi>0.) ],alpha[ (chi>0.) ],fmt='o')
		#plt.show()
		#plt.errorbar(chi[ (chi>0.) ],alphaErr[ (chi>0.) ],fmt='o')
		#plt.show()
		#plt.errorbar(chi[ (chi>0.) ],numpy.abs(alpha[ (chi>0.) ])/alphaErr[ (chi>0.) ],fmt='o')
		#plt.show()
		plt.hist(beta[cut],bins=1000,label='beta', log=True)
		plt.show()
		#plt.errorbar(chi[ (chi>0.) ],beta[ (chi>0.) ],fmt='o')
		#plt.show()
		#plt.errorbar(chi[ (chi>0.) ],betaErr[ (chi>0.) ],fmt='o')
		#plt.show()
		#plt.errorbar(chi[ (chi>0.) ],numpy.abs(beta[ (chi>0.) ])/betaErr[ (chi>0.) ],fmt='o')
		#plt.show()
		plt.hist(chi[cut]/nbPixel[cut],bins=1000,label='chi2', log=True)
		plt.show()
		
		
		### Put the new alpha and beta in fits file
		if (reObs):
			path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits'
		elif (eBOSS):
			path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits'
		else:
			path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_DR12_Guy/DR12_primery/DR12_primery'+method+'.fits'
		print path

		file_cat = pyfits.open(path,mode='update')
		cat      = file_cat[1].data

		
		print
		print
		print '  nb spectra         : ', cat.size
		print '  alpha diff         : ', cat['ALPHA_2']-alpha
		print '  beta diff          : ', cat['BETA_2']-beta
		print

		#cut = numpy.logical_and( numpy.logical_and( numpy.logical_and( chi>0., numpy.logical_and( numpy.logical_and(  numpy.logical_and( alpha!=alphaStart__, beta!=0.),  numpy.abs(alpha)<=39.5 ), numpy.abs(beta)<=0.25 ) ),  numpy.abs(cat['ALPHA_2'])<=39.5 ), numpy.abs(cat['BETA_2'])<=0.25 )
		cut = numpy.logical_and( chi>0., numpy.logical_and( numpy.logical_and(  numpy.logical_and( alpha!=alphaStart__, beta!=0.),  numpy.abs(alpha)<=39.5 ), numpy.abs(beta)<=0.25 ) )

		cut2 = numpy.logical_and( numpy.abs(cat['ALPHA_2']-alpha)>0.0001, cut )
		#cut2 = numpy.abs(cat['ALPHA_2']-alpha)>0.0001
		#print idx[cut2].size
		#print zip(idx[cut2])

		cut2 = numpy.logical_and( numpy.abs(cat['BETA_2']-beta)>0.0001, cut )
		cut2 = numpy.abs(cat['BETA_2']-beta)>0.0001
		#print idx[cut2].size
		#print zip(idx[cut2])

		a = (cat['ALPHA_2']-alpha)[cut]
		print numpy.mean(a), numpy.min(a), numpy.max(a)
		a = (cat['BETA_2']-beta)[cut]
		print numpy.mean(a), numpy.min(a), numpy.max(a)

		plt.hist( (cat['ALPHA_2']-alpha)[cut],bins=1000,label='alpha', log=True)
		plt.show()
		plt.hist( (cat['BETA_2']-beta)[cut],bins=1000,label='beta', log=True)
		plt.show()
		

		'''
		plt.hist( cat['ALPHA_2'][:50000]-alpha,bins=1000,label='alpha', log=True)
		plt.show()
		plt.hist( cat['BETA_2'][:50000]-beta,bins=1000,label='beta', log=True)
		plt.show()
		'''

		cat['ALPHA_2'] = alpha
		cat['BETA_2']  = beta

		
		### Flag for spectra with alpha_error > alpha
		cat['BETA_1'][ (alphaErr>=numpy.abs(saveAlpha)) ] = -600.
		### Keep the chi^{2} if not reObs
		if (not reObs):
			cat['ALPHA_1'][ (nbPixel>0.) ]  = chi[ (nbPixel>0.) ]/nbPixel[ (nbPixel>0.) ]
			cat['ALPHA_1'][ (nbPixel<=0.) ] = 0.
		

		file_cat.close()
		




	print
	print " ------ End ------"
	print

main()
