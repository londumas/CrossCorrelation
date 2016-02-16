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


	
	### Data
	path     = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path)[1].data[:10]
	
	print cat['ALPHA_2']
	print cat['BETA_2']

	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/alphaAndBeta_0_10.txt')
	print data[:,1]
	print data[:,2]

	for i in range(0,10):
		el = cat[i]

		xxx = el['LAMBDA_RF'][ el['LAMBDA_RF']!=0 ]
		plt.errorbar(xxx,el['NORM_FLUX'][ el['LAMBDA_RF']!=0 ],fmt='o')

		yyy = (el['ALPHA_2'] + el['BETA_2']*(xxx-el['MEAN_FOREST_LAMBDA_RF']) )*el['TEMPLATE'][ el['LAMBDA_RF']!=0 ]
		plt.errorbar(xxx,yyy)

		yyy = (data[i,1] + data[i,2]*(xxx-el['MEAN_FOREST_LAMBDA_RF']) )*el['TEMPLATE'][ el['LAMBDA_RF']!=0 ]
		plt.errorbar(xxx,yyy)

		plt.show()


	print
	print " ------ End ------"
	print

def difTemplate():

	''' 
	Look at the differences between the templates

	'''

	path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/"

	data1 = numpy.loadtxt(path + 'template_0_0.txt')
	templateData1 = interpolate.interp1d(data1[:,0],data1[:,1],bounds_error=False,fill_value=0)
	plt.errorbar(data1[:,0], data1[:,1]/templateData1(1150.), fmt='o', label=r'$Simu$',color='red')

	data = numpy.loadtxt(path + 'template.txt')
	templateData = interpolate.interp1d(data[:,0],data[:,1],bounds_error=False,fill_value=0)
	plt.errorbar(data[:,0], data[:,1]/templateData(1150.), fmt='o', label=r'$Data$',color='blue')

	data3 = numpy.loadtxt(path + 'template_0_0_MocksColab.txt')
	templateData3 = interpolate.interp1d(data3[:,0],data3[:,1],bounds_error=False,fill_value=0)
	plt.errorbar(data3[:,0], data3[:,1]/templateData3(1150.), fmt='o', label=r'$Mock \, colab$',color='green')


	plt.title(r'$Template$', fontsize=40)
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$f(\lambda_{R.F.}) / f(1150.)$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	

	plt.errorbar(data[:,0], (data1[:,1]/templateData1(1150.)-data[:,1]/templateData(1150.))/(data1[:,1]/templateData1(1150.)) , fmt='o')
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$( (f(\lambda_{R.F.}) / f(1150.))_{Data} - (f(\lambda_{R.F.}) / f(1150.))_{Simu} ) / (f(\lambda_{R.F.}) / f(1150.))_{Data})$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

def lookNotFittedSpectra():


	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery_test_PDFMocksJMC.fits'

	cat = pyfits.open(path, memmap=True)[1].data
	print '  cat = ', cat.size
	sizeBefore = cat.size
	idx = numpy.arange(cat.size)
	cut = numpy.logical_or( numpy.logical_or(  numpy.logical_or( cat['ALPHA_2']==1., cat['BETA_2']==0.),  numpy.abs(cat['ALPHA_2'])>=39.5 ), numpy.abs(cat['BETA_2'])>=0.25 )
	#cut = numpy.logical_and( numpy.logical_and(  numpy.logical_and( cat['ALPHA_2']!=1., cat['BETA_2']!=0.),  numpy.abs(cat['ALPHA_2'])<=39.5 ), numpy.abs(cat['BETA_2'])<=0.25 )
	cat = cat[cut]
	idx = idx[cut]
	#cut = numpy.abs(cat['ALPHA_2']) == cat['ALPHA_2'][-1]
	#cat = cat[cut]
	#idx = idx[cut]
		
	print '  cat = ', cat.size
	print '  good = ', sizeBefore-cat.size
	print numpy.asarray(zip(cat['ALPHA_2'],cat['BETA_2'],idx))

	#weight = []
	#meanFlux = []
	#SN = []

	for i in range(0,cat.size):
		el = cat[i]
		cut = numpy.logical_and( el['NORM_FLUX_IVAR']>0.,el['DELTA_IVAR']>0.)
		print el['ALPHA_2'], el['BETA_2'], el['LAMBDA_RF'][cut].size, idx[i], el['Z_VI']
		#weight += [ numpy.mean(el['DELTA_WEIGHT'][cut]) ]
		#meanFlux += [  numpy.mean(el['NORM_FLUX'][cut]) ]
		#SN += [  numpy.mean(el['NORM_FLUX'][cut]/numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5)) ]
		#print meanFlux[i]

		print numpy.mean(el['LAMBDA_OBS'][cut])
			
		#if (numpy.mean(el['NORM_FLUX'][cut])<10.): continue
		
		plt.plot(el['LAMBDA_RF'][cut],el['NORM_FLUX'][cut],label='flux')
		#plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='flux err')
		#plt.plot(el['LAMBDA_RF'][cut],el['DELTA'][cut],label='delta')
		#plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['DELTA_IVAR'][cut],-0.5),label='delta err')
		#plt.plot(el['LAMBDA_RF'][cut],el['DELTA_WEIGHT'][cut],label='delta weight')
		#plt.plot(el['LAMBDA_RF'][cut],el['FLUX_DLA'][cut],label='DLA')

		template = (el['ALPHA_2']+el['BETA_2']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]
		plt.plot(el['LAMBDA_RF'][cut],template,label='template')
		myTools.deal_with_plot(False,False,True)
		plt.show()
			

	plt.errorbar(weight,cat['ALPHA_2'],fmt='o',label='weight')
	plt.errorbar(meanFlux,cat['ALPHA_2'],fmt='o',label='meanFlux')
	plt.errorbar(SN,cat['ALPHA_2'],fmt='o',label='SN')
	myTools.deal_with_plot(False,False,True)
	plt.ylim([ -7.,7. ])
	plt.show()
	plt.errorbar(weight,cat['BETA_2'],fmt='o',label='weight')
	plt.errorbar(meanFlux,cat['BETA_2'],fmt='o',label='meanFlux')
	plt.errorbar(SN,cat['BETA_2'],fmt='o',label='SN')
	myTools.deal_with_plot(False,False,True)
	plt.ylim([ -0.4,0.4 ])
	plt.show()

	return

def distribSomething():

	path = ['/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits']
	path = ['/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits']
	#path = ['/home/gpfs/manip/mnt0607/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits']
	cat = pyfits.open(path[0])[1].data[:10000]
	#cat2 = cat[ numpy.logical_or( numpy.logical_or(  numpy.logical_or( cat['ALPHA_2']==1., cat['BETA_2']==0.),  numpy.abs(cat['ALPHA_2'])>=39.5 ), numpy.abs(cat['BETA_2'])>=0.25 ) ]
	#cat = cat[ numpy.logical_and( numpy.logical_and(  numpy.logical_and( cat['ALPHA_2']!=1., cat['BETA_2']!=0.),  numpy.abs(cat['ALPHA_2'])<=39.5 ), numpy.abs(cat['BETA_2'])<=0.25 ) ]
	print cat.size
	#print cat2.size

	'''
	### All spectra
	value = []
	for el in cat:
		cut = numpy.logical_and( el['NORM_FLUX_IVAR']>0.,el['DELTA_IVAR']>0.)
		value += [ (el['ALPHA_2']+el['BETA_2']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut] ]
	#value = numpy.asarray(value).flatten()
	plt.hist(value,bins=100)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	'''

	'''
	value = cat['ALPHA_2']
	plt.hist(value,bins=1000)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	value = cat['BETA_2']
	plt.hist(value,bins=1000)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	cat = cat[ numpy.abs(cat['ALPHA_2'])>=3.5 ]
	print cat.size
	for i in range(0,cat.size):
		el = cat[i]
		cut = el['DELTA_WEIGHT']>0.
		plt.plot(el['LAMBDA_RF'][cut],el['NORM_FLUX'][cut],label='flux')
		plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='flux err')
		plt.plot(el['LAMBDA_RF'][cut],el['DELTA'][cut],label='delta')
		plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['DELTA_IVAR'][cut],-0.5),label='delta err')
		plt.plot(el['LAMBDA_RF'][cut],el['DELTA_WEIGHT'][cut],label='delta weight')
		plt.plot(el['LAMBDA_RF'][cut],el['FLUX_DLA'][cut],label='DLA')
		plt.plot(el['LAMBDA_RF'][cut],el['TEMPLATE'][cut],label='template')
		myTools.deal_with_plot(False,False,True)
		plt.show()
	'''

	return


def meanDelta():

	forest = 'LYA'
	alphaStart__ = 1.3
	if (forest == 'LYB'):
		lambdaRFMin__      = 800.
		lambdaRFMax__      = 1020.
	elif (forest == 'LYA'):
		lambdaRFMin__      = 1040.
		lambdaRFMax__      = 1200.
	elif (forest == 'SIIV'):
		lambdaRFMin__      = 1286.
		lambdaRFMax__      = 1380.
	elif (forest == 'CIV'):
		lambdaRFMin__      = 1410.
		lambdaRFMax__      = 1530.
	elif (forest == 'MGII'):
		lambdaRFMin__      = 1570.
		lambdaRFMax__      = 2790.

	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits'
	#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits'

	#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation_test_soved_shift_withMetals_withError_with_template/Box_000/Simu_000/Data/delta.fits'

	cat = pyfits.open(path, memmap=True)[1].data[:50000]

	cat = cat[ numpy.logical_and( numpy.logical_and(  numpy.logical_and( cat['ALPHA_2']!=alphaStart__, cat['BETA_2']!=0.),  numpy.abs(cat['ALPHA_2'])<=39.5 ), numpy.abs(cat['BETA_2'])<=0.25 ) ]
	#cat = cat[ numpy.logical_or( numpy.logical_or(  numpy.logical_or( cat['ALPHA_2']==alphaStart__, cat['BETA_2']==0.),  numpy.abs(cat['ALPHA_2'])>=39.5 ), numpy.abs(cat['BETA_2'])>=0.25 ) ]
	print cat.size
	print cat[ (cat['BETA_1']==-600.) ].size



	print '  nb |alpha| > 39.5  = ', cat[ (cat['ALPHA_2']>39.5) ].size
	print '  nb |beta|  > 0.25  = ', cat[ (cat['BETA_2']>0.25) ].size

	cat['DELTA'][ (cat['NORM_FLUX_IVAR']<0.) ] = 0.
	cat['DELTA'][ (cat['DELTA_IVAR']<0.) ] = 0.
	cat['DELTA'][ (cat['DELTA_WEIGHT']<0.) ] = 0.
	cat['DELTA'][ (cat['FLUX_DLA']<0.8) ] = 0.
	cat['DELTA'][ (cat['LAMBDA_RF']<lambdaRFMin__) ] = 0.
	cat['DELTA'][ (cat['LAMBDA_RF']>lambdaRFMax__) ] = 0.

	cat['NORM_FLUX'][ (cat['NORM_FLUX_IVAR']<0.) ] = 0.
	cat['NORM_FLUX'][ (cat['DELTA_IVAR']<0.) ] = 0.
	cat['NORM_FLUX'][ (cat['DELTA_WEIGHT']<0.) ] = 0.
	cat['NORM_FLUX'][ (cat['FLUX_DLA']<0.8) ] = 0.
	cat['NORM_FLUX'][ (cat['LAMBDA_RF']<lambdaRFMin__) ] = 0.
	cat['NORM_FLUX'][ (cat['LAMBDA_RF']>lambdaRFMax__) ] = 0.

	cat['DELTA_WEIGHT'][ (cat['NORM_FLUX_IVAR']<0.) ] = 0.
	cat['DELTA_WEIGHT'][ (cat['DELTA_IVAR']<0.) ] = 0.
	cat['DELTA_WEIGHT'][ (cat['DELTA_WEIGHT']<0.) ] = 0.
	cat['DELTA_WEIGHT'][ (cat['FLUX_DLA']<0.8) ] = 0.
	cat['DELTA_WEIGHT'][ (cat['LAMBDA_RF']<lambdaRFMin__) ] = 0.
	cat['DELTA_WEIGHT'][ (cat['LAMBDA_RF']>lambdaRFMax__) ] = 0.

	for lines in skyLines__:
		cut = numpy.logical_and( cat["LAMBDA_OBS"]>lines[0] , cat["LAMBDA_OBS"]<lines[1] )
		cat['DELTA'][ cut ]        = 0.
		cat['NORM_FLUX'][ cut ]    = 0.
		cat['DELTA_WEIGHT'][ cut ] = 0.
	
	print cat.size
	cut_noForestPixel = numpy.logical_and(cat['DELTA_WEIGHT']>0., numpy.logical_and( (cat['LAMBDA_RF']>=lambdaRFMin__) , (cat['LAMBDA_RF']<lambdaRFMax__) )).astype(int)
	len_forest = numpy.sum(cut_noForestPixel,axis=1)
	cat = cat[ (len_forest>=nbBinRFMin__) ]
	print cat.size

	### Get the mean flux
	meanFlux = numpy.average( cat['NORM_FLUX'], weights=cat['DELTA_WEIGHT'],axis=1)
	print meanFlux.size
	print meanFlux[ numpy.isfinite(meanFlux) ].size

	### Get the mean delta
	meanDelta = numpy.average( cat['DELTA'], weights=cat['DELTA_WEIGHT'],axis=1)
	print meanDelta.size
	print meanDelta[ numpy.isfinite(meanDelta) ].size

	### Get mean flux_DLA
	meanDLA = numpy.average( cat['FLUX_DLA'], weights=cat['DELTA_WEIGHT'],axis=1)
	print cat[ meanDLA!=1. ].size

	### Get the mean S/N
	cat['DELTA_WEIGHT'][ (cat['DELTA_WEIGHT']>0.) ] = 1.
	meanSNR = numpy.average( cat['NORM_FLUX']*numpy.sqrt(cat['NORM_FLUX_IVAR']), weights=cat['DELTA_WEIGHT'],axis=1)

	'''
	### \alpha vs. <flux>
	xxx = meanFlux
	yyy = cat['ALPHA_2']
	def chi2(a0,a1,a2):
		return numpy.sum(numpy.power( yyy-(a0*xxx) ,2.))
	m = Minuit(chi2,a0=1.,error_a0=1., print_level=-1, errordef=0.01) 	
	m.migrad()
	a0 = m.values['a0']
	print a0, numpy.mean(meanFlux), numpy.mean(cat['ALPHA_2'])
	plt.errorbar(meanFlux, cat['ALPHA_2'],fmt='o')
	plt.errorbar(meanFlux, a0*meanFlux,fmt='o')
	plt.xlabel(r'$<flux>$', fontsize=40)
	plt.ylabel(r'$\alpha$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	'''
	### < \delta >
	plt.hist(meanDelta,bins=1000)
	plt.xlabel(r'$< \delta >$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	'''
	### \chi^{2}/NDF
	plt.hist(cat['ALPHA_1'],bins=1000)
	plt.xlabel(r'$\chi^{2}/N.D.F.$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	'''
	### <SNR> vs. <delta>
	plt.errorbar(numpy.abs(meanDelta[meanDLA==1.]),meanSNR[meanDLA==1.],fmt='o',label='noDLA')
	plt.errorbar(numpy.abs(meanDelta[meanDLA!=1.]),meanSNR[meanDLA!=1.],fmt='o',color='red',label='withDLA')
	plt.errorbar(numpy.abs(meanDelta[ (cat['BETA_1']==-600) ]),meanSNR[ (cat['BETA_1']==-600) ],fmt='o',color='green')
	plt.xlabel(r'$|<\delta>|$', fontsize=40)
	plt.ylabel(r'$<SNR>$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### <SNR> vs. |<delta>|
	plt.errorbar(numpy.abs(meanDelta[meanDLA==1.]),meanSNR[meanDLA==1.],fmt='o')
	plt.errorbar(numpy.abs(meanDelta[meanDLA!=1.]),meanSNR[meanDLA!=1.],fmt='o',color='red')
	plt.errorbar(numpy.abs(meanDelta[ (cat['BETA_1']==-600) ]),meanSNR[ (cat['BETA_1']==-600) ],fmt='o',color='green')
	plt.xlabel(r'$|<\delta>|$', fontsize=40)
	plt.ylabel(r'$<SNR>$', fontsize=40)
	myTools.deal_with_plot(True,True,True)
	plt.show()
	'''
	### |<delta>| vs. chi^2/NDF
	plt.errorbar(cat['ALPHA_1'],numpy.abs(meanDelta),fmt='o')
	plt.xlabel(r'$\chi^{2}/NDF$', fontsize=40)
	plt.ylabel(r'$|< \delta >|$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### <SNR> vs. chi^2/NDF
	plt.errorbar(cat['ALPHA_1'],meanSNR,fmt='o')
	plt.xlabel(r'$\chi^{2}/NDF$', fontsize=40)
	plt.ylabel(r'$< SNR >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### z vs. chi^2
	plt.errorbar(cat['ALPHA_1'],cat['Z_VI'],fmt='o')
	plt.xlabel(r'$\chi^{2}/NDF$', fontsize=40)
	plt.ylabel(r'$z$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### |<delta>| vs. nbPixel
	plt.errorbar(cat['NB_PIXEL'],numpy.abs(meanDelta),fmt='o')
	plt.xlabel(r'$nb \, pixel$', fontsize=40)
	plt.ylabel(r'$|< \delta >|$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### alpha vs. |<delta>|
	plt.errorbar(numpy.abs(meanDelta),cat['ALPHA_2'],fmt='o')
	plt.xlabel(r'$|< \delta >|$', fontsize=40)
	plt.ylabel(r'$cat[ALPHA_2]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### beta vs. |<delta>|
	plt.errorbar(numpy.abs(meanDelta),cat['BETA_2'],fmt='o')
	plt.xlabel(r'$|< \delta >|$', fontsize=40)
	plt.ylabel(r'$cat[BETA_2]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### alpha vs. <SNR>
	plt.errorbar(meanSNR,cat['ALPHA_2'],fmt='o')
	plt.xlabel(r'$<SNR>$', fontsize=40)
	plt.ylabel(r'$alpha$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### beta vs. <SNR>
	plt.errorbar(meanSNR,cat['BETA_2'],fmt='o')
	plt.xlabel(r'$<SNR>$', fontsize=40)
	plt.ylabel(r'$beta$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### alpha vs. Z
	plt.errorbar(cat['Z_VI'],cat['ALPHA_2'],fmt='o')
	plt.xlabel(r'$Z$', fontsize=40)
	plt.ylabel(r'$alpha$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### beta vs. Z
	plt.errorbar(cat['Z_VI'],cat['BETA_2'],fmt='o')
	plt.xlabel(r'$Z$', fontsize=40)
	plt.ylabel(r'$beta$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### |<delta>| vs. z
	plt.errorbar(cat['Z_VI'],numpy.abs(meanDelta),fmt='o')
	plt.xlabel(r'$z$', fontsize=40)
	plt.ylabel(r'$|< \delta >|$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	integralNumber = 100.*numpy.asarray( [ meanSNR[ (meanSNR<0.01*i) ].size for i in numpy.arange(0.,10000.) ] )/meanSNR.size
	
	plt.errorbar(0.01*numpy.arange(0.,10000.),integralNumber,fmt='o')
	plt.xlabel(r'$SNR$', fontsize=40)
	plt.ylabel(r'$p \, (\%)$', fontsize=40)
	myTools.deal_with_plot(True,True,True)
	plt.show()
	'''
	
	#cat = cat[ meanSNR>10. ]
	#meanDelta = meanDelta[ meanSNR>10. ]
	#meanSNR = meanSNR[ meanSNR>10. ]
	cat = cat[ numpy.abs(meanDelta)>0.5 ]
	meanDelta = meanDelta[ numpy.abs(meanDelta)>0.5 ]
	meanSNR = meanSNR[ numpy.abs(meanDelta)>0.5 ]
	#meanDelta = meanDelta[cat['Z_VI']>4 ]
	#meanSNR = meanSNR[ cat['Z_VI']>4 ]
	#cat = cat[ cat['Z_VI']>4 ]

	print '  size before each spctrum : ',  cat.size
	for i in range(0,cat.size):
		el = cat[i]
		cut = el['DELTA_WEIGHT']>0.
		template = (el['ALPHA_2']+el['BETA_2']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]

		if (template[ template<=0. ].size!=0):
			#print template[ template<=0. ]
			continue
		print meanSNR[i],meanDelta[i], el['Z_VI'], el['PLATE'], el['MJD'], el['FIBERID'], el['ALPHA_2'], el['BETA_2'], el['BETA_1']


		#myTools.plotOnSpectra_plate(el['PLATE'], el['MJD'], el['FIBERID'])
		#myTools.plotOnSpectra_spec(el['PLATE'], el['MJD'], el['FIBERID'],el['Z_VI'])

		plt.errorbar(el['LAMBDA_RF'][cut],el['NORM_FLUX'][cut],yerr=numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='$Data$')
		#plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='flux err')
		#plt.plot(el['LAMBDA_RF'][cut],el['DELTA'][cut],label='delta')
		#plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['DELTA_IVAR'][cut],-0.5),label='delta err')
		#plt.plot(el['LAMBDA_RF'][cut],el['DELTA_WEIGHT'][cut],label='delta weight')
		#plt.plot(el['LAMBDA_RF'][cut],el['FLUX_DLA'][cut],label='DLA')
		
		
		plt.plot(el['LAMBDA_RF'][cut],template,label=r'$QSO \, continuum$',color='red')
		#plt.plot(el['LAMBDA_RF'][cut],el['NORM_FLUX'][cut]/(el['DELTA'][cut]+1.),label=r'$QSO \, continuum$',color='red')

		myTools.deal_with_plot(False,False,True)
		plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
		plt.ylabel(r'$Normalized \, flux$', fontsize=40)
		plt.show()
	

	return

def plot_spectra_i_want():

	a = [(1113, 6880, 56543, 401), (1268, 7135, 56564, 740), (9054, 4371, 55830, 186), (9456, 4225, 55455, 542), (9772, 5706, 56165, 114), (10274, 4660, 56191, 732), (10712, 3736, 55214, 562), (11438, 5699, 55953, 716), (11571, 6598, 56574, 732), (11595, 5131, 55835, 964), (11817, 5126, 55923, 608), (12811, 4669, 55831, 478), (12822, 5126, 55923, 940), (14199, 5127, 55889, 259), (18604, 4271, 55507, 642), (21207, 6369, 56217, 759), (21893, 7235, 56603, 172), (22124, 6781, 56599, 338), (22421, 4386, 55540, 996), (22515, 7337, 56686, 469), (23135, 6782, 56576, 972), (26235, 4506, 55568, 524), (26253, 3753, 55486, 712), (29403, 4457, 55858, 582), (29789, 7301, 56746, 456), (30040, 7300, 56707, 958), (30228, 5151, 56567, 538), (31602, 7328, 56715, 986), (34099, 7310, 56693, 110), (36714, 5714, 56660, 534), (37090, 7294, 56739, 478), (38348, 5177, 56245, 63), (38477, 5293, 55953, 634), (39606, 3817, 55277, 480), (41143, 7304, 56745, 122), (43833, 5729, 56598, 926), (44308, 5771, 56011, 744), (48284, 5316, 55955, 424), (51273, 5806, 56310, 18), (53378, 5784, 56029, 582), (54163, 4637, 55616, 596), (55262, 7284, 56683, 74), (56907, 7384, 56715, 732), (65093, 5350, 56009, 636), (66059, 4635, 55615, 812), (66875, 7096, 56683, 828), (67066, 7096, 56683, 142), (69147, 7105, 56740, 184), (71081, 7380, 56753, 54), (75784, 6430, 56299, 282), (77339, 4849, 55945, 958), (78388, 4648, 55673, 792), (78707, 3841, 55572, 514), (81092, 7115, 56667, 554), (82314, 4654, 55659, 37), (82491, 6407, 56311, 480), (82560, 6644, 56384, 437), (87902, 4749, 55633, 764), (92656, 5982, 56074, 926), (94664, 6636, 56367, 16), (95437, 3793, 55214, 328), (97960, 7417, 56753, 718), (99968, 5989, 56312, 996), (102521, 5864, 56047, 576), (102799, 7117, 56685, 426), (104110, 6817, 56455, 559), (106206, 7408, 56780, 895), (107148, 7406, 56805, 586), (110132, 5862, 56045, 264), (112709, 6006, 56105, 356), (113375, 6629, 56365, 556), (115627, 3857, 55272, 250), (116506, 7027, 56448, 472), (117014, 3855, 55268, 546), (117936, 7339, 56739, 452), (120992, 3859, 55246, 230), (122107, 6015, 56096, 82), (124381, 5904, 56046, 500), (128042, 6046, 56096, 212), (128945, 6024, 56088, 778), (130722, 6795, 56425, 594), (132755, 4013, 55629, 846), (135990, 3959, 55679, 674), (138103, 3943, 55336, 441), (141591, 3931, 55350, 802), (142016, 6042, 56101, 50), (142259, 6043, 56096, 396), (143367, 4961, 55719, 671), (143629, 3940, 55327, 40), (145175, 4727, 55693, 293), (145606, 5203, 56034, 88), (147205, 5005, 55751, 526), (147802, 6314, 56191, 492), (151194, 4983, 55836, 531), (152393, 4984, 55827, 240), (154281, 5143, 55828, 44), (154306, 4085, 55452, 986), (154571, 4195, 55452, 377), (155657, 4375, 55889, 866), (155958, 5147, 55854, 843), (159034, 5024, 55854, 822), (164417, 4288, 55501, 924), (165147, 6590, 56273, 871), (166294, 6148, 56209, 185), (168214, 6305, 56563, 846)]


	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path, memmap=True)[1].data[:50000]

	for ell in a:
		#el = cat[ numpy.logical_and(  numpy.logical_and(cat['PLATE']==ell[1],cat['MJD']==ell[2]),cat['FIBERID']==ell[3]) ]
		#if (el.size==0): break
		#el = el[0]
		el = cat[ell[0]]

		cut = numpy.logical_and( (el['DELTA_WEIGHT']>0.), el['FLUX_DLA']>0.8 )
		for lines in skyLines__:
                                cut[ numpy.logical_and( el["LAMBDA_OBS"]>lines[0] , el["LAMBDA_OBS"]<lines[1] ) ] = False 
		template = (el['ALPHA_2']+el['BETA_2']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]

		test = el['NORM_FLUX'][cut]/(el['DELTA'][cut]*el['FLUX_DLA'][cut]+1.)

		if (test[ test<=0. ].size!=0):
			continue
		print el['ALPHA_2'], el['BETA_2'], el['Z_VI'], test[0], test[-1], numpy.mean(el['NORM_FLUX'][cut]/numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5))

		#myTools.plotOnSpectra_plate(el['PLATE'], el['MJD'], el['FIBERID'])
		#myTools.plotOnSpectra_spec(el['PLATE'], el['MJD'], el['FIBERID'],el['Z_VI'])

		plt.errorbar(el['LAMBDA_RF'][cut],el['NORM_FLUX'][cut],yerr=numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='$Data$')
		#plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='flux err')
		plt.plot(el['LAMBDA_RF'][cut],el['DELTA'][cut],label='delta')
		#plt.plot(el['LAMBDA_RF'][cut],numpy.power(el['DELTA_IVAR'][cut],-0.5),label='delta err')
		#plt.plot(el['LAMBDA_RF'][cut],el['DELTA_WEIGHT'][cut],label='delta weight')
		if (numpy.mean(el['FLUX_DLA'][cut])!=1.): plt.plot(el['LAMBDA_RF'][cut],el['FLUX_DLA'][cut],label='DLA')
		
		
		#plt.plot(el['LAMBDA_RF'][cut],template,label=r'$QSO \, continuum$',color='red')
		plt.plot(el['LAMBDA_RF'][cut],el['NORM_FLUX'][cut]/(el['DELTA'][cut]+1.),label=r'$QSO \, continuum$')

		myTools.deal_with_plot(False,False,False)
		plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
		plt.ylabel(r'$Normalized \, flux$', fontsize=40)
		plt.show()

	return

def distribPixelZ():

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path)[1].data
	lambdaRFLine__ = 1550.77845
	cut = cat['DELTA_WEIGHT']>0.
	z = cat['LAMBDA_OBS'][cut]/lambdaRFLine__-1.
	weight = cat['DELTA_WEIGHT'][cut]
	print numpy.average(z,weights=weight)
	weight /= numpy.mean(weight)
	plt.hist(z,bins=numpy.arange(0.,5.,0.05),histtype='step',label='CIV',weights=weight)


	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/SIIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path)[1].data
	lambdaRFLine__ = 1402.77291
	cut = cat['DELTA_WEIGHT']>0.
	z = cat['LAMBDA_OBS'][cut]/lambdaRFLine__-1.
	weight = cat['DELTA_WEIGHT'][cut]
	print numpy.average(z,weights=weight)
	weight /= numpy.mean(weight)
	plt.hist(z,bins=numpy.arange(0.,5.,0.05),histtype='step',label='SiIV',weights=weight)

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path)[1].data
	lambdaRFLine__ = 1215.67
	cut = cat['DELTA_WEIGHT']>0.
	z = cat['LAMBDA_OBS'][cut]/lambdaRFLine__-1.
	weight = cat['DELTA_WEIGHT'][cut]
	print numpy.average(z,weights=weight)
	weight /= numpy.mean(weight)
	plt.hist(z,bins=numpy.arange(0.,5.,0.05),histtype='step',label='LYA',weights=weight)


	myTools.deal_with_plot(False,False,True)
	plt.show()

def lookTemplate():
	'''
	'''

	#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/MGII/FitsFile_DR12_Guy/DR12_Guy_210000_211186.fits'	
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path)[1].data

	'''
	plt.hist(cat['NORM_FACTOR'],bins=1000)
	myTools.deal_with_plot(False,True,True)
	plt.show()
	'''

	cut = (cat['NORM_FLUX_IVAR'] > 0.)
	flux      = cat['NORM_FLUX'][cut]
	lambdaRF  = cat['LAMBDA_RF'][cut]
	lambdaOBS = cat['LAMBDA_OBS'][cut]
	weight    = cat['NORM_FLUX_IVAR'][cut]

	axisX, mean, err, number = myTools.Get_TProfile(lambdaRF,flux,lambdaRFTemplateBinEdges__,weight)#numpy.ones(flux.size))
	plt.errorbar(axisX,mean,yerr=err,fmt='o')
	plt.show()
	
	
	axisX, mean, err, number = myTools.Get_TProfile(lambdaOBS,flux,100,weight)#numpy.ones(flux.size))
	plt.errorbar(axisX,mean,yerr=err,fmt='o')
	plt.show()

	from scipy import interpolate
	iii = interpolate.interp1d(axisX,mean,bounds_error=False,fill_value=0)


	flux /= iii(lambdaOBS)
	cut = numpy.isfinite(flux)
	flux      = flux[cut]
	lambdaRF  = lambdaRF[cut]
	weight    = weight[cut]

	axisX, mean, err, number = myTools.Get_TProfile(lambdaRF,flux,lambdaRFTemplateBinEdges__,weight)#numpy.ones(flux.size))
	plt.errorbar(axisX,mean,yerr=err,fmt='o')
	plt.show()

	return 
def distribSpaceLambda():
	
	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Correlation/src/out.txt')
	plt.hist(data[ data<20. ],bins=2000)
	myTools.deal_with_plot(False,False,True)
	plt.show()
def correlation_norma_integral():
	'''
	'''

	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery_test_PDFMocksJMC_meanLambda.fits'
	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path)[1].data

	print '  < alpha > = ', numpy.mean(cat['ALPHA_2'])
	print '  < beta  > = ', numpy.mean(cat['BETA_2'])

	plt.hist(cat['ALPHA_2'],bins=500)
	plt.show()
	plt.hist(cat['BETA_2'],bins=500)
	plt.show()
	'''
	cat['FLUX_IVAR'][ (cat['LAMBDA_RF']<lambdaRFMin__) ] = 0.
	cat['FLUX_IVAR'][ (cat['LAMBDA_RF']>lambdaRFMax__) ] = 0.
	cat['FLUX_IVAR'][ (cat['FLUX_IVAR']<=0.) ] = 0.
	cat['FLUX_IVAR'][ (cat['FLUX_IVAR']>0.) ] = 1.
	cat['FLUX'][ (cat['FLUX_IVAR']==0.) ] = 0.

	meanFlux = numpy.average( cat['FLUX'], weights=cat['FLUX_IVAR'],axis=1)

	plt.errorbar(meanFlux,cat['NORM_FACTOR'],fmt='o')
	plt.xlabel(r'$< meanFlux >$', fontsize=40)
	plt.ylabel(r'$NORM \, FACTOR$', fontsize=40)
	myTools.deal_with_plot(True,True,True)
	plt.show()
	'''
	return
def testCosmo():
	'''

	'''

	import cosmolopy.distance as cosmology

	### Constants
	z = 2.4

	### Magneville:
	# omegaM = 0.267804
	# omegaL = 0.73
	# omegaK = 0.002117
	# h      = 0.71

	### Fiduciel BAOFIT
	# omegaM = 0.27
	# omegaL = 0.73
	# omegaK = 0
	# h      = 0.7
	# rd = 149.7

	
	cosmoCMV = {'omega_M_0':0.267804, 'omega_lambda_0':0.73, 'omega_k_0':0.002117, 'h':0.71}
	cosmoFID = {'omega_M_0':0.27, 'omega_lambda_0':0.73, 'omega_k_0':0, 'h':0.70}
	
	DAZ_CMV = cosmology.angular_diameter_distance(z,**cosmoCMV)
	HZ_CMV  = cosmology.hubble_distance_z(z,**cosmoCMV) #cosmology.hubble_z(z,**cosmoCMV)*3.085677581e22/1.e3

	DAZ_FID = cosmology.angular_diameter_distance(z,**cosmoFID)
	HZ_FID  = cosmology.hubble_distance_z(z,**cosmoFID) #cosmology.hubble_z(z,**cosmoFID)*3.085677581e22/1.e3

	print DAZ_CMV, DAZ_FID
	print HZ_CMV, HZ_FID

	print DAZ_FID/DAZ_CMV
	print HZ_FID/HZ_CMV

	return



#plot_spectra_i_want()
#lookNotFittedSpectra()
#distribSomething()
meanDelta()
#distribPixelZ()
#lookTemplate()
#distribSpaceLambda()

#correlation_norma_integral()
#testCosmo()


















