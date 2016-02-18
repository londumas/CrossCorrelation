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

	#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/Box_000/Simu_000/Data/delta.fits'

	cat = pyfits.open(path, memmap=True)[1].data

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

	a = [880,1088,1268,1977,2013,2716,2744,2856,3073,3109,3627,3640,4020,5083,5567,5661,6033,6077,6371,7509,8062,9253,9813,9894,10571,11595,11768,12285,13235,14159,14199,14337,18359,18503,
19604,19724,20371,20710,21074,21134,21528,21895,22115,22136,22163,23053,23057,23196,24527,25217,25531,25556,25771,25820,26253,27008,27051,27475,27788,29135,29403,29417,29714,29975,30476,32151,
32352,32400,32799,33043,33449,34099,34107,34405,34556,34770,35099,35884,36069,36446,36814,37461,37779,38009,38516,38729,39606,39801,40526,40850,41146,41227,41387,41505,41910,42450,42576,42716,
42927,43207,43580,43666,43687,43836,43992,44887,45173,45532,46081,46324,46408,46968,47095,47361,47925,48322,48323,48634,49123,49325,50129,50202,50300,50376,50900,51710,51775,53279,53378,53428,
53536,54034,55826,56273,56460,56707,56829,56872,57710,57971,58703,58757,58905,58938,59640,60326,60540,60741,60880,61184,61937,62002,62258,62265,62720,63282,63813,64580,64664,64939,65345,65615,
66105,67152,67211,67468,67932,68044,68992,69851,70262,71029,71614,71748,71994,72063,72150,72306,73793,73960,74049,74345,75330,75458,75768,76495,76574,77111,77186,77271,77280,77446,77968,79100,
79356,79416,80768,80897,80929,81092,81118,81213,81287,81382,81703,81863,81893,82485,82633,82648,82663,82849,83379,83483,83975,84086,84143,84957,85219,85994,86200,86818,87080,87796,88304,88354,
88382,88488,88850,89788,89790,91228,91646,91771,92635,92730,92970,93943,94085,94335,94785,94888,94895,94900,95299,95464,95522,95576,96706,97020,97141,97285,97341,97343,98488,99234,99244,99542,
99569,100371,100540,101283,101925,102022,102057,102553,103684,103730,104638,105177,105600,105701,105724,106126,106335,106959,107106,107863,107913,107974,108652,108777,109345,109451,109533,109594,
109770,109876,110943,111162,111735,111855,112239,112340,112415,112709,112725,113027,113153,113275,113316,113341,113511,113732,114964,115580,115769,116639,117014,117768,118042,119400,120301,120757,
120992,121598,121808,122209,122231,123103,123312,123647,123660,123762,124252,124322,124960,125177,125504,126650,127040,127185,127746,128029,128190,128195,128583,128945,129570,129581,129604,130131,
130310,130325,130474,131036,131580,131680,132659,132867,133562,133597,133602,133826,134232,134416,134558,134790,134817,135063,135240,137591,137653,138017,138103,138493,138730,138943,139484,139638,
140194,141591,142181,142318,142678,143025,143367,143652,143746,143910,144017,144593,145006,145785,145891,146940,147037,147205,147307,147418,147504,147595,148498,148969,149261,150782,150817,150834,
150931,151123,151335,151683,152587,152686,153307,153648,153796,154235,154275,154438,154571,154717,155024,155209,155352,155657,155958,156251,156933,157865,158098,158756,158895,160213,160327,161023,161801,161986,
162494,162788,162811,164284,164314,164524,165323,165354,165419,165475,166372,166825,166885,167097,167233,168358,168589,169018,169150,169196,169295,169481,169757,170859,171286,171427]

	a = [(1268, 7135, 56564, 740), (2851, 4219, 55480, 462), (3041, 4218, 55479, 83), (4464, 7151, 56598, 266), (5497, 4221, 55443, 634), (5814, 6274, 56550, 980), (9186, 4371, 55830, 842), (11595, 5131, 55835, 964), (11817, 5126, 55923, 608), (13235, 4313, 55857, 300), (14199, 5127, 55889, 259), (19689, 4397, 55921, 837), (20926, 4395, 55828, 806), (20956, 4387, 55534, 488), (22136, 6781, 56274, 697), (26253, 3753, 55486, 712), (26329, 3677, 55205, 736), (27051, 5945, 56213, 468), (27788, 4494, 55569, 808), (28457, 7300, 56707, 261), (29403, 4457, 55858, 582), (31949, 3690, 55182, 808), (32400, 3762, 55507, 612), (33852, 4442, 55532, 960), (34099, 7310, 56693, 110), (34107, 4445, 55869, 730), (35388, 7374, 56751, 161), (38287, 7515, 56781, 574), (39606, 3817, 55277, 480), (39801, 5714, 56660, 123), (42927, 5813, 56363, 468), (43687, 3822, 55544, 838), (47095, 4640, 55927, 844), (47793, 7445, 56720, 38), (48284, 5316, 55955, 424), (50376, 4572, 55622, 698), (50534, 5785, 56269, 486), (51609, 7308, 56709, 842), (52012, 4569, 55631, 924), (53378, 5784, 56029, 582), (56829, 6660, 56370, 10), (58905, 4562, 55570, 849), (58938, 5872, 56027, 760), (59929, 3770, 55234, 624), (60305, 3832, 55289, 472), (60550, 7083, 56722, 570), (61417, 6456, 56339, 286), (62002, 3833, 55290, 512), (65542, 7088, 56657, 994), (67584, 6707, 56383, 672), (71748, 4741, 55704, 720), (72347, 6441, 56364, 382), (72835, 6649, 56364, 94), (75042, 5369, 56272, 824), (76353, 5373, 56010, 531), (77453, 3789, 55269, 202), (77589, 4648, 55673, 716), (79424, 7107, 56740, 800), (81092, 7115, 56667, 554), (81287, 5378, 56011, 2), (82314, 4654, 55659, 37), (82491, 6407, 56311, 480), (82648, 6411, 56331, 247), (82959, 7112, 56666, 780), (83483, 6644, 56384, 700), (86818, 5973, 56067, 166), (91498, 6681, 56419, 192), (96083, 6482, 56358, 103), (96706, 6970, 56444, 678), (99139, 3974, 55320, 812), (102022, 6483, 56341, 947), (103702, 6754, 56414, 632), (104110, 6817, 56455, 559), (105177, 6625, 56386, 328), (108195, 6817, 56455, 62), (108652, 3985, 55320, 660), (108995, 5999, 56089, 836), (109345, 4837, 55707, 810), (109533, 4046, 55605, 200), (109872, 6496, 56363, 346), (110569, 5437, 55973, 800), (112709, 6006, 56105, 356), (113153, 6006, 56105, 262), (114964, 5445, 55987, 891), (115627, 3857, 55272, 250), (118364, 5456, 55980, 574), (123762, 6752, 56366, 922), (126486, 6725, 56390, 184), (128945, 6024, 56088, 778), (130326, 5901, 56039, 308), (132867, 3957, 55664, 990), (132902, 4716, 55738, 538), (133399, 5487, 55982, 116), (133826, 6790, 56430, 655), (133935, 5488, 56013, 92), (138103, 3943, 55336, 441), (141096, 5010, 55748, 515), (141591, 3931, 55350, 802), (142181, 3926, 55327, 948), (142259, 6043, 56096, 396), (143367, 4961, 55719, 671), (143652, 5204, 56036, 678), (144593, 3934, 55336, 986), (145241, 6034, 56103, 942), (151335, 6319, 56452, 101), (152587, 4990, 55743, 569), (152686, 4990, 55743, 642), (153648, 4078, 55358, 8), (154275, 4194, 55450, 296), (154281, 5143, 55828, 44), (154438, 5144, 55829, 359), (154571, 4195, 55452, 377), (155296, 4196, 55478, 944), (155958, 5147, 55854, 843), (158895, 6300, 56449, 620), (161932, 4206, 55471, 283), (163949, 6592, 56535, 722), (165354, 6136, 56206, 494), (166372, 6148, 56209, 129), (166825, 6143, 56267, 832), (166885, 6157, 56238, 678), (168214, 6305, 56563, 846), (169481, 6515, 56536, 152)]

	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path, memmap=True)[1].data

	z = []

	for ell in a:
		#el = cat[ numpy.logical_and(  numpy.logical_and(cat['PLATE']==ell[1],cat['MJD']==ell[2]),cat['FIBERID']==ell[3]) ]
		#if (el.size==0): break
		#el = el[0]
		try:
			el = cat[ell[0]]
		except:
			continue

		cut = numpy.logical_and( (el['DELTA_WEIGHT']>0.), el['FLUX_DLA']>0.8 )
		for lines in skyLines__:
                                cut[ numpy.logical_and( el["LAMBDA_OBS"]>lines[0] , el["LAMBDA_OBS"]<lines[1] ) ] = False 
		template = (el['ALPHA_2']+el['BETA_2']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]

		test = el['NORM_FLUX'][cut]/(el['DELTA'][cut]*el['FLUX_DLA'][cut]+1.)

		### Remove with DLA
		#if (numpy.mean(el['FLUX_DLA'][cut])==1.):
		#	continue

		if (test[ test<=0. ].size!=0):
			continue
		print el['PLATE'], el['MJD'], el['FIBERID']
		print el['ALPHA_2'], el['BETA_2'], el['Z_VI'], test[0], test[-1], numpy.mean(el['NORM_FLUX'][cut]/numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5))
		#print test[0], test[-1]
		print el['ALPHA_2'], el['BETA_2'], el['MEAN_FOREST_LAMBDA_RF']
		z += [el['Z_VI']]

		lambdaDD = el['LAMBDA_OBS'][cut]

		#myTools.plotOnSpectra_plate(el['PLATE'], el['MJD'], el['FIBERID'])
		#myTools.plotOnSpectra_spec(el['PLATE'], el['MJD'], el['FIBERID'],el['Z_VI'])

		plt.errorbar(lambdaDD,el['NORM_FLUX'][cut],yerr=numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='$Data$')
		#plt.plot(lambdaDD,numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='flux err')
		plt.plot(lambdaDD,el['DELTA'][cut],label='delta')
		#plt.plot(lambdaDD,numpy.power(el['DELTA_IVAR'][cut],-0.5),label='delta err')
		plt.plot(lambdaDD,numpy.power(el['DELTA_WEIGHT'][cut],0.5),label='delta weight')		
		
		plt.plot(lambdaDD,template,label=r'$QSO \, continuum$',color='red')
		plt.plot(lambdaDD,el['NORM_FLUX'][cut]/(el['DELTA'][cut]*el['FLUX_DLA'][cut]+1.),label=r'$QSO \, continuum$')

		if (numpy.mean(el['FLUX_DLA'][cut])!=1.): plt.plot(lambdaDD,el['FLUX_DLA'][cut],label='DLA')

		myTools.deal_with_plot(False,False,False)
		plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
		plt.ylabel(r'$Normalized \, flux$', fontsize=40)
		plt.show()

	plt.hist(z)
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



plot_spectra_i_want()
#lookNotFittedSpectra()
#distribSomething()
#meanDelta()
#distribPixelZ()
#lookTemplate()
#distribSpaceLambda()

#correlation_norma_integral()
#testCosmo()


















