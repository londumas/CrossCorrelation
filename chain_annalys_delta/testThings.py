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
	
	print cat['ALPHA']
	print cat['BETA']

	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/alphaAndBeta_0_10.txt')
	print data[:,1]
	print data[:,2]

	for i in range(0,10):
		el = cat[i]

		xxx = el['LAMBDA_RF'][ el['LAMBDA_RF']!=0 ]
		plt.errorbar(xxx,el['NORM_FLUX'][ el['LAMBDA_RF']!=0 ],fmt='o')

		yyy = (el['ALPHA'] + el['BETA']*(xxx-el['MEAN_FOREST_LAMBDA_RF']) )*el['TEMPLATE'][ el['LAMBDA_RF']!=0 ]
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
	cut = numpy.logical_or( numpy.logical_or(  numpy.logical_or( cat['ALPHA']==1., cat['BETA']==0.),  numpy.abs(cat['ALPHA'])>=39.5 ), numpy.abs(cat['BETA'])>=0.25 )
	#cut = numpy.logical_and( numpy.logical_and(  numpy.logical_and( cat['ALPHA']!=1., cat['BETA']!=0.),  numpy.abs(cat['ALPHA'])<=39.5 ), numpy.abs(cat['BETA'])<=0.25 )
	cat = cat[cut]
	idx = idx[cut]
	#cut = numpy.abs(cat['ALPHA']) == cat['ALPHA'][-1]
	#cat = cat[cut]
	#idx = idx[cut]
		
	print '  cat = ', cat.size
	print '  good = ', sizeBefore-cat.size
	print numpy.asarray(zip(cat['ALPHA'],cat['BETA'],idx))

	#weight = []
	#meanFlux = []
	#SN = []

	for i in range(0,cat.size):
		el = cat[i]
		cut = numpy.logical_and( el['NORM_FLUX_IVAR']>0.,el['DELTA_IVAR']>0.)
		print el['ALPHA'], el['BETA'], el['LAMBDA_RF'][cut].size, idx[i], el['Z_VI']
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

		template = (el['ALPHA']+el['BETA']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]
		plt.plot(el['LAMBDA_RF'][cut],template,label='template')
		myTools.deal_with_plot(False,False,True)
		plt.show()
			

	plt.errorbar(weight,cat['ALPHA'],fmt='o',label='weight')
	plt.errorbar(meanFlux,cat['ALPHA'],fmt='o',label='meanFlux')
	plt.errorbar(SN,cat['ALPHA'],fmt='o',label='SN')
	myTools.deal_with_plot(False,False,True)
	plt.ylim([ -7.,7. ])
	plt.show()
	plt.errorbar(weight,cat['BETA'],fmt='o',label='weight')
	plt.errorbar(meanFlux,cat['BETA'],fmt='o',label='meanFlux')
	plt.errorbar(SN,cat['BETA'],fmt='o',label='SN')
	myTools.deal_with_plot(False,False,True)
	plt.ylim([ -0.4,0.4 ])
	plt.show()

	return

def distribSomething():

	path = ['/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits']
	path = ['/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits']
	#path = ['/home/gpfs/manip/mnt0607/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits']
	cat = pyfits.open(path[0])[1].data[:10000]
	#cat2 = cat[ numpy.logical_or( numpy.logical_or(  numpy.logical_or( cat['ALPHA']==1., cat['BETA']==0.),  numpy.abs(cat['ALPHA'])>=39.5 ), numpy.abs(cat['BETA'])>=0.25 ) ]
	#cat = cat[ numpy.logical_and( numpy.logical_and(  numpy.logical_and( cat['ALPHA']!=1., cat['BETA']!=0.),  numpy.abs(cat['ALPHA'])<=39.5 ), numpy.abs(cat['BETA'])<=0.25 ) ]
	print cat.size
	#print cat2.size

	'''
	### All spectra
	value = []
	for el in cat:
		cut = numpy.logical_and( el['NORM_FLUX_IVAR']>0.,el['DELTA_IVAR']>0.)
		value += [ (el['ALPHA']+el['BETA']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut] ]
	#value = numpy.asarray(value).flatten()
	plt.hist(value,bins=100)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	'''

	'''
	value = cat['ALPHA']
	plt.hist(value,bins=1000)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	value = cat['BETA']
	plt.hist(value,bins=1000)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	cat = cat[ numpy.abs(cat['ALPHA'])>=3.5 ]
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

	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/FitsFile_DR12_Guy_Margala/DR12_primery/DR12_primery.fits'
	#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest+'/DR14/DR14_primery/DR14_primery.fits'
	#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/DR12_Nicolas/delta.fits'
	#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits'
	#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits'

	#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_000/Simu_000/Data/delta.fits'

	cat = pyfits.open(path, memmap=True)[1].data

	cat = cat[ numpy.logical_and( numpy.logical_and(  numpy.logical_and( cat['ALPHA']!=alphaStart__, cat['BETA']!=0.),  numpy.abs(cat['ALPHA'])<=99.5 ), numpy.abs(cat['BETA'])<=0.55 ) ]
	#cat = cat[ numpy.logical_or( numpy.logical_or(  numpy.logical_or( cat['ALPHA']==alphaStart__, cat['BETA']==0.),  numpy.abs(cat['ALPHA'])>=99.5 ), numpy.abs(cat['BETA'])>=0.55 ) ]
	print cat.size
	print cat[ (cat['BIT']==-600.) ].size

	'''
	plt.hist(cat['Z_VI'],bins=100)
	plt.show()
	plt.plot(cat['RA'],cat['DEC'],marker='o')
	plt.show()
	'''

	print '  nb |alpha| > 39.5  = ', cat[ (cat['ALPHA']>39.5) ].size
	print '  nb |beta|  > 0.25  = ', cat[ (cat['BETA']>0.25) ].size

	cat['DELTA'][ (cat['NORM_FLUX_IVAR']<0.) ] = 0.
	cat['DELTA'][ (cat['DELTA_WEIGHT']<0.) ] = 0.
	#cat['DELTA'][ (cat['FLUX_DLA']<0.8) ] = 0.
	#cat['DELTA'][ (cat['LAMBDA_RF']<lambdaRFMin__) ] = 0.
	#cat['DELTA'][ (cat['LAMBDA_RF']>lambdaRFMax__) ] = 0.

	cat['NORM_FLUX'][ (cat['NORM_FLUX_IVAR']<0.) ] = 0.
	cat['NORM_FLUX'][ (cat['DELTA_WEIGHT']<0.) ] = 0.
	#cat['NORM_FLUX'][ (cat['FLUX_DLA']<0.8) ] = 0.
	#cat['NORM_FLUX'][ (cat['LAMBDA_RF']<lambdaRFMin__) ] = 0.
	#cat['NORM_FLUX'][ (cat['LAMBDA_RF']>lambdaRFMax__) ] = 0.

	cat['DELTA_WEIGHT'][ (cat['NORM_FLUX_IVAR']<0.) ] = 0.
	cat['DELTA_WEIGHT'][ (cat['DELTA_WEIGHT']<0.) ] = 0.
	#cat['DELTA_WEIGHT'][ (cat['FLUX_DLA']<0.8) ] = 0.
	#cat['DELTA_WEIGHT'][ (cat['LAMBDA_RF']<lambdaRFMin__) ] = 0.
	#cat['DELTA_WEIGHT'][ (cat['LAMBDA_RF']>lambdaRFMax__) ] = 0.

	'''
	for lines in skyLines__:
		cut = numpy.logical_and( cat["LAMBDA_OBS"]>lines[0] , cat["LAMBDA_OBS"]<lines[1] )
		cat['DELTA'][ cut ]        = 0.
		cat['NORM_FLUX'][ cut ]    = 0.
		cat['DELTA_WEIGHT'][ cut ] = 0.
	'''

	### get nb of pixels in forest
	#cut_noForestPixel  = (cat['NORM_FLUX_IVAR']>0.).astype(int)
	#len_forest         = numpy.sum(cut_noForestPixel,axis=1)
	#cat = cat[len_forest<30]


	print cat.size
	#cut_noForestPixel = numpy.logical_and(cat['DELTA_WEIGHT']>0., numpy.logical_and( (cat['LAMBDA_RF']>=lambdaRFMin__) , (cat['LAMBDA_RF']<lambdaRFMax__) )).astype(int)
	cut_noForestPixel = (cat['DELTA_WEIGHT']>0.).astype(int)
	len_forest = numpy.sum(cut_noForestPixel,axis=1)
	cat = cat[ (len_forest>=nbBinRFMin__) ]
	print cat.size

	### Get the mean flux
	meanFlux = numpy.average( cat['NORM_FLUX'], weights=cat['DELTA_WEIGHT'],axis=1)
	print meanFlux.size
	print meanFlux[ numpy.isfinite(meanFlux) ].size

	### Get the mean delta
	meanDelta = numpy.average( cat['DELTA'], weights=cat['DELTA_WEIGHT'],axis=1)
	print numpy.min(meanDelta), numpy.max(meanDelta)
	print meanDelta[ numpy.isfinite(meanDelta) ].size

	### Get mean flux_DLA
	meanDLA = numpy.average( cat['FLUX_DLA'], weights=cat['DELTA_WEIGHT'],axis=1)
	print cat[ meanDLA!=1. ].size

	### Get the mean S/N
	cat['DELTA_WEIGHT'][ (cat['DELTA_WEIGHT']>0.) ] = 1.
	meanSNR = numpy.average( cat['NORM_FLUX']*numpy.sqrt(cat['NORM_FLUX_IVAR']), weights=cat['DELTA_WEIGHT'],axis=1)

	cat = cat[ meanDLA==1. ]
	meanSNR = meanSNR[ meanDLA==1. ]
	meanDelta = meanDelta[ meanDLA==1. ]

	"""
	### \alpha vs. <flux>
	xxx = meanFlux
	yyy = cat['ALPHA']
	def chi2(a0,a1,a2):
		return numpy.sum(numpy.power( yyy-(a0*xxx) ,2.))
	m = Minuit(chi2,a0=1.,error_a0=1., print_level=-1, errordef=0.01) 	
	m.migrad()
	a0 = m.values['a0']
	print a0, numpy.mean(meanFlux), numpy.mean(cat['ALPHA'])
	plt.errorbar(meanFlux, cat['ALPHA'],fmt='o')
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
	
	
	### \chi^{2}/NDF
	plt.hist(cat['ALPHA'],bins=1000)
	plt.xlabel(r'$\chi^{2}/N.D.F.$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	### <SNR> vs. <delta>
	plt.errorbar(numpy.abs(meanDelta[meanDLA==1.]),meanSNR[meanDLA==1.],fmt='o',label='noDLA')
	plt.errorbar(numpy.abs(meanDelta[meanDLA!=1.]),meanSNR[meanDLA!=1.],fmt='o',color='red',label='withDLA')
	plt.errorbar(numpy.abs(meanDelta[ (cat['BIT']==-600) ]),meanSNR[ (cat['BIT']==-600) ],fmt='o',color='green')
	plt.xlabel(r'$|<\delta>|$', fontsize=40)
	plt.ylabel(r'$<SNR>$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### <SNR> vs. |<delta>|
	plt.errorbar(numpy.abs(meanDelta[meanDLA==1.]),meanSNR[meanDLA==1.],fmt='o')
	plt.errorbar(numpy.abs(meanDelta[meanDLA!=1.]),meanSNR[meanDLA!=1.],fmt='o',color='red')
	plt.errorbar(numpy.abs(meanDelta[ (cat['BIT']==-600) ]),meanSNR[ (cat['BIT']==-600) ],fmt='o',color='green')
	plt.xlabel(r'$|<\delta>|$', fontsize=40)
	plt.ylabel(r'$<SNR>$', fontsize=40)
	myTools.deal_with_plot(True,True,True)
	plt.show()
	
	
	### |<delta>| vs. chi^2/NDF
	plt.errorbar(cat['ALPHA'],numpy.abs(meanDelta),fmt='o')
	plt.xlabel(r'$\chi^{2}/NDF$', fontsize=40)
	plt.ylabel(r'$|< \delta >|$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### <SNR> vs. chi^2/NDF
	plt.errorbar(cat['ALPHA'],meanSNR,fmt='o')
	plt.xlabel(r'$\chi^{2}/NDF$', fontsize=40)
	plt.ylabel(r'$< SNR >$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### z vs. chi^2
	plt.errorbar(cat['ALPHA'],cat['Z'],fmt='o')
	plt.xlabel(r'$\chi^{2}/NDF$', fontsize=40)
	plt.ylabel(r'$z$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### |<delta>| vs. nbPixel
	#plt.errorbar(cat['NB_PIXEL'],numpy.abs(meanDelta),fmt='o')
	#plt.xlabel(r'$nb \, pixel$', fontsize=40)
	#plt.ylabel(r'$|< \delta >|$', fontsize=40)
	#myTools.deal_with_plot(False,False,True)
	#plt.show()
	### alpha vs. |<delta>|
	plt.errorbar(numpy.abs(meanDelta),cat['ALPHA'],fmt='o')
	plt.xlabel(r'$|< \delta >|$', fontsize=40)
	plt.ylabel(r'$cat[ALPHA]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### beta vs. |<delta>|
	plt.errorbar(numpy.abs(meanDelta),cat['BETA'],fmt='o')
	plt.xlabel(r'$|< \delta >|$', fontsize=40)
	plt.ylabel(r'$cat[BETA]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### alpha vs. <SNR>
	plt.errorbar(meanSNR,cat['ALPHA'],fmt='o')
	plt.xlabel(r'$<SNR>$', fontsize=40)
	plt.ylabel(r'$alpha$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### beta vs. <SNR>
	plt.errorbar(meanSNR,cat['BETA'],fmt='o')
	plt.xlabel(r'$<SNR>$', fontsize=40)
	plt.ylabel(r'$beta$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### alpha vs. Z
	plt.errorbar(cat['Z'],cat['ALPHA'],fmt='o')
	plt.xlabel(r'$Z$', fontsize=40)
	plt.ylabel(r'$alpha$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### beta vs. Z
	plt.errorbar(cat['Z'],cat['BETA'],fmt='o')
	plt.xlabel(r'$Z$', fontsize=40)
	plt.ylabel(r'$beta$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### |<delta>| vs. z
	plt.errorbar(cat['Z'],numpy.abs(meanDelta),fmt='o')
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
	"""

	
	cat = cat[ meanSNR>30. ]
	meanDelta = meanDelta[ meanSNR>30. ]
	meanSNR = meanSNR[ meanSNR>30. ]
	#cat = cat[ numpy.abs(meanDelta)>0.5 ]
	#meanDelta = meanDelta[ numpy.abs(meanDelta)>0.5 ]
	#meanSNR = meanSNR[ numpy.abs(meanDelta)>0.5 ]
	#meanDelta = meanDelta[cat['Z_VI']>4 ]
	#meanSNR = meanSNR[ cat['Z_VI']>4 ]
	#cat = cat[ cat['Z_VI']>4 ]

	print '  size before each spctrum : ',  cat.size
	for i in range(0,cat.size):
		el = cat[i]
		cut = el['DELTA_WEIGHT']>0.
		lamRF = el['LAMBDA_OBS']/(1.+el['Z'])
		cut = numpy.logical_and( numpy.logical_and( el['DELTA_WEIGHT']>0., lamRF>=lambdaRFMin__ ), lamRF<=lambdaRFMax__)
		template = (el['ALPHA']+el['BETA']*(lamRF[cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]

		#if (template[ template<=0. ].size!=0):
			#print template[ template<=0. ]
		#	continue

		
		#if (template[ numpy.logical_and( lamRF[cut]>=lambdaRFMin__,lamRF[cut]<lambdaRFMax__) ].size<50):
		#	continue
			
		#print meanSNR[i],meanDelta[i], el['Z'], el['PLATE'], el['MJD'], el['FIBERID'], el['ALPHA'], el['BETA'], el['BIT']
		print '[',el['PLATE'],',', el['MJD'],',', el['FIBERID'],']'

		#myTools.plotOnSpectra_plate(el['PLATE'], el['MJD'], el['FIBERID'])
		#myTools.plotOnSpectra_spec(el['PLATE'], el['MJD'], el['FIBERID'],el['Z'], [lamRF[cut],template*el['FLUX_DLA'][cut]])

		plt.errorbar(lamRF[cut],el['NORM_FLUX'][cut],yerr=numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label=r'$Data$', markersize=8,linewidth=2)
		#plt.plot(lamRF[cut],numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5),label='flux err')
		#plt.plot(lamRF[cut],el['DELTA'][cut],label='delta')
		#plt.plot(lamRF[cut],numpy.power(el['DELTA_IVAR'][cut],-0.5),label='delta err')
		#plt.plot(lamRF[cut],el['DELTA_WEIGHT'][cut],label='delta weight')
		#plt.plot(lamRF[cut],template*el['FLUX_DLA'][cut],label='DLA')
		
		
		plt.plot(lamRF[cut],template*el['FLUX_DLA'][cut],label=r'$Quasar \, continuum$',color='red', markersize=8,linewidth=4)
		#plt.plot(lamRF[cut],el['NORM_FLUX'][cut]/(el['DELTA'][cut]+1.),label=r'$QSO \, continuum$',color='red')

		#myTools.deal_with_plot(False,False,True)
		plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
		plt.ylabel(r'$\phi(\lambda_{R.F.})$', fontsize=40)
		#plt.legend(fontsize=40, numpoints=1,ncol=1)
		plt.xlim( [lambdaRFMin__,lambdaRFMax__] )
		#plt.ylim( [numpy.min(el['NORM_FLUX'][cut])*1.01,numpy.max(template*el['FLUX_DLA'][cut])*1.1] )
		plt.show()
		

	return

def plot_spectra_i_want():

	a = [(5814,), (8887,), (10565,), (11871,), (11901,), (13536,), (14923,), (20576,), (20794,), (22106,), (24657,), (28252,), (29292,), (30806,), (32400,), (33330,), (38477,), (39353,), (41529,), (44466,), (47785,), (52274,), (53214,), (61417,), (62002,), (64708,), (66255,), (67584,), (68941,), (69147,), (71081,), (71748,), (73578,), (75042,), (75295,), (77453,), (77482,), (77589,), (77912,), (82491,), (94857,), (95299,), (96661,), (97519,), (97960,), (99968,), (100428,), (102022,), (102521,), (104110,), (105046,), (105177,), (106206,), (106760,), (107148,), (107884,), (110086,), (110132,), (110136,), (111442,), (112250,), (124381,), (125304,), (126571,), (127453,), (127749,), (129062,), (129581,), (130722,), (130883,), (132867,), (133016,), (137418,), (137891,), (138103,), (143775,), (147362,), (151934,), (152586,), (154438,), (154571,), (156854,), (158598,), (161932,), (162795,), (164417,), (165147,), (166910,), (166965,), (168197,)]

	path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
	cat = pyfits.open(path, memmap=True)[1].data

	z = []
	l_pix = []

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
		template = (el['ALPHA']+el['BETA']*(el['LAMBDA_RF'][cut]-el['MEAN_FOREST_LAMBDA_RF']))*el['TEMPLATE'][cut]

		test = el['NORM_FLUX'][cut]/(el['DELTA'][cut]*el['FLUX_DLA'][cut]+1.)

		### Remove with DLA
		#if (numpy.mean(el['FLUX_DLA'][cut])!=1.):
		#	continue

		if (test[ test<=0. ].size!=0):
			continue
		#print el['PLATE'], el['MJD'], el['FIBERID']
		print el['ALPHA'], el['BETA'] #, el['Z_VI'], test[0], test[-1], numpy.mean(el['NORM_FLUX'][cut]/numpy.power(el['NORM_FLUX_IVAR'][cut],-0.5))
		#print test[0], test[-1]
		#print el['ALPHA'], el['BETA'], el['MEAN_FOREST_LAMBDA_RF']
		z += [el['Z_VI']]
		l_pix += [el['LAMBDA_OBS'][cut]]

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

		if (numpy.mean(el['FLUX_DLA'][cut])!=1.):
			print 'DLA'
			plt.plot(lambdaDD,el['FLUX_DLA'][cut],label='DLA')

		myTools.deal_with_plot(False,False,False)
		plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
		plt.ylabel(r'$Normalized \, flux$', fontsize=40)
		plt.show()

	plt.hist(z)
	plt.show()
	plt.hist(numpy.array(l_pix))
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

	#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/MGII/FitsFile_DR12_Guy/DR12_Guy1000011186.fits'	
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

	print '  < alpha > = ', numpy.mean(cat['ALPHA'])
	print '  < beta  > = ', numpy.mean(cat['BETA'])

	plt.hist(cat['ALPHA'],bins=500)
	plt.show()
	plt.hist(cat['BETA'],bins=500)
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


#path     = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_Margala/DR12_coAdd/DR12_coAdd.fits'
#cat = pyfits.open(path)[1].data
	
#print cat['BIT']














