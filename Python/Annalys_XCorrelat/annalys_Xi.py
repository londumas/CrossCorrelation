# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >



import subprocess
import sys
import os

import astropy.io.fits as pyfits
import numpy
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#import decimal
#import profile
from iminuit import Minuit
import root_numpy
import ROOT
#import array
#import re

#import warnings
#warnings.filterwarnings("error")

### My tools
import myTools
### My constants
import XCorrelat_const as const

path__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_bin10Mpc_m2_Prod_eBOSS_2/correlation_notDebug_allSky_eBOSS.root"
path__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m2_Prod_DR12_11_CoAddDelta/source.root"


def get_Data(path,draw=False):

	xi1D__   = numpy.ndarray(shape=(3,const.nbBin1D))
	xi2D__   = numpy.ndarray(shape=(3,const.nbBin2D))

	xi2D_2_distX = numpy.ndarray(shape=(const.nbBinX2D,const.nbBinY2D))
	xi2D_2_distY = numpy.ndarray(shape=(const.nbBinX2D,const.nbBinY2D))
	
	rootFile = ROOT.TFile(path)
	hXi1DMeanDist = rootFile.Get("hXi1DMeanDist")
	hXi1D         = rootFile.Get("hXi1D_")
	hXi2D         = rootFile.Get("hXi2D_")
	hXi2DMeanDistRparal2D = rootFile.Get("hXi2DMeanDistRparal2D")
	hXi2DMeanDistRperp2D  = rootFile.Get("hXi2DMeanDistRperp2D")

	for i in range(0,const.nbBin1D):
		xi1D__[0,i] = hXi1DMeanDist.GetBinContent(i+1)
		xi1D__[1,i] = -1.*hXi1D.GetBinContent(i+1)
		xi1D__[2,i] = hXi1D.GetBinError(i+1)

	for i in range(0,const.nbBin2D):
		xi2D__[0,i] = 0.
		xi2D__[1,i] = -1.*hXi2D.GetBinContent(i/const.nbBinY2D+1,i%const.nbBinY2D+1)
		xi2D__[2,i] = hXi2D.GetBinError(i/const.nbBinY2D+1,i%const.nbBinY2D+1)
		xi2D_2_distX[i/const.nbBinY2D,i%const.nbBinY2D] = hXi2DMeanDistRperp2D.GetBinContent(i/const.nbBinY2D+1,i%const.nbBinY2D+1)
		xi2D_2_distY[i/const.nbBinY2D,i%const.nbBinY2D] = hXi2DMeanDistRparal2D.GetBinContent(i/const.nbBinY2D+1,i%const.nbBinY2D+1)
	
	### Xi2D
	xi2D_2 = numpy.ndarray(shape=(const.nbBinX2D,const.nbBinY2D))
	for i in range(0,const.nbBin2D):
		xi2D_2[i/const.nbBinY2D,i%const.nbBinY2D] = xi2D__[1,i]

	if (draw):
		### Xi1D
		coef = -1.
		plt.errorbar(xi1D__[0,:],coef*xi1D__[1,:],yerr=coef*xi1D__[2,:],marker="o",label='all')
		plt.title(r'', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$\xi(|s|)$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.show()
		### Xi1D*|s|^{2}
		coef = -1.*xi1D__[0,:]*xi1D__[0,:]
		plt.errorbar(xi1D__[0,:],coef*xi1D__[1,:],yerr=coef*xi1D__[2,:],marker="o",label='all')
		plt.title(r'', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s|^{2}.\xi(|s|)$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.show()
		### Xi2D
		coef = -1.
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xticks([0.,50.,100.,150.])
		plt.imshow(coef*xi2D_2, origin='lower',extent=[const.minY2D, const.maxY2D, const.minX2D, const.maxX2D],interpolation='None')
		plt.title(r'', fontsize=40)
		plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()

		### Xi2D*|s|^{2}
		coef = -1.*numpy.sqrt(xi2D_2_distX*xi2D_2_distX + xi2D_2_distY*xi2D_2_distY)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xticks([0.,50.,100.,150.])
		plt.imshow(coef*xi2D_2, origin='lower',extent=[const.minY2D, const.maxY2D, const.minX2D, const.maxX2D],interpolation='None') #interpolation='None'
		plt.title(r'', fontsize=40)
		plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar = plt.colorbar()
	        cbar.set_label(r'$|s|^{2}.\xi(\, \overrightarrow{s} \,)$',size=40)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()

	return xi1D__, xi2D__, xi2D_2, xi2D_2_distX, xi2D_2_distY

def read_data_from_root():
	'''

	'''
	#XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_Flat
	rootFile = ROOT.TFile("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_Upgrade2/correlation_notDebug_allSky_DR12.root")
	#rootFile = ROOT.TFile("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m2_Prod_DR12_13_CoAddDelta/correlation_notDebug_allSky_DR12.root")

	hXi1D_        = rootFile.Get("hXi1D_")
	hXi1DMeanDist = rootFile.Get("hXi1DMeanDist")
	hXi1DMu1      = rootFile.Get("hXi1DMu1")
	hXi1DMu2      = rootFile.Get("hXi1DMu2")
	hXi1DMu3      = rootFile.Get("hXi1DMu3")
	nbBin = hXi1D_.GetNbinsX()

	binContent = numpy.zeros(nbBin)
	binError   = numpy.zeros(nbBin)
	binCenter  = numpy.zeros(nbBin)
	binMu1     = numpy.zeros(nbBin)
	binMu2     = numpy.zeros(nbBin)
	binMu3     = numpy.zeros(nbBin)
	ebinMu1    = numpy.zeros(nbBin)
	ebinMu2    = numpy.zeros(nbBin)
	ebinMu3    = numpy.zeros(nbBin)

	for i in range(0,nbBin):
		binContent[i] = hXi1D_.GetBinContent(i+1)
		binError[i]   = hXi1D_.GetBinError(i+1)
		binCenter[i]  = hXi1DMeanDist.GetBinContent(i+1)
		binMu1[i]     = hXi1DMu1.GetBinContent(i+1)
		binMu2[i]     = hXi1DMu2.GetBinContent(i+1)
		binMu3[i]     = hXi1DMu3.GetBinContent(i+1)
		ebinMu1[i]    = hXi1DMu1.GetBinError(i+1)
		ebinMu2[i]    = hXi1DMu2.GetBinError(i+1)
		ebinMu3[i]    = hXi1DMu3.GetBinError(i+1)

	coef = -1.
	plt.errorbar(binCenter,coef*binMu1,yerr=ebinMu1,marker="o",label='mu1')
	plt.errorbar(binCenter,coef*binMu2,yerr=ebinMu2,marker="o",label='mu2')
	plt.errorbar(binCenter,coef*binMu3,yerr=ebinMu3,marker="o",label='mu3')
	plt.errorbar(binCenter,coef*binContent,yerr=binError,marker="o",label='all')
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$\xi(|s|)$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	plt.errorbar(binCenter,coef*binMu1*binCenter,yerr=ebinMu1*binCenter,marker="o",label='mu1')
	plt.errorbar(binCenter,coef*binMu2*binCenter,yerr=ebinMu2*binCenter,marker="o",label='mu2')
	plt.errorbar(binCenter,coef*binMu3*binCenter,yerr=ebinMu3*binCenter,marker="o",label='mu3')
	plt.errorbar(binCenter,coef*binContent*binCenter,yerr=binError*binCenter,marker="o",label='all')
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$|s|.\xi(|s|)$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	plt.errorbar(binCenter,coef*binMu1*binCenter*binCenter,yerr=ebinMu1*binCenter*binCenter,marker="o",label='mu1')
	plt.errorbar(binCenter,coef*binMu2*binCenter*binCenter,yerr=ebinMu2*binCenter*binCenter,marker="o",label='mu2')
	plt.errorbar(binCenter,coef*binMu3*binCenter*binCenter,yerr=ebinMu3*binCenter*binCenter,marker="o",label='mu3')
	plt.errorbar(binCenter,coef*binContent*binCenter*binCenter,yerr=binError*binCenter*binCenter,marker="o",label='all')
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$|s|^{2}.\xi(|s|)$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	return

def covar_matrix():

	#XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_Upgrade2
	path     = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_allFlat2/correlation_notDebug_oneRegion_"
	nbRegion = 80
	nbBin1D  = 16

	nbBinX   = 32
	nbBinY   = 16
	nbBin2D  = nbBinX*nbBinY
	maxXI    = 160.

	### Get the overall xi1D
	rootFile = ROOT.TFile("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_allFlat2/source.root")
	hXi1D_        = rootFile.Get("hXi1D_")
	hXi1DMeanDist = rootFile.Get("hXi1DMeanDist")
	binCenter  = numpy.zeros(nbBin1D)
	binContent = numpy.zeros(nbBin1D)
	binError   = numpy.zeros(nbBin1D)

	for i in range(0,nbBin1D):
		binCenter[i]  = hXi1DMeanDist.GetBinContent(i+1)
		binContent[i] = hXi1D_.GetBinContent(i+1)
		binError[i]   = hXi1D_.GetBinError(i+1)

	boot_1D = numpy.ndarray(shape=(nbBin1D, nbRegion))
	boot_2D = numpy.ndarray(shape=(nbBin2D, nbRegion))

	for i in range(0,nbRegion):

		tmp_path = path + str(i) + "_over_" + str(nbRegion) + "_DR12.root"
		rootFile = ROOT.TFile(tmp_path)

		hXi1D_        = rootFile.Get("hXi1D_")
		hXi2D_        = rootFile.Get("hXi2D_")

		for j in range(0,nbBin1D):
			boot_1D[j,i] = hXi1D_.GetBinContent(j+1)

		for j in range(0,nbBin2D):
			boot_2D[j,i] = hXi2D_.GetBinContent(j/nbBinY+1,j%nbBinY+1)

	### 1D
	cov = numpy.cov(boot_1D)
	plt.imshow(cov, interpolation='None', origin='lower')
	plt.grid(True)
	plt.show()

	err_1D = [ numpy.sqrt(cov[i,i]/nbRegion) for i in range(0,nbBin1D) ]
	err_1D = numpy.asarray(err_1D)

	cov = numpy.corrcoef(boot_1D)
	plt.imshow(cov, interpolation='None', origin='lower')
	plt.grid(True)
	plt.show()

	### 2D
	cov = numpy.cov(boot_2D)
	plt.imshow(cov, interpolation='None', origin='lower')
	plt.grid(True)
	plt.show()

	cov = numpy.corrcoef(boot_2D)
	plt.imshow(cov, interpolation='None', origin='lower')
	plt.grid(True)
	plt.show()


	### Plot the Xi_1D
	coef = 1.*binCenter*binCenter
	plt.errorbar(binCenter,coef*binContent+1,yerr=coef*err_1D,marker="o",label='all err boot')
	plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$|s|.\xi(|s|)$', fontsize=40)
	myTools.deal_with_plot()
	plt.show()

	plt.errorbar(err_1D,binError,marker="o",label='all err boot')
	plt.xlabel(r'$err from cov$', fontsize=40)
	plt.ylabel(r'$binError$', fontsize=40)
	myTools.deal_with_plot()
	plt.show()

	print zip(err_1D,binError)

def Convert_Fits_To_Ntuples():
	'''

	Convert my fits file to root files

	'''
	

	path     = Get_Data_Fits_File()[1]
	file_cat = pyfits.open(path)
	cat      = file_cat[1].data[:1000]
	
	rootFile = ROOT.TFile("Results0_6999_Mock__DR11.root", "recreate") 
	tree = ROOT.TTree("SDSSQSO", "SDSS QSO")


	### Attributes of forest	
	plate     = array.array("i", [0])
	mjd       = array.array("i", [0])
	fiber     = array.array("i", [0])
	ra_J2000  = array.array("d", [0])
	dec_J2000 = array.array("d", [0])
	RedShift  = array.array("d", [0])
	NbPixel   = array.array("i", [0])

	### Array of forest
	LambdaAbs         = array.array("f", [0.]*nbBinRFMax__)
	DeltaFluxMeth2    = array.array("f", [0.]*nbBinRFMax__)
	DeltaFluxErrMeth2 = array.array("f", [0.]*nbBinRFMax__)
	zAbs              = array.array("f", [0.]*nbBinRFMax__)
	DLACorr           = array.array("f", [1.]*nbBinRFMax__)

	#print LambdaAbs

	tree.Branch("plate",     plate,     "plate/I")
	tree.Branch("mjd",       mjd,       "mjd/I")
	tree.Branch("fiber",     fiber,     "fiber/I")
	tree.Branch("ra_J2000",  ra_J2000,  "ra_J2000/D")
	tree.Branch("dec_J2000", dec_J2000, "dec_J2000/D")
	tree.Branch("RedShift",     RedShift,     "RedShift/D")
	tree.Branch("NbPixel",   NbPixel,   "NbPixel/I")

	tree.Branch("LambdaAbs",          LambdaAbs,         'LambdaAbs['+str(nbBinRFMax__)+']/F')
	tree.Branch("DeltaFluxMeth2",     DeltaFluxMeth2,    'DeltaFluxMeth2['+str(nbBinRFMax__)+']/F') 
	tree.Branch("DeltaFluxErrMeth2",  DeltaFluxErrMeth2, 'DeltaFluxErrMeth2['+str(nbBinRFMax__)+']/F');
	tree.Branch("zAbs",               zAbs,              'zAbs['+str(nbBinRFMax__)+']/F'); 
	tree.Branch("DLACorr",            DLACorr   ,        'DLACorr['+str(nbBinRFMax__)+']/F'); 

	for el in cat:

		plate[0]     = el['PLATE']
		mjd[0]       = el['MJD']
		fiber[0]     = el['FIBERID']
		ra_J2000[0]  = el['RA']
		dec_J2000[0] = el['DEC']
		RedShift[0]  = el['Z_VI']
		NbPixel[0]   = el['NB_PIXEL']

		for i in range(0,nbBinRFMax__):
			LambdaAbs[i]      = el['LAMBDA_OBS'][i]
			DeltaFluxMeth2[i] = el['DELTA'][i]
			zAbs[i]           = el['LAMBDA_OBS'][i]/lambdaRFLya__ -1.
			try:
				DeltaFluxErrMeth2[i] = numpy.power(el['DELTA_IVAR'][i],-0.5)
			except:
				DeltaFluxErrMeth2[i] = 0.

		tree.Fill() 

	rootFile.Write() 
	rootFile.Close()

	return


covar_matrix()






