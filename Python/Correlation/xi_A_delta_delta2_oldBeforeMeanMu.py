# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >



import subprocess
import sys
import os
import numpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import cosmolopy.distance as cosmology ### Cosmology
from iminuit import Minuit


import sys
import myTools
from myTools import Get_TProfile
from const_delta import *


nbRegion = 80

### 1D
min1D   = 0.
max1D   = 200.
nbBin1D = 20
binSize = max1D/nbBin1D

### 2D
minX2D   = min1D
maxX2D   = max1D
nbBinX2D = nbBin1D

minY2D   = -max1D
maxY2D   = max1D
nbBinY2D = 2*nbBin1D

nbBin2D  = nbBinX2D*nbBinY2D

### Mu
nbBinM = 25;

forest1 = sys.argv[1]
forest2 = sys.argv[2]
qso1    = sys.argv[3]
qso2    = sys.argv[4]

path1 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test/'

def loadData(path1D, path2D):
	'''
	
		Resize the data from a bining of 1 Mpc.h-1 to a binSize
	
	'''

	
	### 2D
	data = numpy.loadtxt(path2D)
	print path2D

	save0 = data[:,0]
	save1 = data[:,1]
	save2 = data[:,2]
	save3 = data[:,3]
	save4 = data[:,4]
	save5 = data[:,5]
	save6 = data[:,6]

	tmp_save0 = numpy.zeros(shape=(nbBinX2D,nbBinY2D))
	tmp_save1 = numpy.zeros(shape=(nbBinX2D,nbBinY2D))
	tmp_save2 = numpy.zeros(shape=(nbBinX2D,nbBinY2D))
	tmp_save3 = numpy.zeros(shape=(nbBinX2D,nbBinY2D))
	tmp_save4 = numpy.zeros(shape=(nbBinX2D,nbBinY2D))
	tmp_save5 = numpy.zeros(shape=(nbBinX2D,nbBinY2D))
	tmp_save6 = numpy.zeros(shape=(nbBinX2D,nbBinY2D))

	for i in range(0,len(save0)):
		iX = i/int(maxY2D-minY2D)
		iY = i%int(maxY2D-minY2D)

		idX = iX/int(binSize)
		idY = iY/int(binSize)

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]
		tmp_save2[idX][idY] += save2[i]
		tmp_save3[idX][idY] += save3[i]
		tmp_save4[idX][idY] += save4[i]
		tmp_save5[idX][idY] += save5[i]
		tmp_save6[idX][idY] += save6[i]

	xi2D = numpy.zeros(shape=(nbBinX2D,nbBinY2D,3))

	for i in range(0,nbBinX2D):
		for j in range(0,nbBinY2D):

			if (tmp_save4[i][j] == 0.): continue

			xxx = numpy.sqrt( (tmp_save2[i][j]/tmp_save5[i][j])**2. + (tmp_save3[i][j]/tmp_save5[i][j])**2. )
			yyy = tmp_save0[i][j] / tmp_save5[i][j]
			yer = numpy.sqrt( (tmp_save1[i][j]/tmp_save5[i][j] - yyy*yyy)/tmp_save6[i][j] )
			xi2D[i][j][0] = xxx
			xi2D[i][j][1] = yyy
			xi2D[i][j][2] = yer

	### Mu
	print path1D
	data = numpy.loadtxt(path1D)

	save0 = data[:,0]
	save1 = data[:,1]
	save2 = data[:,2]
	save3 = data[:,3]
	save4 = data[:,4]
	save5 = data[:,5]

	tmp_save0 = numpy.zeros(shape=(nbBin1D,nbBinM))
	tmp_save1 = numpy.zeros(shape=(nbBin1D,nbBinM))
	tmp_save2 = numpy.zeros(shape=(nbBin1D,nbBinM))
	tmp_save3 = numpy.zeros(shape=(nbBin1D,nbBinM))
	tmp_save4 = numpy.zeros(shape=(nbBin1D,nbBinM))
	tmp_save5 = numpy.zeros(shape=(nbBin1D,nbBinM))

	tmp_save00 = numpy.zeros(shape=(nbBin1D,3))
	tmp_save11 = numpy.zeros(shape=(nbBin1D,3))
	tmp_save22 = numpy.zeros(shape=(nbBin1D,3))
	tmp_save33 = numpy.zeros(shape=(nbBin1D,3))
	tmp_save44 = numpy.zeros(shape=(nbBin1D,3))
	tmp_save55 = numpy.zeros(shape=(nbBin1D,3))

	tmp_save000 = numpy.zeros(nbBin1D)
	tmp_save111 = numpy.zeros(nbBin1D)
	tmp_save222 = numpy.zeros(nbBin1D)
	tmp_save333 = numpy.zeros(nbBin1D)
	tmp_save444 = numpy.zeros(nbBin1D)
	tmp_save555 = numpy.zeros(nbBin1D)

	binSizeX = int(max1D)/nbBin1D
	binSizeY = 100/nbBinM

	
	for i in range(0,len(save0)):
		iX = i/100
		iY = i%100

		### for mu
		idX = iX/binSizeX
		idY = iY/binSizeY

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]
		tmp_save2[idX][idY] += save2[i]
		tmp_save3[idX][idY] += save3[i]
		tmp_save4[idX][idY] += save4[i]
		tmp_save5[idX][idY] += save5[i]

		### for wedges
		if (iY<10 or iY>90):
			idY = 0
		elif ( (iY>=10 and iY<25) or (iY<=90 and iY>75) ):
			idY = 1
		else:
			idY = 2

		tmp_save00[idX][idY] += save0[i]
		tmp_save11[idX][idY] += save1[i]
		tmp_save22[idX][idY] += save2[i]
		tmp_save33[idX][idY] += save3[i]
		tmp_save44[idX][idY] += save4[i]
		tmp_save55[idX][idY] += save5[i]

		### for xi1D
		tmp_save000[idX] += save0[i]
		tmp_save111[idX] += save1[i]
		tmp_save222[idX] += save2[i]
		tmp_save333[idX] += save3[i]
		tmp_save444[idX] += save4[i]
		tmp_save555[idX] += save5[i]

	xiMu = numpy.zeros(shape=(nbBin1D,nbBinM,3))
	xiWe = numpy.zeros(shape=(nbBin1D,3,3))
	xi1D = numpy.zeros(shape=(nbBin1D,3))

	for i in range(0,nbBin1D):
		for j in range(0,nbBinM):

			if (tmp_save4[i][j] != 0.):
				xxx = tmp_save2[i][j]/tmp_save4[i][j]
				yyy = tmp_save0[i][j]/tmp_save4[i][j]
				yer = numpy.sqrt( (tmp_save1[i][j]/tmp_save4[i][j] - yyy*yyy)/tmp_save5[i][j] )
				xiMu[i][j][0] = xxx
				xiMu[i][j][1] = yyy
				xiMu[i][j][2] = yer
		for j in range(0,3):

			if (tmp_save44[i][j] != 0.):
				xxx = tmp_save22[i][j]/tmp_save44[i][j]
				yyy = tmp_save00[i][j]/tmp_save44[i][j]
				yer = numpy.sqrt( (tmp_save11[i][j]/tmp_save44[i][j] - yyy*yyy)/tmp_save55[i][j] )
				xiWe[i][j][0] = xxx
				xiWe[i][j][1] = yyy
				xiWe[i][j][2] = yer

		if (tmp_save444[i] != 0.):
			xxx = tmp_save222[i]/tmp_save444[i]
			yyy = tmp_save000[i]/tmp_save444[i]
			yer = numpy.sqrt( (tmp_save111[i]/tmp_save444[i] - yyy*yyy)/tmp_save555[i] )

			xi1D[i][0] = xxx
			xi1D[i][1] = yyy
			xi1D[i][2] = yer


	return xi1D, xi2D, xiMu, xiWe
def plot_Xi_1D(rescale):

	cut = (xi1D_[:,2] != 0.)
	xxx = xi1D_[:,0][cut]
	yyy = xi1D_[:,1][cut]
	yer = xi1D_[:,2][cut]

	if (rescale==0):
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o')
		plt.ylabel(r'$\xi (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, fmt='o')
		plt.ylabel(r'$|s|.\xi (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, fmt='o')
		plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1+'} \, - \, \delta_{'+forest2+'}$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
	
	return
def plot_Xi_2D(rescale):

	xxx = numpy.transpose(xi2D_[:,:,0])
	yyy = numpy.transpose(xi2D_[:,:,1])
	yer = numpy.transpose(xi2D_[:,:,2])
        
	fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xticks([ 0.,50.,100.,150.,200.])
	ax.set_yticks([ -200.,-150.,-100.,-50.,0.,50.,100.,150.,200.])

	if (rescale==0):
                plt.imshow(yyy, origin='lower',extent=[minX2D, maxX2D, minY2D, maxY2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
        if (rescale==1):
                plt.imshow(xxx*yyy, origin='lower',extent=[minX2D, maxX2D, minY2D, maxY2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|.\xi(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
        if (rescale==2):
                plt.imshow(xxx**2.*yyy, origin='lower',extent=[minX2D, maxX2D, minY2D, maxY2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|^{2}.\xi(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

        plt.title(r'$\delta_{'+forest1+'} \, - \, \delta_{'+forest2+'}$', fontsize=40)
        plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
        plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
        plt.grid(True)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
        plt.show()

        return
def plotMu(rescale):

	xxx = xiMu_[:,:,0]
	yyy = xiMu_[:,:,1]
	yer = xiMu_[:,:,2]

	fig = plt.figure()
	ax = fig.add_subplot(111)
	#ax.set_yticks([ 0.,50.,100.,150.])

	if (rescale==0):
		plt.imshow(yyy, origin='lower', interpolation='None',extent=[-1., 1.,min1D, max1D],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		plt.imshow(xxx*yyy, origin='lower', interpolation='None',extent=[-1., 1.,min1D, max1D],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \, [h^{-1}.Mpc])$',size=40)
	if (rescale==2):
		plt.imshow(xxx*xxx*yyy, origin='lower', interpolation='None',extent=[-1., 1.,min1D, max1D],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \, [(h^{-1}.Mpc)^{2}])$',size=40)


        plt.title(r'$\delta_{'+forest1+'} \, - \, \delta_{'+forest2+'}$', fontsize=40)
	plt.xlabel(r'$\mu$', fontsize=40)
	plt.ylabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
	
	return
def plotWe(rescale):

	###
	cut = (xiWe_[:,0,2]!=0.)
	xxx0 = xiWe_[:,0,0][cut]
	yyy0 = xiWe_[:,0,1][cut]
	yer0 = xiWe_[:,0,2][cut]
	#yyy0 -= yyy0[-1]
	###
	cut = (xiWe_[:,1,2]!=0.)
	xxx1 = xiWe_[:,1,0][cut]
	yyy1 = xiWe_[:,1,1][cut]
	yer1 = xiWe_[:,1,2][cut]
	#yyy1 -= yyy1[-1]
	###
	cut = (xiWe_[:,2,2]!=0.)
	xxx2 = xiWe_[:,2,0][cut]
	yyy2 = xiWe_[:,2,1][cut]
	yer2 = xiWe_[:,2,2][cut]
	#yyy2 -= yyy2[-1]

	if (rescale==0):
		plt.errorbar(xxx0, yyy0, yerr=yer0, fmt='o', label=r'$0.8 < |\mu|$')
		plt.errorbar(xxx1, yyy1, yerr=yer1, fmt='o', label=r'$0.5 < |\mu| \leq 0.8$')
		plt.errorbar(xxx2, yyy2, yerr=yer2, fmt='o', label=r'$|\mu| \leq 0.5$')
		plt.ylabel(r'$\xi (|s|)$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	if (rescale==1):
		plt.errorbar(xxx0, xxx0*yyy0, yerr=xxx0*yer0, fmt='o', label=r'$0.8 < |\mu|$')
		plt.errorbar(xxx1, xxx1*yyy1, yerr=xxx1*yer1, fmt='o', label=r'$0.5 < |\mu| \leq 0.8$')
		plt.errorbar(xxx2, xxx2*yyy2, yerr=xxx2*yer2, fmt='o', label=r'$|\mu| \leq 0.5$')
		plt.ylabel(r'$|s|.\xi (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	if (rescale==2):
		plt.errorbar(xxx0, xxx0*xxx0*yyy0, yerr=xxx0*xxx0*yer0, fmt='o', label=r'$0.8 < |\mu|$')
		plt.errorbar(xxx1, xxx1*xxx1*yyy1, yerr=xxx1*xxx1*yer1, fmt='o', label=r'$0.5 < |\mu| \leq 0.8$')
		plt.errorbar(xxx2, xxx2*xxx2*yyy2, yerr=xxx2*xxx2*yer2, fmt='o', label=r'$|\mu| \leq 0.5$')
		plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=2)
	
        plt.title(r'$\delta_{'+forest1+'} \, - \, \delta_{'+forest2+'}$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()


xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_A_delta_delta2_Mu_LYA_CIV.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_A_delta_delta2_2D_LYA_CIV.txt')
plot_Xi_1D(0)
plot_Xi_1D(1)
plot_Xi_1D(2)
plot_Xi_2D(0)
plot_Xi_2D(1)
plot_Xi_2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)

