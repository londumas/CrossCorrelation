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



import myTools
from myTools import Get_TProfile
from const_delta import *


nbRegion = 80

### 1D
min1D   = 0.
max1D   = 160.
nbBin1D = 32

### 2D
minX2D   = min1D
maxX2D   = max1D
nbBinX2D = nbBin1D

minY2D   = -max1D
maxY2D   = max1D
nbBinY2D = 2*nbBin1D

nbBin2D  = nbBinX2D*nbBinY2D


def loadData(binSize=1):
	'''
	
		Resize the data from a bining of 1 Mpc.h-1 to a binSize
	
	'''

	data1D = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/xi_delta_delta_1D_'+forest__+'.txt')

	save0 = data1D[:,3]
	save1 = data1D[:,4]
	save2 = data1D[:,5]
	save3 = data1D[:,6]
	save4 = data1D[:,7]

	tmp_save0 = numpy.zeros(nbBin1D)
	tmp_save1 = numpy.zeros(nbBin1D)
	tmp_save2 = numpy.zeros(nbBin1D)
	tmp_save3 = numpy.zeros(nbBin1D)
	tmp_save4 = numpy.zeros(nbBin1D)

	for i in range(0,len(save0)):
		idx = i/binSize
		tmp_save0[idx] += save0[i]
		tmp_save1[idx] += save1[i]
		tmp_save2[idx] += save2[i]
		tmp_save3[idx] += save3[i]
		tmp_save4[idx] += save4[i]

	xi1D = numpy.zeros(shape=(nbBin1D,3) )

	for i in range(0,nbBin1D):
		xxx = tmp_save2[i]/tmp_save3[i]
		yyy = tmp_save0[i]/tmp_save3[i]
		yer = numpy.sqrt( (tmp_save1[i]/tmp_save3[i] - yyy*yyy)/tmp_save4[i] )

		xi1D[i][0] = xxx
		xi1D[i][1] = yyy
		xi1D[i][2] = yer

	return xi1D
def plot_Xi_1D(rescale):

	xxx = xi1D_[:,0]
	yyy = xi1D_[:,1]
	yer = xi1D_[:,2]

	if (rescale==0):
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o')
		plt.ylabel(r'$\xi_{\delta,\delta \, Auto} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, fmt='o')
		plt.ylabel(r'$|s|.\xi_{\delta,\delta \, Auto} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, fmt='o')
		plt.ylabel(r'$|s|^{2}.\xi_{\delta,\delta \, Auto} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
	
	return
def plot_Xi_2D(rescale):

        ### If from cpp
        data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/xi_delta_delta_2D_'+forest__+'.txt')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xticks([ 0.,50.,100.,150.])
	ax.set_yticks([ 0.,50.,100.,150.])

        ### Xi2D
        xi2D = numpy.ndarray(shape=(nbBinX2D,nbBinX2D))
        xxx  = numpy.ndarray(shape=(nbBinX2D,nbBinX2D))
        for el in data:
                xi2D[int(el[0])/nbBinX2D,int(el[0])%nbBinX2D] = el[3]
                xxx[int(el[0])/nbBinX2D,int(el[0])%nbBinX2D] = numpy.sqrt(el[1]**2.+el[2]**2.)

        if (rescale==0):
                plt.imshow(xi2D, origin='lower',extent=[minX2D, maxX2D, minX2D, maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
        if (rescale==1):
                plt.imshow(xxx*xi2D, origin='lower',extent=[minX2D, maxX2D, minX2D, maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|.\xi(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
        if (rescale==2):
                plt.imshow(xxx**2.*xi2D, origin='lower',extent=[minX2D, maxX2D, minX2D, maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|^{2}.\xi(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

        plt.title(r'', fontsize=40)
        plt.ylabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
        plt.xlabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
        plt.grid(True)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
        plt.show()

        return

xi1D_ = loadData(max1D/nbBin1D)
plot_Xi_1D(0)
plot_Xi_1D(1)
plot_Xi_1D(2)
plot_Xi_2D(0)
plot_Xi_2D(1)
plot_Xi_2D(2)
