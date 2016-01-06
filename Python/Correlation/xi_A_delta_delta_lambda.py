# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import myTools
from const_delta import *

import subprocess
import time
import sys
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit

nbRegion__ = 0

### 1D
min1D__ = 1.
max1D__ = 1.11
nbBin1D__ = 110

### Theta
minTheta_ = 0.;
maxTheta_ = 0.003;
nbBinM__ = 30;


path1__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'

### Parameters
forest__ = sys.argv[1]
qso__    = sys.argv[3]
if (len(sys.argv)>5):
	wickIdx__ = int(sys.argv[5])
if (len(sys.argv)>6):
	box__   = int(sys.argv[6])
	simul__ = int(sys.argv[7])


if (forest__=='LYA'):
	lines = LYA_lines
	names = LYA_lines_names
	lambdaRFLine     = 1215.67
if (forest__=='CIV'):
	lines = CIV_lines
	names = CIV_lines_names
	lambdaRFLine     = 1550.77845
if (forest__=='MGII'):
	lines = MGII_lines
	names = MGII_lines_names
	lambdaRFLine     = 2803.5324
if (forest__=='LYB'):
	lines = LYB_lines
	names = LYB_lines_names
	lambdaRFLine     = 1025.7223
if (forest__=='SIIV'):
	lines = SIIV_lines
	names = SIIV_lines_names
	lambdaRFLine     = 1402.77291

def loadData(path1D):

	print path1D

	### Mu
	data = numpy.loadtxt(path1D)

	save0 = data[:,0]
	save1 = data[:,1]
	save2 = data[:,2]
	save3 = data[:,3]
	save4 = data[:,4]
	save5 = data[:,5]
	save6 = data[:,6]

	nbBin1D__ = data[:,0].size/nbBinM__
	print nbBin1D__

	tmp_save0 = numpy.zeros(shape=(nbBinM__,nbBin1D__))
	tmp_save1 = numpy.zeros(shape=(nbBinM__,nbBin1D__))
	tmp_save2 = numpy.zeros(shape=(nbBinM__,nbBin1D__))
	tmp_save3 = numpy.zeros(shape=(nbBinM__,nbBin1D__))
	tmp_save4 = numpy.zeros(shape=(nbBinM__,nbBin1D__))
	tmp_save5 = numpy.zeros(shape=(nbBinM__,nbBin1D__))
	tmp_save6 = numpy.zeros(shape=(nbBinM__,nbBin1D__))

	tmp_save00 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save11 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save22 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save33 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save44 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save55 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save66 = numpy.zeros(shape=(nbBin1D__,3))

	tmp_save000 = numpy.zeros(nbBin1D__)
	tmp_save111 = numpy.zeros(nbBin1D__)
	tmp_save222 = numpy.zeros(nbBin1D__)
	tmp_save333 = numpy.zeros(nbBin1D__)
	tmp_save444 = numpy.zeros(nbBin1D__)
	tmp_save555 = numpy.zeros(nbBin1D__)
	tmp_save666 = numpy.zeros(nbBin1D__)
	
	for i in range(0,len(save0)):
		iX = i/nbBinM__
		iY = i%nbBinM__
		idY = iY

		tmp_save0[idY][iX] += save0[i]
		tmp_save1[idY][iX] += save1[i]
		tmp_save2[idY][iX] += save2[i]
		tmp_save3[idY][iX] += save3[i]
		tmp_save4[idY][iX] += save4[i]
		tmp_save5[idY][iX] += save5[i]
		tmp_save6[idY][iX] += save6[i]

		if (iY>2): continue

		tmp_save00[iX][idY] += save0[i]
		tmp_save11[iX][idY] += save1[i]
		tmp_save22[iX][idY] += save2[i]
		tmp_save33[iX][idY] += save3[i]
		tmp_save44[iX][idY] += save4[i]
		tmp_save55[iX][idY] += save5[i]
		tmp_save66[iX][idY] += save6[i]

		### for xi1D
		tmp_save000[iX] += save0[i]
		tmp_save111[iX] += save1[i]
		tmp_save222[iX] += save2[i]
		tmp_save333[iX] += save3[i]
		tmp_save444[iX] += save4[i]
		tmp_save555[iX] += save5[i]
		tmp_save666[iX] += save6[i]
	
	xiMu = numpy.zeros(shape=(nbBinM__,nbBin1D__,4))
	xiMu[:,:,0] = tmp_save2/tmp_save5
	xiMu[:,:,1] = tmp_save3/tmp_save5
	xiMu[:,:,2] = tmp_save0/tmp_save5
	xiMu[:,:,3] = numpy.sqrt( (tmp_save1/tmp_save5 - xiMu[:,:,2]*xiMu[:,:,2])/tmp_save6)
	cut = (tmp_save5==0.)
	xiMu[:,:,0][cut] = 0.
	xiMu[:,:,1][cut] = 0.
	xiMu[:,:,2][cut] = 0.
	xiMu[:,:,3][cut] = 0.

	xiWe = numpy.zeros(shape=(nbBin1D__,3,3))
	xiWe[:,:,0] = tmp_save22/tmp_save55
	xiWe[:,:,1] = tmp_save00/tmp_save55
	xiWe[:,:,2] = numpy.sqrt( (tmp_save11/tmp_save55 - xiWe[:,:,1]*xiWe[:,:,1])/tmp_save66 )
	cut = (tmp_save55==0.)
	xiWe[:,:,0][cut] = 0.
	xiWe[:,:,1][cut] = 0.
	xiWe[:,:,2][cut] = 0.

	xi1D = numpy.zeros(shape=(nbBin1D__,3))
	xi1D[:,0] = tmp_save222/tmp_save555
	xi1D[:,1] = tmp_save000/tmp_save555
	xi1D[:,2] = numpy.sqrt( (tmp_save111/tmp_save555 - xi1D[:,1]*xi1D[:,1])/tmp_save666 )
	cut = (tmp_save555==0.)
	xi1D[:,0][cut] = 0.
	xi1D[:,1][cut] = 0.
	xi1D[:,2][cut] = 0.

	return xi1D, xiMu, xiWe
def plotXi():

	xxx = xi1D_[:,0]
	yyy = xi1D_[:,1]
	yer = xi1D_[:,2]

	cut = (yer!=0.)
	xxx = xxx[ cut ]
	yyy = yyy[ cut ]
	yer = yer[ cut ]

	yMin    = numpy.min(yyy)
	yMax    = numpy.max(yyy)
	nbLines = lines.size
	for i in range(0,nbLines):
		for j in range(0,i):
			if ( lines[i]!=lambdaRFLine and lines[j]!=lambdaRFLine ): continue
			line = max( lines[i]/lines[j], lines[j]/lines[i])
			if (line<min1D__ or line>max1D__): continue
			xLi  = [line,line]
			yLi  = [yMin,yMax]
			name = names[i] + ' - ' + names[j]
			plt.plot(xLi,yLi,color='green',linewidth=2)
			plt.text(line, yMax, name, rotation='vertical', fontsize=20)

	plt.errorbar(xxx, yyy, yerr=yer, fmt='o')
	plt.ylabel(r'$\xi^{ff}$', fontsize=40)
	plt.xlabel(r'$\lambda_{Obs.,2}/\lambda_{Obs.,1}$', fontsize=40)
	plt.xlim([ min1D__-0.01, max1D__+0.01 ])
	myTools.deal_with_plot(False,False,False)

	plt.show()
def plotMu():

	import matplotlib as mpl

	xxx = xiMu_[:,:,0]
	muu = xiMu_[:,:,1]
	yyy = xiMu_[:,:,2]
	yer = xiMu_[:,:,3]

	yyy[ (yer==0.) ] = float('nan')
	yer[ (yer==0.) ] = float('nan')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	extent=[min1D__, max1D__, minTheta_, maxTheta_]

	plt.imshow(yyy, origin='lower',aspect='auto',extent=extent,interpolation='None',vmin=-0.002,vmax=0.002)  ##interpolation='None'
	cbar = plt.colorbar()
	plt.xlabel(r'$\lambda_{Obs.,2}/\lambda_{Obs.,1}$', fontsize=40)
	plt.ylabel(r'$\theta \, [rad]$', fontsize=40)
	cbar.set_label(r'$\xi^{ff}$',size=40)
	plt.grid(True)

	plt.xlim([ min1D__, max1D__ ])
	plt.ylim([ minTheta_, maxTheta_ ])


	yMin    = numpy.min(muu)
	yMax    = numpy.max(muu)
	nbLines = lines.size
	for i in range(0,nbLines):
		for j in range(0,i):
			if ( lines[i]!=lambdaRFLine and lines[j]!=lambdaRFLine ): continue
			line = max( lines[i]/lines[j], lines[j]/lines[i])
			if (line<min1D__ or line>max1D__): continue
			xLi  = [line,line]
			yLi  = [yMin,yMax]
			name = names[i] + ' - ' + names[j]
			plt.plot(xLi,yLi,color='green',linewidth=2)
			plt.text(line, yMax, name, rotation='vertical', fontsize=20)


	plt.show()
def plotWe(rescale):

	a = ["theta < 0.0001", "0.0001 < theta < 0.0002"]

	yMin    = 0.
	yMax    = 0.005
	nbLines = lines.size
	for i in range(0,nbLines):
		for j in range(0,i):
			#if ( lines[i]!=1550.77845 and lines[j]!=1550.77845 and lines[i]!=1548.2049 and lines[j]!=1548.2049 ): continue
			#if ( lines[i]!=1402.77291 and lines[j]!=1402.77291 and lines[i]!=1393.76018 and lines[j]!=1393.76018 ): continue
			if ( lines[i]!=1215.67 and lines[j]!=1215.67): continue
			line = max( lines[i]/lines[j], lines[j]/lines[i])
			if (line<min1D__ or line>max1D__): continue
			xLi  = [line,line]
			yLi  = [yMin,yMax]
			name = names[i] + ' - ' + names[j]
			plt.plot(xLi,yLi,color='green',linewidth=2)
			plt.text(line, yMax, name, rotation='vertical', fontsize=20)

	for i in range(0,2):

		#if (i>=1):
		#	xiWe_[:,i,1] += 0.05
		#if (i>=2):
		#	xiWe_[:,i,1] += 0.05	

		###
		cut = (xiWe_[:,i,2]!=0.)

		if (xiWe_[:,i,0][cut].size==0):
			continue

		xxx = xiWe_[:,i,0][cut]
		yyy = xiWe_[:,i,1][cut]
		yer = xiWe_[:,i,2][cut]
		coef = numpy.power(xxx,rescale)
		

		plt.errorbar(xxx, coef*yyy, yerr=coef*yer, marker='o', label=r'$'+a[i]+'$')
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=1)
	
	plt.ylabel(r'$\xi^{ff}$', fontsize=40)
	plt.xlabel(r'$\lambda_{Obs.,2}/\lambda_{Obs.,1}$', fontsize=40)
	plt.xlim([ min1D__-0.01, max1D__+0.01 ])
	myTools.deal_with_plot(False,False,False)

	


	plt.show()

xi1D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/xi_A_delta_delta_lambda_Mu_'+forest__+'.txt')

plotXi()
plotMu()
plotWe(0)



























