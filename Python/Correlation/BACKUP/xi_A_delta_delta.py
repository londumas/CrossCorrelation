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
import const

nbRegion = 80

### 1D
min1D   = 0.
max1D   = 200.
nbBin1D__ = 50
binSize = max1D/nbBin1D__

### 2D
minX2D   = min1D
maxX2D   = max1D
nbBinX2D__ = nbBin1D__

minY2D   = min1D
maxY2D   = max1D
nbBinY2D__ = nbBin1D__

nbBin2D__  = nbBinX2D__*nbBinY2D__

### Mu
nbBinMCalcul__ = 50
nbBinM__ = 25

forest__ = sys.argv[1]

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

	print '  <z> = ', numpy.sum(save4)/numpy.sum(save5)

	tmp_save0 = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save1 = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save2 = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save3 = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save4 = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save5 = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save6 = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))

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

	xi2D = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__,3))

	for i in range(0,nbBinX2D__):
		for j in range(0,nbBinY2D__):

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
	save6 = data[:,6]

	tmp_save0 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save1 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save2 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save3 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save4 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save5 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save6 = numpy.zeros(shape=(nbBin1D__,nbBinM__))

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

	binSizeX = int(max1D)/nbBin1D__
	binSizeY = nbBinMCalcul__/nbBinM__

	
	for i in range(0,len(save0)):
		iX = i/nbBinMCalcul__
		iY = i%nbBinMCalcul__

		### for mu
		idX = iX/binSizeX
		idY = iY/binSizeY

		#print i, idX, idY

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]
		tmp_save2[idX][idY] += save2[i]
		tmp_save3[idX][idY] += save3[i]
		tmp_save4[idX][idY] += save4[i]
		tmp_save5[idX][idY] += save5[i]
		tmp_save6[idX][idY] += save6[i]

		### for wedges
		if (iY>40):
			idY = 0
		elif (iY>25 and iY<=40):
			idY = 1
		else:
			idY = 2

		tmp_save00[idX][idY] += save0[i]
		tmp_save11[idX][idY] += save1[i]
		tmp_save22[idX][idY] += save2[i]
		tmp_save33[idX][idY] += save3[i]
		tmp_save44[idX][idY] += save4[i]
		tmp_save55[idX][idY] += save5[i]
		tmp_save66[idX][idY] += save6[i]

		### for xi1D
		tmp_save000[idX] += save0[i]
		tmp_save111[idX] += save1[i]
		tmp_save222[idX] += save2[i]
		tmp_save333[idX] += save3[i]
		tmp_save444[idX] += save4[i]
		tmp_save555[idX] += save5[i]
		tmp_save666[idX] += save6[i]

	xiMu = numpy.zeros(shape=(nbBin1D__,nbBinM__,3))
	xiWe = numpy.zeros(shape=(nbBin1D__,3,3))
	xi1D = numpy.zeros(shape=(nbBin1D__,3))

	for i in range(0,nbBin1D__):
		for j in range(0,nbBinM__):

			if (tmp_save4[i][j] != 0.):
				xxx = tmp_save2[i][j]/tmp_save5[i][j]
				yyy = tmp_save0[i][j]/tmp_save5[i][j]
				yer = numpy.sqrt( (tmp_save1[i][j]/tmp_save5[i][j] - yyy*yyy)/tmp_save6[i][j] )
				xiMu[i][j][0] = xxx
				xiMu[i][j][1] = yyy
				xiMu[i][j][2] = yer
		for j in range(0,3):

			if (tmp_save44[i][j] != 0.):
				xxx = tmp_save22[i][j]/tmp_save55[i][j]
				yyy = tmp_save00[i][j]/tmp_save55[i][j]
				yer = numpy.sqrt( (tmp_save11[i][j]/tmp_save55[i][j] - yyy*yyy)/tmp_save66[i][j] )
				xiWe[i][j][0] = xxx
				xiWe[i][j][1] = yyy
				xiWe[i][j][2] = yer

		if (tmp_save444[i] != 0.):
			xxx = tmp_save222[i]/tmp_save555[i]
			yyy = tmp_save000[i]/tmp_save555[i]
			yer = numpy.sqrt( (tmp_save111[i]/tmp_save555[i] - yyy*yyy)/tmp_save666[i] )

			xi1D[i][0] = xxx
			xi1D[i][1] = yyy
			xi1D[i][2] = yer


	return xi1D, xi2D, xiMu, xiWe
def plotXi(rescale):

	cut = (xi1D_[:,2] != 0.)
	xxx = xi1D_[:,0][cut]
	yyy = xi1D_[:,1][cut]
	yer = xi1D_[:,2][cut]

	#yyy -= yyy[-1]

	if (rescale==0):
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o')
		plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, fmt='o')
		plt.ylabel(r'$|s|.\xi^{ff} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, fmt='o')
		plt.ylabel(r'$|s|^{2}.\xi^{ff} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
	
	return
def plotXi2D(rescale):

	xxx = numpy.transpose(xi2D_[:,:,0])
	yyy = numpy.transpose(xi2D_[:,:,1])
	yer = numpy.transpose(xi2D_[:,:,2])

	yyy[ (yer==0.) ] = float('nan')
        
	fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xticks([ 0.,50.,100.,150.,200.])
	ax.set_yticks([ 0.,50.,100.,150.,200.])

	if (rescale==0):
                plt.imshow(yyy, origin='lower',extent=[minX2D, maxX2D, minX2D, maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$\xi^{ff}(\, \overrightarrow{s} \,)$',size=40)
        if (rescale==1):
                plt.imshow(xxx*yyy, origin='lower',extent=[minX2D, maxX2D, minX2D, maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|.\xi^{ff}(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
        if (rescale==2):
                plt.imshow(xxx**2.*yyy, origin='lower',extent=[minX2D, maxX2D, minX2D, maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|^{2}.\xi^{ff}(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

        plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
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

	yyy[ (yer==0.) ] = float('nan')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	#ax.set_yticks([ 0.,50.,100.,150.])

	if (rescale==0):
		plt.imshow(yyy, origin='lower', interpolation='None',extent=[0., 1.,min1D, max1D],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi^{ff}(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		plt.imshow(xxx*yyy, origin='lower', interpolation='None',extent=[0., 1.,min1D, max1D],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|.\xi^{ff}(\, \overrightarrow{s} \, [h^{-1}.Mpc])$',size=40)
	if (rescale==2):
		plt.imshow(xxx*xxx*yyy, origin='lower', interpolation='None',extent=[0., 1.,min1D, max1D],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|^{2}.\xi^{ff}(\, \overrightarrow{s} \, [(h^{-1}.Mpc)^{2}])$',size=40)


        plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
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
		plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=1)
	if (rescale==1):
		plt.errorbar(xxx0, xxx0*yyy0, yerr=xxx0*yer0, fmt='o', label=r'$0.8 < |\mu|$')
		plt.errorbar(xxx1, xxx1*yyy1, yerr=xxx1*yer1, fmt='o', label=r'$0.5 < |\mu| \leq 0.8$')
		plt.errorbar(xxx2, xxx2*yyy2, yerr=xxx2*yer2, fmt='o', label=r'$|\mu| \leq 0.5$')
		plt.ylabel(r'$|s|.\xi^{ff} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	if (rescale==2):
		plt.errorbar(xxx0, xxx0*xxx0*yyy0, yerr=xxx0*xxx0*yer0, fmt='o', label=r'$0.8 < |\mu|$')
		plt.errorbar(xxx1, xxx1*xxx1*yyy1, yerr=xxx1*xxx1*yer1, fmt='o', label=r'$0.5 < |\mu| \leq 0.8$')
		plt.errorbar(xxx2, xxx2*xxx2*yyy2, yerr=xxx2*xxx2*yer2, fmt='o', label=r'$|\mu| \leq 0.5$')
		plt.ylabel(r'$|s|^{2}.\xi^{ff} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=2)
	
        plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	plt.show()
def plotMultipol(xi):
	'''

	Plot the multipol of the cross-correlation

	'''

	### Plot or not
	plot = False
	### Get the data
	xxx = xi[:,:,0]
	muu = numpy.arange(0.,1.,1./nbBinM__) + 1./nbBinM__/2.
	yyy = xi[:,:,1]
	yer = xi[:,:,2]

	### Array to keep the results
	result_xi = numpy.zeros(shape=(nbBin1D__,5,3))
	for i in range(0,5):
		result_xi[:,i,0] = numpy.mean(xxx,axis=1)

	### Array with name of variable
	nameArray = [ 'xi0','xi1','xi2','xi3','xi4']

	for i in range(0,nbBin1D__):

		cut         = (xxx[i,:]!=0.)
		tmpyyy      = yyy[i,:][cut]
		tmpyer      = yer[i,:][cut]
		xxxMu       = muu[cut]
		xxxMuPower1 = numpy.power(xxxMu,1.)
		xxxMuPower2 = numpy.power(xxxMu,2.)
		xxxMuPower3 = numpy.power(xxxMu,3.)
		xxxMuPower4 = numpy.power(xxxMu,4.)

		if (xxxMu.size==0):
			continue
		
		
		### Define the fit function
		def chi2(xi0,xi1,xi2,xi3,xi4):
			fit = (xi0 - 0.5*xi2 + 0.375*xi4) + (xi1-1.5*xi3)*xxxMuPower1 + (1.5*xi2 - 3.75*xi4)*xxxMuPower2 + 2.5*xi3*xxxMuPower3 + 4.375*xi4*xxxMuPower4
			return numpy.sum( numpy.power( (tmpyyy-fit)/tmpyer ,2.) )

		### Init ad perform the fit
		m = Minuit(chi2, xi0=0.,error_xi0=0.1,xi1=0.,error_xi1=0.1,xi2=0.,error_xi2=0.1,xi3=0.,error_xi3=0.1,xi4=0.,error_xi4=0.1, print_level=-1, errordef=0.01,fix_xi1=True,fix_xi3=True,fix_xi4=True) 	
		m.migrad()

		### Stock the results
		for j in range(0,5): 
			result_xi[i,j,1] = m.values[ nameArray[j] ]
			result_xi[i,j,2] = m.errors[ nameArray[j] ]

		'''
		xi0 = m.values[ nameArray[0] ]
		xi1 = m.values[ nameArray[1] ]
		xi2 = m.values[ nameArray[2] ]
		xi3 = m.values[ nameArray[3] ]
		xi4 = m.values[ nameArray[4] ]
		fit = (xi0 - 0.5*xi2 + 0.375*xi4) + (xi1-1.5*xi3)*xxxMuPower1 + (1.5*xi2 - 3.75*xi4)*xxxMuPower2 + 2.5*xi3*xxxMuPower3 + 4.375*xi4*xxxMuPower4

		plt.errorbar(xxxMu, tmpyyy, yerr=tmpyer, fmt='o')
		plt.errorbar(xxxMu, fit)
		myTools.deal_with_plot(False,False,False)
		plt.show()
		
		
		### Residuals
		plt.errorbar(xxxMu, (yyy[i,:]-fit)/yer[i,:], fmt='o')
		'''

	if (plot):
		### Get the results
		xxx = result_xi[:,0,0]

		### Show the result
		for i in range(0,3):

			for j in range(0,5):
				if ( result_xi[:,j,1][ (result_xi[:,j,2]!=0.) ].size == 0): continue
			
				tmp_xxx = xxx[ (result_xi[:,j,1]!=0.) ]
				tmp_yyy = result_xi[:,j,1][ (result_xi[:,j,1]!=0.) ]
				tmp_yer = result_xi[:,j,2][ (result_xi[:,j,1]!=0.) ]
				coef    = numpy.power(tmp_xxx,i)
				plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$\xi_{'+str(j)+'}$')

			tmp_xxx = xi1D_[:,0][ (xi1D_[:,2]>0.) ]
			tmp_yyy = xi1D_[:,1][ (xi1D_[:,2]>0.) ]
			tmp_yer = xi1D_[:,2][ (xi1D_[:,2]>0.) ]
			coef    = numpy.power(tmp_xxx,i)
			plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$\xi$')

			plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
			plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
			plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
			myTools.deal_with_plot(False,False,True)
			plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
			plt.show()

	return result_xi
def fitCamb(data,pathToFile,mulpol=0):
	'''

	'''

	### Constants
	startFit   = 30.
	endFit     = 190.
	maxForGues = 60.
	idx=1
	if (mulpol==2):
		idx = 2


	### Get the data
	cut = numpy.logical_and( (data[:,0]>=startFit),(data[:,0]<endFit) )
	xxx = data[:,0][cut]
	yyy = data[:,1][cut]
	yer = data[:,2][cut]

	### Get Camb
	data_camb = numpy.loadtxt(pathToFile)
	yyy_Camb = numpy.interp(xxx,data_camb[1:,0],data_camb[1:,idx])

	b1b2_init = 1. #numpy.mean(yyy[ (xxx<maxForGues) ]/yyy_Camb[ (xxx<maxForGues) ])
	print '  b1b2_init = ', b1b2_init

	'''
	### Show result
	for i in range(0,3):

		coef = numpy.power(xxx,i)
		plt.errorbar(xxx,coef*yyy,yerr=coef*yer,fmt='o')

		coef = numpy.power(xxx,i)
		plt.plot(xxx,coef*b1b2_init*yyy_Camb,marker='o')
		myTools.deal_with_plot(False,False,False)
		plt.show()
	'''

	### Define the fit function
	def chi2(b1b2,roof):
		fit = yyy_Camb*b1b2+roof
		return numpy.sum( numpy.power( (yyy-fit)/yer ,2.) )

	### Init and perform the fit
	
	m = Minuit(chi2, b1b2=b1b2_init,error_b1b2=0.1, roof=0.,error_roof=0.1,print_level=-1, errordef=0.01,fix_roof=True) 	
	m.migrad()

	### Get result
	b1b2 = m.values[ 'b1b2' ]
	roof = m.values[ 'roof' ]
	print '  b1b2 = ', b1b2
	print '  roof = ', roof

	### Print chi^2
	print '  DoF   = ', yyy.size, ' - ', 1
	print '  chi^2 = ', numpy.sum( numpy.power( (yyy-yyy_Camb*b1b2-roof)/yer ,2.) )



	### Get the data
	xxx = data[:,0]
	yyy = data[:,1]
	yer = data[:,2]

	### Get the smooth CAMB
	cut = (data_camb[1:,0] <= numpy.amax(xxx)*1.1)
	size = cut[cut].size
	result_1D_camb = numpy.zeros( shape=(size,3) )
	result_1D_camb[:,0] = data_camb[1:,0][cut]
	result_1D_camb[:,1] = b1b2*data_camb[1:,idx][cut]+roof
	result_1D_camb[:,2] = 0.0000000001

	### Show result
	for i in range(0,3):

		coef = numpy.power(xxx,i)
		plt.errorbar(xxx,coef*yyy,yerr=coef*yer,fmt='o')

		coef = numpy.power(result_1D_camb[:,0],i)
		plt.plot(result_1D_camb[:,0],coef*result_1D_camb[:,1],color='red')

		if (i==0):
			plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
		if (i==1):
			plt.ylabel(r'$|s|^{1}.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (i==2):
			plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
		plt.title(r'', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		myTools.deal_with_plot(False,False,False)
		plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
		plt.show()
	
	return b1b2
def saveListRealMocks(ni,nj):
	'''

	'''

	path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
	pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_A_delta_delta_result_'

	list1D       = numpy.zeros( shape=(nbBin1D__,ni*nj) )
	list2D       = numpy.zeros( shape=(nbBin2D__,ni*nj) )
	listMu       = numpy.zeros( shape=(nbBin1D__*nbBinM__,ni*nj) )
	listWe       = numpy.zeros( shape=(nbBin1D__,3,ni*nj) )
	listMultipol = numpy.zeros( shape=(nbBin1D__,5,ni*nj) )

	for i in numpy.arange(ni):
		for j in numpy.arange(nj):
			tmpPath =  path + str(i)+'/Simu_00'+str(j)+'/Results/xi_A_delta_delta_'
			xi1D, xi2D, xiMu, xiWe = loadData(tmpPath+'Mu_LYA.txt',tmpPath+'2D_LYA.txt')
			list1D[:,i*10+j]         = xi1D[:,1]
			list2D[:,i*10+j]         = xi2D[:,:,1].flatten()
			listMu[:,i*10+j]         = xiMu[:,:,1].flatten()
			listWe[:,:,i*10+j]       = xiWe[:,:,1]
			listMultipol[:,:,i*10+j] = plotMultipol(xiMu)[:,:,1]

	numpy.save(pathToSave+'1D',list1D)
	numpy.save(pathToSave+'2D',list2D)
	numpy.save(pathToSave+'Mu',listMu)
	numpy.save(pathToSave+'We',listWe)
	numpy.save(pathToSave+'Multipol',listMultipol)

	cov1D = numpy.cov(list1D)
	cov2D = numpy.cov(list2D)
	covMu = numpy.cov(listMu)
	numpy.save(pathToSave+'cov_1D',cov1D)
	numpy.save(pathToSave+'cov_2D',cov2D)
	numpy.save(pathToSave+'cov_Mu',covMu)

	return










path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path+'xi_A_delta_delta_Mu_'+forest__+'.txt',path+'xi_A_delta_delta_2D_'+forest__+'.txt')

plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,const.pathToCamb__)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)

plotMu(0)
plotMu(1)
plotMu(2)

plotWe(0)
plotWe(1)
plotWe(2)



xi1DD_, xi2DD_, xiMuD_, xiWeD_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_test_PDFMocksJMC__z_0_2_5/xi_A_delta_delta_Mu_'+forest__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_test_PDFMocksJMC__z_0_2_5/xi_A_delta_delta_2D_'+forest__+'.txt')

for rescale in range(0,3):
	cut = (xi1D_[:,2] != 0.)
	xxx = xi1D_[:,0][cut]
	yyy = xi1D_[:,1][cut]
	yer = xi1D_[:,2][cut]

	#yyy -= yyy[-1]

	if (rescale==0):
		plt.errorbar(xxx, yyy, yerr=yer, marker='o',label='Simu 00')
		plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, marker='o',label='Simu 00')
		plt.ylabel(r'$|s|.\xi^{ff} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, marker='o',label='Simu 00')
		plt.ylabel(r'$|s|^{2}.\xi^{ff} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)

	cut = (xi1DD_[:,2] != 0.)
	xxx = xi1DD_[:,0][cut]
	yyy = xi1DD_[:,1][cut]
	yer = xi1DD_[:,2][cut]

	#yyy -= yyy[-1]

	if (rescale==0):
		plt.errorbar(xxx, yyy, yerr=yer, marker='o',label='Data')
		plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, marker='o',label='Data')
		plt.ylabel(r'$|s|.\xi^{ff} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, marker='o',label='Data')
		plt.ylabel(r'$|s|^{2}.\xi^{ff} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)


	
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()





plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,const.pathToCamb__)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)

plotMu(0)
plotMu(1)
plotMu(2)

plotWe(0)
plotWe(1)
plotWe(2)


listXi1D = []
name = []

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta_Mu_LYA.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta_2D_LYA.txt')
listXi1D += [xi1D_]
name += ['LYA-LYA']

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta_Mu_CIV.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta_2D_CIV.txt')
listXi1D += [xi1D_]
name += ['CIV-CIV']

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta_Mu_SIIV.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta_2D_SIIV.txt')
listXi1D += [xi1D_]
name += ['SiIV-SiIV']


nbBinMCalcul__ = 100
nbBinY2D__ = 2*nbBin1D__
minY2D   = -max1D
xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta2_Mu_LYA_CIV.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta2_2D_LYA_CIV.txt')
listXi1D += [xi1D_]
name += ['LYA-CIV']

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta2_Mu_LYA_SIIV.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta2_2D_LYA_SIIV.txt')
listXi1D += [xi1D_]
name += ['LYA-SiIV']

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta2_Mu_SIIV_CIV.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma/xi_A_delta_delta2_2D_SIIV_CIV.txt')
listXi1D += [xi1D_]
name += ['CIV-SiIV']


for rescale in [0,1,2]:
	for i in range(0, len(name)):

		xxx = listXi1D[i][:,0]
		yyy = listXi1D[i][:,1]
		yer = listXi1D[i][:,2]

		if (i!=0):
			if (rescale==0): print name[i], 1./( listXi1D[0][1,1]/yyy[1])
			yyy = yyy*listXi1D[0][1,1]/yyy[1]
			yer = yer*listXi1D[0][1,1]/yyy[1]
		else:
			if (rescale==0): print name[i], yyy[1]

		if (rescale==0):
			plt.errorbar(xxx, yyy, yerr=yer, marker='o',label=name[i])
			plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
		if (rescale==1):
			plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, marker='o',label=name[i])
			plt.ylabel(r'$|s|.\xi^{ff} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (rescale==2):
			plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, marker='o',label=name[i])
			plt.ylabel(r'$|s|^{2}.\xi^{ff} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()





ni = 10
nj = 10
saveListRealMocks(ni,nj)

'''
cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_A_delta_delta_result_cov_1D.npy')
myTools.plot2D(cov)
myTools.plot2D(myTools.getCorrelationMatrix(cov))
myTools.plotCovar([cov],['a'])

cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_A_delta_delta_result_cov_2D.npy')
myTools.plot2D(cov)
myTools.plot2D(myTools.getCorrelationMatrix(cov))
myTools.plotCovar([cov],['a'],50,50)
'''


#xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_Mu_LYA.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_2D_LYA.txt')
#xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_A_delta_delta_Mu_LYA.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_A_delta_delta_2D_LYA.txt')

'''
xi1DD, xi2DD, xiMuDD, xiWeDD = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_Mu_LYA.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_2D_LYA.txt')
result_Multipol_DD = plotMultipol(xiMuDD)


path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
tmpPath =  path + '0/Simu_000/Results/xi_A_delta_delta_'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(tmpPath+'Mu_LYA.txt',tmpPath+'2D_LYA.txt')
result_Multipol = plotMultipol(xiMu_)

rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/'

xi1D_[:,1] = numpy.mean(numpy.load(rawPath+'xi_A_delta_delta_result_1D.npy'),axis=1)
xi1D_[:,2] = numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_A_delta_delta_result_cov_1D.npy')))/numpy.sqrt(100.)

print numpy.mean(numpy.load(rawPath+'xi_A_delta_delta_result_2D.npy'),axis=1).size
xi2D_[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_A_delta_delta_result_2D.npy'),axis=1), nbBinX2D__,nbBinY2D__)
xi2D_[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_A_delta_delta_result_cov_2D.npy')))/numpy.sqrt(100.), nbBinX2D__,nbBinY2D__)

xiMu_[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_A_delta_delta_result_Mu.npy'),axis=1), nbBin1D__,nbBinM__)
xiMu_[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_A_delta_delta_result_cov_Mu.npy')))/numpy.sqrt(100.), nbBin1D__,nbBinM__)

result_Multipol[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_A_delta_delta_result_Multipol.npy'),axis=2)
result_Multipol[:,:,2] = numpy.var(numpy.load(rawPath+'xi_A_delta_delta_result_Multipol.npy'),axis=2)/numpy.sqrt(100.)

xiWe_[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_A_delta_delta_result_We.npy'),axis=2)
xiWe_[:,:,2] = numpy.var(numpy.load(rawPath+'xi_A_delta_delta_result_We.npy'),axis=2)/numpy.sqrt(100.)
'''

path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path+'xi_A_delta_delta_Mu_SIIV.txt',path+'xi_A_delta_delta_2D_SIIV.txt')

plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')



plotXi2D(0)
plotXi2D(1)
plotXi2D(2)

plotMu(0)
plotMu(1)
plotMu(2)

plotWe(0)
plotWe(1)
plotWe(2)

### Multipol
xi1D_[:,0] = result_Multipol[:,0,0]
xi1D_[:,1] = result_Multipol[:,0,1]
xi1D_[:,2] = result_Multipol[:,0,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
### Multipol
xi1D_[:,0] = result_Multipol[:,2,0]
xi1D_[:,1] = result_Multipol[:,2,1]
xi1D_[:,2] = result_Multipol[:,2,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat',2)
### Multipol

for rescale in range(0,3):
	###
	cut = (result_Multipol[:,0,2]!=0.)
	xxx0 = result_Multipol[:,0,0][cut]
	yyy0 = result_Multipol[:,0,1][cut]
	yer0 = result_Multipol[:,0,2][cut]
	###
	cut = (result_Multipol[:,2,2]!=0.)
	xxx1 = result_Multipol[:,2,0][cut]
	yyy1 = result_Multipol[:,2,1][cut]
	yer1 = result_Multipol[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Simu \, \xi_{0}$',color='red')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Simu \, \xi_{2}$',color='red')

	###
	cut = (result_Multipol_DD[:,0,2]!=0.)
	xxx0 = result_Multipol_DD[:,0,0][cut]
	yyy0 = result_Multipol_DD[:,0,1][cut]
	yer0 = result_Multipol_DD[:,0,2][cut]
	###
	cut = (result_Multipol_DD[:,2,2]!=0.)
	xxx1 = result_Multipol_DD[:,2,0][cut]
	yyy1 = result_Multipol_DD[:,2,1][cut]
	yer1 = result_Multipol_DD[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Data \, \xi_{0}$',color='blue')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Data \, \xi_{2}$',color='blue')

	if (rescale==0):
		plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{ff} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{ff} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()







### Numbers
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_1D_LYA.txt')
data[:,2] /= data[:,5]
plt.errorbar(data[:,2],data[:,6],fmt='o',color='blue',label='Data')

data_SAVE = numpy.array(numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_A_delta_delta_1D_LYA.txt'))
data_SAVE.fill(0.)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_A_delta_delta_1D_LYA.txt'
		data = numpy.loadtxt(path)
		data[:,2] /= data[:,5]
		data_SAVE += data
data_SAVE /= 100.
plt.errorbar(data_SAVE[:,2],data_SAVE[:,6],fmt='o',color='red',label='<Simu>')
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$nb \, pairs \, LYA-LYA$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()



### Errors
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_1D_LYA.txt')
data[:,2] /= data[:,5]
data[:,3] = data[:,1]/data[:,5] - (data[:,0]/data[:,5])**2.
plt.errorbar(data[:,2],data[:,3],fmt='o',color='blue',label='Data')

data_SAVE = numpy.array(numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_A_delta_delta_1D_LYA.txt'))
data_SAVE.fill(0.)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_A_delta_delta_1D_LYA.txt'
		data = numpy.loadtxt(path)
		data[:,2] /= data[:,5]
		data[:,3] = data[:,1]/data[:,5] - (data[:,0]/data[:,5])**2.
		data_SAVE += data
data_SAVE /= 100.
plt.errorbar(data_SAVE[:,2],data_SAVE[:,3],fmt='o',color='red',label='<Simu>')

plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$Var \, LYA-LYA$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

