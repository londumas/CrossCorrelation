# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >


import sys
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit
from scipy import interpolate
import myTools
from const_delta import *
import astropy.io.fits as pyfits
from scipy import interpolate

forest1__ = sys.argv[1]
forest2__ = sys.argv[2]
qso1__    = sys.argv[3]
qso2__    = sys.argv[4]

### 1D
min1D__   = 0.
max1D__   = 200.
nbBin1D__ = 50
binSize__ = max1D__/nbBin1D__

### 2D
minX2D__   = min1D__
maxX2D__   = max1D__
nbBinX2D__ = nbBin1D__

minY2D__   = min1D__
maxY2D__   = max1D__
nbBinY2D__ = nbBin1D__
nbBin2D__  = nbBinX2D__*nbBinY2D__

### Mu
nbBinM__ = 50;


#path1 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/new_2560-pixel-spectra/'



def loadData(path1D,coef):

	#print path1D
	data = numpy.loadtxt(path1D)

	### 1D
	save0 = data[:,0]
	save1 = data[:,1]

	#print '  <z> = ', numpy.sum(data[:,3])/numpy.sum(data[:,0])

	tmp_save0  = numpy.zeros(nbBin1D__)
	tmp_save1  = numpy.zeros(nbBin1D__)

	for i in range(0,len(save0)):
		idx = i/int(binSize__)

		tmp_save0[idx] += save0[i]
		tmp_save1[idx] += save1[i]

	### Fill arrays
	xxx = tmp_save1/tmp_save0
	yyy = tmp_save0/coef
	yer = numpy.sqrt(yyy)

	### Set to zero the empty fields
	cut = (yyy==0.)
	xxx[cut] = 0.
	yer[cut] = 0.


	return xxx, yyy, yer
def loadData2D(path2D,coef):

	data = numpy.loadtxt(path2D)

	### 2D
	save0 = data[:,0]
	save1 = data[:,1]
	save2 = data[:,2]

	xi2D   = numpy.zeros( shape=(nbBinX2D__,nbBinY2D__,3) )
	sPerp  = numpy.zeros( shape=(nbBinX2D__,nbBinY2D__) )
	sPara  = numpy.zeros( shape=(nbBinX2D__,nbBinY2D__) )

	for i in range(0,len(save0)):
		iX = i/int(maxY2D__-minY2D__)
		iY = i%int(maxY2D__-minY2D__)

		idX = iX/int(binSize__)
		idY = iY/int(binSize__)

		xi2D[idX][idY][1] += save0[i]
		sPerp[idX][idY]   += save1[i]
		sPara[idX][idY]   += save2[i]

	### Fill arrays
	xi2D[:,:,0] = numpy.sqrt( sPerp**2. + sPara**2. )/xi2D[:,:,1]
	xi2D[:,:,1] /= coef

	### Set to zero the empty fields
	cut = (xi2D[:,:,1]==0.)
	xi2D[:,:,1][cut] = 0.

	return xi2D
def loadDataMu(pathMu,coef):

	### Mu
	data = numpy.loadtxt(pathMu)

	save0 = data[:,0]
	save1 = data[:,1]

	tmp_save0 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save1 = numpy.zeros(shape=(nbBin1D__,nbBinM__))

	tmp_save00 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save11 = numpy.zeros(shape=(nbBin1D__,3))

	tmp_save000 = numpy.zeros(nbBin1D__)
	tmp_save111 = numpy.zeros(nbBin1D__)

	binSizeX = int(max1D__)/nbBin1D__
	binSizeY = 50/nbBinM__

	for i in range(0,len(save0)):
		iX = i/50
		iY = i%50

		### for mu
		idX = iX/binSizeX
		idY = iY/binSizeY

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]

		### for wedges
		if (iY>40):
			idY = 0
		elif (iY>25 and iY<=40):
			idY = 1
		else:
			idY = 2

		tmp_save00[idX][idY] += save0[i]
		tmp_save11[idX][idY] += save1[i]

		### for xi1D
		tmp_save000[idX] += save0[i]
		tmp_save111[idX] += save1[i]

	### Set arrays
	xiMu = numpy.zeros(shape=(nbBin1D__,nbBinM__,3))
	xiWe = numpy.zeros(shape=(nbBin1D__,3,3))
	xi1D = numpy.zeros(shape=(nbBin1D__,3))

	### Fill arrays
	xiMu[:,:,0] = tmp_save1/tmp_save0
	xiMu[:,:,1] = tmp_save0/coef

	xiWe[:,:,0] = tmp_save11/tmp_save00
	xiWe[:,:,1] = tmp_save00/coef

	xi1D[:,0] = tmp_save111/tmp_save000
	xi1D[:,1] = tmp_save000/coef
	xi1D[:,2] = numpy.sqrt(tmp_save000/coef)

	### Set to zero the empty fields
	cut = (tmp_save0 == 0.)
	xiMu[:,:,0][cut] = 0.
	xiMu[:,:,1][cut] = 0.

	cut = (tmp_save00 == 0.)
	xiWe[:,:,0][cut] = 0.
	xiWe[:,:,1][cut] = 0.

	cut = (tmp_save000 == 0.)
	xi1D[:,0][cut] = 0.
	xi1D[:,1][cut] = 0.

	return xi1D, xiMu, xiWe
def plot_Xi_1D(xi, rescale):

	cut = (xi[:,2]!=0.)
	xxx = xi[:,0][cut]
	yyy = xi[:,1][cut]
	yer = xi[:,2][cut]

	coef = numpy.power(xxx,rescale)
	plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|^{1}.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
	
	return
def plotXi2D(xi, rescale):

	xxx = xi[:,:,0]
	yyy = xi[:,:,1]
	cut = (yyy==0.)
	yyy[cut] = float('nan')

	xxx = numpy.transpose(xxx)
	yyy = numpy.transpose(yyy)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xticks([ 0.,50.,100.,150.,200.])
	ax.set_yticks([ 0.,50.,100.,150.,200.])

	coef = numpy.power(xxx,rescale)
	plt.imshow(coef*yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__], interpolation='None')
	cbar = plt.colorbar()

	if (rescale==0):
		cbar.set_label(r'$\xi^{qq}(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		cbar.set_label(r'$|s|.\xi^{qq}(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
	if (rescale==2):
		cbar.set_label(r'$|s|^{2}.\xi^{qq}(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

	plt.title(r'$'+qso1__+' \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
def plotMu(xiMu, rescale):

	xxx = xiMu[:,:,0]
	yyy = xiMu[:,:,1]

	cut = (yyy==0.)
	yyy[cut] = float('nan')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	#ax.set_yticks([ 0.,50.,100.,150.])

	coef = numpy.power(xxx,rescale)
	plt.imshow(coef*yyy, origin='lower', interpolation='None',extent=[0., 1.,min1D__, max1D__],aspect='auto')
	cbar = plt.colorbar()

	if (rescale==0):
		cbar.set_label(r'$\xi^{qq}(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		cbar.set_label(r'$|s|.\xi^{qq}(\, \overrightarrow{s}) \, [h^{-1}.Mpc]$',size=40)
	if (rescale==2):
		cbar.set_label(r'$|s|^{2}.\xi^{qq}(\, \overrightarrow{s}) \, [(h^{-1}.Mpc)^{2}]$',size=40)

	plt.title(r'$'+qso1__+' \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$\mu$', fontsize=40)
	plt.ylabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
def plotWe(xiWe, rescale):

	###
	cut = (xiWe[:,0,2]!=0.)
	xxx0 = xiWe[:,0,0][cut]
	yyy0 = xiWe[:,0,1][cut]
	yer0 = xiWe[:,0,2][cut]
	###
	cut = (xiWe[:,1,2]!=0.)
	xxx1 = xiWe[:,1,0][cut]
	yyy1 = xiWe[:,1,1][cut]
	yer1 = xiWe[:,1,2][cut]
	###
	cut = (xiWe[:,2,2]!=0.)
	xxx2 = xiWe[:,2,0][cut]
	yyy2 = xiWe[:,2,1][cut]
	yer2 = xiWe[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)
	coef2 = numpy.power(xxx2,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$0.8 < |\mu|$')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$0.5 < |\mu| \leq 0.8$')
	plt.errorbar(xxx2, coef2*yyy2, yerr=coef2*yer2, marker='o', label=r'$|\mu| \leq 0.5$')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=2)
	
	plt.title(r'$'+qso1__+' \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
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

		cut         = (yer[i,:]!=0.)
		tmpyyy      = yyy[i,:][cut]
		tmpyer      = yer[i,:][cut]
		xxxMu       = muu[cut]
		xxxMuPower1 = numpy.power(xxxMu,1.)
		xxxMuPower2 = numpy.power(xxxMu,2.)
		xxxMuPower3 = numpy.power(xxxMu,3.)
		xxxMuPower4 = numpy.power(xxxMu,4.)

		if (xxxMu.size<=1):
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

		color = ['blue','red','red']
		xxx = result_xi[:,0,0]

		### Show the result
		for i in range(0,3):
			coef = numpy.power(xxx,i)

			for j in range(0,5):
				if ( result_xi[:,j,1][ (result_xi[:,j,1]!=0.) ].size == 0): continue
				plt.errorbar(xxx, coef*result_xi[:,j,1],  yerr=coef*result_xi[:,j,2],  fmt='o', label=r'$\xi_{'+str(j)+'}$',color=color[j])
			#plt.errorbar(xxx, coef*xi1D_[:,1],  yerr=coef*xi1D_[:,2],  fmt='o', label=r'$\xi$')

			plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
			plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
			plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			myTools.deal_with_plot(False,False,True)
			plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
			plt.show()

	return result_xi
def fitCamb(data,pathToFile,mulpol=0):
	'''

	'''

	plot=True

	### Constants
	startFit = 30.
	maxForGues = 80.
	idx=1
	if (mulpol==2):
		idx = 2

	### Get the data
	cut = (data[:,0]>=startFit)
	xxx = data[:,0][cut]
	yyy = data[:,1][cut]
	yer = data[:,2][cut]

	### Get Camb
	data_camb = numpy.loadtxt(pathToFile)
	yyy_Camb = numpy.interp(xxx,data_camb[1:,0],data_camb[1:,idx])

	bPow2_init = -1. #numpy.mean(yyy[ (xxx<maxForGues) ]/yyy_Camb[ (xxx<maxForGues) ])
	print '  b_init = ', numpy.sqrt(bPow2_init)
	'''
	### Show result
	for i in range(0,3):

		coef = numpy.power(xxx,i)
		plt.errorbar(xxx,coef*yyy,yerr=coef*yer,fmt='o')

		coef = numpy.power(xxx,i)
		plt.plot(xxx,coef*bPow2_init*yyy_Camb,marker='o')
		myTools.deal_with_plot(False,False,False)
		plt.show()
	'''
	### Define the fit function
	def chi2(bPow2):
		return numpy.sum( numpy.power( (yyy-yyy_Camb*bPow2)/yer ,2.) )

	### Init and perform the fit
	
	m = Minuit(chi2, bPow2=bPow2_init,error_bPow2=0.1,print_level=-1, errordef=0.01) 	
	m.migrad()

	### Get result
	b = numpy.sqrt(m.values[ 'bPow2' ])
	print '  b = ', b

	### Print chi^2
	print '  DoF   = ', yyy.size, ' - ', 1
	print '  chi^2 = ', numpy.sum( numpy.power( (yyy-yyy_Camb*b*b)/yer ,2.) )
	result_Fit = numpy.asarray([b,yyy.size,numpy.sum( numpy.power( (yyy-yyy_Camb*b*b)/yer ,2.) )])

	### Get the data
	xxx = data[:,0]
	yyy = data[:,1]
	yer = data[:,2]

	### Get the smooth CAMB
	cut = (data_camb[1:,0] <= numpy.amax(xxx)*1.1)
	size = cut[cut].size
	result_1D_camb = numpy.zeros( shape=(size,3) )
	result_1D_camb[:,0] = data_camb[1:,0][cut]
	result_1D_camb[:,1] = b*b*data_camb[1:,idx][cut]
	result_1D_camb[:,2] = 0.0000000001

	if (plot):
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
	
	return result_Fit
def getACorrelation(nd,nr,nbRand,rawPath):
	'''

		Given a path to DD, DR, RR, return all the correlations.

		nd           : number of data
		nr           : number of random
		nbRand       : number of realisation of random
		rawPath : path to ASCII files

	'''

	### Which estimator
	LS = True
	### Plot or not
	plot = False

	print rawPath
	
	### Coef to normalize
	coefDD = nd*(nd-1.)/2.
	if (LS): coefDR = nd*nr
	coefRR = nr*(nr-1.)/2.

	### 2D: DD
	dd_2D = loadData2D(rawPath + 'xi_QSO_QSO_2D_QSO_DD.txt',coefDD)
	### Mu: DD
	dd_1D, dd_Mu, dd_We  = loadDataMu(rawPath + 'xi_QSO_QSO_Mu_QSO_DD.txt',coefDD)

	
	### 2D: RR,DR
	rr_2D = loadData2D(rawPath + 'xi_QSO_QSO_2D_QSO_RR_0.txt',coefRR)
	if (LS):
		dr_2D = loadData2D(rawPath + 'xi_QSO_QSO_2D_QSO_DR_0.txt',coefDR)
	### Mu: RR,DR
	rr_1D, rr_Mu, rr_We = loadDataMu(rawPath + 'xi_QSO_QSO_Mu_QSO_RR_0.txt',coefRR)
	if (LS):
		dr_1D, dr_Mu, dr_We = loadDataMu(rawPath + 'xi_QSO_QSO_Mu_QSO_DR_0.txt',coefDR)
	
	for i in range(1,nbRand):
	
		### 2D:
		rr_2D += loadData2D(rawPath + 'xi_QSO_QSO_2D_QSO_RR_'+str(i)+'.txt',coefRR)
	        if (LS):
			dr_2D +=  loadData2D(rawPath + 'xi_QSO_QSO_2D_QSO_DR_'+str(i)+'.txt',coefDR)
		### Mu:
		tmp_rr_1D, tmp_rr_Mu, tmp_rr_We = loadDataMu(rawPath + 'xi_QSO_QSO_Mu_QSO_RR_'+str(i)+'.txt',coefRR)
		rr_1D += tmp_rr_1D
		rr_Mu += tmp_rr_Mu
		rr_We += tmp_rr_We
		if (LS): 
			tmp_dr_1D, tmp_dr_Mu, tmp_dr_We = loadDataMu(rawPath + 'xi_QSO_QSO_Mu_QSO_DR_'+str(i)+'.txt',coefDR)
			dr_1D += tmp_dr_1D
			dr_Mu += tmp_dr_Mu
			dr_We += tmp_dr_We
		
	

	### 2D:
	rr_2D /= nbRand
	if (LS): dr_2D /= nbRand
	
	### Mu:
	rr_1D /= nbRand
	if (LS): dr_1D /= nbRand
	rr_Mu /= nbRand
	if (LS): dr_Mu /= nbRand
	rr_We /= nbRand
	if (LS): dr_We /= nbRand

	if (plot):
		### If plot
		plt.errorbar(rr_1D[:,0],rr_1D[:,1],fmt='o', label='RR')
		plt.errorbar(dr_1D[:,0],dr_1D[:,1],fmt='o', label='DR')
		plt.errorbar(dd_1D[:,0],dd_1D[:,1],fmt='o', label='DD')
		myTools.deal_with_plot(False,False,True)
		plt.show()

		### Plot 1D
		plt.errorbar(dd_1D[:,0], (dd_1D[:,0]-2.*dr_1D[:,0]+rr_1D[:,0])/dd_1D[:,0] )
		myTools.deal_with_plot(False,False,True)
		plt.show()

		### Plot 2D
		myTools.plot2D( (dd_2D[:,:,0]-2.*dr_2D[:,:,0]+rr_2D[:,:,0])/dd_2D[:,:,0] )

		### Plot 2D
		myTools.plot2D( (dd_Mu[:,:,0]-2.*dr_Mu[:,:,0]+rr_Mu[:,:,0])/dd_Mu[:,:,0] )

	### Get the Landy-Saley estimator
	if (LS): 
		### 1D:
		result_1D = numpy.array( dd_1D )
		result_1D[:,1] = (dd_1D[:,1]-2.*dr_1D[:,1]+rr_1D[:,1])/rr_1D[:,1]
		result_1D[:,2] = numpy.sqrt(dd_1D[:,1]*coefDD)/(coefDD*rr_1D[:,1])
		cut = (dd_1D[:,1]==0.)
		result_1D[:,1][cut] = 0.
		result_1D[:,2][cut] = 0.
	
		### 2D:
		result_2D = numpy.array( dd_2D )
		result_2D[:,:,1] = (dd_2D[:,:,1]-2.*dr_2D[:,:,1]+rr_2D[:,:,1])/rr_2D[:,:,1]
		cut = (dd_2D[:,:,1]==0.)
		result_2D[:,:,1][cut] = 0.
	
		### Mu:
		result_Mu = numpy.zeros(shape=(nbBin1D__,nbBinM__,3))
		result_Mu[:,:,0] = dd_Mu[:,:,0]
		result_Mu[:,:,1] = (dd_Mu[:,:,1] -2.*dr_Mu[:,:,1] + rr_Mu[:,:,1])/rr_Mu[:,:,1]
		result_Mu[:,:,2] = numpy.sqrt(dd_Mu[:,:,1]*coefDD)/(coefDD*rr_Mu[:,:,1])
		cut = (dd_Mu[:,:,1]==0.)
		result_Mu[:,:,2][cut] = 0.
		result_Mu[:,:,1][cut] = 0.
	
		### xiWe
		result_We = numpy.zeros(shape=(nbBin1D__,3,3))
		result_We[:,:,0] = dd_We[:,:,0]
		result_We[:,:,1] = (dd_We[:,:,1] -2.*dr_We[:,:,1] + rr_We[:,:,1])/rr_We[:,:,1]
		result_We[:,:,2] = numpy.sqrt(dd_We[:,:,1]*coefDD)/(coefDD*rr_We[:,:,1])
		cut = (dd_We[:,:,1]==0.)
		result_We[:,:,2][cut] = 0.
		result_We[:,:,1][cut] = 0.

	### Save
	numpy.save(rawPath+'xi_QSO_QSO_result_1D',result_1D)
	numpy.save(rawPath+'xi_QSO_QSO_result_2D',result_2D)
	numpy.save(rawPath+'xi_QSO_QSO_result_Mu',result_Mu)
	numpy.save(rawPath+'xi_QSO_QSO_result_We',result_We)

	result_Multipol = plotMultipol(result_Mu)
	numpy.save(rawPath+'xi_QSO_QSO_result_Multipol',result_Multipol)

	return result_1D, result_2D, result_Mu, result_We,result_Multipol
def saveOnetRealMocks(i,j):
	'''

	'''

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_withRSD.fits'
	cat = pyfits.open(path)[1].data
	nd = cat.size
	nr = cat.size
	del cat
	rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
	nbRand = 10
	result_1D,result_2D,result_Mu,result_We,result_Multipol = getACorrelation(nd,nr,nbRand,rawPath)

	return
def saveListRealMocks(ni,nj):
	'''

	'''

	pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_QSO_QSO_result_'

	
	### For fit
	cov = numpy.load(pathToSave+'cov_1D_backup.npy')
	arg_fit = {}
	arg_fit['a'] = 125.08255109750635
	arg_fit['b'] = -1.8611671002360046
	arg_fit['c'] = 0.006109139533639488
	arg_fit['A'] = 0.003731291955288748
	arg_fit['sigma'] = 9.870124119190448
	arg_fit['mean'] = 106.73882608124278
	arg_fit['fix_A'] = False
	arg_fit['fix_sigma'] = True
	arg_fit['fix_mean'] = True
	arg_fit['background'] = 'inv'
	arg_fit['cov'] = 'full'
	arg_fit['x_min'] = 60.
	arg_fit['x_max'] = 160.
	

	list1D       = numpy.zeros( shape=(nbBin1D__,ni*nj) )
	list2D       = numpy.zeros( shape=(nbBin2D__,ni*nj) )
	listMu       = numpy.zeros( shape=(nbBin1D__*nbBinM__,ni*nj) )
	listWe       = numpy.zeros( shape=(nbBin1D__,3,ni*nj) )
	listMultipol = numpy.zeros( shape=(nbBin1D__,5,ni*nj) )
	listBAO     = numpy.zeros( shape=(6,2,ni*nj) )

	for i in numpy.arange(ni):
		for j in numpy.arange(nj):
			print i,j
			xi1D       = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/xi_QSO_QSO_result_1D.npy')
			xi2D       = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/xi_QSO_QSO_result_2D.npy')
			xiMu       = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/xi_QSO_QSO_result_Mu.npy')
			xiWe       = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/xi_QSO_QSO_result_We.npy')
			xiMultipol = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/xi_QSO_QSO_result_Multipol.npy')
			list1D[:,i*10+j]         = xi1D[:,1]
			list2D[:,i*10+j]         = xi2D[:,:,1].flatten()
			listMu[:,i*10+j]         = xiMu[:,:,1].flatten()
			listWe[:,:,i*10+j]       = xiWe[:,:,1]
			listMultipol[:,:,i*10+j] = xiMultipol[:,:,1]

			a,b,c,d,e,f,g =  myTools.fit_BAO(xi1D[:,0],xi1D[:,1],cov,arg_fit)
			#plt.errorbar(xi1D[:,0],xi1D[:,0]**2.*xi1D[:,1],yerr=xi1D[:,0]**2.*xi1D[:,2],fmt='o')
			#plt.plot(c[0],c[0]**2.*c[1])
			#plt.plot(g[0],g[0]**2.*g[1])
			#plt.show()
			#print a['mean']
			print a
			listBAO[:,0,i*10+j] = numpy.array([ el for el in a.values() ])
			listBAO[:,1,i*10+j] = numpy.array([ el for el in b.values() ])

	numpy.save(pathToSave+'1D',list1D)
	numpy.save(pathToSave+'2D',list2D)
	numpy.save(pathToSave+'Mu',listMu)
	numpy.save(pathToSave+'We',listWe)
	numpy.save(pathToSave+'Multipol',listMultipol)
	numpy.save(pathToSave+'BAO',listBAO)

	cov1D = numpy.cov(list1D)
	cov2D = numpy.cov(list2D)
	covMu = numpy.cov(listMu)
	numpy.save(pathToSave+'cov_1D',cov1D)
	numpy.save(pathToSave+'cov_2D',cov2D)
	numpy.save(pathToSave+'cov_Mu',covMu)

	return






i = sys.argv[5]
j = sys.argv[6]
saveOnetRealMocks(i,j)


xi1D = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_QSO_QSO_result_1D.npy')
#plot_Xi_1D(xi1D,0)
#plot_Xi_1D(xi1D,1)
#plot_Xi_1D(xi1D,2)

fitCamb(xi1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Data/CAMB/CAMB_2_4/xi-z2.4.dat',0)
















