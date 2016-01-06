# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# created by Hélion du Mas des Bourboux
# <Helion.du-Mas-des-Bourboux@cea.fr>

import numpy
import matplotlib.pyplot as plt
import pylab
import myTools


type__     = '2D'
nbParam__  = 16 #31
ni__ = 10
nj__ = 10
otherFit__ = True

### 1D
min1D__   = 0.
max1D__   = 200.
nbBin1D__ = 50
binSize__ = max1D__/nbBin1D__

### 2D
minX2D__   = min1D__
maxX2D__   = max1D__
nbBinX2D__ = nbBin1D__
minY2D__   = -max1D__
maxY2D__   = max1D__
nbBinY2D__ = 2*nbBin1D__
nbBin2D__  = nbBinX2D__*nbBinY2D__
nbBin__    = nbBinX2D__*nbBinY2D__

nbMocks__  = ni__*nj__
xxx__     = numpy.zeros(nbBin__)
muu__     = numpy.zeros(nbBin__)
xi_dat__  = numpy.zeros(nbBin__)
xi_fit__  = numpy.zeros(nbBin__)
xi_err__  = numpy.zeros(nbBin__)
xi_res__  = numpy.zeros(nbBin__)
xxxAll__        = numpy.zeros( shape=(nbBin__,nbMocks__) )
muuAll__        = numpy.zeros( shape=(nbBin__,nbMocks__) )
xi_datAll__     = numpy.zeros( shape=(nbBin__,nbMocks__) )
xi_fitAll__     = numpy.zeros( shape=(nbBin__,nbMocks__) )
xi_errAll__     = numpy.zeros( shape=(nbBin__,nbMocks__) )
xi_resAll__     = numpy.zeros( shape=(nbBin__,nbMocks__) )
param__         = numpy.zeros( shape=(nbParam__,2,nbMocks__) )
covar__         = numpy.zeros( shape=(nbParam__,nbParam__,nbMocks__) )
chi2__          = numpy.zeros( shape=(4,nbMocks__) )
chi2OtherFit__  = numpy.zeros( shape=(4,nbMocks__) )
paramName__     = numpy.asarray(['\\beta','b.(1+\\beta)','gamma-bias','gamma-beta','\\Delta v','b_{2}','b_{2}.\\beta_{2}',
				'BAO \, amplitude','\\alpha_{iso}','\\alpha_{\parallel}','\\alpha_{\perp}',
				'gamma-scale','Rad \, strength','Rad \, anisotropy','Rad \, mean \, free \, path',
				'Rad \, quasar \, lifetime','a0','a1','a2','a3','a4','a5','a6','a7','a7','a8','a9','a10','a11','a12','a13'])


def getResults():
	'''
	
		Set the array of results and of data
	
	'''

	for i in numpy.arange(ni__):
		for j in numpy.arange(nj__):
			idx = i*nj__+j

			pathData__     = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/BaoFit_q_f/bao' + type__
			path__         = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/BaoFit_q_f/bao' + type__
			pathOtherFit__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/BaoFit_q_f/bao' + type__

			name2 = '.data'
			name1 = '.residuals.dat'
			name3 = '.save.pars'
			name4 = '.fit.chisq'
			name5 = '.save.pcov'
	
			
			### Get data and fit
			data  = numpy.loadtxt(path__+name1)
			save0 = data[:,0].astype(int)
			save4 = data[:,4]
			save5 = data[:,5]
			save7 = data[:,7]
			save9 = data[:,9]
			xxxAll__[:,idx][save0]    = save4
			muuAll__[:,idx][save0]    = save5
			xi_fitAll__[:,idx][save0] = save7
			xi_errAll__[:,idx][save0] = save9
	
			### Get only data
			data = numpy.loadtxt(pathData__+name2)
			save0 = data[:,0].astype(int)
			save1 = data[:,1]
			xi_datAll__[:,idx][save0] = save1
			
			### Get the parameters of the fit
			data    = numpy.loadtxt(path__+name3)
			param__[:,0,idx] = data[:,1]
			param__[:,1,idx] = data[:,2]
	
			### Get the chi2 of the fit
			data = numpy.loadtxt(path__+name4)
			chi2__[:,idx] = data
			if (otherFit__):
				data = numpy.loadtxt(pathOtherFit__+name4)
				chi2OtherFit__[:,idx] = data


			### Get correlation parameters
			data = numpy.loadtxt(path__+name5)
			for el in data:
				covar__[el[0],el[1],idx] = el[2]

	numpy.save('param',param__)

	mean = numpy.average(param__[9,0,:],weights=numpy.power(param__[9,1,:],-2.))
	err = numpy.power(numpy.sum(numpy.power(param__[9,1,:],-2.)),-0.5)
	print '  alpha_parallel = ', mean, '  +/-  ', err
	mean = numpy.average(param__[10,0,:],weights=numpy.power(param__[10,1,:],-2.))
	err = numpy.power(numpy.sum(numpy.power(param__[10,1,:],-2.)),-0.5)
	print '  alpha_perp = ', mean, '  +/-  ', err

	mean = numpy.mean(param__[9,0,:])
	err = numpy.sqrt( numpy.var(param__[9,0,:])/(ni__*nj__))
	print '  alpha_parallel = ', mean, '  +/-  ', err
	mean = numpy.mean(param__[10,0,:])
	err = numpy.sqrt( numpy.var(param__[10,0,:])/(ni__*nj__))
	print '  alpha_perp = ', mean, '  +/-  ', err

	### Pass from covariance matrix to correlation matrix
	for i in numpy.arange(nbMocks__).astype(int):
		covar__[:,:,i] = myTools.getCorrelationMatrix(covar__[:,:,i])

	### Get the mean of all correlation-functions
	xi_resAll__[:] = (xi_datAll__-xi_fitAll__)/xi_errAll__
	xxx__[:]       = numpy.mean(xxxAll__,axis=1)[:]
	muu__[:]       = numpy.mean(muuAll__,axis=1)[:]
	xi_dat__[:]    = numpy.mean(xi_datAll__,axis=1)[:]
	xi_fit__[:]    = numpy.mean(xi_fitAll__,axis=1)[:]
	xi_err__[:]    = numpy.mean(xi_errAll__,axis=1)[:]/numpy.sqrt(nbMocks__)
	xi_res__[:]    = numpy.mean( (xi_datAll__-xi_fitAll__)/xi_errAll__,axis=1)[:]

	return
def plotParameters():
	'''
	'''

	####################################################################
	### 2D: alpha

	###
	plt.errorbar(param__[9,0,:],param__[10,0,:],xerr=param__[9,1,:],yerr=param__[10,1,:], fmt='o')
	plt.xlabel(r'$\alpha_{\parallel}$')
	plt.ylabel(r'$\alpha_{\perp}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

	###
	plt.errorbar((param__[9,0,:]-1.)/param__[9,1,:],(param__[10,0,:]-1.)/param__[10,1,:], linestyle="", marker="o")
	plt.xlabel(r'$(\alpha_{\parallel}-1.)/\sigma_{\alpha_{\parallel}}$')
	plt.ylabel(r'$(\alpha_{\perp}-1.)/\sigma_{\alpha_{\perp}}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()



	####################################################################
	### \alpha_{\parallel}

	### all:
	plt.errorbar(numpy.arange(nbMocks__),param__[9,0,:],yerr=param__[0,1,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$\alpha_{\parallel}$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### distrib:
	plt.hist(param__[9,0,:], bins=20)
	plt.xlabel(r'$\alpha_{\parallel}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### all: poll
	plt.errorbar(numpy.arange(nbMocks__),(param__[9,0,:]-1.)/param__[9,1,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$(\alpha_{\parallel}-1.)/\sigma$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### distrib: poll
	plt.hist((param__[9,0,:]-1.)/param__[9,1,:], bins=20)
	plt.xlabel(r'$(\alpha_{\parallel}-1.)/\sigma$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()



	####################################################################
	### \alpha_{\perp}

	### all:
	plt.errorbar(numpy.arange(nbMocks__),param__[10,0,:],yerr=param__[10,1,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$\alpha_{\perp}$')
	myTools.deal_with_plot(False,False,True)
	plt.xlim(-1,nbMocks__+1)
	plt.show()
	
	### distrib:
	plt.hist(param__[10,0,:], bins=20)
	plt.xlabel(r'$\alpha_{\perp}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### all: poll
	plt.errorbar(numpy.arange(nbMocks__),(param__[10,0,:]-1.)/param__[10,1,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$(\alpha_{\perp}-1.)/\sigma$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### distrib: poll
	plt.hist((param__[10,0,:]-1.)/param__[10,1,:], bins=20)
	plt.xlabel(r'$(\alpha_{\perp}-1.)/\sigma$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()



	####################################################################
	### Correlation between alpha

	### all:
	plt.errorbar(numpy.arange(nbMocks__),covar__[9,10,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$Correlation$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### distrib:
	plt.hist(covar__[9,10,:], bins=20)
	plt.xlabel(r'$Correlation$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()


	####################################################################
	### Delta-V
	
	plt.errorbar(numpy.arange(nbMocks__),param__[4,0,:],yerr=param__[4,1,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$\Delta \, v$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	####################################################################
	### beta
	
	plt.errorbar(numpy.arange(nbMocks__),param__[0,0,:],yerr=param__[0,1,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$\beta$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	####################################################################
	### bias2
	
	plt.errorbar(numpy.arange(nbMocks__),param__[5,0,:],yerr=param__[5,1,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$b_{2}$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()
def plotAllParameters():
	'''

	'''

	for i in numpy.arange(paramName__.size):

		### If the parameter isn't fitted
		if (param__[i,1,0] == 0.): continue

		### all:
		plt.errorbar(numpy.arange(nbMocks__),param__[i,0,:], yerr=param__[i,1,:], linestyle="", marker="o")
		plt.xlabel(r'$index \, mock$')
		plt.ylabel(r'$'+paramName__[i]+'$')
		plt.xlim(-1,nbMocks__+1)
		myTools.deal_with_plot(False,False,True)
		plt.show()

		### hist:
		yyy = param__[i,0,:]
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.hist(yyy, bins=20)
		myTools.deal_with_plot(False,False,False)
		plt.xlabel(r'$'+paramName__[i]+'$')
		plt.ylabel(r'$\#$')
		textstr = '$Nb=%d$ \n $\mu=%.5e  +/-  %.2e$ \n $\sigma=%.5e$'%(yyy.size,numpy.mean(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size),numpy.std(yyy))
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
		plt.show()

		'''
		### hist pull:
		yyy = (param__[i,0,:]-1.)/param__[i,1,:]
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.hist(yyy, bins=20)
		myTools.deal_with_plot(False,False,False)
		plt.xlabel(r'$('+paramName__[i]+'-1)/err$')
		plt.ylabel(r'$\#$')
		textstr = '$Nb=%d$ \n $\mu=%.5e  +/-  %.2e$ \n $\sigma=%.5e$'%(yyy.size,numpy.mean(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size),numpy.std(yyy))
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
		plt.show()
		'''
	return
def chi2Study():
	'''

	'''

	
	nbDOF = chi2__[1,0]-chi2__[2,0]
	print 
	print '  number of bins       = ', chi2__[1,0].astype(int)
	print '  number of parameters = ', chi2__[2,0].astype(int)
	print
	print

	
	####################################################################
	### chi^2

	### all:
	plt.errorbar(numpy.arange(nbMocks__),chi2__[0,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$\chi^{2}$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### hist:
	yyy = chi2__[0,:]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.hist(yyy, bins=20)
	myTools.deal_with_plot(False,False,False)
	plt.xlabel(r'$\chi^{2}$')
	plt.ylabel(r'$\#$')
	textstr = '$Nb=%d$ \n $\mu=%.5e  +/-  %.2e$ \n $\sigma=%.5e$'%(yyy.size,numpy.mean(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size),numpy.std(yyy))
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
	plt.show()

	
	####################################################################
	### chi^2_{1,1}

	nbDOF = chi2OtherFit__[1,0]-chi2OtherFit__[2,0]

	### all:
	plt.errorbar(numpy.arange(nbMocks__),chi2OtherFit__[0,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$\chi^{2}_{1,1}$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### hist:
	yyy = chi2OtherFit__[0,:]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.hist(yyy, bins=20)
	myTools.deal_with_plot(False,False,False)
	plt.xlabel(r'$\chi^{2}_{1,1}$')
	plt.ylabel(r'$\#$')
	textstr = '$Nb=%d$ \n $\mu=%.5e  +/-  %.2e$ \n $\sigma=%.5e$'%(yyy.size,numpy.mean(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size),numpy.std(yyy))
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
	plt.show()

	####################################################################
	### chi^2_{1,1}-chi^2

	### all:
	plt.errorbar(numpy.arange(nbMocks__),chi2OtherFit__[0,:]-chi2__[0,:], linestyle="", marker="o")
	plt.xlabel(r'$index \, mock$')
	plt.ylabel(r'$\chi^{2}_{1,1}-\chi^{2}$')
	plt.xlim(-1,nbMocks__+1)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### distrib:
	plt.hist(chi2OtherFit__[0,:]-chi2__[0,:], bins=20)
	plt.xlabel(r'$\chi^{2}_{1,1}-\chi^{2}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	

	####################################################################
	### chi^2_{1,1}-chi^2 vs. alpha

	### alpha_{parallel}:
	plt.errorbar(param__[9,0,:],chi2OtherFit__[0,:]-chi2__[0,:], linestyle="", marker="o")
	plt.xlabel(r'$\alpha_{\parallel}$')
	plt.ylabel(r'$\chi^{2}_{1,1}-\chi^{2}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### alpha_{perp}:
	plt.errorbar(param__[10,0,:],chi2OtherFit__[0,:]-chi2__[0,:], linestyle="", marker="o")
	plt.xlabel(r'$\alpha_{\perp}$')
	plt.ylabel(r'$\chi^{2}_{1,1}-\chi^{2}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	a = (chi2OtherFit__[0,:]-chi2__[0,:])
	print a[ a<0. ].size
	print a[ a<0. ]
	
	from scipy.stats import chisqprob
	### Get number under one sigma, two sigma ...
	nbParametrDiff = chi2__[2,0].astype(int)-chi2OtherFit__[2,0].astype(int)
	print '  parameter number difference  = ', nbParametrDiff
	tmp_a = chi2OtherFit__[0,:]-chi2__[0,:]

	print 1.-chisqprob(2.3, 2)
	sigma = [2.30, 6.18, 11.83, 19.33, 28.74]
	
	for i in range(0,5):
		print ' Delta Chi^2 = ', sigma[i], '  sigma = ', i+1, ' nbReal = ', tmp_a[ tmp_a<sigma[i] ].size, '  expected = ', tmp_a.size*(1.-chisqprob(sigma[i], 2)), chisqprob(sigma[i], 2)
	

def plotCorrelationFunction():

	####################################################################
	### Simu

	###
	for i in range(0,nbMocks__):
		plt.errorbar(numpy.arange(nbBin__),xi_dat__[:,i], fmt='o')
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$\xi_{Simu}(s)$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### mean:
	meanXi2D = numpy.mean(xi_dat__,axis=1)
	plt.errorbar(numpy.arange(nbBin__),meanXi2D, fmt='o')
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$<\xi_{Simu}(s)>$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### 2D:
	meanXi2D = myTools.convert1DTo2D(meanXi2D,nbBinY2D__,nbBinX2D__)
	myTools.plot2D(meanXi2D,[minX2D__, maxX2D__, minY2D__, maxY2D__],'s_{\perp} \, [h^{-1} Mpc]','s_{\parallel} \, [h^{-1} Mpc]','\\xi_{Simu}^{qf}')



	####################################################################
	### Fit

	###
	for i in range(0,nbMocks__):
		plt.errorbar(numpy.arange(nbBin__),xi_fit__[:,i], fmt='o')
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$\xi_{Fit}(s)$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### mean:
	meanXi2D = numpy.mean(xi_fit__,axis=1)
	plt.errorbar(numpy.arange(nbBin__),meanXi2D, fmt='o')
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$<\xi_{Fit}(s)>$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### 2D:
	meanXi2D = myTools.convert1DTo2D(meanXi2D,nbBinY2D__,nbBinX2D__)
	myTools.plot2D(meanXi2D,[minX2D__, maxX2D__, minY2D__, maxY2D__],'s_{\perp} \, [h^{-1} Mpc]','s_{\parallel} \, [h^{-1} Mpc]','\\xi_{Fit}^{qf}')

def plotDataAndFit_1D():
	'''
	'''
	
	if (type__=='1D'):
		xxx = numpy.arange(0., nbBin__ )*4. + 2.
		xxx[(xxx__!=0.)] = xxx__[(xxx__!=0.)]
		b = '|s|'
	elif (type__=='2D'):
		xxx = numpy.arange(0., nbBin__ )
		b = 's'
	
	for i in numpy.arange(0,3):
		
		if (i==0):
			a = ''
			c = ''
		elif (i==1):
			a = '|s|.'
			c = ' \, [h^{-1}.Mpc]'
		else:
			a = '|s|^{2}.'
			c = ' \, [(h^{-1}.Mpc)^{2}]'

		coef = numpy.power(xxx__,1.*i)
		
		plt.errorbar(xxx, coef*xi_dat__, yerr=coef*xi_err__, linestyle="", marker="o", color='blue')
		plt.errorbar(xxx, coef*xi_fit__, color='red')		
		plt.xlabel(r'$'+b+' \, [h^{-1}.Mpc] $')
		plt.ylabel(r'$'+a+'\\xi('+b+') '+c+'$')
		myTools.deal_with_plot(False,False,True)
		plt.show()
	
	return
def plotDataAndFit_2D():
	'''
	'''
	
	### but the xi from a flat array to a 2D array
	xxx     = myTools.convert1DTo2D(xxx__,nbBinY2D__,nbBinX2D__)
	yyy_dat = myTools.convert1DTo2D(xi_dat__,nbBinY2D__,nbBinX2D__)
	yyy_fit = myTools.convert1DTo2D(xi_fit__,nbBinY2D__,nbBinX2D__)
	yyy_err = myTools.convert1DTo2D(xi_err__,nbBinY2D__,nbBinX2D__)
	yyy_res = myTools.convert1DTo2D(xi_res__,nbBinY2D__,nbBinX2D__)
	
	yyy_dat[ (yyy_fit==0.) ] = float('nan')
	yyy_err[ (yyy_fit==0.) ] = float('nan')
	yyy_res[ (yyy_fit==0.) ] = float('nan')
	yyy_fit[ (yyy_fit==0.) ] = float('nan')
	edge = [0., 200., -200., 200.]
	
	### Plot the arrays
	for i in numpy.arange(0,3):

		coef = numpy.power(xxx,1.*i)
		
		a = ''
		if (i==1):
			a += '|s|.'
		elif (i==2):
			a += '|s|^{2}.'
			
		### data
		myTools.plot2D(coef*yyy_dat, edge, 's_{\\perp} \\, [h^{-1} Mpc]', 's_{\parallel} \\, [h^{-1} Mpc]', a+'\\xi_{Simu}(\\, \\overrightarrow{s} \\,)', 'Simu')
		### Fit
		myTools.plot2D(coef*yyy_fit, edge, 's_{\\perp} \\, [h^{-1} Mpc]', 's_{\parallel} \\, [h^{-1} Mpc]', a+'\\xi_{Fit}(\\, \\overrightarrow{s} \\,)', 'Fit')
		### residuals
		myTools.plot2D(coef*yyy_res, edge, 's_{\\perp} \\, [h^{-1} Mpc]', 's_{\parallel} \\, [h^{-1} Mpc]', a+'(\\xi_{Simu}-\\xi_{Fit})/\\xi_{error \, Simu}', 'residuals')
		
	return
def plotDataAndFit_2DSlice():
	'''
	'''
	
	### but the xi from a flat array to a 2D array
	xxx     = numpy.zeros( shape=(100,50) )
	yyy_dat = numpy.zeros( shape=(100,50) )
	yyy_fit = numpy.zeros( shape=(100,50) )
	yyy_err = numpy.zeros( shape=(100,50) )
	
	for i in range(0,nbBin__):
		xxx[i/50][i%50]     = xxx__[i]
		yyy_dat[i/50][i%50] = xi_dat__[i]
		yyy_fit[i/50][i%50] = xi_fit__[i]
		yyy_err[i/50][i%50] = xi_err__[i]
	
	yyy_dat[ (yyy_fit==0.) ] = float('nan')
	yyy_err[ (yyy_fit==0.) ] = float('nan')
	yyy_fit[ (yyy_fit==0.) ] = float('nan')
	
	
	
	
	### Slice of constant s_perp
	for i in range(0,50):
		tmpxxx = numpy.arange(0,100)*4.-200.+2.
		#plt.errorbar(tmpxxx, yyy_dat[:,i], yerr=yyy_err[:,i], fmt='o')
		plt.errorbar(tmpxxx, yyy_fit[:,i])
	
	plt.xlabel(r'$s_{\parallel} \, [Mpc.h^{-1}]$')
	plt.ylabel(r'$\xi(|s|)$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	### Slice of constant parallel
	for i in range(0,100):
		tmpxxx = numpy.arange(0,50)*4.+2.
		#plt.errorbar(tmpxxx, yyy_dat[i,:], yerr=yyy_err[i,:], fmt='o')
		plt.errorbar(tmpxxx, yyy_fit[i,:])
	
	plt.xlabel(r'$s_{\perp} \, [Mpc.h^{-1}]$')
	plt.ylabel(r'$\xi(|s|)$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	return
def plotDataAndFit_S():
	'''
	
		Plot \\xi( |s| )
	
	'''
	
	binSize = 4.
	
	sss     = numpy.arange(0.,50.)*4.+2.
	yyy_dat = numpy.zeros( shape=(50,2) )
	yyy_fit = numpy.zeros( shape=(50,2) )
	
	cut = numpy.logical_and( (xxx__!=0.), (xxx__<200.) )
	xxx    = xxx__[ cut ]
	xi_dat = xi_dat__[ cut ]
	xi_fit = xi_fit__[ cut ]
	xi_err = 1./(xi_err__[ cut ]**2.)
	
	for i in range(0,xxx.size):
		sIdx = int(xxx[i]/binSize)
		yyy_dat[sIdx][0] += xi_dat[i]*xi_err[i]
		yyy_dat[sIdx][1] += xi_err[i]
		yyy_fit[sIdx][0] += xi_fit[i]*xi_err[i]
		yyy_fit[sIdx][1] += xi_err[i]
	
	yyy_dat[:,0] /= yyy_dat[:,1]
	yyy_dat[:,1]  = numpy.sqrt(1./yyy_dat[:,1])
	yyy_fit[:,0] /= yyy_fit[:,1]
	yyy_fit[:,1]  = numpy.sqrt(1./yyy_fit[:,1])
	
	### Plot the results
	
	for i in numpy.arange(0,3):
		
		a = ''
		if (i==1):
			a += '|s|.'
		elif (i==2):
			a += '|s|^{2}.'

		coef = numpy.power(sss,1.*i)
		
		plt.errorbar(sss, coef*yyy_dat[:,0], yerr=coef*yyy_dat[:,1], linestyle="", marker="o", color='blue',label=r'$<Simu>$')
		plt.errorbar(sss, coef*yyy_fit[:,0], color='red',label=r'$<Fit>$')
		plt.xlabel(r'$|s| \, [h^{-1} Mpc]$')
		plt.ylabel(r'$'+a+'\\xi(|s|)$')
		myTools.deal_with_plot(False,False,True)
		
		plt.show()
		
	return
def plotDataAndFit_Mu():
	'''
	
		Plot \\xi( |s| )
	
	'''
	
	binSize = 4.
	
	sss     = numpy.arange(0.,50.)*4.+2.
	yyy_dat = numpy.zeros( shape=(50,3,2) )
	yyy_fit = numpy.zeros( shape=(50,3,2) )
	
	cut = numpy.logical_and( (xxx__!=0.), (xxx__<200.) )
	xxx    = xxx__[ cut ]
	muu    = numpy.abs(muu__[ cut ])
	xi_dat = xi_dat__[ cut ]
	xi_fit = xi_fit__[ cut ]
	xi_err = 1./(xi_err__[ cut ]**2.)
	
	for i in range(0,xxx.size):
		sIdx = int(xxx[i]/binSize)
		
		if (muu[i]>0.8):
			yyy_dat[sIdx][0][0] += xi_dat[i]*xi_err[i]
			yyy_dat[sIdx][0][1] += xi_err[i]
			yyy_fit[sIdx][0][0] += xi_fit[i]*xi_err[i]
			yyy_fit[sIdx][0][1] += xi_err[i]
		elif (muu[i]<=0.8 and muu[i]>0.5):
			yyy_dat[sIdx][1][0] += xi_dat[i]*xi_err[i]
			yyy_dat[sIdx][1][1] += xi_err[i]
			yyy_fit[sIdx][1][0] += xi_fit[i]*xi_err[i]
			yyy_fit[sIdx][1][1] += xi_err[i]
		else:
			yyy_dat[sIdx][2][0] += xi_dat[i]*xi_err[i]
			yyy_dat[sIdx][2][1] += xi_err[i]
			yyy_fit[sIdx][2][0] += xi_fit[i]*xi_err[i]
			yyy_fit[sIdx][2][1] += xi_err[i]
	
	for i in range(0,3):
		yyy_dat[:,i,0] /= yyy_dat[:,i,1]
		yyy_dat[:,i,1]  = numpy.sqrt(1./yyy_dat[:,i,1])
		yyy_fit[:,i,0] /= yyy_fit[:,i,1]
		yyy_fit[:,i,1]  = numpy.sqrt(1./yyy_fit[:,i,1])
	
	### Plot the results
	
	for i in range(0,3):
		
		a = ''
		if (i==1):
			a += '|s|.'
		elif (i==2):
			a += '|s|^{2}.'
			
		coef = numpy.power(sss,1.*i)
		
		for j in range(0,3):
			
			plt.errorbar(sss, coef*yyy_dat[:,j,0], yerr=coef*yyy_dat[:,j,1], linestyle="", marker="o",label=r'$<Simu>$')
			plt.errorbar(sss, coef*yyy_fit[:,j,0], color='red',label=r'$<Fit>$')
			plt.xlabel(r'$|s| \, [h^{-1} Mpc]$')
			plt.ylabel(r'$'+a+'\\xi(|s|)$')
			myTools.deal_with_plot(False,False,True)
			
		plt.show()
		
	return

def getHistoResiduals():
	'''
	
	'''

	###
	yyy = xi_resAll__[numpy.isfinite(xi_resAll__)]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.hist(yyy, bins=100)
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.xlabel(r'$(\xi_{simu}-\xi_{fit})/\sigma_{simu}$')
	plt.ylabel(r'$\#$')
	plt.rc('font', **{'size':'30'})
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.grid(True, which='both')
	mng = plt.get_current_fig_manager()
	textstr = '$Nb=%d$ \n $\mu=%.5e  +/-  %.2e$ \n $\sigma=%.5e  +/-  %.2e$'%(yyy.size,numpy.mean(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size),numpy.std(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size))
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
	plt.show()

	### Fit a gaussian
	myTools.fitAGaussian(yyy, 100, [1,0.,1.], False)
	myTools.fitAGaussian(yyy, 100, [1,0.,1.], True)
	
	### Constants
	nbBins = 100
	
	### Plot the residuals as a function of bin index
	xxx = numpy.arange(0., nbBin__ )
	yyy = xi_dat__-xi_fit__
	yyy[(xi_fit__==0.)] = 0.
	yer = xi_err__
	yer[(xi_fit__==0.)] = 0.
	plt.errorbar(xxx, yyy, yerr=yer, linestyle="", marker="o", color='blue')
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$\xi_{data}-\xi_{fit}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	### Plot the residuals/err as a function of bin index
	xxx = numpy.arange(0., nbBin__ )
	yyy = (xi_dat__-xi_fit__)/xi_err__
	yyy[(xi_fit__==0.)] = 0.
	plt.errorbar(xxx, yyy, linestyle="", marker="o", color='blue')
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$\frac{\xi_{data}-\xi_{fit}}{\sigma}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()	
		
	### Plot the distribution
	yyy = ( xi_dat__[ (xi_err__!=0.) ] - xi_fit__[ (xi_err__!=0.) ])/xi_err__[ (xi_err__!=0.) ]

	fig = plt.figure()
	ax = fig.add_subplot(111)

	ax.hist(yyy, bins=nbBins)
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	
	
	plt.xlabel(r'$\frac{data-fit}{\sigma_{data}}$')
	plt.ylabel(r'$\#$')
	
	plt.rc('font', **{'size':'30'})
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.grid(True, which='both')
	
	mng = plt.get_current_fig_manager()
	#mng.resize(*mng.window.maxsize())
	textstr = '$\mu=%.5e$\n$\sigma=%.5e$'%(numpy.mean(yyy), numpy.std(yyy))
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
	
	plt.show()
	
	#fitAGaussian(yyy, 50, [1,0.,1.],'\\frac{data-fit}{\\sigma_{data}}','\\#')
	
	return
def getChiScan():
	'''
		get the \chi^2 scan 
		1st line: 
		2nd line: settings
		3rd line: best fit
		4-end lines: scan
	'''
	
	### Constants
	sizeX = 50
	sizeY = 50
	name = '.scan.dat'

	### Arrays
	alphaPara = numpy.zeros(sizeX*sizeY)
	alphaPerp = numpy.zeros(sizeX*sizeY)
	chi2      = numpy.zeros(sizeX*sizeY)
	listAllFits = numpy.zeros( shape=(2,ni__*nj__) )

	for i in range(0,ni__):
		for j in range(0,nj__):

			path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f_covFromFit/bao' + type__
	
			### Create a file in /tmp with the scan minus the two first lines
			idx = 0
			f     = open(path+name)
			tmp_f = open('/tmp/scan.txt','w')
			for line in f:
				if (idx>=2): tmp_f.write(line)
				idx += 1
			tmp_f.close()
			f.close()
	
			### Data
			data = numpy.loadtxt('/tmp/scan.txt')
	
			### Best Fit
			chi2_bestFit = data[0][-1]
	
			### Scan
			alphaPara += data[1:,9]
			alphaPerp += data[1:,10]
			chi2      += data[1:,-1]-chi2_bestFit

	### All the bestFits
	listAllFits[0,:] =  param__[10,0,:]
	listAllFits[1,:] =  param__[9,0,:]

	alphaPara /= ni__*nj__
	alphaPerp /= ni__*nj__
	chi2      /= ni__*nj__
	
	### Put data in 2D array
	alphaPara_2D = numpy.zeros( shape=(sizeX,sizeY) )
	alphaPerp_2D = numpy.zeros( shape=(sizeX,sizeY) )
	chi2_2D      = numpy.zeros( shape=(sizeX,sizeY) )
	for k in numpy.arange(0,chi2.size):
		i = k/sizeY
		j = k%sizeY
		alphaPara_2D[i][j] = alphaPara[k]
		alphaPerp_2D[i][j] = alphaPerp[k]
		chi2_2D[i][j]      = chi2[k]


	listAllFits 
	
	
	edge = [0.7, 1.4, 0.7, 1.4]
	### Plot the data
	plotChi2Scan(chi2_2D,edge,True, listAllFits,'\\alpha_{\\perp}','\\alpha_{\\parallel}','\\Delta \\chi^{2} = \\chi^{2}-\\chi^{2}_{best \\, fit}')

	
	return
def plotChi2Scan(data, edge=None, contour=False, bestFit=None, xTitle='-',yTitle='-',zTitle='-',title='-'):
	'''
	'''
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	plt.imshow(data, origin='lower',extent=edge, interpolation='None') ##, interpolation='None'
	cbar = plt.colorbar()
	plt.xlabel(r'$'+xTitle+'$', fontsize=40)
	plt.ylabel(r'$'+yTitle+'$', fontsize=40)
	cbar.set_label(r'$'+zTitle+'$',size=40)
	
	### Limit for color
	plt.clim(0.,63.69)

	plt.grid(True)
	#cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	
	### [1s = 2.30, 2s = 6.18, 3s = 11.83, 4s = 19.33]
	### http://www.reid.ai/2012/09/chi-squared-distribution-table-with.html
	if (contour):
		pylab.contour(data, [2.30, 6.18, 11.83, 19.33, 28.74], extent=edge, linewidths=3, colors = 'red', origin='lower', hold='on')
	
	if (bestFit[0,:].size > 0):
		for i in range(0,bestFit[0,:].size):
			plt.scatter(bestFit[0,i],bestFit[1,i],marker='+',color='white',linewidths=5, hold='on')
			plt.xlim(edge[:2])
			plt.ylim(edge[2:])

		plt.scatter(0.90896, 1.04853, marker='+',color='red',linewidths=10, hold='on')
		plt.xlim(edge[:2])
		plt.ylim(edge[2:])

	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	return
def getChiScan1D():
	'''
		get the \chi^2 scan 
		1st line: 
		2nd line: settings
		3rd line: best fit
		4-end lines: scan
	'''

	### idex of the alpha
	idxAlpha = 9

	### keep the true errors (left and right)
	trueError = numpy.zeros( shape=(2,ni__*nj__) ) 
	
	### Constants
	name = '.scan.dat'

	for i in range(0,ni__):
		for j in range(0,nj__):

			path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f_covFromFit_fixedBAO_FixAlphaParal/bao' + type__
	
			### Create a file in /tmp with the scan minus the two first lines
			idx = 0
			f     = open(path+name)
			tmp_f = open('/tmp/scan.txt','w')
			for line in f:
				if (idx>=2): tmp_f.write(line)
				idx += 1
			tmp_f.close()
			f.close()
	
			### Data
			data = numpy.loadtxt('/tmp/scan.txt')
	
			### Best Fit
			alpha_bestFit = data[0][idxAlpha]
			chi2_bestFit  = data[0][-1]
			#print data[0][9], data[0][10], param__[idxAlpha,0,i*10+j], param__[idxAlpha,1,i*10+j]
	
			### Scan
			alpha = data[1:,idxAlpha]
			chi2  = data[1:,-1]-chi2_bestFit

			### Find the error at one sigma
			sigma = 1.
			size = alpha.size
			idxResult = numpy.asarray( [ k for k in range(0,size-1) if (chi2[k]-sigma)*(chi2[k+1]-sigma) <=0. ] )
			idxResult = [ (alpha[k]-alpha[k+1])/(chi2[k]-chi2[k+1])*sigma + alpha[k]-(alpha[k]-alpha[k+1])/(chi2[k]-chi2[k+1])*chi2[k] for k in idxResult ]
			a = idxResult-alpha_bestFit
			errMin = numpy.amax(a[ a<=0. ])
			errMax = numpy.amin(a[ a>0. ])
			print ' M1:  alpha = ', alpha_bestFit, '  +/- ', param__[idxAlpha,1,i*10+j]
			print ' M2:  alpha = ', alpha_bestFit, '  ', errMin, ' + ', errMax
			trueError[0,i*10+j] = errMin
			trueError[1,i*10+j] = errMax
			print errMin, errMax
			#if ((numpy.abs(errMin)/param__[idxAlpha,1,i*10+j] >1.2) or (errMax/param__[idxAlpha,1,i*10+j]> 1.2) or (numpy.abs(errMin)/param__[idxAlpha,1,i*10+j] <0.9) or (errMax/param__[idxAlpha,1,i*10+j]< 0.9)):

			if (True):
				plt.errorbar(alpha,chi2,fmt='o')
				plt.errorbar([alpha_bestFit],[0.],xerr=param__[idxAlpha,1,i*10+j],fmt='o',color='red',label='Best fit')
				plt.plot(alpha,numpy.ones(alpha.size)*1.,color='green')
				plt.plot(alpha,numpy.ones(alpha.size)*4.,color='green')
				plt.plot(alpha,numpy.ones(alpha.size)*9.,color='green')
				plt.xlabel(r'$\alpha_{\parallel}$', fontsize=40)
				plt.ylabel(r'$\Delta \chi^{2} = \chi^{2}-\chi^{2}_{best \, fit}$', fontsize=40)
				plt.ylim([ -2., numpy.amax(chi2)+2. ])
				myTools.deal_with_plot(False,False,True)
				plt.show()
			
	
	plt.hist(numpy.abs(trueError[0,:]),bins=20)
	plt.xlabel(r'$error \, left$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	plt.hist(trueError[1,:],bins=20)
	plt.xlabel(r'$error \, right$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	plt.hist(numpy.abs(trueError[0,:])/param__[idxAlpha,1,:],bins=20)
	plt.xlabel(r'$error \, left \, /  \,error \, BAOFIT$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	plt.hist(trueError[1,:]/param__[idxAlpha,1,:],bins=20)
	plt.xlabel(r'$error \, right \, / \, error \, BAOFIT$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	trueError[0,:] = numpy.abs(trueError[0,:])

	for i in range(0,ni__*nj__):
		if (param__[idxAlpha,0,i]<1.): param__[idxAlpha,1,i] = trueError[0,i]
		else: param__[idxAlpha,1,i] = trueError[1,i]

	return







getResults()
#getChiScan()
#getChiScan1D()

plotParameters()
#plotAllParameters()
#chi2Study()
plotDataAndFit_1D()

if (type__=='2D'):
	plotDataAndFit_2D()
	plotDataAndFit_2DSlice()
	plotDataAndFit_S()
	plotDataAndFit_Mu()

getHistoResiduals()

getChiScan1D()
getChiScan()












