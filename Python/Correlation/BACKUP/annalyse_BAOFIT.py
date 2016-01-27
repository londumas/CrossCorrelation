# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# created by Hélion du Mas des Bourboux
# <Helion.du-Mas-des-Bourboux@cea.fr>

import numpy
import matplotlib.pyplot as plt
import pylab
import myTools
from scipy import interpolate

type__        = '2D'
path__        = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/BaoFit_q_f__LYA__QSO/bao' + type__ + ''  ##_scanAlphaParal  scanAlphaPerp ##SmalScan  #middleScan
pathData__    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/BaoFit_q_f__LYA__QSO/bao' + type__
#path__        = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation_correctedForest_withMoreMetals/Box_000/Simu_000/Results/BaoFit_q_f__LYA__QSO/bao' + type__
#pathData__    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation_correctedForest_withMoreMetals/Box_000/Simu_000/Results/BaoFit_q_f__LYA__QSO/bao' + type__

nbParam__ = 30
nbBin__   = 5000
xxx__    = numpy.zeros( nbBin__ )
muu__    = numpy.zeros( nbBin__ )
xi_dat__ = numpy.zeros( nbBin__ )
xi_fit__ = numpy.zeros( nbBin__ )
xi_err__ = numpy.zeros( nbBin__ )
param__  = numpy.zeros( shape=(nbParam__,2) )
chi2__   = numpy.zeros( 4 )
paramName     = numpy.asarray(['\\beta','b.(1+\\beta)','gamma-bias','gamma-beta','\\Delta v','b_{2}','b_{2}.\\beta_{2}',
				'1+f','SigmaNL-perp',
				'BAO \, amplitude','\\alpha_{iso}','\\alpha_{\parallel}','\\alpha_{\perp}',
				'gamma-scale','Rad \, strength','Rad \, anisotropy','Rad \, mean \, free \, path',
				'Rad \, quasar \, lifetime','a0','a1','a2','a3','a4','a5','a6','a7','a7','a8','a9','a10','a11','a12','a13'])
paramName__ = paramName[:nbParam__]


def getResults():
	'''
	
		Set the array of results and of data
	
	'''
	
	name2 = '.data'
	name1 = '.residuals.dat'
	name3 = '.save.pars'
	name4 = '.fit.chisq'
	
	### Get data and fit
	data  = numpy.loadtxt(path__+name1)
	save0 = data[:,0].astype(int)
	save4 = data[:,4]
	save5 = data[:,5]
	save7 = data[:,7]
	save9 = data[:,9]
	xxx__[save0]    = save4
	muu__[save0]    = save5
	xi_fit__[save0] = save7
	xi_err__[save0] = save9
	
	### Get only data
	data = numpy.loadtxt(pathData__+name2)
	save0 = data[:,0].astype(int)
	save1 = data[:,1]
	xi_dat__[save0] = save1
	
	### Get the parameters of the fit
	data    = numpy.loadtxt(path__+name3)
	param__[:,0] = data[:,1]
	param__[:,1] = data[:,2]
	
	### Get the parameters of the fit
	data = numpy.loadtxt(path__+name4)
	chi2__[:] = data
	
	return
def printResults():
	'''

	'''
	
	### Get the chi^2
	tmp_chi2    = chi2__[0]
	tmp_NBBin   = chi2__[1]
	tmp_NBParam = chi2__[2]
	print tmp_chi2
	
	### Get the fit parameters
	tmp_beta               = param__[0,0]
	tmp_beta_err           = param__[0,1]
	tmp_delta_v            = param__[4,0]
	tmp_delta_v_err        = param__[4,1]
	tmp_bias2              = param__[5,0]
	tmp_bias2_err          = param__[5,1]
	tmp_beta2_bias2        = param__[6,0]
	tmp_beta2_bias2_err    = param__[6,1]
	tmp_alpha_parallel     = param__[9,0]
	tmp_alpha_parallel_err = param__[9,1]
	tmp_alpha_perp         = param__[10,0]
	tmp_alpha_perp_err     = param__[10,1]
	
	string =  """
||  chiSquare / dof      || beta         ||  delta-v      ||  bias2        ||  beta2*bias2      ||  BAO alpha-parallel  ||  BAO alpha-perp  ||
||  %1.3e / (%u - %u)  ||  %1.3e ± %1.3e  ||  %1.3e ±  %1.3e   ||  %1.3e ± %1.3e    ||  %1.3e ± %1.3e       ||  %1.3e ± %1.3e       || %1.3e ± %1.3e  ||
""" % (tmp_chi2,tmp_NBBin,tmp_NBParam,tmp_beta,tmp_beta_err,tmp_delta_v,tmp_delta_v_err,tmp_bias2,tmp_bias2_err,tmp_beta2_bias2,tmp_beta2_bias2_err,tmp_alpha_parallel,tmp_alpha_parallel_err,tmp_alpha_perp,tmp_alpha_perp_err)
	print string

	return
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
		
		plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
		plt.grid()
		plt.xlabel(r'$'+b+' \, [h^{-1}.Mpc] $')
		plt.ylabel(r'$'+a+'\\xi('+b+') '+c+'$')
		plt.rc('font', **{'size':'30'})
		plt.tick_params(axis='both', which='major', labelsize=30)
		plt.grid(True, which='both')
		
		plt.show()
	
	return
def plotDataAndFit_2D():
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
		myTools.plot2D(coef*yyy_dat, edge, 's_{\\perp} \\, [h^{-1} Mpc]', 's_{\parallel} \\, [h^{-1} Mpc]', a+'\\xi(\\, \\overrightarrow{s} \\,)', 'data')
		### Fit
		myTools.plot2D(coef*yyy_fit, edge, 's_{\\perp} \\, [h^{-1} Mpc]', 's_{\parallel} \\, [h^{-1} Mpc]', a+'\\xi(\\, \\overrightarrow{s} \\,)', 'fit')
		### residuals
		myTools.plot2D(coef*(yyy_dat-yyy_fit)/yyy_err, edge, 's_{\\perp} \\, [h^{-1} Mpc]', 's_{\parallel} \\, [h^{-1} Mpc]', a+'\\xi_{data}-\\xi_{fit}', 'residuals')
		
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
	plt.rc('font', **{'size':'30'})
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.grid(True, which='both')
	plt.show()
	
	### Slice of constant parallel
	for i in range(0,100):
		tmpxxx = numpy.arange(0,50)*4.+2.
		#plt.errorbar(tmpxxx, yyy_dat[i,:], yerr=yyy_err[i,:], fmt='o')
		plt.errorbar(tmpxxx, yyy_fit[i,:])
	
	plt.xlabel(r'$s_{\perp} \, [Mpc.h^{-1}]$')
	plt.ylabel(r'$\xi(|s|)$')
	plt.rc('font', **{'size':'30'})
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.grid(True, which='both')
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
		
		plt.errorbar(sss, coef*yyy_dat[:,0], yerr=coef*yyy_dat[:,1], linestyle="", marker="o", color='blue')
		plt.errorbar(sss, coef*yyy_fit[:,0], color='red')
		
		plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
		plt.grid()
		plt.xlabel(r'$|s| \, [h^{-1} Mpc]$')
		plt.ylabel(r'$'+a+'\\xi(|s|)$')
		plt.rc('font', **{'size':'30'})
		plt.tick_params(axis='both', which='major', labelsize=30)
		plt.grid(True, which='both')
		
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
			
			plt.errorbar(sss, coef*yyy_dat[:,j,0], yerr=coef*yyy_dat[:,j,1], linestyle="", marker="o")
			plt.errorbar(sss, coef*yyy_fit[:,j,0], color='red')
			
			plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
			plt.grid()
			plt.xlabel(r'$|s| \, [h^{-1} Mpc]$')
			plt.ylabel(r'$'+a+'\\xi(|s|)$')
			plt.rc('font', **{'size':'30'})
			plt.tick_params(axis='both', which='major', labelsize=30)
			plt.grid(True, which='both')
			
		plt.show()
		
	return

def getHistoResiduals():
	'''
	
	'''
	
	### Constants
	nbBins = 100
	
	### Plot the residuals as a function of bin index
	xxx = numpy.arange(0., nbBin__ )
	yyy = xi_dat__-xi_fit__
	yyy[(xi_fit__==0.)] = 0.
	yer = xi_err__
	yer[(xi_fit__==0.)] = 0.
	plt.errorbar(xxx, yyy, yerr=yer, linestyle="", marker="o", color='blue')
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$\xi_{data}-\xi_{fit}$')
	plt.rc('font', **{'size':'30'})
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.grid(True, which='both')
	plt.show()
	### Plot the residuals/err as a function of bin index
	xxx = numpy.arange(0., nbBin__ )
	yyy = (xi_dat__-xi_fit__)/xi_err__
	yyy[(xi_fit__==0.)] = 0.
	plt.errorbar(xxx, yyy, linestyle="", marker="o", color='blue')
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.xlabel(r'$index \, bin$')
	plt.ylabel(r'$\frac{\xi_{data}-\xi_{fit}}{\sigma}$')
	plt.rc('font', **{'size':'30'})
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.grid(True, which='both')
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
	sizeX = 100
	sizeY = 100
	name = '.scan.dat'
	
	### Create a file in /tmp with the scan minus the two first lines
	idx = 0
	f     = open(path__+name)
	tmp_f = open('/tmp/scan.txt','w')
	for line in f:
		if (idx>=2): tmp_f.write(line)
		idx += 1
	tmp_f.close()
	f.close()
	
	### Data
	data = numpy.loadtxt('/tmp/scan.txt')
	
	### Best Fit
	alphaPara_bestFit = data[0][11]
	alphaPerp_bestFit = data[0][12]
	chi2_bestFit      = data[0][-1]
	print alphaPara_bestFit, alphaPerp_bestFit, chi2_bestFit
	
	### Scan
	alphaPara = data[1:,11]
	alphaPerp = data[1:,12]
	chi2      = data[1:,-1]
	
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
	
	edge = [0.5, 1.5, 0.5, 1.5]

	### Get best fit plus bootstrap
	nbBoot = 10000
	bestFit  = numpy.zeros( shape=(2,nbBoot+1) )
	chi2Boot = numpy.zeros( shape=(nbBoot+1) )
	bestFit[0,0] = alphaPerp_bestFit
	bestFit[1,0] = alphaPara_bestFit
	chi2Boot[0]  = chi2_bestFit
	for i in range(0,nbBoot):
		try:
			data     = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/BACKUP_2015_01_15/FitsFile_DR12_Guy/BaoFit_q_f/Bootstraps/boot_'+str(i).zfill(4)+'/bao2D.save.pars')
			dataChi2 = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/BACKUP_2015_01_15/FitsFile_DR12_Guy/BaoFit_q_f/Bootstraps/boot_'+str(i).zfill(4)+'/bao2D.fit.chisq')
			bestFit[0,i+1] = data[12,1]
			bestFit[1,i+1] = data[11,1]
			chi2Boot[i+1]  = dataChi2[0]
			#print '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/BaoFit_q_f/Bootstraps/boot_'+str(i).zfill(4)+'/bao2D.save.pars'
			#print bestFit[0,i+1], bestFit[1,i+1]
		except Exception:
			print i


	### Get the delta of chi^2
	chi2_2D[ (chi2_2D==0.) ] = float('nan')
	chi2_2D[ (chi2_2D!=0.) ] = chi2_2D[ (chi2_2D!=0.) ]-chi2_bestFit

	
	### Plot the data
	plotChi2Scan(chi2_2D,edge,True, bestFit,'\\alpha_{\\perp}','\\alpha_{\\parallel}','\\Delta \\chi^{2} = \\chi^{2}-\\chi^{2}_{best \\, fit}')

	###
	print numpy.mean(bestFit[0,1:][ bestFit[0,1:]>0. ])
	print numpy.mean(bestFit[1,1:])
	functionChi2 = interpolate.interp2d(alphaPerp_2D, alphaPara_2D, chi2_2D)
	chi2Realisation = numpy.asarray( [ functionChi2(bestFit[0,i],  bestFit[1,i]) for i in range(0,nbBoot+1) ] )

	### Get the Delta chi^2 value for fiducial model
	chi2Fiducial = functionChi2(1.,1.)
	print chi2Fiducial
	print chi2Realisation[ chi2Realisation >= chi2Fiducial ].size

	plt.hist(bestFit[0,1:][ bestFit[0,1:]>0. ], bins=20)
	plt.show()
	plt.hist(bestFit[1,1:][ bestFit[0,1:]>0. ], bins=20)
	plt.show()
	plt.hist(chi2Boot[1:][ bestFit[0,1:]>0. ]-chi2Boot[0], bins=20)
	plt.show()
	plt.hist(chi2Realisation[ bestFit[0,1:]>0. ], bins=20)
	plt.show()

	myTools.fitAGaussian(bestFit[0,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], False)
	myTools.fitAGaussian(bestFit[0,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], True)
	myTools.fitAGaussian(bestFit[1,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], False)
	myTools.fitAGaussian(bestFit[1,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], True)
	myTools.fitAGaussian(chi2Boot[1:][ bestFit[0,1:]>0.6 ]-chi2Boot[0], 50, [1,0.,50.], False)
	myTools.fitAGaussian(chi2Boot[1:][ bestFit[0,1:]>0.6 ]-chi2Boot[0], 50, [1,0.,50.], True)

	from scipy.stats import chisqprob
	### Get number under one sigma, two sigma ...
	nbParametrDiff = 2
	print '  parameter number difference  = ', nbParametrDiff
	tmp_a = chi2Realisation[1:]
	print tmp_a.size

	print 1.-chisqprob(2.3, 2)
	sigma = [2.30, 6.18, 11.83, 19.33, 28.74]
	
	for i in range(0,5):
		print ' Delta Chi^2 = ', sigma[i], '  sigma = ', i+1, ' nbReal = ', tmp_a[ tmp_a<sigma[i] ].size, '  expected = ', tmp_a.size*(1.-chisqprob(sigma[i], 2)), chisqprob(sigma[i], 2)

	return
def getChiScanToyMC():
	'''
		get the \chi^2 scan 
		1st line: 
		2nd line: settings
		3rd line: best fit
		4-end lines: scan
	'''
	
	### Constants
	sizeX = 100
	sizeY = 100
	name  = '.scan.dat'
	name2 = '.toymc.dat'
	
	### Create a file in /tmp with the scan minus the two first lines
	idx = 0
	f     = open(path__+name)
	tmp_f = open('/tmp/scan.txt','w')
	for line in f:
		if (idx>=2): tmp_f.write(line)
		idx += 1
	tmp_f.close()
	f.close()
	
	### Data
	data = numpy.loadtxt('/tmp/scan.txt')
	
	### Best Fit
	alphaPara_bestFit = data[0][11]
	alphaPerp_bestFit = data[0][12]
	chi2_bestFit      = data[0][-1]
	print alphaPara_bestFit, alphaPerp_bestFit, chi2_bestFit
	
	### Scan
	alphaPara = data[1:,11]
	alphaPerp = data[1:,12]
	chi2      = data[1:,-1]
	
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
	
	edge = [0.5, 1.5, 0.5, 1.5]

	### Get the delta of chi^2
	chi2_2D[ (chi2_2D==0.) ] = float('nan')
	chi2_2D[ (chi2_2D!=0.) ] = chi2_2D[ (chi2_2D!=0.) ]-chi2_bestFit

	### Create a file in /tmp with the scan minus the two first lines
	idx = 0
	f     = open('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/BACKUP_2015_01_15/FitsFile_DR12_Guy/BaoFit_q_f/toymc/bao2D.toymc.dat')
	tmp_f = open('/tmp/scan.txt','w')
	for line in f:
		if (idx>=2): tmp_f.write(line)
		idx += 1
	tmp_f.close()
	f.close()

	### Data
	data = numpy.loadtxt('/tmp/scan.txt')

	### Best fit and toy-MC
	data = data[1:]
	print ' nb toy model done = ',  data[:,0].size
	data = data[ data[:,-1]!=0 ]
	nbToyMC = data[:,0].size
	print '  nb toy model done correct = ', nbToyMC

	bestFit  = numpy.zeros( shape=(2,nbToyMC+1) )
	bestFit[0,0] = alphaPerp_bestFit
	bestFit[1,0] = alphaPara_bestFit

	bestFit[1,1:] = data[:,11]
	bestFit[0,1:] = data[:,12]


	
	### Plot the data
	plotChi2Scan(chi2_2D,edge,True, bestFit,'\\alpha_{\\perp}','\\alpha_{\\parallel}','\\Delta \\chi^{2} = \\chi^{2}-\\chi^{2}_{best \\, fit}')



	functionChi2 = interpolate.interp2d(alphaPerp_2D, alphaPara_2D, chi2_2D)
	chi2Realisation = numpy.asarray( [ functionChi2(bestFit[0,i],  bestFit[1,i]) for i in range(0,nbToyMC+1) ] )

	### Get the Delta chi^2 value for fiducial model
	chi2Fiducial = functionChi2(1.,1.)
	print chi2Fiducial
	print chi2Realisation[ chi2Realisation >= chi2Fiducial ].size

	plt.hist(bestFit[0,1:][ bestFit[0,1:]>0. ], bins=20)
	plt.show()
	plt.hist(bestFit[1,1:][ bestFit[0,1:]>0. ], bins=20)
	plt.show()
	#plt.hist(chi2Boot[1:][ bestFit[0,1:]>0. ]-chi2Boot[0], bins=20)
	#plt.show()
	plt.hist(chi2Realisation[ bestFit[0,1:]>0. ], bins=20)
	plt.show()

	#myTools.fitAGaussian(bestFit[0,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], False)
	#myTools.fitAGaussian(bestFit[0,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], True)
	#myTools.fitAGaussian(bestFit[1,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], False)
	#myTools.fitAGaussian(bestFit[1,1:][ bestFit[0,1:]>0.6 ], 50, [1,0.,1.], True)
	#myTools.fitAGaussian(chi2Boot[1:][ bestFit[0,1:]>0.6 ]-chi2Boot[0], 50, [1,0.,50.], False)
	#myTools.fitAGaussian(chi2Boot[1:][ bestFit[0,1:]>0.6 ]-chi2Boot[0], 50, [1,0.,50.], True)

	from scipy.stats import chisqprob
	### Get number under one sigma, two sigma ...
	nbParametrDiff = 2
	print '  parameter number difference  = ', nbParametrDiff
	tmp_a = chi2Realisation[1:]
	print tmp_a.size

	print 1.-chisqprob(2.3, 2)
	sigma = [2.30, 6.18, 11.83, 19.33, 28.74]
	
	for i in range(0,5):
		print ' Delta Chi^2 = ', sigma[i], '  sigma = ', i+1, ' nbReal = ', tmp_a[ tmp_a<sigma[i] ].size, '  expected = ', tmp_a.size*(1.-chisqprob(sigma[i], 2)), chisqprob(sigma[i], 2)

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
	idxAlpha = 11
	
	### Constants
	name = '.scan.dat'
	
	### Create a file in /tmp with the scan minus the two first lines
	idx = 0
	f     = open(path__+name)
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
	print data[0][9], data[0][10], param__[idxAlpha,0], param__[idxAlpha,1]
	
	### Scan
	alpha = data[1:,idxAlpha]
	chi2  = data[1:,-1]-chi2_bestFit

	### Chi^2 scan
	plt.errorbar(alpha,chi2,fmt='o')
	### Best fit
	plt.errorbar([alpha_bestFit],[0.],xerr=param__[idxAlpha,1],fmt='o',color='red',label='Best fit')
	### 1,2,3,4 sigma lines
	plt.plot(alpha,numpy.ones(alpha.size)*1.,color='green')
	plt.plot(alpha,numpy.ones(alpha.size)*4.,color='green')
	plt.plot(alpha,numpy.ones(alpha.size)*9.,color='green')
	plt.plot(alpha,numpy.ones(alpha.size)*16.,color='green')
	plt.xlabel(r'$'+paramName__[idxAlpha]+'$', fontsize=40)
	plt.ylabel(r'$\Delta \chi^{2} = \chi^{2}-\chi^{2}_{best \, fit}$', fontsize=40)
	plt.ylim([ -2., numpy.amax(chi2)+2. ])
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	return
def plotChi2Scan(data, edge=None, contour=False, bestFit=None, xTitle='-',yTitle='-',zTitle='-',title='-'):
	'''
	'''
	
	fig = plt.figure()
	ax = fig.add_subplot(111)

	myTools.deal_with_plot(False,False,False)
	ax.set_xticks([ 0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5 ])
	ax.set_yticks([ 0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5 ])
	
	plt.imshow(data, origin='lower',extent=edge, interpolation='None') ##, interpolation='None'
	cbar = plt.colorbar()
	plt.xlabel(r'$'+xTitle+'$', fontsize=40)
	plt.ylabel(r'$'+yTitle+'$', fontsize=40)
	cbar.set_label(r'$'+zTitle+'$',size=40)

	### Limit for color
	print numpy.amax(data)
	plt.clim(0.,63.69)
	
	plt.grid(True)
	#cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	
	### [1s = 2.30, 2s = 6.18, 3s = 11.83, 4s = 19.33]
	### http://www.reid.ai/2012/09/chi-squared-distribution-table-with.html
	if (contour):
		pylab.contour(data, [2.30, 6.18, 11.83, 19.33, 28.74], extent=edge, linewidths=3, colors = 'red', origin='lower', hold='on') #0.57,1.47, 
	
	if (bestFit[0,:].size > 1):
		for i in range(0,bestFit[0,1:].size):
			plt.scatter(bestFit[0,i+1],bestFit[1,i+1],marker='+',color='white',linewidths=5, hold='on')

	if (bestFit[0,:].size > 0):
		plt.scatter(bestFit[0,0],bestFit[1,0], marker='+',color='red',linewidths=10, hold='on')
		
	plt.xlim(edge[:2])
	plt.ylim(edge[2:])

	plt.show()
	
	return
'''
for i in range(0,10):
	for j in range(0,10):
		path__        = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f_covFromFit/bao' + type__ + ''
		pathData__    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f/bao' + type__
		getResults()
		getChiScan()
'''

getResults()
#getChiScan1D()
#getChiScanToyMC()
printResults()
#getChiScan()
plotDataAndFit_1D()
if (type__=='2D'):
	plotDataAndFit_2D()
	plotDataAndFit_2DSlice()
	plotDataAndFit_S()
	plotDataAndFit_Mu()
getHistoResiduals()
getChiScan()

