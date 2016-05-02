# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >


import astropy.io.fits as pyfits
import numpy
import matplotlib.pyplot as plt
import subprocess
import time
import iminuit as minuit
from iminuit import Minuit
import pylab

### Perso lib
import const_delta

#import warnings
#warnings.filterwarnings("error")

def Rebin(a, shape, weight=None):
	'''

		Function to rebin an array
		from Pierre Laurent ( pierre.laurent@cea.fr )
		/home/usr201/mnt/plaurent/libraries/mylib/mytools.py —> la fonction s’appelle rebin.

		a:     array
		shape: array of bin

	'''

	sh = [[shape[i],a.shape[i]//shape[i]] for i in numpy.arange(len(shape))]
	sh = [item  for sublist in sh for item in sublist]
	a = a.reshape(tuple(sh))

	if weight is None:
		return a.mean(axis=tuple(numpy.arange(1,len(sh),2)))
	else:
		weight = weight.reshape(tuple(sh))
		weight = numpy.where(weight==0,1e-37,weight)
		return numpy.where(numpy.average(weight,axis=tuple(numpy.arange(1,len(sh),2)))==0,0,(numpy.average(a,axis=tuple(numpy.arange(1,len(sh),2)),weights=weight)))
def GetHisto(data,nbBin,weights=None):
	''' 
		Get an histo 
	'''

	if (weights is None):
		hist, axisX = numpy.histogram(data,bins=nbBin)
	else:
		hist, axisX = numpy.histogram(data,bins=nbBin, weights=weights)
	xxx = numpy.array([ axisX[i]+(axisX[i+1]-axisX[i])/2. for i in range(0,axisX.size-1) ])

	return numpy.asarray(zip(xxx,hist))
def Get_TProfile(ar1, ar2, nbBin1, we2):
	'''
		Create a plot similar to a ROOT TProfile
	'''

	number, axisX, axisY = numpy.histogram2d(ar1, ar2, (nbBin1,1))
	weight, axisX, axisY = numpy.histogram2d(ar1, ar2, (nbBin1,1), weights=we2)
        mean,   axisX, axisY = numpy.histogram2d(ar1, ar2, (nbBin1,1), weights=we2*ar2)
        err,    axisX, axisY = numpy.histogram2d(ar1, ar2, (nbBin1,1), weights=we2*(ar2**2.))

	### find the axis X
        axisX = numpy.array([ axisX[i]+(axisX[i+1]-axisX[i])/2. for i in range(0,axisX.size-1) ])
	
	### Get only not empty bins
	bool_number = (number[:,0]>1)
	
	axisX  = axisX[bool_number]
	number = number[:,0][bool_number]
	weight = weight[:,0][bool_number]
	mean   = mean[:,0][bool_number]
	err    = err[:,0][bool_number]

	mean  = mean/weight
	err   = numpy.sqrt((err/weight-mean**2.)/number)

	return axisX, mean, err, number
def Get_2DTProfile(ar1, ar2, ar3, nbBinsX, nbBinsY,we):
	'''
	'''

	d = numpy.array(zip(ar1,ar2,ar3))
	number, axis = numpy.histogramdd( d, (nbBinsX,nbBinsY,1))
	weight, axis = numpy.histogramdd( d, (nbBinsX,nbBinsY,1), weights=we  )
	mean,   axis = numpy.histogramdd( d, (nbBinsX,nbBinsY,1), weights=we*ar3)
	err,    axis = numpy.histogramdd( d, (nbBinsX,nbBinsY,1), weights=we*(ar3**2.))

	mean   /= weight
	err    = numpy.sqrt((err/weight-mean**2.)/number)

	mean   = mean[:,:,0]
	err    = err[:,:,0]
	number = number[:,:,0]

	### find the axis X
	#axisX  = axis[0]
        #axisX = numpy.array([ axisX[i]+(axisX[i+1]-axisX[i])/2. for i in range(0,axisX.size-1) ])
	### find the axis Y
	#axisY  = axis[1]
        #axisY = numpy.array([ axisY[i]+(axisY[i+1]-axisY[i])/2. for i in range(0,axisY.size-1) ])

	### For test look at the histo
	#plt.imshow(mean,origin='lower',extent=[0., 10., 0., 10.],interpolation='None')
	#cbar = plt.colorbar()
	#plt.show()

	return mean, err, number
def deal_with_plot(logx=False, logy=False, legend=False):
	'''
		Usefull function for plots
	'''
	
	#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	#plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.rc('font', **{'size':'40'})
	plt.tick_params(axis='both', which='major', labelsize=40)
	#plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
	#plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)

	if (logx):
		plt.xscale('log')
	if (logy):
		plt.yscale('log')
	#plt.tight_layout()
	plt.grid(True, which='both')
	plt.linewidth = 40

	ax = plt.gca()
	[i.set_linewidth(3) for i in ax.spines.itervalues()]
	plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.12)

	if (legend):
		plt.legend(fontsize=40, numpoints=1,ncol=2)
	return
def isReadyForNewJobs(max_nbJobsMe, max_nbJobsAll,word='echo'):
	'''
		This function ends when there is room in the cluster for a job.
		

	'''

	commandQstatMe  = "qstat | grep "+ word +" | wc -l"
	commandQstatAll = "qstat -u '*' | wc -l"

	## Duration to wait before checking again
	waitTime = 60
	## Number of lines that aren't jobs
	noJobsLine = 2
	
	notDoProd = True
	while (notDoProd):
		
		## Counts the jobs from me
		nbJobsMe = int(subprocess.check_output(commandQstatMe, shell=True))
		
		## Counts jobs from everybody
		nbJobsAll = int(subprocess.check_output(commandQstatAll, shell=True))
		if (nbJobsAll>0):
			nbJobsAll = nbJobsAll-noJobsLine

		## Look if there are too many jobs
		if (nbJobsMe < max_nbJobsMe and nbJobsAll < max_nbJobsAll):
			notDoProd = False
		else:
			time.sleep(waitTime)

	## Return if there is room for a new job
	return
def plot1D(data, xTitle='-',yTitle='-',title='-',label=None):
	'''
	'''

	size = len(data)

	for i in numpy.arange(size):
		if (label is not None): plt.plot(numpy.arange(0,data[i].size),data[i],marker='o',label=label[i])
		else: plt.plot(numpy.arange(0,data[i].size),data[i],marker='o')

	plt.xlabel(r'$'+xTitle+'$', fontsize=40)
	plt.ylabel(r'$'+yTitle+'$', fontsize=40)

	if (label is not None): deal_with_plot(False,False,True)
	else: deal_with_plot()

	plt.show()

	return
def plot2D(data, edge=None, xTitle=None,yTitle=None,zTitle=None,title=None):
	'''
	'''
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	plt.imshow(data, origin='lower',extent=edge, interpolation='None')
	
	cbar = plt.colorbar()
	if(xTitle is not None): plt.xlabel(r'$'+xTitle+'$', fontsize=40)
	if(yTitle is not None): plt.ylabel(r'$'+yTitle+'$', fontsize=40)
	if(zTitle is not None): cbar.set_label(r'$'+zTitle+'$',size=40)
	if(title is not None): plt.title(r'$'+title+'$', fontsize=40)

	deal_with_plot()
	cbar.update_ticks()
	plt.show()
	
	return
def plot_2d_correlation(xi2D, x_power=0):


	origin='upper'
	#extent=[0.,200.,200., -200.]

	xxx = numpy.transpose(xi2D[:,:,0])
	yyy = numpy.transpose(xi2D[:,:,1])
	yer = numpy.transpose(xi2D[:,:,2])

	cut = (yer==0)
	if (xxx[cut].size==xxx.size):
		return
	yyy[ cut ] = float('nan')
	coef = numpy.power(xxx,x_power)

	fig = plt.figure()
	ax = fig.add_subplot(111)

	plt.imshow(coef*yyy, origin=origin, interpolation='None') #,extent=extent)
	cbar = plt.colorbar()

	plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	deal_with_plot(False,False,False)

	plt.show()
	
	return
def getCorrelationMatrix(cov):
	'''
		Get the correlation matrix from a covaraince matrix
	'''

	### Size
	size = cov[:,0].size

	### Get normalisation factor
	invSqrtDiag = 1./numpy.sqrt(numpy.diag(cov))

	### Normalize
	cor = numpy.array(cov)
	for i in range(0,size):
		cor[:,i] *= invSqrtDiag[i]
		cor[i,:] *= invSqrtDiag[i]
		cor[i,i] = 1.

	return cor
def getCovarianceMatrix(cor,diag):
	'''
		Get the covariance matrix from a correlation matrix and a diagonal
	'''

	### Size
	size = diag.size

	### Get normalisation factor
	sqrtDiag = numpy.sqrt(diag)

	### Normalize
	cov = numpy.array(cor)
	for i in range(0,size):
		cov[:,i] *= sqrtDiag[i]
	for i in range(0,size):
		cov[i,:] *= sqrtDiag[i]

	return cov
def convert1DTo2D(array1D,nbX,nbY):
	'''
		convert a 1D array to a 2D array
	'''

	array2D = numpy.zeros(shape=(nbX,nbY))

	for k in range(0,array1D.size):
		i = k/nbY
		j = k%nbY

		array2D[i][j] = array1D[k]

	return array2D
def convert2DTo1D(array2D,nbX,nbY):
	'''
		convert a 2D array to a 1D array
	'''

	array1D = numpy.zeros(nbX*nbY)

	for i in range(0,nbX):
		for j in range(0,nbY):
			k = i*nbY + j
			array1D[k] = array2D[i][j]

	return array1D
def fitAGaussian(data, bins, p0=[1,0.,1.], log=False):
	'''
		from http://stackoverflow.com/questions/11507028/fit-a-gaussian-function
	'''
	
	import numpy
	from scipy.optimize import curve_fit
	import matplotlib.pyplot as plt
	
	#data = numpy.random.normal(size=10000)
	
	hist, bin_edges = numpy.histogram(data, density=True, bins=bins)
	bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
	
	# Define model function to be used to fit to the data above:
	def gauss(x, *p):
		A, mu, sigma = p
		return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))
	
	# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
	#p0 = [1., 50., 50.]
	
	coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
	
	# Get the fitted curve
	minData = numpy.amin(data)
	maxData = numpy.amax(data)
	hist_fit = gauss(numpy.arange(minData,maxData,(maxData-minData)/1000.), *coeff)
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(bin_centres, hist, label=r'$Data$', drawstyle='steps')
	ax.plot(numpy.arange(minData,maxData,(maxData-minData)/1000.), hist_fit, label=r'$Fit$')
	
	# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
	print '\n'
	print 'Fitted mean = ', coeff[1], ' +/- ', numpy.sqrt(var_matrix[1][1])
	print 'Fitted standard deviation = ', coeff[2], ' +/- ', numpy.sqrt(var_matrix[2][2])
	print '\n'
	
	
	#plt.rc('font', **{'size':'30'})
	#ax.tick_params(axis='both', which='major', labelsize=30)
	#ax.grid(True, which='both')
	#ax.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	
	textstr = '$\mu=%.2f$\n$\sigma=%.2f$'%(coeff[1],coeff[2])
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)

	deal_with_plot(False, log, True)

	plt.show()
	
	return coeff, var_matrix
def fit_BAO(dist,obs,cov,arg_fit):
	'''

	from Pierre Laurent ( pierre.laurent@cea.fr )

	arg_fit = {}
        arg_fit['a'] = 125.08255109750635
        arg_fit['b'] = -1.8611671002360046
        arg_fit['c'] = 0.006109139533639488
        arg_fit['A'] = 0.003731291955288748
        arg_fit['sigma'] = 9.870124119190448
        arg_fit['mean'] = 106.73882608124278
        arg_fit['fix_A'] = False
        arg_fit['fix_sigma'] = True
        arg_fit['fix_mean'] = False
        arg_fit['background'] = 'inv'
        arg_fit['cov'] = 'full'
        arg_fit['x_min'] = 60.
        arg_fit['x_max'] = 160.

	'''

	cond_range = (dist > arg_fit['x_min']) & (dist < arg_fit['x_max'])
	x_fit = dist[cond_range]
	y_fit = obs[cond_range]
	
	
	if (arg_fit['background'] == 'poly2'):
		def BAOfunc(x,A,mean,sigma,a,b,c):
			return A*numpy.exp(-(x-mean)**2/(2*sigma**2)) + a*x**2 + b*x + c
	elif (arg_fit['background'] == 'inv'):
		def BAOfunc(x,A,mean,sigma,a,b,c):
			return A*numpy.exp(-(x-mean)**2/(2*sigma**2)) + a/(x**2) + b/x + c
	
		
	if (arg_fit['cov'] == 'var_only'):
		cov_fit = cov[cond_range]
		def chi2(A,mean,sigma,a,b,c):
			return numpy.sum(((BAOfunc(x_fit,A,mean,sigma,a,b,c) - y_fit)**2)/cov_fit)
	elif (arg_fit['cov'] == 'full'):
		cov_fit = cov[cond_range].T[cond_range].T
		def chi2(A,mean,sigma,a,b,c):
			return numpy.dot(numpy.dot((BAOfunc(x_fit,A,mean,sigma,a,b,c) - y_fit).T,numpy.linalg.inv(cov_fit)),(BAOfunc(x_fit,A,mean,sigma,a,b,c) - y_fit))
	
	
	
	m_BAO = minuit.Minuit(chi2,A=arg_fit['A'],mean=arg_fit['mean'],sigma=arg_fit['sigma'],a=arg_fit['a'],b=arg_fit['b'],c=arg_fit['c'],fix_A=arg_fit['fix_A'],fix_mean=arg_fit['fix_mean'],fix_sigma=arg_fit['fix_sigma'],print_level=0)
	m_BAO.migrad()
	m_BAO.hesse()
	
	m_noise = minuit.Minuit(chi2,A=0,mean=0,sigma=0,a=arg_fit['a'],b=arg_fit['b'],c=arg_fit['c'],fix_A=True,fix_mean=True,fix_sigma=True,print_level=0)
	m_noise.migrad()
	m_noise.hesse()
	
	x       = numpy.linspace(numpy.amin(dist),numpy.amax(dist),200)
	y_BAO   = BAOfunc(x,m_BAO.values['A'],m_BAO.values['mean'],m_BAO.values['sigma'],m_BAO.values['a'],m_BAO.values['b'],m_BAO.values['c'])
	y_noise = BAOfunc(x,m_noise.values['A'],m_noise.values['mean'],m_noise.values['sigma'],m_noise.values['a'],m_noise.values['b'],m_noise.values['c'])
	y_noBAO = BAOfunc(x,0.,m_BAO.values['mean'],m_BAO.values['sigma'],m_BAO.values['a'],m_BAO.values['b'],m_BAO.values['c'])
	
	print '  NDF : %i' % x_fit.size
	print '  chi2 for BAO + background : %.2f' % m_BAO.fval
	print '  chi2 for background only : %.2f' % m_noise.fval
	print '  delta chi2 = %.2f '  %(m_noise.fval- m_BAO.fval)
	print '  Amplitude = %f +/- %f' %(m_BAO.values['A'],m_BAO.errors['A'])
	print '  Mean = %f +/- %f' %(m_BAO.values['mean'],m_BAO.errors['mean'])
	print '  Sigma = %f +/- %f' %(m_BAO.values['sigma'],m_BAO.errors['sigma'])
	
	
	return m_BAO.values, m_BAO.errors, [x,y_BAO], m_noise.values, m_noise.errors, [x,y_noise], [x,y_noBAO]
def plotCovar(covList, covNameList,nbBinX=50,nbBinY=100,correlation='q_f'):

	color = ['blue','red','green','orange','black']

	### Diagrams
	###################################################
	### 1D
	min1D__   = 0.
	max1D__   = 200.
	nbBin1D__ = nbBinX
	binSize__ = max1D__/nbBin1D__

	### 2D
	minX2D__   = min1D__
	maxX2D__   = max1D__
	nbBinX2D__ = nbBin1D__

	minY2D__   = -max1D__
	if (correlation=='q_q' or correlation=='f_f'):
		minY2D__ = 0.
	maxY2D__   = max1D__
	nbBinY2D__ = nbBinY

	nbBin2D__  = nbBinX2D__*nbBinY2D__


	nbCov       = len(covList)

	nbBin = covList[0][0,:].size
	print '  nbBin = ', nbBin
	print '  nbBin2D__ = ', nbBin2D__
	
	corList      = [numpy.zeros(shape=(nbBin,nbBin))]*nbCov
	diagList     = [numpy.diagonal(covList[i]) for i in range(0,nbCov)]
	diagSqrtList = [numpy.sqrt(diagList[i]) for i in range(0,nbCov)]

	for i in range(0,nbCov):
		for j in range(0,nbBin):
			for k in range(0,nbBin):
				corList[i][j][k] = covList[i][j][k]/(diagSqrtList[i][j]*diagSqrtList[i][k])
	
		corList[i] = numpy.nan_to_num(corList[i])
	
	### Correlation
	###################################################

	### Show the variogram
	tmp_xxx  = []
	tmp_yyy  = []
	for c in range(0,nbCov):
		if (nbBin==nbBin1D__ or nbBin==161):
			xxx = binSize__*numpy.arange(0,nbBin)
			yyy = numpy.zeros(shape=(1,nbBin,2))
			for k1 in range(0,nbBin):
				for k2 in range(0,k1+1):
					if (corList[c][k1][k2]!=0.):
						yyy[0][k1-k2][0] += corList[c][k1][k2]
						yyy[0][k1-k2][1] += 1.
		else:
			xxx  = binSize__*numpy.arange(0, nbBinY2D__)
			xxx2 = binSize__*numpy.arange(0, nbBinX2D__)
			yyy  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__,2))
			yyy2 = numpy.zeros(shape=(nbBinY2D__,nbBinX2D__,2))
			for k1 in range(0,nbBin):
				i1 = k1/nbBinY2D__
				j1 = k1%nbBinY2D__
				for k2 in range(0,k1+1):
					i2 = k2/nbBinY2D__
					j2 = k2%nbBinY2D__
					if (corList[c][k1][k2]==0.): continue
					
					yyy[abs(i1-i2)][abs(j1-j2)][0]  += corList[c][k1][k2]
					yyy[abs(i1-i2)][abs(j1-j2)][1]  += 1.
					yyy2[abs(j1-j2)][abs(i1-i2)][0] += corList[c][k1][k2]
					yyy2[abs(j1-j2)][abs(i1-i2)][1] += 1.

		
		yyy  = numpy.nan_to_num(yyy[:,:,0]/yyy[:,:,1])
		
		### 0
		tmp_xxx += [xxx]
		tmp_yyy += [yyy]

	parameterFromFit = numpy.zeros(shape=(nbCov,nbBinX2D__,3))

	### Show all the method on the main variogram
	for i in range(0,nbBinX2D__):
		if (i>=1 and (nbBin==nbBin1D__ or nbBin==161 or nbBin==445) ): continue
		
		for c in range(0,nbCov):
			'''
			### fir the covariance by a model
			def chi2(a0,a1,a2):
				model = (a0 + a1*tmp_xxx[c])/(1. + a2*tmp_xxx[c])
				return numpy.sum(numpy.power(tmp_yyy[c][i]-model,2.))

			if (i==0):
				m = Minuit(chi2, a0=1.,error_a0=1., a1=1.,error_a1=1.,  a2=1.,error_a2=1.,fix_a0=True,  print_level=-1, pedantic=False)
			elif (i<=13):
				m = Minuit(chi2, a0=1.,error_a0=1., a1=1.,error_a1=1.,  a2=1.,error_a2=1., print_level=-1, pedantic=False)
			else:
				m = Minuit(chi2, a0=1.,error_a0=1., a1=0.,error_a1=1.,  a2=0.,error_a2=1.,fix_a1=True, fix_a2=True, print_level=-1, pedantic=False)
			m.migrad()

			parameterFromFit[c][i][0] = m.values['a0']
			parameterFromFit[c][i][1] = m.values['a1']
			parameterFromFit[c][i][2] = m.values['a2']
		
			minuit_yyy = (parameterFromFit[c][i][0] + parameterFromFit[c][i][1]*tmp_xxx[c])/(1. + parameterFromFit[c][i][2]*tmp_xxx[c])
			'''
		
			### Plot the results
			if (nbBin==nbBin1D__):
				xTitle = '|\Delta \, |s|_{1}-|s|_{2} | \, [h^{-1}.Mpc]'
			else:
				xTitle = '|\Delta \, s_{\parallel,1}-s_{\parallel,2} | \, [h^{-1}.Mpc], \, for \, |\Delta \, s_{\perp}| = ' + str(i*binSize__)
			plt.errorbar(tmp_xxx[c], tmp_yyy[c][i], label=r'$'+covNameList[c]+'$', marker='o', color=color[c])
			#plt.plot(tmp_xxx[c], minuit_yyy, label=r'$model$', marker='o')

			### Save the results of the fit
			

		plt.xlabel(r'$'+xTitle+'$', fontsize=40)
		plt.ylabel(r'$Cor( s_{1}, s_{2})$', fontsize=40)
		plt.xlim([ numpy.min(tmp_xxx[0])-10., numpy.max(tmp_xxx[0])+10. ])
		deal_with_plot(False,False,True)
		plt.show()
		

	

	### Get the covariance from fit
	covFromFitList = [numpy.zeros(shape=(nbBin,nbBin))]*nbCov
	corFromFitList = [numpy.zeros(shape=(nbBin,nbBin))]*nbCov
	for c in range(0,nbCov):

		meanBasment = numpy.mean( tmp_yyy[c][6:,:])

		if (nbBin==nbBin1D__):
			for k1 in range(0,nbBin):
				for k2 in range(0,k1+1):
					if (corList[c][k1][k2]==0.): continue
					x         = numpy.abs(k1-k2)*binSize__
					#y_fromFit = (parameterFromFit[c][0][0] + parameterFromFit[c][0][1]*x)/(1. + parameterFromFit[c][0][2]*x)
					y_fromFit = tmp_yyy[c][0][abs(k1-k2)]
					covFromFitList[c][k1][k2] = y_fromFit*diagSqrtList[c][k1]*diagSqrtList[c][k2]
					corFromFitList[c][k1][k2] = y_fromFit
					covFromFitList[c][k2][k1] = covFromFitList[c][k1][k2]
					corFromFitList[c][k2][k1] = y_fromFit
		else:
			for k1 in range(0,nbBin):
				i1 = k1/nbBinY2D__
				j1 = k1%nbBinY2D__
				for k2 in range(0,k2+1):
					if (corList[c][k1][k2]==0.): continue
					i2 = k2/nbBinY2D__
					j2 = k2%nbBinY2D__

					x1 = numpy.abs(j1-j2)*binSize__
					x2 = numpy.abs(i1-i2)

					#y_fromFit = (parameterFromFit[c][x2][0] + parameterFromFit[c][x2][1]*x1)/(1. + parameterFromFit[c][x2][2]*x1)
					if (abs(i1-i2)<=5): y_fromFit = tmp_yyy[c][abs(i1-i2)][abs(j1-j2)]
					else: y_fromFit = meanBasment
					covFromFitList[c][k1][k2] = y_fromFit*diagSqrtList[c][k1]*diagSqrtList[c][k2]
					corFromFitList[c][k1][k2] = y_fromFit
					covFromFitList[c][k2][k1] = covFromFitList[c][k1][k2]
					corFromFitList[c][k2][k1] = y_fromFit
	return covFromFitList
def plotOnSpectra_plate(plate,mjd,fiber):
	'''
		Plot a spectra when given in a sp_plate formate

	'''
	
	### Constants
	pathToSpec = '/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_9_1/'
	pathToCat  = '/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_9_1.fits'
	sizeBitsQSO = const_delta.bitsQSO__['NAME'].size
	z = 0. 

	'''
	### Get its redshift
	cat = pyfits.open(pathToCat, memmap=True)[1].data
	cat = cat[ cat['PLATE'] ==  plate]
	cat = cat[ cat['MJD'] == mjd ]
	cat = cat[ cat['FIBERID'] == fiber ]
	cat = cat[0]
	z = cat['Z']
	print '  myTools::plotOnSpectra_plate:: CLASS = ', cat['CLASS']
	print '  myTools::plotOnSpectra_plate:: ZWARNING = ', cat['ZWARNING']
	print '  myTools::plotOnSpectra_plate:: Z_ERR = ', cat['Z_ERR']
	print '  myTools::plotOnSpectra_plate:: BOSS_TARGET1 = ', cat['BOSS_TARGET1']
	print '  myTools::plotOnSpectra_plate:: ANCILLARY_TARGET1 = ', cat['ANCILLARY_TARGET1']
	print '  myTools::plotOnSpectra_plate:: ANCILLARY_TARGET2 = ', cat['ANCILLARY_TARGET2']
	print '  myTools::plotOnSpectra_plate:: EBOSS_TARGET0 = ', cat['EBOSS_TARGET0']
	print '  myTools::plotOnSpectra_plate:: EBOSS_TARGET1 = ', cat['EBOSS_TARGET1']
	print '  myTools::plotOnSpectra_plate:: EBOSS_TARGET2 = ', cat['EBOSS_TARGET2']
	print '  myTools::plotOnSpectra_plate:: BOSS_TARGET1 = ', [ const_delta.bitsQSO__['NAME'][i] for i in range(0,sizeBitsQSO) if ( cat['BOSS_TARGET1'] & 2**const_delta.bitsQSO__['BIT'][i] >0 and const_delta.bitsQSO__['KEY'][i]=='BOSS_TARGET1') ]
	print '  myTools::plotOnSpectra_plate:: ANCILLARY_TARGET1 = ', [ const_delta.bitsQSO__['NAME'][i] for i in range(0,sizeBitsQSO) if ( cat['ANCILLARY_TARGET1'] & 2**const_delta.bitsQSO__['BIT'][i] >0 and const_delta.bitsQSO__['KEY'][i]=='ANCILLARY_TARGET1') ]
	print '  myTools::plotOnSpectra_plate:: ANCILLARY_TARGET2 = ', [ const_delta.bitsQSO__['NAME'][i] for i in range(0,sizeBitsQSO) if ( cat['ANCILLARY_TARGET2'] & 2**const_delta.bitsQSO__['BIT'][i] >0 and const_delta.bitsQSO__['KEY'][i]=='ANCILLARY_TARGET2') ]
	print '  myTools::plotOnSpectra_plate:: EBOSS_TARGET0 = ', [ const_delta.bitsQSO__['NAME'][i] for i in range(0,sizeBitsQSO) if ( cat['EBOSS_TARGET0'] & 2**const_delta.bitsQSO__['BIT'][i] >0 and const_delta.bitsQSO__['KEY'][i]=='EBOSS_TARGET0') ]
	print '  myTools::plotOnSpectra_plate:: EBOSS_TARGET1 = ', [ const_delta.bitsQSO__['NAME'][i] for i in range(0,sizeBitsQSO) if ( cat['EBOSS_TARGET1'] & 2**const_delta.bitsQSO__['BIT'][i] >0 and const_delta.bitsQSO__['KEY'][i]=='EBOSS_TARGET1') ]
	print '  myTools::plotOnSpectra_plate:: EBOSS_TARGET2 = ', [ const_delta.bitsQSO__['NAME'][i] for i in range(0,sizeBitsQSO) if ( cat['EBOSS_TARGET2'] & 2**const_delta.bitsQSO__['BIT'][i] >0 and const_delta.bitsQSO__['KEY'][i]=='EBOSS_TARGET2') ]
	'''

	
	### Get its redshift from reduced cat
	data, plate_list = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Run/list_LYA.npy')
	z = data[:,5][ numpy.logical_and( numpy.logical_and( data[:,0]==plate, data[:,1]==mjd ), data[:,2]==fiber) ][0]
	print '  myTools::plotOnSpectra_plate:: z = ', z
	

	### Get spectrum
	tmpCat = pyfits.open(pathToSpec + 'spPlate-' + str(plate) + '-' + str(mjd) + '.fits', memmap=True)

	### Get attributes
	nbPixels         = tmpCat[0].header['NAXIS1']
	stepPixels       = tmpCat[0].header['CD1_1']
	logLamFirstPixel = tmpCat[0].header['CRVAL1']

	### Create an objects with data
	dtype_out=[('LOGLAM','<f8'),('FLUX','<f8'),('IVAR','<f8'),('AND_MASK',numpy.int16)]
  	cat = numpy.ndarray(shape=nbPixels,dtype=dtype_out)
	cat['LOGLAM']   = numpy.power(10.,numpy.arange(0.,nbPixels)*stepPixels+logLamFirstPixel)
	cat['FLUX']     = tmpCat[0].data[fiber-1]
	cat['IVAR']     = tmpCat[1].data[fiber-1]
	cat['AND_MASK'] = tmpCat[2].data[fiber-1]

	### Plot
	plt.plot( cat['LOGLAM']/(1.+z), cat['IVAR'],label='ivar', color='red')
	plt.plot( cat['LOGLAM']/(1.+z), cat['FLUX'],label='flux',color='blue' )
	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$flux$', fontsize=40)
	deal_with_plot(False,False,True)
	plt.show()
		

	return
def plotOnSpectra_spec(plate,mjd,fiber,z=0.,template=None):
	'''
		Plot a spectra when given in a sp_spec formate

	'''

	from const_delta import *
	import scipy
	lambdaObsMin__ = 3555.
	lambdaObsMax__ = 10290.
	log10lambdaObsMin__ = numpy.log10(lambdaObsMin__)
	log10lambdaObsMax__ = numpy.log10(lambdaObsMax__)

	
	### Constants
	pathToSpec = '/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_8_guy/spectra/'
	pathToCat  = '/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits'
	sizeBitsQSO = const_delta.bitsQSO__['NAME'].size

	if (z!=0.):
		cat = pyfits.open(pathToCat, memmap=True)[1].data
		cat = cat[ numpy.logical_and(numpy.logical_and( cat['PLATE'] ==  plate, cat['MJD'] == mjd, ), cat['FIBERID'] == fiber) ]	
		cat = cat[0]
		z = cat['Z_VI']
		print z

	### Get spectrum
	cat = pyfits.open(pathToSpec + str(plate) + "/spec-" + str(plate) + "-" + str(mjd) + "-" + str(fiber).zfill(4) + ".fits", memmap=True)[1].data

	### Get dispertion \Delta \log(\lambda)
	#difLOGLOBS = numpy.array( [ cat['LOGLAM'][i+1]-cat['LOGLAM'][i] for i in numpy.arange(cat['LOGLAM'].size-1) if ( cat['LOGLAM'][i]!=0. and cat['LOGLAM'][i+1]!=0.) ] )
	#print difLOGLOBS
	#plt.hist(difLOGLOBS)
	#plt.show()

	cat = cat[ (numpy.logical_and( (cat["LOGLAM"]>=log10lambdaObsMin__) , (cat["LOGLAM"]<log10lambdaObsMax__) )) ]
	#cat = cat[ cat["LOGLAM"]>=log10lambdaObsMin__ ]
	for lines in skyLines__:
		cat = cat[ (numpy.logical_or( (cat["LOGLAM"]<=lines[0]) , (cat["LOGLAM"]>=lines[1]) )) ]
		cat = cat[ numpy.logical_and( numpy.logical_and( (cat["IVAR"]>0.), (cat["AND_MASK"]<bit16__)), (numpy.isfinite(cat["FLUX"])) ) ]
	removeSkyLines = scipy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/Calibration/calibration_flux_using_CIV_forest.txt')
	#removeSkyLines = scipy.interpolate.interp1d(numpy.log10(3447.5+removeSkyLines[:,0]),removeSkyLines[:,1],bounds_error=False,fill_value=1)

	#plt.plot(3447.5+removeSkyLines[:,0], removeSkyLines[:,1],color='red')
	#deal_with_plot(False,False,False)
	#plt.show()
	removeSkyLines = scipy.interpolate.interp1d(numpy.log10(3447.5+removeSkyLines[:,0]),removeSkyLines[:,1],bounds_error=False,fill_value=1)

	coef = removeSkyLines(cat["LOGLAM"])
	cat['FLUX'] /= coef


	tmp_logZ = numpy.log10(1.+z)
	normFactor = numpy.mean( cat["FLUX"][ numpy.logical_and( (cat["LOGLAM"]>log10lambdaRFNormaMin__+tmp_logZ), (cat["LOGLAM"]<log10lambdaRFNormaMax__+tmp_logZ) ) ] )
	cat['FLUX'] /= normFactor

	cat['LOGLAM'] = numpy.power(10.,cat['LOGLAM'])#/(1.+z)
	

	#print cat['LOGLAM'][ numpy.logical_and( (cat['LOGLAM']>=const_delta.lambdaRFTemplateMin__) , (cat['LOGLAM']<const_delta.lambdaRFTemplateMax__) ) ].size

	### Plot
	#plt.plot( cat['LOGLAM'], cat['IVAR'],label='ivar', color='red')
	plt.plot( cat['LOGLAM'], cat['FLUX'],label=r'$Data$',color='blue', markersize=8,linewidth=2)
	#if (template is not None): plt.plot( template[0]*(1.+z),template[1],label=r'$Quasar \, continuum$',color='red', markersize=8,linewidth=4)

	
	line = 1215.67*(1.+z)
	yLi = [numpy.min(cat['FLUX'])*1.1,2.*numpy.max(cat['FLUX'])]
	line = 1909*(1.+z)
	plt.plot([line,line],yLi,color='yellow',markersize=8,linewidth=4,label=r'$C-III$')
	line = 1548.2049*(1.+z)
	plt.plot([line,line],yLi,color='orange',markersize=8,linewidth=4,label=r'$C-IV$')
	line = 1393.76018*(1.+z)
        plt.plot([line,line],yLi,color='violet',markersize=8,linewidth=4,label=r'$Si-IV$')
	line = 1215.67*(1.+z)
	plt.plot([line,line],yLi,color='green',markersize=8,linewidth=4,label=r'$Lyman-\alpha$')
	line = 1025.72*(1.+z)
        plt.plot([line,line],yLi,color='pink',markersize=8,linewidth=4,label=r'$Lyman-\beta$')
	

	#plt.text(line, numpy.max(cat['FLUX'])*1.1, r'$Lyman-\alpha$', rotation='vertical', fontsize=30)

	if (template is not None): plt.plot( template[0]*(1.+z),template[1],label=r'$Quasar \, continuum$',color='red', markersize=8,linewidth=4)
	plt.plot([1040.*(1.+z),1200.*(1.+z)],[0.,0.],color='red',markersize=8,linewidth=5,label=r'$Lyman-\alpha forest$')

	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$\phi(\lambda_{Obs.})$', fontsize=40)
	deal_with_plot(False,False,True)
	plt.ylim([ numpy.min(cat['FLUX'])*1.01,numpy.max(cat['FLUX'])*1.15  ])
	plt.legend(fontsize=40, numpoints=1,ncol=1)
	plt.show()
		

	return
def plot_chi2_scan(data, edge=None, contour=False, bestFit=None, xTitle='-',yTitle='-',zTitle='-',title='-'):
	'''
	'''
	
	fig = plt.figure()
	ax = fig.add_subplot(111)

	deal_with_plot(False,False,False)
	plt.ticklabel_format(axis='x',useOffset=False)
	plt.ticklabel_format(axis='y',useOffset=False)
	plt.ticklabel_format(axis='z',useOffset=False)
	ax.set_xticks([ 0.98,0.98,1.0,1.02,1.02 ])
	ax.set_yticks([ 0.98,0.98,1.0,1.02,1.02 ])
	
	plt.imshow(data, origin='lower',extent=edge)#, interpolation='None')
	cbar = plt.colorbar()
	plt.xlabel(r'$'+xTitle+'$', fontsize=40)
	plt.ylabel(r'$'+yTitle+'$', fontsize=40)
	cbar.set_label(r'$'+zTitle+'$',size=40)

	### Limit for color
	#plt.clim(0.,63.69)
	#plt.clim(0.,6.18)
	
	plt.grid(True)
	#cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	
	### [1s = 2.30, 2s = 6.18, 3s = 11.83, 4s = 19.33]
	### http://www.reid.ai/2012/09/chi-squared-distribution-table-with.html
	if (contour):
		pylab.contour(data, [2.30, 6.18, 11.83, 19.33, 28.74], extent=edge, linewidths=3, colors = 'red', origin='lower') #0.57,1.47, 
	
	if (bestFit is not None):
		if (bestFit[0,:].size > 1):
			x = bestFit[0,1:][ bestFit[0,1:]!=0. ]
			y = bestFit[1,1:][ bestFit[1,1:]!=0. ]
			plt.scatter(x,y,marker='+',color='white',linewidths=5, hold='on') #,alpha=0.3)
		if (bestFit[0,:].size > 0):
			plt.scatter(bestFit[0,0],bestFit[1,0], marker='+',color='red',linewidths=10, hold='on')

	### Show center
	plt.scatter([1.],[1.], marker='o',color='green',linewidths=10, hold='on')

	plt.xlim(edge[:2])
	plt.ylim(edge[2:])

	#plt.xlim([0.7,1.4])
	#plt.ylim([0.7,1.4])

	plt.show()
	
	return
def append_files(all_t, path_to_save):
	### Get the FitsFile
	all_t_file  = []
	for el in all_t:
		all_t_file.append( pyfits.open(el,  memmap=True) )
		print el
	
	### Get arrays of size of FitsFile
	all_t_nrows = []
	for el in all_t_file:
		print el[1].data.size
		all_t_nrows.append( el[1].data.shape[0] )
	nrowsTot = numpy.sum(all_t_nrows)
	all_t_nrows = numpy.array(all_t_nrows)

	print 
	print '  ', all_t_nrows.size
	print '  ', nrowsTot

	### Set the Fits_File which will contain all
	hdu = pyfits.BinTableHDU.from_columns(all_t_file[0][1].columns, nrows=nrowsTot)

	### Set the values of each rows
	first = 0
	last  = all_t_nrows[0]
	for i in range(1, all_t_nrows.size ):
		first += all_t_nrows[i-1]
		last  += all_t_nrows[i]
		print first, last
		for colname in all_t_file[0][1].columns.names:
			hdu.data[colname][first:last] = all_t_file[i][1].data[colname]

	### Add some column to the data
	cat     = hdu.data
	sizeMax = cat.size

	print sizeMax

	tbhdu = pyfits.BinTableHDU(data=cat)
	tbhdu.update()
	tbhdu.writeto(path_to_save, clobber=True)

	return
def save_fited_cor_from_cov(path_to_load, path_to_save,nbBinX=50,nbBinY=100,correlation='q_f'):
	"""

	save a fitted version of a correlation matrix from a covariance matrix

	"""

	### Load covariance
	cov = numpy.load(path_to_load)
	plot2D(cov)

	### Ge the correlation matrix
	cor = getCorrelationMatrix(cov)
	plot2D(cor)

	### model the mean and fit
	mean_cor = plotCovar([cor], ['a'], nbBinX,nbBinY,correlation)[0]
	plot2D(mean_cor)

	### Save the mean fitted
	numpy.save(path_to_save, mean_cor)

	### Look at the shape of the mean fitted reduced correlation
	plotCovar([mean_cor], ['a'], nbBinX,nbBinY,correlation)[0]

	return



























