# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
# USE: python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/plot_histos.py
#


import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit
import myTools
from const_delta import *

chunckNb = 1
simulNb  = 1

forest__ = 'SIIV'
#lambdaRFMin__      = 1410.
#lambdaRFMax__      = 1530.

###
#path = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_histos/"
path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/"




### Delta+1 vs. lambda_Obs
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)

		data = numpy.loadtxt('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_histos/hDeltaVsLambdaObs_LYA.txt')
                plt.errorbar(data[:,0]+3600.5, data[:,1], label=r'$DR12$',color='red')

		data = numpy.loadtxt(path+'hDeltaVsLambdaObs_LYA'+mockNumber+'.txt')
		plt.errorbar(data[:,0]+3600.5, data[:,1], label=r'$eBOSS$',color='blue')

		data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/calibration_flux_using_CIV_forest.txt')
		plt.errorbar(data[:,0]+3547.5, data[:,1], label=r'$CIV \, forest \, backup$',color='red')

		data = numpy.loadtxt(path+'hDeltaVsLambdaObs_CIV'+mockNumber+'.txt')
		plt.errorbar(data[:,0]+3547.5, data[:,1], label=r'$CIV \, forest$',color='green')

		#data = numpy.loadtxt(path+'hDeltaVsLambdaObs_SIIV.txt')
                #plt.errorbar(data[:,0]+3547.5, data[:,1], label=r'$CIV \, forest$',color='black')

		#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/hDeltaVsLambdaObs_LYA_JMC.txt')
		#plt.errorbar(data[:,1][ data[:,2]!=0. ], data[:,2][ data[:,2]!=0. ], label=r'$simulation$',color='green')
		#data = numpy.loadtxt('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_histos_test_PDFMocksJMC/hDeltaVsLambdaObs_LYA.txt')
		#plt.errorbar(data[:,0]+3600.5, data[:,1], label=r'$DR12$',color='red')

		'''
		yMin    = numpy.amin(0.)
		yMax    = numpy.amax(1.)
		skyLines__ = numpy.power(10.,skyLines__)
		for i in range(0,skyLinesNames__.size):
			el = skyLines__[i]
			print numpy.mean(el)
			line = numpy.mean(el)
			xLi = [line,line]
			yLi = [yMin,yMax]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, skyLinesNames__[i], rotation='vertical', fontsize=20)
		'''

#plt.title(r'$hDeltaVsLambdaObs$', fontsize=40)
plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$flux$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()




### Template (m2)
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path + 'template_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$'+str(i)+' \, '+str(j)+'$')
		plt.plot([lambdaRFMin__,lambdaRFMin__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
		plt.plot([lambdaRFMax__,lambdaRFMax__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
		print path + 'template_'+forest__+mockNumber+'.txt'
		print data[:,0].size
		if (i==0 and j==0): saveYYY0 = data[:,1]
plt.title(r'$Template$', fontsize=40)
plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$Norm \, flux$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.show()

'''
### Template (m2)
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path + 'template_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0], data[:,1]-saveYYY0, fmt='o', label=r'$'+str(i)+' \, '+str(j)+'$')
		plt.plot([lambdaRFMin__,lambdaRFMin__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
		plt.plot([lambdaRFMax__,lambdaRFMax__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
plt.title(r'$Template$', fontsize=40)
plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$Norm \, flux$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()
'''

### Delta vs. lambda_RF
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'deltaVSLambdaRF_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$'+str(i)+' \, '+str(j)+'$')
		plt.plot([lambdaRFMin__,lambdaRFMin__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
		plt.plot([lambdaRFMax__,lambdaRFMax__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
plt.title(r'$Residual$', fontsize=40)
plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$< \delta >$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.show()


### Delta+1 vs. lambda_RF
data = numpy.loadtxt(path+'hDeltaVsLambdaRF_'+forest__+'.txt')
plt.errorbar(lambdaRFMin__-3.+data[:,0], data[:,1], fmt='o')
plt.plot([lambdaRFMin__,lambdaRFMin__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
plt.plot([lambdaRFMax__,lambdaRFMax__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
plt.title(r'$hDeltaVsLambdaRF$', fontsize=40)
plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$< \delta+1 >$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()


### Delta vs. lambda_Obs
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'deltaVSLambdaObs_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$'+str(i)+' \, '+str(j)+'$')
		print numpy.amin( data[:,0] )

		yMin    = numpy.min(data[:-1,1])
		yMax    = numpy.max(data[:-1,1])
		for i in range(0,skyLinesNames__.size):
			el = skyLines__[i]
			line = numpy.mean(el)
			xLi = [line,line]
			yLi = [yMin,yMax]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, skyLinesNames__[i], rotation='vertical', fontsize=20)

plt.title(r'$Residual$', fontsize=40)
plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$< \delta >$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

### Delta+1 vs. lambda_Obs
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'hDeltaVsLambdaObs_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0]+3600., data[:,1], marker='o', label=r'$'+str(i)+' \, '+str(j)+'$')

		yMin    = 0.
		yMax    = 1.
		for i in range(0,skyLinesNames__.size):
			el = skyLines__[i]
			line = numpy.mean(el)
			xLi = [line,line]
			yLi = [yMin,yMax]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, skyLinesNames__[i], rotation='vertical', fontsize=20)


#plt.title(r'$hDeltaVsLambdaObs$', fontsize=40)
plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$flux$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.show()

### Eta
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'eta_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$'+str(i)+' \, '+str(j)+'$')

plt.title(r'', fontsize=40)
plt.xlabel(r'$z_{pixel}$', fontsize=40)
plt.ylabel(r'$\eta$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.show()

### sigma^2_LSS

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range (0,chunckNb):
	for j in range(0,simulNb):
		mockNumber = '' #'_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'sigma2LSS_'+forest__+mockNumber+'.txt')
		ax.errorbar(data[:,0], data[:,1],yerr=data[:,2], fmt='o', label=r'$Data$')

		
		xxx = data[:,0]
		yyy = data[:,1]
		def chi2(a0,a1,a2):
			return numpy.sum(numpy.power( yyy-(a0+a1*numpy.power(xxx+1.,a2)) ,2.))
		m = Minuit(chi2,a0=1.,error_a0=1.,a1=1.,error_a1=1.,a2=1.,error_a2=1.,print_level=-1, errordef=0.01) 	
		m.migrad()
		a0 = m.values['a0']
		a1 = m.values['a1']
		a2 = m.values['a2']
		
		xxx = numpy.arange(0.,5.,0.01)
		print a0, a1, a2
		ax.plot( xxx, a0+a1*numpy.power(xxx+1.,a2),label=r'$\sigma_{L.S.S.}^{2}=a + b.(1+z)^{\gamma}$' )
		
plt.title(r'', fontsize=40)
plt.xlabel(r'$z_{pixel}$', fontsize=40)
plt.ylabel(r'$\sigma_{L.S.S.}^{2}$', fontsize=40)
myTools.deal_with_plot(False,False,True)


mng = plt.get_current_fig_manager()
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, r"$a=%.2e$" "\n" r"$b=%.2e$" "\n" r"$\gamma=%.2e$"%(a0,a1,a2), transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)

plt.show()




'''
	yyy = xi_resAll__[numpy.isfinite(xi_resAll__)]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.hist(yyy, bins=100)
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.rc('font', **{'size':'30'})
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.grid(True, which='both')
	mng = plt.get_current_fig_manager()
	textstr = '$Nb=%d$ \n $\mu=%.5e  +/-  %.2e$ \n $\sigma=%.5e  +/-  %.2e$'%(yyy.size,numpy.mean(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size),numpy.std(yyy),numpy.std(yyy)/numpy.sqrt(yyy.size))
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
	plt.show()
'''











