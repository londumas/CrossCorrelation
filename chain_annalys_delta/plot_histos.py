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
from scipy import interpolate

chunckNb = 1
simulNb  = 1
mockNumber = ''
isMock_ = False
forest__ = 'LYA'
if (forest__ == 'LYB'):
	lambdaRFMin__      = 800.
	lambdaRFMax__      = 1020.
	shift__            = 3547.5
elif (forest__ == 'LYA'):
	lambdaRFMin__      = 1040.
	lambdaRFMax__      = 1200.
	shift__            = 3500. #3447.5
	
elif (forest__ == 'SIIV'):
	lambdaRFMin__      = 1286.
	lambdaRFMax__      = 1380.
	shift__            = 3447.5
elif (forest__ == 'CIV'):
	lambdaRFMin__      = 1410.
	lambdaRFMax__      = 1530.
	shift__            = 3447.5
elif (forest__ == 'MGII'):
	lambdaRFMin__      = 2100.
	lambdaRFMax__      = 2790.
	shift__            = 3547.5
###
#path = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/FitsFile_DR12_reOBS_eBOSS_Guy/DR12_primery/histos/"
#path = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/FitsFile_DR12_reOBS_Guy/DR12_primery/histos/"
path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest__+'/FitsFile_DR12_Guy_Margala/DR12_primery/histos/' ##_method1
#path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/'+forest__+'/FitsFile_DR14/DR14_primery/histos/'
#path = "/home/gpfs/manip/mnt/bao/hdumasde/Data/MGII/FitsFile_DR12_Guy/DR12_primery/histos/"
#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/'+forest__+'/DR14/DR14_primery/histos/'
#path = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA//FitsFile_eBOSS_Guy/all_eBOSS_primery/histos/"
rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/'


### Delta+1 vs. lambda_Obs
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		if (isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)


		data = numpy.loadtxt(path+'hDeltaVsLambdaObs_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0]+shift__, data[:,1], label=r'$DR14$')
		tmp_data = data[:,1]

		data = numpy.loadtxt('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_Margala/DR12_primery/histos/hDeltaVsLambdaObs_LYA.txt')
                plt.errorbar(data[:,0]+shift__, data[:,1], label=r'$DR12$')
		yyy = data[:,1]
		lll = data[:,0]+shift__
		zzz = lll/1215.67-1.

		data = numpy.loadtxt('/home/gpfs/manip/mnt/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy_Margala/DR12_primery/histos/hDeltaVsLambdaObs_CIV.txt')
		plt.errorbar(data[:,0]+shift__, data[:,1], label=r'$DR12: CIV \, forest \, for \, calibration $')

		#plt.errorbar(data[:,0]+shift__, data[:,1]/tmp_data-1., label=r'$devided$')

		#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/PDF/hDeltaVsLambdaObs_LYA_JMC.txt')
		#plt.errorbar(data[:,1][ data[:,2]!=0. ], data[:,2][ data[:,2]!=0. ], label=r'$v1547$', color='red')
		#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/v1575_hDeltaVsLambdaObs_LYA_JMC.txt')
		#plt.errorbar(data[:,1][ data[:,2]!=0. ], data[:,2][ data[:,2]!=0. ], label=r'$v1575$', color='orange')

		
		def chi2(a0,a1):
			return numpy.sum( numpy.power( yyy-numpy.exp(-a0*numpy.power(1.+zzz,a1) ) ,2.) )
		m = Minuit(chi2,a0=0.0023,error_a0=1.,a1=3.65,error_a1=1.,print_level=-1, errordef=0.01) 	
		m.migrad()
		a0 = m.values['a0']
		a1 = m.values['a1']
		print '  parameter = ', a0, a1
		plt.plot(lll, numpy.exp(-a0*numpy.power(1.+zzz,a1) ) )


		'''
		yMin    = numpy.amin(0.)
		yMax    = numpy.amax(1.)
		skyLines__ = numpy.power(10.,skyLines__)
		for k in range(skyLinesNames__.size):
			el = skyLines__[k]
			print numpy.mean(el)
			line = numpy.mean(el)
			xLi = [line,line]
			yLi = [yMin,yMax]
			plt.plot(xLi,yLi,color='green')
			#plt.text(line, yMax, skyLinesNames__[k], rotation='vertical', fontsize=20)
		'''

#plt.title(r'$hDeltaVsLambdaObs$', fontsize=40)
plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$Normalized \, flux$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()




### Template (m2)
for i in range (0,chunckNb):
	for j in range(0,simulNb):

		#data = numpy.loadtxt('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/histos/template_LYA.txt')
		#plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$DR12$')

		if (isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path + 'template_'+forest__+mockNumber+'.txt')

		plt.errorbar(data[:,0], data[:,1], marker='o') #, label=r'$'+str(i)+' \, '+str(j)+'$')
		plt.plot([lambdaRFMin__,lambdaRFMin__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red', label=r'$forest \, definition$')
		plt.plot([lambdaRFMax__,lambdaRFMax__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='red')
		print path + 'template_'+forest__+mockNumber+'.txt'
		print data[:,0].size
		if (i==0 and j==0): saveYYY0 = data[:,1]

		if (isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'hDeltaVsLambdaRF_'+forest__+mockNumber+'.txt')
		plt.errorbar(lambdaRFMin__-2.5+data[:,0], data[:,1], marker='o')

		plt.errorbar(lambdaRFMin__-2.5+data[:,0], data[:,1]/saveYYY0, marker='o')

#plt.title(r'$Template$', fontsize=40)
plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$Normalized \, flux$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

'''
### Template (m2)
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		if (isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)
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
		if(isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)
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
if(isMock_): 
	path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
	mockNumber = '_'+str(i)+'_'+str(j)
data = numpy.loadtxt(path+'hDeltaVsLambdaRF_'+forest__+mockNumber+'.txt')
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
		if(isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'deltaVSLambdaObs_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$'+str(i)+' \, '+str(j)+'$')
		print numpy.amin( data[:,0] )

		'''
		yMin    = numpy.min(data[:-1,1])
		yMax    = numpy.max(data[:-1,1])
		for i in range(0,skyLinesNames__.size):
			el = skyLines__[i]
			line = numpy.mean(el)
			xLi = [line,line]
			yLi = [yMin,yMax]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, skyLinesNames__[i], rotation='vertical', fontsize=20)
		'''

plt.title(r'$Residual$', fontsize=40)
plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$< \delta >$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

### Delta+1 vs. lambda_Obs
for i in range (0,chunckNb):
	for j in range(0,simulNb):

		if(isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)

		### Get the <delta>
		#data   = numpy.loadtxt(path+'deltaVSLambdaObs_'+forest__+mockNumber+'.txt')
		#interp = interpolate.interp1d(data[:,0],data[:,1],bounds_error=False,fill_value=0)
		

		data = numpy.loadtxt(path+'hDeltaVsLambdaObs_'+forest__+mockNumber+'.txt')
		plt.errorbar(data[:,0]+shift__, data[:,1], label=r'$'+str(i)+' \, '+str(j)+'$')
		#plt.errorbar(data[:,0]+shift__, data[:,1]*(1.+interp(data[:,0]+shift__)), label=r'$'+str(i)+' \, '+str(j)+'$')
		#plt.errorbar(data[:,0]+shift__, interp(data[:,0]+shift__), label=r'$'+str(i)+' \, '+str(j)+'$')

		'''
		yMin    = 0.
		yMax    = 1.
		for i in range(0,skyLinesNames__.size):
			el = skyLines__[i]
			line = numpy.mean(el)
			xLi = [line,line]
			yLi = [yMin,yMax]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, skyLinesNames__[i], rotation='vertical', fontsize=20)
		'''

#plt.title(r'$hDeltaVsLambdaObs$', fontsize=40)
plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
plt.ylabel(r'$flux$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.show()

### Eta
for i in range (0,chunckNb):
	for j in range(0,simulNb):
		if(isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)
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
		if(isMock_): 
			path = rawPath + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			mockNumber = '_'+str(i)+'_'+str(j)
		data = numpy.loadtxt(path+'sigma2LSS_'+forest__+mockNumber+'.txt')
		ax.errorbar(data[:,0], data[:,1],yerr=data[:,2], fmt='o', label=r'$Data$')

		
		xxx = data[:,0]
		yyy = data[:,1]
		def chi2(a0,a1,a2):
			return numpy.sum(numpy.power( yyy-(a0+a1*numpy.power(xxx+1.,a2)) ,2.))
		m = Minuit(chi2,a0=0.,error_a0=1.,a1=1.,error_a1=1.,a2=1.,error_a2=1.,print_level=-1, errordef=0.01, fix_a0=True) 	
		m.migrad()
		a0 = m.values['a0']
		a1 = m.values['a1']
		a2 = m.values['a2']
		
		xxx = numpy.arange(0.,5.,0.01)
		print a0, a1, a2
		ax.plot( xxx, a0+a1*numpy.power(xxx+1.,a2),label=r'$\sigma_{L.S.S.}^{2}=b.(1+z)^{\gamma}$' )
		#ax.plot( xxx, a0+a1*numpy.power(xxx+1.,a2),label=r'$\sigma_{L.S.S.}^{2}=a + b.(1+z)^{\gamma}$' )
		
plt.title(r'', fontsize=40)
plt.xlabel(r'$z_{pixel}$', fontsize=40)
plt.ylabel(r'$\sigma_{L.S.S.}^{2}$', fontsize=40)
myTools.deal_with_plot(False,False,True)

mng = plt.get_current_fig_manager()
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, r"$b=%.2e$" "\n" r"$\gamma=%.2e$"%(a1,a2), transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
#ax.text(0.05, 0.95, r"$a=%.2e$" "\n" r"$b=%.2e$" "\n" r"$\gamma=%.2e$"%(a0,a1,a2), transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)

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











