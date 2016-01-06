# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import sys
import numpy
import matplotlib.pyplot as plt
import myTools
from const_delta import *
from iminuit import Minuit
from scipy import interpolate

forest__ = sys.argv[1]
meanZ__ = 0.
#strongLines_CIV = [1.00166645,1.00522684,1.01407602,1.01576593,1.03718845,1.03891687]


def loadData(path1D):


	print path1D
	data = numpy.loadtxt(path1D)

	cut = (data[:,5]!=0.)
	size = cut[cut].size
	xi1D = numpy.zeros(shape=(size,3))

	xi1D[:,0] = data[:,2][cut]/data[:,4][cut]
	xi1D[:,1] = data[:,0][cut]/data[:,4][cut]
	xi1D[:,2] = numpy.sqrt( (data[:,1][cut]/data[:,4][cut] -xi1D[:,1]**2.)/data[:,5][cut]  )

	meanZ__ = numpy.sum(data[:,3])/numpy.sum(data[:,4])
	print '  <z> = ', meanZ__
	
	print 

	return xi1D

def plot(xi1D):


	if (forest__=='LYA'):
		lines = LYA_lines
		names = LYA_lines_names
	if (forest__=='CIV'):
		lines = CIV_lines
		names = CIV_lines_names
	if (forest__=='MGII'):
		lines = MGII_lines
		names = MGII_lines_names
	if (forest__=='LYB'):
		lines = LYB_lines
		names = LYB_lines_names
	if (forest__=='SIIV'):
		lines = SIIV_lines
		names = SIIV_lines_names

	xxx = xi1D[:,0]
	yyy = xi1D[:,1]
	yer = xi1D[:,2]


	### Show lines in the correlation
	xMax    = numpy.amax(xxx)
	yMin    = numpy.amin(yyy)
	yMax    = numpy.amax(yyy)
	nbLines = lines.size
	for i in range(0,nbLines):
		for j in range(0,i):
			#if (names[i][:3]!=forest__ and names[j][:3]!=forest__): continue
			line = lines[j]/lines[i]

			
			#isPresent = False
			#for k in strongLines_CIV:
			#	if (abs(line-k)<=0.00005): isPresent = True
			#if (not isPresent): continue
			
			line = numpy.abs( find_dist_correlation_lines(meanZ__,lines[j],lines[i]) )
			if (line>xMax): continue
			#print line

			xLi = [line,line]
			yLi = [yMin,yMax]
			name = names[i]+' - '+names[j]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, name, rotation='vertical', fontsize=20)



	plt.errorbar(xxx,yyy, yerr=yer, fmt='o')
	plt.title(r'$1D: \, \delta_{'+forest__+'} \, - \, \delta_{'+forest__+'} $', fontsize=40)
	plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
	plt.ylabel(r'$\xi (s)$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
def saveListRealMocks(ni,nj):

	pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_1D_delta_delta_LYA_result'

	listSimu = numpy.zeros( shape=(445,ni*nj) )

	for i in numpy.arange(ni):
		for j in numpy.arange(nj):
			path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_1D_delta_delta_LYA.txt'
			xi1D = loadData(path)
			listSimu[:,i*10+j] = xi1D[:,1]

	numpy.save(pathToSave,listSimu)
	numpy.save(pathToSave+'_cov',numpy.cov(listSimu))

	return



lines = LYA_lines
for i in range(0,lines.size):
	if (abs(find_dist_correlation_lines(2.35696306629,1215.67,lines[i]))<200.):
		print lines[i], LYA_lines_names[i], find_dist_correlation_lines(2.35696306629,1215.67,lines[i])


xi1D = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_test_PDFMocksJMC/xi_1D_delta_delta_'+forest__+'.txt')
plot(xi1D)


### Data
path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_1D_delta_delta_LYA.txt'
xi1D = loadData(path)
plt.errorbar(xi1D[:,0],xi1D[:,1],yerr=xi1D[:,2],label='Data',fmt='o',color='blue')
### Data
path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_1D_delta_delta_LYA.txt'
xi1D = loadData(path)
plt.errorbar(xi1D[:,0],xi1D[:,1],yerr=xi1D[:,2],label='Data',fmt='o',color='green')
### Simu
path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/xi_1D_delta_delta_LYA.txt'
xi1D = loadData(path)
plt.errorbar(xi1D[:,0],xi1D[:,1],yerr=xi1D[:,2],label='Simu',fmt='o',color='red')
### Simu
path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_withNoNorma_TestAllSkyLines/xi_1D_delta_delta_LYA.txt'
xi1D = loadData(path)
plt.errorbar(xi1D[:,0],xi1D[:,1],yerr=xi1D[:,2],label='Simu',fmt='o',color='black')
### Simu
path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_test_PDFMocksJMC/xi_1D_delta_delta_LYA.txt'
xi1D = loadData(path)
plt.errorbar(xi1D[:,0],xi1D[:,1],yerr=xi1D[:,2],label='Simu',fmt='o',color='cyan')


plt.title(r'$1D: \, \delta_{LYA} \, - \, \delta_{LYA} $', fontsize=40)
plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$\xi (s)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.xlim([ -10., numpy.max(xi1D[:,0])+10. ])
plt.show()



path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_1D_delta_delta_CIV.txt'
xi1D = loadData(path)
plot(xi1D)





### Numbers
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_1D_delta_delta_LYA.txt')
data[:,2] /= data[:,4]
plt.errorbar(data[:,2],data[:,5],fmt='o',color='blue',label='Data')

data_SAVE = numpy.array(numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_1D_delta_delta_LYA.txt'))
data_SAVE.fill(0.)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_1D_delta_delta_LYA.txt'
		data = numpy.loadtxt(path)
		data[:,2] /= data[:,4]
		data_SAVE += data
data_SAVE /= 100.
plt.errorbar(data_SAVE[:,2],data_SAVE[:,5],fmt='o',color='red',label='<Simu>')
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$nb \, pairs \, LYA-LYA \, same \, forest$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()



### Errors
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_1D_delta_delta_LYA.txt')
data[:,2] /= data[:,4]
data[:,3] = (data[:,1]/data[:,4] -(data[:,0]/data[:,4])**2.)
plt.errorbar(data[:,2],data[:,3],fmt='o',color='blue',label='Data')

data_SAVE = numpy.array(numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_1D_delta_delta_LYA.txt'))
data_SAVE.fill(0.)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_1D_delta_delta_LYA.txt'
		data = numpy.loadtxt(path)
		data[:,2] /= data[:,4]
		data[:,3] = (data[:,1]/data[:,4] -(data[:,0]/data[:,4])**2.)
		data_SAVE += data
data_SAVE /= 100.
plt.errorbar(data_SAVE[:,2],data_SAVE[:,3],fmt='o',color='red',label='<Simu>')

plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$Var \, LYA-LYA \, same \, forest$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()








#saveListRealMocks(10,10)
#myTools.plotCovar([numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_1D_delta_delta_LYA_result_cov.npy')],['a'])
### Data
path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_1D_delta_delta_LYA.txt'
xi1D = loadData(path)
plt.errorbar(xi1D[:,0],xi1D[:,1],yerr=xi1D[:,2],label='Data',fmt='o',color='blue')

### Simu
path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_1D_delta_delta_LYA.txt'
xi1D = loadData(path)
xi1D[:,1] = numpy.mean(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_1D_delta_delta_LYA_result.npy'),axis=1)
xi1D[:,2] = numpy.sqrt(numpy.diag(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_1D_delta_delta_LYA_result_cov.npy')))/numpy.sqrt(10.)
plt.errorbar(xi1D[:,0],xi1D[:,1],yerr=xi1D[:,2],label='Simu',fmt='o',color='red')

plt.title(r'$1D: \, \delta_{LYA} \, - \, \delta_{LYA} $', fontsize=40)
plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$\xi (s)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.xlim([ 4., numpy.max(xi1D[:,0])+10. ])
plt.show()
		





cov1D = numpy.cov(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_1D_delta_delta_LYA_result_cov.npy'))
myTools.plot2D(cov1D)
myTools.plot2D(myTools.getCorrelationMatrix(cov1D))























#plot(xi1D)
