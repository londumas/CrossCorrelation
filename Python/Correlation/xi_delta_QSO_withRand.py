# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >


import sys
import numpy
import matplotlib.pyplot as plt

import myTools
from const_delta import *


nbRegion = 80
nbRandom = 10

### 1D
min1D   = 0.
max1D   = 160.
nbBin1D = 160
binSize = max1D/nbBin1D

### 2D
minX2D   = min1D
maxX2D   = max1D
nbBinX2D = nbBin1D

minY2D   = -max1D
maxY2D   = max1D
nbBinY2D = 2*nbBin1D

nbBin2D  = nbBinX2D*nbBinY2D

### Mu
nbBinM = 20;


path1 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/'

forest1 = sys.argv[1]
forest2 = sys.argv[2]
qso1    = sys.argv[3]
qso2    = sys.argv[4]

if (len(sys.argv)>5):
	wickIdx = int(sys.argv[5])



def loadData():

	xi1D     = numpy.zeros(shape=(nbBin1D,3))
	xi1Dboot = numpy.zeros(shape=(nbBin1D,nbRegion))
	binSizeX = 160/nbBin1D

	### 1D
	path  = path1 +'xi_delta_QSO_1D_'+forest1+'_'+qso1+'.txt'
	print path
	data = numpy.loadtxt(path)

	tmp_save0 = numpy.zeros(nbBin1D)
	tmp_save1 = numpy.zeros(nbBin1D)
	tmp_save2 = numpy.zeros(nbBin1D)
	tmp_save3 = numpy.zeros(nbBin1D)
	tmp_save4 = numpy.zeros(nbBin1D)

	for i in range(0,160):
		idX = i/binSizeX

		### for xi1D
		tmp_save0[idX] += data[i][1]
		tmp_save1[idX] += data[i][2]
		tmp_save2[idX] += data[i][3]
		tmp_save3[idX] += data[i][4]
		tmp_save4[idX] += data[i][5]

	for i in range(0,nbBin1D):
		yyy = tmp_save0[i]/tmp_save3[i]
		
		xi1D[i][0] = tmp_save2[i]/tmp_save3[i]
		xi1D[i][1] = yyy
		xi1D[i][2] = numpy.sqrt( (tmp_save1[i]/tmp_save3[i] - yyy*yyy)/tmp_save4[i] )
	
	### deal with random
	path  = numpy.array( [ '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test/' +'xi_delta_QSO_1D_'+forest1+'_QSO_RAND'+str(i)+'.txt' for i in range(0,nbRandom)] )

	data = numpy.loadtxt(path[0])
	for i in range(1,nbRandom):
		data += numpy.loadtxt(path[i])

	tmp_save0 = numpy.zeros(nbBin1D)
	tmp_save3 = numpy.zeros(nbBin1D)

	for i in range(0,160):
		idX = i/binSizeX
		tmp_save0[idX] += data[i][1]
		tmp_save3[idX] += data[i][4]

	for i in range(0,nbBin1D):
		xi1D[i][1] -= tmp_save0[i]/tmp_save3[i]
	
	### deal with bootstraps
	path  = numpy.array( [ path1 +'xi_delta_QSO_1D_'+forest1+'_QSO_bootstrap_'+str(i)+'.txt' for i in range(0,nbRegion)] )

	for r in range(0,nbRegion):
		data = numpy.loadtxt(path[r])

		tmp_save0 = numpy.zeros(nbBin1D)
		tmp_save3 = numpy.zeros(nbBin1D)

		for i in range(0,160):
			idX = i/binSizeX
			tmp_save0[idX] += data[i][1]
			tmp_save3[idX] += data[i][4]

		for i in range(0,nbBin1D):
			xi1Dboot[i][r] = tmp_save0[i]/tmp_save3[i]

		
		### deal with random
		path2  = numpy.array( [ path1 +'xi_delta_QSO_1D_'+forest1+'_QSO_RAND'+str(i)+'_bootstrap_'+str(r)+'.txt' for i in range(0,nbRandom)] )

		data = numpy.loadtxt(path2[0])
		for i in range(1,nbRandom):
			data += numpy.loadtxt(path2[i])

		#mean = numpy.sum(data[:,1])/numpy.sum(data[:,4])
		tmp_save0 = numpy.zeros(nbBin1D)
		tmp_save3 = numpy.zeros(nbBin1D)
	
		for i in range(0,160):
			idX = i/binSizeX
			tmp_save0[idX] += data[i][1]
			tmp_save3[idX] += data[i][4]

		for i in range(0,nbBin1D):
			xi1Dboot[i][r] -= data[i][1]/data[i][4]  #mean
	



	return xi1D, xi1Dboot
def plotXi(rescale):

	xxx = xi1D_[:,0]
	yyy = xi1D_[:,1]
	yer = xi1D_[:,2]

	xxx = xxx[yer!=0.]
	yyy = yyy[yer!=0.]
	yer = yer[yer!=0.]

	if (rescale==0):
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o')
		plt.ylabel(r'$\xi (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, fmt='o')
		plt.ylabel(r'$|s|.\xi (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, fmt='o')
		plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1+'} \, - \, '+qso1+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()


xi1D_, xi1Dboot_ = loadData()
#plotXi(0)
#plotXi(1)
#plotXi(2)

cov = numpy.cov(xi1Dboot_)
cor = numpy.cov(xi1Dboot_)
diag     = numpy.diagonal(cov)
diagSqrt = numpy.sqrt(diag)
for j in range(0,nbBin1D):
	for k in range(0,nbBin1D):
		cor[j][k] = cov[j][k]/(diagSqrt[j]*diagSqrt[k])

### matrix in 2D 
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xticks([ 0.,50.,100.,150.])
ax.set_yticks([ 0.,50.,100.,150.])

plt.imshow(cov, origin='lower',extent=[min1D, max1D, min1D, max1D], interpolation='None')
cbar = plt.colorbar()
cbar.set_label(r'$Cov(|s_{1}|,|s_{2}|)$',size=40)

#plt.title(r'$'+covNameList[i]+'$', fontsize=40)
plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
plt.grid(True)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()
plt.show()


### matrix in 2D 
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xticks([ 0.,50.,100.,150.])
ax.set_yticks([ 0.,50.,100.,150.])

plt.imshow(cor, origin='lower',extent=[min1D, max1D, min1D, max1D], interpolation='None')
cbar = plt.colorbar()
cbar.set_label(r'$Cor(|s_{1}|,|s_{2}|)$',size=40)

#plt.title(r'$'+covNameList[i]+'$', fontsize=40)
plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
plt.grid(True)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()
plt.show()

### Show the variogram
xxx   = numpy.arange(min1D, max1D, binSize)
mean = numpy.zeros(shape=(nbBin1D,2))
for i in range(0,nbBin1D):
	for j in range(0,nbBin1D-i):
		if (cor[i+j][j]!=0. and numpy.isfinite(cor[i+j][j])):
			mean[i][0] += cor[i+j][j]
			mean[i][1] += 1.

mean = numpy.nan_to_num(mean[:,0]/mean[:,1])

plt.plot(xxx, mean)

plt.xlabel(r'$\Delta \, |s|$', fontsize=40)
plt.ylabel(r'$Cor( \Delta \, |s| )$', fontsize=40)

myTools.deal_with_plot(False,False,True)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()

numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test/covBootWithRand_LYA_QSO',cov)

'''
for i in range(0,nbRegion):
	plt.errorbar(xi1D_[:,0], xi1Dboot_[:,i], fmt='o',label=r'$'+str(i)+'$')
myTools.deal_with_plot(False,False,True)
plt.show()
for i in range(0,nbRegion):
	plt.errorbar(xi1D_[:,0], xi1D_[:,0]*xi1Dboot_[:,i], fmt='o',label=r'$'+str(i)+'$')
myTools.deal_with_plot(False,False,True)
plt.show()
for i in range(0,nbRegion):
	plt.errorbar(xi1D_[:,0], xi1D_[:,0]*xi1D_[:,0]*xi1Dboot_[:,i], fmt='o',label=r'$'+str(i)+'$')
myTools.deal_with_plot(False,False,True)
plt.show()
'''

