# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import myTools
from const_delta import *

import sys
import numpy
import matplotlib.pyplot as plt
import pandas
from iminuit import Minuit

nbRegion__ = 0

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

### Mu
nbBinM__ = 20;


path1__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/'

### Parameters
forest1__ = sys.argv[1]
forest2__ = sys.argv[2]
qso1__    = sys.argv[3]
qso2__    = sys.argv[4]
if (len(sys.argv)>5):
	wickIdx__ = int(sys.argv[5])



def loadData(path1D, path2D):

	
	### 2D
	data = numpy.loadtxt(path2D)

	save0 = data[:,5]
	save1 = data[:,6]
	save2 = data[:,7]
	save3 = data[:,8]
	save4 = data[:,9]
	save5 = data[:,10]
	save6 = data[:,0]

	tmp_save0  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save1  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save2  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save3  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save4  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save5  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))

	for i in range(0,len(save0)):
		iX = int(save6[i])/int(maxY2D__-minY2D__)
		iY = int(save6[i])%int(maxY2D__-minY2D__)

		idX = iX/int(binSize__)
		idY = iY/int(binSize__)

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]
		tmp_save2[idX][idY] += save2[i]
		tmp_save3[idX][idY] += save3[i]
		tmp_save4[idX][idY] += save4[i]
		tmp_save5[idX][idY] += save5[i]

	xi2D = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__,3))

	for i in range(0,nbBinX2D__):
		for j in range(0,nbBinY2D__):

			xxx = numpy.sqrt( (tmp_save2[i][j]/tmp_save4[i][j])**2. + (tmp_save3[i][j]/tmp_save4[i][j])**2. )
			yyy = tmp_save0[i][j] / tmp_save4[i][j]
			yer = numpy.sqrt( (tmp_save1[i][j]/tmp_save4[i][j] - yyy*yyy)/tmp_save5[i][j] )
			xi2D[i][j][0] = xxx
			xi2D[i][j][1] = yyy
			xi2D[i][j][2] = yer
	
	### Mu
	data = numpy.loadtxt(path1D)

	save0 = data[:,1]
	save1 = data[:,2]
	save2 = data[:,3]
	save3 = data[:,4]
	save4 = data[:,5]
	save5 = data[:,0].astype(int)

	tmp_save0 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save1 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save2 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save3 = numpy.zeros(shape=(nbBin1D__,nbBinM__))
	tmp_save4 = numpy.zeros(shape=(nbBin1D__,nbBinM__))

	tmp_save00 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save11 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save22 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save33 = numpy.zeros(shape=(nbBin1D__,3))
	tmp_save44 = numpy.zeros(shape=(nbBin1D__,3))

	tmp_save000 = numpy.zeros(nbBin1D__)
	tmp_save111 = numpy.zeros(nbBin1D__)
	tmp_save222 = numpy.zeros(nbBin1D__)
	tmp_save333 = numpy.zeros(nbBin1D__)
	tmp_save444 = numpy.zeros(nbBin1D__)

	binSizeX = int(max1D__)/nbBin1D__
	binSizeY = 100/nbBinM__

	
	for i in range(0,len(save0)):
		iX = save5[i]/100
		iY = save5[i]%100

		### for mu
		idX = iX/binSizeX
		idY = iY/binSizeY

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]
		tmp_save2[idX][idY] += save2[i]
		tmp_save3[idX][idY] += save3[i]
		tmp_save4[idX][idY] += save4[i]

		### for wedges
		if (iY<10 or iY>90):
			idY = 0
		elif ( (iY>=10 and iY<25) or (iY<=90 and iY>75) ):
			idY = 1
		else:
			idY = 2

		tmp_save00[idX][idY] += save0[i]
		tmp_save11[idX][idY] += save1[i]
		tmp_save22[idX][idY] += save2[i]
		tmp_save33[idX][idY] += save3[i]
		tmp_save44[idX][idY] += save4[i]

		### for xi1D
		tmp_save000[idX] += save0[i]
		tmp_save111[idX] += save1[i]
		tmp_save222[idX] += save2[i]
		tmp_save333[idX] += save3[i]
		tmp_save444[idX] += save4[i]

	xiMu = numpy.zeros(shape=(nbBin1D__,nbBinM__,3))
	xiWe = numpy.zeros(shape=(nbBin1D__,3,3))
	xi1D = numpy.zeros(shape=(nbBin1D__,3))

	for i in range(0,nbBin1D__):
		for j in range(0,nbBinM__):
			if (tmp_save4[i][j]==0.): continue
			xxx = tmp_save2[i][j]/tmp_save3[i][j]
			yyy = tmp_save0[i][j]/tmp_save3[i][j]
			yer = numpy.sqrt( (tmp_save1[i][j]/tmp_save3[i][j] - yyy*yyy)/tmp_save4[i][j] )
			xiMu[i][j][0] = xxx
			xiMu[i][j][1] = yyy
			xiMu[i][j][2] = yer
		for j in range(0,3):
			if (tmp_save44[i][j]==0.): continue
			xxx = tmp_save22[i][j]/tmp_save33[i][j]
			yyy = tmp_save00[i][j]/tmp_save33[i][j]
			yer = numpy.sqrt( (tmp_save11[i][j]/tmp_save33[i][j] - yyy*yyy)/tmp_save44[i][j] )
			xiWe[i][j][0] = xxx
			xiWe[i][j][1] = yyy
			xiWe[i][j][2] = yer

		if (tmp_save444[i]==0.): continue
		xxx = tmp_save222[i]/tmp_save333[i]
		yyy = tmp_save000[i]/tmp_save333[i]
		yer = numpy.sqrt( (tmp_save111[i]/tmp_save333[i] - yyy*yyy)/tmp_save444[i] )

		xi1D[i][0] = xxx
		xi1D[i][1] = yyy
		xi1D[i][2] = yer


	return xi1D, xi2D, xiMu, xiWe
def createIni2():
	for ii in range(0,10):
		for jj in range(0,10):
		
			path = str(ii)+'_'+str(jj)
			data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/xi_delta_QSO_2D_LYA_DLA_mocks__'+path+'_.txt')
			
			save2  = data[:,7]
			save3  = data[:,8]
			save4  = data[:,9]
			save11 = data[:,11]
			save6  = data[:,0]
			
			meanRperp  = numpy.zeros(shape=(nbBinX2D__,2))
			meanRparal = numpy.zeros(shape=(nbBinY2D__,2))
			
			for i in range(0,len(save2)):
				iX = int(save6[i])/int(maxY2D__-minY2D__)
				iY = int(save6[i])%int(maxY2D__-minY2D__)
			
				idX = iX/int(binSize__)
				idY = iY/int(binSize__)
			
				meanRperp[idX][0]  += save2[i]
				meanRperp[idX][1]  += save4[i]
				meanRparal[idY][0] += save3[i]
				meanRparal[idY][1] += save4[i]
			
			meanRperp[:,0]  /= meanRperp[:,1]
			meanRparal[:,0] /= meanRparal[:,1]
			
			### Get the s_perp bin center
			stringRperp = ''
			for el in meanRperp[:-1,0]:
				stringRperp += str(el) + ','
			stringRperp += str(meanRperp[-1,0])
			
			### Get the s_paral bin center
			stringRparal = ''
			for el in meanRparal[:-1,0]:
				stringRparal += str(el) + ','
			stringRparal += str(meanRparal[-1,0])

			### Get the redshift center
			stringRedshift = str(numpy.sum(save11)/numpy.sum(save4))
			print ' <z> = ', stringRedshift

			string = """
modelroot =  /home/helion/Documents/Thèse/Code/Libraries/baofit/models/
fiducial  = DR9LyaMocksNLeffPk
nowiggles = DR9LyaMocksSB

cross-correlation = yes
anisotropic       = yes
decoupled         = yes

model-config = value[beta]=1.1;
model-config = fix[(1+beta)*bias]=-0.336;
model-config = fix[gamma-bias]=0.9; fix[gamma-beta]=0;
model-config = value[bias2]=1.5;
model-config = value[beta2*bias2]=0.962524;
model-config = value[delta-v]=0;
model-config = value[BAO amplitude]=1;
model-config = fix[BAO alpha-iso]; fix[BAO alpha-p*]=1;
model-config = fix[gamma-scale]=0;

model-config = binning[beta]={0.5:3.}*50
model-config = binning[BAO alpha-perp]={0.6:1.3}*50

dist-add = rP,rT=0:2,-3:1
data     = /home/hdumasde/Documents/Documents/Thèse/Results/Fit_Bao/DLA_systematiq/bao2D_"""+path+"""

data-format = comoving-cartesian

axis1-bins = {""" + stringRparal   + """}
axis2-bins = {""" + stringRperp    + """}
axis3-bins = {""" + stringRedshift + """}

rmin = 40
rmax = 180

output-prefix = /home/hdumasde/Documents/Documents/Thèse/Results/Fit_Bao/DLA_systematiq/bao2D_atOneOne_"""+path+""".

alt-config = fix[dist*]=0

ndump = 0
"""

		
			text_file = open('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/DLA_systematiq/bao2D_'+path+'.ini', "w")
			text_file.write(string)
			text_file.close()
		
	return
def loadBoot(path1D, path2D, nb, nb2=1, string='', chunk=-1, simul=-1):

	### Number of bin in the file
	nbBinFile1D = int(max1D__)
	nbBinFile2D = int(max1D__*max1D__*2)

	binSizeX = int(max1D__)/nbBin1D__
	aaa      = int(maxY2D__-minY2D__)
	binSize  = int(binSize__)

	### Where to store the bootstrap
	xi1Dboot = numpy.zeros(shape=(nbBin1D__,nb*nb2))
	xi2Dboot = numpy.zeros(shape=(nbBin2D__,nb*nb2))

	pathList1D = []
	pathList2D = []

	if (string==''):
		pathList1D = [ path1D+str(i)+'.txt' for i in range(0,nb) ]
		pathList2D = [ path2D+str(i)+'.txt' for i in range(0,nb) ]
	elif (string=='mock'):
		if (chunk==-1 and simul==-1):
			for i in range(0,nb):
				for j in range(0,nb2):
					pathList1D += [ path1D+str(i)+'_'+str(j)+'_.txt' ]
					pathList2D += [ path2D+str(i)+'_'+str(j)+'_.txt' ]
		else:
			endTxt = '_mocks__'+str(chunk)+'_'+str(simul)+'_.txt'
			pathList1D = [ path1D+str(i)+endTxt for i in range(0,nb) ]
			pathList2D = [ path2D+str(i)+endTxt for i in range(0,nb) ]

	for b in range(0,len(pathList1D)):


		#if (b==1528 or b==1749 or b==1757 or b==1778 or b==1805 or b==1810 or b==1835 or b==1844 or b==1855 or b==1867 or b==1874 or b==1897 or b==1971 or b==1980): continue

		print pathList1D[b]
		print pathList2D[b]

		### 1D
		data = numpy.array(pandas.read_csv(pathList1D[b], header=None, sep=' ', lineterminator='\n', dtype={'a': numpy.dtype('d')}))

		save0 = data[:,1]
		save3 = data[:,4]

		tmp_save000 = numpy.zeros(nbBin1D__)
		tmp_save333 = numpy.zeros(nbBin1D__)
		
		for i in range(0,nbBinFile1D):
			idX = i/binSize
	
			### for xi1D
			tmp_save000[idX] += save0[i]
			tmp_save333[idX] += save3[i]

		xi1Dboot[:,b] = tmp_save000/tmp_save333
		
		### 2D
		data = numpy.array(pandas.read_csv(pathList2D[b], header=None, sep=' ', lineterminator='\n', dtype={'a': numpy.dtype('d')}))

		save0 = data[:,5]
		save4 = data[:,9]
		save6 = data[:,0].astype(int)
		lenSave0 = len(save0)
	
		tmp_save0 = numpy.zeros(nbBin2D__)
		tmp_save4 = numpy.zeros(nbBin2D__)
	
		for i in range(0,lenSave0):
			iX = save6[i]/aaa
			iY = save6[i]%aaa
	
			idX = iX/binSize
			idY = iY/binSize
			idK = idX*nbBinY2D__+idY
	
			tmp_save0[idK] += save0[i]
			tmp_save4[idK] += save4[i]

		xi2Dboot[:,b] = tmp_save0/tmp_save4
		
	return xi1Dboot, xi2Dboot
def loadWick():

	nb = '9995'
	nbBinFile1D = int(max1D__)

	### Matrix with the divident of all wicks
	wick1D_base = numpy.zeros(shape=(nbBin1D__,nbBin1D__))
	### Matrix to get the sum of '1' if all done
	wick1D_base_ones  = numpy.zeros(shape=(nbBin1D__,nbBin1D__))
	### Matrix to get the sum of the diagrams
	wick1D_sum  = numpy.zeros(shape=(nbBin1D__,nbBin1D__))
	### Matrix to get the sum of '1' of the diagrams
	wick1D_sum_ones  = numpy.zeros(shape=(nbBin1D__,nbBin1D__))
	### Create the diagram matrix
	wick1D_Ti  = [numpy.zeros(shape=(nbBin1D__,nbBin1D__))]
	wick1D_Ti += [numpy.zeros(shape=(nbBin1D__,nbBin1D__))]
	wick1D_Ti += [numpy.zeros(shape=(nbBin1D__,nbBin1D__))]
	wick1D_Ti += [numpy.zeros(shape=(nbBin1D__,nbBin1D__))]
	wick1D_Ti += [numpy.zeros(shape=(nbBin1D__,nbBin1D__))]
	wick1D_Ti += [numpy.zeros(shape=(nbBin1D__,nbBin1D__))]
	### Get the name of the diagrams
	wick1D_Ti_Name = ['T1','T12','T123','T1234','T12345','T123456']



	### Get the cross-correlation
	xi1D     = numpy.zeros(nbBin1D__)
	data     = numpy.loadtxt(path1__ + 'xi_delta_QSO_1D_Wick1_LYA_QSO_'+nb+'_crossCorrelation.txt')

	save1 = data[:,1]
	save2 = data[:,2]
	save3 = data[:,3]

	xi1D = save1/save2
	#plt.plot(numpy.arange(0,xi1D.size),xi1D)
	#plt.show()



	### Set the devident of all Wick
	for i in range(0,nbBin1D__):
		for j in range(0,nbBin1D__):
			wick1D_base[i][j] = 1./(save2[i]*save2[j])
			wick1D_base_ones[i][j] = save3[i]*save3[j]


	for i in range(0,4):
		tmp_save2 = numpy.zeros(shape=(nbBin1D__,nbBin1D__))

		path = path1__ +'xi_delta_QSO_1D_Wick'+str(i+1)+'_'+forest1__+'_'+qso1__+'_'+nb+'.txt'
		print path
		data = numpy.loadtxt(path)

		save0 = data[:,0].astype(int)
		save1 = data[:,1].astype(int)
		save2 = data[:,2]
		save3 = data[:,3]
		save4 = data[:,4]
		size  = len(save0)

		for j in range(0,size):
			idX = save0[j]
			idY = save1[j]

			tmp_save2[idX][idY]  += save2[j]
			wick1D_sum[idX][idY] += save3[j]
			wick1D_sum_ones[idX][idY] += save4[j]

		wick1D_Ti[i] = tmp_save2*wick1D_base

	### Add the diagrams
	for i in range(1,6):
		wick1D_Ti[i] += wick1D_Ti[i-1]

	### Show the sum
	for i in range(0,nbBin1D__):
		print i, i, wick1D_sum[i][i]*wick1D_base[i][i], wick1D_sum[i][i], 1./wick1D_base[i][i], wick1D_sum_ones[i][i]/wick1D_base_ones[i][i], wick1D_sum_ones[i][i], wick1D_base_ones[i][i]

	### Set T6
	for i in range(0,nbBin1D__):
		for j in range(0,nbBin1D__):
			wick1D_Ti[5][i][j] = xi1D[i]*xi1D[j]

	
	for k in range(0,6):

		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(wick1D_Ti[k], origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cov(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$'+wick1D_Ti_Name[k]+'$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()



		corList      = numpy.zeros(shape=(nbBin1D__,nbBin1D__))
		diagList     = numpy.diagonal(wick1D_Ti[k])
		diagSqrtList = numpy.sqrt(diagList)

		for i in range(0,nbBin1D__):
			for j in range(0,nbBin1D__):
				corList[i][j] = wick1D_Ti[k][i][j]/(diagSqrtList[i]*diagSqrtList[j])
	
		corList = numpy.nan_to_num(corList)

		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(corList, origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cor(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$'+wick1D_Ti_Name[k]+'$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()
	
	return wick1D_Ti
def loadBootMap():
	'''
	'''

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test/test_map.txt'
	
	data = numpy.loadtxt(path)
	re = data[:,0].astype(int)
	ra = data[:,1]
	de = data[:,2]
	
	for i in range(0,numpy.amax(re)+1):
		cut = (re==i)
		plt.plot(ra[cut], de[cut], linestyle="", marker="o")

	#plt.xlim([0,360.])
	#plt.ylim([-90.,90.])
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.xlabel(r'$R.A. (\degree)$')
	plt.ylabel(r'$Dec. (\degree)$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	return
def convert1DTo2D(array1D):
	'''
		convert a 1D array to a 2D array
	'''

	array2D = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))

	for k in range(0,nbBin2D__):
		i = k/nbBinY2D__
		j = k%nbBinY2D__

		array2D[i][j] = array1D[k]

	return array2D
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
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
def plotXi2D(rescale):

	xxx = numpy.transpose(xi2D_[:,:,0])
	yyy = numpy.transpose(xi2D_[:,:,1])
	yer = numpy.transpose(xi2D_[:,:,2])

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xticks([ 0.,50.,100.,150.,200.])
	ax.set_yticks([ -200.,-150.,-100.,-50.,0.,50.,100.,150.,200.])

	if (rescale==0):
		plt.imshow(yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__], interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		plt.imshow(xxx*yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|.\xi(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
	if (rescale==2):
		plt.imshow(xxx**2.*yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|^{2}.\xi(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
def plotMu(rescale):

	xxx = xiMu_[:,:,0]
	yyy = xiMu_[:,:,1]
	yer = xiMu_[:,:,2]

	fig = plt.figure()
	ax = fig.add_subplot(111)
	#ax.set_yticks([ 0.,50.,100.,150.])

	if (rescale==0):
		plt.imshow(yyy, origin='lower', interpolation='None',extent=[-1., 1.,min1D__, max1D__],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		plt.imshow(xxx*yyy, origin='lower', interpolation='None',extent=[-1., 1.,min1D__, max1D__],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \, [h^{-1}.Mpc])$',size=40)
	if (rescale==2):
		plt.imshow(xxx*xxx*yyy, origin='lower', interpolation='None',extent=[-1., 1.,min1D__, max1D__],aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \, [(h^{-1}.Mpc)^{2}])$',size=40)


	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$\mu$', fontsize=40)
	plt.ylabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
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
		plt.ylabel(r'$\xi (|s|)$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	if (rescale==1):
		plt.errorbar(xxx0, xxx0*yyy0, yerr=xxx0*yer0, fmt='o', label=r'$0.8 < |\mu|$')
		plt.errorbar(xxx1, xxx1*yyy1, yerr=xxx1*yer1, fmt='o', label=r'$0.5 < |\mu| \leq 0.8$')
		plt.errorbar(xxx2, xxx2*yyy2, yerr=xxx2*yer2, fmt='o', label=r'$|\mu| \leq 0.5$')
		plt.ylabel(r'$|s|.\xi (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	if (rescale==2):
		plt.errorbar(xxx0, xxx0*xxx0*yyy0, yerr=xxx0*xxx0*yer0, fmt='o', label=r'$0.8 < |\mu|$')
		plt.errorbar(xxx1, xxx1*xxx1*yyy1, yerr=xxx1*xxx1*yer1, fmt='o', label=r'$0.5 < |\mu| \leq 0.8$')
		plt.errorbar(xxx2, xxx2*xxx2*yyy2, yerr=xxx2*xxx2*yer2, fmt='o', label=r'$|\mu| \leq 0.5$')
		plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=2)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
def plotCovar(covList, covNameList):

	### Diagrams
	###################################################

	nbCov       = len(covList)

	nbBin = int(numpy.sqrt(covList[0].size))
	print '  nbBin = ', nbBin
	
	corList      = [numpy.zeros(shape=(nbBin,nbBin))]*nbCov
	diagList     = [numpy.diagonal(covList[i]) for i in range(0,nbCov)]
	diagSqrtList = [numpy.sqrt(diagList[i]) for i in range(0,nbCov)]

	for i in range(0,nbCov):
		for j in range(0,nbBin):
			for k in range(0,nbBin):
				corList[i][j][k] = covList[i][j][k]/(diagSqrtList[i][j]*diagSqrtList[i][k])
	
		corList[i] = numpy.nan_to_num(corList[i])
		'''
		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(covList[i], origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cov(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$'+covNameList[i]+'$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()


		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(corList[i], origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cor(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$'+covNameList[i]+'$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()
		'''
	'''
	### Covariance
	###################################################
	tmp_xxx = []
	tmp_yyy = []
	for c in range(0,nbCov):
		if (nbBin==nbBin1D__):
			xxx = binSize__*numpy.arange(0,nbBin)
			yyy = numpy.zeros(shape=(1,nbBin,2))
			for k1 in range(0,nbBin):
				for k2 in range(0,k1+1):
					if (covList[c][k1][k2]!=0.):
						yyy[0][k1-k2][0] += covList[c][k1][k2]
						yyy[0][k1-k2][1] += 1.
		else:
			xxx   = binSize__*numpy.arange(0, nbBinY2D__)
			xxx2  = binSize__*numpy.arange(0, nbBinX2D__)
			yyy  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__,2))
			yyy2 = numpy.zeros(shape=(nbBinY2D__,nbBinX2D__,2))
			for k1 in range(0,nbBin):
				i1 = k1/nbBinY2D__
				j1 = k1%nbBinY2D__
				for k2 in range(0,k1+1):
					i2 = k2/nbBinY2D__
					j2 = k2%nbBinY2D__
					if (covList[c][k1][k2]==0.): continue
					
					yyy[abs(i1-i2)][abs(j1-j2)][0]  += covList[c][k1][k2]
					yyy[abs(i1-i2)][abs(j1-j2)][1]  += 1.
					yyy2[abs(j1-j2)][abs(i1-i2)][0] += covList[c][k1][k2]
					yyy2[abs(j1-j2)][abs(i1-i2)][1] += 1.

		
		yyy  = numpy.nan_to_num(yyy[:,:,0]/yyy[:,:,1])

		### 0
		tmp_xxx += [xxx]
		tmp_yyy += [yyy[0]]
		plt.plot(xxx, yyy[0], label=r'$'+covNameList[c]+' \, \Delta \, s_{\perp, idx} = '+str(0)+'$', marker='o')

		plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
		plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
		plt.ylabel(r'$Cov( \Delta \, |s| \, bin \, idx  )$', fontsize=40)

		myTools.deal_with_plot(False,False,True)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.show()

		if (nbBin!=nbBin1D__):

			for i in range(1,nbBinX2D__):
				plt.plot(xxx, yyy[i], label=r'$'+covNameList[c]+' \, \Delta \, s_{\perp, idx} = '+str(i)+'$', marker='o')
	
			plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
			plt.ylabel(r'$Cov( \Delta \, |s| \, bin \, idx  )$', fontsize=40)

			myTools.deal_with_plot(False,False,True)
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.show()

			yyy2 = numpy.nan_to_num(yyy2[:,:,0]/yyy2[:,:,1])

			### 0
			plt.plot(xxx2, yyy2[0], label=r'$'+covNameList[c]+' \, \Delta \, s_{\parallel, idx} = '+str(0)+'$', marker='o')

			plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
			plt.ylabel(r'$Cov( \Delta \, |s|  \, bin \, idx )$', fontsize=40)

			myTools.deal_with_plot(False,False,True)
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.show()
			for i in range(1,nbBinY2D__):
				plt.plot(xxx2, yyy2[i], label=r'$'+covNameList[c]+' \, \Delta \, s_{\parallel, idx} = '+str(i)+'$', marker='o')
			plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
			plt.ylabel(r'$Cov( \Delta \, |s|  \, bin \, idx )$', fontsize=40)

			myTools.deal_with_plot(False,False,True)
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.show()
		
	### Show all the method on the main variogram
	for c in range(0,nbCov):
		if (nbBin==nbBin1D__):
			xTitle = '|\Delta \, |s_{1}|-|s_{2}| | \, [h^{-1}.Mpc]'
		else:
			xTitle = '|\Delta \, s_{\parallel,1}-s_{\parallel,2} | \, [h^{-1}.Mpc] \, for \, |\Delta \, s_{\perp}| = 0'
		plt.plot(tmp_xxx[c], tmp_yyy[c], label=r'$'+covNameList[c]+'$', marker='o')

	plt.xlabel(r'$'+xTitle+'$', fontsize=40)
	plt.ylabel(r'$Cov( s_{1}, s_{2})$', fontsize=40)

	myTools.deal_with_plot(False,False,True)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
	'''
	### Correlation
	###################################################

	### Show the variogram
	tmp_xxx  = []
	tmp_yyy  = []
	for c in range(0,nbCov):
		if (nbBin==nbBin1D__):
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
		'''
		plt.plot(xxx, yyy[0], label=r'$'+covNameList[c]+' \, \Delta \, s_{\perp, idx} = '+str(0)+'$', marker='o')

		plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
		plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
		plt.ylabel(r'$Cor( \Delta \, |s| \, bin \, idx )$', fontsize=40)

		myTools.deal_with_plot(False,False,True)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.show()

		if (nbBin!=nbBin1D__):

			for i in range(1,nbBinX2D__):
				plt.plot(xxx, yyy[i], label=r'$'+covNameList[c]+' \, \Delta \, s_{\perp, idx} = '+str(i)+'$', marker='o')

			plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
			plt.ylabel(r'$Cor( \Delta \, |s| \, bin \, idx  )$', fontsize=40)

			myTools.deal_with_plot(False,False,True)
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.show()

			yyy2 = numpy.nan_to_num(yyy2[:,:,0]/yyy2[:,:,1])

			### 0
			plt.plot(xxx2, yyy2[0], label=r'$'+covNameList[c]+' \, \Delta \, s_{\parallel, idx} = '+str(0)+'$', marker='o')

			plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
			plt.ylabel(r'$Cor( \Delta \, |s| \, bin \, idx  )$', fontsize=40)

			myTools.deal_with_plot(False,False,True)
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.show()
			for i in range(1,nbBinY2D__):
				plt.plot(xxx2, yyy2[i], label=r'$'+covNameList[c]+' \, \Delta \, s_{\parallel, idx} = '+str(i)+'$', marker='o')
			plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
			plt.ylabel(r'$Cor( \Delta \, |s| \, bin \, idx  )$', fontsize=40)

			myTools.deal_with_plot(False,False,True)
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.show()
		'''

	parameterFromFit = numpy.zeros(shape=(nbCov,nbBinX2D__,3))

	### Show all the method on the main variogram
	for i in range(0,nbBinX2D__):
		if (i>=1 and nbBin==nbBin1D__): continue
		
		for c in range(0,nbCov):

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

		
			### Plot the results
			if (nbBin==nbBin1D__):
				xTitle = '|\Delta \, |s_{1}|-|s_{2}| | \, [h^{-1}.Mpc]'
			else:
				xTitle = '|\Delta \, s_{\parallel,1}-s_{\parallel,2} | \, [h^{-1}.Mpc], \, for \, |\Delta \, s_{\perp}| = ' + str(i*binSize__)
			plt.plot(tmp_xxx[c], tmp_yyy[c][i], label=r'$'+covNameList[c]+'$', marker='o')
			plt.plot(tmp_xxx[c], minuit_yyy, label=r'$model$', marker='o')

			### Save the results of the fit
			

		plt.xlabel(r'$'+xTitle+'$', fontsize=40)
		plt.ylabel(r'$Cor( s_{1}, s_{2})$', fontsize=40)

		myTools.deal_with_plot(False,False,True)
		plt.show()
		


	### Get the covaraince from fit
	covFromFitList = [numpy.zeros(shape=(nbBin,nbBin))]*nbCov
	corFromFitList = [numpy.zeros(shape=(nbBin,nbBin))]*nbCov
	for c in range(0,nbCov):
		if (nbBin==nbBin1D__):
			for k1 in range(0,nbBin):
				for k2 in range(0,k1+1):
					if (corList[c][k1][k2]==0.): continue
					x         = numpy.abs(k1-k2)*binSize__
					y_fromFit = (parameterFromFit[c][0][0] + parameterFromFit[c][0][1]*x)/(1. + parameterFromFit[c][0][2]*x)
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

					y_fromFit = (parameterFromFit[c][x2][0] + parameterFromFit[c][x2][1]*x1)/(1. + parameterFromFit[c][x2][2]*x1)
					covFromFitList[c][k1][k2] = y_fromFit*diagSqrtList[c][k1]*diagSqrtList[c][k2]
					corFromFitList[c][k1][k2] = y_fromFit
					covFromFitList[c][k2][k1] = covFromFitList[c][k1][k2]
					corFromFitList[c][k2][k1] = y_fromFit
		
		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(corFromFitList[c], origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cov(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$a$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()
		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(corList[c], origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cov(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$a$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()
		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(covFromFitList[c], origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cov(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$a$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()
		### matrix in 2D 
		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.imshow(covList[c], origin='lower', interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$Cov(|s_{1}|,|s_{2}|)$',size=40)

		plt.title(r'$a$', fontsize=40)
		plt.xlabel(r'$|s_{1}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$|s_{2}| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()
		
	return covFromFitList


createIni2()

'''
for ii in range(0,10):
	for jj in range(0,10):
		path = str(ii)+'_'+str(jj)

		xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_DLA_mocks__'+path+'_.txt', path1__ +'xi_delta_QSO_2D_LYA_DLA_mocks__'+path+'_.txt')
		xi2DCov  = numpy.load(path1__ + 'mockDLA_'+path+'_subsampling_fromFit_2D.npy')[0]

		tmp_xi2D        = numpy.zeros(nbBin2D__)
		tmp_idx         = numpy.arange(0,nbBin2D__)
		nbBin2D = nbBin2D__*(nbBin2D__+1.)/2.
		tmp_xi2DCov     = numpy.zeros(nbBin2D)
		tmp_xi2DCovIdx  = numpy.zeros(nbBin2D)
		tmp_xi2DCovIdx2 = numpy.zeros(nbBin2D)
		idx = 0
		
		for k1 in range(0,nbBin2D__):
			i1       = k1/nbBinY2D__
			j1       = k1%nbBinY2D__
			k11      = j1*nbBinX2D__ + i1
			tmp_xi2D[k11] = xi2D_[i1][j1][1]
		
			for k2 in range(0,k1+1):
				i2       = k2/nbBinY2D__
				j2       = k2%nbBinY2D__
				k22      = j2*nbBinX2D__ + i2
				tmp_xi2DCov[idx]     = xi2DCov[k1][k2]
				tmp_xi2DCovIdx[idx]  = k11
				tmp_xi2DCovIdx2[idx] = k22
				idx += 1
		
		numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/DLA_systematiq/bao2D_'+path+'.data',zip(tmp_idx,tmp_xi2D),fmt='%d %E')
		numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/DLA_systematiq/bao2D_'+path+'.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%d %d %E')
'''	
xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_QSO.txt')
plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)

plt.plot(xi1D_[:,0],xi1D_[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],(xi1D2_[:,1]-xi1D2_[-1,1]), marker='o',label=r'$simu2$')
plt.plot(xi1D[:,0],xi1D[:,1], marker='o',label=r'$Data$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$\xi (|s|)$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plt.plot(xi1D_[:,0],xi1D_[:,0]*xi1D_[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],xi1D2_[:,0]*(xi1D2_[:,1]-xi1D2_[-1,1]), marker='o',label=r'$simu2$')
plt.plot(xi1D[:,0],xi1D[:,0]*xi1D[:,1], marker='o',label=r'$Data$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$|s|.\xi (|s|) \, [(h^{-1}.Mpc)]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plt.plot(xi1D_[:,0],xi1D_[:,0]*xi1D_[:,0]*xi1D_[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],xi1D2_[:,0]*xi1D2_[:,0]*(xi1D2_[:,1]-xi1D2_[-1,1]), marker='o',label=r'$simu2$')
plt.plot(xi1D[:,0],xi1D[:,0]*xi1D[:,0]*xi1D[:,1], marker='o',label=r'$Data$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plt.plot(xi1D_[:,0],xi1D_[:,1]/xi1D[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],(xi1D2_[:,1]-xi1D2_[-1,1])/xi1D[:,1], marker='o',label=r'$simu2$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$\xi (|s|) \, / \, \xi_{data} (|s|)$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)



xi1Dboot_, xi2Dboot_   = loadBoot(path1__ +'xi_delta_QSO_1D_LYA_QSO_randomQSO_', path1__ +'xi_delta_QSO_2D_LYA_QSO_randomQSO_',2000)
numpy.save(path1__ + 'randomQSO_1D',xi1Dboot_)
numpy.save(path1__ + 'randomQSO_2D',xi2Dboot_)
### For Rand
xi1DList = numpy.load(path1__ + 'randomQSO_1D.npy')
xi2DList = numpy.load(path1__ + 'randomQSO_2D.npy')
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_QSO.txt', path1__ +'xi_delta_QSO_2D_LYA_QSO.txt')
xxx  = xi1D_[:,0]
coef = 1. #(xi1D_[:,0]**2.)
mean = numpy.mean(xi1DList,axis=1)
plt.errorbar(xxx, coef*xi1D_[:,1], yerr=coef*xi1D_[:,2], marker='o', color='blue', label=r'$Data$')
plt.errorbar(xxx, coef*mean, marker='o', color='red', label=r'$mean \, over \, random$')

myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$a$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

### For Rand
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_QSO.txt', path1__ +'xi_delta_QSO_2D_LYA_QSO.txt')
xxx  = numpy.arange(0,xi2D_[:,:,1].size)
coef = 1. #(xi1D_[:,0]**2.)
mean = numpy.mean(xi2DList,axis=1)
plt.errorbar(xxx, coef*xi2D_[:,:,1].flatten(), marker='o', color='blue', label=r'$Data$')
plt.errorbar(xxx, coef*mean, marker='o', color='red', label=r'$mean \, over \, random$')

myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$a$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

### plot 2D
array2D = convert1DTo2D(mean)
fig = plt.figure()
ax = fig.add_subplot(111)

plt.imshow(array2D, origin='lower', interpolation='None')
cbar = plt.colorbar()
plt.grid(True)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()
plt.show()
### plot 2D
array2D = xi2D_[:,:,1]
fig = plt.figure()
ax = fig.add_subplot(111)

plt.imshow(array2D, origin='lower', interpolation='None')
cbar = plt.colorbar()
plt.grid(True)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()
plt.show()


xi1DList    = [ xi1DList ]
xi2DList    = [ xi2DList ]
covNameList = ['hello']
covList1D   = [ numpy.cov(el) for el in xi1DList ]
covList2D   = [ numpy.cov(el) for el in xi2DList ]

covFromFit1D = plotCovar(covList1D, covNameList)
covFromFit2D = plotCovar(covList2D, covNameList)

#### Pilot
'''
for ii in range(0,10):
	for jj in range(0,10):

		xi1DList    = [ numpy.load(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_1D.npy') ]
		xi2DList    = [ numpy.load(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_2D.npy') ]
		covNameList = ['Mock DLA: \, '+str(ii)+' / '+str(jj)]

		size = xi1DList[0][0].size
		print ii, jj, xi1DList[0][0].size, xi2DList[0][0].size, size

		covList1D   = [ numpy.cov(el)/size for el in xi1DList ]
		covList2D   = [ numpy.cov(el)/size for el in xi2DList ]

		covFromFit1D = plotCovar(covList1D, covNameList)
		covFromFit2D = plotCovar(covList2D, covNameList)

		numpy.save(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_fromFit_1D.npy',covFromFit1D)
		numpy.save(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_fromFit_2D.npy',covFromFit2D)
	
'''







