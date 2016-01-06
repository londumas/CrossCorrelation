# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import myTools
from const_delta import *

import subprocess
import time
import sys
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit

nbRegion__ = 0

### 1D
min1D__   = 0.
max1D__   = 200. #200.
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
nbBinM__ = 25;


path1__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/'

### Parameters
forest1__ = sys.argv[1]
forest2__ = sys.argv[2]
qso1__    = sys.argv[3]
qso2__    = sys.argv[4]
if (len(sys.argv)>5):
	wickIdx__ = int(sys.argv[5])
if (len(sys.argv)>6):
	box__   = sys.argv[6]
	simul__ = sys.argv[7]


def loadData(path1D, path2D):

	print path2D
	print path1D

	### 2D
	data = numpy.loadtxt(path2D)

	save0 = data[:,0]
	save1 = data[:,1]
	save2 = data[:,2]
	save3 = data[:,3]
	save4 = data[:,4]
	save5 = data[:,5]
	save6 = data[:,6]

	#print '  <z> = ', numpy.sum(save4)/numpy.sum(save5)

	tmp_save0  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save1  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save2  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save3  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save4  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save5  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save6  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))

	for i in range(0,len(save0)):
		iX = i/int(maxY2D__-minY2D__)
		iY = i%int(maxY2D__-minY2D__)

		idX = iX/int(binSize__)
		idY = iY/int(binSize__)

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]
		tmp_save2[idX][idY] += save2[i]
		tmp_save3[idX][idY] += save3[i]
		tmp_save4[idX][idY] += save4[i]
		tmp_save5[idX][idY] += save5[i]
		tmp_save6[idX][idY] += save6[i]

	xi2D = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__,3))

	xi2D[:,:,0] = numpy.sqrt( (tmp_save2/tmp_save5)**2. + (tmp_save3/tmp_save5)**2. )
	xi2D[:,:,1] = tmp_save0 / tmp_save5
	xi2D[:,:,2] = numpy.sqrt( (tmp_save1/tmp_save5 - xi2D[:,:,1]*xi2D[:,:,1])/tmp_save6 )
	cut = (tmp_save5 == 0.)
	xi2D[:,:,0][cut] = 0.
	xi2D[:,:,1][cut] = 0.
	xi2D[:,:,2][cut] = 0.

	### Mu
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

	binSizeX = int(max1D__)/nbBin1D__
	binSizeY = 100/nbBinM__
	
	for i in range(0,len(save0)):
		iX = i/100
		iY = i%100

		### for mu
		idX = iX/binSizeX
		idY = iY/binSizeY

		tmp_save0[idX][idY] += save0[i]
		tmp_save1[idX][idY] += save1[i]
		tmp_save2[idX][idY] += save2[i]
		tmp_save3[idX][idY] += save3[i]
		tmp_save4[idX][idY] += save4[i]
		tmp_save5[idX][idY] += save5[i]
		tmp_save6[idX][idY] += save6[i]
		
		
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
	
	xiMu = numpy.zeros(shape=(nbBin1D__,nbBinM__,4))
	xiMu[:,:,0] = tmp_save2/tmp_save5
	xiMu[:,:,1] = tmp_save3/tmp_save5
	xiMu[:,:,2] = tmp_save0/tmp_save5
	xiMu[:,:,3] = numpy.sqrt( (tmp_save1/tmp_save5 - xiMu[:,:,2]*xiMu[:,:,2])/tmp_save6)
	cut = (tmp_save5==0.)
	xiMu[:,:,0][cut] = 0.
	xiMu[:,:,1][cut] = 0.
	xiMu[:,:,2][cut] = 0.
	xiMu[:,:,3][cut] = 0.

	xiWe = numpy.zeros(shape=(nbBin1D__,3,3))
	xiWe[:,:,0] = tmp_save22/tmp_save55
	xiWe[:,:,1] = tmp_save00/tmp_save55
	xiWe[:,:,2] = numpy.sqrt( (tmp_save11/tmp_save55 - xiWe[:,:,1]*xiWe[:,:,1])/tmp_save66 )
	cut = (tmp_save55==0.)
	xiWe[:,:,0][cut] = 0.
	xiWe[:,:,1][cut] = 0.
	xiWe[:,:,2][cut] = 0.

	xi1D = numpy.zeros(shape=(nbBin1D__,3))
	xi1D[:,0] = tmp_save222/tmp_save555
	xi1D[:,1] = tmp_save000/tmp_save555
	xi1D[:,2] = numpy.sqrt( (tmp_save111/tmp_save555 - xi1D[:,1]*xi1D[:,1])/tmp_save666 )
	cut = (tmp_save555==0.)
	xi1D[:,0][cut] = 0.
	xi1D[:,1][cut] = 0.
	xi1D[:,2][cut] = 0.

	return xi1D, xi2D, xiMu, xiWe
def createIni2(inputFile, outputFile, folderName, dim):
	data = numpy.loadtxt(inputFile)

	save2 = data[:,2]
	save3 = data[:,3]
	save4 = data[:,4]
	save5 = data[:,5]
			
	meanRperp  = numpy.zeros(shape=(nbBinX2D__,2))
	meanRparal = numpy.zeros(shape=(nbBinY2D__,2))
	meanRedshift = numpy.sum(save4)/numpy.sum(save5)
			
	for i in range(0,len(save2)):
		iX = i/int(maxY2D__-minY2D__)
		iY = i%int(maxY2D__-minY2D__)
	
		idX = iX/int(binSize__)
		idY = iY/int(binSize__)
	
		meanRperp[idX][0]  += save2[i]
		meanRperp[idX][1]  += save5[i]
		meanRparal[idY][0] += save3[i]
		meanRparal[idY][1] += save5[i]
			
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
	stringRedshift = str(meanRedshift)
	print ' <z> = ', stringRedshift

	string = """
modelroot =  /home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/models/
fiducial  = DR9LyaMocksNLeffPk
nowiggles = DR9LyaMocksSB

cross-correlation = yes
anisotropic       = yes
decoupled         = yes

model-config = value[beta]=1.1;
model-config = fix[(1+beta)*bias]=-0.336;
model-config = fix[gamma-bias]=0.9; fix[gamma-beta]=0;
model-config = value[bias2]=3.64;
model-config = fix[beta2*bias2]=0.962524;
model-config = value[delta-v]=0;
model-config = value[BAO amplitude]=1;
model-config = fix[BAO alpha-iso]; value[BAO alpha-p*]=1;
model-config = fix[gamma-scale]=0;

model-config = binning[BAO alpha-parallel]={0.6:1.3}*50
model-config = binning[BAO alpha-perp]={0.6:1.3}*50

dist-add = rP,rT=0:2,-3:1
data     = """+folderName+"""bao"""+dim+"""

data-format = comoving-cartesian

axis1-bins = {""" + stringRparal   + """}
axis2-bins = {""" + stringRperp    + """}
axis3-bins = {""" + stringRedshift + """}

rmin = 40
rmax = 180

output-prefix = """+outputFile+"""bao"""+dim+""".

alt-config = fix[dist*]=0

ndump = 0
"""

		
	text_file = open(folderName+'bao'+dim+'.ini', "w")
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
def loadWick(path):

	### Get data
	data = numpy.loadtxt(path)

	### Get the Covariance
	var = numpy.zeros( shape=(nbBin1D__,3) )
	cov = numpy.zeros( shape=(nbBin1D__,nbBin1D__) )

	var[:,0] = data[:nbBin1D__,2]
	var[:,1] = data[:nbBin1D__,3]
	var[:,2] = numpy.sqrt(data[:nbBin1D__,4])

	data0 = data[nbBin1D__:,0].astype(int)
	data1 = data[nbBin1D__:,1].astype(int)
	data2 = data[nbBin1D__:,2]
	data3 = data[nbBin1D__:,3]
	data4 = data[nbBin1D__:,4]

	for i in range(0,data0.size):
		if (data0[i] == data1[i]):
			k = data0[i]
			cov[k][k] = ( data2[i]+var[k][0] )/( var[k][2]**2.*(data3[i]+var[k][1]) )
		else:
			cov[ data0[i] ][ data1[i] ] = data2[i]/( var[data0[i]][2]*var[data1[i]][2]*data3[i] )
			cov[ data1[i] ][ data0[i] ] = data2[i]/( var[data0[i]][2]*var[data1[i]][2]*data3[i] )
	
	#a = [(var[:,0]/(var[:,1]*var[:,2]**2.) )/(xi1D_[:,2]**2.),numpy.diag(cov)/(xi1D_[:,2]**2.),numpy.diag(cov2)/(xi1D_[:,2]**2.),numpy.diag(cov3)/(xi1D_[:,2]**2.)]
	#myTools.plot1D(a)

	return var, cov
def loadBootMap():
	'''
	'''

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/map.txt'
	
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
'''
def convert1DTo2D(array1D):
	
	#convert a 1D array to a 2D array
	

	array2D = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))

	for k in range(0,nbBin2D__):
		i = k/nbBinY2D__
		j = k%nbBinY2D__

		array2D[i][j] = array1D[k]

	return array2D
'''
def plotXi(rescale):

	xxx = xi1D_[:,0]
	yyy = xi1D_[:,1]
	yer = xi1D_[:,2]

	cut = (yer!=0.)
	xxx = xxx[ cut ]
	yyy = yyy[ cut ]
	yer = yer[ cut ]

	if (rescale==0):
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o')
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, fmt='o')
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, fmt='o')
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
def plotXi2D(rescale):


	xxx = numpy.transpose(xi2D_[:,:,0])
	yyy = numpy.transpose(xi2D_[:,:,1])
	yer = numpy.transpose(xi2D_[:,:,2])

	yyy[ (yer==0.) ]  = float('nan')
	yyy[ (xxx<=40.) ] = float('nan')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xticks([ 0.,50.,100.,150.,200.])
	ax.set_yticks([ -200.,-150.,-100.,-50.,0.,50.,100.,150.,200.])

	if (rescale==0):
		plt.imshow(yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__], interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi^{qf}(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		plt.imshow(xxx*yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|.\xi^{qf}(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
	if (rescale==2):
		plt.imshow(xxx**2.*yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|^{2}.\xi^{qf}(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
def plotMu(rescale):

	xxx = xiMu_[:,:,0]
	muu = xiMu_[:,:,1]
	yyy = xiMu_[:,:,2]
	yer = xiMu_[:,:,3]

	yyy[ (yer==0.) ] = float('nan')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	#ax.set_yticks([ 0.,50.,100.,150.])

	if (rescale==0):
		a = ''
		b = ''
	if (rescale==1):
		a = '|s|.'
		b = '[h^{-1}.Mpc]'
	if (rescale==2):
		a = '|s|^{2}.'
		b = '[(h^{-1}.Mpc)^{2}]'
	
	coef = numpy.power(xxx,rescale)
	plt.imshow(coef*yyy, origin='lower', interpolation='None',extent=[-1., 1.,min1D__, max1D__],aspect='auto')
	cbar = plt.colorbar()
	cbar.set_label(r'$'+a+'\\xi^{qf}(\, \overrightarrow{s} \,'+b+')$',size=40)
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$\mu$', fontsize=40)
	plt.ylabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
def plotWe(rescale):

	a = ['0.8 < |\mu|', '0.5 < |\mu| \leq 0.8', '|\mu| \leq 0.5']

	for i in range(0,3):

		###
		cut = (xiWe_[:,i,2]!=0.)

		if (xiWe_[:,i,0][cut].size==0):
			continue

		xxx = xiWe_[:,i,0][cut]
		yyy = xiWe_[:,i,1][cut]
		yer = xiWe_[:,i,2][cut]
		coef = numpy.power(xxx,rescale)
		

		plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$'+a[i]+'$')
	
		if (rescale==0):
			plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
			plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
		if (rescale==1):
			plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
			plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
		if (rescale==2):
			plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
			plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=2)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
def plotMultipol(xi):
	'''

	Plot the multipol of the cross-correlation

	'''

	### Plot or not
	plot = False

	### Get the data
	xxx = xi[:,:,0]
	muu = xi[:,:,1]
	yyy = xi[:,:,2]
	yer = xi[:,:,3]

	### Array to keep the results
	result_xi = numpy.zeros(shape=(nbBin1D__,5,3))
	for i in range(0,5):
		result_xi[:,i,0] = numpy.mean(xxx,axis=1)

	### Array with name of variable
	nameArray = [ 'xi0','xi1','xi2','xi3','xi4']

	for i in range(0,nbBin1D__):

		cut         = (muu[i,:]!=0.)
		tmpyyy      = yyy[i,:][cut]
		tmpyer      = yer[i,:][cut]
		xxxMu       = muu[i,:][cut]
		xxxMuPower1 = numpy.power(xxxMu,1.)
		xxxMuPower2 = numpy.power(xxxMu,2.)
		xxxMuPower3 = numpy.power(xxxMu,3.)
		xxxMuPower4 = numpy.power(xxxMu,4.)
		
		### Define the fit function
		def chi2(xi0,xi1,xi2,xi3,xi4):
			fit = (xi0 - 0.5*xi2 + 0.375*xi4) + (xi1-1.5*xi3)*xxxMuPower1 + (1.5*xi2 - 3.75*xi4)*xxxMuPower2 + 2.5*xi3*xxxMuPower3 + 4.375*xi4*xxxMuPower4
			return numpy.sum( numpy.power( (tmpyyy-fit)/tmpyer ,2.) )

		### Init ad perform the fit
		m = Minuit(chi2, xi0=0.,error_xi0=0.1,xi1=0.,error_xi1=0.1,xi2=0.,error_xi2=0.1,xi3=0.,error_xi3=0.1,xi4=0.,error_xi4=0.1, print_level=-1, errordef=0.01,fix_xi1=True,fix_xi2=False, fix_xi3=True,fix_xi4=True) 	
		m.migrad()
		#m.hesse()

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
		ylabel = ['\\xi^{qf} (|s|)','|s|.\\xi^{qf} (|s|) \, [h^{-1}.Mpc]','|s|^{2}.\\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]']

		### Show the result
		for i in range(0,3):

			for j in range(0,5):
				if ( result_xi[:,j,1][ (result_xi[:,j,1]!=0.) ].size == 0): continue

				tmp_xxx = xxx[ (result_xi[:,j,1]!=0.) ]
				tmp_yyy = result_xi[:,j,1][ (result_xi[:,j,1]!=0.) ]
				tmp_yer = result_xi[:,j,2][ (result_xi[:,j,2]!=0.) ]
				coef    = numpy.power(tmp_xxx,i)
				plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$\xi_{'+str(j)+'}$',color=color[j])

			cut = (xi1D_[:,2]!=0.)
			tmp_xxx = xi1D_[:,0][ cut ]
			tmp_yyy = xi1D_[:,1][ cut ]
			tmp_yer = xi1D_[:,2][ cut ]
			coef    = numpy.power(tmp_xxx,i)
			#plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$\xi$')

			plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
			plt.ylabel(r'$'+ylabel[i]+'$', fontsize=40)
			plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			myTools.deal_with_plot(False,False,True)
			plt.xlim([ numpy.min(tmp_xxx)-10., numpy.max(tmp_xxx)+10. ])
			plt.show()

	return result_xi
def fitCamb(data,pathToFile,mulpol=0):
	'''

	'''

	### Constants
	startFit   = 1.
	endFit     = 200.
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

	b1b2_init = -1. #numpy.mean(yyy[ (xxx<maxForGues) ]/yyy_Camb[ (xxx<maxForGues) ])
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
	
	m = Minuit(chi2, b1b2=b1b2_init,error_b1b2=0.1, roof=0.,error_roof=0.1,print_level=-1, errordef=0.01,fix_roof=False) 	
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
		'''
		plt.plot(xxx, yyy[0], label=r'$'+covNameList[c]+' \, \Delta \, s_{\perp, idx} = '+str(0)+'$', marker='o')

		plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
		plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
		plt.ylabel(r'$Cov( \Delta \, |s| \, bin \, idx  )$', fontsize=40)

		myTools.deal_with_plot(False,False,True)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.show()
		'''
		if (nbBin!=nbBin1D__):

			'''
			for i in range(1,nbBinX2D__):
				plt.plot(xxx, yyy[i], label=r'$'+covNameList[c]+' \, \Delta \, s_{\perp, idx} = '+str(i)+'$', marker='o')
			
			plt.title(r'$Variogram  \, \delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
			plt.xlabel(r'$\Delta \, |s|  \, bin \, idx $', fontsize=40)
			plt.ylabel(r'$Cov( \Delta \, |s| \, bin \, idx  )$', fontsize=40)

			myTools.deal_with_plot(False,False,True)
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.show()
			'''
			yyy2 = numpy.nan_to_num(yyy2[:,:,0]/yyy2[:,:,1])

			### 0
			'''
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
			'''
	### Show all the method on the main variogram
	for c in range(0,nbCov):
		if (nbBin==nbBin1D__):
			xTitle = '|\Delta \, |s_{1}|-|s_{2}| | \, [h^{-1}.Mpc]'
		else:
			xTitle = '|\Delta \, s_{\parallel,1}-s_{\parallel,2} | \, [h^{-1}.Mpc] \, for \, |\Delta \, s_{\perp}| = 0'
		plt.plot(tmp_xxx[c], tmp_yyy[c], label=r'$'+covNameList[c]+'$', marker='o')

	plt.xlabel(r'$'+xTitle+'$', fontsize=40)
	plt.ylabel(r'$Cov( s_{1}, s_{2})$', fontsize=40)
	plt.xlim([ numpy.min(tmp_xxx[0])-10., numpy.max(tmp_xxx[0])+10. ])
	myTools.deal_with_plot(False,False,True)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
	
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
			#plt.plot(tmp_xxx[c], minuit_yyy, label=r'$model$', marker='o')

			### Save the results of the fit
			

		plt.xlabel(r'$'+xTitle+'$', fontsize=40)
		plt.ylabel(r'$Cor( s_{1}, s_{2})$', fontsize=40)
		plt.xlim([ numpy.min(tmp_xxx[0])-10., numpy.max(tmp_xxx[0])+10. ])
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
		'''
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
		'''
	return covFromFitList

def saveListReal(nbReal,path1D,path2D,pathToSave,subsampling=False):
	'''

	'''

	print path1D

	list1D       = numpy.zeros( shape=(nbBin1D__,nbReal) )
	list2D       = numpy.zeros( shape=(nbBin2D__,nbReal) )
	listMu       = numpy.zeros( shape=(nbBin1D__*nbBinM__,nbReal) )
	listWe       = numpy.zeros( shape=(nbBin1D__,3,nbReal) )
	listMultipol = numpy.zeros( shape=(nbBin1D__,5,nbReal) )

	for i in numpy.arange(nbReal):

		xi1D, xi2D, xiMu, xiWe = loadData(path1D+str(i)+'.txt',path2D+str(i)+'.txt')
		list1D[:,i]         = xi1D[:,1]
		list2D[:,i]         = xi2D[:,:,1].flatten()
		listMu[:,i]         = xiMu[:,:,2].flatten()
		listWe[:,:,i]       = xiWe[:,:,1]
		listMultipol[:,:,i] = plotMultipol(xiMu)[:,:,1]


	numpy.save(pathToSave+'1D',list1D)
	numpy.save(pathToSave+'2D',list2D)
	numpy.save(pathToSave+'Mu',listMu)
	numpy.save(pathToSave+'We',listWe)
	numpy.save(pathToSave+'Multipol',listMultipol)

	cov1D = numpy.cov(list1D)
	cov2D = numpy.cov(list2D)
	covMu = numpy.cov(listMu)
	if (subsampling):
		cov1D /= nbReal
		cov2D /= nbReal
		covMu /= nbReal

	numpy.save(pathToSave+'cov_1D',cov1D)
	numpy.save(pathToSave+'cov_2D',cov2D)
	numpy.save(pathToSave+'cov_Mu',covMu)

	return
def saveListRealMocks(ni,nj):
	'''

	'''

	path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
	pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_delta_QSO_1Mpc_result_'

	list1D       = numpy.zeros( shape=(nbBin1D__,ni*nj) )
	list2D       = numpy.zeros( shape=(nbBin2D__,ni*nj) )
	listMu       = numpy.zeros( shape=(nbBin1D__*nbBinM__,ni*nj) )
	listWe       = numpy.zeros( shape=(nbBin1D__,3,ni*nj) )
	listMultipol = numpy.zeros( shape=(nbBin1D__,5,ni*nj) )

	for i in numpy.arange(ni):
		for j in numpy.arange(nj):
			tmpPath =  path + str(i)+'/Simu_00'+str(j)+'/Results/xi_delta_QSO_'
			xi1D, xi2D, xiMu, xiWe = loadData(tmpPath+'Mu_LYA_QSO.txt',tmpPath+'2D_LYA_QSO.txt')
			list1D[:,i*10+j]         = xi1D[:,1]
			list2D[:,i*10+j]         = xi2D[:,:,1].flatten()
			listMu[:,i*10+j]         = xiMu[:,:,2].flatten()
			listWe[:,:,i*10+j]       = xiWe[:,:,1]
			listMultipol[:,:,i*10+j] = plotMultipol(xiMu)[:,:,1]

	numpy.save(pathToSave+'1D',list1D)
	numpy.save(pathToSave+'2D',list2D)
	numpy.save(pathToSave+'Mu',listMu)
	numpy.save(pathToSave+'We',listWe)
	numpy.save(pathToSave+'Multipol',listMultipol)

	'''
	cov1D = numpy.cov(list1D)
	cov2D = numpy.cov(list2D)
	covMu = numpy.cov(listMu)
	numpy.save(pathToSave+'cov_1D',cov1D)
	numpy.save(pathToSave+'cov_2D',cov2D)
	numpy.save(pathToSave+'cov_Mu',covMu)
	'''

	return
def prepareForBAOFIT():

	i = box__
	j = simul__


	### Create the .ini file
	path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_delta_QSO_2D_LYA_QSO.txt'
	path2 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f/'
	path3 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f/'  ###_fixedBAO
	createIni2(path,path3,path2,'2D')

	'''
	### Get the correlation
	path    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_delta_QSO_'
	xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path+'Mu_LYA_QSO.txt',path+'2D_LYA_QSO.txt')

	### Get the covariance matrix
	cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/subSampling_LYA_QSO_result_cov_2D.npy')

	### If correlation matrix from another matrix
	cor = myTools.getCorrelationMatrix(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_delta_QSO_result_cov_2D_meanSubSampling.npy'))
	cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))

	### Ready to save data
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
		tmp_xi2D[k11] = xi2D_[i1,j1,1]
		
		for k2 in range(0,k1+1):
			i2       = k2/nbBinY2D__
			j2       = k2%nbBinY2D__
			k22      = j2*nbBinX2D__ + i2
			#if (k1==k2): tmp_xi2DCov[idx] = cov[k1][k2]
			tmp_xi2DCov[idx] = cov[k1][k2]
			tmp_xi2DCovIdx[idx]  = k11
			tmp_xi2DCovIdx2[idx] = k22
			idx += 1
	
	### Save data
	numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f/bao2D.data',zip(tmp_idx,tmp_xi2D),fmt='%u %1.20e')
	numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/BaoFit_q_f/bao2D.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%u %u %1.20e')
	'''	

	### Send the fit
	command = "/home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/build/baofit -i " + path2 + "bao2D.ini"
	print command
	subprocess.call(command, shell=True)

#prepareForBAOFIT()



'''
path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
sumWeight = numpy.zeros(100)
sumOne    = numpy.zeros(100)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		idx = i*10+j
		data = numpy.loadtxt( path + str(i)+'/Simu_00'+str(j)+'/Results/xi_delta_QSO_2D_LYA_QSO.txt' )
		sumWeight[idx] = numpy.sum(data[:,5])
		sumOne[idx]    = numpy.sum(data[:,6])
'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/'
name = path
path2D = name+'xi_delta_QSO_2D_LYA_QSO_bootstrap_'
path1D = name+'xi_delta_QSO_Mu_LYA_QSO_bootstrap_'

nbReal = 2
sumWeight = numpy.zeros( nbReal )
sumOne    = numpy.zeros( nbReal )
list1D       = numpy.zeros( shape=(nbBin1D__,nbReal) )
list2D       = numpy.zeros( shape=(nbBin2D__,nbReal) )
listMu       = numpy.zeros( shape=(nbBin1D__*nbBinM__,nbReal) )
listWe       = numpy.zeros( shape=(nbBin1D__,3,nbReal) )
listMultipol = numpy.zeros( shape=(nbBin1D__,5,nbReal) )

for i in numpy.arange(nbReal):
	xi1D, xi2D, xiMu, xiWe = loadData(path1D+str(i)+'.txt',path2D+str(i)+'.txt')
	list1D[:,i]         = xi1D[:,1]
	list2D[:,i]         = xi2D[:,:,1].flatten()
	listMu[:,i]         = xiMu[:,:,2].flatten()
	listWe[:,:,i]       = xiWe[:,:,1]
	listMultipol[:,:,i] = plotMultipol(xiMu)[:,:,1]

	data = numpy.loadtxt(path2D+str(i)+'.txt')
	sumWeight[i] = numpy.sum(data[:,5])
	sumOne[i]    = numpy.sum(data[:,6])


### Get the mean
meanXi = numpy.mean(list2D,axis=2)



cov1D = numpy.cov(list1D,aweights=sumWeight)
cov2D = numpy.cov(list2D,aweights=sumWeight)
covMu = numpy.cov(listMu,aweights=sumWeight)
if (subsampling):
	cov1D /= nbReal
	cov2D /= nbReal
	covMu /= nbReal

plotCovar([cov1D], ['a'],nbBinX=50,nbBinY=100)


cov1D = numpy.cov(list1D,aweights=sumOne)
cov2D = numpy.cov(list2D,aweights=sumOne)
covMu = numpy.cov(listMu,aweights=sumOne)
if (subsampling):
	cov1D /= nbReal
	cov2D /= nbReal
	covMu /= nbReal

plotCovar([cov1D], ['a'],nbBinX=50,nbBinY=100)






























