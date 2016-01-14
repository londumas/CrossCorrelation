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


path1__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_test_PDFMocksJMC/'

### Parameters
forest1__ = sys.argv[1]
forest2__ = sys.argv[2]
qso1__    = sys.argv[3]
qso2__    = sys.argv[4]
if (len(sys.argv)>5):
	wickIdx__ = int(sys.argv[5])
if (len(sys.argv)>6):
	box__   = int(sys.argv[6])
	simul__ = int(sys.argv[7])


def loadData(path1D, path2D,selection=0):
	'''

	selection:
		- 0: do all
		- 1: only 1D
		- 2: only 2D

	'''

	print path2D
	print path1D

	xi2D = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__,3))
	xiMu = numpy.zeros(shape=(nbBin1D__,nbBinM__,4))
	xiWe = numpy.zeros(shape=(nbBin1D__,3,3))
	xi1D = numpy.zeros(shape=(nbBin1D__,3))

	if (selection==0 or selection==1):
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
	
		xiMu[:,:,0] = tmp_save2/tmp_save5
		xiMu[:,:,1] = tmp_save3/tmp_save5
		xiMu[:,:,2] = tmp_save0/tmp_save5
		xiMu[:,:,3] = numpy.sqrt( (tmp_save1/tmp_save5 - xiMu[:,:,2]*xiMu[:,:,2])/tmp_save6)
		cut = (tmp_save5==0.)
		xiMu[:,:,0][cut] = 0.
		xiMu[:,:,1][cut] = 0.
		xiMu[:,:,2][cut] = 0.
		xiMu[:,:,3][cut] = 0.

		xiWe[:,:,0] = tmp_save22/tmp_save55
		xiWe[:,:,1] = tmp_save00/tmp_save55
		xiWe[:,:,2] = numpy.sqrt( (tmp_save11/tmp_save55 - xiWe[:,:,1]*xiWe[:,:,1])/tmp_save66 )
		cut = (tmp_save55==0.)
		xiWe[:,:,0][cut] = 0.
		xiWe[:,:,1][cut] = 0.
		xiWe[:,:,2][cut] = 0.

		xi1D[:,0] = tmp_save222/tmp_save555
		xi1D[:,1] = tmp_save000/tmp_save555
		xi1D[:,2] = numpy.sqrt( (tmp_save111/tmp_save555 - xi1D[:,1]*xi1D[:,1])/tmp_save666 )
		cut = (tmp_save555==0.)
		xi1D[:,0][cut] = 0.
		xi1D[:,1][cut] = 0.
		xi1D[:,2][cut] = 0.

	if (selection==0 or selection==2):
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

		xi2D[:,:,0] = numpy.sqrt( (tmp_save2/tmp_save5)**2. + (tmp_save3/tmp_save5)**2. )
		xi2D[:,:,1] = tmp_save0 / tmp_save5
		xi2D[:,:,2] = numpy.sqrt( (tmp_save1/tmp_save5 - xi2D[:,:,1]*xi2D[:,:,1])/tmp_save6 )
		cut = (tmp_save5 == 0.)
		xi2D[:,:,0][cut] = 0.
		xi2D[:,:,1][cut] = 0.
		xi2D[:,:,2][cut] = 0.

	return xi1D, xi2D, xiMu, xiWe
def createIni2(inputFile, pathToOutFit, pathToData,pathToIni, dim, param=None):

	### Get interval for alpha scan
	minAlpha = (param[11] - 0.5*param[11]).astype('str')
	maxAlpha = (param[11] + 0.5*param[11]).astype('str')
	param = param.astype('str')

	data = numpy.loadtxt(inputFile)

	tmp_save2  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save3  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save4  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	tmp_save5  = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__))
	meanRperp  = numpy.zeros(shape=(nbBinX2D__,2))
	meanRparal = numpy.zeros(shape=(nbBinY2D__,2))
	meanRedshift = numpy.sum(data[:,4])/numpy.sum(data[:,5])
			
	for i in range(0,data[:,2].size):
		
		iX = i/int(maxY2D__-minY2D__)
		iY = i%int(maxY2D__-minY2D__)
	
		idX = iX/int(binSize__)
		idY = iY/int(binSize__)
	
		meanRperp[idX][0]  += data[i,2]
		meanRperp[idX][1]  += data[i,5]
		meanRparal[idY][0] += data[i,3]
		meanRparal[idY][1] += data[i,5]
		
		tmp_save2[idX][idY] += data[i,2]
		tmp_save3[idX][idY] += data[i,3]
		tmp_save4[idX][idY] += data[i,4]
		tmp_save5[idX][idY] += data[i,5]
		
	### Get the grid
	rParal = tmp_save3 / tmp_save5
	rPerp  = tmp_save2 / tmp_save5
	z      = tmp_save4 / tmp_save5

	grid = numpy.zeros( shape=(nbBin2D__,4) )
	indexMatrix = numpy.arange(nbBin2D__)
	grid[:,0] = (indexMatrix%nbBinY2D__)*nbBinX2D__ + indexMatrix/nbBinY2D__
	grid[:,1] = rParal.flatten()
	grid[:,2] = rPerp.flatten()
	grid[:,3] = z.flatten()
	numpy.savetxt(pathToData+'bao'+dim+'.grid',zip(grid[:,0],grid[:,1],grid[:,2],grid[:,3]),fmt='%u %1.20e %1.20e %1.20e')
		
	### Get the s_perp bin center
	meanRperp[:,0]  /= meanRperp[:,1]
	stringRperp = ''
	for el in meanRperp[:-1,0]:
		stringRperp += str(el) + ','
	stringRperp += str(meanRperp[-1,0])
			
	### Get the s_paral bin center
	meanRparal[:,0] /= meanRparal[:,1]
	stringRparal = ''
	for el in meanRparal[:-1,0]:
		stringRparal += str(el) + ','
	stringRparal += str(meanRparal[-1,0])

	### Get the redshift center
	stringRedshift = str(meanRedshift)
	print ' <z> = ', stringRedshift

	string = """

## Linear theory P(k) templates with and w/o wiggles
modelroot = /home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/models/
fiducial =  DR9LyaMocksLCDM
nowiggles = DR9LyaMocksLCDMSB

## k-space fit
kspace = true
ell-max = 6

# Model configuration
cross-correlation = yes
anisotropic = yes
decoupled   = yes
custom-grid = yes
pixelize = yes
dist-matrix = yes
dist-matrix-order = 5000

# Parameter setup
model-config = value[beta]=               """+param[0] +""";
model-config = fix[(1+beta)*bias]=        """+param[1] +""";
model-config = fix[gamma-bias]=           """+param[2] +""";
model-config = fix[gamma-beta]=           """+param[3] +""";
model-config = value[delta-v]=            """+param[4] +""";
model-config = value[bias2]=              """+param[5] +""";
model-config = fix[beta2*bias2]=          """+param[6] +""";
model-config = fix[1+f]=                  """+param[7] +""";
model-config = fix[SigmaNL-perp]=         """+param[8] +""";
model-config = fix[BAO amplitude]=        """+param[9] +""";
model-config = fix[BAO alpha-iso]=        """+param[10] +""";
model-config = value[BAO alpha-parallel]= """+param[11] +""";
model-config = value[BAO alpha-perp]=     """+param[12]+""";
model-config = fix[gamma-scale]=          """+param[13]+""";
model-config = fix[pixel scale]=0.7;

## 2D chisq scan in BAO parameters
model-config = binning[BAO alpha-parallel] ={"""+minAlpha+""":"""+maxAlpha+"""}*50
model-config = binning[BAO alpha-perp]     ={"""+minAlpha+""":"""+maxAlpha+"""}*50	

## Reference redshift
zref = 2.3

## Maximum allowed radial dilation (increases the range that model needs to cover)
dilmin = 0.5
dilmax = 2.5

## Non-linear broadening with 1+f = (SigmaNL-par)/(SigmaNL-perp)

# boxprior keeps result positive (since model only depends on squared value)
model-config = boxprior[SigmaNL-perp] @ (0,6);
# un-comment next line to broaden all scales (default is peak only)
#nl-broadband = true

# Broadband distortion model
dist-add = rP,rT=0:2,-3:1

### Data Options #############################################

## Data to analyze
data = """+pathToData+"""bao"""+dim+"""
dist-matrix-name = """+pathToData+"""bao"""+dim+"""

## Data format
data-format = comoving-cartesian
axis1-bins = {""" + stringRparal   + """}
axis2-bins = {""" + stringRperp    + """}
axis3-bins = {""" + stringRedshift + """}

### Analysis Options #########################################

# Cuts to apply before fitting
rmin = 40
rmax = 180

# Generate a second set of outputs with the additive distortion turned off
alt-config = fix[dist*]=0

# Do not dump multipoles (since the distortion model multipole integrals are singular)
ndump = 0

# Prefix to use for all analysis output files
output-prefix = """+pathToOutFit+"""bao"""+dim+""".
"""

	print '  save into ', pathToIni+'bao'+dim+'.ini'
	text_file = open(pathToIni+'bao'+dim+'.ini', "w")
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
	#yyy[ (xxx<=40.) ] = float('nan')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xticks([ 0.,50.,100.,150.,200.])
	ax.set_yticks([ -200.,-150.,-100.,-50.,0.,50.,100.,150.,200.])
	extent=[minX2D__, maxX2D__, maxY2D__,minY2D__]

	coef = numpy.power(xxx,rescale)
	plt.imshow(coef*yyy, origin='upper',extent=extent, interpolation='None')
	cbar = plt.colorbar()

	if (rescale==0):
		cbar.set_label(r'$\xi^{qf}(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		cbar.set_label(r'$|s|.\xi^{qf}(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
	if (rescale==2):
		cbar.set_label(r'$|s|^{2}.\xi^{qf}(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

	#plt.plot( [0.,200.],[0.,4*200.],color='white',linewidth=2 )
	#plt.plot( [0.,200.],[0.,-4*200.],color='white',linewidth=2 )
	#plt.plot( [0.,200.],[0.,200.],color='white',linewidth=2 )
	#plt.plot( [0.,200.],[0.,-200.],color='white',linewidth=2 )

	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()

	plt.xlim([ numpy.min(minX2D__), numpy.max(maxX2D__) ])
	plt.ylim([ numpy.min(maxY2D__), numpy.max(minY2D__) ])

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
	startFit   = 20.
	endFit     = 200.
	maxForGues = 200.
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

	b1b2_init = -0.2 #numpy.mean(yyy[ (xxx<maxForGues) ]/yyy_Camb[ (xxx<maxForGues) ])
	print '  b1b2_init = ', b1b2_init

	
	### Show result
	for i in range(0,3):

		coef = numpy.power(xxx,i)
		plt.errorbar(xxx,coef*yyy,yerr=coef*yer,fmt='o')

		coef = numpy.power(xxx,i)
		plt.plot(xxx,coef*b1b2_init*yyy_Camb,marker='o')
		myTools.deal_with_plot(False,False,False)
		plt.show()
	

	### Define the fit function
	def chi2(b1b2,roof):
		fit = yyy_Camb*b1b2+roof
		return numpy.sum( numpy.power( (yyy-fit)/yer ,2.) )

	### Init and perform the fit
	
	m = Minuit(chi2, b1b2=b1b2_init,error_b1b2=0.1, roof=0.,error_roof=0.1,print_level=-1, errordef=0.01,fix_roof=True) 	
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
def saveListReal(nbReal,path1D,path2D,pathToSave,subsampling=False):
	'''
		nbReal      = number of realisation
		path1D      = path to the xi_Mu correlation
		path2D      = path to the xi_2D correlation
		pathToSave  = path to folder where to save
		subsampling = is it a subsampling?

		### Data
		path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/'
		saveListReal(100,path+'xi_delta_QSO_Mu_LYA_QSO_shuffleQSO_',path+'xi_delta_QSO_2D_LYA_QSO_shuffleQSO_',path+'shuffleQSO_LYA_QSO_',False)

		### Simulation
		path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00' + str(box__)+'/Simu_00'+str(simul__)+'/Results_RandomPosInCell/'
		saveListReal(80,path+'xi_delta_QSO_Mu_LYA_QSO_subsampling_',path+'xi_delta_QSO_2D_LYA_QSO_subsampling_',path+'subSampling_LYA_QSO_',True)

		cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/subsampling_LYA_QSO_cov_1D.npy')
		cor = myTools.getCorrelationMatrix(cov)
		myTools.plot2D(cor)

	'''

	print path1D

	list1D       = numpy.zeros( shape=(nbBin1D__,nbReal) )
	list2D       = numpy.zeros( shape=(nbBin2D__,nbReal) )
	listMu       = numpy.zeros( shape=(nbBin1D__*nbBinM__,nbReal) )
	listWe       = numpy.zeros( shape=(nbBin1D__,3,nbReal) )
	listMultipol = numpy.zeros( shape=(nbBin1D__,5,nbReal) )

	for i in numpy.arange(nbReal):

		try:
			xi1D, xi2D, xiMu, xiWe = loadData(path1D+str(i)+'.txt',path2D+str(i)+'.txt')
			list1D[:,i]         = xi1D[:,1]
			list2D[:,i]         = xi2D[:,:,1].flatten()
			listMu[:,i]         = xiMu[:,:,2].flatten()
			listWe[:,:,i]       = xiWe[:,:,1]
			listMultipol[:,:,i] = plotMultipol(xiMu)[:,:,1]
		except:
			print '   ERROR:: ', path1D+str(i)+'.txt',path2D+str(i)+'.txt'
			commandProd = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Correlation/bin/main.exe"
			tmp_command = "clubatch \"time ; hostname ; " + commandProd + ' 5 0 ' + str(i) + ' 0 0 0 ' + "\""
			subprocess.call(tmp_command, shell=True)

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
def saveListRealMocks(ni,nj,distortion=False):
	'''
		- ni: Box
		- nj: Simu
		- distortion: Flag to use or not the distortion matrix (defalut=False)

		Usage example:
			cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_NicolasDistortion/xi_delta_QSO_result_cov_1D.npy')
			cor = myTools.getCorrelationMatrix(cov)
			myTools.plot2D(cor)
			a = myTools.plotCovar([cor],['a'])

	'''

	### Where to get correlation
	path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
	### Where to save the results
	pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_NicolasDistortionWithDistortion/xi_delta_QSO_result_'

	list1D       = numpy.zeros( shape=(nbBin1D__,ni*nj) )
	list2D       = numpy.zeros( shape=(nbBin2D__,ni*nj) )
	listMu       = numpy.zeros( shape=(nbBin1D__*nbBinM__,ni*nj) )
	listWe       = numpy.zeros( shape=(nbBin1D__,3,ni*nj) )
	listMultipol = numpy.zeros( shape=(nbBin1D__,5,ni*nj) )

	for i in numpy.arange(ni):
		for j in numpy.arange(nj):
			tmpPath =  path + str(i)+'/Simu_00'+str(j)+'/Results_NicolasDistortion/xi_delta_QSO_'

			try:
				xi1D, xi2D, xiMu, xiWe = loadData(tmpPath+'Mu_LYA_QSO.txt',tmpPath+'2D_LYA_QSO.txt')
				list1D[:,i*10+j]         = xi1D[:,1]
				list2D[:,i*10+j]         = xi2D[:,:,1].flatten()
				listMu[:,i*10+j]         = xiMu[:,:,2].flatten()
				listWe[:,:,i*10+j]       = xiWe[:,:,1]
				listMultipol[:,:,i*10+j] = plotMultipol(xiMu)[:,:,1]
			except:
				print '   ERROR:: ', tmpPath

			if (distortion):
				tmp_command = " echo " + str(i) + " " + str(j)
				subprocess.call(tmp_command, shell=True)

				### distortion matrix 1D
				data = numpy.loadtxt(tmpPath+'distortionMatrix_1D_LYA_QSO.txt')
				list1D[:,i*10+j] = numpy.dot(data,xi1D[:,1])
				#myTools.plot2D(data)
				#myTools.plot1D([xi1D[:,1],list1D[:,i*10+j]],'-','-','-',['before','after'])
				### distortion matrix 2D
				data = numpy.loadtxt(tmpPath+'distortionMatrix_2D_LYA_QSO.txt')
				xi1D = myTools.convert2DTo1D(xi2D[:,:,1], nbBinX2D__,nbBinY2D__)
				list2D[:,i*10+j] = numpy.dot(data,xi1D)

				

	numpy.save(pathToSave+'1D',list1D)
	numpy.save(pathToSave+'2D',list2D)
	numpy.save(pathToSave+'Mu',listMu)
	numpy.save(pathToSave+'We',listWe)
	numpy.save(pathToSave+'Multipol',listMultipol)

	
	cov1D = numpy.cov(list1D)
	cov2D = numpy.cov(list2D)
	covMu = numpy.cov(listMu)
	numpy.save(pathToSave+'cov_1D',cov1D)
	numpy.save(pathToSave+'cov_2D',cov2D)
	numpy.save(pathToSave+'cov_Mu',covMu)
	

	return
def saveMeanCorMatrix(ni,nj):
	'''

	Save the mean correlation matrix of the 100 correaltion matrix

	'''

	path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_00'
	pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Results/xi_delta_QSO_result_'

	### Get the 0_0 simu
	cor1D = myTools.getCorrelationMatrix(numpy.load(path + '0/Simu_000/Results/subSampling_LYA_QSO_cov_1D.npy'))
	cor2D = myTools.getCorrelationMatrix(numpy.load(path + '0/Simu_000/Results/subSampling_LYA_QSO_cov_2D.npy'))

	for i in numpy.arange(ni):
		for j in numpy.arange(nj):

			if (i==0 and j==0): continue

			cor1D += myTools.getCorrelationMatrix(numpy.load(path + str(i)+'/Simu_00'+str(j)+'/Results/subSampling_LYA_QSO_cov_1D.npy'))
			cor2D += myTools.getCorrelationMatrix(numpy.load(path + str(i)+'/Simu_00'+str(j)+'/Results/subSampling_LYA_QSO_cov_2D.npy'))
	cor1D /= 100.
	cor2D /= 100.

	numpy.save(pathToSave+'cor_mean_1D', cor1D)
	numpy.save(pathToSave+'cor_mean_2D', cor2D)

	return
def prepareForBAOFIT():

	### Constants
	param = numpy.asarray( [1.6,-0.336,0.9,0.,0.,3.25,0.962524,3.26,1.966,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.] )
	doBootstraps = False
	nbRegions    = 80

	'''
	### Create the .ini file (for simulation)
	i = box__
	j = simul__
	#param = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Annalyse_BAOFIT/param.npy')
	#param = param[:,0,i*10+j]
	rowPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
	path   = rowPath + '/xi_delta_QSO_'
	cov = numpy.load(rowPath + '/subSampling_LYA_QSO_cov_2D.npy')
	cor = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Results/xi_delta_QSO_result_cor_meanFromFit_2D.npy')
	inputFile    = rowPath + 'xi_delta_QSO_2D_LYA_QSO.txt'
	pathToOutFit = rowPath + 'BaoFit_q_f/'
	pathToData   = rowPath + 'BaoFit_q_f/'
	pathToIni    = rowPath + 'BaoFit_q_f/'
	createIni2(inputFile, pathToOutFit, pathToData,pathToIni,'2D',param)
	'''
	
	### For data
	rowPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_DR12_nicolas/'
	path    = rowPath + '/xi_delta_QSO_'
	cov     = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_test_PDFMocksJMC_meanLambda_testNoCap/shuffleForest_LYA_QSO_cov_2D.npy')
	cor     = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy')
	#inputFile    = rowPath + '/xi_delta_QSO_2D_LYA_QSO.txt'
	inputFile    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_DR12_nicolas/xi_delta_QSO_2D_LYA_QSO.txt'
	pathToDistortionMatrix = rowPath + '/xi_delta_QSO_distortionMatrix_2D_LYA_QSO.txt'
	pathToOutFit = rowPath + 'BaoFit_q_f/'
	pathToData   = rowPath + 'BaoFit_q_f/'
	pathToIni    = rowPath + 'BaoFit_q_f/'

	if (doBootstraps):
		pathToOutFit += 'Bootstraps/boot_'+str(wickIdx__).zfill(4)+'/'
		pathToData   += 'Bootstraps/boot_'+str(wickIdx__).zfill(4)+'/'
		pathToIni    += 'Bootstraps/boot_'+str(wickIdx__).zfill(4)+'/'
		subprocess.call('mkdir ' + pathToIni, shell=True)
		#subprocess.call('cp '+rowPath+'/BaoFit_q_f/bao2D.grid ' + pathToIni, shell=True)
	
		### Set the seed
		numpy.random.seed(seed=42)
		array = numpy.random.choice(numpy.arange(1000000).astype(int), size=10000, replace=False)
		numpy.random.seed(seed=array[wickIdx__])
		randomList = numpy.random.randint( 0, high=nbRegions, size=nbRegions)
		print array[wickIdx__]

		subsampling = numpy.load(rowPath+'/subsampling_LYA_QSO_2D.npy')
		bootstrap   = numpy.zeros( shape=(nbBin2D__,nbRegions) )
	
		xi2D1D = numpy.zeros(nbBin2D__)
		for i in randomList:
			bootstrap[:,i] = subsampling[:,i]
		xi2D1D = numpy.mean( bootstrap, axis=1)
		cov = numpy.cov(bootstrap)/nbRegions
		xi2D_ = numpy.zeros(shape=(nbBinX2D__,nbBinY2D__,3))
		xi2D_[:,:,1] = myTools.convert1DTo2D(xi2D1D,nbBinX2D__,nbBinY2D__)
	else:
		### Get the correlation (data)
		xi1D_, xi2D_, xiMu_, xiWe_ = loadData('',inputFile,2)

	### Create the 'bao2D.grid' file
	createIni2(inputFile, pathToOutFit, pathToData,pathToIni,'2D',param)

	### Correlation
	correlation = numpy.zeros( shape=(nbBin2D__,2) )
	indexMatrix = numpy.arange(nbBin2D__)
	correlation[:,0] = (indexMatrix%nbBinY2D__)*nbBinX2D__ + indexMatrix/nbBinY2D__
	correlation[:,1] = xi2D_[:,:,1].flatten()
	cutCorrelation = (correlation[:,1]!=0.)
	numpy.savetxt( pathToData + '/bao2D.data',zip(correlation[:,0][cutCorrelation],correlation[:,1][cutCorrelation]),fmt='%u %1.20e')

	'''
	### Covariance matrix
	cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))
	covarianceMatrix = numpy.zeros( shape=(nbBin2D__*nbBin2D__,3) )
	indexMatrix1 = numpy.arange(nbBin2D__*nbBin2D__).reshape(nbBin2D__,nbBin2D__)/nbBin2D__
	indexMatrix2 = numpy.arange(nbBin2D__*nbBin2D__).reshape(nbBin2D__,nbBin2D__)%nbBin2D__
	indexMatrix1 = (indexMatrix1%nbBinY2D__)*nbBinX2D__ + indexMatrix1/nbBinY2D__
	indexMatrix2 = (indexMatrix2%nbBinY2D__)*nbBinX2D__ + indexMatrix2/nbBinY2D__
	covarianceMatrix[:,0] = numpy.triu(indexMatrix1,k=0).flatten()
	covarianceMatrix[:,1] = numpy.triu(indexMatrix2,k=0).flatten()
	covarianceMatrix[:,2] = numpy.triu(cov,k=0).flatten()
	cutCovarMatrix = (covarianceMatrix[:,2]!=0.)
	if (covarianceMatrix[:,2][cutCovarMatrix].size != int(nbBin2D__*(nbBin2D__+1.)/2.) ):
		print '  xi_delta_QSO.py::prepareForBAOFIT()  size of covariance matrix is incorrect'
		print '  size covariance matrix', cutCovarMatrix[:,2][cutCovarMatrix].size
		print '  size it should have', nbBin2D__*(nbBin2D__+1.)/2.
		return
	numpy.savetxt( pathToData + '/bao2D.cov',zip(covarianceMatrix[:,0][cutCovarMatrix],covarianceMatrix[:,1][cutCovarMatrix],covarianceMatrix[:,2][cutCovarMatrix]),fmt='%u %u %1.20e')

	### Distortion matrix
	dmatData = numpy.loadtxt(pathToDistortionMatrix)
	distortionMatrix = numpy.zeros( shape=(nbBin2D__*nbBin2D__,3) )
	distortionMatrix[:,0] = indexMatrix1.flatten()
	distortionMatrix[:,1] = indexMatrix2.flatten()
	distortionMatrix[:,2] = dmatData.flatten()
	cutDistortionMatrix = (distortionMatrix[:,2]!=0.)
	numpy.savetxt( pathToData + '/bao2D.dmat',zip(distortionMatrix[:,0][cutDistortionMatrix], distortionMatrix[:,1][cutDistortionMatrix], distortionMatrix[:,2][cutDistortionMatrix]),fmt='%u %u %1.20e')
	'''

	### Send the fit
	command = '/home/gpfs/manip/mnt0607/bao/hdumasde/Program/deepzot/bin/baofit -i ' + pathToIni + 'bao2D.ini' #--parameter-scan'   ### --toymc-samples 10000
	print command
	subprocess.call(command, shell=True)
	
	if (doBootstraps):
		subprocess.call('rm '+ pathToData + '/bao2D.ini', shell=True)
		subprocess.call('rm '+ pathToData + '/bao2D.grid', shell=True)
		subprocess.call('rm '+ pathToData + '/bao2D.data', shell=True)
		subprocess.call('rm '+ pathToData + '/bao2D.cov', shell=True)

	return
def replaceValueByMean():
	'''
		
		Replace the values of the correlation by the one of the mean

	'''
	
	rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_NicolasDistortionWithDistortion/xi_delta_QSO_result_'

	list1D = numpy.load(rawPath+'1D.npy')
	list2D = numpy.load(rawPath+'2D.npy')
	listMu = numpy.load(rawPath+'Mu.npy')
	listWe = numpy.load(rawPath+'We.npy')
	listMultipol = numpy.load(rawPath+'Multipol.npy')

	cov1D = numpy.load(rawPath+'cov_1D.npy')
	cov2D = numpy.load(rawPath+'cov_2D.npy')
	covMu = numpy.load(rawPath+'cov_Mu.npy')

	nbReal = list1D[0,:].size

	xi1D_[:,1] = numpy.mean(list1D,axis=1)
	xi1D_[:,2] = numpy.sqrt(numpy.diag(cov1D))/numpy.sqrt(nbReal)

	xi2D_[:,:,1] = myTools.convert1DTo2D( numpy.mean(list2D,axis=1), nbBinX2D__,nbBinY2D__)
	xi2D_[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(cov2D))/numpy.sqrt(nbReal), nbBinX2D__,nbBinY2D__)

	xiMu_[:,:,2] = myTools.convert1DTo2D( numpy.mean(listMu,axis=1), nbBin1D__,nbBinM__)
	xiMu_[:,:,3] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(covMu))/numpy.sqrt(nbReal), nbBin1D__,nbBinM__)

	result_Multipol_[:,:,1] = numpy.mean(listMultipol,axis=2)
	result_Multipol_[:,:,2] = numpy.var(listMultipol,axis=2)/numpy.sqrt(nbReal)

	xiWe_[:,:,1] = numpy.mean(listWe,axis=2)
	xiWe_[:,:,2] = numpy.var(listWe,axis=2)/numpy.sqrt(nbReal)


	return
def plotCovarDifferentMethod():
	'''

	'''

	dim = '2D'

	listPath = [    '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/subsampling_LYA_QSO_cov_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/shuffleQSO_LYA_QSO_cov_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/randomQSO_LYA_QSO_cov_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/shuffleForest_LYA_QSO_cov_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_'+dim+'_meanSubSampling.npy',
			]
	listPath2 = [   '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/subsampling_LYA_QSO_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/shuffleQSO_LYA_QSO_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/randomQSO_LYA_QSO_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/shuffleForest_LYA_QSO_'+dim+'.npy',
			'/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_'+dim+'.npy',
			]
	listName = ['Data \, subsampling',
			'Data \, shuffle \, QSO',
			'Data \, random \, QSO',
			'Data \, shuffle \, forest',
			'Mocks',
			'< Mock \, subsampling >',
			]

	real = [ numpy.load(i) for i in listPath2 ]
	cov  = [ numpy.load(i) for i in listPath ]

	### Plot the realisation
	for i in numpy.arange(len(real)):
		print listName[i]
		for j in numpy.arange(real[i][0,:].size):
			plt.errorbar(numpy.arange(real[i][:,j].size), real[i][:,j],fmt='o',color='blue',alpha=0.1)
		plt.errorbar(numpy.arange(real[i][:,j].size), numpy.mean(real[i],axis=1),fmt='o',color='red',label=r'$Mean$')
		plt.xlabel(r'$bin \, index$', fontsize=40)
		plt.ylabel(r'$\xi(|s|)$', fontsize=40)
		plt.title(r'$'+listName[i]+'$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.xlim([ -1., cov[i][0,:].size+1 ])
		plt.show()


	### Plot diagonal
	for i in numpy.arange(len(cov)):
		plt.errorbar(numpy.arange(cov[i][0,:].size), numpy.diag(cov[i]),fmt='o',label=r'$'+listName[i]+'$')
	plt.xlabel(r'$bin \, index$', fontsize=40)
	plt.ylabel(r'$Var(|s|)$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ -1., cov[i][0,:].size+1 ])
	plt.show()

	myTools.plotCovar(cov,listName)

	return
def saveMeanCov():
	'''

	'''


	tmpcov1D = numpy.zeros( shape=(5000,5000) )
	for i in range(0,10):
		for j in range(0,10):
			print i, j
			tmpcov1D += numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/subSampling_LYA_QSO_cov_2D.npy')

	tmpcov1D /= 100
	numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_meanSubSampling.npy',tmpcov1D)

'''
saveListRealMocks(10,10,True)
cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_NicolasDistortionWithDistortion/xi_delta_QSO_result_cov_1D.npy')
cor = myTools.getCorrelationMatrix(cov)
myTools.plot2D(cor)
a = myTools.plotCovar([cor],['a'])
cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_NicolasDistortionWithDistortion/xi_delta_QSO_result_cov_2D.npy')
cor = myTools.getCorrelationMatrix(cov)
myTools.plot2D(cor)
#a = myTools.plotCovar([cor],['a'])
'''

'''
### Data
#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_test_PDFMocksJMC_meanLambda_testNoCap/'
#saveListReal(5000,path+'xi_delta_QSO_Mu_LYA_QSO_shuffleForest_',path+'xi_delta_QSO_2D_LYA_QSO_shuffleForest_',path+'shuffleForest_LYA_QSO_',False)

cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_test_PDFMocksJMC_meanLambda_testNoCap/shuffleForest_LYA_QSO_cov_2D.npy')
cor = myTools.getCorrelationMatrix(cov)
myTools.plot2D(cor)
a = myTools.plotCovar([cor],['a'])
myTools.plot2D(a)
'''

#prepareForBAOFIT()
xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests6/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests6/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')
'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests5/xi_delta_QSO_distortionMatrix_1D_LYA_QSO.txt'
data = numpy.loadtxt(path)
xi1D_[:,1] = numpy.dot(data,xi1D_[:,1])
myTools.plot2D(data)
'''
'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests5/xi_delta_QSO_distortionMatrix_2D_LYA_QSO.txt'
data = numpy.loadtxt(path)
xi2D_[:,:,1] = myTools.convert1DTo2D(numpy.dot(data,myTools.convert2DTo1D(xi2D_[:,:,1], 50,100)),50,100)
myTools.plot2D(data)
'''

'''
xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results_NicolasDistortion/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results_NicolasDistortion/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')
result_Multipol_ = plotMultipol(xiMu_)
replaceValueByMean()
'''

#plotXi(0)
#plotXi(1)
plotXi(2)

'''
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,0)

plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
'''
plotWe(2)


#prepareForBAOFIT()





