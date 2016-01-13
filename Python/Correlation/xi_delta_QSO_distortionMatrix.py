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
nbBinM__ = 100;




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


	'''
	### --------------------------------------------
	### If needed to save the grid
	rParal = tmp_save3 / tmp_save5
	rPerp  = tmp_save2 / tmp_save5
	z      = tmp_save4 / tmp_save5

	tmp_rParal = numpy.zeros(nbBin2D__)
	tmp_rPerp  = numpy.zeros(nbBin2D__)
	tmp_z      = numpy.zeros(nbBin2D__)
	tmp_idx    = numpy.arange(0,nbBin2D__)
	for k1 in range(0,nbBin2D__):
		i1       = k1/nbBinY2D__
		j1       = k1%nbBinY2D__
		k11      = j1*nbBinX2D__ + i1
		tmp_rParal[k11] = rParal[i1,j1]
		tmp_rPerp[k11]  = rPerp[i1,j1]
		tmp_z[k11]      = z[i1,j1]
	
	numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/BaoFit_q_f/bao2D.grid',zip(tmp_idx,tmp_rParal,tmp_rPerp,tmp_z),fmt='%u %1.20e %1.20e %1.20e')
	### --------------------------------------------
	'''

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
def fitCamb(data,pathToFile,mulpol=0):
	'''

	'''

	### Constants
	startFit   = 25.
	endFit     = 60.
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

#nbPixel = nbBin1D__
nbPixel = nbBin2D__

#xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_DR12_nicolas/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_DR12_nicolas/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')
#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_test_PDFMocksJMC_meanLambda_testNoCap/xi_delta_QSO_distortionMatrix_1D_LYA_QSO.txt'


xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_003/Results_NicolasDistortion/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_003/Results_NicolasDistortion/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_003/Results_NicolasDistortion/xi_delta_QSO_distortionMatrix_2D_LYA_QSO.txt'


print path
data = numpy.loadtxt(path)
myTools.plot2D(data)
plt.hist(data[data!=0.],bins=100)
plt.hist(numpy.diag(data),bins=100)
plt.show()

'''
cov = numpy.array(data)
### Test if symetric
for i in numpy.arange(nbPixel):
	for j in numpy.arange(nbPixel):
		data[i][j] -= cov[j][i]
data[ data==0. ] = numpy.float('nan')
myTools.plot2D(data)
'''

plt.hist(data[data!=0.],bins=100)
plt.hist(numpy.diag(data),bins=100)
plt.show()



if (nbPixel==nbBin1D__):
	xi1D = xi1D_[:,1]
else:
	xi1D = myTools.convert2DTo1D(xi2D_[:,:,1], 50,100)

xi1D_2 = numpy.dot(data,xi1D)

plt.errorbar(numpy.arange(nbPixel),xi1D,fmt='o',label='Before correction')
plt.errorbar(numpy.arange(nbPixel),xi1D_2,fmt='o',label='After correction')
myTools.deal_with_plot(False,False,True)
plt.show()

if (nbPixel==nbBin1D__):
	xi1D_[:,1] = xi1D_2
	plotXi(0)
	plotXi(1)
	plotXi(2)
	pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
	fitCamb(xi1D_,pathToCamb,0)
else:
	xi2D_[:,:,1] = myTools.convert1DTo2D(xi1D_2,50,100)
	plotXi2D(0)
	plotXi2D(1)
	plotXi2D(2)






#path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_test_PDFMocksJMC_meanLambda_testNoCap/'
#path = path__ + 'xi_delta_QSO_distortionMatrix_1D_LYA_QSO.txt'
#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_test_PDFMocksJMC_meanLambda_testNoCap//xi_delta_QSO_distortionMatrix_1D_LYA_QSO.txt'
#print path
#data = numpy.loadtxt(path)
#numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests_test_PDFMocksJMC_meanLambda_testNoCap//xi_delta_QSO_distortionMatrix_1D_LYA_QSO',data)








print data.size
print data[0].size
print data

data[ data==0. ] = numpy.float('nan')
myTools.plot2D(data)

print '  is diag <0.   : ', numpy.diag(data)[ numpy.diag(data)<0. ].size
print '  is diag ==0.  : ', numpy.diag(data)[ numpy.diag(data)==0. ].size


cor = myTools.getCorrelationMatrix(data)
print cor
print cor[ cor==1. ].size
print cor[ cor>1. ].size
print cor[ cor==0. ].size

#cor[ cor==0. ] = numpy.float('nan')
#cor[ cor==1. ] = numpy.float('nan')
#myTools.plot2D(cor)

plt.hist( numpy.diag(data), bins=100  )
plt.show()

plt.hist( cor[ numpy.logical_and( cor!=0., cor!=1.) ], bins=100  )
plt.show()

cor2 = numpy.array(cor)
### Test if symetric
size = cor[:,0].size
for i in range(0,size):
	for j in range(0,size):
		cor2[i,j] -= cor[i,j]

cor2[ cor2==0. ] = numpy.float('nan')
cor2[ cor2==1. ] = numpy.float('nan')
myTools.plot2D(cor2)
		

plt.hist( cor2[ numpy.logical_and( cor!=0., cor!=1.) ], bins=100  )
plt.show()
















