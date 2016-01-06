# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >



import subprocess
import sys
import os

import numpy
import astropy.io.fits as pyfits
from iminuit import Minuit
import matplotlib.pyplot as plt
import cosmolopy.distance as cosmology  ### Cosmology



import myTools
from myTools import Get_TProfile
from const_delta import *


nbRegion = 80

### 1D
min1D   = 0.
max1D   = 160.
nbBin1D = 16

### 2D
minX2D   = min1D
maxX2D   = max1D
nbBinX2D = nbBin1D

minY2D   = -max1D
maxY2D   = max1D
nbBinY2D = 2*nbBin1D

nbBin2D  = nbBinX2D*nbBinY2D


def loadData(binSize=1):
	'''
	
		Resize the data from a bining of 1 Mpc.h-1 to a binSize
	
	'''

	data1D = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/xi_delta_QSO_1D_'+forest__+'.txt')
	#data1D = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/xi_delta_QSO1D.txt_backup)
	save0 = data1D[:,3]
	save1 = data1D[:,4]
	save2 = data1D[:,5]
	save3 = data1D[:,6]
	save4 = data1D[:,7]
	save5 = data1D[:,8]

	tmp_save0 = numpy.zeros(nbBin1D)
	tmp_save1 = numpy.zeros(nbBin1D)
	tmp_save2 = numpy.zeros(nbBin1D)
	tmp_save3 = numpy.zeros(nbBin1D)
	tmp_save4 = numpy.zeros(nbBin1D)

	for i in range(0,len(save0)):
		idx = save5[i]/binSize
		tmp_save0[idx] += save0[i]
		tmp_save1[idx] += save1[i]
		tmp_save2[idx] += save2[i]
		tmp_save3[idx] += save3[i]
		tmp_save4[idx] += save4[i]
		

	xi1D = numpy.zeros(shape=(nbBin1D,3) )

	for i in range(0,nbBin1D):

		if (tmp_save4[i]==0.):
			continue

		xxx = tmp_save2[i]/tmp_save3[i]
		yyy = tmp_save0[i]/tmp_save3[i]
		yer = numpy.sqrt( (tmp_save1[i]/tmp_save3[i] - yyy*yyy)/tmp_save4[i] )

		xi1D[i][0] = xxx
		xi1D[i][1] = yyy
		xi1D[i][2] = yer

	return xi1D


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
	
	plt.title(r'$\delta_{'+forest__+'} \, - \, qso$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
	
	return
def plotXi2D(rescale):

	### If from cpp
	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/xi_delta_QSO_2D_'+forest__+'.txt')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_yticks([ 0.,50.,100.,150.])

	### Xi2D
	xi2D = numpy.ndarray(shape=(nbBinX2D,nbBinY2D))
	xxx  = numpy.ndarray(shape=(nbBinX2D,nbBinY2D))
	for el in data:
		xi2D[int(el[0])/nbBinY2D,int(el[0])%nbBinY2D] = el[3]
		xxx[int(el[0])/nbBinY2D,int(el[0])%nbBinY2D] = numpy.sqrt(el[1]**2.+el[2]**2.)

	if (rescale==0):
		plt.imshow(xi2D, origin='lower',extent=[minY2D, maxY2D, minX2D, maxX2D],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		plt.imshow(xxx*xi2D, origin='lower',extent=[minY2D, maxY2D, minX2D, maxX2D],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|.\xi(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
	if (rescale==2):
		plt.imshow(xxx**2.*xi2D, origin='lower',extent=[minY2D, maxY2D, minX2D, maxX2D],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|^{2}.\xi(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

	plt.title(r'$\delta_{'+forest__+'} \, - \, qso$', fontsize=40)
	plt.ylabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
	plt.xlabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()
	
	return






def get_xiX_same_Forest():
	path = '../test2.fits'
	print path
	file_cat = pyfits.open(path, memmap=True)
	cat = file_cat[1].data
	sizeMax = cat.size
	print sizeMax
	
	### Cosmology
	omegaM0__      = 0.27
	omegaLambda0__ = 0.73
	h__            = 0.71
	cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}
	distfunc, redfunc = cosmology.quick_distance_function(cosmology.comoving_distance, zmax=20.0, zmin=0.0, zstep=0.001, return_inverse=True, **cosmo)
	
	distBinEdges = numpy.arange(0.,1000.+10.,10.)

	tmp_dist   = []
	tmp_delta  = []
	tmp_weight = []
	tmp_var    = []
	tmp_number = []
	
	#ar_distQSO = numpy.dot(numpy.diag(cat['Z_VI']), numpy.ones((sizeMax,nbBinRFMax__)))
	#ar_distQSO = (numpy.diag(cat['Z_VI'])*numpy.ones((sizeMax,nbBinRFMax__)).T).T	
	step = 500
        iStart = 0
        iEnd   = step
	ar_distQSO = numpy.ones((sizeMax,nbBinRFMax__))
        while (iStart<=sizeMax):
                print iStart, iEnd
                ar_distQSO[iStart:iEnd] = numpy.dot(numpy.diag(cat['Z_VI'][iStart:iEnd]), ar_distQSO[iStart:iEnd])
                iStart = iEnd+1
                iEnd   += step
                if (iEnd>sizeMax):
                        iEnd = sizeMax


	ar_cut     = (cat['DELTA_WEIGHT']>0.)
	ar_delta   = cat['DELTA'][ ar_cut ]
	ar_weight  = cat['DELTA_WEIGHT'][ ar_cut ]
	ar_distQSO = distfunc(ar_distQSO[ ar_cut ])*h__
	ar_dist    = ar_distQSO-cat['R_SPHERICAL'][ ar_cut ]
	print ar_dist.size
	print ar_dist.nbytes
	print numpy.average(ar_delta,weights=ar_weight)
	
	number, axisX, axisY = numpy.histogram2d(ar_dist, ar_delta, (distBinEdges,1))
	mean,   axisX, axisY = numpy.histogram2d(ar_dist, ar_delta, (distBinEdges,1), weights=ar_weight*ar_delta)
	weight, axisX, axisY = numpy.histogram2d(ar_dist, ar_delta, (distBinEdges,1), weights=ar_weight)
	var,    axisX, axisY = numpy.histogram2d(ar_dist, ar_delta, (distBinEdges,1), weights=ar_weight*ar_delta*ar_delta)
	keep_delta  = mean[:,0]
	keep_weight = weight[:,0]
	keep_var    = var[:,0]
	keep_number = number[:,0]
		
	print '  Finished the loop'
	
	
	distBinCenter = numpy.arange(5.,1000.,10.)
	keep_delta /= keep_weight
	keep_var = numpy.sqrt((keep_var/keep_weight-keep_delta*keep_delta)/keep_number)
	
	plt.errorbar(distBinCenter, keep_delta, yerr=keep_var)
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	plt.ylabel(r'$\xi_{\delta,q \, 1D} (|s|)$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.show()
	
	return


def xi_delta_QSO():
	'''
	'''

	cat_QSO   = pyfits.open('../QSO/test_QSO.fits', memmap=True)[1].data[:10000]
	cat_delta = pyfits.open('../test2.fits', memmap=True)[1].data[:1000]

	### Binning of the cross-correlation
	distBinEdges = numpy.arange(0.,1000.+10.,10.)

	### keep only needed pixels
	ar_cut   = (cat_delta['DELTA_WEIGHT']>0.)
	ar_delta = cat_delta['DELTA'][ ar_cut ]	
	sizeDelta = ar_delta.size

	ones_delta = numpy.ones((sizeDelta,1))
	delta   = numpy.ones((sizeDelta,1))
	weight  = numpy.ones((sizeDelta,1))
	delta_x = numpy.ones((sizeDelta,1))
	delta_y = numpy.ones((sizeDelta,1))
	delta_z = numpy.ones((sizeDelta,1))

	delta[:,0]   = ar_delta
	weight[:,0]  = cat_delta['DELTA_WEIGHT'][ ar_cut ]
	delta_x[:,0] = cat_delta['X_CARTESIAN'][ ar_cut ]
	delta_y[:,0] = cat_delta['Y_CARTESIAN'][ ar_cut ]
	delta_z[:,0] = cat_delta['Z_CARTESIAN'][ ar_cut ]

	del ar_delta, cat_delta, ar_cut
	print '  Finished preparing delta'

	### Deal with QSO
	sizeQSO = cat_QSO.size

	qso    = numpy.ones((1,sizeQSO))
        qso_x  = numpy.ones((1,sizeQSO))
        qso_y  = numpy.ones((1,sizeQSO))
        qso_z  = numpy.ones((1,sizeQSO))

	qso_x[0,:] = cat_QSO['X_CARTESIAN']
	qso_y[0,:] = cat_QSO['Y_CARTESIAN']
	qso_z[0,:] = cat_QSO['Z_CARTESIAN']

	print '  Finish preparing qso'

	del cat_QSO

	delta   = numpy.dot(delta,qso).ravel()
	weight  = numpy.dot(weight,qso).ravel()
	delta_x = (numpy.dot(delta_x,qso)-numpy.dot(ones_delta,qso_x))**2.
	delta_x += (numpy.dot(delta_y,qso)-numpy.dot(ones_delta,qso_y))**2.
	delta_x += (numpy.dot(delta_z,qso)-numpy.dot(ones_delta,qso_z))**2.
	delta_x = (delta_x**0.5).ravel()

	print '  Finished multiplication'

	number, axisX, axisY = numpy.histogram2d(delta_x, delta, (distBinEdges,1))
        mean,   axisX, axisY = numpy.histogram2d(delta_x, delta, (distBinEdges,1), weights=weight*delta)
        weight, axisX, axisY = numpy.histogram2d(delta_x, delta, (distBinEdges,1), weights=delta)
        #var,    axisX, axisY = numpy.histogram2d(delta_x, delta, (distBinEdges,1), weights=weight*delta**2.)

	print '  Finished'
	
	

	axisX = numpy.arange(5.,995.+10.,10.)
	cut   = (number[:,0]>0.)
	mean  = mean[:,0][cut]/weight[:,0][cut]
	axisX = axisX[cut]
	#var   = ((var[:,0][cut]/weight[:,0][cut]-mean**2.)/number[:,0][cut])**0.5


	print mean
	print mean.size
	print axisX
	print axisX.size
 	

	numpy.save('data_xi_delta_QSO', zip( axisX, mean))
	#return

        plt.errorbar(axisX, mean)
        plt.title(r'', fontsize=40)
        plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
        plt.ylabel(r'$\xi_{\delta,q \, 1D} (|s|)$', fontsize=40)
        myTools.deal_with_plot(False,False,False)
        plt.show()

	cut   = (axisX<200.)
        mean  = mean[cut]
        axisX = axisX[cut]

	plt.errorbar(axisX, mean)
        plt.title(r'', fontsize=40)
        plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
        plt.ylabel(r'$\xi_{\delta,q \, 1D} (|s|)$', fontsize=40)
        myTools.deal_with_plot(False,False,False)
	plt.show()

	return

xi1D_ = loadData(max1D/nbBin1D)
plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
#get_xiX_same_Foar_cut     = (cat['DELTA_WEIGHT']>0.)
#xi_delta_QSO()
