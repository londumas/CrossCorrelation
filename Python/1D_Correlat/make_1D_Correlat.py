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
import myTools
from myTools import Get_TProfile
from const_delta import *
## Cosmology
import cosmolopy.distance as cosmology


def get_xi1D():
	path = '../test2.fits'
	print path
	file_cat = pyfits.open(path, memmap=True)
	cat = file_cat[1].data[:10000]
	sizeMax = cat.size
	print sizeMax
	
	step   = 100
	iStart = 0
	iEnd   = step
	
	distBinEdges = numpy.arange(0.,1000.+1.,1.)
	weight, axisX, axisY = numpy.histogram2d([-1.], [0.], (distBinEdges,1))
	
	keep_delta  = weight[:,0]
	keep_weight = weight[:,0]
	keep_var    = weight[:,0]
	keep_number = weight[:,0]

	while (iStart<=sizeMax):
		
		tmp_dist   = []
		tmp_delta  = []
		tmp_weight = []
		tmp_var    = []
		tmp_number = []
	
		for el in cat[iStart:iEnd]:
		
			ar_cut    = (el['DELTA_WEIGHT']>0.)
			ar_delta  = el['DELTA'][ ar_cut ]
			ar_weight = el['DELTA_WEIGHT'][ ar_cut ]
			ar_dist   = el['R_SPHERICAL'][ ar_cut ]
			
			sizeMax = ar_delta.size
			
			delta1 = numpy.ones((sizeMax,1))
			delta2 = numpy.ones((1,sizeMax))
			
			weight1 = delta1
			weight2 = delta2
			weight1[:,0] = ar_weight
			weight2[0]   = ar_weight
			weight3 = numpy.dot(weight1,weight2)
			
			delta1[:,0] = ar_delta
			delta2[0]   = ar_delta
			delta3 = numpy.dot(delta1,delta2)
			
			dist      = numpy.ones((sizeMax,1))
			dist[:,0] = ar_dist
			dist1     = numpy.dot( dist, numpy.ones((1,sizeMax)))
			dist2     = numpy.dot( dist, numpy.ones((1,sizeMax)))
			dist      = numpy.transpose(dist2)-dist1		
			
			dist    = dist.flat
			delta3  = delta3.flat
                        weight3 = weight3.flat

			cut_moreZero = (dist>0.)
			
			tmp_dist   += list(dist[cut_moreZero])
			tmp_delta  += list(delta3[cut_moreZero])
			tmp_weight += list(weight3[cut_moreZero])

		number, axisX, axisY = numpy.histogram2d(tmp_dist, tmp_delta, (distBinEdges,1))
		mean,   axisX, axisY = numpy.histogram2d(tmp_dist, tmp_delta, (distBinEdges,1), weights=numpy.array(tmp_weight)*tmp_delta)
		weight, axisX, axisY = numpy.histogram2d(tmp_dist, tmp_delta, (distBinEdges,1), weights=tmp_weight)
		var,    axisX, axisY = numpy.histogram2d(tmp_dist, tmp_delta, (distBinEdges,1), weights=numpy.array(tmp_weight)*tmp_delta*tmp_delta)
		keep_delta  = keep_delta  + mean[:,0]
		keep_weight = keep_weight + weight[:,0]
		keep_var    = keep_var    + var[:,0]
		keep_number = keep_number + number[:,0]
		
		iStart  = iEnd+1
		iEnd   += step
		
	print '  Finished the loop'

	cut = (keep_number>0.)
	keep_delta  = keep_delta[cut]
	keep_weight = keep_weight[cut]
	keep_var    = keep_var[cut]
	keep_number = keep_number[cut]
	
	distBinCenter = numpy.arange(0.5,1000.,1.)
	keep_delta /= keep_weight
	keep_var = numpy.sqrt((keep_var/keep_weight-keep_delta*keep_delta)/keep_number)

	### Get for the xi1D(0)
	ar_cut    = (cat['DELTA_WEIGHT']>0.)
	ar_delta  = cat['DELTA'][ ar_cut ]
	ar_weight = cat['DELTA_WEIGHT'][ ar_cut ]
	ar_delta  *= ar_delta
	ar_weight *= ar_weight
	meanDelta = [ numpy.average(ar_delta,weights=ar_weight ) ]
	stdDelta  = [ numpy.average( (ar_delta-meanDelta[0])**2., weights=ar_weight )/numpy.sqrt(ar_delta.size) ]

	distBinCenter = numpy.append( [0.], distBinCenter )
	keep_delta    = numpy.append( meanDelta, keep_delta )
	keep_var      = numpy.append( stdDelta, keep_var)

	numpy.save('data_xi1D', zip( distBinCenter, keep_delta, keep_var))

	return

def plotXi1D():

	### If from Python
	#data = numpy.load('data_xi1D.npy')

	### If from cpp
	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/xi_delta_delta_'+forest__+'.txt')
	
	plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt='o')
	plt.title(r'$1D: \, \delta_{'+forest__+'} \, - \, \delta_{'+forest__+'} $', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	plt.ylabel(r'$\xi (|s|)$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(data[:,0])-10., numpy.max(data[:,0])+10. ])
	plt.show()
	
	return

	
#get_xi1D()
plotXi1D()
