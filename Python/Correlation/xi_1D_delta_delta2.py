# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import sys
import numpy
import matplotlib.pyplot as plt
import myTools
from const_delta import *

forest1 = sys.argv[1]
forest2 = sys.argv[2]
qso1    = sys.argv[3]
qso2    = sys.argv[4]

def plot():

	path1 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/'
	path  = path1 +'xi_1D_delta_delta2_'+forest1+'_'+forest2+'.txt'	
	data = numpy.loadtxt(path)
	print path

	xxx  = data[:,0]
	yyy  = data[:,1]
	yer  = data[:,2]

	plt.errorbar(xxx, yyy, yerr=yer, fmt='o')
	plt.title(r'$1D: \, \delta_{'+forest1+'} \, - \, \delta_{'+forest2+'} $', fontsize=40)
	plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
	plt.ylabel(r'$\xi (s)$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(yyy)+10. ])
	plt.show()
	
	return

plot()
