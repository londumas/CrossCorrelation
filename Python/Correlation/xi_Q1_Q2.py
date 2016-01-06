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

def plot_Xi_1D(rescale,log=False):

	path1 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation/'
	path = path1 +'xi_Q1_Q2_1D_'+qso1+'_'+qso2+'.txt'
	data = numpy.loadtxt(path)

	xxx = data[:,0]
	yyy = data[:,1]
	yer = numpy.sqrt(data[:,1])

	if (rescale==0):
		plt.errorbar(xxx, yyy,                 yerr=yer, fmt='o')
		plt.ylabel(r'$\xi_{qso \, qso} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xxx, yyy*xxx,       yerr=yer*xxx, fmt='o')
		plt.ylabel(r'$|s|.\xi_{qso \, qso} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xxx, yyy*xxx**2.,   yerr=yer*xxx**2., fmt='o')
		plt.ylabel(r'$|s|^{2}.\xi_{qso \, qso} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	if (rescale==3):
		plt.errorbar(xxx, yyy/xxx,       yerr=yer/xxx, fmt='o')
		plt.ylabel(r'$|s|^{-1}.\xi_{qso \, qso} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	if (rescale==4):
		plt.errorbar(xxx, yyy/(xxx**1.7), yerr=yer/(xxx**1.7), fmt='o')
		plt.ylabel(r'$|s|^{-1.7}.\xi_{qso \, qso} (|s|) \, [(h^{-1}.Mpc)^{1.7}]$', fontsize=40)
		
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	if (log):
		plt.xscale('log')
		plt.yscale('log')
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
	
	return
def plot_Xi_2D(rescale):

        ### If from cpp
        data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/xi_delta_delta_2D_'+forest__+'.txt')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xticks([ 0.,50.,100.,150.])
	ax.set_yticks([ 0.,50.,100.,150.])

        ### Xi2D
        xi2D = numpy.ndarray(shape=(const.nbBinX2D,const.nbBinX2D))
        xxx  = numpy.ndarray(shape=(const.nbBinX2D,const.nbBinX2D))
        for el in data:
                xi2D[int(el[0])/const.nbBinX2D,int(el[0])%const.nbBinX2D] = el[3]
                xxx[int(el[0])/const.nbBinX2D,int(el[0])%const.nbBinX2D] = numpy.sqrt(el[1]**2.+el[2]**2.)

        if (rescale==0):
                plt.imshow(xi2D, origin='lower',extent=[const.minX2D, const.maxX2D, const.minX2D, const.maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$\xi(\, \overrightarrow{s} \,)$',size=40)
        if (rescale==1):
                plt.imshow(xxx*xi2D, origin='lower',extent=[const.minX2D, const.maxX2D, const.minX2D, const.maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|.\xi(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
        if (rescale==2):
                plt.imshow(xxx**2.*xi2D, origin='lower',extent=[const.minX2D, const.maxX2D, const.minX2D, const.maxX2D],interpolation='None')
                cbar = plt.colorbar()
                cbar.set_label(r'$|s|^{2}.\xi(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

        plt.title(r'', fontsize=40)
        plt.ylabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
        plt.xlabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
        plt.grid(True)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
        plt.show()

        return

plot_Xi_1D(0)
#plot_Xi_1D(1)
#plot_Xi_1D(2)
plot_Xi_1D(3)
plot_Xi_1D(4)

plot_Xi_1D(0,True)
#plot_Xi_2D(0)
#plot_Xi_2D(1)
#plot_Xi_2D(2)
