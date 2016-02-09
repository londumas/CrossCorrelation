# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_3D.py
#

### Python lib
import subprocess
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit

### Perso lib
import const
import myTools

class CAMB:
	
	def __init__(self):

		self._xi0 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.0.dat')
		self._xi2 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.2.dat')
		self._xi4 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.4.dat')

		return
	def plot_1d(self,x_power=0):

		xxx = self._xi0[:,0]
		yyy = self._xi0[:,1]
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,fmt='o')

		xxx = self._xi2[:,0]
		yyy = self._xi2[:,1]
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,fmt='o')

		xxx = self._xi4[:,0]
		yyy = self._xi4[:,1]
		coef = numpy.power(xxx,x_power)

		plt.errorbar(xxx,coef*yyy,fmt='o')
		if (x_power==0):
			plt.ylabel(r'$ \xi (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|.\xi (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,True)
		plt.show()

		return
	def save_in_one_file(self, path_to_save):

		numpy.savetxt(path_to_save,zip(self._xi0[:,0], self._xi0[:,1], self._xi2[:,1], self._xi4[:,1]),fmt='%1.20e %1.20e %1.20e %1.20e')
		

		return

camb = CAMB()
camb.plot_1d(0)
camb.plot_1d(1)
camb.plot_1d(2)
camb.save_in_one_file('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/CAMB/CAMB_2_25/camb.txt')












