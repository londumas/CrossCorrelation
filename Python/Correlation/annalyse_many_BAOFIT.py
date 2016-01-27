# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_3D.py
#

### Python lib
import numpy
import matplotlib.pyplot as plt

### Perso lib
import myTools
import correlation_3D
import annalyse_BAOFIT


class AnnalyseManyBAOFIT:

	def __init__(self, dic=None, index_parameter=None, path_to_simu='', nbBox=1, nbSimu=1):

		self._nbBox  = nbBox
		self._nbSimu = nbSimu

		if (dic is None):
			dic = correlation_3D.raw_dic_class

		### Get all the fits
		self._listFit = []
		for i in range(0,self._nbBox):
			for j in range(0,self._nbSimu):
				path = path_to_simu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/'
				dic['path_to_txt_file_folder'] = path+'Results/'
				dic['name'] = str(i)+'-'+str(j)
				self._listFit += [annalyse_BAOFIT.AnnalyseBAOFIT(dic, index_parameter, path+'BaoFit_q_f_covFromFit/bao2D.')]

		### Set attributes set after
		self._chi2_scan = None


		return
	def print_results(self):

		for el in self._listFit:
			el.print_results()

		return
	def plot_histo_residuals(self, nbBins=100):

		yyy = numpy.array([])
		for el in self._listFit:
			yyy = numpy.append(yyy,el.get_residuals())

		### histo
		fig = plt.figure()
		ax = fig.add_subplot(111)
	
		ax.hist(yyy, bins=nbBins)
		plt.xlabel(r'$\frac{data-fit}{\sigma_{data}}$')
		plt.ylabel(r'$\#$')
		myTools.deal_with_plot(False,False,False)
		
		mng = plt.get_current_fig_manager()
		textstr = '$nb=%u$\n$\mu=%.5e$\n$\sigma=%.5e$'%(yyy.size, numpy.mean(yyy), numpy.std(yyy))
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
		
		plt.show()
	def plot_chi2_scan(self, sizeX=100, sizeY=100, edge=None, toyMC=False, bootstrap=False):

		self._chi2_scan = numpy.zeros( shape=(sizeX,sizeY) )
		nb = numpy.zeros( shape=(sizeX,sizeY) )

		for el in self._listFit:
			el.read_chi2_scan(sizeX, sizeY)
			cut = numpy.isfinite(el._chi2_scan)
			self._chi2_scan[cut] += el._chi2_scan[cut]
			nb[cut] += 1.

		self._chi2_scan /= nb

		### Plot the data
		myTools.plot_chi2_scan(self._chi2_scan,self._listFit[0]._chi2_scan_edge,True, None,'\\alpha_{\\perp}','\\alpha_{\\parallel}','\\Delta \\chi^{2} = \\chi^{2}-\\chi^{2}_{best \\, fit}')

		return






































