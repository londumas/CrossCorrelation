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
import copy

### Perso lib
import myTools
import correlation_3D
import annalyse_BAOFIT


class AnnalyseManyBAOFIT:

	def __init__(self, dic=None, index_parameter=None, path_to_simu='', nbBox=1, nbSimu=1):

		self._nbBox  = nbBox
		self._nbSimu = nbSimu
		self._nbtot  = self._nbBox*self._nbSimu
		self._path_to_simu = path_to_simu

		if (dic is None):
			dic = correlation_3D.raw_dic_class

		### Get all the fits
		self._listFit = []
		for i in range(0,self._nbBox):
			for j in range(0,self._nbSimu):
				path = self._path_to_simu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/'
				dic['path_to_txt_file_folder'] = path+'Results/'
				dic['name'] = str(i)+'-'+str(j)
				self._listFit += [annalyse_BAOFIT.AnnalyseBAOFIT(dic, index_parameter, path+'BaoFit_q_f_covFromFit/bao2D.')]

		### Init parameters
		self._nbBin1D  = self._listFit[0]._nbBin1D
		self._nbBin2D  = self._listFit[0]._nbBin2D
		self._nbBinX2D = self._listFit[0]._nbBinX2D
		self._nbBinY2D = self._listFit[0]._nbBinY2D
		self._nbBinM   = self._listFit[0]._nbBinM
		### Set attributes set after
		self._chi2_scan = None

		### Create a 'AnnalyseBAOFIT' object for mean_fit
		#dic['path_to_txt_file_folder'] = self._path_to_simu + 'Box_000/Simu_000/Results/'
		#dic['name'] = '<fit \, simu>'
		#self._mean_fit = annalyse_BAOFIT.AnnalyseBAOFIT(dic, index_parameter, path+'BaoFit_q_f_covFromFit/bao2D.')

		### Create a 'Correlation3D' object for mean_data
		dic['path_to_txt_file_folder'] = self._path_to_simu + 'Box_000/Simu_000/Results/'
		dic['name'] = '<simulation>'
		self._mean_data = correlation_3D.Correlation3D(dic)
		self._mean_data._xiMul = self._mean_data.get_multipol(self._mean_data._xiMu)

		return
	def save_list_realisation(self):

		list1D       = numpy.zeros( shape=(self._nbBin1D,self._nbtot) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,self._nbtot) )
		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,self._nbtot) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,3,self._nbtot) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,5,self._nbtot) )

		pathToSave = self._path_to_simu + 'Results_16_01_28/'

		for i in numpy.arange(self._nbtot):
			print i

			list1D[:,i]         = self._listFit[i]._xi1D[:,1]
			list2D[:,i]         = self._listFit[i]._xi2D[:,:,1].flatten()
			listMu[:,i]         = self._listFit[i]._xiMu[:,:,2].flatten()
			listWe[:,:,i]       = self._listFit[i]._xiWe[:,:,1]
			listMultipol[:,:,i] = self._listFit[i].get_multipol( self._listFit[i]._xiMu)[:,:,1]

		numpy.save(pathToSave+'list_Mu',listMu)
		numpy.save(pathToSave+'list_We',listWe)
		numpy.save(pathToSave+'list_1D',list1D)
		numpy.save(pathToSave+'list_2D',list2D)
		numpy.save(pathToSave+'list_Multipol',listMultipol)

		covMu = numpy.cov(listMu)
		cov1D = numpy.cov(list1D)
		cov2D = numpy.cov(list2D)

		numpy.save(pathToSave+'cov_Mu',covMu)
		numpy.save(pathToSave+'cov_1D',cov1D)
		numpy.save(pathToSave+'cov_2D',cov2D)
	
		return
	def set_mean_fit(self):

		list1D       = numpy.zeros( shape=(self._nbBin1D,self._nbtot) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,self._nbtot) )
		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,self._nbtot) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,3,self._nbtot) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,5,self._nbtot) )

		return
	def set_mean_data(self):

		raw_path = self._path_to_simu + 'Results_16_01_28/'
		### we
		path = raw_path + 'list_We.npy'
		listWe = numpy.load(path)
		nbTot = listWe[0,0,:].size
		for i in numpy.arange(3):
			self._mean_data._xiWe[:,i,1] = numpy.mean(listWe[:,i,:],axis=1)
			self._mean_data._xiWe[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listWe[:,i,:] ))/nbTot)
		### mu
		path = raw_path + 'cov_Mu.npy'
		self._mean_data._xiMu[:,:,2] = myTools.convert1DTo2D(numpy.mean(numpy.load(raw_path + 'list_Mu.npy'),axis=1),self._nbBin1D,self._nbBinM)
		self._mean_data._xiMu[:,:,3] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path))/nbTot ),self._nbBin1D,self._nbBinM)
		### 1D
		path = raw_path + 'cov_1D.npy'
		self._mean_data._xi1D[:,1]   = numpy.mean(numpy.load(raw_path + 'list_1D.npy'),axis=1)
		self._mean_data._xi1D[:,2]   = numpy.sqrt( numpy.diag(numpy.load(path))/nbTot )
		### 2D
		path = raw_path + 'cov_2D.npy'
		self._mean_data._xi2D[:,:,1] = myTools.convert1DTo2D(numpy.mean(numpy.load(raw_path + 'list_2D.npy'),axis=1),self._nbBinX2D, self._nbBinY2D)
		self._mean_data._xi2D[:,:,2] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path))/nbTot ),self._nbBinX2D, self._nbBinY2D)
		### Multipol
		path = raw_path + 'list_Multipol.npy'
		listMultipol = numpy.load(path)
		for i in numpy.arange(5):
			self._mean_data._xiMul[:,i,1] = numpy.mean(listMultipol[:,i,:],axis=1)
			self._mean_data._xiMul[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listMultipol[:,i,:] ))/nbTot)

		return
	def print_results(self):

		for el in self._listFit:
			el.print_results()

		return
	def plot_1d(self, x_power=0, other=[], title=True):
		self._mean_data.plot_1d(x_power, other, title)
		
		return
	def plot_we(self, x_power=0, other=[], title=True):
		self._mean_data.plot_we(x_power, other, title)
		
		return
	def plot_2d(self, x_power=0):
		self._mean_data.plot_2d(x_power)
		
		return
	def plot_mu(self, x_power=0):
		self._mean_data.plot_mu(x_power)
	
		return
	def plot_multipol(self, x_power=0):
		self._mean_data.plot_multipol(x_power)
	
		return
	def plot_histo_residuals(self, nbBins=100):

		yyy = numpy.array([])
		for el in self._listFit:
			xi2D = el.get_residuals()
			tmp_yyy  = (xi2D[:,:,1][ (xi2D[:,:,2]>0.) ]).flatten()
			yyy = numpy.append(yyy,tmp_yyy)

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
	def plot_chi2_scan(self, sizeX=100, sizeY=100, edge=None, bestFitData=None):

		self._chi2_scan = numpy.zeros( shape=(sizeX,sizeY) )
		nb = numpy.zeros( shape=(sizeX,sizeY) )

		if (bestFitData is None): bestFit = numpy.zeros( shape=(2,1) )
		else: bestFit = bestFitData

		for el in self._listFit:
			bestFit = numpy.append(bestFit,el.get_best_fit(),axis=1)
			el.read_chi2_scan(sizeX, sizeY)
			cut = numpy.isfinite(el._chi2_scan)
			self._chi2_scan[cut] += el._chi2_scan[cut]
			nb[cut] += 1.

		self._chi2_scan /= nb

		### Plot the data
		myTools.plot_chi2_scan(self._chi2_scan,self._listFit[0]._chi2_scan_edge,True, bestFit,'\\alpha_{\\perp}','\\alpha_{\\parallel}','\\Delta \\chi^{2} = \\chi^{2}-\\chi^{2}_{best \\, fit}')

		return






































