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
import copy

### Perso lib
import correlation_3D
import myTools


raw_dic_Q = {
	'nb_random' : 10
	'estimator' : 'LS'
	'path_to_cat' : 'NOTHING'
	'load_from_txt' : True
}


class Correlation3DQ(correlation_3D.Correlation3D):

	def __init__(self, dic=None, dic2=None):

		if (dic2 is None): return

		cat = pyfits.open(dic2['path_to_cat'])[1].data
		self._nd = cat.size
		self._nr = cat.size
		del cat
		self._nbRand = dic2['nb_random']
		self._load_txt = dic2['load_from_txt']
		self._estimator = 'LS'

		correlation_3D.Correlation3D.__init__(self,dic)

		return
	def read_data(self, path1D, path2D, selection=0, init=False):

		if (self._load_txt):
			xiMu, xiWe, xi1D, xi2D = getACorrelation()
		else:
			xiMu = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_Mu.npy')
			xiWe = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_We.npy')
			xi1D = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_1D.npy')
			xi2D = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_2D.npy')

		return xiMu, xiWe, xi1D, xi2D
	def getACorrelation():

		### Which estimator
		LS = False
		if (self._estimator=='LS'): LS = True
		
		### Coef to normalize
		coefDD = self._nd*(self._nd-1.)/2.
		if (LS): coefDR = self._nd*self._nr
		coefRR = self._nr*(self._nr-1.)/2.
	
		### 2D: DD
		dd_2D = loadData2D(self._path_to_txt_file_folder+'xi_QSO_QSO_2D_QSO_DD.txt',coefDD)
		### Mu: DD
		dd_1D, dd_Mu, dd_We  = loadDataMu(self._path_to_txt_file_folder+'xi_QSO_QSO_Mu_QSO_DD.txt',coefDD)
	
		
		### 2D: RR,DR
		rr_2D = loadData2D(self._path_to_txt_file_folder+'xi_QSO_QSO_2D_QSO_RR_0.txt',coefRR)
		if (LS):
			dr_2D = loadData2D(self._path_to_txt_file_folder+'xi_QSO_QSO_2D_QSO_DR_0.txt',coefDR)
		### Mu: RR,DR
		rr_1D, rr_Mu, rr_We = loadDataMu(self._path_to_txt_file_folder+'xi_QSO_QSO_Mu_QSO_RR_0.txt',coefRR)
		if (LS):
			dr_1D, dr_Mu, dr_We = loadDataMu(self._path_to_txt_file_folder+'xi_QSO_QSO_Mu_QSO_DR_0.txt',coefDR)
		
		for i in range(1,self._nbRand):
		
			### 2D:
			rr_2D += loadData2D(self._path_to_txt_file_folder+'xi_QSO_QSO_2D_QSO_RR_'+str(i)+'.txt',coefRR)
		        if (LS):
				dr_2D +=  loadData2D(self._path_to_txt_file_folder+'xi_QSO_QSO_2D_QSO_DR_'+str(i)+'.txt',coefDR)
			### Mu:
			tmp_rr_1D, tmp_rr_Mu, tmp_rr_We = loadDataMu(self._path_to_txt_file_folder+'xi_QSO_QSO_Mu_QSO_RR_'+str(i)+'.txt',coefRR)
			rr_1D += tmp_rr_1D
			rr_Mu += tmp_rr_Mu
			rr_We += tmp_rr_We
			if (LS): 
				tmp_dr_1D, tmp_dr_Mu, tmp_dr_We = loadDataMu(self._path_to_txt_file_folder+'xi_QSO_QSO_Mu_QSO_DR_'+str(i)+'.txt',coefDR)
				dr_1D += tmp_dr_1D
				dr_Mu += tmp_dr_Mu
				dr_We += tmp_dr_We
			
		
	
		### 2D:
		rr_2D /= self._nbRand
		if (LS): dr_2D /= self._nbRand
		
		### Mu:
		rr_1D /= self._nbRand
		if (LS): dr_1D /= self._nbRand
		rr_Mu /= self._nbRand
		if (LS): dr_Mu /= self._nbRand
		rr_We /= self._nbRand
		if (LS): dr_We /= self._nbRand
	
		### Get the Landy-Saley estimator
		if (LS):		
			### Mu:
			result_Mu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,3))
			cut = (dd_Mu[:,:,1]>1.)
			result_Mu[:,:,0][cut] = dd_Mu[:,:,0][cut]
			result_Mu[:,:,1][cut] = (dd_Mu[:,:,1][cut] -2.*dr_Mu[:,:,1][cut] + rr_Mu[:,:,1][cut])/rr_Mu[:,:,1][cut]
			result_Mu[:,:,2][cut] = numpy.sqrt(dd_Mu[:,:,1][cut]*coefDD)/(coefDD*rr_Mu[:,:,1][cut])
		
			### xiWe
			result_We = numpy.zeros(shape=(self._nbBin1D,3,3))
			cut = (dd_We[:,:,1]>1.)
			result_We[:,:,0][cut] = dd_We[:,:,0][cut]
			result_We[:,:,1][cut] = (dd_We[:,:,1][cut] -2.*dr_We[:,:,1][cut] + rr_We[:,:,1][cut])/rr_We[:,:,1][cut]
			result_We[:,:,2][cut] = numpy.sqrt(dd_We[:,:,1][cut]*coefDD)/(coefDD*rr_We[:,:,1][cut])

			### 1D:
			result_1D = numpy.array( dd_1D )
			result_1D[:,1] = (dd_1D[:,1]-2.*dr_1D[:,1]+rr_1D[:,1])/rr_1D[:,1]
			result_1D[:,2] = numpy.sqrt(dd_1D[:,1]*coefDD)/(coefDD*rr_1D[:,1])
			cut = (dd_1D[:,1]==0.)
			result_1D[:,1][cut] = 0.
			result_1D[:,2][cut] = 0.

			### 2D:
			result_2D = numpy.array( dd_2D )
			result_2D[:,:,1] = (dd_2D[:,:,1]-2.*dr_2D[:,:,1]+rr_2D[:,:,1])/rr_2D[:,:,1]
			cut = (dd_2D[:,:,1]==0.)
			result_2D[:,:,1][cut] = 0.
	
		### Save
		#numpy.save(rawPath+'xi_QSO_QSO_result_1D',result_1D)
		#numpy.save(rawPath+'xi_QSO_QSO_result_2D',result_2D)
		#numpy.save(rawPath+'xi_QSO_QSO_result_Mu',result_Mu)
		#numpy.save(rawPath+'xi_QSO_QSO_result_We',result_We)
	
		#result_Multipol = plotMultipol(result_Mu)
		#numpy.save(rawPath+'xi_QSO_QSO_result_Multipol',result_Multipol)
	
		return result_Mu, result_We, result_1D, result_2D

corr = Correlation3DQ()

















