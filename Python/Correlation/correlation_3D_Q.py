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
import astropy.io.fits as pyfits

### Perso lib
import correlation_3D
import myTools


raw_dic_Q = {
	'nb_random' : 10,
	'estimator' : 'LS',
	'path_to_cat' : 'NOTHING',
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
	def read_data_QSO(self,path1D, path2D, coef, init=False):
		
		### Set arrays
		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
		xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4))
		xiWe = numpy.zeros(shape=(self._nbBin1D,3,3))
		xi1D = numpy.zeros(shape=(self._nbBin1D,3))

		int_binSize = int(self._binSize)
		binSizeY = self._nbBinM_calcul/self._nbBinM

		### Mu
		data = numpy.loadtxt(path1D)	
		save0 = data[:,0]
	
		for i in range(0,save0.size):

			iX = i/self._nbBinM_calcul
			iY = i%self._nbBinM_calcul
	
			### for mu
			idX = iX/int_binSize
			idY = iY/binSizeY
	
			xiMu[idX,idY,0] += data[i,1]
			xiMu[idX,idY,1] += data[i,2]
			xiMu[idX,idY,2] += data[i,0]
	
			### for we
			if (iY>40):
				idY = 0
			elif (iY>25 and iY<=40):
				idY = 1
			else:
				idY = 2
	
			xiWe[idX,idY,0] += data[i,1]
			xiWe[idX,idY,1] += data[i,0]
	
			### for 1D
			xi1D[idX,0] += data[i,1]
			xi1D[idX,1] += data[i,0]
	
		### Mu
		cut = (xiMu[:,:,1]>0.)
		xiMu[:,:,0][cut] /= xiMu[:,:,2][cut]
		xiMu[:,:,1][cut] /= xiMu[:,:,2][cut]
		xiMu[:,:,2][cut] /= coef
		xiMu[:,:,3][cut]  = numpy.sqrt(xiMu[:,:,2][cut])
		### we
		cut = (xiWe[:,:,1]>0.)
		xiWe[:,:,0][cut] /= xiWe[:,:,1][cut]
		xiWe[:,:,1][cut] /= coef
		xiWe[:,:,2][cut]  = numpy.sqrt(xiWe[:,:,1][cut])
		### 1D
		cut = (xi1D[:,1]>0.)
		xi1D[:,0][cut] /= xi1D[:,1][cut]
		xi1D[:,1][cut] /= coef
		xi1D[:,2][cut]  = numpy.sqrt(xi1D[:,1][cut])
		
		
		### 2D
		data = numpy.loadtxt(path2D)
		save0 = data[:,0]
	
		sPerp = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		sPara = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		sZ    = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
	
		for i in range(0,save0.size):
			iX = i/self._nbBinY2D_calcul
			iY = i%self._nbBinY2D_calcul
	
			idX = iX/int_binSize
			idY = iY/int_binSize
	
			xi2D[idX,idY,1] += data[i,0]
			sPerp[idX,idY]  += data[i,1]
			sPara[idX,idY]  += data[i,2]
			sZ[idX,idY]     += data[i,3]
	
		if (init):
			self._meanZ = numpy.sum(data[:,3])/numpy.sum(data[:,0])
			self._xi2D_grid = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,3) )
			cut = (xi2D[:,:,1]>0.)
			self._xi2D_grid[:,:,0][cut] = sPerp[cut] / xi2D[:,:,1][cut]
			self._xi2D_grid[:,:,1][cut] = sPara[cut] / xi2D[:,:,1][cut]
			self._xi2D_grid[:,:,2][cut] = sZ[cut] / xi2D[:,:,1][cut]

		### Fill arrays
		cut = (xi2D[:,:,1]>0.)
		xi2D[:,:,0][cut] = numpy.sqrt( sPerp[cut]**2. + sPara[cut]**2. )/xi2D[:,:,1][cut]
		xi2D[:,:,1][cut] /= coef
		xi2D[:,:,2][cut] = numpy.sqrt(xi2D[:,:,1][cut])
	
		return xiMu, xiWe, xi1D, xi2D
	def get_correlations_QSO(self):

		### Which estimator
		LS = False
		if (self._estimator=='LS'): LS = True
		
		### Coef to normalize
		coefDD = self._nd*(self._nd-1.)/2.
		if (LS): coefDR = self._nd*self._nr
		coefRR = self._nr*(self._nr-1.)/2.
	
		### Path
		path1D = self._path_to_txt_file_folder+'xi_QSO_QSO_Mu_QSO_'
		path2D = self._path_to_txt_file_folder+'xi_QSO_QSO_2D_QSO_'

		### DD
		dd_Mu, dd_We, dd_1D, dd_2D = self.read_data_QSO(path1D+'DD.txt', path2D+'DD.txt', coefDD, True)
		### RR
		rr_Mu, rr_We, rr_1D, rr_2D = self.read_data_QSO(path1D+'RR_0.txt', path2D+'RR_0.txt', coefRR)
		### DR
		if (LS): dr_Mu, dr_We, dr_1D, dr_2D = self.read_data_QSO(path1D+'DR_0.txt', path2D+'DR_0.txt', coefDR)
		
		for i in numpy.arange(1,self._nbRand):

			### RR
			tmp_rr_Mu, tmp_rr_We, tmp_rr_1D, tmp_rr_2D = self.read_data_QSO(path1D+'RR_'+str(i)+'.txt', path2D+'RR_'+str(i)+'.txt', coefRR)
			rr_Mu += tmp_rr_Mu
			rr_We += tmp_rr_We
			rr_1D += tmp_rr_1D
			rr_2D += tmp_rr_2D

			### DR
			if (LS):
				tmp_dr_Mu, tmp_dr_We, tmp_dr_1D, tmp_dr_2D = self.read_data_QSO(path1D+'DR_'+str(i)+'.txt', path2D+'DR_'+str(i)+'.txt', coefDR)
				dr_Mu += tmp_dr_Mu
				dr_We += tmp_dr_We
				dr_1D += tmp_dr_1D
				dr_2D += tmp_dr_2D
			
		
	
		### RR
		rr_Mu /= self._nbRand
		rr_We /= self._nbRand
		rr_1D /= self._nbRand
		rr_2D /= self._nbRand
		### DR
		if (LS):
			dr_Mu /= self._nbRand
			dr_We /= self._nbRand
			dr_1D /= self._nbRand
			dr_2D /= self._nbRand
	
		### Get the Landy-Saley estimator
		if (LS):	
			### Mu:
			xiMu = numpy.array(dd_Mu)
			cut = (xiMu[:,:,1]>0.)
			xiMu[:,:,2][cut] = (dd_Mu[:,:,2][cut] -2.*dr_Mu[:,:,2][cut] + rr_Mu[:,:,2][cut])/rr_Mu[:,:,2][cut]
			xiMu[:,:,3][cut] = numpy.sqrt(dd_Mu[:,:,2][cut]*coefDD)/(coefDD*rr_Mu[:,:,2][cut])

			### xiWe
			xiWe = numpy.array(dd_We)
			cut = (xiWe[:,:,1]>0.)
			xiWe[:,:,1][cut] = (dd_We[:,:,1][cut] -2.*dr_We[:,:,1][cut] + rr_We[:,:,1][cut])/rr_We[:,:,1][cut]
			xiWe[:,:,2][cut] = numpy.sqrt(dd_We[:,:,1][cut]*coefDD)/(coefDD*rr_We[:,:,1][cut])

			### 1D:
			xi1D = numpy.array(dd_1D)
			cut = (xi1D[:,1]>0.)
			xi1D[:,1][cut] = (dd_1D[:,1][cut]-2.*dr_1D[:,1][cut]+rr_1D[:,1][cut])/rr_1D[:,1][cut]
			xi1D[:,2][cut] = numpy.sqrt(dd_1D[:,1][cut]*coefDD)/(coefDD*rr_1D[:,1][cut])

			### 2D:
			xi2D = numpy.array(dd_2D)
			cut = (xi2D[:,:,1]>0.)
			xi2D[:,:,1][cut] = (dd_2D[:,:,1][cut]-2.*dr_2D[:,:,1][cut]+rr_2D[:,:,1][cut])/rr_2D[:,:,1][cut]
			xi2D[:,:,2][cut] = numpy.sqrt(dd_2D[:,:,1][cut]*coefDD)/(coefDD*rr_2D[:,:,1][cut])
	
		return xiMu, xiWe, xi1D, xi2D
	def read_data(self, path1D, path2D, selection=0, init=False):

		if (self._load_txt):
			xiMu, xiWe, xi1D, xi2D = self.get_correlations_QSO()
			xiMul = self.get_multipol(xiMu)

			if (init):
				### Save
				numpy.save(self._path_to_txt_file_folder+'xi_QSO_QSO_result_Mu',xiMu)
				numpy.save(self._path_to_txt_file_folder+'xi_QSO_QSO_result_We',xiWe)
				numpy.save(self._path_to_txt_file_folder+'xi_QSO_QSO_result_1D',xi1D)
				numpy.save(self._path_to_txt_file_folder+'xi_QSO_QSO_result_2D',xi2D)
				numpy.save(self._path_to_txt_file_folder+'xi_QSO_QSO_result_Multipol',xiMul)
		else:
			### init some objects
			path1D = self._path_to_txt_file_folder+'xi_QSO_QSO_Mu_QSO_'
			path2D = self._path_to_txt_file_folder+'xi_QSO_QSO_2D_QSO_'
			self.read_data_QSO(path1D+'DD.txt', path2D+'DD.txt', 1., True)

			xiMu = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_Mu.npy')
			xiWe = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_We.npy')
			xi1D = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_1D.npy')
			xi2D = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_2D.npy')

		return xiMu, xiWe, xi1D, xi2D

dic_class = {
	'minXi': 0.,
	'maxXi': 200.,
	'nbBin': 50,
	'nbBinM': 25,
	'nb_Sub_Sampling': 80,
	'size_bin_calcul_s': 1.,
	'size_bin_calcul_m': 0.02,
	'correlation': 'q_q',
	'path_to_txt_file_folder': '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/Box_000/Simu_000/Results_NicolasDistortion/',
	'f1': 'LYA',
	'f2': 'LYA',
	'q1': 'QSO',
	'q2': 'QSO',
	'name' : 'Data'
}
dic_Q = {
	'nb_random' : 10,
	'estimator' : 'LS',
	'path_to_cat' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/Box_000/Simu_000/Data/QSO_withRSD.fits',
	'load_from_txt' : False
}
dic_CAMB = {
	'mulpol_index' : 0,
	'start_fit'   : 40.,
	'end_fit'     : 180.,
	'b' : -1.,
	'roof' : 0.,
	'fix_roof_nul' : True,
	'guess_b' : False,
	'min_for_guess' : 20.,
	'max_for_guess' : 50.,
}
corr = Correlation3DQ(dic_class,dic_Q)

print corr._meanZ

corr.plot_slice_2d(0)
corr.plot_slice_2d(1)
corr.plot_slice_2d(None,0)

dic_CAMB = corr.fit_CAMB(None,dic_CAMB,False)
corr.plot_CAMB(None,dic_CAMB,0,False)
corr.plot_CAMB(None,dic_CAMB,1,False)
corr.plot_CAMB(None,dic_CAMB,2,False)
corr.plot_1d(0)
corr.plot_1d(1)
corr.plot_1d(2)
corr.plot_we(0)
corr.plot_we(1)
corr.plot_we(2)
corr.plot_2d(0)
corr.plot_2d(1)
corr.plot_2d(2)
corr.plot_mu(0)
corr.plot_mu(1)
corr.plot_mu(2)
corr.plot_multipol(0)
corr.plot_multipol(1)
corr.plot_multipol(2)








