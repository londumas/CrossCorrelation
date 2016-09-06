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
import copy
import sys
import scipy
import CAMB
import cosmolopy.perturbation as cp
from scipy.stats import chisqprob

### Perso lib
import myTools
import const_delta
import const


raw_dic_class = {
	'minXi'                   : 0.,
	'maxXi'                   : 200.,
	'nbBin'                   : 50,
	'nbBinM'                  : 25,
	'nb_Sub_Sampling'         : 80,
	'size_bin_calcul_s'       : 1.,
	'size_bin_calcul_m'       : 0.02,
	'nb_wedges'               : 3,
	'nb_multipole_max'        : 4,
	'remove_residuales'       : 'lambda_OBS',
	'path_to_txt_file_folder' : 'NOTHING',
	'correlation'             : 'q_f',
	'f1'                      : 'LYA',
	'f2'                      : '',
	'q1'                      : 'QSO',
	'q2'                      : '',
	'name'                    : 'SOMETHING'
}
raw_dic_CAMB = {
	'mulpol_index' : 0,
	'start_fit'   : 40.,
	'end_fit'     : 180.,
	'b' : -1.,
	'roof' : 0.,
	'fix_roof_nul' : True,
	'guess_b' : False,
	'min_for_guess' : 20.,
	'max_for_guess' : 50.,
	'chi^{2}' : -1.
}
raw_dic_simu = {
	'path_to_simu' : 'NOTHING',
	'nb_box' : 1,
	'nb_simu' : 1,
	'projected' : False,
	'with_metals_templates' : False,
	'raw' : True
}


class Correlation3D:

	def __init__(self, dic=None):
		"""
		
		correlationType:
			- 'q_q'
			- 'q_f' (default)
			- 'f_f'
			- 'f_f2'
	
		"""
	
		verbose__ = False
	
		if (dic is None):
			dic = copy.deepcopy(raw_dic_class)

		### folder where data are
		self._path_to_txt_file_folder = dic['path_to_txt_file_folder']

		### Name of the correlation
		self._name = dic['name']

		### Should we remove residuales and if yes, which one:
		self._remove_residuals = dic['remove_residuales']

		### forest and QSO name
		self._f1 = dic['f1']
		self._f2 = dic['f2']
		self._q1 = dic['q1']
		self._q2 = dic['q2']
		
		### Correlation type
		self._correlation = dic['correlation']
		if (self._correlation=='q_q'):
			self._label = '\\xi^{qq}'
			self._title = self._q1
			self._prefix = 'xi_QSO_QSO'
			self._middlefix = self._q1
		elif (self._correlation=='q_f'):
			self._label = '\\xi^{qf}'
			self._title = '\delta_{'+self._f1+'} \, - \, '+self._q1
			self._prefix = 'xi_delta_QSO'
			self._middlefix = self._f1+'_'+self._q1
		elif (self._correlation=='f_f'):
			self._label = '\\xi^{ff}'
			self._title = '\delta_{'+self._f1+'}'
			self._prefix = 'xi_A_delta_delta'
			self._middlefix = self._f1
		elif (self._correlation=='f_f2'):
			self._label = '\\xi^{ff}'
			self._title = '\delta_{'+self._f1+'} \, - \, \delta_{'+self._f2+'}'
			self._prefix = 'xi_A_delta_delta2'
			self._middlefix = self._f1+'_'+self._f2
		else:
			if (verbose__): print "  Correlation_3D::__init__::  ERROR:  'correlation' is incorrect "
			return
		path1D = self._path_to_txt_file_folder + self._prefix + '_Mu_' + self._middlefix + '.txt'
		path2D = self._path_to_txt_file_folder + self._prefix + '_2D_' + self._middlefix + '.txt'

		if (verbose__): print '  path1D = ',path1D

		### 1D
		self._min1D   = float(dic['minXi'])
		self._max1D   = float(dic['maxXi'])
		self._nbBin1D = int(dic['nbBin'])
		self._binSize = (self._max1D-self._min1D)/self._nbBin1D
		self._binSize_calcul = int( (self._max1D-self._min1D)/dic['size_bin_calcul_s'] )
		
		### 2D (X)
		self._minX2D   = self._min1D
		self._maxX2D   = self._max1D
		self._nbBinX2D = self._nbBin1D
		self._binSizeX2D = (self._maxX2D-self._minX2D)/self._nbBinX2D
		### 2D (Y)
		self._minY2D   = self._min1D
		self._maxY2D   = self._max1D
		self._nbBinY2D = self._nbBin1D
		self._nbBinY2D_calcul = self._binSize_calcul
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			self._minY2D    = -self._maxY2D
			self._nbBinY2D *= 2
			self._nbBinY2D_calcul *= 2
		self._nbBin2D  = self._nbBinX2D*self._nbBinY2D
		self._binSizeY2D = (self._maxY2D-self._minY2D)/self._nbBinY2D
		
		### Mu
		self._minM = 0.
		self._maxM = 1.
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			self._minM = -self._maxM
		self._nbBinM = dic['nbBinM']
		self._nbBinM_calcul = int( (self._maxM-self._minM)/dic['size_bin_calcul_m'] )

		### Wedges
		self._nb_wedges = dic['nb_wedges']
		if (self._nb_wedges==3):
			if (self._correlation=='q_q' or self._correlation=='f_f'):
				self._label_wedge = ['0.80 < \mu \leq 1.00',   '0.50 < \mu \leq 0.80',   '0.00 \leq \mu \leq 0.50']
			elif (self._correlation=='q_f' or self._correlation=='f_f2'):
				self._label_wedge = ['0.80 < |\mu| \leq 1.00', '0.50 < |\mu| \leq 0.80', '0.00 \leq |\mu| \leq 0.50']
			self.from_mu_index_to_wedges_index = self.from_mu_index_to_3_wedges_index
			self.from_mu_to_wedges_index = self.from_mu_to_3_wedges_index
		elif (self._nb_wedges==4):
			if (self._correlation=='q_q' or self._correlation=='f_f'):
				self._label_wedge = ['0.96 < \mu \leq 1.00', '0.80 < \mu \leq 0.96',   '0.50 < \mu \leq 0.80',   '0.00 \leq \mu \leq 0.50']
			elif (self._correlation=='q_f' or self._correlation=='f_f2'):
				self._label_wedge = ['0.96 < |\mu| \leq 1.00', '0.80 < |\mu| \leq 0.96', '0.50 < |\mu| \leq 0.80', '0.00 \leq |\mu| \leq 0.50']
			self.from_mu_index_to_wedges_index = self.from_mu_index_to_4_wedges_index
			self.from_mu_to_wedges_index = self.from_mu_to_4_wedges_index
			
		### Sub sampling
		self.nb_Sub_Sampling = dic['nb_Sub_Sampling']
		### Number maximum of multipols
		self._nb_multipole_max = dic['index_multipole_max']+1

		### Set attributes set after
		self._xi2D_grid = None
		self._xiMu_grid = None
		self._xiWe_grid = None
		self._xi1D_grid = None
		self._meanZ = None
		self._xiMul = None

		### Correlations data
		self._xiMu, self._xiWe, self._xi1D, self._xi2D = self.read_data(path1D, path2D,0,True)
		
		return
	
	def read_data(self, path1D, path2D, selection=0, init=False):
		'''
	
		selection:
			- 0: do all
			- 1: only 1D
			- 2: only 2D
	
		'''

		verbose__ = False

		if (init):
			if (verbose__): print "  ||                                  ||              ||                ||              ||"
			if (verbose__): print "  ||  selection                       ||  nb pairs    ||  sum weight    ||  <z>         ||"
			if (verbose__): print "  ||                                  ||              ||                ||              ||"

		remove_residuals = False
		if (self._remove_residuals == 'lambda_RF' or self._remove_residuals == 'lambda_OBS'):
			remove_residuals = True

		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
		xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4))
		xiWe = numpy.zeros(shape=(self._nbBin1D,self._nb_wedges,3))
		xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		
		int_binSize_R = int(self._binSize)
	
		if (selection==0 or selection==1):
			### Mu
			data = numpy.loadtxt(path1D)
	
			save0 = data[:,0]
			save1 = data[:,1]
			save2 = data[:,2]
			save3 = data[:,3]
			save4 = data[:,4]
			save5 = data[:,5]
			save6 = data[:,6]
			if (remove_residuals):
				save7 = data[:,7]
				save8 = data[:,8]
			else:
				save7 = numpy.zeros_like(data[:,6])
                                save8 = numpy.zeros_like(data[:,6])

			tmp_save0 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save1 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save2 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save3 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save4 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save5 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save6 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save7 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save8 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )

	
			tmp_save00 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save11 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save22 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save33 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save44 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save55 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save66 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save77 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save88 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
	
			tmp_save000 = numpy.zeros(self._nbBin1D)
			tmp_save111 = numpy.zeros(self._nbBin1D)
			tmp_save222 = numpy.zeros(self._nbBin1D)
			tmp_save333 = numpy.zeros(self._nbBin1D)
			tmp_save444 = numpy.zeros(self._nbBin1D)
			tmp_save555 = numpy.zeros(self._nbBin1D)
			tmp_save666 = numpy.zeros(self._nbBin1D)
			tmp_save777 = numpy.zeros(self._nbBin1D)
			tmp_save888 = numpy.zeros(self._nbBin1D)
	
			binSizeY = self._nbBinM_calcul/self._nbBinM
		
			for i in range(0,save0.size ):
				iX = i/self._nbBinM_calcul
				iY = i%self._nbBinM_calcul
	
				### for mu
				idX = iX/int_binSize_R
				idY = iY/binSizeY
				if (idX>=self._nbBin1D):
					print '  OUT OF BOUNDS'
					continue
	
				tmp_save0[idX][idY] += save0[i]
				tmp_save1[idX][idY] += save1[i]
				tmp_save2[idX][idY] += save2[i]
				tmp_save3[idX][idY] += save3[i]
				tmp_save4[idX][idY] += save4[i]
				tmp_save5[idX][idY] += save5[i]
				tmp_save6[idX][idY] += save6[i]
				tmp_save7[idX][idY] += save7[i]
				tmp_save8[idX][idY] += save8[i]
				
				### for wedges
				idY = self.from_mu_index_to_wedges_index(iY)
	
				tmp_save00[idX][idY] += save0[i]
				tmp_save11[idX][idY] += save1[i]
				tmp_save22[idX][idY] += save2[i]
				tmp_save33[idX][idY] += save3[i]
				tmp_save44[idX][idY] += save4[i]
				tmp_save55[idX][idY] += save5[i]
				tmp_save66[idX][idY] += save6[i]
				tmp_save77[idX][idY] += save7[i]
				tmp_save88[idX][idY] += save8[i]

				### for xi1D
				tmp_save000[idX] += save0[i]
				tmp_save111[idX] += save1[i]
				tmp_save222[idX] += save2[i]
				tmp_save333[idX] += save3[i]
				tmp_save444[idX] += save4[i]
				tmp_save555[idX] += save5[i]
				tmp_save666[idX] += save6[i]
				tmp_save777[idX] += save7[i]
				tmp_save888[idX] += save8[i]

			cut = (tmp_save6>1.)
			xiMu[:,:,0][cut] = tmp_save2[cut]/tmp_save5[cut]
			xiMu[:,:,1][cut] = tmp_save3[cut]/tmp_save5[cut]
			xiMu[:,:,2][cut] = tmp_save0[cut]/tmp_save5[cut]
			xiMu[:,:,3][cut] = numpy.sqrt( (tmp_save1[cut]/tmp_save5[cut] - xiMu[:,:,2][cut]*xiMu[:,:,2][cut])/tmp_save6[cut])
			if (self._remove_residuals == 'lambda_RF'):
				xiMu[:,:,2][cut]-= tmp_save7[cut]/tmp_save5[cut]
			elif (self._remove_residuals == 'lambda_OBS'):
				xiMu[:,:,2][cut]-= tmp_save8[cut]/tmp_save5[cut]
			if (init):
				self._xiMu_grid = numpy.zeros( shape=(self._nbBin1D,self._nbBinM,3) )
				self._xiMu_grid[:,:,0][cut] = xiMu[:,:,0][cut]
				self._xiMu_grid[:,:,1][cut] = xiMu[:,:,1][cut]
				self._xiMu_grid[:,:,2][cut] = tmp_save4[cut]/tmp_save5[cut]

			cut = (tmp_save66>1.)
			xiWe[:,:,0][cut] = tmp_save22[cut]/tmp_save55[cut]
			xiWe[:,:,1][cut] = tmp_save00[cut]/tmp_save55[cut]
			xiWe[:,:,2][cut] = numpy.sqrt( (tmp_save11[cut]/tmp_save55[cut] - xiWe[:,:,1][cut]*xiWe[:,:,1][cut])/tmp_save66[cut] )
			if (self._remove_residuals == 'lambda_RF'):
				xiWe[:,:,1][cut]-= tmp_save77[cut]/tmp_save55[cut]
			elif (self._remove_residuals == 'lambda_OBS'):
				xiWe[:,:,1][cut]-= tmp_save88[cut]/tmp_save55[cut]
			if (init):
				self._xiWe_grid = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges,2) )
				self._xiWe_grid[:,:,0][cut] = xiWe[:,:,0][cut]
				self._xiWe_grid[:,:,1][cut] = tmp_save44[cut]/tmp_save55[cut]

			cut = (tmp_save666>1.)
			xi1D[:,0][cut] = tmp_save222[cut]/tmp_save555[cut]
			xi1D[:,1][cut] = tmp_save000[cut]/tmp_save555[cut]
			xi1D[:,2][cut] = numpy.sqrt( (tmp_save111[cut]/tmp_save555[cut] - xi1D[:,1][cut]*xi1D[:,1][cut])/tmp_save666[cut] )
			if (self._remove_residuals == 'lambda_RF'):
				xi1D[:,1][cut]-= tmp_save777[cut]/tmp_save555[cut]
			elif (self._remove_residuals == 'lambda_OBS'):
				xi1D[:,1][cut]-= tmp_save888[cut]/tmp_save555[cut]
			if (init):
				self._xi1D_grid = numpy.zeros( shape=(self._nbBin1D,2) )
				self._xi1D_grid[:,0][cut] = xi1D[:,0][cut]
				self._xi1D_grid[:,1][cut] = tmp_save444[cut]/tmp_save555[cut]

			if (init):
				if (verbose__): print "  ||  |s| < %u                       ||  %1.4e  ||  %1.4e    ||  %1.4e  ||" % (int(self._max1D), numpy.sum(tmp_save666), numpy.sum(tmp_save555), numpy.sum(data[:,4])/numpy.sum(data[:,5]))
			
		int_binSize_X = int(self._binSizeX2D)
		int_binSize_Y = int(self._binSizeY2D)

		if (selection==0 or selection==2):
			### 2D
			data = numpy.loadtxt(path2D)
	
			save0 = data[:,0]
			save1 = data[:,1]
			save2 = data[:,2]
			save3 = data[:,3]
			save4 = data[:,4]
			save5 = data[:,5]
			save6 = data[:,6]
			if (remove_residuals):
				save7 = data[:,7]
				save8 = data[:,8]
			else:
				save7 = numpy.zeros_like(data[:,6])
                                save8 = numpy.zeros_like(data[:,6])

			tmp_save0  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save1  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save2  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save3  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save4  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save5  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save6  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save7  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save8  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )

	
			for i in range( 0,save0.size ):
				iX = i/self._nbBinY2D_calcul
				iY = i%self._nbBinY2D_calcul
	
				idX = iX/int_binSize_X
				idY = iY/int_binSize_Y
				if (idX>=self._nbBinX2D):
					print '  OUT OF BOUNDS'
					continue
				if (idY>=self._nbBinY2D):
					print '  OUT OF BOUNDS'
					continue

				tmp_save0[idX][idY] += save0[i]
				tmp_save1[idX][idY] += save1[i]
				tmp_save2[idX][idY] += save2[i]
				tmp_save3[idX][idY] += save3[i]
				tmp_save4[idX][idY] += save4[i]
				tmp_save5[idX][idY] += save5[i]
				tmp_save6[idX][idY] += save6[i]
				tmp_save7[idX][idY] += save7[i]
				tmp_save8[idX][idY] += save8[i]

			cut = (tmp_save6>1.)
			xi2D[:,:,0][cut] = numpy.sqrt( (tmp_save2[cut]/tmp_save5[cut])**2. + (tmp_save3[cut]/tmp_save5[cut])**2. )
			xi2D[:,:,1][cut] = tmp_save0[cut] / tmp_save5[cut]
			xi2D[:,:,2][cut] = numpy.sqrt( (tmp_save1[cut]/tmp_save5[cut] - xi2D[:,:,1][cut]*xi2D[:,:,1][cut])/tmp_save6[cut] )
			if (self._remove_residuals == 'lambda_RF'):
				xi2D[:,:,1][cut]-= tmp_save7[cut]/tmp_save5[cut]
			elif (self._remove_residuals == 'lambda_OBS'):
				xi2D[:,:,1][cut]-= tmp_save8[cut]/tmp_save5[cut]

			if (init):
				self._meanZ = numpy.sum(data[:,4])/numpy.sum(data[:,5])
				self._xi2D_grid = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,3) )
				self._xi2D_grid[:,:,0][cut] = tmp_save2[cut] / tmp_save5[cut]
				self._xi2D_grid[:,:,1][cut] = tmp_save3[cut] / tmp_save5[cut]
				self._xi2D_grid[:,:,2][cut] = tmp_save4[cut] / tmp_save5[cut]

				if (verbose__): print "  ||  |s_perp| and |s_paral| < %u    ||  %1.4e  ||  %1.4e    ||  %1.4e  ||" % (int(self._max1D),numpy.sum(tmp_save6), numpy.sum(tmp_save5), self._meanZ)
				if (verbose__): print "  ||                                  ||              ||                ||              ||"
				if (verbose__): print

		return xiMu, xiWe, xi1D, xi2D
	def read_data_from_BAOFIT_data_file(self, path2D):

		### 2D
		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
		data  = numpy.loadtxt(path2D+'.data')
		for el in data:
			i = int(el[0])
			ix = i%self._nbBinX2D
			iy = i/self._nbBinX2D
			xi2D[ix,iy,1] = el[1]

		'''
		### cov
		data = numpy.loadtxt(path2D+'.cov')
		cut  = (data[:,0]==data[:,1])
		idx1 = data[:,0][cut].astype(int)
		idx2 = data[:,1][cut].astype(int)
		var  = data[:,2][cut]
		print numpy.array( zip(idx1,idx2,var) )
		xi2D[idx1,idx1,2] = numpy.sqrt(var)
		'''

		self._xi2D[:,:,1] = xi2D[:,:,1]

		### 1D
		xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		for i in numpy.arange(self._nbBinX2D):
			for j in numpy.arange(self._nbBinY2D):
				if (self._xi2D[i,j,1]==0.): continue
				else:
					if (self._xi2D[i,j,0]>=self._max1D): continue
					idx = int( self._xi2D[i,j,0]/self._binSize )
					ivar = 1./numpy.power(self._xi2D[i,j,2],2.)
					xi1D[idx,0] += ivar*self._xi2D[i,j,0]
					xi1D[idx,1] += ivar*self._xi2D[i,j,1]
					xi1D[idx,2] += ivar
		cut = (xi1D[:,2]>0.)
		xi1D[:,0][cut] /= xi1D[:,2][cut]
		xi1D[:,1][cut] /= xi1D[:,2][cut]
		xi1D[:,2][cut]  = 1./numpy.sqrt(xi1D[:,2][cut])

		self._xi1D[:,:2] = xi1D[:,:2]
		

		return
	def read_metal_model(self, lineName,selection=0):

		path1D = self._path_to_txt_file_folder + self._prefix + '_Metals_Models_Mu_' + self._middlefix + '_' + lineName + '.txt'
		path2D = self._path_to_txt_file_folder + self._prefix + '_Metals_Models_2D_' + self._middlefix + '_' + lineName + '.txt'

		nb_multipol_calculated_templates = 3

		'''
	
		selection:
			- 0: do all
			- 1: only 1D
			- 2: only 2D
	
		'''
		print path1D

		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3,nb_multipol_calculated_templates))
		xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4,nb_multipol_calculated_templates))
		xiWe = numpy.zeros(shape=(self._nbBin1D,self._nb_wedges,3,nb_multipol_calculated_templates))
		xi1D = numpy.zeros(shape=(self._nbBin1D,3,nb_multipol_calculated_templates))
		
		int_binSize = int(self._binSize)
	
		if (selection==0 or selection==1):
			### Mu
			data = numpy.loadtxt(path1D)
			save0 = data[:,0]
	
			tmp_save0 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM,nb_multipol_calculated_templates) )
			tmp_save1 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM,nb_multipol_calculated_templates) )
			tmp_save2 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save3 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save5 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save6 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
	
			tmp_save00 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges,nb_multipol_calculated_templates) )
			tmp_save11 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges,nb_multipol_calculated_templates) )
			tmp_save22 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save33 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save55 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
			tmp_save66 = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges) )
	
			tmp_save000 = numpy.zeros( shape=(self._nbBin1D,nb_multipol_calculated_templates) )
			tmp_save111 = numpy.zeros( shape=(self._nbBin1D,nb_multipol_calculated_templates) )
			tmp_save222 = numpy.zeros(self._nbBin1D)
			tmp_save333 = numpy.zeros(self._nbBin1D)
			tmp_save555 = numpy.zeros(self._nbBin1D)
			tmp_save666 = numpy.zeros(self._nbBin1D)
	
			binSizeY = self._nbBinM_calcul/self._nbBinM
		
			for i in range(0,save0.size ):
				iX = i/self._nbBinM_calcul
				iY = i%self._nbBinM_calcul
	
				### for mu
				idX = iX/int_binSize
				idY = iY/binSizeY
	
				tmp_save0[idX,idY,0] += data[i,0]
				tmp_save0[idX,idY,1] += data[i,1]
				tmp_save0[idX,idY,2] += data[i,2]
				tmp_save2[idX,idY]   += data[i,3]
				tmp_save3[idX,idY]   += data[i,4]
				tmp_save5[idX,idY]   += data[i,6]
				tmp_save6[idX,idY]   += data[i,7]

				### for wedges
				idY = self.from_mu_index_to_wedges_index(iY)
	
				tmp_save00[idX,idY,0] += data[i,0]
				tmp_save00[idX,idY,1] += data[i,1]
				tmp_save00[idX,idY,2] += data[i,2]
				tmp_save22[idX,idY]   += data[i,3]
				tmp_save33[idX,idY]   += data[i,4]
				tmp_save55[idX,idY]   += data[i,6]
				tmp_save66[idX,idY]   += data[i,7]
	
				### for xi1D
				tmp_save000[idX,0] += data[i,0]
				tmp_save000[idX,1] += data[i,1]
				tmp_save000[idX,2] += data[i,2]
				tmp_save222[idX]   += data[i,3]
				tmp_save333[idX]   += data[i,4]
				tmp_save555[idX]   += data[i,6]
				tmp_save666[idX]   += data[i,7]

			
			cut = (tmp_save6>1.)
			for i in numpy.arange(0,nb_multipol_calculated_templates):
				xiMu[:,:,0,i][cut] = tmp_save2[cut]/tmp_save5[cut]
				xiMu[:,:,1,i][cut] = tmp_save3[cut]/tmp_save5[cut]
				xiMu[:,:,2,i][cut] = tmp_save0[:,:,i][cut]/tmp_save5[cut]
				xiMu[:,:,3,i][cut] = numpy.abs(xiMu[:,:,2,i][cut])/numpy.sqrt(tmp_save6[cut])

			cut = (tmp_save66>1.)
			for i in numpy.arange(0,nb_multipol_calculated_templates):
				xiWe[:,:,0,i][cut] = tmp_save22[cut]/tmp_save55[cut]
				xiWe[:,:,1,i][cut] = tmp_save00[:,:,i][cut]/tmp_save55[cut]
				xiWe[:,:,2,i][cut] = numpy.abs(xiWe[:,:,1,i][cut])/numpy.sqrt(tmp_save66[cut])
			
			cut = (tmp_save666>1.)
			for i in numpy.arange(0,nb_multipol_calculated_templates):
				xi1D[:,0,i][cut] = tmp_save222[cut]/tmp_save555[cut]
				xi1D[:,1,i][cut] = tmp_save000[:,i][cut]/tmp_save555[cut]
				xi1D[:,2,i][cut] = numpy.abs(xi1D[:,1,i][cut])/numpy.sqrt(tmp_save666[cut])

			#plt.plot(xi1D[:,0,0][cut], -xi1D[:,1,0][cut])
			#plt.grid()
			#plt.show()
		
		if (selection==0 or selection==2):
			### 2D
			data = numpy.loadtxt(path2D)
			save0 = data[:,0]

			tmp_save0  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,nb_multipol_calculated_templates) )
			tmp_save1  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,nb_multipol_calculated_templates) )
			tmp_save2  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save3  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save5  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save6  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
	
			for i in range( 0,save0.size ):
				iX = i/self._nbBinY2D_calcul
				iY = i%self._nbBinY2D_calcul
	
				idX = iX/int_binSize
				idY = iY/int_binSize
	
				tmp_save0[idX,idY,0] += data[i,0]
				tmp_save0[idX,idY,1] += data[i,1]
				tmp_save0[idX,idY,2] += data[i,2]
				tmp_save2[idX,idY] += data[i,3]
				tmp_save3[idX,idY] += data[i,4]
				tmp_save5[idX,idY] += data[i,6]
				tmp_save6[idX,idY] += data[i,7]

			cut = (tmp_save6>1.)
			for i in numpy.arange(0,nb_multipol_calculated_templates):
				xi2D[:,:,0,i][cut] = numpy.sqrt( (tmp_save2[cut]/tmp_save5[cut])**2. + (tmp_save3[cut]/tmp_save5[cut])**2. )
				xi2D[:,:,1,i][cut] = tmp_save0[:,:,i][cut] / tmp_save5[cut]
				xi2D[:,:,2,i][cut] = numpy.abs(xi2D[:,:,1,i][cut])/numpy.sqrt(tmp_save6[cut])

		return xiMu, xiWe, xi1D, xi2D
	def from_mu_index_to_3_wedges_index(self,iY):

		### for wedges
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			if (iY<10 or iY>90):
				idY = 0
			elif ( (iY>=10 and iY<25) or (iY<=90 and iY>75) ):
				idY = 1
			else:
				idY = 2
		elif (self._correlation=='q_q' or self._correlation=='f_f'):
			if (iY>40):
				idY = 0
			elif (iY>25 and iY<=40):
				idY = 1
			else:
				idY = 2

		return idY
	def from_mu_to_3_wedges_index(self,mu):

		mu = numpy.abs(mu)
		we_index = numpy.ones(mu.size)
		we_index[ (mu>0.8) ]  = 0
		we_index[ (mu<=0.5) ] = 2
 
		return we_index
	def from_mu_index_to_4_wedges_index(self,iY):

		### for wedges
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			if (iY<3 or iY>97):
				idY = 0
			elif ( (iY>=3 and iY<10) or (iY<=97 and iY>90) ):
				idY = 1
			elif ( (iY>=10 and iY<25) or (iY<=90 and iY>75) ):
				idY = 2
			else:
				idY = 3
		elif (self._correlation=='q_q' or self._correlation=='f_f'):
			if (iY>47):
				idY = 0
			elif (iY>40):
				idY = 1
			elif (iY>25 and iY<=40):
				idY = 2
			else:
				idY = 3

		return idY
	def from_mu_to_4_wedges_index(self,mu):
		
		mu = numpy.abs(mu)
		we_index = numpy.zeros(mu.size)
		we_index[ numpy.logical_and(mu>0.8,mu<=0.96) ] = 1
		we_index[ numpy.logical_and(mu>0.5,mu<=0.8) ] = 2
		we_index[ (mu<=0.5) ] = 3
 
		return we_index
	def empty_data(self):

		### Set attributes set after
		self._meanZ = 0.
		self._xiMul = numpy.zeros(shape=(self._nbBin1D,self._nb_multipole_max,3))

		### Correlations data
		self._xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4))
		self._xiWe = numpy.zeros(shape=(self._nbBin1D,self._nb_wedges,3))
		self._xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		self._xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))

		return
	def save_list_realisation(self, realisation_type, nb_realisation):

		if ( realisation_type=='subsampling' and (nb_realisation != self.nb_Sub_Sampling) ):
			print "  Correlation_3D::saveListReal::  WARNING:  nb_realisation != self.nb_Sub_Sampling"

		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,nb_realisation) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges,nb_realisation) )
		list1D       = numpy.zeros( shape=(self._nbBin1D,nb_realisation) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,nb_realisation) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,self._nb_multipole_max,nb_realisation) )

		path1D = self._path_to_txt_file_folder + self._prefix + '_Mu_' + self._middlefix + '_' + realisation_type + '_'
		path2D = self._path_to_txt_file_folder + self._prefix + '_2D_' + self._middlefix + '_' + realisation_type + '_'
		pathToSave = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_'

		nb = 0
		for i in numpy.arange(nb_realisation):

			try:
				xiMu, xiWe, xi1D, xi2D = self.read_data(path1D+str(i)+'.txt',path2D+str(i)+'.txt')
				list1D[:,nb]         = xi1D[:,1]
				list2D[:,nb]         = xi2D[:,:,1].flatten()
				listMu[:,nb]         = xiMu[:,:,2].flatten()
				listWe[:,:,nb]       = xiWe[:,:,1]
				listMultipol[:,:,nb] = self.get_multipol(xiMu)[:,:,1]

				nb += 1
			except:
				print '   ERROR:: ', path1D+str(i)+'.txt',path2D+str(i)+'.txt'
			

                listMu       = listMu[:,:nb]
                listWe       = listWe[:,:,:nb]
                list1D       = list1D[:,:nb]
                list2D       = list2D[:,:nb]
                listMultipol = listMultipol[:,:,:nb]

		numpy.save(pathToSave+'list_Mu',listMu)
		numpy.save(pathToSave+'list_We',listWe)
		numpy.save(pathToSave+'list_1D',list1D)
		numpy.save(pathToSave+'list_2D',list2D)
		numpy.save(pathToSave+'list_Multipol',listMultipol)

		covMu = numpy.cov(listMu)
		cov1D = numpy.cov(list1D)
		cov2D = numpy.cov(list2D)
		if (realisation_type=='subsampling'):
			covMu /= nb_realisation
			cov1D /= nb_realisation
			cov2D /= nb_realisation

		numpy.save(pathToSave+'cov_Mu',covMu)
		numpy.save(pathToSave+'cov_1D',cov1D)
		numpy.save(pathToSave+'cov_2D',cov2D)
	
		return
	def save_list_realisation_simulation(self, dic_class, dic_simu, distortion_matrix=False):

		'''
			if distortion_matrix is True, it will save the mean, not the list.
		'''

		pathToSave = dic_simu['path_to_simu'] + 'Results'
		pathToSave += dic_simu['prefix']
		pathToSave += '/'
		pathToSave += self._prefix + '_' + self._middlefix + '_result_'

		nb_realisation = dic_simu['nb_box']*dic_simu['nb_simu']

		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,nb_realisation) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges,nb_realisation) )
		list1D       = numpy.zeros( shape=(self._nbBin1D,nb_realisation) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,nb_realisation) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,self._nb_multipole_max,nb_realisation) )
		listGrid     = numpy.zeros( shape=(self._nbBin2D,3,nb_realisation) )
		list_mean_z  = numpy.zeros( shape=(nb_realisation) )
		if (distortion_matrix): dist = numpy.zeros( shape=(self._nbBin2D,self._nbBin2D))

		nb = 0
		for i in numpy.arange(dic_simu['nb_box']):
			for j in numpy.arange(dic_simu['nb_simu']):

				raw = dic_simu['path_to_simu'] + 'Box_00' + str(i) + '/Simu_00' + str(j) +'/Results'
				raw += dic_simu['prefix']
				raw += '/'
				dic_class['path_to_txt_file_folder'] = raw
				path_to_distortion_matrix = dic_class['path_to_txt_file_folder'] + self._prefix + '_distortionMatrix_2D_'+ self._middlefix + '.txt'

				try:
					print i, j
					corr = Correlation3D(dic_class)
					if (distortion_matrix): tmp_dist = numpy.loadtxt(path_to_distortion_matrix)
				except:
					print 'ERROR: can not read file : ', i, j
					continue

				listMu[:,nb]         = corr._xiMu[:,:,2].flatten()
				listWe[:,:,nb]       = corr._xiWe[:,:,1]
				list1D[:,nb]         = corr._xi1D[:,1]
				list2D[:,nb]         = corr._xi2D[:,:,1].flatten()
				listMultipol[:,:,nb] = corr.get_multipol(corr._xiMu)[:,:,1]

				listGrid[:,0,nb] = corr._xi2D_grid[:,:,0].flatten()
				listGrid[:,1,nb] = corr._xi2D_grid[:,:,1].flatten()
				listGrid[:,2,nb] = corr._xi2D_grid[:,:,2].flatten()
				list_mean_z[nb]  = corr._meanZ
				if (distortion_matrix): dist += tmp_dist

				nb += 1

		listMu       = listMu[:,:nb]
		listWe       = listWe[:,:,:nb]
		list1D       = list1D[:,:nb]
		list2D       = list2D[:,:nb]
		listMultipol = listMultipol[:,:,:nb]
		listGrid     = listGrid[:,:,:nb]
		list_mean_z  = list_mean_z[:nb]
		if (distortion_matrix): dist /= nb

		numpy.save(pathToSave+'list_Mu',listMu)
		numpy.save(pathToSave+'list_We',listWe)
		numpy.save(pathToSave+'list_1D',list1D)
		numpy.save(pathToSave+'list_2D',list2D)
		numpy.save(pathToSave+'list_Multipol',listMultipol)
		numpy.save(pathToSave+'list_Grid',listGrid)
		numpy.save(pathToSave+'list_MeanZ',list_mean_z)

		covMu = numpy.cov(listMu)
		cov1D = numpy.cov(list1D)
		cov2D = numpy.cov(list2D)

		numpy.save(pathToSave+'cov_Mu',covMu)
		numpy.save(pathToSave+'cov_1D',cov1D)
		numpy.save(pathToSave+'cov_2D',cov2D)
		if (distortion_matrix): numpy.save(pathToSave+'dist_2D',dist)

		return
	def save_mean_realisation_simulation_correlation_matrix(self, dic_simu):

		### Where to save the resulted mean
		pathToSave = dic_simu['path_to_simu'] + 'Results'
		pathToSave += dic_simu['prefix']
		pathToSave += '/'
		pathToSave += self._prefix + '_' + self._middlefix + '_result_cor'
		
		cor_1D = numpy.zeros( shape=(self._nbBin1D,self._nbBin1D) )
		cor_2D = numpy.zeros( shape=(self._nbBin2D,self._nbBin2D) )

		### Stack the correlation matrix
		nb = 0
		for i in range(0,dic_simu['nb_box']):
			for j in range(0,dic_simu['nb_simu']):

				raw  = dic_simu['path_to_simu'] + 'Box_00' + str(i) + '/Simu_00' + str(j) +'/Results'
				raw += dic_simu['prefix']
				raw += '/'
				raw += self._prefix + '_' + self._middlefix
				raw += '_subsampling_cov_1D.npy'

				try:
					tmp_cov  = numpy.load(raw)
					tmp_diag = diag = numpy.diag(tmp_cov)
					if ( tmp_diag[ (tmp_diag<=0.) ].size != 0  ):
						print '  BUG IN : ', i, j, tmp_diag[ (tmp_diag<=0.) ].size, tmp_diag[ (tmp_diag<=0.) ]
						continue
					cor_1D += myTools.getCorrelationMatrix( tmp_cov )
				except:
					print '  BUG IN : ', i, j, raw
					continue

				raw  = dic_simu['path_to_simu'] + 'Box_00' + str(i) + '/Simu_00' + str(j) +'/Results'
                                raw += dic_simu['prefix']
                                raw += '/'
                                raw += self._prefix + '_' + self._middlefix
                                raw += '_subsampling_cov_2D.npy'

                                try:
                                        tmp_cov  = numpy.load(raw)
                                        tmp_diag = diag = numpy.diag(tmp_cov)
                                        if ( tmp_diag[ (tmp_diag<=0.) ].size != 0  ):
                                                print '  BUG IN : ', i, j, tmp_diag[ (tmp_diag<=0.) ].size, tmp_diag[ (tmp_diag<=0.) ]
                                                continue
                                        cor_2D += myTools.getCorrelationMatrix( tmp_cov )
                                except:
                                        print '  BUG IN : ', i, j, raw
                                        continue

				nb += 1

		print nb

		cor_1D /= nb
		cor_2D /= nb
		myTools.plot2D(cor_1D)
                myTools.plot2D(cor_2D)

		### Save the mean
		print '  saving to = ', pathToSave
		numpy.save( pathToSave+'_mean_subsampling_1D', cor_1D)
		numpy.save( pathToSave+'_mean_subsampling_2D', cor_2D)

		### model the mean and fit
                mean_cor_1D = myTools.plotCovar([cor_1D], ['a'])[0]
                myTools.plot2D(mean_cor_1D)
                ### Save the mean fitted
                print '  saving to = ', pathToSave
                numpy.save( pathToSave+'_mean_subsampling_from_fit_1D', mean_cor_1D)
                ### Look at the shape of the mean fitted reduced correlation
                myTools.plotCovar([mean_cor_1D], ['a'])[0]


		### model the mean and fit
		mean_cor_2D = myTools.plotCovar([cor_2D], ['a'],self._nbBinX2D,self._nbBinY2D)[0]
		myTools.plot2D(mean_cor_2D)
		### Save the mean fitted
		print '  saving to = ', pathToSave
		numpy.save( pathToSave+'_mean_subsampling_from_fit_2D', mean_cor_2D)
		### Look at the shape of the mean fitted reduced correlation
		myTools.plotCovar([mean_cor_2D], ['a'],self._nbBinX2D,self._nbBinY2D)[0]

		return 
	def set_values_on_mean_simulation_or_realisation_type(self, dic_simu, realisation_type=None):

		if (realisation_type is None):
			path_to_load = dic_simu['path_to_simu'] + 'Results'
			path_to_load += dic_simu['prefix']
			path_to_load += '/'
			if (self._prefix=='xi_QSO_QSO'):
				path_to_load += self._prefix + '_result_'
			else:
				path_to_load += self._prefix + '_' + self._middlefix + '_result_'
			print path_to_load
		else:
                	if (dic_simu is None):
                	        path_to_load = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_'
                	else:
                	        path_to_load = dic_simu['path_to_simu'] + 'Results/' + self._prefix + '_result_'
                	        realisation_type = 'simulation'

		### mu
		listMu = numpy.load(path_to_load+'list_Mu.npy')
		nb_realisation = listMu[0,:].size
		self._xiMu[:,:,2] = myTools.convert1DTo2D( numpy.mean(listMu,axis=1) ,self._nbBin1D,self._nbBinM)
		self._xiMu[:,:,3] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.cov(listMu))/nb_realisation ),self._nbBin1D,self._nbBinM)
		### we
		listWe = numpy.load(path_to_load+'list_We.npy')
		for i in numpy.arange(self._nb_wedges):
			self._xiWe[:,i,1] = numpy.mean(listWe[:,i,:],axis=1)
			self._xiWe[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listWe[:,i,:] )/nb_realisation))
		### 1D
		list1D = numpy.load(path_to_load+'list_1D.npy')
		self._xi1D[:,1] = numpy.mean(list1D,axis=1)
		self._xi1D[:,2] = numpy.sqrt( numpy.diag(numpy.cov(list1D))/nb_realisation )
		### 2D
		list2D = numpy.load(path_to_load+'list_2D.npy')
		self._xi2D[:,:,1] = myTools.convert1DTo2D( numpy.mean(list2D,axis=1) ,self._nbBinX2D, self._nbBinY2D)
		self._xi2D[:,:,2] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.cov(list2D))/nb_realisation ),self._nbBinX2D, self._nbBinY2D)
		### Multipol
		listMultipol = numpy.load(path_to_load+'list_Multipol.npy')
		self._xiMul = self.get_multipol(self._xiMu)
		for i in numpy.arange(self._nb_multipole_max):
			self._xiMul[:,i,1] = numpy.mean(listMultipol[:,i,:],axis=1)
			self._xiMul[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listMultipol[:,i,:] )/nb_realisation))
		
		if (realisation_type is None):
			### Grid
			listGrid = numpy.load(path_to_load+'list_Grid.npy')
			self._xi2D_grid[:,:,0] = myTools.convert1DTo2D( numpy.mean(listGrid[:,0,:],axis=1),self._nbBinX2D,self._nbBinY2D )
			self._xi2D_grid[:,:,1] = myTools.convert1DTo2D( numpy.mean(listGrid[:,1,:],axis=1),self._nbBinX2D,self._nbBinY2D )
			self._xi2D_grid[:,:,2] = myTools.convert1DTo2D( numpy.mean(listGrid[:,2,:],axis=1),self._nbBinX2D,self._nbBinY2D )
			self._meanZ = numpy.mean(numpy.load(path_to_load+'list_MeanZ.npy'))

		print '  nb_realisation = ', nb_realisation

		return
	def set_error_on_covar_matrix(self, realisation_type, dic_simu=None):

		if (dic_simu is None):
			raw_path = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_'
		else:
			raw_path = dic_simu['path_to_simu'] + 'Results/' + self._prefix + '_result_'
			realisation_type = 'simulation'

		### mu
		path = raw_path + 'cov_Mu.npy'
		self._xiMu[:,:,3] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path)) ),self._nbBin1D,self._nbBinM)
		### we
		path = raw_path + 'list_We.npy'
		listWe = numpy.load(path)
		if (realisation_type=='subsampling'): coef = listWe[0,0,:].size
		else: coef = 1.
		for i in numpy.arange(self._nb_wedges):
			self._xiWe[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listWe[:,i,:] )/coef))
		### 1D
		path = raw_path + 'cov_1D.npy'
		self._xi1D[:,2]   = numpy.sqrt( numpy.diag(numpy.load(path)) )
		### 2D
		path = raw_path + 'cov_2D.npy'
		self._xi2D[:,:,2] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path)) ),self._nbBinX2D, self._nbBinY2D)
		### Multipol
		if (self._xiMul is None): self._xiMul = self.get_multipol(self._xiMu)
		path = raw_path + 'list_Multipol.npy'
		listMultipol = numpy.load(path)
		for i in numpy.arange(self._nb_multipole_max):
			self._xiMul[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listMultipol[:,i,:] )/coef))

		return
	def multiply_by_constant(self, const):

		### Data
		self._xiMu[:,:,2] *= const
		self._xiWe[:,:,1] *= const
		self._xi1D[:,1]   *= const
		self._xi2D[:,:,1] *= const

		### Errors
		self._xiMu[:,:,3] *= const
		self._xiWe[:,:,2] *= const
		self._xi1D[:,2]   *= const
		self._xi2D[:,:,2] *= const

		return
	def apply_distortion_matrix(self):

		### 1d
		path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_1D_'+ self._middlefix + '.txt'
		matrix = numpy.loadtxt(path)
		self._xi1D[:,1] = numpy.dot(matrix,self._xi1D[:,1])
		### 2d
		path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_2D_'+ self._middlefix + '.txt'
		matrix = numpy.loadtxt(path)
		corr = myTools.convert2DTo1D(self._xi2D[:,:,1], self._nbBinX2D, self._nbBinY2D)
		corr = numpy.dot(matrix,corr)
		self._xi2D[:,:,1] = myTools.convert1DTo2D(corr, self._nbBinX2D, self._nbBinY2D)
	
		return
	def get_multipol(self, xiMu, dic=None):
		"""
			https://fr.wikipedia.org/wiki/Polyn%C3%B4me_de_Legendre
		"""
	
		### Array with name of variable
		xi         = numpy.array(['xi']*10)
		index      = numpy.arange(0,10,1).astype('int').astype('str')
		nameArray  = numpy.core.defchararray.add(xi,index)
		if (dic is None):
			dic = {
				'xi0' : True,
				'xi1' : True,
				'xi2' : True,
				'xi3' : True,
				'xi4' : True,
				'xi5' : True,
				'xi6' : True,
				'xi7' : True,
				'xi8' : True,
				'xi9' : True,
				'plot' : False
			}
			if (self._correlation=='q_q' or self._correlation=='f_f'):
				for i in range(0,10):
					if (i%2==1): dic[ nameArray[i] ] = False
			for i in range(self._nb_multipole_max,10):
				dic[ nameArray[i] ] = False
		
		### Get the data
		xxx = xiMu[:,:,0]
		muu = xiMu[:,:,1]
		yyy = xiMu[:,:,2]
		yer = xiMu[:,:,3]

		### Keep the results
		result_xi = numpy.zeros(shape=(self._nbBin1D,self._nb_multipole_max,3))
		meanXXX   = numpy.mean(xxx,axis=1)
		for i in numpy.arange(0,self._nb_multipole_max):
			result_xi[:,i,0] = meanXXX
	
		for i in numpy.arange(self._nbBin1D):
	
			cut         = (muu[i,:]!=0.)
			tmpyyy      = yyy[i,:][cut]
			tmpyer      = yer[i,:][cut]
			xxxMu       = muu[i,:][cut]

			legendre = [ scipy.special.eval_legendre(j,xxxMu) for j in range(0,10) ]
			
			def chi2(xi0,xi1,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9):

				fit = 0.
				if (dic[ nameArray[0] ]): fit += xi0*legendre[0]
				if (dic[ nameArray[1] ]): fit += xi1*legendre[1]
				if (dic[ nameArray[2] ]): fit += xi2*legendre[2]
				if (dic[ nameArray[3] ]): fit += xi3*legendre[3]
				if (dic[ nameArray[4] ]): fit += xi4*legendre[4]
				if (dic[ nameArray[5] ]): fit += xi5*legendre[5]
				if (dic[ nameArray[6] ]): fit += xi6*legendre[6]
				if (dic[ nameArray[7] ]): fit += xi7*legendre[7]
				if (dic[ nameArray[8] ]): fit += xi8*legendre[8]
				if (dic[ nameArray[9] ]): fit += xi9*legendre[9]
				
				return numpy.sum( numpy.power( (tmpyyy-fit)/tmpyer ,2.) )

			### Init ad perform the fit
			m = Minuit(chi2, print_level=-1, pedantic=False)
			m.migrad()
	
			### Keep the results
			for j in numpy.arange(0,self._nb_multipole_max):
				result_xi[i,j,1] = m.values[ nameArray[j] ]
				result_xi[i,j,2] = m.errors[ nameArray[j] ]

			if (dic['plot']):
				plt.errorbar(xxxMu,tmpyyy,yerr=tmpyer, markersize=10,linewidth=2, marker='o',alpha=0.6)
				fit  = 0.
				for j in numpy.arange(0,self._nb_multipole_max):
					fit += m.values[ nameArray[j] ]*scipy.special.eval_legendre(j,xxxMu)
					print m.values[ nameArray[j] ]
				print
				plt.errorbar(xxxMu,fit)
				plt.xlabel(r'$\mu$', fontsize=40)
				plt.ylabel(r'$\xi(\mu,|s|)$', fontsize=40)
				myTools.deal_with_plot(False,False,False)
				plt.show()

				### Residuals
				plt.errorbar( xxxMu, (tmpyyy-fit)/tmpyer, markersize=10,linewidth=2, marker='o',alpha=0.6)
				plt.xlabel(r'$\mu$', fontsize=40)
				plt.ylabel(r'$Residulas \, \xi(\mu,|s|)$', fontsize=40)
				myTools.deal_with_plot(False,False,False)
				plt.show()

		cut = (result_xi[:,:,1]==0.)
		result_xi[cut] = 0.

		return result_xi
	def get_multipol_index(self, index):

		if (self._xiMul is None): self._xiMul = self.get_multipol(self._xiMu)

		xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		xi1D[:,:] = self._xiMul[:,index,:]

		return xi1D
	def get_CAMB(self,dic=None,dim='2D',distortion=False):

		### Get Camb
		if (dic is not None): dic['z'] = 0.
		camb = CAMB.CAMB(dic)._xi0

		if (dim=='2D'):
			xi = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
			xi[:,:,0] = self._xi2D[:,:,0]
			xi[:,:,1] = numpy.interp(xi[:,:,0],camb[:,0],camb[:,1])
			xi[:,:,2] = 0.00000001
			cut = self._xi2D[:,:,2]<=0.
			xi[:,:,1][cut] = 0.
			xi[:,:,2][cut] = 0.
			if (self._correlation=='q_f'): xi[:,:,1] *= -1.
			f = cp.fgrowth(self._xi2D_grid[:,:,2],dic['omega_matter_0'])
			xi[:,:,1] *= f*f
		elif (dim=='Mu'):
			xi = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4))
			xi[:,:,0] = self._xiMu[:,:,0]
			xi[:,:,1] = self._xiMu[:,:,1]
			xi[:,:,2] = numpy.interp(xi[:,:,0],camb[:,0],camb[:,1])
			xi[:,:,3] = 0.00000001
			cut = (self._xiMu[:,:,3]<=0.)
			xi[:,:,2][cut] = 0.
			xi[:,:,3][cut] = 0.
			if (self._correlation=='q_f'): xi[:,:,2] *= -1.
			f = cp.fgrowth(self._xiMu_grid[:,:,2],dic['omega_matter_0'])
			xi[:,:,2] *= f*f
		elif (dim=='We'):
			xi = numpy.zeros(shape=(self._nbBin1D,self._nb_wedges,3))
			xi[:,:,0] = self._xiWe[:,:,0]
			xi[:,:,1] = numpy.interp(xi[:,:,0],camb[:,0],camb[:,1])
			xi[:,:,2] = 0.00000001
			cut = (self._xiWe[:,:,2]<=0.)
			xi[:,:,1][cut] = 0.
			xi[:,:,2][cut] = 0.
			if (self._correlation=='q_f'): xi[:,:,1] *= -1.
			f = cp.fgrowth(self._xiWe_grid[:,:,1],dic['omega_matter_0'])
			xi[:,:,1] *= f*f
		elif (dim=='1D'):
			xi = numpy.zeros(shape=(self._nbBin1D,3))
			xi[:,0] = self._xi1D[:,0]
			xi[:,1] = numpy.interp(xi[:,0],camb[:,0],camb[:,1])
			xi[:,2] = 0.00000001
			cut = (self._xi1D[:,2]<=0.)
			xi[:,1][cut] = 0.
			xi[:,2][cut] = 0.
			if (self._correlation=='q_f'): xi[:,1] *= -1.
			f = cp.fgrowth(self._xi1D_grid[:,1],dic['omega_matter_0'])
			xi[:,1] *= f*f

		if (distortion and (dim=='2D' or dim=='1D') ):
			path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_'+dim+'_'+ self._middlefix + '.txt'
			matrix = numpy.loadtxt(path)

			if (dim=='2D'):
				tmp = myTools.convert2DTo1D(xi[:,:,1], self._nbBinX2D, self._nbBinY2D)
				tmp = numpy.dot(matrix,tmp)
				xi[:,:,1] = myTools.convert1DTo2D(tmp, self._nbBinX2D, self._nbBinY2D)
				cut = (xi[:,:,1]==0.)
				xi[:,:,2][cut] = 0.
			elif (dim=='1D'):
				xi[:,1] = numpy.dot(matrix,xi[:,1])
				cut = (xi[:,1]==0.)
				xi[:,2][cut] = 0.

		return xi
	def fit_CAMB(self,xi1D=None,dic=None,distortion=True):
		
		### Constants
		if (xi1D is None): xi1D = self._xi1D
		if (dic is None):
			dic = raw_dic_CAMB
		
		### Get the data
		cut = numpy.logical_and( (xi1D[:,0]>=dic['start_fit']),(xi1D[:,0]<dic['end_fit']) )
		xxx = xi1D[:,0][cut]
		yyy = xi1D[:,1][cut]
		yer = xi1D[:,2][cut]

		### Get Camb
		camb = dic['CAMB']
		if (dic['mulpol_index']==0):
			yyy_Camb = numpy.interp(xi1D[:,0],camb._xi0[:,0],camb._xi0[:,1])
		elif (dic['mulpol_index']==2):
			yyy_Camb = numpy.interp(xi1D[:,0],camb._xi2[:,0],camb._xi2[:,1])
		elif (dic['mulpol_index']==4):
			yyy_Camb = numpy.interp(xi1D[:,0],camb._xi4[:,0],camb._xi4[:,1])

		if (distortion):
			path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_1D_'+ self._middlefix + '.txt'
			matrix = numpy.loadtxt(path)
			yyy_Camb = numpy.dot(matrix,yyy_Camb)

		yyy_Camb = yyy_Camb[cut]
	
		if (dic['guess_b']):
			b_init = numpy.mean(yyy[ (xxx<dic['min_for_guess']) ]/yyy_Camb[ (xxx<dic['max_for_guess']) ])
			print '  b_init = ', b_init
		else:
			b_init = dic['b']
	
		### Define the chi^{2} for the fit
		def chi2(b,roof):
			fit = yyy_Camb*b+roof
			return numpy.sum( numpy.power( (yyy-fit)/yer ,2.) )
	
		### Init and perform the fit
		
		m = Minuit(chi2, b=b_init,error_b=0.1, roof=dic['roof'],error_roof=0.1,print_level=-1, errordef=0.01,fix_roof=dic['fix_roof_nul'] ) 	
		m.migrad()
	
		### Get result
		dic['b']       = m.values[ 'b' ]
		dic['roof']    = m.values[ 'roof' ]
		dic['chi^{2}'] = numpy.sum( numpy.power( (yyy - (yyy_Camb*dic['b']+dic['roof']) )/yer ,2.) )
		print '  b = ', dic['b']
		print '  roof = ', dic['roof']
	
		### Print chi^2
		if (dic['fix_roof_nul']):
			print '  DoF   = ', yyy.size, ' - 1'
		else:
			print '  DoF   = ', yyy.size, ' - 2'
		print '  chi^{2} = ', dic['chi^{2}']
		
		return dic
	def fit_CAMB_2d(self, realisation_type, correlation_matrix_path=None):

		### Data
		data = self._xi2D[:,:,1].flatten()
		### Cov
		path_to_cov = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_2D.npy'
		if (correlation_matrix_path is None): correlation_matrix_path = path_to_cov
		print '  correlation matrix path = ', correlation_matrix_path
		cov = numpy.load(path_to_cov)
		cor = myTools.getCorrelationMatrix( numpy.load(correlation_matrix_path) )
		cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))
		icov = inverse = numpy.linalg.inv(cov)
		### LYA
		xiMu, xiWe, xi1D, xi2D = self.read_metal_model('LYA')
		xi0_LYA = xi2D[:,:,1,0].flatten()
		xi2_LYA = xi2D[:,:,1,1].flatten()
		xi4_LYA = xi2D[:,:,1,2].flatten()
		### Si-II (1190)
		xiMu, xiWe, xi1D, xi2D = self.read_metal_model('SiII(1190)')
		xi0_Si2a = xi2D[:,:,1,0].flatten()
		xi2_Si2a = xi2D[:,:,1,1].flatten()
		xi4_Si2a = xi2D[:,:,1,2].flatten()
		### Si-II (1193)
		xiMu, xiWe, xi1D, xi2D = self.read_metal_model('SiII(1193)')
		xi0_Si2b = xi2D[:,:,1,0].flatten()
		xi2_Si2b = xi2D[:,:,1,1].flatten()
		xi4_Si2b = xi2D[:,:,1,2].flatten()
		### Si-III (1207)
		xiMu, xiWe, xi1D, xi2D = self.read_metal_model('SiIII(1207)')
		xi0_Si3 = xi2D[:,:,1,0].flatten()
		xi2_Si3 = xi2D[:,:,1,1].flatten()
		xi4_Si3 = xi2D[:,:,1,2].flatten()

		### Chi^{2}
		def chi2(a0,a1,a2,a3,a4,a5):
			model = a0*xi0_LYA + a1*xi2_LYA + a2*xi4_LYA + a3*xi0_Si3 + a4*xi0_Si2a + a5*xi0_Si2b
			return numpy.dot( numpy.dot( (data - model).T,icov),(data - model) )

		m = Minuit(chi2,a0=-1.,error_a0=0.1,a1=1.,error_a1=0.1,a2=0.,error_a2=0.1,a3=0.,error_a3=0.1,a4=0.,error_a4=0.1,a5=0.,error_a5=0.1,print_level=-1, errordef=0.1, fix_a2=True) 	

		print 'Starting fit'
		m.migrad()

		a0 = m.values[ 'a0' ]
		a1 = m.values[ 'a1' ]
		a2 = m.values[ 'a2' ]
		a3 = m.values[ 'a3' ]
		a4 = m.values[ 'a4' ]
		a5 = m.values[ 'a5' ]
		print "  a0 = ", a0
		print "  a1 = ", a1
		print "  a2 = ", a2
		print "  a3 = ", a3
		print "  a4 = ", a4
		print "  a5 = ", a5


		model    = a0*xi0_LYA + a1*xi2_LYA + a2*xi4_LYA + a3*xi0_Si3 + a4*xi0_Si2a + a5*xi0_Si2b
		print '  chi^{2}  = ', numpy.dot( numpy.dot( (data - model).T,icov),(data - model) )
		print '  nb bin   = ', data.size
		print '  nb param = ', 6

		xxx      = self._xi2D[:,:,0]
		var      = self._xi2D[:,:,2].flatten()
		residual = myTools.convert1DTo2D( (data-model)/var,self._nbBinX2D,self._nbBinY2D)
		model    = myTools.convert1DTo2D(model,self._nbBinX2D,self._nbBinY2D)
		###
		myTools.plot2D(residual)
		myTools.plot2D(model)
		###
		myTools.plot2D(residual*xxx)
		myTools.plot2D(model*xxx)
		###
		myTools.plot2D(residual*xxx*xxx)
		myTools.plot2D(model*xxx*xxx)

		return
	def write_metal_model(self, with_metals_templates=False):

		nb_multipol_calculated_templates = 3

		if (self._correlation=='q_f'):
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1
			name_metal_1 = ['']
			name_metal_2 = ['SiII(1260)','SiIII(1207)','SiII(1193)','SiII(1190)']
			name_metal_to_save_1 = ['QSO_']
			name_metal_to_save_2 = ['Si2c','Si3','Si2b','Si2a']
		elif (self._correlation=='f_f'):
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1
			name_metal_1 = ['LYA_','SiII(1260)_','SiIII(1207)_','SiII(1193)_','SiII(1190)_']
			name_metal_2 = ['SiII(1260)','SiIII(1207)','SiII(1193)','SiII(1190)']
			name_metal_to_save_1 = ['Lya_','Si2c_','Si3_','Si2b_','Si2a_']
			name_metal_to_save_2 = ['Si2c','Si3','Si2b','Si2a']

		if (with_metals_templates):
			path_to_BAOFIT += '__withMetalsTemplates'
		#path_to_BAOFIT += '/bao2D_'
		path_to_BAOFIT += '/metTemp_'
		suffix = ['.0.dat','.2.dat','.4.dat']




		for i in numpy.arange(len(name_metal_1)):

			for j in numpy.arange(len(name_metal_2)):
				xiMu, xiWe, xi1D, xi2D = self.read_metal_model( name_metal_1[i]+name_metal_2[j] )

				#myTools.plot_2d_correlation(xi2D[:,:,:,0],0)

				for k in range(0,nb_multipol_calculated_templates):

					correlation = numpy.zeros( shape=(self._nbBin2D,2) )
					indexMatrix = numpy.arange(self._nbBin2D)
					correlation[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
					correlation[:,1] = xi2D[:,:,1,k].flatten()
					print path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k]
					numpy.savetxt( path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k],zip(correlation[:,0],correlation[:,1]),fmt='%u %1.20e')
					#numpy.savetxt( "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results/Test_fitter_Nicolas/bao2D_" + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k],zip(correlation[:,0],correlation[:,1]),fmt='%u %1.20e')
					"""
					cutCorrelation = (correlation[:,1]!=0.)
					if (correlation[:,0][cutCorrelation].size==0): continue

					### Get the array and sort it
					correlation2 = numpy.zeros( shape=(self._nbBin2D,4) )
					for l in numpy.arange(self._nbBin2D):
						idx = correlation[l,0].astype(int)
						correlation2[idx,0] = idx
						correlation2[idx,1] = correlation[l,1]
						correlation2[idx,2] = self._minX2D+0.5*self._binSize + (idx%self._nbBinX2D)*self._binSize
						correlation2[idx,3] = self._minY2D+0.5*self._binSize + (idx/self._nbBinX2D)*self._binSize

					print path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k]
					numpy.savetxt( path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k],zip(correlation2[:,2],correlation2[:,3],correlation2[:,1]),fmt='%u %u %1.20e')
					"""
		return
	def write_metal_model_mean_simulation(self, dic_class, dic_simu):

		nb_multipol_calculated_templates = 3

		if (self._correlation=='q_f'):
			path_to_BAOFIT = dic_simu['path_to_simu']+'Results'+dic_simu['prefix']+'/BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1
			name_metal_1 = ['']
			name_metal_2 = ['SiII(1260)','SiIII(1207)','SiII(1193)','SiII(1190)']
			name_metal_to_save_1 = ['QSO_']
			name_metal_to_save_2 = ['Si2c','Si3','Si2b','Si2a']
		elif (self._correlation=='f_f'):
			path_to_BAOFIT = dic_simu['path_to_simu']+'Results'+dic_simu['prefix']+'/BaoFit_'+self._correlation+'__'+self._f1
			name_metal_1 = ['LYA_','SiII(1260)_','SiIII(1207)_','SiII(1193)_','SiII(1190)_']
			name_metal_2 = ['SiII(1260)','SiIII(1207)','SiII(1193)','SiII(1190)']
			name_metal_to_save_1 = ['Lya_','Si2c_','Si3_','Si2b_','Si2a_']
			name_metal_to_save_2 = ['Si2c','Si3','Si2b','Si2a']

		print path_to_BAOFIT
		path_to_BAOFIT += '/bao2D_'
		suffix = ['.0.dat','.2.dat','.4.dat']


		for i in numpy.arange(len(name_metal_1)):
			for j in numpy.arange(len(name_metal_2)):

				xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3,nb_multipol_calculated_templates))

				for k in numpy.arange(dic_simu['nb_box']):
					for l in numpy.arange(dic_simu['nb_simu']):

						raw = dic_simu['path_to_simu'] + 'Box_00' + str(k) + '/Simu_00' + str(l) +'/Results'
						raw += dic_simu['prefix']
						raw += '/'
						dic_class['path_to_txt_file_folder'] = raw
						corr = Correlation3D(dic_class)

						tmp_xiMu, tmp_xiWe, tmp_xi1D, tmp_xi2D = corr.read_metal_model( name_metal_1[i]+name_metal_2[j] )
						xi2D[:,:,:,:,] += tmp_xi2D

				xi2D[:,:,:,:,] /= dic_simu['nb_box']*dic_simu['nb_simu']

				for k in numpy.arange(0,nb_multipol_calculated_templates):

					correlation = numpy.zeros( shape=(self._nbBin2D,4) )
					indexMatrix = numpy.arange(self._nbBin2D)
					correlation[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
					correlation[:,1] = xi2D[:,:,1,k].flatten()
					cutCorrelation = (correlation[:,1]!=0.)
					if (correlation[:,0][cutCorrelation].size==0): continue
					print path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k]
					numpy.savetxt( path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k],zip(correlation[:,0],correlation[:,1]),fmt='%u %1.20e')
					"""
					### Get the array and sort it
					correlation2 = numpy.zeros( shape=(self._nbBin2D,4) )
					for l in numpy.arange(self._nbBin2D):
						idx = correlation[l,0].astype(int)
						correlation2[idx,0] = idx
						correlation2[idx,1] = correlation[l,1]
						correlation2[idx,2] = self._minX2D+0.5*self._binSize + (idx%self._nbBinX2D)*self._binSize
						correlation2[idx,3] = self._minY2D+0.5*self._binSize + (idx/self._nbBinX2D)*self._binSize

					print path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k]
					numpy.savetxt( path_to_BAOFIT + name_metal_to_save_1[i]+name_metal_to_save_2[j] + suffix[k],zip(correlation2[:,2],correlation2[:,3],correlation2[:,1]),fmt='%u %u %1.20e')
					"""
		return
	def write_pyLyA_config_file(self, with_metals_templates=False, prefix2=''):

		prefix = self._path_to_txt_file_folder + 'pyLyA'
		if (with_metals_templates):
			prefix += '_withMetalsTemplates'
		prefix += prefix2
		prefix += '/'
		print '  SAVING INTO: ', prefix
		subprocess.call('mkdir ' + prefix, shell=True)

		data         = ''
		data_autoQSO = 'None'
		data_cross   = 'None'
		data_auto    = 'None'
		fix          = 'None'
		free         = 'None'
		metals       = 'None'
		if (self._correlation=='q_q'):
			data = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._q1 +'/bao2D'
			data_autoQSO = data+'-exp.fits'
			prefix += 'autoQSO_alone.'
		elif (self._correlation=='q_f'):
			data= self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1 +'/bao2D'
			data_cross = data+'-exp.fits'
			fix = "growth_factor_qso Lpar_cross bias_lya*(1+beta_lya) qso_metal_boost"
			free = "bias_qso drp"
			if (with_metals_templates): metals = "Si3 Si2a Si2b Si2c"
			prefix += 'cross_alone.'
		elif (self._correlation=='f_f'):
			data = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1 +'/bao2D'
			data_auto = data+'-exp.fits'
			fix = "alpha SigmaNL_perp 1+f Lpar_auto"
			if (with_metals_templates): metals = "Si3 Si2a Si2b Si2c"
			free = "Lpar_auto"
			prefix += 'auto_alone.'

		string_ini = """
[PARAMETER]

### String
model         = /home/gpfs/manip/mnt0607/bao/hdumasde/Program/pyLyA/models/DR9LyaMocks/DR9LyaMocks.fits
metal_prefix  = """ +data+"""
data_auto     = """ +data_auto+"""
data_cross    = """ +data_cross+"""
output_prefix = """ +prefix+"""

### String list
metals        = """ +metals+"""
fix           = """ +fix+"""
free          = """ +free+"""

### General
bias_lya*(1+beta_lya) = -0.39536
beta_lya          = 1.352
ap                = 1.
at                = 1.
rmin              = 10.
rmax              = 180.
ell_max           = 6

### Auto-correlation
alpha             = 3.8
Lpar_auto         = 8.
SigmaNL_perp      = 0.
1+f               = 1.

## Cross-correlation
bias_qso          = 3.125
growth_factor_qso = 0.962524
drp               = 0.
Lpar_cross        = 9.

### LLS
bias_lls          = -0.01
beta_lls          = 0.6
L0_lls            = 20

### UV
bias_gamma        = 0.13
bias_prim         = -0.6666666666666666
kappa0            = 300.

### Metals
bias_Si2a         = -0.01
beta_Si2a         = 1.4
bias_Si2b         = -0.01
beta_Si2b         = 1.4
bias_Si2c         = -0.01
beta_Si2c         = 1.4
bias_Si3          = -0.01
beta_Si3          = 1.4
qso_metal_boost   = 1.
"""

		text_file = open(prefix+'ini', "w")
		text_file.write(string_ini)
		text_file.close()

		return data, prefix
	def send_pyLyA(self, with_metals_templates=False, prefix2=''):

		### Create config file
		data, config_file = self.write_pyLyA_config_file(with_metals_templates=with_metals_templates, prefix2=prefix2)

		### Create the fits file
		'''
		command = 'rm ' + data + '-exp.fits'
		print command
		subprocess.call(command, shell=True)
		command = '/home/gpfs/manip/mnt0607/bao/hdumasde/Program/pyLyA/bin/export2fits ' + data
		print command
		subprocess.call(command, shell=True)
		'''

		### Send the fit
		command = '/home/gpfs/manip/mnt0607/bao/hdumasde/Program/pyLyA/bin/fit -c ' + config_file+'ini --verbose --ap 1. --at 1.'
		#command = '/home/gpfs/manip/mnt0607/bao/hdumasde/Program/pyLyA/bin/fit -c ' + config_file+'fit.config --verbose --ap 1. --at 1.'
		print command
		subprocess.call(command, shell=True)

		return
	def write_BAOFIT_ini_and_grid(self, param, dist_matrix=False, with_metals_templates=False, prefix2=''):

		precision = 1000.
		param = param.astype('str')
		### Is there a distortion matrix
		str_dist_matrix = '#'
		if (dist_matrix):
			str_dist_matrix = ''
		### Are there metals
		str_metals = '#'
		if (with_metals_templates):
			str_metals = ''
		### Is it a q_f correlation
		str_corr  = '#'
		str_corr2 = 'value'
		str_corr3 = ''
		if (self._correlation=='q_f'):
			str_corr  = ''
			str_corr2 = 'fix'
			str_corr3 = '#'

		if (self._correlation=='q_q'):
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._q1
		elif (self._correlation=='q_f'):
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1
		elif (self._correlation=='f_f'):
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1

		path_to_data = path_to_BAOFIT + '/bao2D'
		if (with_metals_templates):
			path_to_BAOFIT += '__withMetalsTemplates'
		if (not dist_matrix):
			path_to_BAOFIT += '__noDistortionMatrix'
		path_to_BAOFIT += prefix2
		path_to_BAOFIT += '/bao2D'

		string_ini = """

## Linear theory P(k) templates with and w/o wiggles
modelroot = """ + const.path_to_BAOFIT_model__ + """
#fiducial =  PlanckDR12
#nowiggles = PlanckDR12SB
#omega-matter    = 0.3145695148634867
#hubble-constant = 0.6731
#sigma8 = 0.829
#zref = 2.3

fiducial =  DR9LyaMocks
nowiggles = DR9LyaMocksSB

omega-matter    = 0.27
hubble-constant = 0.7
sigma8 = 0.829
zref = 2.25

## k-space fit
kspace = true
ell-max = 4

# Model configuration
"""+str_corr+"""cross-correlation = yes
anisotropic = yes
decoupled   = yes
custom-grid = yes
pixelize    = yes
#pixelize-alt = yes
"""+str_dist_matrix+"""dist-matrix = yes
dist-matrix-order = """+str(self._nbBin2D)+"""

# Parameter setup
model-config = value[beta]=                      """+param[0] +""";
#model-config = """+str_corr2+"""[(1+beta)*bias]= """+param[1] +""";
model-config = value[(1+beta)*bias]= """+param[1] +""";
model-config = fix[gamma-bias]=                  """+param[2] +""";
model-config = fix[gamma-beta]=                  """+param[3] +""";
"""+str_corr+"""model-config = value[delta-v]=   """+param[4] +""";
"""+str_corr+"""model-config = value[bias2]=     """+param[5] +""";
"""+str_corr+"""model-config = fix[beta2*bias2]= """+param[6] +""";
model-config = fix[1+f]=                         """+param[7] +""";
model-config = fix[SigmaNL-perp]=                """+param[8] +""";
model-config = fix[BAO amplitude]=               """+param[9] +""";
model-config = fix[BAO alpha-iso]=               """+param[10]+""";
model-config = value[BAO alpha-parallel]=        """+param[11]+""";
model-config = value[BAO alpha-perp]=            """+param[12]+""";
model-config = fix[gamma-scale]=                 """+param[13]+""";
#model-config = fix[pixel scale]=0.7;
#model-config = fix[pixel scale2]=0.7;
model-config = fix[pixel scale]=3.15;

## Metal correlations
"""+str_metals+"""metal-model-interpolate = true
"""+str_metals+"""metal-model-name = """+path_to_data+"""
"""+str_metals+"""model-config = value[beta Si2a] = """+param[15]+""";
"""+str_metals+"""model-config = value[bias Si2a] = """+param[16]+""";
"""+str_metals+"""model-config = value[beta Si2b] = """+param[17]+""";
"""+str_metals+"""model-config = value[bias Si2b] = """+param[18]+""";
"""+str_metals+"""model-config = value[beta Si2c] = """+param[19]+""";
"""+str_metals+"""model-config = value[bias Si2c] = """+param[20]+""";
"""+str_metals+"""model-config = fix[beta Si3]  = """+param[21]+""";
"""+str_metals+"""model-config = fix[bias Si3]  = """+param[22]+""";
"""+str_metals+"""model-config = gaussprior[beta Si2a] @ (0,2.8);
"""+str_metals+"""model-config = gaussprior[beta Si2b] @ (0,2.8);
"""+str_metals+"""model-config = gaussprior[beta Si2c] @ (0,2.8);
"""+str_metals+"""model-config = gaussprior[beta Si3]  @ (0,2.8);


## 2D chisq scan in BAO parameters
#model-config = binning[BAO alpha-parallel] ={0.7:1.4}*50
#model-config = binning[BAO alpha-perp]     ={0.7:1.4}*50	
model-config = binning[BAO alpha-parallel] ={0.98:1.02}*50
model-config = binning[BAO alpha-perp]     ={0.98:1.02}*50

## Maximum allowed radial dilation (increases the range that model needs to cover)
dilmin = 0.01
dilmax = 1000.

## Non-linear broadening with 1+f = (SigmaNL-par)/(SigmaNL-perp)

# boxprior keeps result positive (since model only depends on squared value)
model-config = boxprior[SigmaNL-perp] @ (0,6);
# un-comment next line to broaden all scales (default is peak only)
#nl-broadband = true

# Broadband distortion model
#"""+str_corr+"""dist-add = rP,rT=0:2,-3:1
#dist-add = rP,rT=0:2,-3:1
#"""+str_corr3+"""dist-add = -2:0,0:4:2,0

### Data Options #############################################

## Data to analyze
data = """+path_to_data+"""
dist-matrix-name = """+path_to_data+"""

## Data format
data-format = comoving-cartesian
axis1-bins = ["""+ str(self._minY2D) + ':' + str(self._maxY2D) + """]*"""+ str(self._nbBinY2D) +"""
axis2-bins = ["""+ str(self._minX2D) + ':' + str(self._maxX2D) + """]*"""+ str(self._nbBinX2D) +"""
axis3-bins = {"""+ str(self._meanZ) +"""}

### Analysis Options #########################################

# Cuts to apply before fitting
rmin = 40
rmax = 180
#rperp-min = 10

# Generate a second set of outputs with the additive distortion turned off
#alt-config = fix[dist*]=0

# Do not dump multipoles (since the distortion model multipole integrals are singular)
ndump = 0

# Prefix to use for all analysis output files
output-prefix = """ + path_to_BAOFIT + """.
"""

		text_file = open(path_to_BAOFIT+'.ini', "w")
		text_file.write(string_ini)
		text_file.close()
		

		return
	def send_BAOFIT(self, dic_send_BAOFIT=None):
		"""

			Send the fit of the BAO with BAOFIT
			realisation_type:
			cor:=None

		"""

		if (dic_send_BAOFIT is None):
			dic_send_BAOFIT = {
				'realisation_type'        : None,
				'correlation_matrix_path' : None,
				'saving_data' : True,
				'scan' : False,
				'toyMC' : 0,
				'dist_matrix' : False,
				'only_diagonal' : False,
				'with_metals_templates'  : False,
				'get_param_from_file'    : False,
				'dic_simu'               : None,
				'prefix'                 : ''
			}
		realisation_type        = dic_send_BAOFIT['realisation_type']
		correlation_matrix_path = dic_send_BAOFIT['correlation_matrix_path']
		saving_data             = dic_send_BAOFIT['saving_data']
		scan                    = dic_send_BAOFIT['scan']
		toyMC                   = dic_send_BAOFIT['toyMC']
		dist_matrix             = dic_send_BAOFIT['dist_matrix']
		only_diagonal           = dic_send_BAOFIT['only_diagonal']
		with_metals_templates   = dic_send_BAOFIT['with_metals_templates']
		get_param_from_file     = dic_send_BAOFIT['get_param_from_file']
		dic_simu                = dic_send_BAOFIT['dic_simu']
		fit_prefix              = dic_send_BAOFIT['prefix']

		pathToSave = dic_simu['path_to_simu'] + 'Results'
		pathToSave += dic_simu['prefix']
		pathToSave += '/'
		pathToSave += self._prefix + '_' + self._middlefix + '_result_'

		if ( (dic_simu is not None) and realisation_type=='Simu_stack' ):
			path_to_distortion_matrix = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix'] + '/' + self._prefix + '_'+ self._middlefix + '_result_dist_2D.npy'
		else:
			path_to_distortion_matrix = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_2D_'+ self._middlefix + '.txt'
		print '   Path to dmat = ', path_to_distortion_matrix
		
		if (dic_simu is not None):
			raw_path_to_cov = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix'] + '/'
			if (realisation_type=='Simu_stack'):
				self._path_to_txt_file_folder = raw_path_to_cov
		if (self._correlation=='q_q'):
			param = numpy.asarray( [0.268344,0.962524,-2.,0.,0.,3.6,0.,0.,1.966,1.,1.,1.,1.,0.,0.,1.4,-0.01,1.4,-0.01,1.4,-0.01,1.4,-0.01] )
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._q1
		elif (self._correlation=='q_f'):
			param = numpy.asarray( [1.1,-0.351,0.9,0.,0.,3.64,0.962524,1.966,3.26,1.,1.,1.,1.,0.,0.,1.4,-0.01,1.4,-0.01,1.4,-0.01,1.4,-0.01] )
			if (dic_simu is not None):
				param[0] = 1.329723
				param[1] = -0.400
				param[5] = 3.554
				param[8] = 0.
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1
		elif (self._correlation=='f_f'):
			param = numpy.asarray( [1.2,-0.351,3.8,0.,0.,3.6,0.962524,1.966,3.26,1.,1.,1.,1.,0.,0.,1.4,-0.01,1.4,-0.01,1.4,-0.01,1.4,0.] )
			if (dic_simu is not None):
				param[0] = 1.329723
				param[1] = -0.400
				param[8] = 0.
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1

		path_to_save = path_to_BAOFIT + '/bao2D.'
		print '  saving data in ', path_to_save
		if (with_metals_templates):
			path_to_BAOFIT += '__withMetalsTemplates'
		if (not dist_matrix):
			path_to_BAOFIT += '__noDistortionMatrix'
		prefix2 = fit_prefix
		if (dic_simu is not None): prefix2 += dic_simu['prefix2']
		path_to_BAOFIT += prefix2 + '/'
		print '  path_to_BAOFIT = ', path_to_BAOFIT


		try:
			subprocess.call('mkdir ' + path_to_BAOFIT, shell=True)
		except:
			print '  Folder exists'
		path_to_BAOFIT += 'bao2D.'

		### Get parameter from file
		if (get_param_from_file):
			data  = numpy.loadtxt(path_to_BAOFIT+'save.pars')
			param = data[:,1]

		self.write_BAOFIT_ini_and_grid(param,dist_matrix, with_metals_templates, prefix2)

		if (saving_data):
			
			### Save .grid
			print '  Saving .grid'
			grid = numpy.zeros( shape=(self._nbBin2D,4) )
			indexMatrix = numpy.arange(self._nbBin2D)
			grid[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
			grid[:,1] = self._xi2D_grid[:,:,1].flatten()
			grid[:,2] = self._xi2D_grid[:,:,0].flatten()
			grid[:,3] = self._xi2D_grid[:,:,2].flatten()
			numpy.savetxt(path_to_save+'grid',zip(grid[:,0],grid[:,1],grid[:,2],grid[:,3]),fmt='%u %1.20e %1.20e %1.20e')
			del grid

			### Save .data
			print '  Saving .data'
			correlation = numpy.zeros( shape=(self._nbBin2D,2) )
			indexMatrix = numpy.arange(self._nbBin2D)
			correlation[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
			correlation[:,1] = self._xi2D[:,:,1].flatten()
			cutCorrelation = (correlation[:,1]!=0.)
			numpy.savetxt( path_to_save + 'data',zip(correlation[:,0][cutCorrelation],correlation[:,1][cutCorrelation]),fmt='%u %1.20e')
			del correlation, indexMatrix, cutCorrelation
			
			### Save .cov
			if (self._correlation=='q_q'):
				if (dic_simu is not None):
					path_to_cov       = raw_path_to_cov + 'xi_QSO_QSO_result_cov_2D.npy'
					nb_of_realisation = numpy.load(raw_path_to_cov + 'xi_QSO_QSO_result_list_1D.npy')[0,:].size
			elif (self._correlation=='q_f'):
				if (dic_simu is not None):
					path_to_cov       = raw_path_to_cov + 'xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
					nb_of_realisation = numpy.load(raw_path_to_cov + 'xi_delta_QSO_LYA_QSO_result_list_1D.npy')[0,:].size
				else: path_to_cov = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_2D.npy'
			elif (self._correlation=='f_f'):
				if (dic_simu is not None):
					path_to_cov       = raw_path_to_cov + 'xi_A_delta_delta_LYA_result_cov_2D.npy'
					nb_of_realisation = numpy.load(raw_path_to_cov + 'xi_A_delta_delta_LYA_result_list_1D.npy')[0,:].size
				else: path_to_cov = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_2D.npy'

			if (correlation_matrix_path is None): correlation_matrix_path = path_to_cov
			print '  correlation matrix path = ', correlation_matrix_path

			print '  covariance matrix path = ', path_to_cov
			cov = numpy.load(path_to_cov)
			cor = myTools.getCorrelationMatrix( numpy.load(correlation_matrix_path) )

			print '  Saving .cov'
			cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))
			if ( (dic_simu is not None) and realisation_type=='Simu_stack' ):
				print '  nb of realisation = ', nb_of_realisation
				cov /= nb_of_realisation
			if (only_diagonal):
				'''
				mean_cor_off_diag = self._nbBin2D/(self._nbBin2D-1.)*numpy.mean( cor-numpy.diag( numpy.diag(cor)) )
				print '  The mean off diagonal term is = ', mean_cor_off_diag
				mean_cor_matrix = numpy.ones( shape=(self._nbBin2D,self._nbBin2D) )*mean_cor_off_diag
				mean_cor_matrix -= numpy.diag(numpy.diag(mean_cor_matrix))
				cor = numpy.diag( numpy.diag(cor) ) + mean_cor_matrix
				cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))
				'''
				cov = numpy.diag( numpy.diag(cov) )

			#print cov[0,0], self._xi2D[0,0,2]*self._xi2D[0,0,2]
			covarianceMatrix = numpy.zeros( shape=(self._nbBin2D*self._nbBin2D,3) )
			indexMatrix1 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)/self._nbBin2D
			indexMatrix2 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)%self._nbBin2D
			indexMatrix1 = (indexMatrix1%self._nbBinY2D)*self._nbBinX2D + indexMatrix1/self._nbBinY2D
			indexMatrix2 = (indexMatrix2%self._nbBinY2D)*self._nbBinX2D + indexMatrix2/self._nbBinY2D
			covarianceMatrix[:,0] = numpy.triu(indexMatrix1,k=0).flatten()
			covarianceMatrix[:,1] = numpy.triu(indexMatrix2,k=0).flatten()
			covarianceMatrix[:,2] = numpy.triu(cov,k=0).flatten()
			cutCovarMatrix = (covarianceMatrix[:,2]!=0.)

			numpy.savetxt( path_to_save + 'cov',zip(covarianceMatrix[:,0][cutCovarMatrix],covarianceMatrix[:,1][cutCovarMatrix],covarianceMatrix[:,2][cutCovarMatrix]),fmt='%u %u %1.20e')
			del cov, covarianceMatrix, cutCovarMatrix
			
			### Save .dmat
			if (dist_matrix):
				print '  Saving .dmat'
				if (dic_simu['prefix']=='_raw_from_JeanMarc' or dic_simu['prefix']=='_delta_gaussian' or dic_simu['prefix']=='_only_LR' or dic_simu['prefix']=='_only_LR_noRand' or dic_simu['prefix']=='_only_LR_noRand_noRSD'):
					dmatData = numpy.eye(self._nbBin2D)
				elif ( (dic_simu is not None) and realisation_type=='Simu_stack' ):
					print '  LOAD DMAT FROM: ', path_to_distortion_matrix
					dmatData = numpy.load(path_to_distortion_matrix)
				else:
					dmatData = numpy.loadtxt(path_to_distortion_matrix)
				indexMatrix1 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)/self._nbBin2D
				indexMatrix2 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)%self._nbBin2D
				indexMatrix1 = (indexMatrix1%self._nbBinY2D)*self._nbBinX2D + indexMatrix1/self._nbBinY2D
				indexMatrix2 = (indexMatrix2%self._nbBinY2D)*self._nbBinX2D + indexMatrix2/self._nbBinY2D
				distortionMatrix = numpy.zeros( shape=(self._nbBin2D*self._nbBin2D,3) )
				distortionMatrix[:,0] = indexMatrix1.flatten()
				distortionMatrix[:,1] = indexMatrix2.flatten()
				distortionMatrix[:,2] = dmatData.flatten()
				cutDistortionMatrix = (distortionMatrix[:,2]!=0.)
				numpy.savetxt( path_to_save + 'dmat',zip(distortionMatrix[:,0][cutDistortionMatrix], distortionMatrix[:,1][cutDistortionMatrix], distortionMatrix[:,2][cutDistortionMatrix]),fmt='%u %u %1.20e')
				del dmatData, distortionMatrix, cutDistortionMatrix
			

		print '\n\n'
		### Send the fit (#--parameter-scan'   ### --toymc-samples 10000)
		command = const.path_to_BAOFIT_bin__ + ' -i ' + path_to_BAOFIT + 'ini'
		if (scan):
			command += ' --parameter-scan'
		if (toyMC!=0):
			command += ' --toymc-samples ' + str(toyMC)
		print command
		subprocess.call(command, shell=True)

		print '\n\n'

		return
	def get_mapping_2D_to_1D(self):

		### Constants
		nb_random = 1000000

		### to return
		mapping_2D_to_1D = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,self._nbBin1D) )

		for i in numpy.arange(0.,self._nbBinX2D):
			for j in numpy.arange(0.,self._nbBinY2D):

				### to get a line probability
				rand = numpy.random.random(nb_random)
				r_perp = self._minX2D + self._binSize*numpy.sqrt(i*i +rand*( 2.*i+1. ) )

				### to get a flat probability
				rand = numpy.random.random(nb_random)
				r_paral = self._minY2D + self._binSize*(j + rand)

				### Find bin of the wedge
				r        = numpy.sqrt( r_perp*r_perp + r_paral*r_paral )
				r_index  = (r/self._binSize).astype('int')

				for k in numpy.arange(0,nb_random):
					if (r[k]>=self._max1D): continue
					mapping_2D_to_1D[i,j,r_index[k]] += 1.
		mapping_2D_to_1D /= nb_random

		return mapping_2D_to_1D
	def get_mapping_2D_to_we(self):

		### Constants
		nb_random = 1000000

		### Usefull values
		mu_bin_size = 1./self._nbBinM

		### to return
		mapping_2D_to_we = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,self._nbBin1D,self._nb_wedges) )

		for i in numpy.arange(0.,self._nbBinX2D):
			for j in numpy.arange(0.,self._nbBinY2D):

				### to get a line probability
				rand = numpy.random.random(nb_random)
				r_perp = self._minX2D + self._binSize*numpy.sqrt(i*i +rand*( 2.*i+1. ) )

				### to get a flat probability
				rand = numpy.random.random(nb_random)
				r_paral = self._minY2D + self._binSize*(j + rand)

				### Find bin of the wedge
				r        = numpy.sqrt( r_perp*r_perp + r_paral*r_paral )
				r_index  = (r/self._binSize).astype('int')
				mu       = numpy.abs(r_paral)/r
				we_index = self.from_mu_to_wedges_index(mu)

				for k in numpy.arange(0,nb_random):
					if (r[k]>=self._max1D): continue
					mapping_2D_to_we[i,j,r_index[k],we_index[k]] += 1.
		mapping_2D_to_we /= nb_random

		return mapping_2D_to_we
	def plot_1d(self, x_power=0, other=[], title=True, color=None):

		with_camb = False
		list_corr_to_plot = numpy.append( [self],other )
		if (color is None):
			color = ['blue','red','green','orange','black','blue']
		nb = 0

		for el in list_corr_to_plot:
			TMP_xxx = el._xi1D[:,0]
			TMP_yyy = el._xi1D[:,1]
			TMP_yer = el._xi1D[:,2]
			cut = (TMP_yer>0.)
		
			TMP_cut = (TMP_yer>0.)
			if (TMP_xxx[cut].size==0):
				continue
			TMP_xxx = TMP_xxx[ TMP_cut ]
			TMP_yyy = TMP_yyy[ TMP_cut ]
			TMP_yer = TMP_yer[ TMP_cut ]
			TMP_coef = numpy.power(TMP_xxx,x_power)
			plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=5,linewidth=2, fmt='o',label=r'$'+el._name+'$', color=color[nb])
			nb += 1

		if (with_camb):
			camb = CAMB.CAMB()
			camb = camb._xi0
			coef2 = numpy.power(camb[:,0],x_power)
			plt.errorbar(camb[:,0],-0.6*coef2*camb[:,1],linewidth=2,label=r'$CAMB, \, b=-0.6$')

		### 
		xxx = self._xi1D[:,0]
		xxx = xxx[ self._xi1D[:,2]>0. ]

		if (title): plt.title(r'$'+self._title+'$', fontsize=40)
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s| \cdot '+self._label+' (|s|)$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2} \cdot '+self._label+' (|s|)$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=40, numpoints=1,ncol=1,loc=0)
		plt.show()
		
		return
	def plot_1d_with_residuals(self, x_power=0, other=[], title=True):

		fig1 = plt.figure(1)
		frame1=fig1.add_axes((.1,.3,.8,.6))

		xxx = self._xi1D[:,0]
		yyy = self._xi1D[:,1]
		yer = self._xi1D[:,2]
		
		cut = (yer>0.)
		if (xxx[cut].size==0):
			return

		xxx = xxx[ cut ]
		yyy = yyy[ cut ]
		yer = yer[ cut ]
		coef = numpy.power(xxx,x_power)
		frame1.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$'+self._name+'$', markersize=10,linewidth=2)

		for el in other:
			TMP_xxx = el._xi1D[:,0]
			TMP_yyy = el._xi1D[:,1]
			TMP_yer = el._xi1D[:,2]
		
			TMP_cut = (TMP_yer>0.)
			if (TMP_xxx[cut].size==0):
				return
			TMP_xxx = TMP_xxx[ TMP_cut ]
			TMP_yyy = TMP_yyy[ TMP_cut ]
			TMP_yer = TMP_yer[ TMP_cut ]
			TMP_coef = numpy.power(TMP_xxx,x_power)
			frame1.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, label=r'$'+el._name+'$', markersize=10,linewidth=2, fmt='o')

		if (title): plt.title(r'$'+self._title+'$', fontsize=40)
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1} \, \rm{Mpc})^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=40, numpoints=1,loc=4)
		plt.tick_params(axis='both', which='major', labelsize=20)

		###
		frame2=fig1.add_axes((.1,.1,.8,.2))
		frame2.errorbar([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ], [0.,0.],color='red')

		for el in other:
			TMP_xxx = el._xi1D[:,0]
			TMP_yyy = el._xi1D[:,1]
			TMP_yer = el._xi1D[:,2]
		
			TMP_cut = (TMP_yer>0.)
			if (TMP_xxx[cut].size==0):
				return
			TMP_xxx = TMP_xxx[ TMP_cut ]
			TMP_yyy = TMP_yyy[ TMP_cut ]
			TMP_yer = TMP_yer[ TMP_cut ]
			TMP_coef = numpy.power(TMP_xxx,x_power)

			frame2.errorbar(TMP_xxx, TMP_coef*(TMP_yyy-yyy)*numpy.sqrt(2.)/numpy.sqrt( TMP_yer*TMP_yer + yer*yer), markersize=10,linewidth=2, fmt='o',label=r'$('+el._name+'-'+self._name+')/\sigma$')

		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,True)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.legend(fontsize=30, numpoints=1,loc=1)
		plt.tick_params(axis='both', which='major', labelsize=20)
		plt.show()
		
		return
	def plot_we(self, we_index=0,x_power=0, other=[], title=True):

		list_corr_to_plot = numpy.append( [self],other )
		
		for el in list_corr_to_plot:
			cut = (el._xiWe[:,we_index,2]>0.)
			if ( el._xiWe[:,we_index,1][cut].size == 0): continue

			xxx = el._xiWe[:,we_index,0][cut]
			yyy = el._xiWe[:,we_index,1][cut]
			yer = el._xiWe[:,we_index,2][cut]
			coef = numpy.power(xxx,x_power)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$'+el._name+'$', markersize=10,linewidth=2)
		
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
			plt.legend(fontsize=30, numpoints=1,ncol=2, loc=4)
		if (x_power==1):
			plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
			plt.legend(fontsize=30, numpoints=1,ncol=2, loc=4)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1} \, \rm{Mpc})^{2}]$', fontsize=40)
			plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		
		if (title): plt.title(r'$'+self._label_wedge[we_index]+'$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,False)
		plt.show()
		
		return
	def plot_we_all(self, x_power=0, other=[], title=True):

		list_corr_to_plot = numpy.append( [self],other )
		
		for i in numpy.arange(0,self._nb_wedges):
		
			for el in list_corr_to_plot:
				cut = (el._xiWe[:,i,2]>0.)
				if ( el._xiWe[:,i,1][cut].size == 0): continue

				xxx = el._xiWe[:,i,0][cut]
				yyy = el._xiWe[:,i,1][cut]
				yer = el._xiWe[:,i,2][cut]
				coef = numpy.power(xxx,x_power)
				plt.errorbar(xxx, coef*yyy, yerr=coef*yer, marker='o', label=r'$'+self._label_wedge[i]+'$', markersize=10,linewidth=2)
		
			if (x_power==0):
				plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
				plt.legend(fontsize=30, numpoints=1,ncol=2, loc=4)
			if (x_power==1):
				plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
				plt.legend(fontsize=30, numpoints=1,ncol=2, loc=4)
			if (x_power==2):
				plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1} \, \rm{Mpc})^{2}]$', fontsize=40)
				plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		
		if (title): plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,False)
		plt.show()
		
		return
	def plot_2d(self, x_power=0):

		origin='lower'
		extent=[self._minX2D, self._maxX2D, self._minY2D, self._maxY2D]
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			origin='upper'
			extent=[self._minX2D, self._maxX2D, self._maxY2D, self._minY2D]
	
		xxx = numpy.transpose(self._xi2D[:,:,0])
		yyy = numpy.transpose(self._xi2D[:,:,1])
		yer = numpy.transpose(self._xi2D[:,:,2])

		#yer[ xxx<=40. ] = 0.
		#yer[ xxx>180. ] = 0.
	
		cut = (yer<=0.)
		if (xxx[cut].size==xxx.size):
			return
		yyy[ cut ] = float('nan')
		coef = numpy.power(xxx,x_power)
	
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xticks([ i for i in numpy.arange(self._minX2D-50., self._maxX2D+50., 50.) ])
		ax.set_yticks([ i for i in numpy.arange(self._minY2D-50., self._maxY2D+50., 50.) ])

		plt.imshow(coef*yyy, origin=origin,extent=extent, interpolation='None')
		cbar = plt.colorbar()
	
		if (x_power==0):
			cbar.set_label(r'$'+self._label+'(\, r_{\parallel},r_{\perp} \,)$',size=40)
		if (x_power==1):
			cbar.set_label(r'$|r|.'+self._label+'(\, r_{\parallel},r_{\perp} \,)$',size=40)
		if (x_power==2):
			cbar.set_label(r'$|r|^{2}.'+self._label+'(\, r_{\parallel},r_{\perp} \,)$',size=40)
	
		'''
		plt.plot( [0.,200.],[0.,4*200.],color='white',linewidth=2 )
		plt.plot( [0.,200.],[0.,-4*200.],color='white',linewidth=2 )
		plt.plot( [0.,200.],[0.,200.],color='white',linewidth=2 )
		plt.plot( [0.,200.],[0.,-200.],color='white',linewidth=2 )
		plt.xlim( [0.,200.] )
		plt.ylim( [-200.,200.] )
		'''

		#plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$r_{\perp} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.ylabel(r'$r_{\parallel} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		myTools.deal_with_plot(False,False,False)

		plt.show()
		
		return
	def plot_slice_2d(self,sliceX=None,sliceY=None, other=[]):

		list_corr = [self] + other
		i = 0
		for el in list_corr:
			i += 1

			if (sliceX is not None):
				mean = numpy.mean(el._xi2D_grid[sliceX,:,0])
				xxx  = el._xi2D_grid[sliceX,:,1]
				yyy  = el._xi2D[sliceX,:,1]
				yer  = el._xi2D[sliceX,:,2]
			elif (sliceY is not None):
				mean = numpy.mean(el._xi2D_grid[:,sliceY,1])
				xxx  = el._xi2D_grid[:,sliceY,0]
				yyy  = el._xi2D[:,sliceY,1]
                	        yer  = el._xi2D[:,sliceY,2]

			cut = (yer>0.)
			xxx = xxx[cut]
			yyy = yyy[cut]
			yer = yer[cut]
			if (xxx.size==0): continue
			#if (i==1): plt.errorbar(xxx, yyy, yerr=yer, fmt='o', markersize=10,linewidth=2,marker='o')
			plt.errorbar(xxx, yyy, yerr=yer,linewidth=2,marker='o',alpha=0.6, label=r'$'+el._name+'$')

		if (sliceX is not None):
			plt.xlabel(r'$r_{\parallel} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
			plt.title(r'$\overline{r}_{\perp} = %.2f \, [h^{-1} \, \rm{Mpc}]$' % mean)
		else:
			plt.xlabel(r'$r_{\perp} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
			plt.title(r'$\overline{r}_{\parallel} = %.2f \, [h^{-1} \, \rm{Mpc}]$' % mean)
		plt.ylabel(r'$'+self._label+'(r_{\parallel},r_{\perp})$',fontsize=40)

                myTools.deal_with_plot(False,False,True)
                plt.show()

		return
	def plot_mu(self, x_power=0):
	
		xxx = self._xiMu[:,:,0]
		muu = self._xiMu[:,:,1]
		yyy = self._xiMu[:,:,2]
		yer = self._xiMu[:,:,3]
		coef = numpy.power(xxx,x_power)
		
		cut = (yer==0.)
		if (xxx[cut].size==xxx.size):
			return
		yyy[ cut ] = float('nan')
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xticks([ i for i in numpy.arange(self._minM-1., self._maxM+1., 0.5) ])
		ax.set_yticks([ i for i in numpy.arange(self._minY2D-50., self._maxY2D+50., 50.) ])
		extent=[self._minM, self._maxM,self._min1D, self._max1D]
		
		if (x_power==0):
			a = ''
			b = ''
		if (x_power==1):
			a = '|s|.'
			b = '[h^{-1} \, \rm{Mpc}]'
		if (x_power==2):
			a = '|s|^{2}.'
			b = '[(h^{-1} \, \rm{Mpc})^{2}]'
		
		
		plt.imshow(coef*yyy, origin='lower', interpolation='None',extent=extent,aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$'+a+self._label+'\, (|s|, \mu) \, '+b+'$',size=40)
		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$\mu$', fontsize=40)
		plt.ylabel(r'$|s| \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		myTools.deal_with_plot(False,False,False)

		plt.show()
	
		return
	def plot_CAMB(self, xi1D=None, dic=None, x_power=0,distortion=True):

		if (xi1D is None): xi1D = self._xi1D
		if (dic is None):
			dic = {
				'mulpol_index' : 0,
				'start_fit'   : 40.,
				'end_fit'     : 180.,
				'b' : -1.,
				'roof' : 0.,
				'fix_roof_nul' : True,
				'guess_b' : False,
				'min_for_guess' : 20.,
				'max_for_guess' : 50.
			}

		### Get the data
		xxx = xi1D[:,0]
		yyy = xi1D[:,1]
		yer = xi1D[:,2]
	
		### Get the smooth CAMB
		camb = dic['CAMB']
		if (dic['mulpol_index']==0):   camb = camb._xi0
		elif (dic['mulpol_index']==2): camb = camb._xi2
		elif (dic['mulpol_index']==4): camb = camb._xi4

		if (distortion):
			result_1D_camb = numpy.zeros( shape=(self._nbBin1D,3) )
			result_1D_camb[:,0] = xxx
			result_1D_camb[:,1] = numpy.interp(xxx,camb[:,0],camb[:,1])

			#coef2 = numpy.power(result_1D_camb[:,0],x_power)
			#plt.plot(result_1D_camb[:,0],coef2*result_1D_camb[:,1]*yyy[4]/result_1D_camb[4,1],color='green', markersize=10,linewidth=2)

			path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_1D_'+ self._middlefix + '.txt'
			matrix = numpy.loadtxt(path)
			result_1D_camb[:,1] = numpy.dot(matrix,result_1D_camb[:,1])
			### if get the biais manually
			#result_1D_camb[:,1] *= yyy[4]/result_1D_camb[4,1]
			result_1D_camb[:,1] = dic['b']*result_1D_camb[:,1] + dic['roof']
			result_1D_camb[:,2] = 0.0000000001
		else:
			cut = (camb[:,0] <= numpy.amax(xxx)*1.1)
			size = cut[cut].size
			result_1D_camb = numpy.zeros( shape=(size,3) )
			result_1D_camb[:,0] = camb[:,0][cut]
			result_1D_camb[:,1] = dic['b']*camb[:,1][cut]+dic['roof']
			result_1D_camb[:,2] = 0.0000000001
			### if get the biais manually
			#result_1D_camb[:,1] = camb[:,1][cut]*yyy[8]/numpy.interp(xxx[8],camb[:,0],camb[:,1])
	
		### Show result
		coef = numpy.power(xxx,x_power)
		#yyy = (yyy-numpy.interp(xxx,result_1D_camb[:,0],result_1D_camb[:,1]))/yer
		#yer = 0.00000001
		plt.errorbar(xxx,coef*yyy,yerr=coef*yer,fmt='o', markersize=10,linewidth=2)
	
		coef2 = numpy.power(result_1D_camb[:,0],x_power)
		plt.plot(result_1D_camb[:,0],coef2*result_1D_camb[:,1],color='red', markersize=10,linewidth=2)
	
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|^{1}.'+self._label+' (|s|) \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1} \, \rm{Mpc})^{2}]$', fontsize=40)
			
		plt.xlabel(r'$|s| \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		#plt.ylim([ min( min( numpy.amin(coef*yyy)*1.1,numpy.amin(coef*yyy)*0.9),min( numpy.amin(coef2*result_1D_camb[:,1])*1.1,numpy.amin(coef2*result_1D_camb[:,1])*0.9) ) , max( max( numpy.amax(coef*yyy)*1.1,numpy.amax(coef*yyy)*0.9),max( numpy.amax(coef2*result_1D_camb[:,1])*1.1,numpy.amax(coef2*result_1D_camb[:,1])*0.9) ) ])
		myTools.deal_with_plot(False,False,False)
		plt.show()

		return
	def plot_multipol(self, x_power=0, other=[]):

		if (self._xiMul is None):
			self._xiMul = self.get_multipol(self._xiMu)
		for el in other:
			if (el._xiMul is None):
				el._xiMul = el.get_multipol(el._xiMu)
	
		color = ['blue','green','red','orange','black']
		ylabel = ['\\xi^{qf} (|s|)','|s|.\\xi^{qf} (|s|) \, [h^{-1} \, \rm{Mpc}]','|s|^{2}.\\xi^{qf} (|s|) \, [(h^{-1} \, \rm{Mpc})^{2}]']
	
		### Show the result
		for i in numpy.arange(0, self._nb_multipole_max ):

			##if (i!=0): continue

			cut = (self._xiMul[:,i,2]>0.)
			if ( self._xiMul[:,i,1][cut].size == 0): continue

			xxx = self._xiMul[:,i,0][cut]
			yyy = self._xiMul[:,i,1][cut]
			yer = self._xiMul[:,i,2][cut]
			coef    = numpy.power(xxx,x_power)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', markersize=10,linewidth=2, label=r'$\xi_{'+str(i)+'}$', alpha=0.6)

			for el in other:
				cut = (el._xiMul[:,i,2]>0.)
				if ( el._xiMul[:,i,1][cut].size == 0): continue

				xxx = el._xiMul[:,i,0][cut]
				yyy = el._xiMul[:,i,1][cut]
				yer = el._xiMul[:,i,2][cut]
				coef    = numpy.power(xxx,x_power)
				#plt.errorbar(xxx, coef*yyy, yerr=coef*yer,color='red', markersize=10,linewidth=2, label=r'$Mean \, simu \, \xi_{'+str(i)+'}$')
				plt.errorbar(xxx, coef*yyy, markersize=10,linewidth=2, label=r'$Mean \, simu \, \xi_{'+str(i)+'}$', alpha=0.6)
	
		plt.title(r'$'+self._title+'$', fontsize=40)
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|r|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|r|.'+self._label+' (|r|)$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|r|^{2}.'+self._label+' (|r|)$', fontsize=40)
		plt.xlabel(r'$|r| \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,True)
		plt.show()

		return
	def plot_cov_cor_matrix(self, realisation_type, dim='1D'):

		pathToLoad = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_'+dim+'.npy'

		cov = numpy.load(pathToLoad)
		myTools.plot2D(cov)
		cor = myTools.getCorrelationMatrix(cov)
		myTools.plot2D(cor)
		myTools.plotCovar([cov],[realisation_type])

		return
	def plot_cov_cor_matrix_different_method(self, realisation_type, dim='1D', path_to_other_cov=None, path_to_other_list=None):

		path_to_load = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_'
		real = [ numpy.load(path_to_load + el + '_list_'+dim+'.npy') for el in realisation_type ]
		cov  = [ numpy.load(path_to_load + el + '_cov_'+dim+'.npy') for el in realisation_type ]
		if (path_to_other_list is not None):
			for el in path_to_other_list:
				realisation_type += [ el[0] ]
				real             += [ numpy.load(el[1] + '_list_'+dim+'.npy') ]
				cov              += [ numpy.load(el[1] + '_cov_'+dim+'.npy')  ]
		if (path_to_other_cov is not None):
			for el in path_to_other_cov:
				cov  += [ numpy.load(el[1]) ]
				realisation_type = numpy.append( realisation_type, el[0] )

		### Plot the realisation
		for i in numpy.arange( len(real) ):
			print realisation_type[i]
			for j in numpy.arange(real[i][0,:].size):
				plt.errorbar(numpy.arange(real[i][:,j].size), real[i][:,j],color='blue',alpha=0.1)
			xxx = numpy.arange(real[i][:,j].size)
			yyy = numpy.mean(real[i],axis=1)
			yer = numpy.sqrt( numpy.var(real[i],axis=1)/real[i][:,j].size )
			plt.errorbar(xxx, yyy, yerr=yer,marker='o',color='red',label=r'$Mean$')
			plt.xlabel(r'$bin \, index$', fontsize=40)
			plt.ylabel(r'$\xi(|s|)$', fontsize=40)
			plt.title(r'$'+realisation_type[i]+'$', fontsize=40)
			myTools.deal_with_plot(False,False,True)
			plt.xlim([ -1., cov[i][0,:].size+1 ])
			plt.show()

			if (dim=='Mu'):
                                yyy = numpy.mean(real[i],axis=1)
                                myTools.plot2D( myTools.convert1DTo2D(yyy,self._nbBin1D,self._nbBinM))
			if (dim=='1D'):
				xxx = self._xi1D[:,0]
                        	yyy = numpy.mean(real[i],axis=1)
                        	yer = numpy.sqrt( numpy.var(real[i],axis=1)/real[i][:,j].size )
                        	plt.errorbar(xxx, yyy, yerr=yer,marker='o',color='red',label=r'$Mean$')
                        	plt.xlabel(r'$bin \, index$', fontsize=40)
                        	plt.ylabel(r'$\xi(|s|)$', fontsize=40)
                        	plt.title(r'$'+realisation_type[i]+'$', fontsize=40)
                        	myTools.deal_with_plot(False,False,True)
                        	plt.xlim([ self._min1D, self._max1D ])
                        	plt.show()
			if (dim=='2D'):
				yyy = numpy.mean(real[i],axis=1)
				myTools.plot2D( myTools.convert1DTo2D(yyy,self._nbBinX2D,self._nbBinY2D))

		### Plot diagonal
		for i in numpy.arange(len(cov)):
			plt.errorbar(numpy.arange(cov[i][0,:].size), numpy.diag(cov[i]),label=r'$'+realisation_type[i]+'$')
		plt.xlabel(r'$bin \, index$', fontsize=40)
		plt.ylabel(r'$Var(|s|)$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.xlim([ -1., cov[i][0,:].size+1 ])
		plt.show()

		realisation_type = numpy.array(realisation_type)

		### Plot Matrix
		for i in numpy.arange( len(cov) ):
			myTools.plot2D(myTools.getCorrelationMatrix(cov[i]), None,None,None,None,realisation_type[i])

		realisation_type[0] = 'Subsampling'
		myTools.plotCovar(cov,realisation_type)

		return
	def plot_distortion_matrix(self, dim='1D'):
	
		path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_'+dim+'_'+ self._middlefix + '.txt'
		matrix = numpy.loadtxt(path)

		#print matrix[ (matrix!=0.) ].size
		plt.hist( matrix[ matrix>=0.5 ], bins=100 )
		plt.show()
		plt.hist( matrix[ numpy.logical_and((matrix!=0.),matrix<=0.5) ], bins=100 )
                plt.show()

		### Blind the zero terms
		matrix[ (matrix==0.) ] = numpy.float('nan')

		myTools.plot2D(matrix)
		#myTools.plotCovar([matrix],['Distortion matrix'])

		return
	def plot_grid(self,index,minus=False):

		'''
			index :
				- 0 = s_perp
				- 1 = s_parallel
				- 2 = redshift of the pairs
		'''

		z_axis = [ '<s_{\perp}>','<s_{\parallel}>','<z_{pairs}>' ]
		if (minus):
			z_axis = [ '(<s_{\perp}> - center)/center \, [\%]','(<s_{\parallel}> - center)/center \, [\%]','(<z_{pairs}> - mean)/mean \, [\%]' ]

		origin='lower'
		extent=[self._minX2D, self._maxX2D, self._minY2D, self._maxY2D]
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			origin='upper'
			extent=[self._minX2D, self._maxX2D, self._maxY2D, self._minY2D]

		yyy = numpy.array( self._xi2D_grid[:,:,index])
		if (minus and index==0):
			for i in numpy.arange(self._nbBinX2D):
				yyy[i,:] -= self._minX2D+(i+0.5)*self._binSize
				yyy[i,:] *= 100./(self._minX2D+(i+0.5)*self._binSize)
		if (minus and index==1):
			for i in numpy.arange(self._nbBinY2D):
				yyy[:,i] -= self._minY2D+(i+0.5)*self._binSize
				yyy[:,i] *= 100./(self._minY2D+(i+0.5)*self._binSize)
		if (minus and index==2):
			mean = numpy.mean(yyy)
			yyy = 100.*(yyy-mean)/mean

		yyy = numpy.transpose(yyy)
		yer = numpy.transpose(self._xi2D[:,:,2])

		cut = (yer<=0)
		if (yyy[cut].size==yyy.size):
			return
		yyy[ cut ] = float('nan')

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_xticks([ i for i in numpy.arange(self._minX2D-50., self._maxX2D+50., 50.) ])
		ax.set_yticks([ i for i in numpy.arange(self._minY2D-50., self._maxY2D+50., 50.) ])

		plt.imshow(yyy, origin=origin,extent=extent, interpolation='None')
		cbar = plt.colorbar()
	
		#plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$s_{\perp} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		plt.ylabel(r'$s_{\parallel} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
		cbar.set_label(r'$'+z_axis[index]+'$',size=40)
		plt.grid(True)
		#cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		myTools.deal_with_plot(False,False,False)

		plt.show()

		return
	def plot_map_sub_sampling(self, euclidean=False):
	
		pathToLoad = self._path_to_txt_file_folder + self._prefix + '_map_'+ self._middlefix + '.txt'
	
		data = numpy.loadtxt(pathToLoad)
		re = data[:,1].astype(int)
		ra = data[:,2]
		de = data[:,3]
		pa = data[:,4]

		### Test if number region == self.nb_Sub_Sampling
		if (numpy.amax(re)+1 != self.nb_Sub_Sampling):
			print "  Correlation_3D::plot_map_sub_sampling::  ERROR:  numpy.amax(re)+1 != self.nb_Sub_Sampling"
			return
	
		for i in numpy.arange(0,self.nb_Sub_Sampling):
			cut = (re==i)
			plt.errorbar(ra[cut], de[cut], fmt="o")
	
		#plt.xlim([0,360.])
		#plt.ylim([-90.,90.])
		plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
		if (euclidean):
			plt.xlabel(r'$x \, [h^{-1} \, \rm{Mpc}]$')
			plt.ylabel(r'$y \, [h^{-1} \, \rm{Mpc}]$')
		else:
			plt.xlabel(r'$R.A. \, [\degree]$')
			plt.ylabel(r'$Dec. \, [\degree]$')
		myTools.deal_with_plot(False,False,False)
		plt.show()
	

		### Plot the sky with color as function of nb of pairs
		from matplotlib  import cm

		fig = plt.figure(figsize=(6,6))
		ax = fig.add_subplot(111)
		ax.grid(True,linestyle='-',color='0.75')

		if (euclidean):
			plt.xlabel(r'$x \, [h^{-1} \, \rm{Mpc}]$')
			plt.ylabel(r'$y \, [h^{-1} \, \rm{Mpc}]$')
		else:
			plt.xlabel(r'$R.A. \, [\degree]$')
			plt.ylabel(r'$Dec. \, [\degree]$')

		ax.scatter(ra,de,s=20,c=numpy.log10(pa), marker = 'o', cmap = cm.jet )

		#cbar = plt.colorbar()
		#cbar.set_label(r'$nb \, pairs \, [\#]$',size=40)

		plt.show()

		### Distribution
		plt.hist(pa, bins=1000)
		plt.xlabel(r'$nb \, pairs \, [\#]$')
		plt.ylabel(r'$\#$')
		plt.show()
		### Distribution
		plt.hist(numpy.log10(pa), bins=1000)
		plt.xlabel(r'$log( nb \, pairs ) \, [\#]$')
		plt.ylabel(r'$\#$')
		plt.show()

		return
	def plot_list_pairs_sub_sampling(self, euclidean=False):
	
		pathToLoad = self._path_to_txt_file_folder + self._prefix + '_list_pairs_'+ self._middlefix + '.txt'
	
		data = numpy.loadtxt(pathToLoad)
		re = data[:,1].astype(int)
		ra = data[:,2]
		de = data[:,3]
		pa = data[:,4]
	
		plt.errorbar(ra, de, fmt="o")

		plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
		if (euclidean):
			plt.xlabel(r'$x \, [h^{-1} \, \rm{Mpc}]$')
			plt.ylabel(r'$y \, [h^{-1} \, \rm{Mpc}]$')
		else:
			plt.xlabel(r'$R.A. \, [\degree]$')
			plt.ylabel(r'$Dec. \, [\degree]$')
		myTools.deal_with_plot(False,False,False)
		plt.show()
	

		### Plot the sky with color as function of nb of pairs
		from matplotlib  import cm

		fig = plt.figure(figsize=(6,6))
		ax = fig.add_subplot(111)
		ax.grid(True,linestyle='-',color='0.75')

		if (euclidean):
			plt.xlabel(r'$x \, [h^{-1} \, \rm{Mpc}]$')
			plt.ylabel(r'$y \, [h^{-1} \, \rm{Mpc}]$')
		else:
			plt.xlabel(r'$R.A. \, [\degree]$')
			plt.ylabel(r'$Dec. \, [\degree]$')

		ax.scatter(ra,de,s=20,c=numpy.log10(pa), marker = 'o', cmap = cm.jet )

		#cbar = plt.colorbar()
		#cbar.set_label(r'$nb \, pairs \, [\#]$',size=40)

		plt.show()

		### Distribution
		plt.hist(pa, bins=1000)
		plt.xlabel(r'$nb \, pairs \, [\#]$')
		plt.ylabel(r'$\#$')
		plt.show()
		### Distribution
		plt.hist(numpy.log10(pa), bins=1000)
		plt.xlabel(r'$log( nb \, pairs ) \, [\#]$')
		plt.ylabel(r'$\#$')
		plt.show()

		return
	def print_null_test(self,realisation_type=None,sufix=''):

		print '  corr._name'
		fix_to_zero_=False

		def fit_a_constant(yyy,cov,fix_to_zero=False):

			ones = numpy.ones(yyy.size)
			def get_chi2_null_test(a0):
                        	result = numpy.dot(numpy.dot( (yyy-a0*ones).T,numpy.linalg.inv(cov)),yyy-a0*ones)
                        	return result

			m = Minuit(get_chi2_null_test,a0=0.,error_a0=1.,print_level=-1,errordef=0.01, fix_a0=fix_to_zero)
			m.migrad()
			a0 = m.values['a0']
			chi2 = get_chi2_null_test(a0)
			print '  chi2 = ', chi2, ' nb DOF = ', yyy.size, '  p = ', chisqprob(chi2, yyy.size), '  a0 = ', a0

			return

		if ( realisation_type is None):

			### Mu
			yyy  = self._xiMu[:,:,2].flatten()
			cov  = numpy.diag( numpy.power(self._xiMu[:,:,3].flatten(),2.) )
			print '  Mu  ', fit_a_constant(yyy,cov,fix_to_zero_)
			### We
			for i in numpy.arange(self._nb_wedges):
				yyy  = self._xiWe[:,i,1]
                        	cov  = numpy.diag( numpy.power(self._xiWe[:,i,2],2.) )
				print '  We  ', fit_a_constant(yyy,cov,fix_to_zero_)
			### 1D
			yyy  = self._xi1D[:,1]
                        cov  = numpy.diag( numpy.power(self._xi1D[:,2],2.) )
			print '  1D   ', fit_a_constant(yyy,cov,fix_to_zero_)
			### 2D			
			yyy  = self._xi2D[:,:,1].flatten()
                        cov  = numpy.diag( numpy.power(self._xi2D[:,:,2].flatten(),2.) )
			print '  2D    ', fit_a_constant(yyy,cov,fix_to_zero_)
		else:
			type_cov_cor = 'cov'
			if (sufix=='_from_fit'): type_cov_cor = 'cor'

			### Mu
			#yyy  = self._xiMu[:,:,2].flatten()
			#type_cor = 'Mu'
			#path_to_load = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type +'_'+type_cov_cor+'_'+type_cor+sufix+'.npy'
			#cov  = numpy.load(path_to_load)
			#chi2 = get_chi2_null_test(yyy,cov)
			#print '  Mu    chi2 = ', chi2, ' nb DOF = ', yyy.size, '  p = ', chisqprob(chi2, yyy.size)

			### We
			#for i in numpy.arange(self._nb_wedges):
			#	yyy  = self._xiWe[:,i,1]
                        #	cov  = numpy.diag( numpy.power(self._xiWe[:,i,2],2.) )
                        #	chi2 = get_chi2_null_test(yyy,cov)
			#	print '  We',str(i),' chi2 = ', chi2, ' nb DOF = ', yyy.size, '  p = ', chisqprob(chi2, yyy.size)

			### 1D
			yyy  = self._xi1D[:,1]
			type_cor = '1D'
                        path_to_load = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type +'_'+type_cov_cor+'_'+type_cor+sufix+'.npy'
                        cov  = numpy.load(path_to_load)
			if (sufix=='_from_fit'): cov = myTools.getCovarianceMatrix(cov,numpy.power(self._xi1D[:,2],2.))
			print '  1D    ', fit_a_constant(yyy,cov,fix_to_zero_)

			### 2D			
			yyy  = self._xi2D[:,:,1].flatten()
			type_cor = '2D'
                        path_to_load = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type +'_'+type_cov_cor+'_'+type_cor+sufix+'.npy'
                        cov  = numpy.load(path_to_load)
			if (sufix=='_from_fit'): cov = myTools.getCovarianceMatrix(cov,numpy.power(self._xi2D[:,:,2].flatten(),2.))
			print '  2D    ', fit_a_constant(yyy,cov,fix_to_zero_)

		return
	def plot_Wick_calculation(self,dim,diagram_list,diagram_name):

		diagram_list = numpy.array(diagram_list)
		cov = []
		nb  = []
		we  = []
		for diagram in diagram_list:
			path_to_load = self._path_to_txt_file_folder+'xi_delta_QSO_'+dim+'_Wick_T'+diagram+'_'+self._f1+'_'+self._q1 + '.txt'
			print path_to_load
			cov += [numpy.loadtxt(path_to_load)]
			nb += [numpy.loadtxt(self._path_to_txt_file_folder+'xi_delta_QSO_'+dim+'_Wick_T'+diagram+'_number_'+self._f1+'_'+self._q1 + '.txt')]
			we += [numpy.loadtxt(self._path_to_txt_file_folder+'xi_delta_QSO_'+dim+'_Wick_T'+diagram+'_weight_'+self._f1+'_'+self._q1 + '.txt')]


		### Add other technics
		cov += [numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26/xi_delta_QSO_LYA_QSO_subsampling_cov_'+dim+'.npy')]
		diagram_name = numpy.append(diagram_name,['subSampling'])

		"""

		### Get the factor <xi_A>*<xi_B>
		#xi = numpy.matrix( self._xi2D[:,:,1].flatten() )
		#expected_value_xi = xi.T * xi
		#cor_expected_value_xi = myTools.getCorrelationMatrix(expected_value_xi)
		#myTools.plot2D(expected_value_xi)
		#myTools.plot2D(cor_expected_value_xi)

		### Get the result: Cov_A_B
		#cov -= expected_value_xi

		### Look at the proportion of calculated data
		for i in numpy.arange(diagram_list.size):
			nb[i][ (nb[i]==0.) ] = numpy.float('nan')
			myTools.plot2D(nb[i])

		### Look at the proportion of calculated data
		for i in numpy.arange(diagram_list.size):
			plt.plot( numpy.diag(nb[i]), label=r'$'+diagram_name[i]+'$')
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
                plt.show()

		### Look at the proportion of calculated data
                for i in numpy.arange(diagram_list.size):
                        plt.plot( numpy.diag(we[i]), label=r'$'+diagram_name[i]+'$')
                plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
                plt.show()
		

		"""		
		### look at the variance
		#for i in numpy.arange(diagram_name.size):
		#	plt.plot( numpy.diag(cov[i]), label=r'$'+diagram_name[i]+'$',alpha=0.7,linewidth=2)
		plt.plot( numpy.diag(cov[2]), label=r'$'+diagram_name[2]+'$',alpha=0.7,linewidth=2)
		plt.plot( numpy.diag(cov[0]), label=r'$'+diagram_name[0]+'$',alpha=0.7,linewidth=2)
		plt.plot( numpy.diag(cov[1]), label=r'$'+diagram_name[1]+'$',alpha=0.7,linewidth=2)
		#if (dim=='1D'): plt.plot( numpy.power(self._xi1D[:,2].flatten(),2.),label=r'$only \, error$' )
		#elif (dim=='2D'): plt.plot( numpy.power(self._xi2D[:,:,2].flatten(),2.),label=r'$only \, error$' )
		#plt.plot( numpy.power(self._xi2D[:,:,1].flatten(),2.) ,label=r'$<\xi_{B}> \cdot <\xi_{A}>$' )
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlabel(r'$bin \, index$', fontsize=40)
		plt.ylabel(r'$Variance$', fontsize=40)
		plt.show()

		
		### look at the variance
		if (dim=='1D'):
			coef_dist = numpy.power(self._xi1D[:,0].flatten(),-2.)
                	for i in numpy.arange(diagram_name.size):
                	        plt.errorbar( self._xi1D[:,0].flatten(), numpy.diag(cov[i])/coef_dist, label=r'$'+diagram_name[i]+'$', fmt='o')
			plt.errorbar( self._xi1D[:,0].flatten(), numpy.power(self._xi1D[:,2].flatten(),2.)/coef_dist,label=r'$only \, error$', fmt='o')
                	plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
                	plt.show()
		elif (dim=='2D'):
			coef_dist = numpy.power(self._xi2D[:,:,0].flatten(),-2.)
                        for i in numpy.arange(diagram_name.size):
                                plt.errorbar( self._xi2D[:,:,0].flatten(), numpy.diag(cov[i])/coef_dist, label=r'$'+diagram_name[i]+'$', fmt='o')
                        plt.errorbar( self._xi2D[:,:,0].flatten(), numpy.power(self._xi2D[:,:,2].flatten(),2.)/coef_dist,label=r'$only \, error$', fmt='o')
                        plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
                        plt.show()
		

		### Find the correlation function
		for i in numpy.arange(diagram_name.size):
			cor = myTools.getCorrelationMatrix(cov[i])
			cov[i][ (cov[i]==0.) ] = numpy.float('nan')
			#myTools.plot2D(cov[i])

			cor[ (cor==0.) ] = numpy.float('nan')
                	myTools.plot2D(cor)

			#if ( cor[ numpy.logical_and(cor!=1.,cor!=0.) ].size != 0 ):
			#	plt.hist( cor[numpy.logical_and(cor!=1.,cor!=0.) ], bins=500 )
			#	plt.show()

		myTools.plotCovar(cov, diagram_name)

		return
	def save_value_xi_in_txt(self):

		### 1D
		path = self._path_to_txt_file_folder + self._prefix + '_txt_1D_' + self._middlefix + '.txt'
		numpy.savetxt(path,zip(numpy.arange(self._xi1D[:,0].size),self._xi1D[:,0],self._xi1D[:,1],self._xi1D[:,2]))
		### 2D
		path = self._path_to_txt_file_folder + self._prefix + '_txt_2D_' + self._middlefix + '.txt'
		numpy.savetxt(path,zip(numpy.arange(self._xi2D[:,:,0].flatten().size),self._xi2D[:,:,0].flatten(),self._xi2D[:,:,1].flatten(),self._xi2D[:,:,2].flatten()))


		return


























































