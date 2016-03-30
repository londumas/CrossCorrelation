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

### Perso lib
import myTools
import const_delta
import const
import CAMB

raw_dic_class = {
	'minXi': 0.,
	'maxXi': 200.,
	'nbBin': 50,
	'nbBinM': 25,
	'nb_Sub_Sampling': 80,
	'size_bin_calcul_s': 1.,
	'size_bin_calcul_m': 0.02,
	'correlation': 'q_f',
	'path_to_txt_file_folder': 'NOTHING',
	'f1': 'LYA',
	'f2': 'LYA',
	'q1': 'QSO',
	'q2': 'QSO',
	'name' : 'Data'
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
	
		verbose__ = True	
	
		if (dic is None):
			dic = copy.deepcopy(raw_dic_class)

		### folder where data are
		self._path_to_txt_file_folder = dic['path_to_txt_file_folder']

		### Name of the correlation
		self._name = dic['name']

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
		
		### Mu
		self._minM = 0.
		self._maxM = 1.
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			self._minM = -self._maxM
		self._nbBinM = dic['nbBinM']
		self._nbBinM_calcul = int( (self._maxM-self._minM)/dic['size_bin_calcul_m'] )
			
		### Sub sampling
		self.nb_Sub_Sampling = int(dic['nb_Sub_Sampling'])

		### Set attributes set after
		self._xi2D_grid = None
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

		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
		xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4))
		xiWe = numpy.zeros(shape=(self._nbBin1D,3,3))
		xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		
		int_binSize = int(self._binSize)
	
		if (selection==0 or selection==1):
			### Mu
			data = numpy.loadtxt(path1D)
	
			save0 = data[:,0]
			save1 = data[:,1]
			save2 = data[:,2]
			save3 = data[:,3]
			save5 = data[:,5]
			save6 = data[:,6]
	
			tmp_save0 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save1 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save2 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save3 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save5 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save6 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
	
			tmp_save00 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save11 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save22 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save33 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save55 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save66 = numpy.zeros( shape=(self._nbBin1D,3) )
	
			tmp_save000 = numpy.zeros(self._nbBin1D)
			tmp_save111 = numpy.zeros(self._nbBin1D)
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
	
				tmp_save0[idX][idY] += save0[i]
				tmp_save1[idX][idY] += save1[i]
				tmp_save2[idX][idY] += save2[i]
				tmp_save3[idX][idY] += save3[i]
				tmp_save5[idX][idY] += save5[i]
				tmp_save6[idX][idY] += save6[i]
				
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
	
				tmp_save00[idX][idY] += save0[i]
				tmp_save11[idX][idY] += save1[i]
				tmp_save22[idX][idY] += save2[i]
				tmp_save33[idX][idY] += save3[i]
				tmp_save55[idX][idY] += save5[i]
				tmp_save66[idX][idY] += save6[i]
	
				### for xi1D
				tmp_save000[idX] += save0[i]
				tmp_save111[idX] += save1[i]
				tmp_save222[idX] += save2[i]
				tmp_save333[idX] += save3[i]
				tmp_save555[idX] += save5[i]
				tmp_save666[idX] += save6[i]

			cut = (tmp_save5!=0.)
			xiMu[:,:,0][cut] = tmp_save2[cut]/tmp_save5[cut]
			xiMu[:,:,1][cut] = tmp_save3[cut]/tmp_save5[cut]
			xiMu[:,:,2][cut] = tmp_save0[cut]/tmp_save5[cut]
			xiMu[:,:,3][cut] = numpy.sqrt( (tmp_save1[cut]/tmp_save5[cut] - xiMu[:,:,2][cut]*xiMu[:,:,2][cut])/tmp_save6[cut])

			cut = (tmp_save55!=0.)
			xiWe[:,:,0][cut] = tmp_save22[cut]/tmp_save55[cut]
			xiWe[:,:,1][cut] = tmp_save00[cut]/tmp_save55[cut]
			xiWe[:,:,2][cut] = numpy.sqrt( (tmp_save11[cut]/tmp_save55[cut] - xiWe[:,:,1][cut]*xiWe[:,:,1][cut])/tmp_save66[cut] )
	
			cut = (tmp_save555!=0.)
			xi1D[:,0][cut] = tmp_save222[cut]/tmp_save555[cut]
			xi1D[:,1][cut] = tmp_save000[cut]/tmp_save555[cut]
			xi1D[:,2][cut] = numpy.sqrt( (tmp_save111[cut]/tmp_save555[cut] - xi1D[:,1][cut]*xi1D[:,1][cut])/tmp_save666[cut] )

			if (init):
				if (verbose__): print "  ||  |s| < %u                       ||  %1.4e  ||  %1.4e    ||  %1.4e  ||" % (int(self._max1D), numpy.sum(tmp_save666), numpy.sum(tmp_save555), numpy.sum(data[:,4])/numpy.sum(data[:,5]))
			
	
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
	
			tmp_save0  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save1  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save2  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save3  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save4  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save5  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save6  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
	
			for i in range( 0,save0.size ):
				iX = i/self._nbBinY2D_calcul
				iY = i%self._nbBinY2D_calcul
	
				idX = iX/int_binSize
				idY = iY/int_binSize
	
				tmp_save0[idX][idY] += save0[i]
				tmp_save1[idX][idY] += save1[i]
				tmp_save2[idX][idY] += save2[i]
				tmp_save3[idX][idY] += save3[i]
				tmp_save4[idX][idY] += save4[i]
				tmp_save5[idX][idY] += save5[i]
				tmp_save6[idX][idY] += save6[i]

			cut = (tmp_save5!=0.)
			xi2D[:,:,0][cut] = numpy.sqrt( (tmp_save2[cut]/tmp_save5[cut])**2. + (tmp_save3[cut]/tmp_save5[cut])**2. )
			xi2D[:,:,1][cut] = tmp_save0[cut] / tmp_save5[cut]
			xi2D[:,:,2][cut] = numpy.sqrt( (tmp_save1[cut]/tmp_save5[cut] - xi2D[:,:,1][cut]*xi2D[:,:,1][cut])/tmp_save6[cut] )

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

		'''
	
		selection:
			- 0: do all
			- 1: only 1D
			- 2: only 2D
	
		'''
		print path1D

		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3,3))
		xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4,3))
		xiWe = numpy.zeros(shape=(self._nbBin1D,3,3,3))
		xi1D = numpy.zeros(shape=(self._nbBin1D,3,3))
		
		int_binSize = int(self._binSize)
	
		if (selection==0 or selection==1):
			### Mu
			data = numpy.loadtxt(path1D)
			save0 = data[:,0]
	
			tmp_save0 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM,3) )
			tmp_save1 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM,3) )
			tmp_save2 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save3 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save5 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
			tmp_save6 = numpy.zeros( shape=(self._nbBin1D,self._nbBinM) )
	
			tmp_save00 = numpy.zeros( shape=(self._nbBin1D,3,3) )
			tmp_save11 = numpy.zeros( shape=(self._nbBin1D,3,3) )
			tmp_save22 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save33 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save55 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save66 = numpy.zeros( shape=(self._nbBin1D,3) )
	
			tmp_save000 = numpy.zeros( shape=(self._nbBin1D,3) )
			tmp_save111 = numpy.zeros( shape=(self._nbBin1D,3) )
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
			for i in numpy.arange(0,3):
				xiMu[:,:,0,i][cut] = tmp_save2[cut]/tmp_save5[cut]
				xiMu[:,:,1,i][cut] = tmp_save3[cut]/tmp_save5[cut]
				xiMu[:,:,2,i][cut] = tmp_save0[:,:,i][cut]/tmp_save5[cut]
				xiMu[:,:,3,i][cut] = numpy.abs(xiMu[:,:,2,i][cut])/numpy.sqrt(tmp_save6[cut])

			cut = (tmp_save66>1.)
			for i in numpy.arange(0,3):
				xiWe[:,:,0,i][cut] = tmp_save22[cut]/tmp_save55[cut]
				xiWe[:,:,1,i][cut] = tmp_save00[:,:,i][cut]/tmp_save55[cut]
				xiWe[:,:,2,i][cut] = numpy.abs(xiWe[:,:,1,i][cut])/numpy.sqrt(tmp_save66[cut])
			
			cut = (tmp_save666>1.)
			for i in numpy.arange(0,3):
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
	
			tmp_save0  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,3) )
			tmp_save1  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D,3) )
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
			for i in numpy.arange(0,3):
				xi2D[:,:,0,i][cut] = numpy.sqrt( (tmp_save2[cut]/tmp_save5[cut])**2. + (tmp_save3[cut]/tmp_save5[cut])**2. )
				xi2D[:,:,1,i][cut] = tmp_save0[:,:,i][cut] / tmp_save5[cut]
				xi2D[:,:,2,i][cut] = numpy.abs(xi2D[:,:,1,i][cut])/numpy.sqrt(tmp_save6[cut])
		
		return xiMu, xiWe, xi1D, xi2D
	def empty_data(self):

		### Set attributes set after
		self._meanZ = 0.
		self._xiMul = numpy.zeros(shape=(self._nbBin1D,5,3))

		### Correlations data
		self._xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4))
		self._xiWe = numpy.zeros(shape=(self._nbBin1D,3,3))
		self._xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		self._xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))

		return
	def save_list_realisation(self, realisation_type, nb_realisation):

		if ( realisation_type=='subsampling' and (nb_realisation != self.nb_Sub_Sampling) ):
			print "  Correlation_3D::saveListReal::  WARNING:  nb_realisation != self.nb_Sub_Sampling"

		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,nb_realisation) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,3,nb_realisation) )
		list1D       = numpy.zeros( shape=(self._nbBin1D,nb_realisation) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,nb_realisation) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,5,nb_realisation) )

		path1D = self._path_to_txt_file_folder + self._prefix + '_Mu_' + self._middlefix + '_' + realisation_type + '_'
		path2D = self._path_to_txt_file_folder + self._prefix + '_2D_' + self._middlefix + '_' + realisation_type + '_'
		pathToSave = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_'

		for i in numpy.arange(nb_realisation):
			print i

			xiMu, xiWe, xi1D, xi2D = self.read_data(path1D+str(i)+'.txt',path2D+str(i)+'.txt')
			list1D[:,i]         = xi1D[:,1]
			list2D[:,i]         = xi2D[:,:,1].flatten()
			listMu[:,i]         = xiMu[:,:,2].flatten()
			listWe[:,:,i]       = xiWe[:,:,1]
			listMultipol[:,:,i] = self.get_multipol(xiMu)[:,:,1]
			'''
			except:
				print '   ERROR:: ', path1D+str(i)+'.txt',path2D+str(i)+'.txt'
				commandProd = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Correlation/bin/main.exe"
				tmp_command = "clubatch \"time ; hostname ; " + commandProd + ' 5 0 ' + str(i) + ' 0 0 0 ' + "\""
				subprocess.call(tmp_command, shell=True)
			'''

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
	def save_list_realisation_simulation(self, dic_class, dic_simu):

		pathToSave = dic_simu['path_to_simu'] + 'Results'
		pathToSave += dic_simu['prefix']
		pathToSave += '/'
		pathToSave += self._prefix + '_' + self._middlefix + '_result_'

		nb_realisation = dic_simu['nb_box']*dic_simu['nb_simu']

		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,nb_realisation) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,3,nb_realisation) )
		list1D       = numpy.zeros( shape=(self._nbBin1D,nb_realisation) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,nb_realisation) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,5,nb_realisation) )
		listGrid     = numpy.zeros( shape=(self._nbBin2D,3,nb_realisation) )

		nb = 0
		for i in range(0,dic_simu['nb_box']):
			for j in range(0,dic_simu['nb_simu']):

				raw = dic_simu['path_to_simu'] + 'Box_00' + str(i) + '/Simu_00' + str(j) +'/Results'
				raw += dic_simu['prefix']
				raw += '/'
				dic_class['path_to_txt_file_folder'] = raw

				try:
					corr = Correlation3D(dic_class)
				except:
					print i, j
					continue

				listMu[:,nb]         = corr._xiMu[:,:,2].flatten()
				listWe[:,:,nb]       = corr._xiWe[:,:,1]
				list1D[:,nb]         = corr._xi1D[:,1]
				list2D[:,nb]         = corr._xi2D[:,:,1].flatten()
				listMultipol[:,:,nb] = corr.get_multipol(corr._xiMu)[:,:,1]

				listGrid[:,0,nb] = corr._xi2D_grid[:,:,0].flatten()
				listGrid[:,1,nb] = corr._xi2D_grid[:,:,1].flatten()
				listGrid[:,2,nb] = corr._xi2D_grid[:,:,2].flatten()

				nb += 1

		listMu       = listMu[:,:nb]
		listWe       = listWe[:,:,:nb]
		list1D       = list1D[:,:nb]
		list2D       = list2D[:,:nb]
		listMultipol = listMultipol[:,:,:nb]
		listGrid     = listGrid[:,:,:nb]

		print listMu[0,:].size, nb

		numpy.save(pathToSave+'list_Mu',listMu)
		numpy.save(pathToSave+'list_We',listWe)
		numpy.save(pathToSave+'list_1D',list1D)
		numpy.save(pathToSave+'list_2D',list2D)
		numpy.save(pathToSave+'list_Multipol',listMultipol)
		numpy.save(pathToSave+'list_Grid',listGrid)

		covMu = numpy.cov(listMu)
		cov1D = numpy.cov(list1D)
		cov2D = numpy.cov(list2D)

		numpy.save(pathToSave+'cov_Mu',covMu)
		numpy.save(pathToSave+'cov_1D',cov1D)
		numpy.save(pathToSave+'cov_2D',cov2D)

		return
	def set_values_on_mean_simulation(self, dic_simu):

		path_to_load = dic_simu['path_to_simu'] + 'Results'
		path_to_load += dic_simu['prefix']
		path_to_load += '/'
		if (self._prefix=='xi_QSO_QSO'):
			path_to_load += self._prefix + '_result_'
		else:
			path_to_load += self._prefix + '_' + self._middlefix + '_result_'

		### mu
		listMu = numpy.load(path_to_load+'list_Mu.npy')
		nb_realisation = listMu[0,:].size
		self._xiMu[:,:,2] = myTools.convert1DTo2D( numpy.mean(listMu,axis=1) ,self._nbBin1D,self._nbBinM)
		self._xiMu[:,:,3] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.cov(listMu))/nb_realisation ),self._nbBin1D,self._nbBinM)
		### we
		listWe = numpy.load(path_to_load+'list_We.npy')
		for i in numpy.arange(3):
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
		for i in numpy.arange(5):
			self._xiMul[:,i,1] = numpy.mean(listMultipol[:,i,:],axis=1)
			self._xiMul[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listMultipol[:,i,:] )/nb_realisation))
		### Grid
		listGrid = numpy.load(path_to_load+'list_Grid.npy')
		self._xi2D_grid[:,:,0] = myTools.convert1DTo2D( numpy.mean(listGrid[:,0,:],axis=1),self._nbBinX2D,self._nbBinY2D )
		self._xi2D_grid[:,:,1] = myTools.convert1DTo2D( numpy.mean(listGrid[:,1,:],axis=1),self._nbBinX2D,self._nbBinY2D )
		self._xi2D_grid[:,:,2] = myTools.convert1DTo2D( numpy.mean(listGrid[:,2,:],axis=1),self._nbBinX2D,self._nbBinY2D )
		self._meanZ = numpy.mean(self._xi2D_grid[:,:,2])

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
		for i in numpy.arange(3):
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
		for i in numpy.arange(5):
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
	
		### Array with name of variable
		if (dic is None):
			dic = {
				'xi0' : True,
				'xi1' : False,
				'xi2' : True,
				'xi3' : False,
				'xi4' : False,
				'plot' : False
			}
		nameArray = [ 'xi0','xi1','xi2','xi3','xi4']
		nbXi = len(nameArray)

		### Get the data
		xxx = xiMu[:,:,0]
		muu = xiMu[:,:,1]
		yyy = xiMu[:,:,2]
		yer = xiMu[:,:,3]

		### Keep the results
		result_xi = numpy.zeros(shape=(self._nbBin1D,nbXi,3))
		meanXXX   = numpy.mean(xxx,axis=1)
		for i in numpy.arange(0,nbXi):
			result_xi[:,i,0] = meanXXX
	
		for i in numpy.arange(self._nbBin1D):
	
			cut         = (muu[i,:]!=0.)
			tmpyyy      = yyy[i,:][cut]
			tmpyer      = yer[i,:][cut]
			xxxMu       = muu[i,:][cut]
			xxxMuPower1 = numpy.power(xxxMu,1.)
			xxxMuPower2 = numpy.power(xxxMu,2.)
			xxxMuPower3 = numpy.power(xxxMu,3.)
			xxxMuPower4 = numpy.power(xxxMu,4.)
			
			### Define the fit function
			def chi2(xi0,xi1,xi2,xi3,xi4):
				fit = (xi0 - 0.5*xi2 + 0.375*xi4) + (xi1-1.5*xi3)*xxxMuPower1 + (1.5*xi2 - 3.75*xi4)*xxxMuPower2 + 2.5*xi3*xxxMuPower3 + 4.375*xi4*xxxMuPower4
				return numpy.sum( numpy.power( (tmpyyy-fit)/tmpyer ,2.) )
	
			### Init ad perform the fit
			m = Minuit(chi2, xi0=0.,error_xi0=0.1,xi1=0.,error_xi1=0.1,xi2=0.,error_xi2=0.1,xi3=0.,error_xi3=0.1,xi4=0.,error_xi4=0.1, print_level=-1, errordef=0.01,
				fix_xi0=not dic['xi0'],fix_xi1=not dic['xi1'],fix_xi2=not dic['xi2'], fix_xi3=not dic['xi3'],fix_xi4=not dic['xi4'])
			m.migrad()
	
			### Keep the results
			for j in numpy.arange(0,nbXi):
				result_xi[i,j,1] = m.values[ nameArray[j] ]
				result_xi[i,j,2] = m.errors[ nameArray[j] ]

			if (dic['plot']):
				xi0 = m.values[ nameArray[0] ]
				xi1 = m.values[ nameArray[1] ]
				xi2 = m.values[ nameArray[2] ]
				xi3 = m.values[ nameArray[3] ]
				xi4 = m.values[ nameArray[4] ]
				print xi0,xi1,xi2,xi3,xi4
				plt.errorbar(xxxMu,tmpyyy,yerr=tmpyer, fmt='o')
				plt.errorbar(xxxMu,(xi0 - 0.5*xi2 + 0.375*xi4) + (xi1-1.5*xi3)*xxxMuPower1 + (1.5*xi2 - 3.75*xi4)*xxxMuPower2 + 2.5*xi3*xxxMuPower3 + 4.375*xi4*xxxMuPower4)
				plt.show()

		cut = (result_xi[:,:,2]==0.1)
		result_xi[cut] = 0.

		return result_xi
	def get_multipol_index(self, index):

		if (self._xiMul is None): self._xiMul = self.get_multipol(self._xiMu)

		xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		xi1D[:,:] = self._xiMul[:,index,:]

		return xi1D
	def get_CAMB(self,distortion=True):

		### Get Camb
		xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		data_camb = numpy.loadtxt(const.pathToCamb__)
		xi1D[:,0] = self._xi1D[:,0]
		xi1D[:,1] = numpy.interp(xi1D[:,0],data_camb[1:,0],data_camb[1:,1])
		xi1D[:,2] = 0.00000001

		if (distortion):
			path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_1D_'+ self._middlefix + '.txt'
			matrix = numpy.loadtxt(path)
			xi1D[:,1] = numpy.dot(matrix,xi1D[:,1])

		return xi1D
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
		camb = CAMB.CAMB('CHRISTOPHE')
		if (dic['mulpol_index']==0):
			yyy_Camb = numpy.interp(xi1D[:,0],camb._xi0[:,0],camb._xi0[:,1])
		elif (dic['mulpol_index']==2):
			yyy_Camb = numpy.interp(xi1D[:,0],camb._xi2[:,0],camb._xi2[:,1])

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
	def write_metal_model(self, with_metals_templates=False, plot=False):

		path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_q_f__'+self._f1+'__'+self._q1
		if (with_metals_templates):
			path_to_BAOFIT += '__withMetalsTemplates'
		path_to_BAOFIT += '/bao2D_QSO_'

		suffix = ['.0.dat','.2.dat','.4.dat']
		name   = ['SiII(1260)','SiIII(1207)','SiII(1193)','SiII(1190)']
		name_to_save = ['Si2c','Si3','Si2b','Si2a']

		for j in numpy.arange(len(name)):
			xiMu, xiWe, xi1D, xi2D = self.read_metal_model(name[j])

			if (plot): myTools.plot_2d_correlation(xi2D[:,:,:,0],0)

			for i in range(0,3):
				correlation = numpy.zeros( shape=(self._nbBin2D,4) )
        	        	indexMatrix = numpy.arange(self._nbBin2D)
        	        	correlation[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
        	        	correlation[:,1] = xi2D[:,:,1,i].flatten()
        	        	cutCorrelation = (correlation[:,1]!=0.)
				if (correlation[:,0][cutCorrelation].size==0): continue

				### Get the array and sort it
				correlation2 = numpy.zeros( shape=(self._nbBin2D,4) )
				for k in numpy.arange(self._nbBin2D):
					idx = correlation[k,0].astype(int)
					correlation2[idx,0] = idx
					correlation2[idx,1] = correlation[k,1]
					correlation2[idx,2] = self._minX2D+0.5*self._binSize + (idx%self._nbBinX2D)*self._binSize
					correlation2[idx,3] = self._minY2D+0.5*self._binSize + (idx/self._nbBinX2D)*self._binSize

				print path_to_BAOFIT + name_to_save[j] + suffix[i]
				numpy.savetxt( path_to_BAOFIT + name_to_save[j] + suffix[i],zip(correlation2[:,2],correlation2[:,3],correlation2[:,1]),fmt='%u %u %1.20e')

		return
	def write_BAOFIT_ini_and_grid(self, param, dist_matrix=False, with_metals_templates=False):

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

		if (with_metals_templates):
			path_to_BAOFIT += '__withMetalsTemplates'
		path_to_BAOFIT += '/bao2D'

		string_ini = """

## Linear theory P(k) templates with and w/o wiggles
modelroot = """ + const.path_to_BAOFIT_model__ + """
#fiducial =  DR9LyaMocksLCDM
#nowiggles = DR9LyaMocksLCDMSB

fiducial =  DR9LyaMocks
nowiggles = DR9LyaMocksSB

#modelroot = /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Mock_JMLG/Produce_CAMB/
#fiducial  = christophe_Eisentien_Hu_Simu
#nowiggles = christophe_Eisentien_Hu_Simu_SB

omega-matter    = 0.27
hubble-constant = 0.7
sigma8 = 0.794961
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
"""+str_dist_matrix+"""dist-matrix = yes
dist-matrix-order = """+str(self._nbBin2D)+"""

# Parameter setup
model-config = value[beta]=                      """+param[0] +""";
model-config = """+str_corr2+"""[(1+beta)*bias]= """+param[1] +""";
#model-config = value[(1+beta)*bias]= """+param[1] +""";
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
model-config = fix[pixel scale]=0.7;

## Metal correlations
"""+str_metals+"""metal-model-interpolate = true
"""+str_metals+"""metal-model-name = """+path_to_BAOFIT+"""
"""+str_metals+"""model-config = value[beta Si2a]=1.4;
"""+str_metals+"""model-config = value[bias Si2a]=-0.01;
"""+str_metals+"""model-config = value[beta Si2b]=1.4;
"""+str_metals+"""model-config = value[bias Si2b]=-0.01;
"""+str_metals+"""model-config = value[beta Si2c]=1.4;
"""+str_metals+"""model-config = value[bias Si2c]=-0.01;
"""+str_metals+"""model-config = fix[beta Si3]=1.4;
"""+str_metals+"""model-config = fix[bias Si3]=-0.01;
"""+str_metals+"""model-config = gaussprior[beta Si2a] @ (0,2.8);
"""+str_metals+"""model-config = gaussprior[beta Si2b] @ (0,2.8);
"""+str_metals+"""model-config = gaussprior[beta Si2c] @ (0,2.8);
#model-config = boxprior[beta Si2c] @ (0.00001,2.8);
#model-config = gaussprior[beta Si3] @ (0,2.8);


## 2D chisq scan in BAO parameters
#model-config = binning[BAO alpha-parallel] ={0.7:1.4}*50
#model-config = binning[BAO alpha-perp]     ={0.7:1.4}*50	
model-config = binning[BAO alpha-parallel] ={0.98:1.02}*50
model-config = binning[BAO alpha-perp]     ={0.98:1.02}*50

## Maximum allowed radial dilation (increases the range that model needs to cover)
dilmin = 0.01
#dilmax = 2.5
dilmax = 3.5

## Non-linear broadening with 1+f = (SigmaNL-par)/(SigmaNL-perp)

# boxprior keeps result positive (since model only depends on squared value)
model-config = boxprior[SigmaNL-perp] @ (0,6);
# un-comment next line to broaden all scales (default is peak only)
#nl-broadband = true

# Broadband distortion model
"""+str_corr+"""dist-add = rP,rT=0:2,-3:1
"""+str_corr3+"""dist-add = -2:0,0:4:2,0

### Data Options #############################################

## Data to analyze
data = """+path_to_BAOFIT+"""
dist-matrix-name = """+path_to_BAOFIT+"""

## Data format
data-format = comoving-cartesian
axis1-bins = ["""+ str(self._minY2D) + ':' + str(self._maxY2D) + """]*"""+ str(self._nbBinY2D) +"""
axis2-bins = ["""+ str(self._minX2D) + ':' + str(self._maxX2D) + """]*"""+ str(self._nbBinX2D) +"""
axis3-bins = {"""+ str(self._meanZ) +"""}

### Analysis Options #########################################

# Cuts to apply before fitting
rmin = 40
rmax = 180
#rperp-min = 5

# Generate a second set of outputs with the additive distortion turned off
alt-config = fix[dist*]=0

# Do not dump multipoles (since the distortion model multipole integrals are singular)
ndump = 0

# Prefix to use for all analysis output files
output-prefix = """ + path_to_BAOFIT + """.
"""

		text_file = open(path_to_BAOFIT+'.ini', "w")
		text_file.write(string_ini)
		text_file.close()
		

		return
	def send_BAOFIT(self, realisation_type, correlation_matrix_path=None, saving_data=True, scan=False, toyMC=0, dist_matrix=False, only_diagonal=False, with_metals_templates=False):
		"""

			Send the fit of the BAO with BAOFIT
			realisation_type:
			cor:=None

		"""

		if (self._correlation=='q_q'):
			param = numpy.asarray( [0.268344,0.962524,-2.,0.,0.,3.6,0.,0.,1.966,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.] )
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._q1
		elif (self._correlation=='q_f'):
			param = numpy.asarray( [1.2,-0.351,0.9,0.,0.,3.6,0.962524,1.966,3.26,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.] )
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1
		elif (self._correlation=='f_f'):
			param = numpy.asarray( [1.4,-0.351,3.8,0.,0.,3.6,0.962524,1.966,3.26,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.] )
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1

		if (with_metals_templates):
			path_to_BAOFIT += '__withMetalsTemplates'
		path_to_BAOFIT += '/'



		subprocess.call('mkdir ' + path_to_BAOFIT, shell=True)
		path_to_BAOFIT += 'bao2D.'

		self.write_BAOFIT_ini_and_grid(param,dist_matrix, with_metals_templates)

		if (saving_data):

			### Save .grid
			print '  Saving .grid'
			grid = numpy.zeros( shape=(self._nbBin2D,4) )
			indexMatrix = numpy.arange(self._nbBin2D)
			grid[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
			grid[:,1] = self._xi2D_grid[:,:,1].flatten()
			grid[:,2] = self._xi2D_grid[:,:,0].flatten()
			grid[:,3] = self._xi2D_grid[:,:,2].flatten()
			numpy.savetxt(path_to_BAOFIT+'grid',zip(grid[:,0],grid[:,1],grid[:,2],grid[:,3]),fmt='%u %1.20e %1.20e %1.20e')
			del grid

			### Save .data
			print '  Saving .data'
			correlation = numpy.zeros( shape=(self._nbBin2D,2) )
			indexMatrix = numpy.arange(self._nbBin2D)
			correlation[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
			correlation[:,1] = self._xi2D[:,:,1].flatten()
			cutCorrelation = (correlation[:,1]!=0.)
			numpy.savetxt( path_to_BAOFIT + 'data',zip(correlation[:,0][cutCorrelation],correlation[:,1][cutCorrelation]),fmt='%u %1.20e')
			del correlation, indexMatrix, cutCorrelation

			### Save .cov
			if (self._correlation=='q_q'):
				path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results/xi_QSO_QSO_result_cov_2D.npy'
			elif (self._correlation=='q_f'):
				#path_to_cov = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_2D.npy'
				#path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_PureRaw/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
				#path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_nicolasEstimator/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
				#path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Results_PureRaw/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
				#path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_Raw/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
				#path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_noRSD_PureRaw/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
				path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
			elif (self._correlation=='f_f'):
				path_to_cov = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_2D.npy'
				path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results/xi_A_delta_delta_LYA_result_cov_2D.npy'

			if (correlation_matrix_path is None): correlation_matrix_path = path_to_cov
			print '  correlation matrix path = ', correlation_matrix_path

			cov = numpy.load(path_to_cov)
			cor = myTools.getCorrelationMatrix( numpy.load(correlation_matrix_path) )

			print '  Saving .cov'
			cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))
			#cov /= 88.
			cov /= 100.
			#cov /= 30.
			#cov /= 86.
			if (only_diagonal):

				mean_cor_off_diag = self._nbBin2D/(self._nbBin2D-1.)*numpy.mean( cor-numpy.diag( numpy.diag(cor)) )
				print '  The mean off diagonal term is = ', mean_cor_off_diag
				mean_cor_matrix = numpy.ones( shape=(self._nbBin2D,self._nbBin2D) )*mean_cor_off_diag
				mean_cor_matrix -= numpy.diag(numpy.diag(mean_cor_matrix))
				cor = numpy.diag( numpy.diag(cor) ) + mean_cor_matrix
				cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))
				#cov = numpy.diag( numpy.diag(cov) )

			covarianceMatrix = numpy.zeros( shape=(self._nbBin2D*self._nbBin2D,3) )
			indexMatrix1 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)/self._nbBin2D
			indexMatrix2 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)%self._nbBin2D
			indexMatrix1 = (indexMatrix1%self._nbBinY2D)*self._nbBinX2D + indexMatrix1/self._nbBinY2D
			indexMatrix2 = (indexMatrix2%self._nbBinY2D)*self._nbBinX2D + indexMatrix2/self._nbBinY2D
			covarianceMatrix[:,0] = numpy.triu(indexMatrix1,k=0).flatten()
			covarianceMatrix[:,1] = numpy.triu(indexMatrix2,k=0).flatten()
			covarianceMatrix[:,2] = numpy.triu(cov,k=0).flatten()
			cutCovarMatrix = (covarianceMatrix[:,2]!=0.)
			'''
			if (covarianceMatrix[:,2][cutCovarMatrix].size != int(self._nbBin2D*(self._nbBin2D+1.)/2.) ):
				print '  xi_delta_QSO.py::prepareForBAOFIT()  size of covariance matrix is incorrect'
				print '  size covariance matrix', cutCovarMatrix[:,2][cutCovarMatrix].size
				print '  size it should have', self._nbBin2D*(self._nbBin2D+1.)/2.
				return
			'''
			numpy.savetxt( path_to_BAOFIT + 'cov',zip(covarianceMatrix[:,0][cutCovarMatrix],covarianceMatrix[:,1][cutCovarMatrix],covarianceMatrix[:,2][cutCovarMatrix]),fmt='%u %u %1.20e')
			del cov, covarianceMatrix, cutCovarMatrix

			### Save .dmat
			if (dist_matrix):
				print '  Saving .dmat'
				path_to_distortion_matrix = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_2D_'+ self._middlefix + '.txt'
				dmatData = numpy.loadtxt(path_to_distortion_matrix)
				distortionMatrix = numpy.zeros( shape=(self._nbBin2D*self._nbBin2D,3) )
				distortionMatrix[:,0] = indexMatrix1.flatten()
				distortionMatrix[:,1] = indexMatrix2.flatten()
				distortionMatrix[:,2] = dmatData.flatten()
				cutDistortionMatrix = (distortionMatrix[:,2]!=0.)
				numpy.savetxt( path_to_BAOFIT + 'dmat',zip(distortionMatrix[:,0][cutDistortionMatrix], distortionMatrix[:,1][cutDistortionMatrix], distortionMatrix[:,2][cutDistortionMatrix]),fmt='%u %u %1.20e')
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
	def plot_1d(self, x_power=0, other=[], title=True):

		with_camb = False

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
		plt.errorbar(xxx, coef*yyy, yerr=coef*yer, marker='o', markersize=10,linewidth=2, label=r'$'+self._name+'$')

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
			plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, marker='o',label=r'$'+el._name+'$',color='red')

		if (with_camb):
			camb = CAMB.CAMB()
			camb = camb._xi0
			coef2 = numpy.power(camb[:,0],x_power)
			plt.errorbar(camb[:,0],-0.6*coef2*camb[:,1],linewidth=2,label=r'$CAMB, \, b=-0.6$')

		if (title): plt.title(r'$'+self._title+'$', fontsize=40)
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,True)
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
			plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
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
	def plot_we(self, x_power=0, other=[], title=True):
	
		label = ['0.8 < |\mu|', '0.5 < |\mu| \leq 0.8', '|\mu| \leq 0.5']
		
		
		for i in numpy.arange(0,3):
		
			cut = (self._xiWe[:,i,2]>0.)
			if (self._xiWe[:,i,0][cut].size==0):
				continue
		
			xxx = self._xiWe[:,i,0][cut]
			yyy = self._xiWe[:,i,1][cut]
			yer = self._xiWe[:,i,2][cut]
			coef = numpy.power(xxx,x_power)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$Data \, '+label[i]+'$', markersize=10,linewidth=2)

			for el in other:
				cut = (el._xiWe[:,i,2]>0.)
				if ( el._xiWe[:,i,1][cut].size == 0): continue

				xxx = el._xiWe[:,i,0][cut]
				yyy = el._xiWe[:,i,1][cut]
				yer = el._xiWe[:,i,2][cut]
				coef = numpy.power(xxx,x_power)
				plt.errorbar(xxx, coef*yyy, yerr=coef*yer, label=r'$Simu \,'+label[i]+'$',color='red', markersize=10,linewidth=2)
		
			if (x_power==0):
				plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
				plt.legend(fontsize=30, numpoints=1,ncol=2, loc=4)
			if (x_power==1):
				plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
				plt.legend(fontsize=30, numpoints=1,ncol=2, loc=4)
			if (x_power==2):
				plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
				plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		
		if (title): plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
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
	
		cut = (yer==0)
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
			cbar.set_label(r'$'+self._label+'(\, \overrightarrow{s} \,)$',size=40)
		if (x_power==1):
			cbar.set_label(r'$|s|.'+self._label+'(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
		if (x_power==2):
			cbar.set_label(r'$|s|^{2}.'+self._label+'(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)
	
		'''
		plt.plot( [0.,200.],[0.,4*200.],color='white',linewidth=2 )
		plt.plot( [0.,200.],[0.,-4*200.],color='white',linewidth=2 )
		plt.plot( [0.,200.],[0.,200.],color='white',linewidth=2 )
		plt.plot( [0.,200.],[0.,-200.],color='white',linewidth=2 )
		plt.xlim( [0.,200.] )
		plt.ylim( [-200.,200.] )
		'''

		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
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
			plt.errorbar(xxx, yyy, yerr=yer,linewidth=2,marker='o')

		if (sliceX is not None):
			plt.xlabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
			plt.title(r'$<s_{\perp}> = %2f$' % mean)
		else:
			plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
			plt.title(r'$<s_{\parallel}> = %2f$' % mean)
		plt.ylabel(r'$'+self._label+'(\, \overrightarrow{s} \,)$',fontsize=40)

                myTools.deal_with_plot(False,False,False)
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
			b = '[h^{-1}.Mpc]'
		if (x_power==2):
			a = '|s|^{2}.'
			b = '[(h^{-1}.Mpc)^{2}]'
		
		
		plt.imshow(coef*yyy, origin='lower', interpolation='None',extent=extent,aspect='auto')
		cbar = plt.colorbar()
		cbar.set_label(r'$'+a+self._label+'\, (|s|, \mu) \, '+b+'$',size=40)
		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$\mu$', fontsize=40)
		plt.ylabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
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
		camb = CAMB.CAMB('CHRISTOPHE')
		if (dic['mulpol_index']==0): camb = camb._xi0
		elif (dic['mulpol_index']==2): camb = camb._xi2

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
			#result_1D_camb[:,1] *= yyy[8]/result_1D_camb[4,1]
			result_1D_camb[:,1] = dic['b']*result_1D_camb[:,1] + dic['roof']
		else:
			cut = (camb[:,0] <= numpy.amax(xxx)*1.1)
			size = cut[cut].size
			result_1D_camb = numpy.zeros( shape=(size,3) )
			result_1D_camb[:,0] = camb[:,0][cut]
			result_1D_camb[:,1] = dic['b']*camb[:,1][cut]+dic['roof']
			result_1D_camb[:,2] = 0.0000000001
			### if get the biais manually
			#result_1D_camb[:,1] = camb[:,1][cut]*yyy[4]/numpy.interp(xxx[4],camb[:,0],camb[:,1])
	
		### Show result
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,yerr=coef*yer,fmt='o', markersize=10,linewidth=2)
	
		coef2 = numpy.power(result_1D_camb[:,0],x_power)
		plt.plot(result_1D_camb[:,0],coef2*result_1D_camb[:,1],color='red', markersize=10,linewidth=2)
	
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|^{1}.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
			
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		plt.ylim([ min( min( numpy.amin(coef*yyy)*1.1,numpy.amin(coef*yyy)*0.9),min( numpy.amin(coef2*result_1D_camb[:,1])*1.1,numpy.amin(coef2*result_1D_camb[:,1])*0.9) ) , max( max( numpy.amax(coef*yyy)*1.1,numpy.amax(coef*yyy)*0.9),max( numpy.amax(coef2*result_1D_camb[:,1])*1.1,numpy.amax(coef2*result_1D_camb[:,1])*0.9) ) ])
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
		ylabel = ['\\xi^{qf} (|s|)','|s|.\\xi^{qf} (|s|) \, [h^{-1}.Mpc]','|s|^{2}.\\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]']
	
		### Show the result
		for i in numpy.arange(0, self._xiMul[0,:,0].size ):

			cut = (self._xiMul[:,i,2]>0.)
			if ( self._xiMul[:,i,1][cut].size == 0): continue

			xxx = self._xiMul[:,i,0][cut]
			yyy = self._xiMul[:,i,1][cut]
			yer = self._xiMul[:,i,2][cut]
			coef    = numpy.power(xxx,x_power)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o',color='blue', markersize=10,linewidth=2, label=r'$Data \, \xi_{'+str(i)+'}$')

			for el in other:
				cut = (el._xiMul[:,i,2]>0.)
				if ( el._xiMul[:,i,1][cut].size == 0): continue

				xxx = el._xiMul[:,i,0][cut]
				yyy = el._xiMul[:,i,1][cut]
				yer = el._xiMul[:,i,2][cut]
				coef    = numpy.power(xxx,x_power)
				#plt.errorbar(xxx, coef*yyy, yerr=coef*yer,color='red', markersize=10,linewidth=2, label=r'$Mean \, simu \, \xi_{'+str(i)+'}$')
				plt.errorbar(xxx, coef*yyy,color='red', markersize=10,linewidth=2, label=r'$Mean \, simu \, \xi_{'+str(i)+'}$')
	
		plt.title(r'$'+self._title+'$', fontsize=40)
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
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
	def plot_cov_cor_matrix_different_method(self, realisation_type, dim='1D', path=None):

		realisation_type = numpy.array(realisation_type)

		path_to_load = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_'
		real = [ numpy.load(path_to_load + el + '_list_'+dim+'.npy') for el in realisation_type ]
		cov  = [ numpy.load(path_to_load + el + '_cov_'+dim+'.npy') for el in realisation_type ]

		if (path is not None):
			for el in path:
				cov  += [ numpy.load(el[1]) ]
				realisation_type = numpy.append( realisation_type, el[0] )

		### Plot the realisation
		for i in numpy.arange( len(real) ):
			print realisation_type[i]
			for j in numpy.arange(real[i][0,:].size):
				plt.errorbar(numpy.arange(real[i][:,j].size), real[i][:,j],color='blue',alpha=0.1)
			plt.errorbar(numpy.arange(real[i][:,j].size), numpy.mean(real[i],axis=1),marker='o',color='red',label=r'$Mean$')
			plt.xlabel(r'$bin \, index$', fontsize=40)
			plt.ylabel(r'$\xi(|s|)$', fontsize=40)
			plt.title(r'$'+realisation_type[i]+'$', fontsize=40)
			myTools.deal_with_plot(False,False,True)
			plt.xlim([ -1., cov[i][0,:].size+1 ])
			plt.show()

		### Plot diagonal
		for i in numpy.arange(len(cov)):
			plt.errorbar(numpy.arange(cov[i][0,:].size), numpy.diag(cov[i]),label=r'$'+realisation_type[i]+'$')
		plt.xlabel(r'$bin \, index$', fontsize=40)
		plt.ylabel(r'$Var(|s|)$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.xlim([ -1., cov[i][0,:].size+1 ])
		plt.show()

		### Plot Matrix
		for i in numpy.arange( len(real) ):
			myTools.plot2D(myTools.getCorrelationMatrix(cov[i]), None,None,None,None,realisation_type[i])

		myTools.plotCovar(cov,realisation_type)

		return
	def plot_distortion_matrix(self, dim='1D'):
	
		path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_'+dim+'_'+ self._middlefix + '.txt'
		matrix = numpy.loadtxt(path)
		myTools.plot2D(matrix)
		myTools.plotCovar([matrix],['Distortion matrix'])

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
			z_axis = [ '<s_{\perp}> - mean \, [\%]','<s_{\parallel}> - mean \, [\%]','<z_{pairs}> - mean \, [\%]' ]

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
		plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
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
			plt.xlabel(r'$x \, [h^{-1.}.Mpc]$')
			plt.ylabel(r'$y \, [h^{-1.}.Mpc]$')
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
			plt.xlabel(r'$x \, [h^{-1.}.Mpc]$')
			plt.ylabel(r'$y \, [h^{-1.}.Mpc]$')
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


'''

dic_CAMB = {
	'mulpol_index' : 0,
	'start_fit'   : 20.,
	'end_fit'     : 60.,
	'b' : -1.,
	'roof' : 0.,
	'fix_roof_nul' : True,
	'guess_b' : False,
	'min_for_guess' : 20.,
	'max_for_guess' : 50.,
	}
dic_class = {
	'minXi': 0.,
	'maxXi': 200.,
	'nbBin': 50,
	'nbBinM': 25,
	'nb_Sub_Sampling': 80,
	'size_bin_calcul_s': 1.,
	'size_bin_calcul_m': 0.02,
	'correlation': 'q_f',
	'path_to_txt_file_folder': 'NOTHING',
	'f1': 'LYA',
	'f2': 'a',
	'q1': 'QSO',
	'q2': 'a',
	'name' : 'Mean \, simu'
	}
dic_simu = {
	'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/',
	'nb_box' : 1,
	'nb_simu' : 10,
	'projected' : False
}
i = sys.argv[1]
j = sys.argv[2]	
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
corr = Correlation3D(dic_class)
corr.save_list_realisation('subsampling', 80)
correlation_matrix_path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy'
corr.send_BAOFIT('subsampling',correlation_matrix_path,True,False,False,False)
'''

















