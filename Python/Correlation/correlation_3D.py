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
import myTools
import const


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
	'start_fit'   : 20.,
	'end_fit'     : 70.,
	'b' : -1.,
	'roof' : 0.,
	'fix_roof_nul' : True,
	'guess_b' : False,
	'min_for_guess' : 20.,
	'max_for_guess' : 50.,
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
		
		if (dic is None):
			dic = raw_dic_class

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
			print "  Correlation_3D::__init__::  ERROR:  'correlation' is incorrect "
			return
		path1D = self._path_to_txt_file_folder + self._prefix + '_Mu_' + self._middlefix + '.txt'
		path2D = self._path_to_txt_file_folder + self._prefix + '_2D_' + self._middlefix + '.txt'

		print
		print '  path1D = ', path1D

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
			
			
	
		if (selection==0 or selection==2):
			### 2D
			data = numpy.loadtxt(path2D)
	
			save0 = data[:,0]
			save1 = data[:,1]
			save2 = data[:,2]
			save3 = data[:,3]
			save5 = data[:,5]
			save6 = data[:,6]
	
			tmp_save0  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save1  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save2  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
			tmp_save3  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
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
				tmp_save5[idX][idY] += save5[i]
				tmp_save6[idX][idY] += save6[i]

			cut = (tmp_save5!=0.)
			xi2D[:,:,0][cut] = numpy.sqrt( (tmp_save2[cut]/tmp_save5[cut])**2. + (tmp_save3[cut]/tmp_save5[cut])**2. )
			xi2D[:,:,1][cut] = tmp_save0[cut] / tmp_save5[cut]
			xi2D[:,:,2][cut] = numpy.sqrt( (tmp_save1[cut]/tmp_save5[cut] - xi2D[:,:,1][cut]*xi2D[:,:,1][cut])/tmp_save6[cut] )

			if (init):
				self._meanZ = numpy.sum(data[:,4])/numpy.sum(data[:,5])

		return xiMu, xiWe, xi1D, xi2D
	def read_grid(self):

		path = self._path_to_txt_file_folder + self._prefix + '_2D_' + self._middlefix + '.txt'
		data = numpy.loadtxt(path)
		int_binSize = int(self._binSize)
	
		tmp_save2  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		tmp_save3  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		tmp_save4  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		tmp_save5  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		meanRperp  = numpy.zeros( shape=(self._nbBinX2D,2))
		meanRparal = numpy.zeros( shape=(self._nbBinY2D,2))
	
		for i in range( 0,data[:,2].size ):
			iX = i/self._nbBinY2D_calcul
			iY = i%self._nbBinY2D_calcul
	
			idX = iX/int_binSize
			idY = iY/int_binSize
	
			tmp_save2[idX][idY] += data[i,2]
			tmp_save3[idX][idY] += data[i,3]
			tmp_save4[idX][idY] += data[i,4]
			tmp_save5[idX][idY] += data[i,5]

			meanRperp[idX][0]  += data[i,2]
			meanRperp[idX][1]  += data[i,5]
			meanRparal[idY][0] += data[i,3]
			meanRparal[idY][1] += data[i,5]
		
		### Get the grid
		cut = (tmp_save5!=0.)
		rPerp  = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		rParal = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		z      = numpy.zeros( shape=(self._nbBinX2D,self._nbBinY2D) )
		rPerp[cut]  = tmp_save2[cut]/tmp_save5[cut]
		rParal[cut] = tmp_save3[cut]/tmp_save5[cut]
		z[cut]      = tmp_save4[cut]/tmp_save5[cut]

		grid = numpy.zeros( shape=(self._nbBin2D,4) )
		indexMatrix = numpy.arange(self._nbBin2D)
		grid[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
		grid[:,1] = rParal.flatten()
		grid[:,2] = rPerp.flatten()
		grid[:,3] = z.flatten()
		
		### Get the s_perp bin center
		meanRperp[:,0]  /= meanRperp[:,1]
			
		### Get the s_paral bin center
		meanRparal[:,0] /= meanRparal[:,1]

		return grid, meanRperp[:,0], meanRparal[:,0]
	def save_list_realisation(self, realisation_type, nb_realisation):

		if ( realisation_type=='subsampling' and (nb_realisation != self.nb_Sub_Sampling) ):
			print "  Correlation_3D::saveListReal::  WARNING:  nb_realisation != self.nb_Sub_Sampling"

		list1D       = numpy.zeros( shape=(self._nbBin1D,nb_realisation) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,nb_realisation) )
		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,nb_realisation) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,3,nb_realisation) )
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
			#listMultipol[:,:,i] = plotMultipol(xiMu)[:,:,1]
			'''
			except:
				print '   ERROR:: ', path1D+str(i)+'.txt',path2D+str(i)+'.txt'
				commandProd = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Correlation/bin/main.exe"
				tmp_command = "clubatch \"time ; hostname ; " + commandProd + ' 5 0 ' + str(i) + ' 0 0 0 ' + "\""
				subprocess.call(tmp_command, shell=True)
			'''
		
		numpy.save(pathToSave+'list_1D',list1D)
		numpy.save(pathToSave+'list_2D',list2D)
		numpy.save(pathToSave+'list_Mu',listMu)
		numpy.save(pathToSave+'list_We',listWe)
		#numpy.save(pathToSave+'_list_Multipol',listMultipol)
	
		cov1D = numpy.cov(list1D)
		cov2D = numpy.cov(list2D)
		covMu = numpy.cov(listMu)
		if (realisation_type=='subsampling'):
			cov1D /= nb_realisation
			cov2D /= nb_realisation
			covMu /= nb_realisation
	
		numpy.save(pathToSave+'cov_1D',cov1D)
		numpy.save(pathToSave+'cov_2D',cov2D)
		numpy.save(pathToSave+'cov_Mu',covMu)
	
		return
	def set_error_on_covar_matrix(self, realisation_type):

		### mu
		#path = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_Mu.npy'
		#self._xiMu[:,:,3] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path)) ),self._nbBin1D,self._nbBinM)
		### we
		path = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_We.npy'
		print numpy.load(path)
		### 1D
		path = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_1D.npy'
		self._xi1D[:,2]   = numpy.sqrt( numpy.diag(numpy.load(path)) )
		### 2D
		path = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_2D.npy'
		self._xi2D[:,:,2] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path)) ),self._nbBinX2D, self._nbBinY2D)

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
				'xi4' : False
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

		cut = (result_xi[:,:,2]==0.1)
		result_xi[cut] = 0.

		return result_xi
	def fit_CAMB(self,dic=None):
		
		### Constants
		if (dic is None):
			dic = raw_dic_CAMB
		
		### Get the data
		cut = numpy.logical_and( (self._xi1D[:,0]>=dic['start_fit']),(self._xi1D[:,0]<dic['end_fit']) )
		xxx = self._xi1D[:,0][cut]
		yyy = self._xi1D[:,1][cut]
		yer = self._xi1D[:,2][cut]
	
		### Get Camb
		idx=1
		if (dic['mulpol_index']==2):
			idx = 2
		data_camb = numpy.loadtxt(const.pathToCamb__)
		yyy_Camb  = numpy.interp(xxx,data_camb[1:,0],data_camb[1:,idx])
	
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
		dic['b'] = m.values[ 'b' ]
		dic['roof'] = m.values[ 'roof' ]
		print '  b = ', dic['b']
		print '  roof = ', dic['roof']
	
		### Print chi^2
		if (dic['fix_roof_nul']):
			print '  DoF   = ', yyy.size, ' - 1'
		else:
			print '  DoF   = ', yyy.size, ' - 2'
		print '  chi^{2} = ', numpy.sum( numpy.power( (yyy - (yyy_Camb*dic['b']+dic['roof']) )/yer ,2.) )
		
		return dic
	def write_BAOFIT_ini_and_grid(self, param):

		precision = 1000.
		param = param.astype('str')
		path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_q_f__'+self._f1+'__'+self._q1+'/bao2D'

		grid, meanRperp, meanRparal = self.read_grid()

		### Saving .grid
		numpy.savetxt(path_to_BAOFIT+'.grid',zip(grid[:,0],grid[:,1],grid[:,2],grid[:,3]),fmt='%u %1.20e %1.20e %1.20e')
		del grid

		### Get the s_perp bin center
		meanRperp = (precision*meanRperp).astype(int)/precision
		stringRperp = ''
		for el in meanRperp[:-1]:
			stringRperp += str(el) + ','
		stringRperp += str( meanRperp[-1] )

		### Get the s_paral bin center
		meanRparal = (precision*meanRparal).astype(int)/precision
		stringRparal = ''
		for el in meanRparal[:-1]:
			stringRparal += str(el) + ','
		stringRparal += str( meanRparal[-1])

		### Get the redshift center
		meanZ = (precision*self._meanZ).astype(int)/precision
		stringRedshift = str(meanZ)

		string_ini = """

## Linear theory P(k) templates with and w/o wiggles
modelroot = """ + const.path_to_BAOFIT_model__ + """
fiducial =  DR9LyaMocksLCDM
nowiggles = DR9LyaMocksLCDMSB

## k-space fit
kspace = true
ell-max = 6

# Model configuration
cross-correlation = yes
anisotropic = yes
decoupled   = yes
custom-grid = yes
pixelize = yes
#dist-matrix = yes
dist-matrix-order = """+str(self._nbBin2D)+"""

# Parameter setup
model-config = value[beta]=               """+param[0] +""";
model-config = fix[(1+beta)*bias]=        """+param[1] +""";
model-config = fix[gamma-bias]=           """+param[2] +""";
model-config = fix[gamma-beta]=           """+param[3] +""";
model-config = value[delta-v]=            """+param[4] +""";
model-config = value[bias2]=              """+param[5] +""";
model-config = fix[beta2*bias2]=          """+param[6] +""";
model-config = fix[1+f]=                  """+param[7] +""";
model-config = fix[SigmaNL-perp]=         """+param[8] +""";
model-config = fix[BAO amplitude]=        """+param[9] +""";
model-config = fix[BAO alpha-iso]=        """+param[10] +""";
model-config = value[BAO alpha-parallel]= """+param[11] +""";
model-config = value[BAO alpha-perp]=     """+param[12]+""";
model-config = fix[gamma-scale]=          """+param[13]+""";
model-config = fix[pixel scale]=0.7;

## 2D chisq scan in BAO parameters
model-config = binning[BAO alpha-parallel] ={0.5:1.5}*100
model-config = binning[BAO alpha-perp]     ={0.5:1.5}*100	

## Reference redshift
zref = 2.3

## Maximum allowed radial dilation (increases the range that model needs to cover)
dilmin = 0.5
dilmax = 2.5

## Non-linear broadening with 1+f = (SigmaNL-par)/(SigmaNL-perp)

# boxprior keeps result positive (since model only depends on squared value)
model-config = boxprior[SigmaNL-perp] @ (0,6);
# un-comment next line to broaden all scales (default is peak only)
#nl-broadband = true

# Broadband distortion model
dist-add = rP,rT=0:2,-3:1

### Data Options #############################################

## Data to analyze
data = """+path_to_BAOFIT+"""
dist-matrix-name = """+path_to_BAOFIT+"""

## Data format
data-format = comoving-cartesian
axis1-bins = {""" + stringRparal   + """}
axis2-bins = {""" + stringRperp    + """}
axis3-bins = {""" + stringRedshift + """}

### Analysis Options #########################################

# Cuts to apply before fitting
rmin = 30
rmax = 180

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
	def send_BAOFIT(self, realisation_type, correlation_matrix_path=None, saving_data=True, scan=False, toyMC=0):
		"""

			Send the fit of the BAO with BAOFIT
			realisation_type:
			cor:=None

		"""

		### Constants
		param = numpy.asarray( [1.6,-0.336,0.9,0.,0.,3.25,0.962524,3.26,1.966,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.] )

		path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_q_f__'+self._f1+'__'+self._q1 + '/'
		subprocess.call('mkdir ' + path_to_BAOFIT, shell=True)
		path_to_BAOFIT += 'bao2D.'

		self.write_BAOFIT_ini_and_grid(param)

		if (saving_data):

			path_to_cov = self._path_to_txt_file_folder + self._prefix + '_'+ self._middlefix + '_' + realisation_type + '_cov_2D.npy'
			path_to_distortion_matrix = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_2D_'+ self._middlefix + '.txt'
			if (correlation_matrix_path is None): correlation_matrix_path = path_to_cov
			print '  correlation matrix path = ', correlation_matrix_path

			### For data
			cov = numpy.load(path_to_cov)
			cor = myTools.getCorrelationMatrix( numpy.load(correlation_matrix_path) )
	
			### Save .data
			correlation = numpy.zeros( shape=(self._nbBin2D,2) )
			indexMatrix = numpy.arange(self._nbBin2D)
			correlation[:,0] = (indexMatrix%self._nbBinY2D)*self._nbBinX2D + indexMatrix/self._nbBinY2D
			correlation[:,1] = self._xi2D[:,:,1].flatten()
			cutCorrelation = (correlation[:,1]!=0.)
			numpy.savetxt( path_to_BAOFIT + 'data',zip(correlation[:,0][cutCorrelation],correlation[:,1][cutCorrelation]),fmt='%u %1.20e')
			del correlation, indexMatrix, cutCorrelation
			
			### Save .cov
			print '  Saving .cov'
			cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))
			covarianceMatrix = numpy.zeros( shape=(self._nbBin2D*self._nbBin2D,3) )
			indexMatrix1 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)/self._nbBin2D
			indexMatrix2 = numpy.arange(self._nbBin2D*self._nbBin2D).reshape(self._nbBin2D,self._nbBin2D)%self._nbBin2D
			indexMatrix1 = (indexMatrix1%self._nbBinY2D)*self._nbBinX2D + indexMatrix1/self._nbBinY2D
			indexMatrix2 = (indexMatrix2%self._nbBinY2D)*self._nbBinX2D + indexMatrix2/self._nbBinY2D
			covarianceMatrix[:,0] = numpy.triu(indexMatrix1,k=0).flatten()
			covarianceMatrix[:,1] = numpy.triu(indexMatrix2,k=0).flatten()
			covarianceMatrix[:,2] = numpy.triu(cov,k=0).flatten()
			cutCovarMatrix = (covarianceMatrix[:,2]!=0.)
			if (covarianceMatrix[:,2][cutCovarMatrix].size != int(self._nbBin2D*(self._nbBin2D+1.)/2.) ):
				print '  xi_delta_QSO.py::prepareForBAOFIT()  size of covariance matrix is incorrect'
				print '  size covariance matrix', cutCovarMatrix[:,2][cutCovarMatrix].size
				print '  size it should have', self._nbBin2D*(self._nbBin2D+1.)/2.
				return
			numpy.savetxt( path_to_BAOFIT + 'cov',zip(covarianceMatrix[:,0][cutCovarMatrix],covarianceMatrix[:,1][cutCovarMatrix],covarianceMatrix[:,2][cutCovarMatrix]),fmt='%u %u %1.20e')
			del cov, covarianceMatrix, cutCovarMatrix

			'''
			### Save .dmat
			print '  Saving .dmat'
			dmatData = numpy.loadtxt(path_to_distortion_matrix)
			distortionMatrix = numpy.zeros( shape=(self._nbBin2D*self._nbBin2D,3) )
			distortionMatrix[:,0] = indexMatrix1.flatten()
			distortionMatrix[:,1] = indexMatrix2.flatten()
			distortionMatrix[:,2] = dmatData.flatten()
			cutDistortionMatrix = (distortionMatrix[:,2]!=0.)
			numpy.savetxt( path_to_BAOFIT + 'dmat',zip(distortionMatrix[:,0][cutDistortionMatrix], distortionMatrix[:,1][cutDistortionMatrix], distortionMatrix[:,2][cutDistortionMatrix]),fmt='%u %u %1.20e')
			del dmatData, distortionMatrix, cutDistortionMatrix
			'''

		### Send the fit (#--parameter-scan'   ### --toymc-samples 10000)
		command = const.path_to_BAOFIT_bin__ + ' -i ' + path_to_BAOFIT + 'ini'
		if (scan):
			command += ' --parameter-scan'
		if (toyMC!=0):
			command += ' --toymc-samples ' + str(toyMC)
		print command
		subprocess.call(command, shell=True)

		return
	def plot_1d(self, x_power=0, other=[], title=True):

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
		plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$'+self._name+'$')

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
			plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, fmt='o', label=r'$'+el._name+'$')

		if (title): plt.title(r'$'+self._title+'$', fontsize=40)
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
	def plot_we(self, x_power=0, title=True):
	
		label = ['0.8 < |\mu|', '0.5 < |\mu| \leq 0.8', '|\mu| \leq 0.5']
		
		
		for i in numpy.arange(0,3):
		
			cut = (self._xiWe[:,i,2]>0.)
			if (self._xiWe[:,i,0][cut].size==0):
				continue
		
			xxx = self._xiWe[:,i,0][cut]
			yyy = self._xiWe[:,i,1][cut]
			yer = self._xiWe[:,i,2][cut]
			coef = numpy.power(xxx,x_power)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$'+label[i]+'$')
		
			if (x_power==0):
				plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
				plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
			if (x_power==1):
				plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
				plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
			if (x_power==2):
				plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
				plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=2)
		
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
	
		#plt.plot( [0.,200.],[0.,4*200.],color='white',linewidth=2 )
		#plt.plot( [0.,200.],[0.,-4*200.],color='white',linewidth=2 )
		#plt.plot( [0.,200.],[0.,200.],color='white',linewidth=2 )
		#plt.plot( [0.,200.],[0.,-200.],color='white',linewidth=2 )
	
		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
	
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
		plt.show()
	
		return
	def plot_CAMB(self, dic=None, x_power=0):
		
		if (dic is None):
			dic = {
				'mulpol_index' : 0,
				'start_fit'   : 20.,
				'end_fit'     : 70.,
				'b' : -1.,
				'roof' : 0.,
				'fix_roof_nul' : True,
				'guess_b' : False,
				'min_for_guess' : 20.,
				'max_for_guess' : 50.
			}

		### Get the data
		xxx = self._xi1D[:,0]
		yyy = self._xi1D[:,1]
		yer = self._xi1D[:,2]
	
		### Get the smooth CAMB
		idx=1
		if (dic['mulpol_index']==2):
			idx = 2
		data_camb = numpy.loadtxt(const.pathToCamb__)
		yyy_Camb  = numpy.interp(xxx,data_camb[1:,0],data_camb[1:,idx])
		cut = (data_camb[1:,0] <= numpy.amax(xxx)*1.1)
		size = cut[cut].size
		result_1D_camb = numpy.zeros( shape=(size,3) )
		result_1D_camb[:,0] = data_camb[1:,0][cut]
		result_1D_camb[:,1] = dic['b']*data_camb[1:,idx][cut]+dic['roof']
		result_1D_camb[:,2] = 0.0000000001
	
		### Show result
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,yerr=coef*yer,fmt='o')
	
		coef2 = numpy.power(result_1D_camb[:,0],x_power)
		plt.plot(result_1D_camb[:,0],coef2*result_1D_camb[:,1],color='red')
	
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|^{1}.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
			
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		plt.ylim([ min( numpy.amin(coef*yyy)*1.1,numpy.amin(coef*yyy)*0.9) , max( numpy.amax(coef*yyy)*1.1,numpy.amax(coef*yyy)*0.9) ])
		myTools.deal_with_plot(False,False,False)
		plt.show()

		return
	def plot_multipol(self, x_power=0):

		if (self._xiMul is None):
			self._xiMul = self.get_multipol(self._xiMu)
	
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
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$\xi_{'+str(i)+'}$',color=color[i])
	
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
				plt.errorbar(numpy.arange(real[i][:,j].size), real[i][:,j],fmt='o',color='blue',alpha=0.1)
			plt.errorbar(numpy.arange(real[i][:,j].size), numpy.mean(real[i],axis=1),fmt='o',color='red',label=r'$Mean$')
			plt.xlabel(r'$bin \, index$', fontsize=40)
			plt.ylabel(r'$\xi(|s|)$', fontsize=40)
			plt.title(r'$'+realisation_type[i]+'$', fontsize=40)
			myTools.deal_with_plot(False,False,True)
			plt.xlim([ -1., cov[i][0,:].size+1 ])
			plt.show()

		### Plot diagonal
		for i in numpy.arange(len(cov)):
			plt.errorbar(numpy.arange(cov[i][0,:].size), numpy.diag(cov[i]),fmt='o',label=r'$'+realisation_type[i]+'$')
		plt.xlabel(r'$bin \, index$', fontsize=40)
		plt.ylabel(r'$Var(|s|)$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.xlim([ -1., cov[i][0,:].size+1 ])
		plt.show()

		myTools.plotCovar(cov,realisation_type)

		return
	def plot_distortion_matrix(self, dim='1D'):
	
		path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_'+dim+'_'+ self._middlefix + '.txt'
		matrix = numpy.loadtxt(path)
		myTools.plot2D(matrix)
		myTools.plotCovar([matrix],['Distortion matrix'])

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









