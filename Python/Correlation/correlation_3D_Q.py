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
import sys

### Perso lib
import correlation_3D
import myTools
import const

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
		if (not self._load_txt):
			self._xiMul = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_Multipol.npy')

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
		cut = (xiMu[:,:,2]!=0.)
		xiMu[:,:,0][cut] /= xiMu[:,:,2][cut]
		xiMu[:,:,1][cut] /= xiMu[:,:,2][cut]
		xiMu[:,:,2][cut] /= coef
		xiMu[:,:,3][cut]  = numpy.sqrt(xiMu[:,:,2][cut])

		### we
		cut = (xiWe[:,:,1]!=0.)
		xiWe[:,:,0][cut] /= xiWe[:,:,1][cut]
		xiWe[:,:,1][cut] /= coef
		xiWe[:,:,2][cut]  = numpy.sqrt(xiWe[:,:,1][cut])
		### 1D
		cut = (xi1D[:,1]!=0.)
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

			xiMu  = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_Mu.npy')
			xiWe  = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_We.npy')
			xi1D  = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_1D.npy')
			xi2D  = numpy.load(self._path_to_txt_file_folder+'xi_QSO_QSO_result_2D.npy')

		return xiMu, xiWe, xi1D, xi2D
	def save_list_realisation_simulation(self, dic_class, dic_Q, dic_simu):

		nb_realisation = dic_simu['nb_box']*dic_simu['nb_simu']
		pathToSave = dic_simu['path_to_simu'] + 'Results'+dic_simu['prefix']+'/' + self._prefix + '_result_'

		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,nb_realisation) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,3,nb_realisation) )
		list1D       = numpy.zeros( shape=(self._nbBin1D,nb_realisation) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,nb_realisation) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,self._nb_multipole_max,nb_realisation) )
		listGrid     = numpy.zeros( shape=(self._nbBin2D,3,nb_realisation) )
		list_mean_z  = numpy.zeros( shape=(nb_realisation) )

		nb = 0
		for i in range(0,dic_simu['nb_box']):
			for j in range(0,dic_simu['nb_simu']):

				print i, j

				raw = dic_simu['path_to_simu'] + 'Box_00' + str(i) + '/Simu_00' + str(j) +'/'
				dic_class['path_to_txt_file_folder'] = raw+'Results'+dic_simu['prefix']+'/'
				dic_Q['path_to_cat']                 = raw+'Data/'+dic_Q['sufix']+'.fits'

				corr = Correlation3DQ(dic_class,dic_Q)
				list1D[:,nb]         = corr._xi1D[:,1]
				list2D[:,nb]         = corr._xi2D[:,:,1].flatten()
				listMu[:,nb]         = corr._xiMu[:,:,2].flatten()
				listWe[:,:,nb]       = corr._xiWe[:,:,1]
				listMultipol[:,:,nb] = corr.get_multipol(corr._xiMu)[:,:,1]

				listGrid[:,0,nb] = corr._xi2D_grid[:,:,0].flatten()
				listGrid[:,1,nb] = corr._xi2D_grid[:,:,1].flatten()
				listGrid[:,2,nb] = corr._xi2D_grid[:,:,2].flatten()
				list_mean_z[nb]  = corr._meanZ

				nb += 1

		listMu       = listMu[:,:nb]
		listWe       = listWe[:,:,:nb]
		list1D       = list1D[:,:nb]
		list2D       = list2D[:,:nb]
		listMultipol = listMultipol[:,:,:nb]
		listGrid     = listGrid[:,:,:nb]
		list_mean_z  = list_mean_z[:nb]

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

		return
	def write_BAOFIT_ini_and_grid(self, param, dist_matrix=False, with_metals_templates=False, prefix2=''):

		precision = 1000.
		param = param.astype('str')

		if (self._correlation=='q_q'):
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._q1
		elif (self._correlation=='q_f'):
			path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1

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

### For data
#fiducial =  DR9LyaMocksLCDM
#nowiggles = DR9LyaMocksLCDMSB

### For Mocks
fiducial =  DR9LyaMocks
nowiggles = DR9LyaMocksSB

omega-matter    = 0.27
hubble-constant = 0.7
sigma8          = 0.794961
zref            = 2.25

## k-space fit
kspace = true
ell-max = 4

# Model configuration
anisotropic = yes
decoupled   = yes
custom-grid = yes
#combined-bias = yes
pixelize    = yes

# Parameter setup

model-config = value[beta]=               """+param[0] +""";
model-config = value[(1+beta)*bias]=      4.;
model-config = fix[gamma-bias]=           """+param[2] +""";
model-config = fix[gamma-beta]=           """+param[3] +""";
model-config = fix[SigmaNL-perp]=         """+param[7] +""";
model-config = fix[1+f]=                  """+param[8] +""";
model-config = fix[BAO amplitude]=        """+param[9] +""";
model-config = fix[BAO alpha-iso]=        """+param[10] +""";
model-config = value[BAO alpha-parallel]= """+param[11] +""";
model-config = value[BAO alpha-perp]=     """+param[12]+""";
model-config = fix[gamma-scale]=          """+param[13]+""";
model-config = fix[pixel scale]=3.15;
#model-config = value[beta*bias]=            """+param[1] +""";

# Broadband distortion model
#dist-add = rP,rT=0:2,-3:1
#dist-add = -2:0,0:4:2,0

## 2D chisq scan in BAO parameters
model-config = binning[BAO alpha-parallel] ={0.98:1.02}*50
model-config = binning[BAO alpha-perp]     ={0.98:1.02}*50	

## Maximum allowed radial dilation (increases the range that model needs to cover)
dilmin = 0.2
dilmax = 20.

# boxprior keeps result positive (since model only depends on squared value)
#model-config = boxprior[SigmaNL-perp] @ (0,6);
# un-comment next line to broaden all scales (default is peak only)
#nl-broadband = true



### Data Options #############################################

## Data to analyze
data             = """+path_to_data+"""

## Data format
data-format = comoving-cartesian
axis1-bins = ["""+ str(self._minY2D) + ':' + str(self._maxY2D) + """]*"""+ str(self._nbBinY2D) +"""
axis2-bins = ["""+ str(self._minX2D) + ':' + str(self._maxX2D) + """]*"""+ str(self._nbBinX2D) +"""
axis3-bins = {"""+ str(self._meanZ) +"""}

### Analysis Options #########################################

# Cuts to apply before fitting
rmin = 40
rmax = 180

# Generate a second set of outputs with the additive distortion turned off
#alt-config = fix[dist*]=0

# Do not dump multipoles (since the distortion model multipole integrals are singular)
zdump = """+ str(self._meanZ) +"""
ndump = 1000

# Prefix to use for all analysis output files
output-prefix = """ + path_to_BAOFIT + """.
"""

		text_file = open(path_to_BAOFIT+'.ini', "w")
		text_file.write(string_ini)
		text_file.close()

		return









