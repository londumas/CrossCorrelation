# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_1D.py
#

### Python lib
import subprocess
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit
import copy

### Perso lib
import myTools
import const
import const_delta
import CAMB

raw_dic_simu = {
	'path_to_simu' : 'NOTHING',
	'nb_box' : 1,
	'nb_simu' : 1,
	'prefix' : ''
}

class Correlation_1D:
	
	def __init__(self, dic=None):
		"""
		correlationType:
			- 'f_f_r'
			- 'f_f_lRF'
			- 'f_f_lRF_devide' (default)
			- 'f_f2_r'
			- 'f_f2_lRF'
			- 'f_f2_lRF_devide'
		"""
		
		if (dic==None):
			dic = {
			'correlation': 'f_f_lRF_devide',
			'path_to_txt_file_folder': 'NOTHING',
			'f1': 'LYA',
			'f2': 'LYA',
			'nb_Sub_Sampling': 80,
			'name' : 'Data'
			}

		verbose__ = False	

		### folder wher data are
		self._path_to_txt_file_folder = copy.deepcopy(dic['path_to_txt_file_folder'])

		### Name of the correlation
		self._name = dic['name']

		### forest and QSO name
		self._f1 = dic['f1']
		self._f2 = dic['f2']

		### Get lines
		if   (self._f1=='LYB'):
			self._lines1 = const_delta.LYB_lines
			self._name_line1 = const_delta.LYB_lines_names
		elif (self._f1=='LYA'):
			self._lines1 = const_delta.LYA_lines
			self._name_line1 = const_delta.LYA_lines_names
		elif (self._f1=='SIIV'):
			self._lines1 = const_delta.SIIV_lines
			self._name_line1 = const_delta.SIIV_lines_names
		elif (self._f1=='CIV'):
			self._lines1 = const_delta.CIV_lines
			self._name_line1 = const_delta.CIV_lines_names
		elif (self._f1=='MGII'):
			self._lines1 = const_delta.MGII_lines
			self._name_line1 = const_delta.MGII_lines_names
	
		if   (self._f2=='LYB'):
			self._lines2 = const_delta.LYB_lines
			self._name_line2 = const_delta.LYB_lines_names
		elif (self._f2=='LYA'):
			self._lines2 = const_delta.LYA_lines
			self._name_line2 = const_delta.LYA_lines_names
		elif (self._f2=='SIIV'):
			self._lines2 = const_delta.SIIV_lines
			self._name_line2 = const_delta.SIIV_lines_names
		elif (self._f2=='CIV'):
			self._lines2 = const_delta.CIV_lines
			self._name_line2 = const_delta.CIV_lines_names
		elif (self._f2=='MGII'):
			self._lines2 = const_delta.MGII_lines
			self._name_line2 = const_delta.MGII_lines_names
		
		### Correlation type
		self._correlation = dic['correlation']
		if (self._correlation=='f_f_r'):
			self._xTitle = '|s| \, [h^{-1}.Mpc]'
			self._yTitle = '\\xi (|s|)'
			self._title = '\delta_{'+self._f1+'}'
			self._prefix = 'xi_1D_delta_delta'
			self._middlefix = self._f1
			self._lines2 = self._lines1
			self._name_line2 = self._name_line1
		elif (self._correlation=='f_f_lRF'):
			self._xTitle = '\Delta \lambda_{R.F.} \, [\AA]'
			self._yTitle = '\\xi(\Delta \lambda_{R.F.})'
			self._title = '\delta_{'+self._f1+'}'
			self._prefix = 'xi_1DlRF_delta_delta'
			self._middlefix = self._f1
			self._lines2 = self._lines1
			self._name_line2 = self._name_line1
		elif (self._correlation=='f_f_lRF_devide'):
			self._xTitle = '\lambda_{1}/\lambda_{2}'
			self._yTitle = '\\xi(\lambda_{1}/\lambda_{2})'
			self._title = '\delta_{'+self._f1+'}'
			self._prefix = 'xi_1DlRFDevide_delta_delta'
			self._middlefix = self._f1
			self._lines2 = self._lines1
			self._name_line2 = self._name_line1
		elif (self._correlation=='f_f2_r'):
			self._xTitle = '|s| \, [h^{-1}.Mpc]'
			self._yTitle = '\\xi (|s|)'
			self._title = '\delta_{'+self._f1+'} \, - \, \delta_{'+self._f2+'}'
			self._prefix = 'xi_1D_delta_delta2'
			self._middlefix = self._f1 + '_' + self._f2
		elif (self._correlation=='f_f2_lRF'):
			self._xTitle = '\Delta \lambda_{R.F.} \, [\AA]'
			self._yTitle = '\\xi(\Delta \lambda_{R.F.})'
			self._title = '\delta_{'+self._f1+'} \, - \, \delta_{'+self._f2+'}'
			self._prefix = 'xi_1DlRF_delta_delta2'
			self._middlefix = self._f1 + '_' + self._f2
		elif (self._correlation=='f_f2_lRF_devide'):
			self._xTitle = '\lambda_{1}/\lambda_{2}'
			self._yTitle = '\\xi(\lambda_{1}/\lambda_{2})'
			self._title = '\delta_{'+self._f1+'} \, - \, \delta_{'+self._f2+'}'
			self._prefix = 'xi_1DlRFDevide_delta_delta2'
			self._middlefix = self._f1 + '_' + self._f2

		path = self._path_to_txt_file_folder + self._prefix + '_' + self._middlefix + '.txt'
		if (verbose__): print '  Correlation_1D::__init__::  path = ', path

		### Set attributes set after
		self._grid  = None
		self._meanZ = None

		### Get data
		self._xi = self.fill_data(path, True)

		self._nbBin = self._xi[:,0].size

		return

	def fill_data(self, path, init=False):

		verbose__ = False	

		data = numpy.loadtxt(path)
		cut = (data[:,5]>1.)

		nbBin = data[:,2][ cut ].size
		if (nbBin==0):
			print "  Correlation_1D::fill_data::  WARNING:  file is empty: path = ", path
			return numpy.zeros(shape=(1,3))

		xi = numpy.zeros(shape=(nbBin,3))

		xi[:,0] = data[:,2][cut]/data[:,4][cut]
		xi[:,1] = data[:,0][cut]/data[:,4][cut]
		xi[:,2] = numpy.sqrt( (data[:,1][cut]/data[:,4][cut] -xi[:,1]**2. )/data[:,5][cut]  )
		

		if (init):
			self._nbBin = nbBin
			self._meanZ = numpy.sum(data[:,3][cut])/numpy.sum(data[:,4][cut])
			self._grid           = numpy.zeros( shape=(self._nbBin,2) )
			self._grid[:,0][cut] = data[:,2][cut]/data[:,4][cut]
			self._grid[:,1][cut] = data[:,3][cut]/data[:,4][cut]

			if (verbose__): 
				print "  ||              ||                ||              ||"
				print "  || nb pairs     ||  sum weight    ||  <z>         ||"
				print "  ||              ||                ||              ||"
				print "  ||  %1.4e  ||  %1.4e    ||  %1.4e  ||" % (numpy.sum(data[:,5][cut]), numpy.sum(data[:,4][cut]), self._meanZ)
				print "  ||              ||                ||              ||"

		return xi
	def save_list_realisation_simulation(self, dic_class, dic_simu):

		pathToSave = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix'] + '/'
		pathToSave += self._prefix + '_' + self._middlefix + '_result_'

		nb_realisation = dic_simu['nb_box']*dic_simu['nb_simu']
		list1D         = numpy.zeros( shape=(self._nbBin,nb_realisation) )
		listGrid       = numpy.zeros( shape=(self._nbBin,2,nb_realisation) )
		list_mean_z    = numpy.zeros( shape=(nb_realisation) )

		nb = 0
		for i in range(0,dic_simu['nb_box']):
			for j in range(0,dic_simu['nb_simu']):

				raw = dic_simu['path_to_simu'] + 'Box_00' + str(i) + '/Simu_00' + str(j) +'/'
				dic_class['path_to_txt_file_folder'] = raw + 'Results' + dic_simu['prefix'] + '/'

				try:
					corr = Correlation_1D(dic_class)
					list1D[:,nb]     = corr._xi[:,1]
					listGrid[:,0,nb] = corr._grid[:,0]
					listGrid[:,1,nb] = corr._grid[:,1]
					list_mean_z[nb]  = corr._meanZ
					nb += 1
				except:
					print '  ERROR: ', i, j
		
		list1D      = list1D[:,:nb]
		listGrid    = listGrid[:,:,:nb]
		list_mean_z = list_mean_z[:nb]

		numpy.save(pathToSave+'list_1D',list1D)
		numpy.save(pathToSave+'list_Grid',listGrid)
		numpy.save(pathToSave+'list_MeanZ',list_mean_z)

		return
	def set_values_on_mean_simulation(self, dic_simu):

		path_to_load = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix'] + '/'
		path_to_load += self._prefix + '_' + self._middlefix + '_result_'
		print path_to_load

		list1 = numpy.load(path_to_load+'list_1D.npy')
		nb_realisation = list1[0,:].size
		self._xi[:,1] = numpy.mean(list1,axis=1)
		self._xi[:,2] = numpy.sqrt( numpy.diag(numpy.cov(list1))/nb_realisation )

		### Grid
		listGrid = numpy.load(path_to_load+'list_Grid.npy')
		self._grid[:,0] = numpy.mean(listGrid[:,0,:],axis=1)
		self._grid[:,1] = numpy.mean(listGrid[:,1,:],axis=1)
		self._meanZ = numpy.mean(numpy.load(path_to_load+'list_MeanZ.npy'))

		print '  nb_realisation = ', nb_realisation

		return
	def fit_CAMB(self,distortion=False,dic=None):

		### Get the data
		xxx = self._xi[:,0]
		yyy = self._xi[:,1]
		yer = self._xi[:,2]
		xMin = numpy.amin(xxx)
		xMax = numpy.amax(xxx)

		### Get CAMB
		camb = CAMB.CAMB(dic)._xi0
		xxx_Camb = copy.deepcopy(xxx)
		yyy_Camb = numpy.interp(xxx,camb[:,0],camb[:,1])

		if (distortion):
			path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_'+ self._middlefix + '.txt'
			matrix = numpy.loadtxt(path)
			xxx_Camb = numpy.append(xxx_Camb,numpy.zeros(2100-xxx.size) )
			yyy_Camb = numpy.append(yyy_Camb,numpy.zeros(2100-yyy_Camb.size) )
			yyy_Camb = numpy.dot(matrix,yyy_Camb)
			yyy_Camb = yyy_Camb[(xxx_Camb!=0.)]
			xxx_Camb = xxx_Camb[(xxx_Camb!=0.)]
		b = 0.15
		#b = 0.05

		'''
		### Chi^{2}
		def chi2(b):
			model = yyy_Camb*b
			return numpy.sum( numpy.power( (yyy-model)/yer ,2.) )
	
		### Init and perform the fit
		m = Minuit(chi2, b=1.,error_b=0.1,print_level=-1, errordef=0.01) 	
		m.migrad()

		b = m.values['b']
		#b = 0.01
		print b
		'''

		for i in numpy.arange(1):
			coef = numpy.power(xxx,i)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o',color='blue')

			coef = numpy.power(xxx_Camb,i)
			plt.errorbar(xxx_Camb,coef*b*yyy_Camb,color='red')

			plt.title(r'$'+self._title+'$', fontsize=40)
			plt.xlabel(r'$'+self._xTitle+'$', fontsize=40)
			plt.ylabel(r'$'+self._yTitle+'$', fontsize=40)

			if (self._correlation=='f_f_r' or self._correlation=='f_f2_r' or self._correlation=='f_f_lRF' or self._correlation=='f_f_lRF'):
				plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
			if (self._correlation=='f_f_lRF_devide' or self._correlation=='f_f2_lRF_devide'):
				plt.xlim([ 0.99*xMin, 1.01*xMax ])

			myTools.deal_with_plot(False,False,False)
			plt.show()

		return
	def plot(self, with_lines=False, verbose=False, other=[]):

		cut_l = 0
		cut_h = -1
		list_corr_to_plot = numpy.append( [self],other )

		xMin = numpy.amin(self._xi[cut_l:cut_h,0])
		xMax = numpy.amax(self._xi[cut_l:cut_h,0])
		yMin = numpy.amin(self._xi[cut_l:cut_h,1])
		yMax = numpy.amax(self._xi[cut_l:cut_h,1])

		for el in list_corr_to_plot:

			print el._name, el._xi[0,1]

			if (el._correlation == 'f_f2_lRF_devide'): 
				TMP_xxx = 1./el._xi[cut_l:cut_h,0]
			else:
				TMP_xxx = el._xi[cut_l:cut_h,0]
			TMP_yyy = el._xi[cut_l:cut_h,1]
			TMP_yer = el._xi[cut_l:cut_h,2]
			xxx = TMP_xxx

			plt.errorbar(TMP_xxx, TMP_yer, label=r'$'+el._name+'$', markersize=8,linewidth=2)

			xMin = min(xMin, numpy.amin(TMP_xxx) )
			xMax = max(xMax, numpy.amax(TMP_xxx) )
			yMin = min(yMin, numpy.amin(TMP_yyy) )
			yMax = max(yMax, numpy.amax(TMP_yyy) )

		if (with_lines):

			if (verbose): print ' ||  name_1 - name_2 || line || lambda_rf_1 || lamnda_rf_2 || '

			yLi = [yMin,yMax]
			nbLines1 = self._lines1.size

			for i in range(0,nbLines1):

				if (self._correlation=='f_f_r' or self._correlation=='f_f_lRF' or self._correlation=='f_f_lRF_devide'):
					nbLines2 = i
				else:
					nbLines2 = self._lines2.size

				for j in range(0,nbLines2):

					if ( self._lines1[i]!=1215.67 and self._lines2[j]!=1215.67 ): continue

					if (self._correlation=='f_f_r' or self._correlation=='f_f2_r'):
						line = numpy.abs( const_delta.find_dist_correlation_lines(self._meanZ,self._lines2[j], self._lines1[i]) )
						if (line>xMax): continue
					elif (self._correlation=='f_f_lRF' or self._correlation=='f_f2_lRF'):
						line = abs(self._lines1[i]-self._lines2[j])
						if (line==0. or line>xMax): continue
					elif (self._correlation=='f_f_lRF_devide' or self._correlation=='f_f2_lRF_devide'):
						line = self._lines2[j]/self._lines1[i]
						if (line<xMin or line>xMax): line = 1./line
						if (line<xMin or line>xMax): continue
					if (verbose): print ' ||  ', self._name_line1[i] ,' - ', self._name_line2[j], ' || ', line, ' || ', self._lines1[i], ' || ', self._lines2[j], ' || '

					xLi = [line,line]
					name = self._name_line1[i]+' - ' + self._name_line2[j]
					plt.plot(xLi,yLi,color='green')
					plt.text(line, yMax, name, rotation='vertical', fontsize=20)

		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$'+self._xTitle+'$', fontsize=40)
		plt.ylabel(r'$'+self._yTitle+'$', fontsize=40)

		if (self._correlation=='f_f_r' or self._correlation=='f_f2_r' or self._correlation=='f_f_lRF' or self._correlation=='f_f_lRF'):
			plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
		if (self._correlation=='f_f_lRF_devide' or self._correlation=='f_f2_lRF_devide'):
			plt.xlim([ 0.99*xMin, 1.01*xMax ])

		myTools.deal_with_plot(False,False,True)
		plt.show()

		return
	def plot_distortion_matrix(self):
	
		path = self._path_to_txt_file_folder + self._prefix + '_distortionMatrix_'+ self._middlefix + '.txt'
		matrix = numpy.loadtxt(path)
		
		### Blind the zero terms
		matrix[ (matrix==0.) ] = numpy.float('nan')

		myTools.plot2D(matrix)
		#myTools.plotCovar([matrix],['Distortion matrix'])

		return



"""
correlationType:
	- 'f_f_r'
	- 'f_f_lRF'
	- 'f_f_lRF_devide' (default)
	- 'f_f2_r'
	- 'f_f2_lRF'
	- 'f_f2_lRF_devide'
"""
dic_class = {
	'correlation': 'f_f_lRF_devide',
	'path_to_txt_file_folder': 'NONE',
	'f1': 'LYA',
	'f2': 'a',
	'nb_Sub_Sampling': 80,
	'name' : 'NAME'
}
dic_simu = {
	'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
	'nb_box' : 10,
	'nb_simu' : 10,
	'prefix' : ''
}
dic_CAMB_corr = {
	'z' : 0.,
	'h_0' : 0.70,
	'omega_matter_0' : 0.27,
	'omega_lambda_0' : 0.73,
	'source'         : 'CAMB_Mocks_me'
}

list_corr = []


### Data
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_primery_sky/'
dic_class['name'] = "Data"
dic_class['correlation'] = "f_f_lRF"
corr = Correlation_1D(dic_class)
list_corr += [corr]
#corr.plot_distortion_matrix()
#dic_CAMB_corr['z'] = corr._meanZ
#corr.fit_CAMB()
#corr.fit_CAMB(True,dic_CAMB_corr)

'''
### Data
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
dic_class['name'] = "Data"
dic_class['correlation'] = "f_f_lRF"
corr = Correlation_1D(dic_class)
list_corr += [corr]
#corr.plot_distortion_matrix()
#dic_CAMB_corr['z'] = corr._meanZ
#corr.fit_CAMB()
#corr.fit_CAMB(True,dic_CAMB_corr)
list_corr[0].plot(False,False,list_corr[1:])
'''

'''
#
dic_simu['path_to_simu'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/'
dic_class['f1'] = 'LYA'
dic_class['name'] = 'Mean \, Simulation'
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_000/Simu_000/Results/'
dic_class['correlation'] = "f_f_r"
corr = Correlation_1D(dic_class)
corr.save_list_realisation_simulation(dic_class, dic_simu)
corr.set_values_on_mean_simulation(dic_simu)
list_corr += [corr]
corr.plot_distortion_matrix()
dic_CAMB_corr['z'] = corr._meanZ
corr.fit_CAMB(True,dic_CAMB_corr)
'''

"""
#
dic_simu['prefix'] = '_no_metals'
dic_simu['path_to_simu'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/'
dic_class['f1'] = 'LYA'
dic_class['name'] = 'Mean \, Simulation \, no \, metals'
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_000/Simu_000/Results_no_metals/'
dic_class['correlation'] = "f_f_r"
corr = Correlation_1D(dic_class)
corr.save_list_realisation_simulation(dic_class, dic_simu)
corr.set_values_on_mean_simulation(dic_simu)
list_corr += [corr]



"""
#list_corr[0].plot(True,True)
list_corr[0].plot(False,False,list_corr[1:])






























































