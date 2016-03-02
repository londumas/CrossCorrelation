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
	'projected' : False
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
		print '  Correlation_1D::__init__::  path = ', path

		### Get data
		self._xi = self.fill_data(path, True)

		self._nbBin = self._xi[:,0].size


		
		return

	def fill_data(self, path, init=False):

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

			print "  ||              ||                ||              ||"
			print "  || nb pairs     ||  sum weight    ||  <z>         ||"
			print "  ||              ||                ||              ||"
			print "  ||  %1.4e  ||  %1.4e    ||  %1.4e  ||" % (numpy.sum(data[:,5][cut]), numpy.sum(data[:,4][cut]), self._meanZ)
			print "  ||              ||                ||              ||"

		return xi
	def save_list_realisation_simulation(self, dic_class, dic_simu):

		pathToSave = dic_simu['path_to_simu']
		if (dic_simu['projected']):
			pathToSave += 'Results_nicolasEstimator/'
		else:
			pathToSave += 'Results/'
		pathToSave += self._prefix + '_' + self._middlefix + '_result_list'

		nb_realisation = dic_simu['nb_box']*dic_simu['nb_simu']
		list1          = numpy.zeros( shape=(self._nbBin,nb_realisation) )

		for i in range(0,dic_simu['nb_box']):
			for j in range(0,dic_simu['nb_simu']):

				raw = dic_simu['path_to_simu'] + 'Box_00' + str(i) + '/Simu_00' + str(j) +'/'
				if (dic_simu['projected']):
					dic_class['path_to_txt_file_folder'] = raw+'Results_nicolasEstimator/'
				else:
					dic_class['path_to_txt_file_folder'] = raw+'Results/'

				corr = Correlation_1D(dic_class)
				list1[:,i*10+j] = corr._xi[:,1]

		numpy.save(pathToSave,list1)

		return
	def set_values_on_mean_simulation(self, dic_simu):

		path_to_load = dic_simu['path_to_simu']
		if (dic_simu['projected']):
			path_to_load += 'Results_nicolasEstimator/'
		else:
			path_to_load += 'Results/'
		path_to_load += self._prefix + '_' + self._middlefix + '_result_list.npy'

		nb_realisation = dic_simu['nb_box']*dic_simu['nb_simu']

		list1 = numpy.load(path_to_load)
		self._xi[:,1] = numpy.mean(list1,axis=1)
		self._xi[:,2] = numpy.sqrt( numpy.diag(numpy.cov(list1))/nb_realisation )

		return
	def fit_CAMB(self):

		### Get the data
		cut = (self._xi[:,0]>10.)
		xxx = self._xi[:,0][cut]
		yyy = self._xi[:,1][cut]
		yer = self._xi[:,2][cut]

		### Get CAMB
		camb = CAMB.CAMB()._xi0
		yyy_Camb = numpy.interp(xxx,camb[:,0],camb[:,1])

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

		for i in [0,1,2]:
			coef = numpy.power(xxx,i)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o')

			coef = numpy.power(camb[:,0],i)
			plt.errorbar(camb[:,0],coef*b*camb[:,1])
			myTools.deal_with_plot(False,False,False)
			plt.xlim([ 0., 200. ])

			plt.show()

		return
	def plot(self, with_lines=False, verbose=False, other=[]):

		cut_l = 0 #7
		cut_h = -1 #300

		if (self._correlation == 'f_f2_lRF_devide'): 
			xxx = 1./self._xi[cut_l:cut_h,0]
		else:
			xxx = self._xi[cut_l:cut_h,0]
		yyy = self._xi[cut_l:cut_h,1]
		yer = self._xi[cut_l:cut_h,2]
		plt.errorbar(xxx, yyy, yerr=yer, marker='o', label=r'$'+self._name+'$', markersize=8,linewidth=2)

		xMin = numpy.amin(xxx)
		xMax = numpy.amax(xxx)
		yMin = numpy.amin(yyy)
		yMax = numpy.amax(yyy)

		for el in other:
			if (el._correlation == 'f_f2_lRF_devide'): 
				TMP_xxx = 1./el._xi[cut_l:cut_h,0]
			else:
				TMP_xxx = el._xi[cut_l:cut_h,0]
			TMP_yyy = el._xi[cut_l:cut_h,1]
			TMP_yer = el._xi[cut_l:cut_h,2]
			plt.errorbar(TMP_xxx, TMP_yyy, yerr=TMP_yer, marker='o', label=r'$'+el._name+'$', markersize=8,linewidth=2,color='red')

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

					#if ( self._lines1[i]!=1215.67 and self._lines2[j]!=1215.67 ): continue
					#if ( self._lines1[i]==1238.821 or self._lines2[j]==1238.821 ): continue
					#if ( self._name_line1[i][:4]!='SiII' and self._name_line2[j][:4]!='SiII' ): continue
					#if ( self._name_line1[i]=='SiII(1304)' or self._name_line2[j]=='SiII(1304)' ): continue

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
	'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/',
	'nb_box' : 1,
	'nb_simu' : 10,
	'projected' : True
}


list_corr = []

### Data
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
dic_class['name'] = "Data"
dic_class['correlation'] = "f_f_lRF_devide"
corr = Correlation_1D(dic_class)
list_corr += [corr]
#
dic_class['name'] = 'Mean \, simu'
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_000/Simu_000/Results_nicolasEstimator/'
dic_class['correlation'] = "f_f_lRF_devide"
corr = Correlation_1D(dic_class)
corr.save_list_realisation_simulation(dic_class,dic_simu)
corr.set_values_on_mean_simulation(dic_simu)
list_corr += [corr]

list_corr[0].plot(False,False,list_corr[1:])
































































