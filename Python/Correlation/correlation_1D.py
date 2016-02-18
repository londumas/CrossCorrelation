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

		cut = 0

		xxx = self._xi[cut:,0]
		yyy = self._xi[cut:,1]
		yer = self._xi[cut:,2]
		plt.errorbar(xxx, yyy, yerr=yer, marker='o', label=r'$'+self._name+'$', markersize=8,linewidth=2)

		xMin = numpy.amin(xxx)
		xMax = numpy.amax(xxx)
		yMin = numpy.amin(yyy)
		yMax = numpy.amax(yyy)

		for el in other:
			if (el._correlation == 'f_f2_lRF_devide'): 
				TMP_xxx = 1./el._xi[cut:,0]
			else:
				TMP_xxx = el._xi[cut:,0]
			TMP_yyy = el._xi[cut:,1]
			TMP_yer = el._xi[cut:,2]
			plt.errorbar(TMP_xxx, TMP_yyy, yerr=TMP_yer, marker='o', label=r'$'+el._name+'$', markersize=8,linewidth=2)

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
	'f2': 'LYA',
	'nb_Sub_Sampling': 80,
	'name' : 'NAME'
}



### Data
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
dic_class['name'] = "data"
dic_class['correlation'] = "f_f_r"
dic_class['f1'] = "LYA"
corrD = Correlation_1D(dic_class)
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
dic_class['name'] = "dataM1"
dic_class['correlation'] = "f_f_r"
dic_class['f1'] = "LYA"
corrDD = Correlation_1D(dic_class)

corrD.plot(False,False,[corrDD])
































































