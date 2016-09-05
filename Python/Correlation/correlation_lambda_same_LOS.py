# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#

### Python lib
import subprocess
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
from iminuit import Minuit

### Perso lib
import myTools
import const
import const_delta




# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#

### Python lib
import subprocess
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
from iminuit import Minuit

### Perso lib
import myTools
import const
import const_delta

raw_dic_class = {
	'min_l1_over_l2': 0.40,
	'max_l1_over_l2': 2.30,
	'min_theta': 0.,
	'max_theta': 0.003,
	'min_visual_z': -0.04,
	'max_visual_z': 0.04,
	'size_bin_calcul_l1_overl2': 1.e-03,
	'size_bin_calcul_theta': 1.e-04,
	'correlation': 'q_f',
	'path_to_txt_file_folder': 'NOTHING',
	'f1': 'LYA',
	'f2': 'LYA',
	'q': 'QSO',
	'name' : 'Data'
}

class CorrelationLambda:
	
	def __init__(self, dic=None):

		if (dic is None):
			dic = raw_dic_class

		### l1_over_l2
		self._minL1oL2 = dic['min_l1_over_l2']
		self._maxL1oL2 = dic['max_l1_over_l2']
		self._nbBinL1oL2 = int( (self._maxL1oL2-self._minL1oL2)/dic['size_bin_calcul_l1_overl2'] )

		### forest and QSO name
		self._f1 = dic['f1']
		self._q1 = dic['q']

		### Get lines
		if   (self._f1=='LYB'):
			self._line_RF = const_delta.lambda_RF_line_LYB
			self._lines1 = const_delta.LYB_lines
			self._name_line1 = const_delta.LYB_lines_names
		elif (self._f1=='LYA'):
			self._line_RF = const_delta.lambda_RF_line_LYA
			self._lines1 = const_delta.LYA_lines
			self._name_line1 = const_delta.LYA_lines_names
		elif (self._f1=='SIIV'):
			self._line_RF = const_delta.lambda_RF_line_SIIV
			self._lines1 = const_delta.SIIV_lines
			self._name_line1 = const_delta.SIIV_lines_names
		elif (self._f1=='CIV'):
			self._line_RF = const_delta.lambda_RF_line_CIV
			self._lines1 = const_delta.CIV_lines
			self._name_line1 = const_delta.CIV_lines_names
		elif (self._f1=='MGII'):
			self._line_RF = const_delta.lambda_RF_line_MGII
			self._lines1 = const_delta.MGII_lines
			self._name_line1 = const_delta.MGII_lines_names

		### folder where data are
		self._path_to_txt_file_folder = dic['path_to_txt_file_folder']

		### Correlation type
		self._correlation = dic['correlation']
		if (self._correlation=='q_f'):
			self._label = '\\xi^{qf}'
			self._xTitle = '\lambda_{Obs., pix}/\lambda_{Obs., QSO}'
			self._yTitle = self._label
			self._title = '\delta_{'+self._f1+'} \, - \, '+self._q1
			self._prefix = 'xi_delta_QSO_lambda_same_LOS'
			self._middlefix = self._f1+'_'+self._q1

		path = self._path_to_txt_file_folder + self._prefix + '_Mu_' + self._middlefix + '.txt'

		print
		print '  path = ', path

		### Name of the correlation
		self._name = dic['name']

		### Read data
		self._xi1D = self.read_data(path)

		return
	def read_data(self, path):

		xi1D = numpy.zeros(shape=(self._nbBinL1oL2,3))

		tmp_save000 = numpy.zeros(self._nbBinL1oL2)
		tmp_save111 = numpy.zeros(self._nbBinL1oL2)
		tmp_save222 = numpy.zeros(self._nbBinL1oL2)
		tmp_save333 = numpy.zeros(self._nbBinL1oL2)
		tmp_save444 = numpy.zeros(self._nbBinL1oL2)
		tmp_save555 = numpy.zeros(self._nbBinL1oL2)
		tmp_save666 = numpy.zeros(self._nbBinL1oL2)

		data = numpy.loadtxt(path)
		save0 = data[:,0]
		save1 = data[:,1]
		save2 = data[:,2]
		save3 = data[:,3]
		save4 = data[:,4]
		save5 = data[:,5]

		for i in range(0,len(save0)):
			### for xi1D
			tmp_save000[i] += save0[i]
			tmp_save111[i] += save1[i]
			tmp_save222[i] += save2[i]
			tmp_save333[i] += save3[i]
			tmp_save444[i] += save4[i]
			tmp_save555[i] += save5[i]
		
		cut = (tmp_save555>1.)
		xi1D[:,0][cut] = tmp_save222[cut]/tmp_save444[cut]
		xi1D[:,1][cut] = tmp_save000[cut]/tmp_save444[cut]
		xi1D[:,2][cut] = numpy.sqrt( (tmp_save111[cut]/tmp_save444[cut] - xi1D[:,1][cut]*xi1D[:,1][cut])/tmp_save555[cut] )

		return xi1D
	def plot_1d(self, with_lines=False, other=[], verbose=False):
	
		list_corr_to_plot = numpy.append( [self],other )
		
                xMin = numpy.amin(self._xi1D[:,0])
                xMax = numpy.amax(self._xi1D[:,0])
                yMin = numpy.amin(self._xi1D[:,1])
                yMax = numpy.amax(self._xi1D[:,1])

                for el in list_corr_to_plot:

			cut = (el._xi1D[:,2]>0.)
			xxx = el._xi1D[:,0][cut]
			yyy = el._xi1D[:,1][cut]
			yer = el._xi1D[:,2][cut]
			plt.errorbar(self._line_RF*xxx, yyy, yerr=yer, fmt='o', label=r'$'+el._name+'$')
			plt.errorbar(self._line_RF*xxx, yyy, color='blue')

			xMin = min(xMin, numpy.amin(xxx) )
                        xMax = max(xMax, numpy.amax(xxx) )
                        yMin = min(yMin, numpy.amin(yyy) )
                        yMax = max(yMax, numpy.amax(yyy) )

		if (with_lines):

			if (verbose): print ' ||  name_1 - name_2 || line || lambda_rf_1 || lamnda_rf_2 || '

			yLi = [yMin,yMax]
			nbLines1 = self._lines1.size

			for i in range(0,nbLines1):
				nbLines2 = 1

				for j in range(0,nbLines2):
					line = self._lines1[i] #/self._line_RF

					if (line<self._line_RF*xMin or line>self._line_RF*xMax): continue
					xLi = [line,line]
					plt.plot(xLi,yLi,color='green')

					if (verbose): print ' ||  ', self._name_line1[i] ,' - QSO || ', line, ' || ', self._lines1[i], ' || ', self._line_RF, ' || '
					name = self._name_line1[i]+' - QSO '
		
					plt.text(line, yMax, name, rotation='vertical', fontsize=20)



		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$'+self._xTitle+'$', fontsize=40)
		plt.ylabel(r'$'+self._yTitle+'$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.show()

		return


dic = {
	'min_l1_over_l2': 0.40,
	'max_l1_over_l2': 2.30,
	'size_bin_calcul_l1_overl2': 1.e-03,
	'correlation': 'q_f',
	'path_to_txt_file_folder': 'NOTHING',
	'f1': 'CIV',
	'q': 'QSO_DR14_v1_0',
	'name' : 'Data'
}


corr_list = []
dic['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/'
a = CorrelationLambda(dic)
corr_list += [a]
a.plot_1d(True,corr_list[1:])



