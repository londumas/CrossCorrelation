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

		### theta
		self._minTheta = dic['min_theta']
		self._maxTheta = dic['max_theta']
		self._nbBinTheta = int( (self._maxTheta-self._minTheta)/dic['size_bin_calcul_theta'] )

		### visualisation in mu
		self._minVisu = dic['min_visual_z']
		self._maxVisu = dic['max_visual_z']

		### forest and QSO name
		self._f1 = dic['f1']
		self._f2 = dic['f2']
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

		### folder where data are
		self._path_to_txt_file_folder = dic['path_to_txt_file_folder']

		### Correlation type
		self._correlation = dic['correlation']
		if (self._correlation=='q_f'):
			self._label = '\\xi^{qf}'
			self._xTitle = '\lambda_{Obs., pix}/\lambda_{Obs., QSO}'
			self._yTitle = self._label+' \, (\\theta<'+str(self._maxTheta)+' \, rad)'
			self._title = '\delta_{'+self._f1+'} \, - \, '+self._q1
			self._prefix = 'xi_delta_QSO_lambda'
			self._middlefix = self._f1+'_'+self._q1
		elif (self._correlation=='f_f'):
			self._label = '\\xi^{ff}'
			self._xTitle = '\lambda_{Obs.,2}/\lambda_{Obs.,1}'
			self._yTitle = self._label+' \, (\\theta<'+str(self._maxTheta)+' \, rad)'
			self._title = '\delta_{'+self._f1+'}'
			self._prefix = 'xi_A_delta_delta_lambda'
			self._middlefix = self._f1
		elif (self._correlation=='f_f2'):
			self._label = '\\xi^{ff}'
			self._xTitle = '\lambda_{Obs.,2}/\lambda_{Obs.,1}'
			self._yTitle = self._label+' \, (\\theta<'+str(self._maxTheta)+' \, rad)'
			self._title = '\delta_{'+self._f1+'} \, - \, \delta_{'+self._f2+'}'
			self._prefix = 'xi_A_delta_delta2_lambda'
			self._middlefix = self._f1+'_'+self._f2

		path = self._path_to_txt_file_folder + self._prefix + '_Mu_' + self._middlefix + '.txt'

		print
		print '  path = ', path

		### Name of the correlation
		self._name = dic['name']

		### Read data
		self._xiMu, self._xiWe, self._xi1D = self.read_data(path)

		return
	def read_data(self, path):

		xiMu = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2,4))
		xiWe = numpy.zeros(shape=(self._nbBinL1oL2,3,3))
		xi1D = numpy.zeros(shape=(self._nbBinL1oL2,3))

		tmp_save0 = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2))
		tmp_save1 = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2))
		tmp_save2 = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2))
		tmp_save3 = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2))
		tmp_save4 = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2))
		tmp_save5 = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2))
		tmp_save6 = numpy.zeros(shape=(self._nbBinTheta,self._nbBinL1oL2))
	
		tmp_save00 = numpy.zeros(shape=(self._nbBinL1oL2,3))
		tmp_save11 = numpy.zeros(shape=(self._nbBinL1oL2,3))
		tmp_save22 = numpy.zeros(shape=(self._nbBinL1oL2,3))
		tmp_save33 = numpy.zeros(shape=(self._nbBinL1oL2,3))
		tmp_save44 = numpy.zeros(shape=(self._nbBinL1oL2,3))
		tmp_save55 = numpy.zeros(shape=(self._nbBinL1oL2,3))
		tmp_save66 = numpy.zeros(shape=(self._nbBinL1oL2,3))
	
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
		save6 = data[:,6]

		for i in range(0,len(save0)):
			iX = i/self._nbBinTheta
			iY = i%self._nbBinTheta
			idY = iY
	
			tmp_save0[idY][iX] += save0[i]
			tmp_save1[idY][iX] += save1[i]
			tmp_save2[idY][iX] += save2[i]
			tmp_save3[idY][iX] += save3[i]
			tmp_save4[idY][iX] += save4[i]
			tmp_save5[idY][iX] += save5[i]
			tmp_save6[idY][iX] += save6[i]
	
			if (iY>2): continue
	
			tmp_save00[iX][idY] += save0[i]
			tmp_save11[iX][idY] += save1[i]
			tmp_save22[iX][idY] += save2[i]
			tmp_save33[iX][idY] += save3[i]
			tmp_save44[iX][idY] += save4[i]
			tmp_save55[iX][idY] += save5[i]
			tmp_save66[iX][idY] += save6[i]
	
			### for xi1D
			tmp_save000[iX] += save0[i]
			tmp_save111[iX] += save1[i]
			tmp_save222[iX] += save2[i]
			tmp_save333[iX] += save3[i]
			tmp_save444[iX] += save4[i]
			tmp_save555[iX] += save5[i]
			tmp_save666[iX] += save6[i]
		
		cut = (tmp_save6>1.)
		xiMu[:,:,0][cut] = tmp_save2[cut]/tmp_save5[cut]
		xiMu[:,:,1][cut] = tmp_save3[cut]/tmp_save5[cut]
		xiMu[:,:,2][cut] = tmp_save0[cut]/tmp_save5[cut]
		xiMu[:,:,3][cut] = numpy.sqrt( (tmp_save1[cut]/tmp_save5[cut] - xiMu[:,:,2][cut]*xiMu[:,:,2][cut])/tmp_save6[cut])
	
		cut = (tmp_save66>1.)
		xiWe[:,:,0][cut] = tmp_save22[cut]/tmp_save55[cut]
		xiWe[:,:,1][cut] = tmp_save00[cut]/tmp_save55[cut]
		xiWe[:,:,2][cut] = numpy.sqrt( (tmp_save11[cut]/tmp_save55[cut] - xiWe[:,:,1][cut]*xiWe[:,:,1][cut])/tmp_save66[cut] )
	
		cut = (tmp_save666>1.)
		xi1D[:,0][cut] = tmp_save222[cut]/tmp_save555[cut]
		xi1D[:,1][cut] = tmp_save000[cut]/tmp_save555[cut]
		xi1D[:,2][cut] = numpy.sqrt( (tmp_save111[cut]/tmp_save555[cut] - xi1D[:,1][cut]*xi1D[:,1][cut])/tmp_save666[cut] )

		return xiMu, xiWe, xi1D
	def plot_1d(self, with_lines=False, verbose=False):
	
		cut = (self._xi1D[:,2]>0.)
		xxx = self._xi1D[:,0][cut]
		yyy = self._xi1D[:,1][cut]
		yer = self._xi1D[:,2][cut]
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o', label=r'$'+self._name+'$')

		if (with_lines):

			if (verbose): print ' ||  name_1 - name_2 || line || lambda_rf_1 || lamnda_rf_2 || '

			xMin = numpy.min(xxx)
			xMax = numpy.max(xxx)
			yMin = numpy.min(yyy)
			yMax = numpy.max(yyy)
			yLi = [yMin,yMax]
			nbLines1 = self._lines1.size

			for i in range(0,nbLines1):

				if (self._correlation=='f_f'):
					nbLines2 = i
				if (self._correlation=='q_f'):
					nbLines2 = 1
				else:
					nbLines2 = self._lines2.size

				for j in range(0,nbLines2):

					if (self._correlation=='q_f'):
						line = self._lines1[i]/self._line_RF
					else:
						line = self._lines2[j]/self._lines1[i]
						if (line<xMin or line>xMax): line = 1./line
					if (line<xMin or line>xMax): continue
					xLi = [line,line]
					plt.plot(xLi,yLi,color='green')

					if (self._correlation=='q_f'):
						if (verbose): print ' ||  ', self._name_line1[i] ,' - QSO || ', line, ' || ', self._lines1[i], ' || ', self._line_RF, ' || '
						name = self._name_line1[i]+' - QSO'
					else:
						if (verbose): print ' ||  ', self._name_line1[i] ,' - ', self._name_line2[j], ' || ', line, ' || ', self._lines1[i], ' || ', self._lines2[j], ' || '
						name = self._name_line1[i]+' - ' + self._name_line2[j]
		
					plt.text(line, yMax, name, rotation='vertical', fontsize=20)



		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$'+self._xTitle+'$', fontsize=40)
		plt.ylabel(r'$'+self._yTitle+'$', fontsize=40)
		myTools.deal_with_plot(False,False,True)
		plt.show()

		return
	def plot_we(self, with_lines=False, verbose=False):
	
		label = ["\\theta < 0.0001", "0.0001 < \\theta < 0.0002", "0.0002 < \\theta < 0.0003"]
	
		for i in numpy.arange(3):

			cut = (self._xiWe[:,i,2]>0.)
			if (self._xiWe[:,i,0][cut].size==0):
				continue
	
			xxx = self._xiWe[:,i,0][cut]
			yyy = self._xiWe[:,i,1][cut]
			yer = self._xiWe[:,i,2][cut]
			
			plt.errorbar(xxx, yyy, yerr=yer, fmt='o', label=r'$'+label[i]+'$')

		if (with_lines):

			if (verbose): print ' ||  name_1 - name_2 || line || lambda_rf_1 || lamnda_rf_2 || '

			xMin = numpy.min(xxx)
			xMax = numpy.max(xxx)
			yMin = numpy.min(yyy)
			yMax = numpy.max(yyy)
			yLi = [yMin,yMax]
			nbLines1 = self._lines1.size

			for i in range(0,nbLines1):

				if (self._correlation=='f_f'):
					nbLines2 = i
				if (self._correlation=='q_f'):
					nbLines2 = 1
				else:
					nbLines2 = self._lines2.size

				for j in range(0,nbLines2):

					if (self._correlation=='q_f'):
						line = self._lines1[i]/self._line_RF
					else:
						line = self._lines2[j]/self._lines1[i]
						if (line<xMin or line>xMax): line = 1./line
					if (line<xMin or line>xMax): continue
					xLi = [line,line]
					plt.plot(xLi,yLi,color='green')

					if (self._correlation=='q_f'):
						if (verbose): print ' ||  ', self._name_line1[i] ,' - QSO || ', line, ' || ', self._lines1[i], ' || ', self._line_RF, ' || '
						name = self._name_line1[i]+' - QSO'
					else:
						if (verbose): print ' ||  ', self._name_line1[i] ,' - ', self._name_line2[j], ' || ', line, ' || ', self._lines1[i], ' || ', self._lines2[j], ' || '
						name = self._name_line1[i]+' - ' + self._name_line2[j]
		
					plt.text(line, yMax, name, rotation='vertical', fontsize=20)

		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$'+self._xTitle+'$', fontsize=40)
		plt.ylabel(r'$'+self._yTitle+'$', fontsize=40)
		plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=1)
		plt.xlim([ self._minL1oL2-0.01, self._maxL1oL2+0.01 ])
		myTools.deal_with_plot(False,False,True)
		plt.show()

		return
	def plot_mu(self, with_lines=False, verbose=False):
	
		xxx = self._xiMu[:,:,0]
		muu = self._xiMu[:,:,1]
		yyy = self._xiMu[:,:,2]
		yer = self._xiMu[:,:,3]

		cut = (yer<=0.)
		yyy[cut] = float('nan')
		yer[cut] = float('nan')
	
		fig = plt.figure()
		ax = fig.add_subplot(111)
		extent=[ self._minL1oL2, self._maxL1oL2, self._minTheta, self._maxTheta]
	
		plt.imshow(yyy, origin='lower',aspect='auto',extent=extent,interpolation='None',vmin=self._minVisu, vmax=self._maxVisu)
		cbar = plt.colorbar()

		

		'''
		plt.xlim([ min1D__, max1D__ ])
		plt.ylim([ minTheta_, maxTheta_ ])
	
		yMin    = numpy.min(muu)
		yMax    = 0.95*numpy.max(muu)
		nbLines = lines.size
		for i in range(0,nbLines):
			line = lines[i]/lambdaRFLine
			if (line<min1D__ or line>max1D__): continue
			xLi  = [line,line]
			yLi  = [yMin,yMax]
			name = 'QSO - ' + names[i]
			plt.plot(xLi,yLi,color='green',linewidth=2)
			plt.text(line, yMax, name, rotation='vertical', fontsize=20)
		'''

		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$'+self._xTitle+'$', fontsize=40)
		plt.ylabel(r'$\theta \, [rad]$', fontsize=40)
		cbar.set_label(r'$'+self._yTitle+'$',size=40)
		myTools.deal_with_plot(False,False,False)
		plt.show()

		return

dic = {
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
	'f1': 'CIV',
	'f2': 'a',
	'q': 'QSO_DR7_DR12_EBOSS',
	'name' : 'Data'
}

### For f_f
'''
dic['min_l1_over_l2'] = 1.
dic['max_l1_over_l2'] = 1.11
dic['correlation']    = 'f_f'
dic['min_visual_z'] = -0.002
dic['max_visual_z'] = 0.002
'''
### For f_f2
'''
dic['min_l1_over_l2'] = 0.89
dic['max_l1_over_l2'] = 1.11
dic['correlation']    = 'f_f2'
dic['min_visual_z'] = -0.002
dic['max_visual_z'] = 0.002
'''

dic['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
a = CorrelationLambda(dic)
a.plot_1d()
a.plot_we()
a.plot_mu()


























