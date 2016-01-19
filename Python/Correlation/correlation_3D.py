# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#
#

import numpy
import matplotlib.pyplot as plt


dic_constants = {
	'minXi': 0.,
	'maxXi': 200.,
	'nbBin': 50,
	'nbBinM': 25,
	'nb_Sub_Sampling': 80,
	'size_bin_calcul_s': 1.,
	'size_bin_calcul_m': 0.02,
	'correlation': 'q_f',
	'path_To_Txt_File': ['NOTHING','NOTHING'] }


class Correlation_3D:
	
	def __init__(self, dic_constants):
		"""
		
		correlationType:
			- 'q_q'
			- 'q_f' (default)
			- 'f_f'
			- 'f_f2'
	
		"""
		
		### Correlation type
		self._correlation = dic_constants['correlation']
		if (self._correlation=='q_q'):
			self._label = '\\xi^{qq}'
		elif (self._correlation=='q_f'):
			self._label = '\\xi^{qf}'
		elif (self._correlation=='f_f'):
			self._label = '\\xi^{ff}'
		elif (self._correlation=='f_f2'):
			self._label = '\\xi^{ff}'
		else:
			print "  Correlation_3D::__init__::  ERROR:  'correlation' is incorrect "
			return

		### 1D
		self._min1D   = float(dic_constants['minXi'])
		self._max1D   = float(dic_constants['maxXi'])
		self._nbBin1D = int(dic_constants['nbBin'])
		self._binSize = (self._max1D-self._min1D)/self._nbBin1D
		self._binSize_calcul = int( (self._max1D-self._min1D)/dic_constants['size_bin_calcul_s'] )
		
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
		self._nbBinM = dic_constants['nbBinM']
		self._nbBinM_calcul = int( (self._maxM-self._minM)/dic_constants['size_bin_calcul_m'] )
			
		### Sub sampling
		self.nb_Sub_Sampling = int(dic_constants['nb_Sub_Sampling'])

		### Correlations data
		self._xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
		self._xiMu = numpy.zeros(shape=(self._nbBin1D,self._nbBinM,4))
		self._xiWe = numpy.zeros(shape=(self._nbBin1D,3,3))
		self._xi1D = numpy.zeros(shape=(self._nbBin1D,3))
		self.fill_data(dic_constants['path_To_Txt_File'][0],dic_constants['path_To_Txt_File'][1])
		
		return
	
	def fill_data(self, path1D, path2D, selection=0):
		'''
	
		selection:
			- 0: do all
			- 1: only 1D
			- 2: only 2D
	
		'''
		
		int_binSize = int(self._binSize)
	
		if (selection==0 or selection==1):
			### Mu
			data = numpy.loadtxt(path1D)
			print path1D
	
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
		
			self._xiMu[:,:,0] = tmp_save2/tmp_save5
			self._xiMu[:,:,1] = tmp_save3/tmp_save5
			self._xiMu[:,:,2] = tmp_save0/tmp_save5
			self._xiMu[:,:,3] = numpy.sqrt( (tmp_save1/tmp_save5 - self._xiMu[:,:,2]*self._xiMu[:,:,2])/tmp_save6)
			cut = (tmp_save5==0.)
			self._xiMu[:,:,0][cut] = 0.
			self._xiMu[:,:,1][cut] = 0.
			self._xiMu[:,:,2][cut] = 0.
			self._xiMu[:,:,3][cut] = 0.
	
			self._xiWe[:,:,0] = tmp_save22/tmp_save55
			self._xiWe[:,:,1] = tmp_save00/tmp_save55
			self._xiWe[:,:,2] = numpy.sqrt( (tmp_save11/tmp_save55 - self._xiWe[:,:,1]*self._xiWe[:,:,1])/tmp_save66 )
			cut = (tmp_save55==0.)
			self._xiWe[:,:,0][cut] = 0.
			self._xiWe[:,:,1][cut] = 0.
			self._xiWe[:,:,2][cut] = 0.
	
			self._xi1D[:,0] = tmp_save222/tmp_save555
			self._xi1D[:,1] = tmp_save000/tmp_save555
			self._xi1D[:,2] = numpy.sqrt( (tmp_save111/tmp_save555 - self._xi1D[:,1]*self._xi1D[:,1])/tmp_save666 )
			cut = (tmp_save555==0.)
			self._xi1D[:,0][cut] = 0.
			self._xi1D[:,1][cut] = 0.
			self._xi1D[:,2][cut] = 0.
			
			
	
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
				
			self._xi2D[:,:,0] = numpy.sqrt( (tmp_save2/tmp_save5)**2. + (tmp_save3/tmp_save5)**2. )
			self._xi2D[:,:,1] = tmp_save0 / tmp_save5
			self._xi2D[:,:,2] = numpy.sqrt( (tmp_save1/tmp_save5 - self._xi2D[:,:,1]*self._xi2D[:,:,1])/tmp_save6 )
			cut = (tmp_save5 == 0.)
			self._xi2D[:,:,0][cut] = 0.
			self._xi2D[:,:,1][cut] = 0.
			self._xi2D[:,:,2][cut] = 0.

		return
	def plot_1d(self, x_power=0):

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
		plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o')
		
		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
		#plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		#myTools.deal_with_plot(False,False,False)
		plt.show()
		
		return
	def plot_2d(self, x_power=0):

		origin='lower'
		extent=[self._minX2D, self._maxX2D, self._minY2D, self._maxY2D]
		if (self._correlation=='q_f'):
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
	
		#plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
		plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
		plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
	
		#plt.xlim([ self._minX2D, self._maxX2D ])
		#plt.ylim([ self._maxY2D, self._minY2D ])
	
		plt.show()
		
		return
	def plot_Mu(self, x_power=0):
	
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
		#plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
		plt.xlabel(r'$\mu$', fontsize=40)
		plt.ylabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
		plt.grid(True)
		cbar.formatter.set_powerlimits((0, 0))
		cbar.update_ticks()
		plt.show()
	
		return
	def plot_We(self, x_power=0):
	
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
		
		#plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		#myTools.deal_with_plot(False,False,False)
		plt.show()
		
		return
		
		
'''	
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
path1D = path + 'xi_A_delta_delta_Mu_LYA.txt'
path2D = path + 'xi_A_delta_delta_2D_LYA.txt'
dic_constants['path_To_Txt_File'][0] = path1D
dic_constants['path_To_Txt_File'][1] = path2D
dic_constants['correlation'] = 'f_f'
'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
path1D = path + 'xi_delta_QSO_Mu_LYA_QSO.txt'
path2D = path + 'xi_delta_QSO_2D_LYA_QSO.txt'
dic_constants['path_To_Txt_File'][0] = path1D
dic_constants['path_To_Txt_File'][1] = path2D
dic_constants['correlation'] = 'q_f'

corr = Correlation_3D(dic_constants) 	

corr.plot_1d(0)
corr.plot_1d(1)
corr.plot_1d(2)
corr.plot_2d(0)
corr.plot_2d(1)
corr.plot_2d(2)
corr.plot_We(0)
corr.plot_We(1)
corr.plot_We(2)
corr.plot_Mu(0)
corr.plot_Mu(1)
corr.plot_Mu(2)
		
		
		
	
