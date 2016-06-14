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

### Perso lib
import correlation_3D
import myTools

par_name = numpy.asarray(['\\beta_{f}','b_{f} \cdot (1+\\beta_{f})','\\gamma_{b_{f}}','\\gamma_{\\beta_{f}}','\\Delta v','b_{q}','b_{q} \cdot \\beta_{q}',
	'1+f','SigmaNL-perp',
	'BAO \, amplitude','\\alpha_{iso}','\\alpha_{\parallel}','\\alpha_{\perp}',
	'gamma-scale','Rad \, strength','Rad \, anisotropy','Rad \, mean \, free \, path',
	'Rad \, quasar \, lifetime','a0','a1','a2','a3','a4','a5','a6','a7','a7','a8','a9','a10','a11','a12','a13'])
raw_index_parameter = {
	'beta'                        : 0,
	'b.(1+beta)'                  : 1,
	'gamma-bias'                  : 2,
	'gamma-beta'                  : 3,
	'Delta v'                     : 4,
	'b_{2}'                       : 5,
	'b_{2}.beta_{2}'              : 6,
	'SigmaNL-perp'                : 7,
	'1+f'                         : 8,
	'BAO \, amplitude'            : 9,
	'alpha_{iso}'                 : 10,
	'alpha_paral'                 : 11, 
	'alpha_perp'                  : 12,
	'gamma-scale'                 : 13,
	'pixel-scale'                 : 14,
	'beta Si2a'                   : 15,
	'bias Si2a'                   : 16,
	'beta Si2b'                   : 17,
	'bias Si2b'                   : 18,
	'beta Si2c'                   : 19,
	'bias Si2c'                   : 20,
	'beta Si3'                    : 21,
	'bias Si3'                    : 22
}
par_name_f = numpy.asarray(['\\beta_{f}','b_{f} \cdot  (1+\\beta_{f})','\\gamma_{b_{f}}','\\gamma_{\\beta}','SigmaNL-perp','1+f',
	'BAO \, amplitude','\\alpha_{iso}','\\alpha_{\parallel}','\\alpha_{\perp}', 'gamma-scale',
	'pixel-scale'])
raw_index_parameter_f = {
	'beta'                        : 0,
	'b.(1+beta)'                  : 1,
	'gamma-bias'                  : 2,
	'gamma-beta'                  : 3,
	'SigmaNL-perp'                : 4,
	'1+f'                         : 5,
	'BAO \, amplitude'            : 6,
	'alpha_{iso}'                 : 7,
	'alpha_paral'                 : 8, 
	'alpha_perp'                  : 9,
	'gamma-scale'                 : 10,
	'pixel-scale'                 : 11
}

class AnnalyseBAOFIT(correlation_3D.Correlation3D):

	def __init__(self, dic=None, index_parameter=None, path_to_BAOFIT=None):

		correlation_3D.Correlation3D.__init__(self,dic)

		### index of parameters
		if (index_parameter is None):
			if (self._correlation=='q_f'):
				self._index_parameter = copy.deepcopy(raw_index_parameter)
			elif (self._correlation=='f_f'):
				self._index_parameter = copy.deepcopy(raw_index_parameter_f)
		else:
			self._index_parameter = copy.deepcopy(index_parameter)

		if (path_to_BAOFIT is None):
			if (self._correlation=='q_q'):
				self._path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._q1 + '/bao2D.'
			elif (self._correlation=='q_f'):
				self._path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1+'__'+self._q1 + '/bao2D.'
			elif (self._correlation=='f_f'):
				self._path_to_BAOFIT = self._path_to_txt_file_folder + 'BaoFit_'+self._correlation+'__'+self._f1 + '/bao2D.'
		else:
			self._path_to_BAOFIT = path_to_BAOFIT

		print self._path_to_BAOFIT

		### Set attributes set after
		if (self._correlation=='q_f'):
			self._par_name       = copy.deepcopy(par_name)
		elif (self._correlation=='f_f'):
			self._par_name       = copy.deepcopy(par_name_f)
		self._xi2D_fit       = None
		self._nbParam        = None
		self._param          = None
		self._chi2           = None
		self._chi2_scan      = None
		self._chi2_scan_edge = None

		### Get BAOFIT results
		self.read_fit_data()

		return

	def read_fit_data(self):

		### String names
		residuals = 'residuals.dat'
		param     = 'save.pars'
		chi2      = 'fit.chisq'	
		
		### Get data and fit
		self._xi2D_fit = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,9))
		data  = numpy.loadtxt(self._path_to_BAOFIT+residuals)
		index = data[:,0].astype(int)
		for el in data:
			i = int(el[0])
			ix = i%self._nbBinX2D
			iy = i/self._nbBinX2D
			self._xi2D_fit[ix,iy,:] = el[1:]

	
		### Get the parameters of the fit
		data    = numpy.loadtxt(self._path_to_BAOFIT+param)
		self._nbParam = data[:,1].size
		self._param   = numpy.zeros( shape=(self._nbParam,2) )
		self._param[:,0] = data[:,1]
		self._param[:,1] = data[:,2]
	
		### Get the parameters of the fit
		data = numpy.loadtxt(self._path_to_BAOFIT+chi2)
		self._chi2 = data
	
		return
	def read_chi2_scan(self, sizeX=100, sizeY=100):

		### Constants
		name = 'scan.dat'
		i_alpha_paral = self._index_parameter['alpha_paral']
		i_alpha_perp  = self._index_parameter['alpha_perp']

		path_to_BAOFIT = self._path_to_BAOFIT
		#path_to_BAOFIT = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/BACKUP_2015_01_15/FitsFile_DR12_Guy/BaoFit_q_f/bao2D.'
		
		### Create a file a tmp with the scan minus the two first lines
		idx = 0
		f     = open(path_to_BAOFIT+name)
		tmp_f = open(self._path_to_BAOFIT+'tmp_scan.txt','w')
		for line in f:
			if (idx>=2): tmp_f.write(line)
			idx += 1
		tmp_f.close()
		f.close()
		
		### Data
		data = numpy.loadtxt(self._path_to_BAOFIT+'tmp_scan.txt')

		### delete tmp file
		subprocess.call('rm '+self._path_to_BAOFIT+'tmp_scan.txt', shell=True)
		
		### Scan
		bestFit = numpy.zeros( shape=(2,1) )
		bestFit[0,0] = data[0,i_alpha_perp]
		bestFit[1,0] = data[0,i_alpha_paral]
		chi2_bestFit = data[0][-1]
		chi2 = data[1:,-1]
		
		### Put data in 2D array
		self._chi2_scan = numpy.zeros( shape=(sizeX,sizeY) )
		for k in numpy.arange(0,chi2.size):
			i = k/sizeY
			j = k%sizeY
			self._chi2_scan[i][j] = chi2[k]
		
		### Get the delta of chi^2
		self._chi2_scan[ (self._chi2_scan==0.) ] = float('nan')
		self._chi2_scan[ (self._chi2_scan!=0.) ] = self._chi2_scan[ (self._chi2_scan!=0.) ]-chi2_bestFit

		self._chi2_scan_edge = [numpy.amin(data[1:,i_alpha_perp]), numpy.amax(data[1:,i_alpha_perp]), numpy.amin(data[1:,i_alpha_paral]), numpy.amax(data[1:,i_alpha_paral])]

		return bestFit
	def read_toyMC(self):

		### Constants
		name = 'toymc.dat'
		i_alpha_paral = self._index_parameter['alpha_paral']
		i_alpha_perp  = self._index_parameter['alpha_perp']

		#path_to_BAOFIT = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/BACKUP_2015_01_15/FitsFile_DR12_Guy/BaoFit_q_f/toymc/bao2D.'

		## Create a file a tmp with the toyMc minus the two first lines
		idx = 0
		f     = open(self._path_to_BAOFIT+name)
		tmp_f = open(self._path_to_BAOFIT+'tmp_toymc.txt','w')
		for line in f:
			if (idx>=2): tmp_f.write(line)
			idx += 1
		tmp_f.close()
		f.close()

		### Data
		data = numpy.loadtxt(self._path_to_BAOFIT+'tmp_toymc.txt')

		### delete tmp file
		subprocess.call('  rm '+self._path_to_BAOFIT+'tmp_toymc.txt', shell=True)

		### Best fit and toy-MC
		data = data[1:]

		print '  nb toy model done = ',  data[:,0].size
		data = data[ data[:,-1]!=0 ]
		nbToyMC = data[:,0].size
		print '  nb toy model done correct = ', nbToyMC

		bestFit  = numpy.zeros( shape=(2,nbToyMC) )
		bestFit[1,:] = data[:,i_alpha_paral]
		bestFit[0,:] = data[:,i_alpha_perp]

		return bestFit
	def read_bootStrap(self):

		### String names
		param     = 'save.pars'
		i_alpha_paral = self._index_parameter['alpha_paral']
		i_alpha_perp  = self._index_parameter['alpha_perp']

		path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/BACKUP_2015_01_15/FitsFile_DR12_Guy/BaoFit_q_f/Bootstraps/'

		### Get best fit plus bootstrap
		nbBoot = 10000
		bestFit  = numpy.zeros( shape=(2,nbBoot) )
		nb = 0
		for i in range(0,nbBoot):
			try:
				data     = numpy.loadtxt(path + 'boot_'+str(i).zfill(4)+'/bao2D.'+param)
				bestFit[0,nb] = data[i_alpha_perp,1]
				bestFit[1,nb] = data[i_alpha_paral,1]
				nb += 1
			except Exception:
				print '  bootstrap missing : ', i

		return bestFit
	def get_best_fit(self):

		bestFit = numpy.zeros( shape=(2,1) )
		bestFit[0,0] = self._param[self._index_parameter['alpha_perp'],0]
		bestFit[1,0] = self._param[self._index_parameter['alpha_paral'],0]

		return bestFit
	def get_residuals(self):
		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
		xi2D[:,:,0] = self._xi2D_fit[:,:,3]
		xi2D[:,:,1] = self._xi2D_fit[:,:,6]
		xi2D[:,:,2] = self._xi2D_fit[:,:,8]
		cut = (xi2D[:,:,2]>0.)
		xi2D[:,:,1][cut] = (self._xi2D[:,:,1][cut]-xi2D[:,:,1][cut])/xi2D[:,:,2][cut]

		return xi2D
	def print_results(self):

		### Get the chi^2
		tmp_chi2    = self._chi2[0]
		tmp_NBBin   = self._chi2[1]
		tmp_NBParam = self._chi2[2]
		print
		print '  chi^{2} = ', tmp_chi2
	
		### Get the fit parameters
		tmp_beta               = self._param[ self._index_parameter['beta'],0]
		tmp_beta_err           = self._param[ self._index_parameter['beta'],1]
		tmp_delta_v            = self._param[ self._index_parameter['Delta v'],0]
		tmp_delta_v_err        = self._param[ self._index_parameter['Delta v'],1]
		tmp_bias2              = self._param[ self._index_parameter['b_{2}'],0]
		tmp_bias2_err          = self._param[ self._index_parameter['b_{2}'],1]
		tmp_beta2_bias2        = self._param[ self._index_parameter['b_{2}.beta_{2}'],0]
		tmp_beta2_bias2_err    = self._param[ self._index_parameter['b_{2}.beta_{2}'],1]
		tmp_alpha_parallel     = self._param[ self._index_parameter['alpha_paral'],0]
		tmp_alpha_parallel_err = self._param[ self._index_parameter['alpha_paral'],1]
		tmp_alpha_perp         = self._param[ self._index_parameter['alpha_perp'],0]
		tmp_alpha_perp_err     = self._param[ self._index_parameter['alpha_perp'],1]

		string =  """
|| chi2 / dof              || beta                  || delta-v                 || bias2                 || beta2*bias2           ||  BAO alpha-parallel   ||  BAO alpha-perp       ||
|| %1.3e / (%u - %u) || %1.3e ± %1.3e || %1.3e ±  %1.3e || %1.3e ± %1.3e || %1.3e ± %1.3e || %1.3e ± %1.3e || %1.3e ± %1.3e ||
""" % (tmp_chi2,tmp_NBBin,tmp_NBParam,tmp_beta,tmp_beta_err,tmp_delta_v,tmp_delta_v_err,tmp_bias2,tmp_bias2_err,tmp_beta2_bias2,tmp_beta2_bias2_err,tmp_alpha_parallel,tmp_alpha_parallel_err,tmp_alpha_perp,tmp_alpha_perp_err)
		print string

		return
	def get_data_and_fit_1d(self,path_to_mapping_1D):

		### Load the mapping from 2D to 1D
		mapping_2D_to_1D = numpy.load( path_to_mapping_1D )

		### Fit
		xi1D_data = numpy.zeros(shape=(self._nbBin1D,3))
		xi1D_fit  = numpy.zeros(shape=(self._nbBin1D,3))
		for i in numpy.arange(self._nbBinX2D):
			for j in numpy.arange(self._nbBinY2D):

				if (self._xi2D_fit[i,j,8]<=0.): continue
				ivar = 1./( self._xi2D_fit[i,j,8]*self._xi2D_fit[i,j,8] )
				for k in numpy.arange(self._nbBin1D):

					coef = mapping_2D_to_1D[i,j,k]
					if (coef==0.): continue

					xi1D_data[k,0] += coef*ivar*self._xi2D_fit[i,j,3]
					xi1D_data[k,1] += coef*ivar*self._xi2D_fit[i,j,7]
					xi1D_data[k,2] += coef*ivar

					xi1D_fit[k,0] += coef*ivar*self._xi2D_fit[i,j,3]
					xi1D_fit[k,1] += coef*ivar*self._xi2D_fit[i,j,6]
					xi1D_fit[k,2] += coef*ivar

		### Data
		cut = (xi1D_data[:,2]>0.)
		xi1D_data[:,0][cut] /= xi1D_data[:,2][cut]
		xi1D_data[:,1][cut] /= xi1D_data[:,2][cut]
		xi1D_data[:,2][cut]  = 1./numpy.sqrt(xi1D_data[:,2][cut])
		### Fit
		cut = (xi1D_fit[:,2]>0.)
		xi1D_fit[:,0][cut] /= xi1D_fit[:,2][cut]
		xi1D_fit[:,1][cut] /= xi1D_fit[:,2][cut]
		xi1D_fit[:,2][cut]  = 1./numpy.sqrt(xi1D_fit[:,2][cut])

		return xi1D_data, xi1D_fit
	def get_data_and_fit_we(self,path_to_mapping):

		### Load the mapping from 2D to we
		mapping_2D_to_we = numpy.load( path_to_mapping )

		### Fit
		xi1D_data = numpy.zeros(shape=(self._nbBin1D,self._nb_wedges,3))
		xi1D_fit  = numpy.zeros(shape=(self._nbBin1D,self._nb_wedges,3))
		for i in numpy.arange(self._nbBinX2D):
			for j in numpy.arange(self._nbBinY2D):

				if (self._xi2D_fit[i,j,8]<=0.): continue
				ivar = 1./( self._xi2D_fit[i,j,8]*self._xi2D_fit[i,j,8] )
				for k in numpy.arange(self._nbBin1D):
					for l in numpy.arange(self._nb_wedges):

						coef = mapping_2D_to_we[i,j,k,l]
						if (coef==0.): continue
						xi1D_data[k,l,0] += coef*ivar*self._xi2D_fit[i,j,3]
						xi1D_data[k,l,1] += coef*ivar*self._xi2D_fit[i,j,7]
						xi1D_data[k,l,2] += coef*ivar

						xi1D_fit[k,l,0] += coef*ivar*self._xi2D_fit[i,j,3]
						xi1D_fit[k,l,1] += coef*ivar*self._xi2D_fit[i,j,6]
						xi1D_fit[k,l,2] += coef*ivar

		cut = (xi1D_data[:,:,2]>0.)
		xi1D_data[:,:,0][cut] /= xi1D_data[:,:,2][cut]
		xi1D_data[:,:,1][cut] /= xi1D_data[:,:,2][cut]
		xi1D_data[:,:,2][cut]  = 1./numpy.sqrt(xi1D_data[:,:,2][cut])
		cut = (xi1D_fit[:,:,2]>0.)
		xi1D_fit[:,:,0][cut] /= xi1D_fit[:,:,2][cut]
		xi1D_fit[:,:,1][cut] /= xi1D_fit[:,:,2][cut]
		xi1D_fit[:,:,2][cut]  = 1./numpy.sqrt(xi1D_fit[:,:,2][cut])

		return xi1D_data, xi1D_fit
	def plot_data_and_fit_1d(self,x_power,path_to_mapping_1D):

		xi1D_data, xi1D_fit = self.get_data_and_fit_1d(path_to_mapping_1D)

		### Data
		cut = (xi1D_data[:,2]>0.)
		coef = numpy.power(xi1D_data[:,0][cut],x_power)
		xxx = xi1D_data[:,0][cut]
		plt.errorbar(xi1D_data[:,0][cut], coef*xi1D_data[:,1][cut], yerr=coef*xi1D_data[:,2][cut], fmt='o', label=r'$'+self._name+'$', markersize=10,linewidth=2)
		### Fit
		cut = (xi1D_fit[:,2]>0.)
		coef = numpy.power(xi1D_fit[:,0][cut],x_power)
		plt.errorbar(xi1D_fit[:,0][cut], coef*xi1D_fit[:,1][cut], label=r'$Fit$', color='red',linewidth=2)

		if (x_power==0):
			plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,True)
		plt.show()
		
		return
	def plot_data_and_fit_we(self,x_power,path_to_mapping):

		color = ['blue', 'green', 'orange','cyan']

		xi1D_data, xi1D_fit = self.get_data_and_fit_we(path_to_mapping)
		'''
		cut = (xi1D_data[:,:,2]>0.)
		xi1D_data[:,:,0][cut] /= xi1D_data[:,:,2][cut]
		xi1D_data[:,:,1][cut] /= xi1D_data[:,:,2][cut]
		xi1D_data[:,:,2][cut]  = 1./numpy.sqrt(xi1D_data[:,:,2][cut])
		cut = (xi1D_fit[:,:,2]>0.)
		xi1D_fit[:,:,0][cut] /= xi1D_fit[:,:,2][cut]
		xi1D_fit[:,:,1][cut] /= xi1D_fit[:,:,2][cut]
		xi1D_fit[:,:,2][cut]  = 1./numpy.sqrt(xi1D_fit[:,:,2][cut])
		'''
		for i in numpy.arange(0,self._nb_wedges):
		
			cut = (self._xiWe[:,i,2]>0.)
			if (self._xiWe[:,i,0][cut].size==0):
				continue
		
			cut = (xi1D_data[:,i,2]>0.)
			xxx = xi1D_data[:,i,0][cut]
			yyy = xi1D_data[:,i,1][cut]
			yer = xi1D_data[:,i,2][cut]
			coef = numpy.power(xxx,x_power)
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=r'$'+self._label_wedge[i]+'$', color=color[i], markersize=10,linewidth=2)

			cut = (xi1D_fit[:,i,2]>0.)
			xxxF = xi1D_fit[:,i,0][cut]
			yyyF = xi1D_fit[:,i,1][cut]
			yerF = xi1D_fit[:,i,2][cut]
			coefF = numpy.power(xxxF,x_power)
			plt.errorbar(xxxF, coefF*yyyF, color='red',linewidth=2)

			#plt.errorbar(xxx[cut], coef[cut]*(yyy[cut]-yyyF)/yer[cut], fmt='o', markersize=10,linewidth=2, label=r'$'+label[i]+'$')

			if (x_power==0):
				plt.ylabel(r'$'+self._label+' (|s|)$', fontsize=40)
				plt.legend(fontsize=40, numpoints=1,ncol=2, loc=0) #4
			if (x_power==1):
				plt.ylabel(r'$|s|.'+self._label+' (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
				plt.legend(fontsize=40, numpoints=1,ncol=2, loc=0) #4
			if (x_power==2):
				plt.ylabel(r'$|s|^{2}.'+self._label+' (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
				plt.legend(fontsize=40, numpoints=1,ncol=2, loc=0) #2

		#plt.ylabel(r'$Residuals (|s|)$', fontsize=40)
		
		plt.title(r'$'+self._title+'$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(xxx)+10. ])
		myTools.deal_with_plot(False,False,False)
		plt.show()

		
		return
	def plot_given_2d(self, x_power=0, xi2D=None, label=None):

		if (xi2D is None):
			print '  annalyseBAOFIT::plot_2d::  xi2D==None'
			return
		if (label is None):
			label = self._label

		origin='lower'
		extent=[self._minX2D, self._maxX2D, self._minY2D, self._maxY2D]
		if (self._correlation=='q_f' or self._correlation=='f_f2'):
			origin='upper'
			extent=[self._minX2D, self._maxX2D, self._maxY2D, self._minY2D]
	
		xxx = numpy.transpose(xi2D[:,:,0])
		yyy = numpy.transpose(xi2D[:,:,1])
		yer = numpy.transpose(xi2D[:,:,2])
	
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
			cbar.set_label(r'$'+label+'(\, \overrightarrow{s} \,)$',size=40)
		if (x_power==1):
			cbar.set_label(r'$|s|.'+label+'(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
		if (x_power==2):
			cbar.set_label(r'$|s|^{2}.'+label+'(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)
	
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
		myTools.deal_with_plot(False,False,False)

		plt.show()
		
		return
	def plot_fit_2d(self,x_power=0):

		xi2D = numpy.zeros(shape=(self._nbBinX2D,self._nbBinY2D,3))
		xi2D[:,:,0] = self._xi2D_fit[:,:,3]
		xi2D[:,:,1] = self._xi2D_fit[:,:,6]
		xi2D[:,:,2] = self._xi2D_fit[:,:,8]
		self.plot_given_2d(x_power, xi2D, self._label+'_{fit}')

		return
	def plot_slice_fit_2d(self,sliceX=None,sliceY=None, other=[]):

		list_corr = [self] + other
		i = 0
		fit = True
		for el in list_corr:
			i += 1

			if (sliceX is not None):
				mean = numpy.mean(self._xi2D_grid[sliceX,:,0])
				xxx = el._xi2D_grid[sliceX,:,1]
				yyy = el._xi2D[sliceX,:,1]
				yer = el._xi2D[sliceX,:,2]
				try:
					xxx2 = el._xi2D_grid[sliceX,:,1]
					yyy2 = el._xi2D_fit[sliceX,:,6]
					yer2 = el._xi2D_fit[sliceX,:,8]
				except:
					print 'no fit'
					fit = False
			elif (sliceY is not None):
				mean = numpy.mean(self._xi2D_grid[:,sliceY,1])
				xxx = self._xi2D_grid[:,sliceY,0]
				yyy = self._xi2D[:,sliceY,1]
                	        yer = self._xi2D[:,sliceY,2]
				try:
					xxx2 = el._xi2D_grid[:,sliceY,0]
					yyy2 = el._xi2D_fit[:,sliceY,6]
					yer2 = el._xi2D_fit[:,sliceY,8]
				except:
					print 'no fit'
					fit = False
			else:
				return

			cut = (yer>0.)
			xxx = xxx[cut]
			yyy = yyy[cut]
			yer = yer[cut]
			if (xxx.size==0): return
			plt.errorbar(xxx, yyy, yerr=yer, fmt='o', markersize=10,linewidth=2)

			if (fit):
				cut = (yer2>0.)
				xxx2 = xxx2[cut]
				yyy2 = yyy2[cut]
				yer2 = yer2[cut]
				if (xxx2.size!=0): plt.errorbar(xxx2, yyy2, linewidth=2,color='red')
		
		if (sliceX is not None):
			plt.xlabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
			plt.title(r'$<s_{\perp}> = %.2f \, Mpc.h^{-1}$' % mean)
		else:
			plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
			plt.title(r'$<s_{\parallel}> = %.2f \, Mpc.h^{-1}$' % mean)
		plt.ylabel(r'$'+self._label+'(\, \overrightarrow{s} \,)$',fontsize=40)	
                myTools.deal_with_plot(False,False,False)
                plt.show()

		return
	def plot_residuals_2d(self,x_power=0):

		self.plot_given_2d(x_power, self.get_residuals(), self._label+'_{residuals}')

		return
	def plot_histo_residuals(self, nbBins=100):

		xi2D = self.get_residuals()
		yyy  = (xi2D[:,:,1][ (xi2D[:,:,2]>0.) ]).flatten()
	
		fig = plt.figure()
		ax = fig.add_subplot(111)
	
		ax.hist(yyy, bins=nbBins)
		plt.xlabel(r'$\frac{data-fit}{\sigma_{data}}$')
		plt.ylabel(r'$\#$')
		myTools.deal_with_plot(False,False,False)
		
		mng = plt.get_current_fig_manager()
		textstr = '$nb=%u$\n$\mu=%.5e$\n$\sigma=%.5e$'%(yyy.size, numpy.mean(yyy), numpy.std(yyy))
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=40, verticalalignment='top', bbox=props)
		
		plt.show()
	def plot_chi2_scan(self, sizeX=100, sizeY=100, edge=None, toyMC=False, bootstrap=False, simulation=None):
		
		### Read chi2 scan
		bestFit = self.read_chi2_scan(sizeX, sizeY)

		if (edge is not None):
			self._chi2_scan_edge = edge

		if (toyMC):
			bestFit = numpy.append( bestFit,self.read_toyMC(),axis=1 )
		if (bootstrap):
			bestFit = numpy.append( bestFit,self.read_bootStrap(),axis=1 )
		if (simulation is not None):
			bestFit = numpy.append( bestFit,simulation,axis=1 )

		### Plot the data
		myTools.plot_chi2_scan(self._chi2_scan,self._chi2_scan_edge,True, bestFit,'\\alpha_{\\perp}','\\alpha_{\\parallel}','\\Delta \\chi^{2} = \\chi^{2}-\\chi^{2}_{best \\, fit}')

		return





















































