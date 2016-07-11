# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_3D.py
#

### Python lib
import numpy
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import copy
import subprocess

### Perso lib
import myTools
import correlation_3D
import correlation_3D_Q
import annalyse_BAOFIT
import annalyse_BAOFIT_Q

class AnnalyseManyBAOFIT:

	def __init__(self, dic=None, index_parameter=None, dic_Q=None, dic_simu=None, load_correlation=True, isPyLyA=False):

		if (dic_simu==None):
			dic_simu = correlation_3D.raw_dic_simu

		self._isPyLyA      = isPyLyA
		self._nbBox        = dic_simu['nb_box']
		self._nbSimu       = dic_simu['nb_simu']
		self._nbtot        = self._nbBox*self._nbSimu
		self._path_to_simu = dic_simu['path_to_simu']
		self._name         = dic['name']
		self._param        = None
		self._list_chi2    = None
		self._minos        = None
		self._redshift     = numpy.zeros( self._nbtot )
		self._dic          = dic
		if (self._dic is None):
			self._dic = correlation_3D.raw_dic_class
		
		### Get all the fits
		self._listFit = []
		for i in range(0,self._nbBox):
			for j in range(0,self._nbSimu):

				path = self._path_to_simu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/'
				dic['path_to_txt_file_folder'] = path+'Results' + dic_simu['prefix']+'/'
				dic['name'] = str(i)+'-'+str(j)

				if (load_correlation or (i==0 and j==0)):

					if (True):
						if (dic['correlation']=='q_q'):
							self._listFit += [annalyse_BAOFIT_Q.AnnalyseBAOFIT(dic, index_parameter, None,dic_Q)]
						else:
							if (dic['correlation']=='q_f'):
								path_to_BAOFIT = dic['path_to_txt_file_folder'] + 'BaoFit_'+dic['correlation']+'__'+dic['f1']+'__'+dic['q1']+dic_simu['prefix2']+'/bao2D.'
								if (self._isPyLyA): path_to_BAOFIT = dic['path_to_txt_file_folder'] + 'pyLyA'+dic_simu['prefix2']+'/cross_alone.'
							elif (dic['correlation']=='f_f'):
								path_to_BAOFIT = dic['path_to_txt_file_folder'] + 'BaoFit_'+dic['correlation']+'__'+dic['f1']+dic_simu['prefix2']+'/bao2D.'
								if (self._isPyLyA): path_to_BAOFIT = dic['path_to_txt_file_folder'] + 'pyLyA'+dic_simu['prefix2']+'/auto_alone.'
							self._listFit += [annalyse_BAOFIT.AnnalyseBAOFIT(dic=dic, dic_simu=dic_simu, index_parameter=index_parameter, path_to_BAOFIT=path_to_BAOFIT, isPyLyA=self._isPyLyA)]

						if (load_correlation):
							### Get all the resulted parameters
							nb = i*10+j
							if (i==0 and j==0):
								nb_param = self._listFit[nb]._param[:,0].size
								self._param = numpy.zeros( shape=(self._nbtot,nb_param,2) )
								self._list_chi2 = numpy.zeros( shape=(self._nbtot,4) )
							self._param[nb,:,:]   = self._listFit[nb]._param[:,:]
							self._list_chi2[nb,:] = self._listFit[nb]._chi2
					#except:
					#	print '  ERROR n°0 : ', i, j, dic['path_to_txt_file_folder']
				if (not load_correlation):

					if (dic['correlation']=='q_q'):
						path_to_data = dic['path_to_txt_file_folder'] + 'BaoFit_'+dic['correlation']+'__'+dic['q1']
						path_to_BAOFIT = path_to_data + dic_simu['prefix2']+'/bao2D.'
					elif (dic['correlation']=='q_f'):
						path_to_data = dic['path_to_txt_file_folder'] + 'BaoFit_'+dic['correlation']+'__'+dic['f1']+'__'+dic['q1']
						path_to_BAOFIT = path_to_data + dic_simu['prefix2']+'/bao2D.'
						if (self._isPyLyA): path_to_BAOFIT = dic['path_to_txt_file_folder'] + 'pyLyA'+dic_simu['prefix2']+'/cross_alone.'
					elif (dic['correlation']=='f_f'):
						path_to_data = dic['path_to_txt_file_folder'] + 'BaoFit_'+dic['correlation']+'__'+dic['f1']
						path_to_BAOFIT = path_to_data + dic_simu['prefix2']+'/bao2D.'
						if (self._isPyLyA): path_to_BAOFIT = dic['path_to_txt_file_folder'] + 'pyLyA'+dic_simu['prefix2']+'/auto_alone.'
					path_to_data += '/bao2D.'
					try:
					#if (True):
						nbParam, param, chi2, minos = self.read_fit_data_only_parameter(path_to_BAOFIT)
						nb = i*10+j
						if (i==0 and j==0):
							self._param     = numpy.zeros( shape=(self._nbtot,nbParam,2) )
							self._list_chi2 = numpy.zeros( shape=(self._nbtot,4) )
							self._minos     = numpy.zeros( shape=(self._nbtot,2,4) )
						self._param[nb,:,:]   = param
						self._list_chi2[nb,:] = chi2
						self._minos[nb,:]     = minos

						### Read redshift
						#self._redshift[nb] = numpy.mean( numpy.loadtxt(path_to_data+'grid')[:,3])
					except:
						print '  ERROR n°1: ', i, j
						#print path_to_BAOFIT
						#command = 'clubatch \" time ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Program/pyLyA/bin/fit --config_file /home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_no_metals/pyLyA/cross_alone.ini \"'
						#command = 'clubatch \" time ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Program/pyLyA/bin/fit --config_file /home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/pyLyA/cross_alone.ini \"'
						#command = 'ls /home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/pyLyA_withMetalsTemplates/ '
						#subprocess.call(command, shell=True)
						#command = 'clubatch \" time ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Program/pyLyA/bin/fit --config_file /home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/pyLyA_withMetalsTemplates/cross_alone.ini \"'
						#print command
						#subprocess.call(command, shell=True)

		if (load_correlation):
			self._nbtot = len(self._listFit)

		### Init parameters
		self._nbBin1D  = self._listFit[0]._nbBin1D
		self._nbBin2D  = self._listFit[0]._nbBin2D
		self._nbBinX2D = self._listFit[0]._nbBinX2D
		self._nbBinY2D = self._listFit[0]._nbBinY2D
		self._nbBinM   = self._listFit[0]._nbBinM
		### Set attributes set after
		self._chi2_scan = None

		### Create a 'AnnalyseBAOFIT' object for mean_fit
		#dic['path_to_txt_file_folder'] = self._path_to_simu + 'Box_000/Simu_000/Results/'
		#dic['name'] = '<fit \, simu>'
		#self._mean_fit = annalyse_BAOFIT.AnnalyseBAOFIT(dic, index_parameter, path+'BaoFit_q_f_covFromFit/bao2D.')

		### Create a 'Correlation3D' object for mean_data
		dic['path_to_txt_file_folder'] = self._path_to_simu + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
		if (dic['correlation']=='q_q'):
			self._mean_data = correlation_3D_Q.Correlation3DQ(dic,dic_Q)
			self._mean_data._xiMul = self._mean_data.get_multipol(self._mean_data._xiMu)
		else:
			self._mean_data = correlation_3D.Correlation3D(dic)
			self._mean_data._xiMul = self._mean_data.get_multipol(self._mean_data._xiMu)

		return
	def read_fit_data_only_parameter(self, path_to_BAOFIT):

		### String names
		if not self._isPyLyA:
			residuals = 'residuals.dat'
			param     = 'save.pars'
			chi2      = 'fit.chisq'
		else:
			if (self._dic['correlation']=='q_q'):
				prefix = 'autoQSO'
			elif (self._dic['correlation']=='q_f'):
				prefix = 'cross'
			elif (self._dic['correlation']=='f_f'):
				prefix = 'auto'
			residuals = prefix+'_all_residuals.dat'
			param     = 'save.pars'
			chi2      = prefix+'_fit.chisq'

		### Get the parameters covariance of the fit
		#data    = numpy.loadtxt(path_to_BAOFIT+param_covar)

		### Get the parameters of the fit
		#data    = numpy.loadtxt(self._path_to_BAOFIT+param)
		if not self._isPyLyA:
			data    = numpy.loadtxt(path_to_BAOFIT+param)
		else:
			data    = numpy.loadtxt(path_to_BAOFIT+param,usecols=(0,2,3))

		nbParam = data[:,1].size
		param   = numpy.zeros( shape=(nbParam,2) )
		param[:,0] = data[:,1]
		param[:,1] = data[:,2]

		### Get the parameters of the fit
		data = numpy.loadtxt(path_to_BAOFIT+chi2)
		chi2 = data

		### Read minos value
		minos = numpy.zeros(shape=(2,4))
		try:
			data = numpy.loadtxt(path_to_BAOFIT+'minos_save.pars',skiprows=1,usecols=(5,1,2,3,4,6))
			minos[:,0] = data[:,0]
			minos[:,1] = data[:,1]
			minos[:,2] = data[:,2]
			minos[:,3] = data[:,3]*data[:,4]*data[:,5]
		except:
			print '  No minos errors', 

		return nbParam, param, chi2, minos
	def save_list_realisation(self):

		list1D       = numpy.zeros( shape=(self._nbBin1D,self._nbtot) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,self._nbtot) )
		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,self._nbtot) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges,self._nbtot) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,5,self._nbtot) )

		pathToSave = self._path_to_simu + 'Results/'

		for i in numpy.arange(self._nbtot):
			print i

			list1D[:,i]         = self._listFit[i]._xi1D[:,1]
			list2D[:,i]         = self._listFit[i]._xi2D[:,:,1].flatten()
			listMu[:,i]         = self._listFit[i]._xiMu[:,:,2].flatten()
			listWe[:,:,i]       = self._listFit[i]._xiWe[:,:,1]
			listMultipol[:,:,i] = self._listFit[i].get_multipol( self._listFit[i]._xiMu)[:,:,1]

		numpy.save(pathToSave+'list_Mu',listMu)
		numpy.save(pathToSave+'list_We',listWe)
		numpy.save(pathToSave+'list_1D',list1D)
		numpy.save(pathToSave+'list_2D',list2D)
		numpy.save(pathToSave+'list_Multipol',listMultipol)

		covMu = numpy.cov(listMu)
		cov1D = numpy.cov(list1D)
		cov2D = numpy.cov(list2D)

		numpy.save(pathToSave+'cov_Mu',covMu)
		numpy.save(pathToSave+'cov_1D',cov1D)
		numpy.save(pathToSave+'cov_2D',cov2D)
	
		return
	def set_mean_fit(self):

		list1D       = numpy.zeros( shape=(self._nbBin1D,self._nbtot) )
		list2D       = numpy.zeros( shape=(self._nbBin2D,self._nbtot) )
		listMu       = numpy.zeros( shape=(self._nbBin1D*self._nbBinM,self._nbtot) )
		listWe       = numpy.zeros( shape=(self._nbBin1D,self._nb_wedges,self._nbtot) )
		listMultipol = numpy.zeros( shape=(self._nbBin1D,5,self._nbtot) )

		return
	def set_mean_data(self):

		raw_path = self._path_to_simu + 'Results/'
		### we
		path = raw_path + 'list_We.npy'
		listWe = numpy.load(path)
		nbTot = listWe[0,0,:].size
		for i in numpy.arange(self._nb_wedges):
			self._mean_data._xiWe[:,i,1] = numpy.mean(listWe[:,i,:],axis=1)
			self._mean_data._xiWe[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listWe[:,i,:] ))/nbTot)
		### mu
		path = raw_path + 'cov_Mu.npy'
		self._mean_data._xiMu[:,:,2] = myTools.convert1DTo2D(numpy.mean(numpy.load(raw_path + 'list_Mu.npy'),axis=1),self._nbBin1D,self._nbBinM)
		self._mean_data._xiMu[:,:,3] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path))/nbTot ),self._nbBin1D,self._nbBinM)
		### 1D
		path = raw_path + 'cov_1D.npy'
		self._mean_data._xi1D[:,1]   = numpy.mean(numpy.load(raw_path + 'list_1D.npy'),axis=1)
		self._mean_data._xi1D[:,2]   = numpy.sqrt( numpy.diag(numpy.load(path))/nbTot )
		### 2D
		path = raw_path + 'cov_2D.npy'
		self._mean_data._xi2D[:,:,1] = myTools.convert1DTo2D(numpy.mean(numpy.load(raw_path + 'list_2D.npy'),axis=1),self._nbBinX2D, self._nbBinY2D)
		self._mean_data._xi2D[:,:,2] = myTools.convert1DTo2D(numpy.sqrt( numpy.diag(numpy.load(path))/nbTot ),self._nbBinX2D, self._nbBinY2D)
		### Multipol
		path = raw_path + 'list_Multipol.npy'
		listMultipol = numpy.load(path)
		for i in numpy.arange(5):
			self._mean_data._xiMul[:,i,1] = numpy.mean(listMultipol[:,i,:],axis=1)
			self._mean_data._xiMul[:,i,2] = numpy.sqrt( numpy.diag(numpy.cov( listMultipol[:,i,:] ))/nbTot)

		return
	def print_list_parameter(self,index):

		for i in numpy.arange(self._param[:,index,0].size):

                        yyy = self._param[i,index,0]
                        yer = self._param[i,index,1]
			print i, yyy, yer
		return
	def print_results(self,index, other=[]):

		list_of_fit = numpy.append( [self],other )
		color = ['blue', 'red', 'green', 'cyan', 'orange', 'black']

		for el in list_of_fit:

			yyy = el._param[:,index,0][ (el._param[:,index,1]>0.) ]
			yer = el._param[:,index,1][ (el._param[:,index,1]>0.) ]

			print 
			print '      ', el._name
			print
			print '  --  ', el._listFit[0]._par_name[index], '  --  '
			print '      mean = ', numpy.mean(yyy)
			print '      std  = ', numpy.std(yyy)
			print '      err  = ', numpy.std(yyy)/numpy.sqrt(el._nbtot)
			print
			values  = yyy
			weights = numpy.power(yer,-2.)
			average = numpy.average(values, weights=weights)
			print '    w mean = ', average
			print '    w std  = ', numpy.sqrt( numpy.average((values-average)**2, weights=weights))
			print '    w err  = ', numpy.sqrt( numpy.average((values-average)**2, weights=weights)/el._nbtot)
			print '  --  '
			print

		### mock index evolution
		nb = 0
		for el in list_of_fit:
			plt.errorbar(numpy.arange(el._nbtot), el._param[:,index,0], yerr=el._param[:,index,1], marker='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
			nb += 1
		plt.xlabel(r'$simulation \, index$')
		plt.ylabel(r'$'+self._listFit[0]._par_name[index]+'$')
		plt.xlim([-1, el._nbtot+1])
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()

		### histogram
		nb = 0
		max_hist = 0.
		for el in list_of_fit:
			yyy = el._param[:,index,0][ (el._param[:,index,1]>0.) ]
			a, b, c = plt.hist(yyy,bins=10,linewidth=2, histtype='step',alpha=0.6, color=color[nb])
			max_hist = numpy.max( [numpy.max(a), max_hist])
			nb += 1
		plt.ylim([0.,1.2*max_hist])
		nb = 0
		for el in list_of_fit:
			### xxx
			x       = el._param[:,index,0][ (el._param[:,index,1]>0.) ]
			weights = numpy.ones(el._param[:,index,1][ (el._param[:,index,1]>0.) ].size) #numpy.power(el._param[:,index,1][ (el._param[:,index,1]>0.) ], -2.)
			mean_x  = numpy.average(x, weights=weights)
			err_x   = numpy.sqrt( numpy.average((x-mean_x)**2, weights=weights)/x.size) 
			### xxx
			label = el._name+" \, <"+el._listFit[0]._par_name[index]+"> = %1.4f \pm %1.4f" % ( mean_x,err_x )
			plt.plot( [mean_x,mean_x], [0.,1.2*max_hist], linestyle='dashed',linewidth=2, label=r'$'+label+'$',alpha=0.6, color=color[nb])
			nb += 1
		plt.xlabel(r'$'+self._listFit[0]._par_name[index]+'$')
		plt.ylabel(r'$\#$')
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()

		### histogram of poll plot
		nb = 0
		max_hist = 0.
		for el in list_of_fit:
			y       = el._param[:,index,0][ (el._param[:,index,1]>0.) ]
			weights = numpy.power(el._param[:,index,1][ (el._param[:,index,1]>0.) ], -2.)
			mean_y  = numpy.average(y, weights=weights)
			yyy     = (y-mean_y)/el._param[:,index,1][ (el._param[:,index,1]>0.) ]
			sigma   = numpy.average((yyy)**2, weights=weights)
			label   = el._name+" \, : \, \sigma = %1.2f" % ( sigma )
			a, b, c = plt.hist(yyy,bins=10,linewidth=2, histtype='step',alpha=0.6, color=color[nb], label=r'$'+label+'$')
			max_hist = numpy.max( [numpy.max(a), max_hist])
			nb += 1
		plt.ylim([0.,1.2*max_hist])
		plt.xlabel(r'$('+self._listFit[0]._par_name[index]+' - <'+self._listFit[0]._par_name[index]+'>)/\sigma$')
		plt.ylabel(r'$\#$')
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()

		
		### correlation between the versions
		nb = 0
		if (list_of_fit.size!=1):
			for el in list_of_fit[1:]:
				yyy1 = self._param[:,index,0][ (self._param[:,index,1]>0.) ]
				yyy2 = el._param[:,index,0][ (el._param[:,index,1]>0.) ]
				if (yyy1.size != yyy2.size): continue
				plt.errorbar(yyy1,yyy2, fmt='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
				nb += 1
			plt.xlabel(r'$'+self._listFit[0]._par_name[index]+' \, ' +self._name+'$')
			plt.ylabel(r'$'+self._listFit[0]._par_name[index]+' \, version \, i$')
			myTools.deal_with_plot(False,False,True)
			plt.legend(fontsize=20, numpoints=2,ncol=1)
			plt.show()
		'''
		### correlation with redshift
		nb = 0
		for el in list_of_fit:
			yyy1 = el._redshift[ (el._redshift!=0.) ]
			yyy2 = el._param[:,index,0][ (el._param[:,index,1]>0.) ]
			plt.errorbar(yyy1,yyy2, fmt='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
			nb += 1
		plt.xlabel(r'$<z>$')
		plt.ylabel(r'$'+self._listFit[0]._par_name[index]+'$')
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()
		'''

		return
	def print_results_for_chi2(self,other=[]):

		list_of_fit = numpy.append( [self],other )
		color = ['blue', 'red', 'green', 'cyan', 'orange', 'black']

		for el in list_of_fit:

			yyy = el._list_chi2[:,0][ (el._list_chi2[:,0]>0.) ]

			print 
			print '      chi^2'
			print
			print '      mean = ', numpy.mean(yyy)
			print '      std  = ', numpy.std(yyy)
			print '      err  = ', numpy.std(yyy)/numpy.sqrt(el._nbtot)
			print

		### mock index evolution
		nb = 0
		for el in list_of_fit:
			plt.errorbar(numpy.arange(el._nbtot), el._list_chi2[:,0], marker='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
			nb += 1
		plt.xlabel(r'$simulation \, index$')
		plt.ylabel(r'$\chi^{2} \, ('+str(int(el._list_chi2[0,1]))+'-'+str(int(el._list_chi2[0,2]))+')$')
		plt.xlim([-1, el._nbtot+1])
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()

		### histogram
		nb = 0
		max_hist = 0.
		for el in list_of_fit:
			yyy = el._list_chi2[:,0][ (el._list_chi2[:,0]>0.) ]
			a, b, c = plt.hist(yyy,bins=10,linewidth=2, histtype='step',alpha=0.6, color=color[nb])
			max_hist = numpy.max( [numpy.max(a), max_hist])
			nb += 1
		plt.ylim([0.,1.2*max_hist])
		nb = 0
		for el in list_of_fit:
			### xxx
			x       = el._list_chi2[:,0][ (el._list_chi2[:,0]>0.) ]
			weights = numpy.ones(x.size)
			mean_x  = numpy.average(x, weights=weights)
			err_x   = numpy.sqrt( numpy.average((x-mean_x)**2, weights=weights)/x.size) 
			### xxx
			label = '$' + el._name+" $ \n $ <\chi^{2}> / NDF = %1.4f \, (%u - %u) $ \n $ \sigma = %1.2f $" % ( mean_x,int(el._list_chi2[0,1]),int(el._list_chi2[0,2]), err_x*numpy.sqrt(x.size)  )
			plt.plot( [mean_x,mean_x], [0.,1.2*max_hist], linestyle='dashed',linewidth=2, label=r''+label,alpha=0.6, color=color[nb])
			nb += 1
		#plt.xlabel(r'$\chi^{2} \, ('+str(int(el._list_chi2[0,1]))+'-'+str(int(el._list_chi2[0,2]))+')$')
		plt.xlabel(r'$\chi^{2} $')
		plt.ylabel(r'$\#$')
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()
		
		### correlation between the versions
		nb = 0
		if (list_of_fit.size!=1):
			for el in list_of_fit[1:]:
				yyy1 = self._list_chi2[:,0][ (self._list_chi2[:,0]>0.) ]
				yyy2 = el._list_chi2[:,0][ (el._list_chi2[:,0]>0.) ]
				if (yyy1.size != yyy2.size): continue
				plt.errorbar(yyy1,yyy2, fmt='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
				nb += 1
			plt.xlabel(r'$\chi^{2} \, ('+str(int(el._list_chi2[0,1]))+'-'+str(int(el._list_chi2[0,2]))+') \, ' + self._name+'$')
			plt.ylabel(r'$\chi^{2} \, ('+str(int(el._list_chi2[0,1]))+'-'+str(int(el._list_chi2[0,2]))+') \, version \, i$')
			myTools.deal_with_plot(False,False,True)
			plt.legend(fontsize=20, numpoints=2,ncol=1)
			plt.show()

		### correlation with redshift
		nb = 0
		for el in list_of_fit:
			yyy1 = el._redshift[ (el._redshift!=0.) ]
			yyy2 = el._list_chi2[:,0][ (el._list_chi2[:,0]>0.) ]
			plt.errorbar(yyy1,yyy2, fmt='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
			nb += 1
		plt.xlabel(r'$<z>$')
		plt.ylabel(r'$\chi^{2} \, ('+str(int(el._list_chi2[0,1]))+'-'+str(int(el._list_chi2[0,2]))+')$')
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()
		

		return
	def print_poll_plots(self,index, other=[]):

		list_of_fit = numpy.append( [self],other )
		color = ['blue', 'red', 'green', 'cyan', 'orange', 'black']	

		### mock index evolution
		nb = 0
		for el in list_of_fit:
			xxx = numpy.arange(el._nbtot)
			yyy = el._minos[:,index,0]
			yer1 = el._minos[:,index,1]
			yer2 = el._minos[:,index,2]
			yyy[(el._minos[:,index,3]==0.)] = 0.
			yer1[(el._minos[:,index,3]==0.)] = 0.
			yer2[(el._minos[:,index,3]==0.)] = 0.
			yer = [ numpy.abs(yer1),yer2 ]
			plt.errorbar(xxx, yyy, yerr=yer, marker='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
			nb += 1
		plt.xlabel(r'$simulation \, index$')
		plt.ylabel(r'$'+self._listFit[0]._par_name[index]+'$')
		plt.xlim([-1, el._nbtot+1])
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()

		### mock index evolution for poll
		nb = 0
		for el in list_of_fit:
			xxx = numpy.arange(el._nbtot)
			yyy = el._minos[:,index,0]
			yer = numpy.abs(el._minos[:,index,1])
			yer[yyy-1.>0.] = el._minos[:,index,2][yyy-1.>0.]
			yer[(el._minos[:,index,3]==0.)] = 0.
			yyy   = (yyy-1.)/yer
			plt.errorbar(xxx[(yer!=0.)], yyy[(yer!=0.)], marker='o', markersize=10,linewidth=2, label=r'$'+el._name+'$',alpha=0.6, color=color[nb])
			nb += 1
		plt.xlabel(r'$simulation \, index$')
		plt.ylabel(r'$'+self._listFit[0]._par_name[index]+'$')
		plt.xlim([-1, el._nbtot+1])
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()

		### histogram of poll plot
		nb = 0
		max_hist = 0.
		for el in list_of_fit:

			yyy = el._minos[:,index,0]
			yer = numpy.abs(el._minos[:,index,1])
			yer[yyy-1.>0.] = el._minos[:,index,2][yyy-1.>0.]

			cut = numpy.abs(el._minos[:,index,1])>0.
			yyy = yyy[cut]
			yer = yer[cut]

			www   = numpy.power(yer, -2.)
			mmy   = numpy.average(yyy, weights=www)
			yyy   = (yyy-1.)/yer
			sigma = numpy.average((yyy)**2, weights=www)

			label   = el._name+" \, : \, \sigma = %1.2f" % ( sigma )
			a, b, c = plt.hist(yyy,bins=10,linewidth=2, histtype='step',alpha=0.6, color=color[nb], label=r'$'+label+'$')
			max_hist = numpy.max( [numpy.max(a), max_hist])
			nb += 1

		plt.ylim([0.,1.2*max_hist])
		plt.xlabel(r'$('+self._listFit[0]._par_name[index]+' - <'+self._listFit[0]._par_name[index]+'>)/\sigma$')
		plt.ylabel(r'$\#$')
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=20, numpoints=2,ncol=1)
		plt.show()
		

		return
	def get_best_fit(self):

		bestFit = numpy.zeros( shape=(2,1) )
		for el in self._listFit:
			bestFit = numpy.append(bestFit,el.get_best_fit(),axis=1)

		return bestFit
	def plot_1d(self, x_power=0, other=[], title=True):
		self._mean_data.plot_1d(x_power, other, title)
		
		return
	def plot_we(self, x_power=0, other=[], title=True):
		self._mean_data.plot_we(x_power, other, title)
		
		return
	def plot_2d(self, x_power=0):
		self._mean_data.plot_2d(x_power)
		
		return
	def plot_mu(self, x_power=0):
		self._mean_data.plot_mu(x_power)
	
		return
	def plot_multipol(self, x_power=0):
		self._mean_data.plot_multipol(x_power)
	
		return
	def plot_histo_residuals(self, other=[]):

		### Constants
		nbBins=100

		list_of_fit = numpy.append( [self],other )

		### histo
		fig = plt.figure()
		ax = fig.add_subplot(111)

		for el in list_of_fit:
			yyy = numpy.array([])
			for el in self._listFit:
				xi2D = el.get_residuals()
				tmp_yyy  = (xi2D[:,:,1][ (xi2D[:,:,2]>0.) ]).flatten()
				yyy = numpy.append(yyy,tmp_yyy)

			ax.hist(yyy, bins=nbBins, histtype='step', linewidth=2,alpha=0.6)

		plt.xlabel(r'$\frac{data-fit}{\sigma_{data}}$')
		plt.ylabel(r'$\#$')
		myTools.deal_with_plot(False,False,False)
		
		mng = plt.get_current_fig_manager()
		textstr = '$nb=%u$\n$\mu=%.5e$\n$\sigma=%.5e$'%(yyy.size, numpy.mean(yyy), numpy.std(yyy))
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)
		
		plt.show()

		return
	def plot_scatter_hist(self,index_1, index_2, other=None, single_point=None):

		list_of_fit = numpy.array( [self] )
		if (other is not None): list_of_fit = numpy.append( list_of_fit,other )

		### Constants
		nb_bins = 10
		color = ['blue', 'red', 'green', 'cyan', 'orange', 'black']

		nullfmt = NullFormatter()         # no labels
		
		# definitions for the axes
		left, width = 0.1, 0.65
		bottom, height = 0.1, 0.65
		bottom_h = left_h = left + width + 0.02
		
		rect_scatter = [left, bottom, width, height]
		rect_histx = [left, bottom_h, width, 0.2]
		rect_histy = [left_h, bottom, 0.2, height]
		
		# start with a rectangular Figure
		plt.figure(1, figsize=(8, 8))
		
		###
		axScatter = plt.axes(rect_scatter)
		plt.xlabel(r'$'+self._listFit[0]._par_name[index_1]+'$')
		plt.ylabel(r'$'+self._listFit[0]._par_name[index_2]+'$')
		plt.rc('font', **{'size':'40'})
		plt.tick_params(axis='both', which='major', labelsize=40)
		plt.grid(True, which='both')
		plt.linewidth = 40
		###
		axHistx = plt.axes(rect_histx)
		plt.rc('font', **{'size':'40'})
		plt.tick_params(axis='both', which='major', labelsize=40)
		plt.grid(True, which='both')
		plt.linewidth = 40
		###
		axHisty = plt.axes(rect_histy)
		plt.rc('font', **{'size':'40'})
		plt.tick_params(axis='both', which='major', labelsize=40)
		plt.grid(True, which='both')
		plt.linewidth = 40

		# no labels
		axHistx.xaxis.set_major_formatter(nullfmt)
		axHisty.yaxis.set_major_formatter(nullfmt)		


		nb = 0
		max_hist_x = 0.
		max_hist_y = 0.
		for el in list_of_fit:

			### xxx
			cut = numpy.logical_and(el._param[:,index_1,1]>0., el._param[:,index_2,1]>0.)
			#cut = numpy.logical_and(el._minos[:,1,3]>0., numpy.logical_and(el._minos[:,0,3]>0.,numpy.logical_and(el._param[:,index_1,1]>0., el._param[:,index_2,1]>0.)))
			x     = el._param[:,index_1,0][cut]
			if (x.size==0): return
			x_err = el._param[:,index_1,1][cut]
			weights = numpy.ones(x_err.size) #numpy.power(x_err,-2.)
			mean_x = numpy.average(x, weights=weights)
			### yyy
			y     = el._param[:,index_2,0][cut]
			if (y.size==0): return
			y_err = el._param[:,index_2,1][cut]
			weights = numpy.ones(y_err.size) #numpy.power(y_err,-2.)
			mean_y = numpy.average(y, weights=weights)
			if (x.size==0 or y.size==0 or x.size!=y.size):
				print ' ERROR ! x.size==0 or y.size==0 or x.size!=y.size'
				plt.clf()
				return
			else: print ' size is = ', x.size

			# the scatter plot:
			axScatter.errorbar(x, y, fmt='o',linewidth=2, markersize=10, color=color[nb],alpha=0.6)

			# now determine nice limits by hand:
			limit = numpy.zeros(4)
			mean  = numpy.zeros(2)
			mean[0]  = numpy.mean(x)
			limit[0] = numpy.min(x)
			limit[1] = numpy.max(x)
			mean[1]  = numpy.mean(y)
			limit[2] = numpy.min(y)
			limit[3] = numpy.max(y)

			limit_for_hist = numpy.array(limit)

			
			for i in [0,2]:
				if ( limit_for_hist[i]<0. ): limit_for_hist[i]*= 1.001
				elif ( limit_for_hist[i]>0. ): limit_for_hist[i]*= 0.999
			for i in [1,3]:
				if ( limit_for_hist[i]<0. ): limit_for_hist[i]*= 0.999
				elif ( limit_for_hist[i]>0. ): limit_for_hist[i]*= 1.001
			

			###
			bins = numpy.arange( limit_for_hist[0],limit_for_hist[1]+(limit_for_hist[1]-limit_for_hist[0])/nb_bins, (limit_for_hist[1]-limit_for_hist[0])/nb_bins )
			a, b, c = axHistx.hist(x, bins=bins, histtype='step',linewidth=2, color=color[nb],alpha=0.6)
			max_hist_x = numpy.max( [numpy.max(a), max_hist_x])

			###
			bins = numpy.arange( limit_for_hist[2],limit_for_hist[3]+(limit_for_hist[3]-limit_for_hist[2])/nb_bins, (limit_for_hist[3]-limit_for_hist[2])/nb_bins )
			a, b, c = axHisty.hist(y, bins=bins, orientation='horizontal', histtype='step',linewidth=2, color=color[nb],alpha=0.6)
			max_hist_y = numpy.max( [numpy.max(a), max_hist_y])

			# now determine nice limits by hand:
			for i in [0,2]:
				limit[i] -= (mean[i/2]-limit[i])*0.1
			for i in [1,3]:
				limit[i] += (limit[i]-mean[i/2])*0.1

			if (nb==0):
				full_limit = numpy.array(limit)
			else:
				for i in [0,2]:
					full_limit[i] = numpy.min( [full_limit[i],limit[i]] )
				for i in [1,3]:
					full_limit[i] = numpy.max( [full_limit[i],limit[i]] )

			nb += 1

		### Plot single points
		if (single_point is not None):
			for el in single_point:
				axScatter.errorbar(el[0], el[1], fmt='*',linewidth=2, markersize=30, color=color[nb],alpha=0.9)
			nb += 1


		### Set mean line in histogram
		nb = 0
		for el in list_of_fit:

			### xxx
			cut = numpy.logical_and(el._param[:,index_1,1]>0., el._param[:,index_2,1]>0.)
			#cut = numpy.logical_and(el._minos[:,1,3]>0., numpy.logical_and(el._minos[:,0,3]>0.,numpy.logical_and(el._param[:,index_1,1]>0., el._param[:,index_2,1]>0.)))
			x       = el._param[:,index_1,0][cut]
			x_err   = el._param[:,index_1,1][cut]
			weights = numpy.ones(x_err.size) #numpy.power(x_err,-2.)
			mean_x  = numpy.average(x, weights=weights)
			err_x   = numpy.sqrt( numpy.average((x-mean_x)**2, weights=weights)/x.size) 
			### yyy
			y     = el._param[:,index_2,0][cut]
			y_err = el._param[:,index_2,1][cut]
			weights = numpy.ones(y_err.size) #numpy.power(y_err,-2.)
			mean_y  = numpy.average(y, weights=weights)
			err_y   = numpy.sqrt( numpy.average((y-mean_y)**2, weights=weights)/x.size) 

			### scatter
			axScatter.set_xlim( [full_limit[0],full_limit[1]] )
			axScatter.set_ylim( [full_limit[2],full_limit[3]] )
			axScatter.plot( [mean_x,mean_x], [full_limit[2],full_limit[3]],color=color[nb], linestyle='dashed',linewidth=2,alpha=0.6)
			axScatter.plot( [full_limit[0],full_limit[1]], [mean_y,mean_y],color=color[nb], linestyle='dashed',linewidth=2,alpha=0.6)
			### xxx
			axHistx.set_xlim( [full_limit[0],full_limit[1]] )
			axHistx.set_ylim([0.,1.1*max_hist_x])
			val = myTools.format_number_with_precision(mean_x,err_x)
                        err = myTools.format_number_with_precision(err_x,err_x)
			label = el._name+' \, <'+self._listFit[0]._par_name[index_1]+"> = "+val+" \pm" +err
			axHistx.plot( [mean_x,mean_x], [0.,1.1*max_hist_x],color=color[nb], linestyle='dashed',linewidth=2, label=r'$'+label+'$',alpha=0.6)
			axHistx.legend(fontsize=20, numpoints=2,ncol=1)
			
			### yyy
			axHisty.set_ylim( [full_limit[2],full_limit[3]] )
			axHisty.set_xlim([0.,1.1*max_hist_y])
			val = myTools.format_number_with_precision(mean_y,err_y)
			err = myTools.format_number_with_precision(err_y,err_y)
			label = el._name+' \, <'+self._listFit[0]._par_name[index_2]+"> = "+val+" \pm" +err
			axHisty.plot( [0.,1.1*max_hist_y], [mean_y,mean_y],color=color[nb], linestyle='dashed',linewidth=2, label=r'$'+label+'$',alpha=0.6)
			axHisty.legend(fontsize=20, numpoints=2,ncol=1)

			nb += 1

		plt.show()

		return
	def plot_scatter_hist_two_correlation(self,other,index_1, index_2):

		print self._listFit[0]._par_name[index_1], " x " , other._listFit[0]._par_name[index_2]
		### Constants
		nb_bins = 10
		color = ['blue', 'red', 'green', 'cyan', 'orange', 'black']

		nullfmt = NullFormatter() # no labels
		
		# definitions for the axes
		left, width = 0.1, 0.65
		bottom, height = 0.1, 0.65
		bottom_h = left_h = left + width + 0.02
		
		rect_scatter = [left, bottom, width, height]
		rect_histx = [left, bottom_h, width, 0.2]
		rect_histy = [left_h, bottom, 0.2, height]
		
		# start with a rectangular Figure
		plt.figure(1, figsize=(8, 8))
		
		###
		axScatter = plt.axes(rect_scatter)
		plt.xlabel(r'$'+self._listFit[0]._par_name[index_1]+'$')
		plt.ylabel(r'$'+other._listFit[0]._par_name[index_2]+'$')
		plt.rc('font', **{'size':'40'})
		plt.tick_params(axis='both', which='major', labelsize=40)
		plt.grid(True, which='both')
		plt.linewidth = 40
		###
		axHistx = plt.axes(rect_histx)
		plt.rc('font', **{'size':'40'})
		plt.tick_params(axis='both', which='major', labelsize=40)
		plt.grid(True, which='both')
		plt.linewidth = 40
		###
		axHisty = plt.axes(rect_histy)
		plt.rc('font', **{'size':'40'})
		plt.tick_params(axis='both', which='major', labelsize=40)
		plt.grid(True, which='both')
		plt.linewidth = 40

		# no labels
		axHistx.xaxis.set_major_formatter(nullfmt)
		axHisty.yaxis.set_major_formatter(nullfmt)		


		nb = 0
		max_hist_x = 0.
		max_hist_y = 0.

		### xxx
		x     = self._param[:,index_1,0][ (self._param[:,index_1,1]>0.) ]
		x_err = self._param[:,index_1,1][ (self._param[:,index_1,1]>0.) ]
		weights = numpy.ones(x_err.size) #numpy.power(x_err,-2.)
		mean_x = numpy.average(x, weights=weights)
		### yyy
		y     = other._param[:,index_2,0][ (other._param[:,index_2,1]>0.) ]
		y_err = other._param[:,index_2,1][ (other._param[:,index_2,1]>0.) ]
		weights = numpy.ones(y_err.size) #numpy.power(y_err,-2.)
		mean_y = numpy.average(y, weights=weights)

		print '  Correlation coef = ', numpy.corrcoef(x,y)[1,0]

		# the scatter plot:
		axScatter.errorbar(x, y, fmt='o',linewidth=2, markersize=10, color=color[nb],alpha=0.6)

		# now determine nice limits by hand:
		limit = numpy.zeros(4)
		mean  = numpy.zeros(2)
		mean[0]  = numpy.mean(x)
		limit[0] = numpy.min(x)
		limit[1] = numpy.max(x)
		mean[1]  = numpy.mean(y)
		limit[2] = numpy.min(y)
		limit[3] = numpy.max(y)

		limit_for_hist = numpy.array(limit)

		
		for i in [0,2]:
			if ( limit_for_hist[i]<0. ): limit_for_hist[i]*= 1.001
			elif ( limit_for_hist[i]>0. ): limit_for_hist[i]*= 0.999
		for i in [1,3]:
			if ( limit_for_hist[i]<0. ): limit_for_hist[i]*= 0.999
			elif ( limit_for_hist[i]>0. ): limit_for_hist[i]*= 1.001
		

		###
		bins = numpy.arange( limit_for_hist[0],limit_for_hist[1]+(limit_for_hist[1]-limit_for_hist[0])/nb_bins, (limit_for_hist[1]-limit_for_hist[0])/nb_bins )
		a, b, c = axHistx.hist(x, bins=bins, histtype='step',linewidth=2, color=color[nb],alpha=0.6)
		max_hist_x = numpy.max( [numpy.max(a), max_hist_x])

		###
		bins = numpy.arange( limit_for_hist[2],limit_for_hist[3]+(limit_for_hist[3]-limit_for_hist[2])/nb_bins, (limit_for_hist[3]-limit_for_hist[2])/nb_bins )
		a, b, c = axHisty.hist(y, bins=bins, orientation='horizontal', histtype='step',linewidth=2, color=color[nb],alpha=0.6)
		max_hist_y = numpy.max( [numpy.max(a), max_hist_y])

		# now determine nice limits by hand:
		for i in [0,2]:
			limit[i] -= (mean[i/2]-limit[i])*0.1
		for i in [1,3]:
			limit[i] += (limit[i]-mean[i/2])*0.1

		if (nb==0):
			full_limit = numpy.array(limit)
		else:
			for i in [0,2]:
				full_limit[i] = numpy.min( [full_limit[i],limit[i]] )
			for i in [1,3]:
				full_limit[i] = numpy.max( [full_limit[i],limit[i]] )

		### Set mean line in histogram

		### xxx
		x       = self._param[:,index_1,0][ (self._param[:,index_1,1]>0.) ]
		x_err   = self._param[:,index_1,1][ (self._param[:,index_1,1]>0.) ]
		weights = numpy.ones(x_err.size) #numpy.power(x_err,-2.)
		mean_x  = numpy.average(x, weights=weights)
		err_x   = numpy.sqrt( numpy.average((x-mean_x)**2, weights=weights)/x.size) 
		### yyy
		y     = other._param[:,index_2,0][ (other._param[:,index_2,1]>0.) ]
		y_err = other._param[:,index_2,1][ (other._param[:,index_2,1]>0.) ]
		weights = numpy.ones(y_err.size) #numpy.power(y_err,-2.)
		mean_y  = numpy.average(y, weights=weights)
		err_y   = numpy.sqrt( numpy.average((y-mean_y)**2, weights=weights)/x.size) 

		### scatter
		axScatter.set_xlim( [full_limit[0],full_limit[1]] )
		axScatter.set_ylim( [full_limit[2],full_limit[3]] )
		axScatter.plot( [mean_x,mean_x], [full_limit[2],full_limit[3]],color=color[nb], linestyle='dashed',linewidth=2,alpha=0.6)
		axScatter.plot( [full_limit[0],full_limit[1]], [mean_y,mean_y],color=color[nb], linestyle='dashed',linewidth=2,alpha=0.6)
		### xxx
		axHistx.set_xlim( [full_limit[0],full_limit[1]] )
		axHistx.set_ylim([0.,1.1*max_hist_x])
		val = myTools.format_number_with_precision(mean_x,err_x)
		err = myTools.format_number_with_precision(err_x,err_x)
		label = self._name+' \, <'+self._listFit[0]._par_name[index_1]+"> = "+val+" \pm" +err
		axHistx.plot( [mean_x,mean_x], [0.,1.1*max_hist_x],color=color[nb], linestyle='dashed',linewidth=2, label=r'$'+label+'$',alpha=0.6)
		axHistx.legend(fontsize=20, numpoints=2,ncol=1)
		
		### yyy
		axHisty.set_ylim( [full_limit[2],full_limit[3]] )
		axHisty.set_xlim([0.,1.1*max_hist_y])
		val = myTools.format_number_with_precision(mean_y,err_y)
		err = myTools.format_number_with_precision(err_y,err_y)
		label = other._name+' \, <'+other._listFit[0]._par_name[index_2]+"> = "+val+" \pm" +err
		axHisty.plot( [0.,1.1*max_hist_y], [mean_y,mean_y],color=color[nb], linestyle='dashed',linewidth=2, label=r'$'+label+'$',alpha=0.6)
		axHisty.legend(fontsize=20, numpoints=2,ncol=1)

		plt.show()

		return
	def plot_chi2_scan(self, sizeX=100, sizeY=100, edge=None, bestFitData=None):

		self._chi2_scan = numpy.zeros( shape=(sizeX,sizeY) )
		nb = numpy.zeros( shape=(sizeX,sizeY) )

		if (bestFitData is None): bestFit = numpy.zeros( shape=(2,1) )
		else: bestFit = bestFitData

		for el in self._listFit:
			bestFit = numpy.append(bestFit,el.get_best_fit(),axis=1)
			el.read_chi2_scan(sizeX, sizeY)
			cut = numpy.isfinite(el._chi2_scan)
			self._chi2_scan[cut] += el._chi2_scan[cut]
			nb[cut] += 1.

		self._chi2_scan /= nb

		### Plot the data
		myTools.plot_chi2_scan(self._chi2_scan,self._listFit[0]._chi2_scan_edge,True, bestFit,'\\alpha_{\\perp}','\\alpha_{\\parallel}','\\Delta \\chi^{2} = \\chi^{2}-\\chi^{2}_{best \\, fit}')

		return






































