# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#

### Python lib
import subprocess
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit
import sys
import copy

### Perso lib
import const_delta
import myTools
import correlation_3D
import correlation_3D_Q
import annalyse_BAOFIT
import annalyse_many_BAOFIT
import CAMB
import matplotlib.pyplot as plt
import matplotlib
import const

path_to_where_to_plot = const.path_to_where_to_plot__ + '/data_and_simulation/all_realisations/'
path_to_where_to_plot = const.path_to_where_to_plot__ + '/Data/Various/'


def plot_data_simulation():
	"""
		get similar plots as Delubac et al. for the wedges
	"""

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'nb_wedges'               : 4,
		'index_multipole_max' : 4,
		'remove_residuales' : 'lambda_OBS',
		'path_to_txt_file_folder': '',
		'correlation': 'f_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'NOTHING'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '',
		'prefix2'   : '',
	}


	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
	data = correlation_3D.Correlation3D(dic_class)
	#data.save_list_realisation('subsampling', 80)
	data.set_error_on_covar_matrix('subsampling')

	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
	simulation = correlation_3D.Correlation3D(dic_class)
	#data.save_list_realisation_simulation(dic_class, dic_simu)
	simulation.set_values_on_mean_simulation(dic_simu)
	### Get path to the list of simulations results
	path_to_load = dic_simu['path_to_simu'] + 'Results'
	path_to_load += dic_simu['prefix']
	path_to_load += '/'
	if (simulation._prefix=='xi_QSO_QSO'):
		path_to_load += simulation._prefix + '_result_'
	else:
		path_to_load += simulation._prefix + '_' + simulation._middlefix + '_result_'
	print path_to_load


	def plot_1d(x_power):

		list_simulation = numpy.load(path_to_load+'list_1D.npy')
		nb_realisation = list_simulation[0,:].size
		print nb_realisation

		for i in numpy.arange(0,nb_realisation):
			TMP_xxx = simulation._xi1D[:,0]
			TMP_yyy = list_simulation[:,i]
			TMP_yer = simulation._xi1D[:,2]	
			TMP_cut = (TMP_yer>0.)
			TMP_xxx = TMP_xxx[ TMP_cut ]
			TMP_yyy = TMP_yyy[ TMP_cut ]
			TMP_yer = TMP_yer[ TMP_cut ]
			TMP_coef = numpy.power(TMP_xxx,x_power)
			plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy, linewidth=1, color='gray',alpha=0.3)

		### Mean simulation and 1 sigma variance
		TMP_xxx = simulation._xi1D[:,0]
		TMP_yyy = simulation._xi1D[:,1]
		TMP_yer = simulation._xi1D[:,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_yyy_minus_one_sigma = TMP_yyy - TMP_yer*numpy.sqrt(nb_realisation)
		TMP_yyy_plus_one_sigma  = TMP_yyy + TMP_yer*numpy.sqrt(nb_realisation)
		TMP_coef = numpy.power(TMP_xxx,x_power)
		plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy, markersize=10,linewidth=2, color='blue')
		plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy_minus_one_sigma, linewidth=2, color='blue',linestyle='dashed')
		plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy_plus_one_sigma, linewidth=2, color='blue',linestyle='dashed')

		### Data
		TMP_xxx = data._xi1D[:,0]
		TMP_yyy = data._xi1D[:,1]
		TMP_yer = data._xi1D[:,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='red', ecolor='black',linestyle='None')

		###
		xxx = data._xi1D[:,0]
		xxx = xxx[ data._xi1D[:,2]>0. ]

		if (x_power==0):
			plt.ylabel(r'$'+data._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s| \cdot '+data._label+' (|s|) \, [$h$^{-1}$Mpc$]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2} \cdot '+data._label+' (|s|) \, [($h$^{-1}$Mpc$)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [$h$^{-1}$Mpc$]$', fontsize=40)
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlim([ data._min1D,data._max1D ])
		myTools.deal_with_plot(False,False,True)

		name = path_to_where_to_plot + 'xi_' + data._correlation + '_1D_rescale_index_' + str(x_power) + '.png'
		plt.savefig(name, bbox_inches='tight')
		plt.clf()
		#plt.show()
		
		return
	def plot_we(we_index,x_power):

		fig = plt.figure()
		ax = fig.add_subplot(111)

		list_simulation = numpy.load(path_to_load+'list_We.npy')
		nb_realisation = list_simulation[0,0,:].size
		print nb_realisation
		
		for i in numpy.arange(0,nb_realisation):
			TMP_xxx = simulation._xiWe[:,we_index,0]
			TMP_yyy = list_simulation[:,we_index,i]
			TMP_yer = simulation._xiWe[:,we_index,2]	
			TMP_cut = (TMP_yer>0.)
			TMP_xxx = TMP_xxx[ TMP_cut ]
			TMP_yyy = TMP_yyy[ TMP_cut ]
			TMP_yer = TMP_yer[ TMP_cut ]
			TMP_coef = numpy.power(TMP_xxx,x_power)
			ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, linewidth=1, color='gray',alpha=0.3)
		
		### Mean simulation and 1 sigma variance
		TMP_xxx = simulation._xiWe[:,we_index,0]
		TMP_yyy = simulation._xiWe[:,we_index,1]
		TMP_yer = simulation._xiWe[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_yyy_minus_one_sigma = TMP_yyy - TMP_yer*numpy.sqrt(nb_realisation)
		TMP_yyy_plus_one_sigma  = TMP_yyy + TMP_yer*numpy.sqrt(nb_realisation)
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, markersize=10,linewidth=2, color='blue')
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy_minus_one_sigma, linewidth=2, color='blue',linestyle='dashed')
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy_plus_one_sigma, linewidth=2, color='blue',linestyle='dashed')

		### Data
		TMP_xxx = data._xiWe[:,we_index,0]
		TMP_yyy = data._xiWe[:,we_index,1]
		TMP_yer = data._xiWe[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='red', ecolor='black',linestyle='None')

		###
		xxx = data._xiWe[:,we_index,0]
		xxx = xxx[ data._xiWe[:,we_index,2]>0. ]

		if (x_power==0):
			plt.ylabel(r'$'+data._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s| \cdot '+data._label+' (|s|) \, [$h$^{-1}$Mpc$]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2} \cdot '+data._label+' (|s|) \, [($h$^{-1}$Mpc$)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [$h$^{-1}$Mpc$]$', fontsize=40)
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlim([ data._min1D,data._max1D ])

		mng = plt.get_current_fig_manager()
		props = dict(boxstyle='round', facecolor=None, alpha=0.)
		ax.text(0.05, 0.95, r'$'+data._label_wedge[we_index]+'$', transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)

		myTools.deal_with_plot(False,False,True)
		name = path_to_where_to_plot + 'xi_' + data._correlation + '_We_index_' + str(we_index) + '_rescale_index_' + str(x_power) + '.png'
		plt.savefig(name, bbox_inches='tight')
		plt.clf()
		#plt.show()

		return

	
	plot_1d(0)
	plot_1d(1)
	plot_1d(2)
	
	plot_we(0,2)
	plot_we(1,2)
	plot_we(2,2)
	plot_we(3,2)


	return
def plot_dataLya_dataCIV():
	"""
		get similar plots as Delubac et al. for the wedges
	"""

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'nb_wedges'               : 3,
		'index_multipole_max' : 4,
		'remove_residuales' : 'lambda_OBS',
		'path_to_txt_file_folder': '',
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'NOTHING'
	}

	dic_class['name'] = 'Ly\\alpha \, forest'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
	data_lya = correlation_3D.Correlation3D(dic_class)
	data_lya.set_error_on_covar_matrix('subsampling')
	print data_lya._meanZ

	dic_class['name'] = 'CIV \, forest'
	dic_class['f1'] = 'CIV'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
	data_civ = correlation_3D.Correlation3D(dic_class)
	data_civ.set_error_on_covar_matrix('subsampling')
	print data_civ._meanZ

	
	def plot_1d(x_power=0, other=[]):

		list_corr_to_plot = other
		color = ['blue','red','green','orange']
		nb = 0

		for el in list_corr_to_plot:
			TMP_xxx = el._xi1D[:,0]
			TMP_yyy = el._xi1D[:,1]
			TMP_yer = el._xi1D[:,2]
			cut = (TMP_yer>0.)
		
			TMP_cut = (TMP_yer>0.)
			if (TMP_xxx[cut].size==0):
				continue
			TMP_xxx = TMP_xxx[ TMP_cut ]
			TMP_yyy = TMP_yyy[ TMP_cut ]
			TMP_yer = TMP_yer[ TMP_cut ]
			TMP_coef = numpy.power(TMP_xxx,x_power)
			plt.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, marker='o',alpha=0.6,label=r'$'+el._name+'$', color=color[nb])
			nb += 1

		###
		xxx = list_corr_to_plot[0]._xi1D[:,0]
		xxx = xxx[ list_corr_to_plot[0]._xi1D[:,2]>0. ]

		if (x_power==0):
			plt.ylabel(r'$'+list_corr_to_plot[0]._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s| \cdot '+list_corr_to_plot[0]._label+' (|s|) \, [$h$^{-1}$Mpc$]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2} \cdot '+list_corr_to_plot[0]._label+' (|s|) \, [($h$^{-1}$Mpc$)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [$h$^{-1}$Mpc$]$', fontsize=40)
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlim([ list_corr_to_plot[0]._min1D,list_corr_to_plot[0]._max1D ])
		myTools.deal_with_plot(False,False,True)
		plt.legend(fontsize=40, numpoints=1,ncol=1, loc=0)
		name = path_to_where_to_plot + 'xi_' + list_corr_to_plot[0]._correlation + '_1D_lya_and_civ_rescale_index_' + str(x_power) + '.png'
		#plt.savefig(name, bbox_inches='tight')
		#plt.clf()
		plt.show()

		return
	def plot_we(we_index,x_power):

		fig = plt.figure()
		ax = fig.add_subplot(111)

		### Data
		TMP_xxx = data_lya._xiWe[:,we_index,0]
		TMP_yyy = data_lya._xiWe[:,we_index,1]
		TMP_yer = data_lya._xiWe[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='blue', ecolor='black',linestyle='None')

		### Data
		TMP_xxx = data_civ._xiWe[:,we_index,0]
		TMP_yyy = data_civ._xiWe[:,we_index,1]
		TMP_yer = data_civ._xiWe[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='red', ecolor='black',linestyle='None')

		###
		xxx = data_lya._xi1D[:,0]
		xxx = xxx[ data_lya._xi1D[:,2]>0. ]

		if (x_power==0):
			plt.ylabel(r'$'+data_lya._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s| \cdot '+data_lya._label+' (|s|) \, [$h$^{-1}$Mpc$]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2} \cdot '+data_lya._label+' (|s|) \, [($h$^{-1}$Mpc$)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [$h$^{-1}$Mpc$]$', fontsize=40)
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlim([ data_lya._min1D,data_lya._max1D ])

		mng = plt.get_current_fig_manager()
		props = dict(boxstyle='round', facecolor=None, alpha=0.)
		ax.text(0.05, 0.95, r'$'+data_lya._label_wedge[we_index]+'$', transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)

		myTools.deal_with_plot(False,False,True)
		plt.show()

		return

	#plot_1d(0,[data_lya,data_civ])
	#plot_1d(1,[data_lya,data_civ])
	#plot_1d(2,[data_lya,data_civ])
	plot_we(0,2)
	plot_we(1,2)
	plot_we(2,2)

	return
def plot_simulation_raw_cooked():
	
	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'nb_wedges'               : 3,
		'index_multipole_max' : 4,
		'remove_residuales' : 'lambda_OBS',
		'path_to_txt_file_folder': '',
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'NOTHING'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '',
		'prefix2'   : '',
	}
	if (dic_class['correlation']=='q_f'):
		path_to_BAOFIT     = dic_simu['path_to_simu'] + 'Results'+dic_simu['prefix']+'/BaoFit_q_f__LYA__QSO'+dic_simu['prefix2']+'/bao2D.'
		path_to_mapping_we =  '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_we.npy'
		path_to_mapping_1D =  '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_1D.npy'
	elif (dic_class['correlation']=='f_f'):
		path_to_BAOFIT     = dic_simu['path_to_simu'] + 'Results'+dic_simu['prefix']+'/BaoFit_f_f__LYA'+dic_simu['prefix2']+'/bao2D.'
		path_to_mapping_we = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_we.npy'
		path_to_mapping_1D = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_1D.npy'

	### Raw
	dic_class['name']  = 'Raw'
	dic_simu['prefix'] = '_raw_from_JeanMarc'
	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
	if (dic_class['correlation']=='q_f'):
		path_to_BAOFIT     = dic_simu['path_to_simu'] + 'Results'+dic_simu['prefix']+'/BaoFit_q_f__LYA__QSO'+dic_simu['prefix2']+'/bao2D.'
	elif (dic_class['correlation']=='f_f'):
		path_to_BAOFIT     = dic_simu['path_to_simu'] + 'Results'+dic_simu['prefix']+'/BaoFit_f_f__LYA'+dic_simu['prefix2']+'/bao2D.'
	raw = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, None, path_to_BAOFIT)
	raw.set_values_on_mean_simulation(dic_simu)
	xi1D_data, xi1D_fit = raw.get_data_and_fit_1d(path_to_mapping_1D)
	raw._xi1D = xi1D_data
	raw_xi1D_fit = xi1D_fit
	xi1D_data, xi1D_fit = raw.get_data_and_fit_we(path_to_mapping_we)
	raw._xiWe = xi1D_data
	raw_xiWe_fit = xi1D_fit

	### Cooked
	dic_simu['prefix2'] = '__withMetalsTemplates'
	dic_simu['prefix'] = ''
	dic_class['name'] = 'Cooked'
	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
	if (dic_class['correlation']=='q_f'):
		path_to_BAOFIT     = dic_simu['path_to_simu'] + 'Results'+dic_simu['prefix']+'/BaoFit_q_f__LYA__QSO'+dic_simu['prefix2']+'/bao2D.'
	elif (dic_class['correlation']=='f_f'):
		path_to_BAOFIT     = dic_simu['path_to_simu'] + 'Results'+dic_simu['prefix']+'/BaoFit_f_f__LYA'+dic_simu['prefix2']+'/bao2D.'
	cooked = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, None, path_to_BAOFIT)
	cooked.set_values_on_mean_simulation(dic_simu)
	xi1D_data, xi1D_fit = cooked.get_data_and_fit_1d(path_to_mapping_1D)
	cooked._xi1D = xi1D_data
	cooked_xi1D_fit = xi1D_fit
	xi1D_data, xi1D_fit = cooked.get_data_and_fit_we(path_to_mapping_we)
	cooked._xiWe = xi1D_data
	cooked_xiWe_fit = xi1D_fit

	### Range
	cooked._min1D = 40.
	cooked._max1D = 180.

	def plot_1d(x_power):

		fig = plt.figure()
		ax = fig.add_subplot(111)

		### Raw
		TMP_xxx = raw._xi1D[:,0]
		TMP_yyy = raw._xi1D[:,1]
		TMP_yer = raw._xi1D[:,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='blue',linestyle='None')

		### Raw fit
		TMP_xxx = raw_xi1D_fit[:,0]
		TMP_yyy = raw_xi1D_fit[:,1]
		TMP_yer = raw_xi1D_fit[:,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, linewidth=2, color='blue',linestyle='dashed')

		### Cooked
		TMP_xxx = cooked._xi1D[:,0]
		TMP_yyy = cooked._xi1D[:,1]
		TMP_yer = cooked._xi1D[:,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='green',linestyle='None')

		### Cooked fit
		TMP_xxx = cooked_xi1D_fit[:,0]
		TMP_yyy = cooked_xi1D_fit[:,1]
		TMP_yer = cooked_xi1D_fit[:,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, linewidth=2, color='green',linestyle='dashed')

		###
		xxx = cooked._xi1D[:,0]
		xxx = xxx[ cooked._xi1D[:,2]>0. ]

		if (x_power==0):
			plt.ylabel(r'$'+cooked._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s| \cdot '+cooked._label+' (|s|) \, [$h$^{-1}$Mpc$]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2} \cdot '+cooked._label+' (|s|) \, [($h$^{-1}$Mpc$)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [$h$^{-1}$Mpc$]$', fontsize=40)
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlim([ cooked._min1D,cooked._max1D ])
		ax.set_xticks([ i for i in numpy.arange(cooked._min1D, cooked._max1D+20., 20.) ])
		myTools.deal_with_plot(False,False,True)
		

		name = path_to_where_to_plot + 'xi_' + cooked._correlation + '_1D_rescale_index_' + str(x_power) + '.png'
		#plt.savefig(name, bbox_inches='tight')
		#plt.clf()
		plt.show()
		
		return
	def plot_we(we_index,x_power):

		fig = plt.figure()
		ax = fig.add_subplot(111)

		### Raw
		TMP_xxx = raw._xiWe[:,we_index,0]
		TMP_yyy = raw._xiWe[:,we_index,1]
		TMP_yer = raw._xiWe[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='blue',linestyle='None')

		### Raw fit
		TMP_xxx = raw_xiWe_fit[:,we_index,0]
		TMP_yyy = raw_xiWe_fit[:,we_index,1]
		TMP_yer = raw_xiWe_fit[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, linewidth=2, color='blue',linestyle='dashed')
		
		### Cooked
		TMP_xxx = cooked._xiWe[:,we_index,0]
		TMP_yyy = cooked._xiWe[:,we_index,1]
		TMP_yer = cooked._xiWe[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, yerr=TMP_coef*TMP_yer, markersize=10,linewidth=2, fmt='o', color='green',linestyle='None')

		### Cooked fit
		TMP_xxx = cooked_xiWe_fit[:,we_index,0]
		TMP_yyy = cooked_xiWe_fit[:,we_index,1]
		TMP_yer = cooked_xiWe_fit[:,we_index,2]
		TMP_cut = (TMP_yer>0.)
		TMP_xxx = TMP_xxx[ TMP_cut ]
		TMP_yyy = TMP_yyy[ TMP_cut ]
		TMP_yer = TMP_yer[ TMP_cut ]
		TMP_coef = numpy.power(TMP_xxx,x_power)
		ax.errorbar(TMP_xxx, TMP_coef*TMP_yyy, linewidth=2, color='green',linestyle='dashed')

		###
		xxx = cooked._xi1D[:,0]
		xxx = xxx[ cooked._xi1D[:,2]>0. ]

		if (x_power==0):
			plt.ylabel(r'$'+cooked._label+' (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s| \cdot '+cooked._label+' (|s|) \, [$h$^{-1}$Mpc$]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2} \cdot '+cooked._label+' (|s|) \, [($h$^{-1}$Mpc$)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [$h$^{-1}$Mpc$]$', fontsize=40)
		plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
		plt.xlim([ cooked._min1D,cooked._max1D ])
		ax.set_xticks([ i for i in numpy.arange(cooked._min1D, cooked._max1D+20., 20.) ])

		mng = plt.get_current_fig_manager()
		props = dict(boxstyle='round', facecolor=None, alpha=0.)
		ax.text(0.05, 0.95, r'$'+cooked._label_wedge[we_index]+'$', transform=ax.transAxes, fontsize=30, verticalalignment='top', bbox=props)


		myTools.deal_with_plot(False,False,True)
		name = path_to_where_to_plot + 'xi_' + cooked._correlation + '_We_index_' + str(we_index) + '_rescale_index_' + str(x_power) + '.png'
		#plt.savefig(name, bbox_inches='tight')
		#plt.clf()
		plt.show()

		return

	
	plot_1d(0)
	plot_1d(1)
	plot_1d(2)
	
	plot_we(0,2)
	plot_we(1,2)
	plot_we(2,2)


	return

plot_data_simulation()
#plot_dataLya_dataCIV()
#plot_simulation_raw_cooked()







