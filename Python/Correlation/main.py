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


#all_t = [ '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits',
#'/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits',
#'/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits']
#myTools.append_files(all_t, '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/FitsFile_DR12_reOBS_eBOSS_noCoADD_Guy/DR12_primery/DR12_primery_reOBS_eBOSS_noCoADD.fits')

i = sys.argv[1]
j = sys.argv[2]

"""
for i in range(4490534, 4490634):
	command = 'qdel ' + str(i)
	print command
	subprocess.call(command, shell=True)
"""

def plotOne():


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
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'no'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '_only_LR_noRand_noRSD', #
		'prefix2'   : '', #'_sigmaNL_0_pixelSize_3_15', #'_fit_from_10_to_180_sigmaNL_0_pixelSize_3_15', #'_fit_from_10_to_180'
	}
	dic_send_BAOFIT = {
		'realisation_type'        : 'Simu_stack', #'subsampling', #
		#'correlation_matrix_path' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_cor_2D_from_fit.npy',
		#'correlation_matrix_path' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results/xi_A_delta_delta_LYA_result_cor_2D_from_fit.npy',
		'correlation_matrix_path' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_cor_2D_from_fit.npy',
		#'correlation_matrix_path' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results/xi_delta_QSO_LYA_QSO_result_cor_2D_from_fit.npy',
		#'correlation_matrix_path' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_no_metals/xi_delta_QSO_LYA_QSO_result_cor_2D_from_fit.npy',
		#'correlation_matrix_path' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/xi_delta_QSO_LYA_QSO_subsampling_cor_2D_from_fit.npy',
		'saving_data'             : True,
		'scan'                    : False,
		'toyMC'                   : 0,
		'dist_matrix'             : True,
		'only_diagonal'           : False,
		'with_metals_templates'   : False,
		'get_param_from_file'     : False,
		'dic_simu'                : dic_simu, #None, #dic_simu
		'prefix'                  : ''
	}
	#dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Results'+dic_simu['prefix']+'/'
	#print dic_class['path_to_txt_file_folder']
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_Margala__projected_PlanckCosmo_correctedDistance/'
	print dic_class['path_to_txt_file_folder']
	corr = correlation_3D.Correlation3D(dic_class)
	#corr.write_metal_model(False)
	#corr.send_BAOFIT(dic_send_BAOFIT)
	#corr.send_pyLyA(with_metals_templates=False,prefix2='_rmin10')
	#corr.save_list_realisation_simulation(dic_class, dic_simu,False)
	#corr.set_values_on_mean_simulation_or_realisation_type(dic_simu)
	#corr.send_BAOFIT(dic_send_BAOFIT)
	#corr.send_pyLyA(with_metals_templates=False,prefix2='_rmin10')
	#corr.set_values_on_mean_simulation_or_realisation_type(None,'shuffleForest')
        #corr.set_values_on_mean_simulation_or_realisation_type(None,'shuffleQSO')
        #corr.set_values_on_mean_simulation_or_realisation_type(None,'randomQSO')
	#corr.set_values_on_mean_simulation_or_realisation_type(None,'randomQSOXYZ')
	#corr.save_mean_realisation_simulation_correlation_matrix(dic_simu)
	#corr.write_metal_model(False)
	#corr.write_metal_model_mean_simulation(dic_class, dic_simu)
	#corr.read_data_from_BAOFIT_data_file('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/TESTS/Nicolas_Correlation/xcf_dr12_v5_8_guy_c2_baseline_projected')

	### Mapping from 2D to wedge
	'''
	#path_to_save = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_we_with_4_wedges'
	path_to_save = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_A_delta_delta_LYA_result_mapping_2D_to_we_with_4_wedges'
	print path_to_save
	mapping_2D_to_we = corr.get_mapping_2D_to_we()
	numpy.save( path_to_save, mapping_2D_to_we )
	mapping_2D_to_we[ mapping_2D_to_we == 0. ] = numpy.float('nan')
	myTools.plot2D( mapping_2D_to_we[:,:,30,0] )
	myTools.plot2D( mapping_2D_to_we[:,:,30,1] )
	myTools.plot2D( mapping_2D_to_we[:,:,30,2] )
	myTools.plot2D( mapping_2D_to_we[:,:,30,3] )
	'''
	### Mapping from 2D to 1D
	'''
	#path_to_save = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_1D'
	path_to_save = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_A_delta_delta_LYA_result_mapping_2D_to_1D'
	mapping_2D_to_1D = corr.get_mapping_2D_to_1D()
	numpy.save( path_to_save, mapping_2D_to_1D )
	#mapping_2D_to_1D[ mapping_2D_to_1D == 0. ] = numpy.float('nan')
	#myTools.plot2D( mapping_2D_to_1D[:,:,20] )
	'''

	### Look at the covariance matrix
	'''
	path_to_load = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_A_delta_delta_LYA_result_cov_2D.npy'
	path_to_save = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_A_delta_delta_LYA_result_cor_2D_from_fit'
	print path_to_load
	print path_to_save
	myTools.save_fited_cor_from_cov(path_to_load, path_to_save,50,50,'f_f')
	'''	
	'''
	#path_to_load = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'
        #path_to_save = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_delta_QSO_LYA_QSO_result_cor_2D_from_fit'
	path_to_load = dic_class['path_to_txt_file_folder']+'/xi_delta_QSO_LYA_QSO_subsampling_cov_1D.npy'
        path_to_save = dic_class['path_to_txt_file_folder']+'/xi_delta_QSO_LYA_QSO_subsampling_cor_1D_from_fit'
        print path_to_load
        print path_to_save
        myTools.save_fited_cor_from_cov(path_to_load, path_to_save,50,100,'q_f')
	path_to_load = dic_class['path_to_txt_file_folder']+'/xi_delta_QSO_LYA_QSO_subsampling_cov_2D.npy'
        path_to_save = dic_class['path_to_txt_file_folder']+'/xi_delta_QSO_LYA_QSO_subsampling_cor_2D_from_fit'
        print path_to_load
        print path_to_save
        myTools.save_fited_cor_from_cov(path_to_load, path_to_save,50,100,'q_f')
	'''
	
	### Wick decomposition
	#corr.plot_Wick_calculation('1')
	#corr.plot_Wick_calculation('12')
	#corr.plot_Wick_calculation('2D',['1_withDelta','12_withDelta'],['T1','T12'])
	#corr.plot_Wick_calculation(['1','12'],['T1','T12'])
	#corr.plot_Wick_calculation(['1_withDelta'],['T1'])

	print corr._meanZ
	dic_CAMB_corr = {
	'z' : corr._meanZ,
	'h_0' : 0.70,
	'omega_matter_0' : 0.27,
	'omega_lambda_0' : 0.73,
	'source'         : 'CAMB_Mocks_me'
	}
	CAMB1 = CAMB.CAMB(dic_CAMB_corr)
	dic_CAMB = {
		'mulpol_index' : 0,
		'start_fit'   : 40.,
		'end_fit'     : 180.,
		'b' : -1.,
		'roof' : 0.,
		'fix_roof_nul' : True,
		'guess_b' : False,
		'min_for_guess' : 20.,
		'max_for_guess' : 50.,
		'CAMB' : CAMB1
	}

	### Look at the distortion matrix
	#corr.plot_distortion_matrix('1D')
	#corr.plot_distortion_matrix('2D')
	#xi = corr.get_CAMB(dic_CAMB_corr,'2D',True)
	#corr._xi2D = xi
	'''
	xi = corr.get_CAMB(dic_CAMB_corr,'Mu',False)
	corr._xiMu = xi
	xi = corr.get_CAMB(dic_CAMB_corr,'We',False)
	corr._xiWe = xi
	xi = corr.get_CAMB(dic_CAMB_corr,'1D',False)
	corr._xi1D = xi
	'''

	'''
	dic_simu = {
                'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/',
                'nb_box' : 10,
                'nb_simu' : 10,
                'prefix' : '_test_mean_stack_lambda_OBS'
        }
        dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Results_test_mean_stack_lambda_OBS/'
        dic_class['name'] = 'q-f'
        corr2 = correlation_3D.Correlation3D(dic_class)

	corr._xi1D[:,1]   -= corr2._xi1D[:,1]
        corr._xi2D[:,:,1] -= corr2._xi2D[:,:,1]
	'''

	"""
	xiMu, xiWe, xi1D, xi2D = corr.read_metal_model( 'LYA' )	
	path = corr._path_to_txt_file_folder + corr._prefix + '_distortionMatrix_2D_'+ corr._middlefix + '.txt'
	print path
	matrix = numpy.loadtxt(path)

	corr._xi2D[:,:,0] = xi2D[:,:,0,2]
	corr._xi2D[:,:,1] = xi2D[:,:,1,2]

	corr.plot_2d(0)
        corr.plot_2d(1)
        corr.plot_2d(2)
	"""

	#matrix[matrix==0.] = numpy.float('nan')
	#for l in numpy.arange(corr._nbBin2D):
	#	matrix[l,l] = numpy.float('nan')
	#myTools.plot2D(matrix)

	"""
	yyy_Camb = numpy.dot(matrix,xi2D[:,:,1,2].flatten())
	corr._xi2D[:,:,1] = myTools.convert1DTo2D(yyy_Camb,corr._nbBinX2D, corr._nbBinY2D)
	corr.plot_2d(0)
        corr.plot_2d(1)
        corr.plot_2d(2)
	"""

	#corr.save_list_realisation('subsampling', 80)
	#corr.save_list_realisation('shuffleQSO', 100)
	#corr.save_list_realisation('shuffleForest', 100)
	#corr.save_list_realisation('randomQSO', 1000)
	#corr.save_list_realisation('randomQSOXYZ', 1000)
	#corr.save_list_realisation('randomQSOFromFile', 10)
        #corr.set_values_on_mean_simulation_or_realisation_type(None,'shuffleForest')
        #corr.set_values_on_mean_simulation_or_realisation_type(None,'shuffleQSO')
	#corr.set_values_on_mean_simulation_or_realisation_type(None,'randomQSO')
        #corr.set_values_on_mean_simulation_or_realisation_type(None,'randomQSOXYZ')
	#corr.set_values_on_mean_simulation_or_realisation_type(None,'randomQSOFromFile')

	#corr.plot_cov_cor_matrix('subsampling','1D')
	#corr.plot_cov_cor_matrix('subsampling','2D')
	### Fit
	#corr.set_error_on_covar_matrix('subsampling')
	#corr.save_value_xi_in_txt()
	#corr.send_BAOFIT(dic_send_BAOFIT)

	"""

		

	
	c0, c2, c4 = CAMB1.get_multipol_coef(-0.168777,1.364)
	print c0, c2, c4
	
	'''
	dic_CAMB['mulpol_index'] = 0
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,0,:],dic_CAMB,False)
	print dic_CAMB['b']
	#dic_CAMB['b'] = c0
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB,0,False)
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB,1,False)
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB,2,False)

	dic_CAMB['mulpol_index'] = 2
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,2,:],dic_CAMB,False)
	print dic_CAMB['b']
	#dic_CAMB['b'] = c2
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB,0,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB,1,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB,2,False)

	dic_CAMB['mulpol_index'] = 4
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,4,:],dic_CAMB,False)
	print dic_CAMB['b']
	#dic_CAMB['b'] = c4
	corr.plot_CAMB(corr._xiMul[:,4,:],dic_CAMB,0,False)
	corr.plot_CAMB(corr._xiMul[:,4,:],dic_CAMB,1,False)
	corr.plot_CAMB(corr._xiMul[:,4,:],dic_CAMB,2,False)
	'''

	#
	#for ii in numpy.arange(0,50):
	#	corr.plot_slice_2d(ii)
	#for ii in numpy.arange(0,100):
	#	corr.plot_slice_2d(None,ii)
	

	
	for ii in numpy.arange(0,50):
		corr.plot_slice_2d(ii)
	"""	
	print corr._meanZ
	corr.plot_grid(0,True)
	corr.plot_grid(1,True)
	corr.plot_grid(2,True)
	
	#plt.hist(corr._xi2D[:,:,1].flatten(),bins=30)
	#plt.show()
	#plt.hist(corr._xi2D[:,:,1].flatten()/corr._xi2D[:,:,2].flatten(),bins=30)
        #plt.show()
	
	corr.plot_1d(0)
	corr.plot_1d(1)
	corr.plot_1d(2)
	#for ii in range(0,4):
	#	for jj in range(0,3):
	#		corr.plot_we(ii,jj)
	for ii in range(0,3):
		corr.plot_we_all(ii)
	
	corr._xi2D[:,:,1] /= corr._xi2D[:,:,2]	
	corr._xiMu[:,:,2] /= corr._xiMu[:,:,3]
	corr.plot_slice_2d(0)
	
	corr.plot_2d(0)
	corr.plot_2d(1)
	corr.plot_2d(2)
	corr.plot_mu(0)
	corr.plot_mu(1)
	corr.plot_mu(2)
	
	corr.plot_multipol(0)
	corr.plot_multipol(1)
	corr.plot_multipol(2)
	dic_CAMB = corr.fit_CAMB(dic_CAMB)
	corr.plot_CAMB(dic_CAMB,0)
	corr.plot_CAMB(dic_CAMB,1)
	corr.plot_CAMB(dic_CAMB,2)
	corr.plot_map_sub_sampling()
	

	#list_simulation = [ ['Simulation', '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results/xi_delta_QSO_LYA_QSO_result'] ]
	#corr.plot_cov_cor_matrix_different_method(['subsampling'],'1D',None,list_simulation )
	#corr.plot_cov_cor_matrix_different_method(['subsampling'],'2D',[ ['Simulation','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_withRand_XYZ_nicolasEstimator/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'] ] )
	
	#corr.plot_cov_cor_matrix_different_method(['subsampling'],'2D',[ ['simu','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'], [ 'not \, removing', '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_LYA_QSO_subsampling_cov_2D.npy'] ] )
	
	#path = dic_class['path_to_txt_file_folder']
	#corr.plot_cov_cor_matrix_different_method(['shuffleForest','shuffleQSO','randomQSO','subsampling'],'We')
	
	return
def plotMany():

	list_correlation = ['LYA','SIIV','CIV']
	list_corr = []
	list_value = []

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'correlation': 'q_f',
		'path_to_txt_file_folder': 'NOTHING',
		'f1': 'LYA',
		'f2': 'LYA',
		'q1': 'QSO',
		'q2': 'QSO',
		'name' : 'Data'
	}
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
	dic_class['name'] = 'Data'

	for el in list_correlation:
		dic_class['f1']   = el
		dic_class['name'] = el+' \, \\times \, QSO'
		corr = correlation_3D.Correlation3D(dic_class)
		corr.set_error_on_covar_matrix('subsampling')
		list_corr += [corr]
		list_value += [ corr._xi1D[2,1] ]
		print el, corr._meanZ

	'''
	print list_correlation[0], list_value[0]
	for i in numpy.arange( 1, len(list_correlation) ):
		list_corr[i].multiply_by_constant( list_value[0]/list_value[i] )
		print list_correlation[i],'/',list_correlation[0], '   ',  list_value[i]/list_value[0]
	'''

	list_corr[0].plot_1d(0, list_corr[1:], False )
	list_corr[0].plot_1d(1, list_corr[1:], False )
	list_corr[0].plot_1d(2, list_corr[1:], False )

	return
def plotBosseBOSS():

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
                'nb_wedges'               : 4,
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '',
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'Mean \, Simulation'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '_withRand_XYZ_nicolasEstimator'
	}
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Results_withRand_XYZ_nicolasEstimator/'
	#corr4 = correlation_3D.Correlation3D(dic_class)
	#corr4.save_list_realisation_simulation(dic_class, dic_simu)
        #corr4.set_values_on_mean_simulation(dic_simu)

	dic_class['q1'] = 'QSO'
        dic_class['remove_residuales'] = 'lambda_OBS'
        #dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_coAdd_2016_05_26/'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26/' 
        dic_class['name'] = 'LYA - LYA'
        corr = correlation_3D.Correlation3D(dic_class)
	#corr.save_list_realisation('subsampling', 80)
	corr.set_error_on_covar_matrix('subsampling')
        #print corr._meanZ

	#corr4.print_null_test()
	#print ''
	#corr4.print_null_test('subsampling')
	#print ''
        #corr4.print_null_test('subsampling','_from_fit')

	###
	#dic_class['q1'] = 'QSO'
	#dic_class['remove_residuales'] = 'lambda_OBS'
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_coAdd_2016_05_26/'
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26/'
	#dic_class['name'] = 'residuals \, LYA - QSO \, (shuffle \, Forest)'
	#dic_class['name'] = 'Sky - Sky'
	#corr5 = correlation_3D.Correlation3D(dic_class)
	#corr5.set_values_on_mean_simulation_or_realisation_type(None,'randomQSOXYZ')
	#corr5.save_list_realisation('shuffleForest', 1000)
	#corr5.save_list_realisation('subsampling', 80)
	#corr5.set_error_on_covar_matrix('subsampling')
	#print corr5._meanZ

        #corr5.print_null_test()
        #print ''
        #corr5.print_null_test('randomQSOXYZ')
	#print ''
        #corr5.print_null_test('randomQSOXYZ','_from_fit')
	#return

	#corr5.plot_cov_cor_matrix('subsampling','1D')
	#corr5._xi1D[:,1] /= corr5._xi1D[:,2]
	#corr5._xi2D[:,:,1] /= corr5._xi2D[:,:,2]
	#corr5.plot_1d(0)
	#corr5.plot_2d(0)
	#corr5.plot_2d(1)
	#corr4.plot_2d(2)


	###
	corr_1 = copy.deepcopy(corr)
	corr_1._name = 'randomQSO \, XYZ'
	corr_1.set_values_on_mean_simulation_or_realisation_type(None,'randomQSOXYZ')
	###
	corr_2 = copy.deepcopy(corr)
	corr_2._name = 'randomQSO'
        corr_2.set_values_on_mean_simulation_or_realisation_type(None,'randomQSO')
	###
	corr_3 = copy.deepcopy(corr)
	corr_3._name = 'shuffleForest'
        corr_3.set_values_on_mean_simulation_or_realisation_type(None,'shuffleForest')
	###
	corr_4 = copy.deepcopy(corr)
	corr_4._name = 'shuffleQSO'
        corr_4.set_values_on_mean_simulation_or_realisation_type(None,'shuffleQSO')
	
	corr_1.print_null_test()
        print ''
        corr_1.print_null_test('randomQSOXYZ')
        print ''
        corr_1.print_null_test('randomQSOXYZ','_from_fit')

	color = ['blue','red','green','orange','black','blue']
	plt.hist(corr_1._xi2D[:,:,1].flatten()/corr_1._xi2D[:,:,2].flatten(),histtype='step',bins=30,label=r'$'+corr_1._name+'$', color=color[0])
        plt.hist(corr_2._xi2D[:,:,1].flatten()/corr_2._xi2D[:,:,2].flatten(),histtype='step',bins=30,label=r'$'+corr_2._name+'$', color=color[1])
        plt.hist(corr_3._xi2D[:,:,1].flatten()/corr_3._xi2D[:,:,2].flatten(),histtype='step',bins=30,label=r'$'+corr_3._name+'$', color=color[2])
        plt.hist(corr_4._xi2D[:,:,1].flatten()/corr_4._xi2D[:,:,2].flatten(),histtype='step',bins=30,label=r'$'+corr_4._name+'$', color=color[3])
	#plt.hist(corr._xi2D[:,:,1].flatten()/corr._xi2D[:,:,2].flatten(),histtype='step',bins=30,label=r'$'+corr._name+'$',color=color[4])
	plt.show()

	myTools.fitAGaussian( (corr._xi2D[:,:,1]/corr._xi2D[:,:,2]).flatten(), 50,p0=[1,0.,1.],log=True)
	myTools.fitAGaussian( (corr_1._xi2D[:,:,1]/corr_1._xi2D[:,:,2]).flatten(), 50, p0=[1,0.,1.],log=True)
	myTools.fitAGaussian( (corr_2._xi2D[:,:,1]/corr_2._xi2D[:,:,2]).flatten(), 50,p0=[1,0.,1.],log=True)
	myTools.fitAGaussian( (corr_3._xi2D[:,:,1]/corr_3._xi2D[:,:,2]).flatten(), 50, p0=[1,0.,1.],log=True)
	myTools.fitAGaussian( (corr_4._xi2D[:,:,1]/corr_4._xi2D[:,:,2]).flatten(), 50, p0=[1,0.,1.],log=True)




	corr_1.plot_1d(0,[corr_2,corr_3,corr_4],False)
	corr_1.plot_1d(1,[corr_2,corr_3,corr_4],False)
	corr_1.plot_1d(2,[corr_2,corr_3,corr_4],False)
	corr_1.plot_we(0,[corr_2,corr_3,corr_4],False)
	corr_1.plot_we(1,[corr_2,corr_3,corr_4],False)
	corr_1.plot_we(2,[corr_2,corr_3,corr_4],False)


	for i in numpy.arange(0,50):
		corr.plot_slice_2d(i,None,[corr6])


	return
### BAOFIT + Many simulations
def BAOFIT_Many_simulations():

	isPyLyA = True

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'nb_wedges'          : 4,
		'correlation': 'q_f',
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'q-f \, raw'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '', #'_no_metals', #'_raw_from_JeanMarc', #
		'prefix2'   : '_withMetalsTemplates_rmin10', #'_withMetalsTemplates', #'_withMetalsTemplates_rmin10', #
	}

	### First is raw
	aMB = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic=dic_class, index_parameter=None, dic_Q=None, dic_simu=dic_simu, load_correlation=False, isPyLyA=isPyLyA)

	
	### Second is other
	other = []
	'''
	dic_class['name']  = 'q-f \, full'
	dic_simu['prefix'] = ''
	dic_simu['prefix2'] = '__withMetalsTemplates'
	other += [ annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,None,dic_simu,False) ]
	#aMB = other[0]
	#other = []
	'''
	
	### If f_f
	if (dic_class['correlation']=='f_f'):

		if (isPyLyA==True):
			#aMB.print_results_for_chi2()
			
			aMB.plot_scatter_hist(0,1)
			aMB.plot_scatter_hist(0,2)
			aMB.plot_scatter_hist(0,3)
			aMB.plot_scatter_hist(1,2)
			aMB.plot_scatter_hist(1,3)
			aMB.plot_scatter_hist(2,3)
			
			aMB.print_results(0)
			aMB.print_results(1)
			aMB.print_results(2)
			aMB.print_results(3)
			aMB.print_results(5)

		aMB.print_list_parameter(8)
		print
		aMB.print_list_parameter(9)

		aMB.print_results_for_chi2()

		aMB.plot_scatter_hist(0,1)
		aMB.plot_scatter_hist(0,8)
		aMB.plot_scatter_hist(0,9)
		aMB.plot_scatter_hist(1,8)
		aMB.plot_scatter_hist(1,9)
		aMB.plot_scatter_hist(8,9)
		
		aMB.print_results(0)
		aMB.print_results(1)
		aMB.print_results(8)
		aMB.print_results(9)
	### If q_f
	elif (dic_class['correlation']=='q_f'):

		aMB.print_poll_plots(0)
		aMB.print_poll_plots(1)

		#aMB.print_list_parameter(11)
		print
		#aMB.print_list_parameter(12)

		#aMB.print_results_for_chi2(other)

		if (isPyLyA==True):
                        #aMB.print_results_for_chi2()

			'''
			aMB.plot_scatter_hist(1,5)
                        aMB.plot_scatter_hist(2,3)
                        aMB.plot_scatter_hist(2,7)
                        aMB.plot_scatter_hist(3,7)
			#aMB.print_results(9)
                        #aMB.print_results(11)
                        #aMB.print_results(13)
                        #aMB.print_results(15)
			'''
			for i in range(15):
				for j in range(i+1,15):
					aMB.plot_scatter_hist(i,j)
                        aMB.print_results(1)
                        aMB.print_results(2)
                        aMB.print_results(3)
                        aMB.print_results(5)
                        aMB.print_results(7)
		
		aMB.plot_scatter_hist(1,0,other)
		aMB.plot_scatter_hist(1,4,other)
		aMB.plot_scatter_hist(1,5,other)
		aMB.plot_scatter_hist(1,11,other)
		aMB.plot_scatter_hist(1,12,other)
		
		
		aMB.plot_scatter_hist(0,4,other)
		aMB.plot_scatter_hist(0,5,other)
		aMB.plot_scatter_hist(0,11,other)
		aMB.plot_scatter_hist(0,12,other)
		aMB.plot_scatter_hist(4,5,other)
		aMB.plot_scatter_hist(4,11,other)
		aMB.plot_scatter_hist(4,12,other)
		aMB.plot_scatter_hist(5,11,other)
		aMB.plot_scatter_hist(5,12,other)
		aMB.plot_scatter_hist(11,12,other)
		
		aMB.print_results(0,other)
		aMB.print_results(1,other)
		aMB.print_results(4,other)
		aMB.print_results(5,other)
		
		aMB.print_results(11,other)
		aMB.print_results(12,other)
	
	#aMB.set_mean_data()
	#print aMB.get_best_fit()
	#corr.plot_chi2_scan(100,100, [0.5,1.5,0.5,1.5], False, False, aMB.get_best_fit())

	return
def BAOFIT_One_fit_data():

	prefix2 = ''
	prefix2 = '_withMetalsTemplates_rmin10'
	#prefix2 += '_fiducial'
	#prefix2 += '_fit_down_to_10'
	#prefix2 += '_fit_up_to_300'
	
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
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'NOTHING'
	}
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/'

	if (dic_class['correlation']=='q_f'):
		path_to_BAOFIT     = dic_class['path_to_txt_file_folder'] + 'BaoFit_q_f__LYA__QSO'+prefix2+'/bao2D.'
		path_to_BAOFIT     = dic_class['path_to_txt_file_folder'] + '/pyLyA'+prefix2+'/cross_alone.'
		path_to_mapping_we =  '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_we_with_4_wedges.npy'
		path_to_mapping_1D =  '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_1D.npy'
	elif (dic_class['correlation']=='f_f'):
		path_to_BAOFIT     = dic_class['path_to_txt_file_folder'] + 'BaoFit_f_f__LYA'+prefix2+'/bao2D.'

		path_to_mapping_we = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_we_with_4_wedges.npy'
		path_to_mapping_1D = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_1D.npy'

	corr = annalyse_BAOFIT.AnnalyseBAOFIT(dic=dic_class,dic_simu=None, index_parameter=None, path_to_BAOFIT=path_to_BAOFIT, isPyLyA=True)

	'''
	for i in numpy.arange(0,50):
		corr.plot_slice_fit_2d(i)
	for i in numpy.arange(0,100):
		corr.plot_slice_fit_2d(None,i)
	corr.print_results()
	'''
	
	#corr.print_results()
	corr.plot_data_and_fit_1d(0,path_to_mapping_1D)
	corr.plot_data_and_fit_1d(1,path_to_mapping_1D)
	corr.plot_data_and_fit_1d(2,path_to_mapping_1D)
	#corr.plot_data_and_fit_we(0,path_to_mapping_we)
	#corr.plot_data_and_fit_we(1,path_to_mapping_we)
	#corr.plot_data_and_fit_we(2,path_to_mapping_we)
	
	corr.plot_fit_2d(0)
	corr.plot_fit_2d(1)
	corr.plot_fit_2d(2)
	corr.plot_residuals_2d(0)
	corr.plot_residuals_2d(1)
	corr.plot_residuals_2d(2)
	corr.plot_histo_residuals(nbBins=100,rmin=10, rmax=180)

	return
def BAOFIT_One_fit_simulation():

	on_stack = True

	dic_class = {
			'minXi': 0.,
			'maxXi': 200.,
			'nbBin': 50,
			'nbBinM': 25,
			'nb_Sub_Sampling': 80,
			'size_bin_calcul_s': 1.,
			'size_bin_calcul_m': 0.02,
			'nb_wedges'               : 4,
			'remove_residuales' : 'lambda_OBS',
			'index_multipole_max' : 4,
			'path_to_txt_file_folder': '',
			'correlation': 'f_f',
			'f1': 'LYA',
			'f2': '',
			'q1': 'QSO',
			'q2': '',
			'name' : 'simu'
	}
	dic_simu = {
			'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
			'nb_box' : 10,
			'nb_simu' : 10,
			'prefix' : '_raw_from_JeanMarc',
			'prefix2'   : '' #'__noDistortionMatrix_sigmaNL_0_pixelSize_3_15', #'__noDistortionMatrix_fit_from_10_to_180'
	}
	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
	path_to_BAOFIT = None

	if (dic_class['correlation']=='q_f'):
		if (on_stack): path_to_BAOFIT = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix']+'/BaoFit_q_f__LYA__QSO'+dic_simu['prefix2']+'/bao2D.'
		else: path_to_BAOFIT = dic_class['path_to_txt_file_folder'] +'/BaoFit_q_f__LYA__QSO'+dic_simu['prefix2']+'/bao2D.'
		path_to_mapping_we =  '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_we_with_4_wedges.npy'
		path_to_mapping_1D =  '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_1D.npy'
	elif (dic_class['correlation']=='f_f'):
		if (on_stack): path_to_BAOFIT = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix']+'/BaoFit_f_f__LYA'+dic_simu['prefix2']+'/bao2D.'
		else: path_to_BAOFIT = dic_class['path_to_txt_file_folder'] +'/BaoFit_f_f__LYA'+dic_simu['prefix2']+'/bao2D.'
		path_to_mapping_we = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_we_with_4_wedges.npy'
		path_to_mapping_1D = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_1D.npy'
	
	path_to_BAOFIT = dic_class['path_to_txt_file_folder'] +'pyLyA/auto_alone.'
	corr = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, None, path_to_BAOFIT,isPyLyA=True)
	#corr.set_values_on_mean_simulation(dic_simu)
	
	#corr.print_results()
	corr.plot_data_and_fit_1d(0,path_to_mapping_1D)
	corr.plot_data_and_fit_1d(1,path_to_mapping_1D)
	corr.plot_data_and_fit_1d(2,path_to_mapping_1D)
	corr.plot_data_and_fit_we(0,2,path_to_mapping_we)
	corr.plot_data_and_fit_we(1,2,path_to_mapping_we)
	corr.plot_data_and_fit_we(2,2,path_to_mapping_we)
	
	corr.plot_fit_2d(0)
	corr.plot_fit_2d(1)
	corr.plot_fit_2d(2)
	corr.plot_residuals_2d(0)
	corr.plot_residuals_2d(1)
	corr.plot_residuals_2d(2)
	corr.plot_histo_residuals()

	return
def BAOFIT_Compare_simulations():

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'nb_wedges'          : 4,
		'correlation': 'q_f',
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'q-f \, raw'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '_raw_from_JeanMarc',
		'prefix2'   : '' #'__withMetalsTemplates'
	}

	### qso_qso
	corr_q_f = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,None,dic_simu,False)

	### forest_forest
	dic_simu['prefix2'] = ''
	dic_class['correlation'] = 'f_f'
	dic_class['name']        = 'f-f \, raw'
	corr_f_f = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,None,dic_simu,False)

        corr_q_f.plot_scatter_hist_two_correlation(corr_q_f,11,12)
        corr_f_f.plot_scatter_hist_two_correlation(corr_f_f,8,9)

	corr_q_f.plot_scatter_hist_two_correlation(corr_f_f,11,8)
        corr_q_f.plot_scatter_hist_two_correlation(corr_f_f,11,9)
        corr_q_f.plot_scatter_hist_two_correlation(corr_f_f,12,8)
        corr_q_f.plot_scatter_hist_two_correlation(corr_f_f,12,9)

	return
def BAOFIT_Compare_simulations_pyLya():

	isPyLyA = True
	
	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'nb_wedges'          : 4,
		'correlation': 'q_f',
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'q-f \, raw'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '_raw_from_JeanMarc',
		'prefix2'   : '' #'_withMetalsTemplates_rmin10', #'_withMetalsTemplates', #'_withMetalsTemplates_rmin10', #
	}

	### First is raw
	aMB = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,None,dic_simu,False,isPyLyA=isPyLyA)

	### Second is other
	other = []
	dic_class['name']  = 'q-f \, cooked \, no \, metals'
	dic_simu['prefix'] = '_no_metals'
	dic_simu['prefix2'] = ''
	other += [ annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,None,dic_simu,False,isPyLyA=isPyLyA) ]

	aMB.plot_scatter_hist(1,5,other=other)
	aMB.plot_scatter_hist(2,3,other=other)
	

	return
def plot_compare_some_correlations():

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
                'nb_wedges'               : 3,
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '',
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'LYA - QSO'
	}

	### nb pairs LYA = 
	### nb pairs CIV = 
        dic_class['name'] = 'CAMB'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/' 
        corr_1 = correlation_3D.Correlation3D(dic_class)
	#corr_1.save_list_realisation('subsampling', 80)
	#corr_1.set_error_on_covar_matrix('subsampling')

	### Correlation 2
	### nb pairs LYA = 2.660747916700000e+10
	### nb pairs CIV = 
	dic_class['q1'] = 'QSO'
        dic_class['name'] = 'CAMB+topHat'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo_with_PDF_nicolas/'
	corr_2 = correlation_3D.Correlation3D(dic_class)
	#corr_2.save_list_realisation('subsampling', 80)
	#corr_2.set_error_on_covar_matrix('subsampling')
	
	xiMu, xiWe, xi1D, xi2D = corr_1.read_metal_model('SiIII(1207)')
	corr_1._xi2D = xi2D[:,:,:,2]

	xiMu, xiWe, xi1D, xi2D = corr_2.read_metal_model('SiIII(1207)')
	corr_2._xi2D = xi2D[:,:,:,2]
	
	corr_1.plot_slice_2d(0,None,[corr_2])
	corr_1.plot_slice_2d(1,None,[corr_2])

	"""
	print numpy.mean( corr_2._xi2D[:,:,2].flatten()/corr_1._xi2D[:,:,2].flatten() )
	plt.plot( corr_2._xi2D[:,:,2].flatten()/corr_1._xi2D[:,:,2].flatten(),alpha=0.7,linewidth=2)
	plt.legend(fontsize=30, numpoints=1,ncol=2, loc=2)
	plt.xlabel(r'$bin \, index$', fontsize=40)
	plt.ylabel(r'$\sigma_{DR14} / \sigma_{DR12}$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	corr_2.plot_we(0,2,[corr_1],True)
	corr_2.plot_1d(0,[corr_1],False)
	corr_2.plot_1d(1,[corr_1],False)
	corr_2.plot_1d(2,[corr_1],False)
	### Wedge
	for i in range(0,4):
		for j in range(0,3):
			corr_2.plot_we(i,j,[corr_1],True)
	"""

	return
def plot_compare_errors_correlations():

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
                'nb_wedges'               : 3,
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '',
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'LYA - QSO'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '', #
		'prefix2'   : '',
	}

	### Data
        dic_class['name'] = 'Data'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/' 
        corr_1 = correlation_3D.Correlation3D(dic_class)
	corr_1.set_error_on_covar_matrix('subsampling')
	xxx = numpy.arange(5000)
	error_data = corr_1._xi2D[:,:,2].flatten()
	plt.plot(xxx,corr_1._xi2D[:,:,2].flatten()/error_data,label='data subSampling',linewidth=2,alpha=0.7)

	### Raw
	dic_class['name'] = 'Raw'
	dic_simu['prefix'] = '_raw_from_JeanMarc'
	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Results'+dic_simu['prefix']+'/'
        corr_2 = correlation_3D.Correlation3D(dic_class) 
	corr_2.set_values_on_mean_simulation_or_realisation_type(dic_simu)
	plt.plot(xxx,10.*corr_2._xi2D[:,:,2].flatten()/error_data,label='Var of 100 raw mocks',linewidth=2,alpha=0.7)

	### Cooked
	for a in range(1):
		for b in range(10):
			dic_class['name'] = 'Cooked'
			dic_simu['prefix'] = ''
			dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_00'+str(a)+'/Simu_00'+str(b)+'/Results'+dic_simu['prefix']+'/'
        		corr_3 = correlation_3D.Correlation3D(dic_class)
			corr_3.set_error_on_covar_matrix('subsampling')
			if (a==0 and b==0): error = corr_3._xi2D[:,:,2]
			else: error += corr_3._xi2D[:,:,2]
	plt.plot(xxx,error.flatten()/10./error_data,label='One cooked mock subsampling',linewidth=2,alpha=0.7)

	corr_3.set_values_on_mean_simulation_or_realisation_type(dic_simu)
	plt.plot(xxx,10.*corr_3._xi2D[:,:,2].flatten()/error_data,label='Var of 100 cooked mocks',linewidth=2,alpha=0.7)

	plt.legend(fontsize=40, numpoints=1,ncol=1, loc=0)
	plt.xlabel('bin index')
	plt.ylabel('error / error data')
	plt.grid(True, which='both')
	plt.rc('font', **{'size':'40'})
	plt.tick_params(axis='both', which='major', labelsize=40)
	plt.linewidth = 40
	plt.show()

	return
def plot_compare_nb_correlations():

	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
                'nb_wedges'               : 3,
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '',
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'LYA - QSO'
	}
	dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '', #
		'prefix2'   : '',
	}

	### Data
        dic_class['name'] = 'Data'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/' 
        corr_1 = correlation_3D.Correlation3D(dic_class)
	xxx = numpy.arange(5000)
	error_data = corr_1._xi2D[:,:,2].flatten()
	plt.plot(xxx,corr_1._xi2D[:,:,2].flatten()/error_data,label='data subSampling',linewidth=2,alpha=0.7)

	### Cooked
	dic_class['name'] = 'Cooked'
	dic_simu['prefix'] = ''
	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Results'+dic_simu['prefix']+'/'
        corr_3 = correlation_3D.Correlation3D(dic_class)
	plt.plot(xxx,corr_3._xi2D[:,:,2].flatten()/error_data,label='One cooked mock',linewidth=2,alpha=0.7)

	plt.legend(fontsize=40, numpoints=1,ncol=1, loc=0)
	plt.xlabel('bin index')
	#plt.ylabel('variance (similar to Wick T0)')
	plt.ylabel('nb pairs / nb pairs data')
	plt.grid(True, which='both')
	plt.rc('font', **{'size':'40'})
	plt.tick_params(axis='both', which='major', labelsize=40)
	plt.linewidth = 40
	plt.show()

	return

#plot_compare_some_correlations()
#plotBosseBOSS()
plotOne()
#plotMany()
#BAOFIT_One_fit_data()
#BAOFIT_One_fit_simulation()
#BAOFIT_Many_simulations()
#BAOFIT_Compare_simulations()
#BAOFIT_Compare_simulations_pyLya()
#plot_compare_errors_correlations()
#plot_compare_nb_correlations()


### Send parameter scan
'''
dic_class = correlation_3D.raw_dic_class
dic_CAMB  = correlation_3D.raw_dic_CAMB
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
dic_class['name'] = 'Data'
corr = correlation_3D.Correlation3D(dic_class)
correlation_matrix_path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy'
corr.send_BAOFIT('subsampling',correlation_matrix_path, False, False, 10000)
'''

"""


### BAOFIT + Data

dic_class = correlation_3D.raw_dic_class
dic_CAMB  = correlation_3D.raw_dic_CAMB
dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
		'remove_residuales' : 'lambda_OBS',
		'index_multipole_max' : 4,
		'path_to_txt_file_folder': '',
		'correlation': 'q_f',
		'f1': 'LYA',
		'f2': '',
		'q1': 'QSO',
		'q2': '',
		'name' : 'simu'
}
dic_simu = {
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '_raw_from_JeanMarc',
		'prefix2'   : '__noDistortionMatrix_sigmaNL_0_pixelSize_3_15', #'__noDistortionMatrix_fit_from_10_to_180'
}
dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
path_to_BAOFIT = None
path_to_BAOFIT = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix']+'/BaoFit_q_f__LYA__QSO'+dic_simu['prefix2']+'/bao2D.'
print path_to_BAOFIT
#path_to_BAOFIT = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix']+'/BaoFit_f_f__LYA'+dic_simu['prefix2']+'/bao2D.'
path_to_mapping_we = dic_simu['path_to_simu']+'Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_we.npy'
path_to_mapping_1D = dic_simu['path_to_simu']+'Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_mapping_2D_to_1D.npy'
#path_to_mapping_we = dic_simu['path_to_simu']+'Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_we.npy'
#path_to_mapping_1D = dic_simu['path_to_simu']+'Results_raw_from_JeanMarc/xi_A_delta_delta_LYA_result_mapping_2D_to_1D.npy'

corr = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, None, path_to_BAOFIT)
corr.set_values_on_mean_simulation(dic_simu)


'''
dic_class = {
                'minXi': 0.,
                'maxXi': 200.,
                'nbBin': 50,
                'nbBinM': 25,
                'nb_Sub_Sampling': 80,
                'size_bin_calcul_s': 1.,
                'size_bin_calcul_m': 0.02,
                'remove_residuales' : 'lambda_OBS',
                'index_multipole_max' : 4,
                'path_to_txt_file_folder': '',
                'correlation': 'q_f',
                'f1': 'LYA',
                'f2': '',
                'q1': 'QSO',
                'q2': '',
                'name' : 'raw'
        }
dic_simu = {
                'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
                'nb_box' : 10,
                'nb_simu' : 10,
                'prefix' : '_raw_from_JeanMarc',
		'prefix2' : '__noDistortionMatrix_fit_from_10_to_180'
}
dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
path_to_BAOFIT = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix']+'/BaoFit_q_f__LYA__QSO'+dic_simu['prefix2']+'/bao2D.'
corr4 = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, None, path_to_BAOFIT)
#corr4.save_list_realisation_simulation(dic_class, dic_simu)
corr4.set_values_on_mean_simulation(dic_simu)



#corr = corr4
#for i in range(0,50):
#	corr.plot_slice_fit_2d(i,None,[corr4])
'''

'''
for i in numpy.arange(0,50):
	corr.plot_slice_fit_2d(i)
for i in numpy.arange(0,100):
	corr.plot_slice_fit_2d(None,i)
corr.print_results()
'''

#corr.print_results()
corr.plot_data_and_fit_1d(0,path_to_mapping_1D)
corr.plot_data_and_fit_1d(1,path_to_mapping_1D)
corr.plot_data_and_fit_1d(2,path_to_mapping_1D)
corr.plot_data_and_fit_we(0,path_to_mapping_we)
corr.plot_data_and_fit_we(1,path_to_mapping_we)
corr.plot_data_and_fit_we(2,path_to_mapping_we)

corr.plot_fit_2d(0)
corr.plot_fit_2d(1)
corr.plot_fit_2d(2)
corr.plot_residuals_2d(0)
corr.plot_residuals_2d(1)
corr.plot_residuals_2d(2)
corr.plot_histo_residuals()



#corr.plot_chi2_scan(50,50, [0.98,1.02,0.98,1.02], False, False)
#corr.plot_chi2_scan(50,50, [0.998,1.002,0.998,1.002], False, False)



### BAOFIT + One simulations
'''
dic_class = correlation_3D.raw_dic_class
dic_CAMB  = correlation_3D.raw_dic_CAMB
index_parameter = annalyse_BAOFIT.raw_index_parameter
index_parameter['alpha_paral'] = 9
index_parameter['alpha_perp']  = 10
dic_class['name'] = 'Simulation'
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_009/Results/'
path_to_BAOFIT = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_009/BaoFit_q_f_covFromFit/bao2D.'
corr = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, index_parameter, path_to_BAOFIT)
corr.print_results()
corr.plot_chi2_scan(50,50)
'''

"""
