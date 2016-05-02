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


def plotOne():


	dic_class = {
		'minXi': 0.,
		'maxXi': 200.,
		'nbBin': 50,
		'nbBinM': 25,
		'nb_Sub_Sampling': 80,
		'size_bin_calcul_s': 1.,
		'size_bin_calcul_m': 0.02,
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
		'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/',
		'nb_box' : 10,
		'nb_simu' : 10,
		'prefix' : '_devide_instead_of_removing_PureRaw'
	}
	dic_send_BAOFIT = {
		'realisation_type'        : 'Simu_stack',
		#'correlation_matrix_path' : None, #'/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results/xi_delta_QSO_LYA_QSO_result_cor_mean_subsampling_from_fit_2D.npy',
		'correlation_matrix_path' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Results_raw_from_JeanMarc/xi_delta_QSO_LYA_QSO_result_cor_2D_from_fit.npy',
		'saving_data'             : True,
		'scan'                    : False,
		'toyMC'                   : 0,
		'dist_matrix'             : False,
		'only_diagonal'           : False,
		'with_metals_templates'   : False,
		'get_param_from_file'     : False,
		'dic_simu'                : dic_simu
	}
	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Results'+dic_simu['prefix']+'/'
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
	corr = correlation_3D.Correlation3D(dic_class)
	#corr.save_list_realisation_simulation(dic_class, dic_simu)
	corr.set_values_on_mean_simulation(dic_simu)
	#corr.save_mean_realisation_simulation_correlation_matrix(dic_simu)
	#corr.write_metal_model(True)
	#corr.read_data_from_BAOFIT_data_file('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/TESTS/Nicolas_Correlation/xcf_dr12_v5_8_guy_c2_baseline_projected')
	

	path_to_load = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_devide_instead_of_removing_PureRaw/xi_A_delta_delta_LYA_result_cov_2D.npy'
	path_to_save = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_devide_instead_of_removing_PureRaw/xi_A_delta_delta_LYA_result_cor_2D_from_fit'
	myTools.save_fited_cor_from_cov(path_to_load, path_to_save,50,50,'f_f')
	

	"""
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
	"""

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

	'''
	xiMu, xiWe, xi1D, xi2D = corr.read_metal_model( 'LYA' )	
	path = corr._path_to_txt_file_folder + corr._prefix + '_distortionMatrix_2D_'+ corr._middlefix + '.txt'
	print path
	matrix = numpy.loadtxt(path)

	corr._xi2D[:,:,0] = xi2D[:,:,0,2]
	corr._xi2D[:,:,1] = xi2D[:,:,1,2]

	corr.plot_2d(0)
        corr.plot_2d(1)
        corr.plot_2d(2)

	#matrix[matrix==0.] = numpy.float('nan')
	#for l in numpy.arange(corr._nbBin2D):
	#	matrix[l,l] = numpy.float('nan')
	#myTools.plot2D(matrix)

	yyy_Camb = numpy.dot(matrix,xi2D[:,:,1,2].flatten())
	corr._xi2D[:,:,1] = myTools.convert1DTo2D(yyy_Camb,corr._nbBinX2D, corr._nbBinY2D)
	corr.plot_2d(0)
        corr.plot_2d(1)
        corr.plot_2d(2)
	'''

	#corr.save_list_realisation('subsampling', 80)
	#corr.plot_cov_cor_matrix('subsampling','1D')
	### Fit
	#corr.set_error_on_covar_matrix('subsampling')

	#corr.send_BAOFIT(dic_send_BAOFIT)

	"""
	print corr._meanZ
	corr.plot_grid(0,True)
	corr.plot_grid(1,True)
	corr.plot_grid(2,True)
	

	c0, c2, c4 = CAMB1.get_multipol_coef(-0.168777,1.364)
	print c0, c2, c4
	
	
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
	"""

	'''
	for ii in numpy.arange(0,50):
		corr.plot_slice_2d(ii,[corr2])
	for ii in numpy.arange(0,100):
		corr.plot_slice_2d(None,ii,[corr2])
	'''
	
	corr.plot_1d(0)
	corr.plot_1d(1)
	corr.plot_1d(2)
	corr.plot_we(0)
	corr.plot_we(1)
	corr.plot_we(2)
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
	

	'''
	corr.plot_cov_cor_matrix_different_method(['subsampling'],'1D',[ ['Simulation','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_withRand_XYZ_nicolasEstimator/xi_delta_QSO_LYA_QSO_result_cov_1D.npy'] ] )
	corr.plot_cov_cor_matrix_different_method(['subsampling'],'2D',[ ['Simulation','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_withRand_XYZ_nicolasEstimator/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'] ] )
	'''
	#corr.plot_cov_cor_matrix_different_method(['subsampling'],'2D',[ ['simu','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results/xi_delta_QSO_LYA_QSO_result_cov_2D.npy'], [ 'not \, removing', '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_LYA_QSO_subsampling_cov_2D.npy'] ] )
	

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
		'remove_residuales' : 'lambda_OBS',
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
		'prefix' : '_raw_from_JeanMarc'
	}
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_000/Simu_000/Results_raw_from_JeanMarc/'
	corr4 = correlation_3D.Correlation3D(dic_class)
	#corr4.save_list_realisation_simulation(dic_class, dic_simu)
        corr4.set_values_on_mean_simulation(dic_simu)

	dic_simu['path_to_simu'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/'
	dic_simu['prefix'] = '_devide_instead_of_removing_PureRaw'
	dic_class['name'] = 'raw \, before'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Results_devide_instead_of_removing_PureRaw/'
	corr3 = correlation_3D.Correlation3D(dic_class)
        #corr.save_list_realisation_simulation(dic_class, dic_simu)
        corr3.set_values_on_mean_simulation(dic_simu)

	#dic_class['name'] = 'simu'
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1597/Box_000/Simu_000/Results/'
	#corr6 = correlation_3D.Correlation3D(dic_class)

	###
	dic_class['q1'] = 'QSO'
	dic_class['remove_residuales'] = 'lambda_OBS'
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
	dic_class['name'] = 'Data'
	corr5 = correlation_3D.Correlation3D(dic_class)
	corr5.set_error_on_covar_matrix('subsampling')
	print corr5._meanZ

	corr5.plot_1d(0,[corr3,corr4],False)
	corr5.plot_1d(1,[corr3,corr4],False)
	corr5.plot_1d(2,[corr3,corr4],False)
	corr5.plot_we(0,[corr4],False)
	corr5.plot_we(1,[corr4],False)
	corr5.plot_we(2,[corr4],False)


	for i in numpy.arange(0,50):
		corr.plot_slice_2d(i,None,[corr6])


	return
def look_result_data():

	

	return

#plotBosseBOSS()
plotOne()
#plotMany()


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
		'prefix' : '_raw_from_JeanMarc'
}
dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
path_to_BAOFIT = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix']+'/BaoFit_q_f__LYA__QSO/bao2D.'

corr = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, None,path_to_BAOFIT)
corr.set_values_on_mean_simulation(dic_simu)
'''
for i in numpy.arange(0,50):
	corr.plot_slice_fit_2d(i)
for i in numpy.arange(0,100):
	corr.plot_slice_fit_2d(None,i)
corr.print_results()
'''

corr.plot_data_and_fit_1d(0)
corr.plot_data_and_fit_1d(1)
corr.plot_data_and_fit_1d(2)
corr.plot_data_and_fit_we(0)
corr.plot_data_and_fit_we(1)
corr.plot_data_and_fit_we(2)
corr.plot_fit_2d(0)
corr.plot_fit_2d(1)
corr.plot_fit_2d(2)
corr.plot_residuals_2d(0)
corr.plot_residuals_2d(1)
corr.plot_residuals_2d(2)
corr.plot_histo_residuals()

corr.plot_chi2_scan(50,50, [0.98,1.02,0.98,1.02], False, False)
#corr.plot_chi2_scan(50,50, [0.998,1.002,0.998,1.002], False, False)
"""
"""

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



### BAOFIT + Many simulations

dic_class = {
	'minXi': 0.,
	'maxXi': 200.,
	'nbBin': 50,
	'nbBinM': 25,
	'nb_Sub_Sampling': 80,
	'size_bin_calcul_s': 1.,
	'size_bin_calcul_m': 0.02,
	'correlation': 'q_f',
	'remove_residuales' : '',
	'path_to_txt_file_folder': '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/',
	'f1': 'LYA',
	'f2': '',
	'q1': 'QSO',
	'q2': '',
	'name' : 'raw'
}
dic_simu = {
	'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/',
	'nb_box' : 10,
	'nb_simu' : 10,
	'prefix' : '_devide_instead_of_removing_PureRaw',
	'prefix2' : ''
}
aMB = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,None,dic_simu,False)

###
dic_class['remove_residuales'] = 'lambda_OBS'
dic_simu['prefix'] = '_nicolasEstimator'
dic_simu['prefix2'] = '__withMetalsTemplates'
dic_class['name']   = 'full'
aMB3 = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,None,dic_simu,False)

aMB.print_results_for_chi2([aMB3])

aMB.plot_scatter_hist(0,4,[aMB3])
aMB.plot_scatter_hist(0,5,[aMB3])
aMB.plot_scatter_hist(0,11,[aMB3])
aMB.plot_scatter_hist(0,12,[aMB3])
aMB.plot_scatter_hist(4,5,[aMB3])
aMB.plot_scatter_hist(4,11,[aMB3])
aMB.plot_scatter_hist(4,12,[aMB3])
aMB.plot_scatter_hist(5,11,[aMB3])
aMB.plot_scatter_hist(5,12,[aMB3])
aMB.plot_scatter_hist(11,12,[aMB3])
#aMB.plot_histo_residuals([aMB3])

aMB.print_results(0,[aMB3])
aMB.print_results(4,[aMB3])
aMB.print_results(5,[aMB3])
aMB.print_results(11,[aMB3])
aMB.print_results(12,[aMB3])

#aMB.set_mean_data()
#print aMB.get_best_fit()
#corr.plot_chi2_scan(100,100, [0.5,1.5,0.5,1.5], False, False, aMB.get_best_fit())


aMB.plot_1d(0)
aMB.plot_1d(1)
aMB.plot_1d(2)
aMB.plot_we(0)
aMB.plot_we(1)
aMB.plot_we(2)
aMB.plot_2d(0)
aMB.plot_2d(1)
aMB.plot_2d(2)
aMB.plot_mu(0)
aMB.plot_mu(1)
aMB.plot_mu(2)
aMB.plot_multipol(0)
aMB.plot_multipol(1)
aMB.plot_multipol(2)
aMB.plot_chi2_scan(50,50,None,corr.get_best_fit())

'''
corr.plot_multipol(0,[aMB._mean_data])
corr.plot_multipol(1,[aMB._mean_data])
corr.plot_multipol(2,[aMB._mean_data])
corr.plot_1d(0,[aMB._mean_data])
corr.plot_1d(1,[aMB._mean_data])
corr.plot_1d(2,[aMB._mean_data])
corr.plot_we(0,[aMB._mean_data])
corr.plot_we(1,[aMB._mean_data])
corr.plot_we(2,[aMB._mean_data])
'''








'''
dic_CAMB = {
	'mulpol_index' : 0,
	'start_fit'   : 20.,
	'end_fit'     : 70.,
	'b' : -1.,
	'roof' : 0.,
	'fix_roof_nul' : True,
	'guess_b' : False,
	'min_for_guess' : 20.,
	'max_for_guess' : 50.,
}
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
	'q1': 'ALL_OBJECTS',
	'q2': 'QSO',
	'name' : 'Data'
}
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
dic_class['name'] = 'Data'
corr = correlation_3D.Correlation3D(dic_class)

#path = [ [ 'Mocks', '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D.npy'],
#	['< Mock \, subsampling >','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_meanSubSampling.npy'] ]
#corr.plot_cov_cor_matrix_different_method( [ 'subsampling', 'shuffleForest', 'shuffleQSO', 'randomQSO'], '2D', path)

corr.set_error_on_covar_matrix('subsampling')

dic_CAMB = corr.fit_CAMB()
corr.plot_CAMB(None,None,0)
corr.plot_CAMB(None,None,1)
corr.plot_CAMB(None,None,2)
corr.plot_multipol(0)
corr.plot_multipol(1)
corr.plot_multipol(2)
corr.plot_we(0)
corr.plot_we(1)
corr.plot_we(2)
corr.plot_1d(0)
corr.plot_1d(1)
corr.plot_1d(2)
corr.plot_we(0)
corr.plot_we(1)
corr.plot_we(2)
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
#corr.save_list_realisation('subsampling', 80)
#corr.plot_cov_cor_matrix('subsampling','1D')
#correlation_matrix_path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy'
#corr.send_BAOFIT('subsampling',correlation_matrix_path,False,True)
'''


### List function Correlation3D
'''
corr.plot_1d(0)
corr.plot_1d(1)
corr.plot_1d(2)
corr.plot_We(0)
corr.plot_We(1)
corr.plot_We(2)
corr.plot_2d(0)
corr.plot_2d(1)
corr.plot_2d(2)
corr.plot_Mu(0)
corr.plot_Mu(1)
corr.plot_Mu(2)
corr.plot_map_sub_sampling()
dic_CAMB = corr.fit_CAMB(dic_CAMB)
corr.plot_CAMB(dic_CAMB)

corr.plot_distortion_matrix()
corr.apply_distortion_matrix()

corr.plot_map_sub_sampling()
corr.plot_distortion_matrix()
corr.apply_distortion_matrix()
corr.plot_cov_cor_matrix('subsampling','1D')
corr.save_list_realisation('subsampling', 80)
corr.plot_cov_cor_matrix('subsampling','1D')
corr.set_error_on_covar_matrix('subsampling')


### Fit
correlation_matrix_path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy'
corr.send_BAOFIT('subsampling',correlation_matrix_path)
'''


### List function AnnalyseBAOFIT
'''
corr.print_results()
corr.plot_data_and_fit_1d(0)
corr.plot_data_and_fit_1d(1)
corr.plot_data_and_fit_1d(2)
corr.plot_data_and_fit_we(0)
corr.plot_data_and_fit_we(1)
corr.plot_data_and_fit_we(2)
corr.plot_fit_2d(0)
corr.plot_fit_2d(1)
corr.plot_fit_2d(2)
corr.plot_residuals_2d(0)
corr.plot_residuals_2d(1)
corr.plot_residuals_2d(2)
corr.plot_histo_residuals()
corr.plot_chi2_scan(100,100, [0.5,1.5,0.5,1.5], True, False)
'''



"""












