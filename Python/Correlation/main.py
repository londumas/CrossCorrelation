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

### Perso lib
import const_delta
import myTools
import correlation_3D
import annalyse_BAOFIT
import annalyse_many_BAOFIT

import matplotlib.pyplot as plt


#all_t = [ '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits',
#'/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits',
#'/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits']
#myTools.append_files(all_t, '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/FitsFile_DR12_reOBS_eBOSS_noCoADD_Guy/DR12_primery/DR12_primery_reOBS_eBOSS_noCoADD.fits')


def plotOne():

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
		'f2': 'a',
		'q1': 'QSO',
		'q2': 'a',
		'name' : 'Data'
	}
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/Box_000/Simu_000/Results/'
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/'
	#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/TESTS/'
	dic_class['name'] = 'Data'
	corr = correlation_3D.Correlation3D(dic_class)
	#corr.write_metal_model()
	#corr.read_data_from_BAOFIT_data_file('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/TESTS/Nicolas_Correlation/xcf_dr12_v5_8_guy_c2_baseline_projected')

	'''
	a = 0
	xiMu, xiWe, xi1D, xi2D = corr.read_metal_model(const_delta.LYB_lines_names[-2])
	
	for i in numpy.arange(2,const_delta.LYB_lines_names.size):
		tmpxiMu, tmpxiWe, tmpxi1D, tmpxi2D = corr.read_metal_model(const_delta.LYB_lines_names[-1-i])
		if (tmpxi2D[:,:,1,0][ (tmpxi2D[:,:,1,0]!=0.) ].size == 0): break

		coef = 1.
		if ( const_delta.LYB_lines_names[-1-i]=="LYA" ): coef = 100.
		tmpxi2D[:,:,1,:] *= coef

		a += 1
		xi2D[:,:,1,:] += tmpxi2D[:,:,1,:]
		corr._xi2D[:,:,:] = tmpxi2D[:,:,:,0]
		corr.plot_2d(0)
	
	corr._xi2D[:,:,:] = xi2D[:,:,:,0]
	corr.plot_2d(0)
	corr.plot_2d(1)
	corr.plot_2d(2)
	'''
	#corr.save_list_realisation('subsampling', 80)
	#corr.plot_cov_cor_matrix('subsampling','1D')
	### Fit
	#correlation_matrix_path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy'
	#corr.send_BAOFIT('subsampling',correlation_matrix_path,False)
	#corr.set_error_on_covar_matrix('subsampling')
	#corr.fit_CAMB_2d('subsampling', correlation_matrix_path)

	dic_CAMB = corr.fit_CAMB(None,dic_CAMB,False)
	corr.plot_CAMB(None,dic_CAMB,0,False)
	corr.plot_CAMB(None,dic_CAMB,1,False)
	corr.plot_CAMB(None,dic_CAMB,2,False)

	
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

	list_corr = []
	#ALL_OBJECTS
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
	### Boss
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt/bao/hdumasde/Mock_JMLG/v1563/Box_000/Simu_000/Results/'
	dic_class['name'] = 'no \, projected'
	corr = correlation_3D.Correlation3D(dic_class)
	#corr.set_error_on_covar_matrix('subsampling')
	list_corr += [corr]
	### Boss
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
	dic_class['name'] = 'projected'
	corr = correlation_3D.Correlation3D(dic_class)
	#corr.set_error_on_covar_matrix('subsampling')
	list_corr += [corr]

	list_corr[0].plot_1d(0, list_corr[1:])
	list_corr[0].plot_1d(1, list_corr[1:])
	list_corr[0].plot_1d(2, list_corr[1:])

	'''
	list_corr[0].plot_we(0, list_corr[1:])
	list_corr[0].plot_we(1, list_corr[1:])
	list_corr[0].plot_we(2, list_corr[1:])
	list_corr[0].plot_1d_with_residuals(0, list_corr[1:])
	list_corr[0].plot_1d_with_residuals(1, list_corr[1:])
	list_corr[0].plot_1d_with_residuals(2, list_corr[1:])
	list_corr[0].plot_1d(0, list_corr[1:])
	list_corr[0].plot_1d(1, list_corr[1:])
	list_corr[0].plot_1d(2, list_corr[1:])
	list_corr[0].plot_we(0, list_corr[1:])
	list_corr[0].plot_we(1, list_corr[1:])
	list_corr[0].plot_we(2, list_corr[1:])
	list_corr[0].plot_2d(0)
	list_corr[0].plot_2d(1)
	list_corr[0].plot_2d(2)
	'''
	list_corr[0]._xi2D[:,:,1] = numpy.sqrt(2.)*(list_corr[1]._xi2D[:,:,1]-list_corr[0]._xi2D[:,:,1])/(numpy.sqrt( list_corr[1]._xi2D[:,:,2]*list_corr[1]._xi2D[:,:,2] + list_corr[0]._xi2D[:,:,2]*list_corr[0]._xi2D[:,:,2]))
	#list_corr[0]._xi2D[:,:,1] = (list_corr[1]._xi2D[:,:,1]-list_corr[0]._xi2D[:,:,1])

	plt.hist( list_corr[0]._xi2D[:,:,1].flatten(), bins=100 )
	plt.show()

	list_corr[0].plot_2d(0)
	list_corr[0].plot_2d(1)
	list_corr[0].plot_2d(2)
	list_corr[0].plot_slice_2d(1)
	list_corr[0].plot_slice_2d(2)
	list_corr[0].plot_slice_2d(3)

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

### BAOFIT + Data

dic_class = correlation_3D.raw_dic_class
dic_CAMB  = correlation_3D.raw_dic_CAMB
index_parameter = annalyse_BAOFIT.raw_index_parameter
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
dic_class['name'] = 'Data'
#dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation_test_soved_shift_withMetals_withError_with_template/Box_000/Simu_000/Results_NicolasDistortion/'
corr = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, index_parameter)
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
corr.plot_chi2_scan(100,100, [0.5,1.5,0.5,1.5], False, False)





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

if (True):
	index_parameter = annalyse_BAOFIT.raw_index_parameter
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
	### Boss
	dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
	dic_class['name'] = 'Data'
	corr = annalyse_BAOFIT.AnnalyseBAOFIT(dic_class, index_parameter)
	corr.set_error_on_covar_matrix('subsampling')

dic_class = correlation_3D.raw_dic_class
index_parameter = annalyse_BAOFIT.raw_index_parameter
index_parameter['alpha_paral'] = 9
index_parameter['alpha_perp']  = 10
aMB = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class, index_parameter, '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/', 1,1)
#aMB.plot_histo_residuals()
#aMB.save_list_realisation()
#aMB.print_results()
aMB.set_mean_data()
'''
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
'''
aMB.plot_chi2_scan(50,50,None,corr.get_best_fit())


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






















