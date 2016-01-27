# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_3D.py
#

### Perso lib
import correlation_3D
import annalyse_BAOFIT
import annalyse_many_BAOFIT






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
'''
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
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
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
corr.plot_chi2_scan(100,100, [0.5,1.5,0.5,1.5], True, False)
'''




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
'''
dic_class = correlation_3D.raw_dic_class
index_parameter = annalyse_BAOFIT.raw_index_parameter
index_parameter['alpha_paral'] = 9
index_parameter['alpha_perp']  = 10
aMB = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class, index_parameter, '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/', 1,10)
aMB.plot_histo_residuals()
#aMB.print_results()
aMB.plot_chi2_scan(50,50)
'''






'''
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
dic_CAMB  = correlation_3D.raw_dic_CAMB
dic_class['path_to_txt_file_folder'] = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
dic_class['name'] = 'Data'
corr = correlation_3D.Correlation3D(dic_class)

#path = [ [ 'Mocks', '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D.npy'],
#	['< Mock \, subsampling >','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_meanSubSampling.npy'] ]
#	corr.plot_cov_cor_matrix_different_method( [ 'subsampling', 'shuffleForest', 'shuffleQSO', 'randomQSO'], '2D', path)
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

#corr.save_list_realisation('subsampling', 80)
#corr.plot_cov_cor_matrix('subsampling','1D')
correlation_matrix_path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy'
corr.send_BAOFIT('subsampling',correlation_matrix_path,False,True)
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






















