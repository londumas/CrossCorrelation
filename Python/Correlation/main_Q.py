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
import cosmolopy.perturbation as cp

### Perso lib
import const_delta
import myTools
import correlation_3D_Q
import annalyse_BAOFIT_Q
import annalyse_many_BAOFIT
import CAMB
import transform_CAMB_by_Bias



dic_class = {
	'minXi': 0.,
	'maxXi': 200.,
	'nbBin': 50,
	'nbBinM': 25,
	'nb_Sub_Sampling': 80,
	'size_bin_calcul_s': 1.,
	'size_bin_calcul_m': 0.02,
	'correlation': 'q_q',
	'path_to_txt_file_folder': '',
	'remove_residuales' : '',
	'f1': 'LYA',
	'f2': '',
	'q1': 'QSO',
	'q2': '',
	'name' : 'Mean \, simulation'
}
dic_Q = {
	'nb_random' : 100,
	'estimator' : 'LS',
	'path_to_cat' : 'NOTHING',
	'load_from_txt' : False,
	'sufix' : 'QSO_withRSD'
}
dic_CAMB = {
	'z' : 2.25,
	'h_0' : 0.70,
	'omega_matter_0' : 0.27,
	'omega_lambda_0' : 0.73,
	'source'         : 'CAMB_Mocks_me'
}
dic_CAMB_fit = {
	'mulpol_index' : 0,
	'start_fit'   : 40.,
	'end_fit'     : 180.,
	'b' : 1.,
	'roof' : 0.,
	'fix_roof_nul' : True,
	'guess_b' : False,
	'min_for_guess' : 20.,
	'max_for_guess' : 50.,
	'CAMB' : CAMB.CAMB(dic_CAMB)
}
dic_simu = {
	'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/',
	'nb_box' : 10,
	'nb_simu' : 10,
	'prefix'    : '_raw_from_JeanMarc',
	'prefix2'   : ''
}
dic_send_BAOFIT = {
	'realisation_type'        : 'Simu_stack',
	'correlation_matrix_path' : None,
	'saving_data' : False,
	'scan' : False,
	'toyMC' : 0,
	'dist_matrix' : False,
	'only_diagonal' : True,
	'with_metals_templates' : False,
	'get_param_from_file' : False,
	'dic_simu'              : dic_simu
}

if (True):

	
	
	i = sys.argv[1]
	j = sys.argv[2]
	
	#i = 0
	#j = 0
	raw = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/'
	dic_class['path_to_txt_file_folder'] = raw+'Results'+dic_simu['prefix']+'/'
	dic_class['name'] = 'random \, XYZ'
	dic_Q['path_to_cat']                 = raw+'Data/'+dic_Q['sufix']+'.fits'
	corr = correlation_3D_Q.Correlation3DQ(dic_class,dic_Q)
	#corr.save_list_realisation_simulation(dic_class, dic_Q, dic_simu)
	corr.set_values_on_mean_simulation(dic_simu)

	'''
	path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results_q_q__noSameCell_randXYZ_withRSD/xi_QSO_QSO_result_cov_1D.npy'
	cov = numpy.load(path_to_cov)
	myTools.plot2D(cov)
	myTools.plot2D(myTools.getCorrelationMatrix(cov))
	myTools.plotCovar( [cov], ['a'] )
	'''

	"""
	### Look at the covariance matrix
	path_to_load = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_QSO_QSO_result_cov_2D.npy'
	path_to_save = dic_simu['path_to_simu']+'Results' + dic_simu['prefix'] + '/xi_QSO_QSO_result_cor_2D_from_fit'
	print path_to_load
	print path_to_save
	myTools.save_fited_cor_from_cov(path_to_load, path_to_save,50,50,'q_q')
	"""

	print corr._meanZ

	### Test CAMB and CAMB transformed
	z_ref = corr._meanZ
	f  = cp.fgrowth(z_ref,0.27)
	f2 = cp.fgrowth(2.5,0.27)

	dic = {
	'sigma_gauss'     : 1.96769,
	'z'               : 0.,
	'a'               : numpy.inf,
	'gamma'           : numpy.inf,
	'cut_QSO'         : 2.3793413839*1.96769,
	'h_0'             : 0.70,
	'omega_matter_0'  : 0.27,
	'omega_lambda_0'  : 0.73,
	'correlation'     : 'q_q',
	'source_for_CAMB' : 'CAMB_Mocks_me',
	'mock_type'       : 'JMC',
	}
	a = transform_CAMB_by_Bias.TRANSFORM_CAMB_BIAS(dic)

	### CAMB
	xxx_CAMB = a._CAMB._xi0[:,0]
	yyy_CAMB = a._CAMB._xi0[:,1]*f*f
	b_camb = 3.5728433556050487
	b_camb_effectif = a._CAMB.get_bias_effectif_knowing_bias_QSO(b_camb,z_ref)

	
	print '  b_camb          = ', b_camb
	print '  b_camb_effectif = ', b_camb_effectif
	print '  c0 = ', b_camb_effectif*b_camb_effectif

	### Transformed CAMB
	r_max = 200.
	xi_r = numpy.arange(0.0001,r_max,0.1)
	transformed_xi_0 = a.get_xi0_transformed(xi_r,0.)
	xxx_CAMB2 = a._CAMB._xi0[:,0]
	yyy_CAMB2 = a._CAMB._xi0[:,1]*f2*f2
	b_trans = numpy.sqrt(numpy.interp(180.,xi_r,transformed_xi_0) / numpy.interp(180.,xxx_CAMB2,yyy_CAMB2))
	b_trans_effectif = a._CAMB.get_bias_effectif_knowing_bias_QSO(b_trans,z_ref)
	print '  b_trans          = ', b_trans
	print '  b_trans_effectif = ', b_trans_effectif
	#b_trans_effectif *= 1.03

	###
	
	for i in range(0,3):
		###
		coef = numpy.power(xxx_CAMB,i)
		plt.plot(xxx_CAMB,coef*yyy_CAMB*b_camb_effectif*b_camb_effectif, label=r'$CAMB$', markersize=10,linewidth=2)
		###
		coef = numpy.power(xi_r,i)*numpy.power(b_trans_effectif/b_trans,2.)
		plt.errorbar(xi_r, coef*transformed_xi_0, label=r'$CAMB \, transformed$', markersize=10,linewidth=2)
		###
		coef = numpy.power(corr._xiMul[:,0,0],i)
		plt.errorbar(corr._xiMul[:,0,0],coef*corr._xiMul[:,0,1], yerr=coef*corr._xiMul[:,0,2], fmt='o', label=r'$Simulation$', markersize=10,linewidth=2)
		###
		plt.xlim([0.,200.])
		myTools.deal_with_plot(False,False,True)
		plt.show()
	
	size = xi_r.size
	xi0 = numpy.zeros( shape=(size,2) )
	xi0[:,0] = xi_r
	xi0[:,1] = transformed_xi_0
	xi2, xi4 = a._CAMB.get_xi2_and_xi4_from_xi0(xi0)
	c0, c2, c4 = a._CAMB.get_multipol_coef_knowing_bias_QSO(b_camb,z_ref)
	print c0, c2, c4



	for i in range(0,3):
		###
		coef = numpy.power(a._CAMB._xi2[:,0],i)
		plt.plot(a._CAMB._xi2[:,0],coef*a._CAMB._xi2[:,1]*c2*f*f, label=r'$CAMB$', markersize=10,linewidth=2)
		###
		coef = numpy.power(xi2[:,0],i)
		plt.plot(xi2[:,0],coef*xi2[:,1]/3., label=r'$CAMB \, transformed$', markersize=10,linewidth=2)
		###
		coef = numpy.power(corr._xiMul[:,2,0],i)
		plt.errorbar(corr._xiMul[:,2,0],coef*corr._xiMul[:,2,1], yerr=coef*corr._xiMul[:,2,2], fmt='o', label=r'$Simulation$', markersize=10,linewidth=2)

		###
		plt.xlim([0.,200.])
		myTools.deal_with_plot(False,False,True)
		plt.show()

	for i in range(0,3):
		###
		coef = numpy.power(a._CAMB._xi4[:,0],i)
		plt.plot(a._CAMB._xi4[:,0],coef*a._CAMB._xi4[:,1]*c4*f*f, label=r'$CAMB$', markersize=10,linewidth=2)
		###
		coef = numpy.power(xi4[:,0],i)
		plt.plot(xi4[:,0],coef*xi4[:,1]/33., label=r'$CAMB \, transformed$', markersize=10,linewidth=2)
		###
		coef = numpy.power(corr._xiMul[:,4,0],i)
		plt.errorbar(corr._xiMul[:,4,0],coef*corr._xiMul[:,4,1], yerr=coef*corr._xiMul[:,4,2], fmt='o', label=r'$Simulation$', markersize=10,linewidth=2)

		###
		plt.xlim([0.,200.])
		myTools.deal_with_plot(False,False,True)
		plt.show()
	


	"""
	dic_CAMB['mulpol_index'] = 0
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,0,:],dic_CAMB_fit,False)
	print numpy.sqrt(dic_CAMB_fit['b'])
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB_fit,0,False)
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB_fit,1,False)
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB_fit,2,False)
	
	dic_CAMB['mulpol_index'] = 2
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,2,:],dic_CAMB_fit,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB_fit,0,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB_fit,1,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB_fit,2,False)

	dic_CAMB['mulpol_index'] = 4
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,4,:],dic_CAMB_fit,False)
	corr.plot_CAMB(corr._xiMul[:,4,:],dic_CAMB_fit,0,False)
	corr.plot_CAMB(corr._xiMul[:,4,:],dic_CAMB_fit,1,False)
	corr.plot_CAMB(corr._xiMul[:,4,:],dic_CAMB_fit,2,False)
	
	
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
	
	corr.plot_grid(0,True)
	corr.plot_grid(1,True)
	corr.plot_grid(2,True)
	"""
	#corr.send_BAOFIT(dic_send_BAOFIT)
	

'''
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Results/BaoFit_q_q__QSO/bao2D.fit.dat')
for i in [0,1,2]:
	coef = numpy.power(data[:,0],i)
	plt.plot( data[:,0], coef*data[:,1])
	plt.plot( data[:,0], coef*data[:,2])
	plt.plot( data[:,0], coef*data[:,3])
	plt.show()
'''

# ---------------------------------
### For Fit
# ---------------------------------
if (False):
	i = sys.argv[1]
	j = sys.argv[2]

	raw = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/'
	dic_class['path_to_txt_file_folder'] = raw+'Results'+dic_simu['prefix']+'/'
	dic_class['name'] = 'random \, XYZ'
	dic_Q['path_to_cat']                 = raw+'Data/'+dic_Q['sufix']+'.fits'

	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results'+dic_simu['prefix']+'/'
	path_to_BAOFIT = dic_simu['path_to_simu'] + 'Results' + dic_simu['prefix']+'/BaoFit_q_q__QSO/bao2D.'
	corr = annalyse_BAOFIT_Q.AnnalyseBAOFIT(dic_class, None,path_to_BAOFIT,dic_Q)
	corr.set_values_on_mean_simulation(dic_simu)
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
	'''
	corr.plot_chi2_scan(50,50, [0.98,1.02,0.98,1.02], False, False, None)





# ---------------------------------
### For Many Fit
# ---------------------------------
if (False):

	i = sys.argv[1]
	j = sys.argv[2]

	raw = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/'
	dic_class['path_to_txt_file_folder'] = raw+'Results'+dic_simu['prefix']+'/'
	dic_class['name'] = 'qso-qso'
	dic_Q['path_to_cat']                 = raw+'Data/'+dic_Q['sufix']+'.fits'
	aMB = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class,None,dic_Q,dic_simu,False)

	#aMB.print_results_for_chi2()
	#aMB.print_results(0)
	#aMB.print_results(1)
	#aMB.print_results(8)
	#aMB.print_results(9)

	aMB.plot_scatter_hist(0,1,None,[ [0.261,3.70*(1.+0.261)] ])
	aMB.plot_scatter_hist(0,8)
	aMB.plot_scatter_hist(0,9)
	aMB.plot_scatter_hist(1,8)
	aMB.plot_scatter_hist(1,9)
	aMB.plot_scatter_hist(8,9)
	aMB.plot_histo_residuals()















