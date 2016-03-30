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

### Perso lib
import const_delta
import myTools
import correlation_3D_Q
import annalyse_BAOFIT_Q
import annalyse_many_BAOFIT

import matplotlib.pyplot as plt



dic_class = {
	'minXi': 0.,
	'maxXi': 200.,
	'nbBin': 50,
	'nbBinM': 25,
	'nb_Sub_Sampling': 80,
	'size_bin_calcul_s': 1.,
	'size_bin_calcul_m': 0.02,
	'correlation': 'q_q',
	'path_to_txt_file_folder': 'NOTHING',
	'f1': 'LYA',
	'f2': 'LYA',
	'q1': 'QSO',
	'q2': 'QSO',
	'name' : 'Mean \, simulation'
}
dic_Q = {
	'nb_random' : 100,
	'estimator' : 'LS',
	'path_to_cat' : 'NOTHING',
	'load_from_txt' : False
}
dic_CAMB = {
	'mulpol_index' : 0,
	'start_fit'   : 60.,
	'end_fit'     : 180.,
	'b' : -1.,
	'roof' : 0.,
	'fix_roof_nul' : True,
	'guess_b' : False,
	'min_for_guess' : 20.,
	'max_for_guess' : 50.,
	'chi^{2}' : -1.
}
dic_simu = {
	'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/',
	'nb_box' : 10,
	'nb_simu' : 10,
	'projected' : False,
	'with_metals_templates' : False,
	'raw' : False,
	'pure_raw' : False,
}


if (True):

	
	'''
	i = sys.argv[1]
	j = sys.argv[2]
	'''
	i = 0
	j = 0
	raw = dic_simu['path_to_simu'] + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/'
	dic_class['path_to_txt_file_folder'] = raw+'Results/'
	dic_Q['path_to_cat']                 = raw+'Data/QSO_withRSD.fits'
	corr = correlation_3D_Q.Correlation3DQ(dic_class,dic_Q)
	#corr.save_list_realisation_simulation(dic_class, dic_Q, dic_simu)
	corr.set_values_on_mean_simulation(dic_simu)

	'''
	path_to_cov = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Results/xi_QSO_QSO_result_cov_2D.npy'
	cov = numpy.load(path_to_cov)
	myTools.plot2D(cov)
	myTools.plot2D(myTools.getCorrelationMatrix(cov))
	myTools.plotCovar( [cov], ['a'] )
	'''

	print corr._meanZ

	'''
	
	dic_CAMB['mulpol_index'] = 0
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,0,:],dic_CAMB,False)
	print numpy.sqrt(dic_CAMB['b'])
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB,0,False)
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB,1,False)
	corr.plot_CAMB(corr._xiMul[:,0,:],dic_CAMB,2,False)
	
	dic_CAMB['mulpol_index'] = 2
	dic_CAMB = corr.fit_CAMB(corr._xiMul[:,2,:],dic_CAMB,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB,0,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB,1,False)
	corr.plot_CAMB(corr._xiMul[:,2,:],dic_CAMB,2,False)
	
	dic_CAMB = corr.fit_CAMB(None,dic_CAMB,False)
	print numpy.sqrt(dic_CAMB['b'])
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
	
	corr.plot_grid(0,True)
	corr.plot_grid(1,True)
	corr.plot_grid(2,True)
	'''	

	corr.send_BAOFIT('',None,False,False, 0, False, True)

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
if (True):
	dic_class['path_to_txt_file_folder'] = dic_simu['path_to_simu'] + 'Box_000/Simu_000/Results/'
	corr = annalyse_BAOFIT_Q.AnnalyseBAOFIT(dic_class, None,None,dic_Q)
	corr.set_values_on_mean_simulation(dic_simu)
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
	corr.plot_chi2_scan(50,50, [0.98,1.02,0.98,1.02], False, False, None)






















