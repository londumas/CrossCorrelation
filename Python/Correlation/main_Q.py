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
	'name' : 'Mean simu'
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
	'path_to_simu' : '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/',
	'nb_box' : 10,
	'nb_simu' : 10,
	'projected' : False,
	'with_metals_templates' : False,
	'raw' : False
}

'''
i=0
j=0
raw = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_00'+str(i)+'/Simu_00'+str(j)+'/'
dic_class['path_to_txt_file_folder'] = raw+'Results/'
dic_Q['path_to_cat']                 = raw+'Data/QSO_withRSD.fits'
corr = correlation_3D_Q.Correlation3DQ(dic_class,dic_Q)
aMB = annalyse_many_BAOFIT.AnnalyseManyBAOFIT(dic_class, None,dic_Q,dic_simu)
#aMB.plot_histo_residuals()
aMB.print_results(0)
aMB.print_results(1)
aMB.print_results(8)
aMB.print_results(9)
'''

### Fit all

### For data
dic_class['name'] = 'Mean \, simu'
dic_Q['load_from_txt'] = False

i = 0
j = 0
raw = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_00'+str(i)+'/Simu_00'+str(j)+'/'
dic_class['path_to_txt_file_folder'] = raw+'Results/'
dic_Q['path_to_cat']                 = raw+'Data/QSO_withRSD.fits'
corr = correlation_3D_Q.Correlation3DQ(dic_class,dic_Q)
#corr.save_list_realisation_simulation(dic_class, dic_Q, dic_simu)
corr.set_values_on_mean_simulation(dic_simu)
corr.send_BAOFIT('',None,False,False, 0, False, True)



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
corr.plot_chi2_scan(100,100, [0.5,1.5,0.5,1.5], False, False)





arg_fit = {}
arg_fit['a'] = 125.08255109750635
arg_fit['b'] = -1.8611671002360046
arg_fit['c'] = 0.006109139533639488
arg_fit['A'] = 0.003731291955288748
arg_fit['sigma'] = 9.870124119190448
arg_fit['mean']  = 106.73882608124278
arg_fit['fix_A'] = False
arg_fit['fix_sigma']  = False
arg_fit['fix_mean']   = False
arg_fit['background'] = 'inv'
arg_fit['cov']   = 'full'
arg_fit['x_min'] = 60.
arg_fit['x_max'] = 160.

listMultipol = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Results/xi_QSO_QSO_result_list_Multipol.npy')
cov0 = numpy.cov( listMultipol[:,0,:] )
cov2 = numpy.cov( listMultipol[:,2,:] )
cov0 = cov0/100.
a,b,c,d,e,f,g =  myTools.fit_BAO(corr._xiMul[:,0,0],corr._xiMul[:,0,1],cov0,arg_fit)
plt.errorbar(corr._xiMul[:,0,0],corr._xiMul[:,0,1]*corr._xiMul[:,0,0]*corr._xiMul[:,0,0],yerr=numpy.sqrt(numpy.diag(cov0))*corr._xiMul[:,0,0]*corr._xiMul[:,0,0],fmt='o')
plt.errorbar(c[0],c[1]*c[0]*c[0])
plt.show()

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


'''
arg_fit = {}
arg_fit['a'] = 125.08255109750635
arg_fit['b'] = -1.8611671002360046
arg_fit['c'] = 0.006109139533639488
arg_fit['A'] = 0.003731291955288748
arg_fit['sigma'] = 9.870124119190448
arg_fit['mean']  = 106.73882608124278
arg_fit['fix_A'] = False
arg_fit['fix_sigma']  = False
arg_fit['fix_mean']   = False
arg_fit['background'] = 'inv'
arg_fit['cov']   = 'full'
arg_fit['x_min'] = 60.
arg_fit['x_max'] = 160.


listMultipol = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Results/xi_QSO_QSO_result_list_Multipol.npy')
cov0 = numpy.cov( listMultipol[:,0,:] )
cov2 = numpy.cov( listMultipol[:,2,:] )
l_corrL  = []
l_b      = []
l_chi2   = []
l_beta   = []
l_chi2_2 = []

nb = 10
mean = numpy.zeros(nb)

for i in range(0,1):
	for j in range(0,1):
		raw = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_00' + str(i) + '/Simu_00' + str(j) +'/'
		dic_class['path_to_txt_file_folder'] = raw+'Results/'
		dic_Q['path_to_cat']                 = raw+'Data/QSO_withRSD.fits'
		dic_Q['load_from_txt']               = False
		corr = Correlation3DQ(dic_class,dic_Q)
		corr.set_error_on_covar_matrix('simulation',dic_simu)
		l_corrL += [corr]
		#
		#dic_CAMB['mulpol_index'] = 0
		#dic_CAMB = corr.fit_CAMB(corr._xiMul[:,0,:],dic_CAMB,False)
		#l_b    += [ numpy.sqrt(dic_CAMB['b']) ]
		#l_chi2 += [dic_CAMB['chi^{2}']]
		#
		#dic_CAMB['mulpol_index'] = 2
		#dic_CAMB = corr.fit_CAMB(corr._xiMul[:,2,:],dic_CAMB,False)
		#l_beta    += [ numpy.sqrt(dic_CAMB['b']) ]
		#l_chi2_2  += [dic_CAMB['chi^{2}']]

		a,b,c,d,e,f,g =  myTools.fit_BAO(corr._xiMul[:,0,0],corr._xiMul[:,0,1],cov0,arg_fit)
		#plt.errorbar(corr._xiMul[:,0,0],corr._xiMul[:,0,1]*corr._xiMul[:,0,0]*corr._xiMul[:,0,0],yerr=corr._xiMul[:,0,2]*corr._xiMul[:,0,0]*corr._xiMul[:,0,0],fmt='o')
		#plt.errorbar(c[0],c[1]*c[0]*c[0])
		#plt.show()
		mean[i*10+j] = a['mean']

plt.hist(mean,bins=10)
plt.show()
'''






'''
#
nb = len(l_corrL)
plt.errorbar( numpy.arange(0,nb),l_b,marker='o' )
plt.errorbar( numpy.arange(0,nb),numpy.ones(nb)*numpy.mean(l_b))
plt.show()
plt.errorbar( numpy.arange(0,nb),l_chi2,marker='o' )
plt.errorbar( numpy.arange(0,nb),numpy.ones(nb)*numpy.mean(l_chi2))
plt.show()
plt.hist(l_b,bins=10)
plt.show()
plt.hist(l_chi2,bins=10)
plt.show()
#
nb = len(corrL)
plt.errorbar( numpy.arange(0,nb),l_beta,marker='o' )
plt.errorbar( numpy.arange(0,nb),numpy.ones(nb)*numpy.mean(l_beta))
plt.show()
plt.errorbar( numpy.arange(0,nb),l_chi2_2,marker='o' )
plt.errorbar( numpy.arange(0,nb),numpy.ones(nb)*numpy.mean(l_chi2_2))
plt.show()
plt.hist(l_beta,bins=10)
plt.show()
plt.hist(l_chi2_2,bins=10)
plt.show()

print numpy.mean(l_b)
print numpy.mean(l_chi2)
print numpy.mean(l_beta)
print numpy.mean(l_chi2_2)
'''



