# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

### Python lib
import numpy
import scipy
import copy
import matplotlib.pyplot as plt
import cosmolopy.perturbation as cp

### Perso lib
import CAMB
import myTools


raw_dic = {
	'sigma_gauss'     : 1.,
	'z'               : 2.25,
	'a'               : 0.,
	'gamma'           : 1.6,
	'cut_QSO'         : 3.36,
	'h_0'             : 0.70,
	'omega_matter_0'  : 0.27,
	'omega_lambda_0'  : 0.73,
	'correlation'     : 'q_f',
	'source_for_CAMB' : 'CAMB',
	'mock_type'       : 'JMC', ## Andreu
}

class TRANSFORM_CAMB_BIAS:

	"""

	Apply the transformation of the linear correlation-function by the bias selection function
		http://arxiv.org/pdf/1205.2018v1.pdf (page 24++)

	/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1575/Pk.cc

	"""
	
	def __init__(self, dic=None):
		
		if (dic is None): dic = copy.deepcopy(raw_dic)
		self._sigma_gauss     = dic['sigma_gauss']
		self._z               = dic['z']
		self._a               = dic['a']
		self._gamma           = dic['gamma']
		self._cut_QSO         = dic['cut_QSO']
		self._h_0             = dic['h_0']
		self._omega_matter_0  = dic['omega_matter_0']
		self._omega_lambda_0  = dic['omega_lambda_0']
		self._source_for_CAMB = dic['source_for_CAMB']

		### Get correlation function
		dic_CAMB = {
			'z' : 0.,
			'h_0' : self._h_0,
			'omega_matter_0' : self._omega_matter_0,
			'omega_lambda_0' : self._omega_lambda_0,
			'source'         : self._source_for_CAMB
		}
		self._CAMB = CAMB.CAMB(dic_CAMB)

		### Mock type
		self.get_parameter_a     = None 
		self.get_parameter_gamma = None 
		if (dic['mock_type'] == 'JMC'):
			self.get_parameter_a     =  self.get_parameter_a_JMC
			self.get_parameter_gamma =  self.get_parameter_gamma_JMC
		elif (dic['mock_type'] == 'Andreu'):
			self.get_parameter_a     =  self.get_parameter_a_Andreu
			self.get_parameter_gamma =  self.get_parameter_gamma_Andreu
		self._a     = self.get_parameter_a(self._z)
		self._gamma = self.get_parameter_gamma(self._z)

		self._correlation = dic['correlation']
		self._F1 = None
		self._F2 = None
		self.lim_min_mean_F1 = None
		self.lim_min_mean_F2 = None

		### Auto-correlation
		if (self._correlation=='i_i'):
			self._F1 = self.f_id
			self._F2 = self.f_id
			self.lim_min_mean_F1 = -numpy.inf
			self.lim_min_mean_F2 = -numpy.inf
		### Auto-correlation
		if (self._correlation=='f_f'):
			self._F1 = self.f_Gunn_Peterson
			self._F2 = self.f_Gunn_Peterson
			self.lim_min_mean_F1 = -numpy.inf
			self.lim_min_mean_F2 = -numpy.inf
		### Cross-correlation
		elif (self._correlation=='q_f'):
			self._F1 = self.f_Gunn_Peterson
			self._F2 = self.f_QSO_Selection
			self.lim_min_mean_F1 = -numpy.inf
			self.lim_min_mean_F2 = self._cut_QSO
		### Auto-correlation QSO
		elif (self._correlation=='q_q'):
			self._F1 = self.f_QSO_Selection
			self._F2 = self.f_QSO_Selection
			self.lim_min_mean_F1 = self._cut_QSO
			self.lim_min_mean_F2 = self._cut_QSO

		return
	### Get the parameters functions
	def get_parameter_a_JMC(self, z):
		return numpy.power(10., -7.04712 + 2.88486*z -0.300941*z*z)
	def get_parameter_gamma_JMC(self, z):
		return self._gamma + z*0.
	def get_parameter_a_Andreu(self, z):
		z_coord = numpy.array([1.96, 2.44, 2.91, 3.39])
		a_coord = numpy.array([0.065, 0.141, 0.275, 0.487])
		return numpy.interp(z,z_coord,a_coord)
	def get_parameter_gamma_Andreu(self, z):
		z_coord     = numpy.array([1.96, 2.44, 2.91, 3.39])
		gamma_coord = numpy.array([1.70, 1.53, 1.38, 1.24])
		return numpy.interp(z,z_coord,gamma_coord)
	### Bias functions
	def f_id(self,x):
		return 1.+x
	def f_Gunn_Peterson(self,x):
		return numpy.exp(-self._a*numpy.exp(self._gamma*x) )
	def f_QSO_Selection(self,x):
		x = numpy.array(x)
		cut = x<self._cut_QSO
		x[cut] = 0.
		cut = x>=self._cut_QSO
		x[cut] = 1.
		return x
	def get_mean_function(self,function_index):
		"""
			Get the mean of the given transforming functions

		"""
		
		if (function_index==0):
			func = self._F1
			lim_min = self.lim_min_mean_F1
		elif (function_index==1):
			func = self._F2
			lim_min = self.lim_min_mean_F2
					
		def func2(x):
			return func(x)*numpy.exp(-0.5*numpy.power(x/self._sigma_gauss,2.))/(numpy.sqrt(2.*numpy.pi)*self._sigma_gauss)
		
		return scipy.integrate.quad(func2,lim_min,numpy.inf)[0]
	def get_mean_function_vs_redshift(self,function_index,z_min, z_max, z_step=0.1):
		"""
			Get the mean of the given transforming functions versus the redshift

		"""
		
		z    = numpy.arange(z_min, z_max, z_step)
		mean = numpy.zeros(z.size)
		
		for i in numpy.arange(mean.size):
			d_growth = cp.fgrowth(z[i],self._omega_matter_0)
			self._a     = self.get_parameter_a(z[i])
			self._gamma = 1.6*d_growth #self.get_parameter_gamma(z[i])*d_growth
			mean[i] = self.get_mean_function(function_index)
			
		self._a     = self.get_parameter_a(self._z)
		self._gamma = self.get_parameter_gamma(self._z)
		
		return numpy.array( zip(z,mean) )
	def get_variance_function_vs_redshift(self,function_index,z_min, z_max, z_step=0.1):
		"""
			Get the mean of the given transforming functions versus the redshift

		"""
		
		z    = numpy.arange(z_min, z_max, z_step)
		mean = numpy.zeros(z.size)
		variance = numpy.zeros(z.size)
		
		for i in numpy.arange(mean.size):
			d_growth = cp.fgrowth(z[i],self._omega_matter_0)
			self._a     = self.get_parameter_a(z[i])
			self._gamma = 1.6*d_growth
			mean[i]     = self.get_mean_function(function_index)

			self._a     *= 2.
			variance[i] = self.get_mean_function(function_index)
			
		self._a     = self.get_parameter_a(self._z)
		self._gamma = self.get_parameter_gamma(self._z)
		
		var_delta = variance/(mean*mean) - 1.

		return numpy.array( zip(z,mean,variance,var_delta) )
	def get_xi0_transformed(self, xi_r, xi_z):
		"""
			Get the transformed correlation function:
			xi_r: coordinate of the points where to calculate
				the correlation function.

		"""
		
		verbose = True

		d_growth    = cp.fgrowth(xi_z,self._omega_matter_0)
		xi0         = numpy.interp(xi_r,self._CAMB._xi0[:,0],self._CAMB._xi0[:,1])*d_growth*d_growth
		self._a     = self.get_parameter_a(xi_z)
		self._gamma = self.get_parameter_gamma(xi_z)
		mean_F1     = self.get_mean_function(0)
		mean_F2     = self.get_mean_function(1)

		if (verbose):
			print '\n'
			print '  Producing the transformation at z = ', xi_z
			print '  D_growth_factor = ', d_growth
			if ( self._correlation=='q_f' or self._correlation=='f_f' ):
				print '  a               = ', self._a
				print '  gamma           = ', self._gamma
			if (self._correlation=='q_f' or self._correlation=='q_q'):
				print '  cut_QSO         = ', self._cut_QSO
			print '  < F1 >          = ', mean_F1
			print '  < F2 >          = ', mean_F2
			print '\n'

		reduced_xi = xi0*xi0/numpy.power(self._sigma_gauss,4.)
		if ( reduced_xi[ reduced_xi>1. ].size  ):
			print '  ERROR:: reduced_xi>1.'
			return None

		def get_transformation(xi):

			sqrt_xi = numpy.sqrt(1.-xi*xi/numpy.power(self._sigma_gauss,4.) )
	
			### Function to integrate
			def F1_times_F2(y,x):
				t  = xi*x/numpy.power(self._sigma_gauss,2.)+sqrt_xi*y
				f1 = self._F1(x)*numpy.exp(-0.5*numpy.power(x/self._sigma_gauss,2.))
				f2 = self._F2(t)*numpy.exp(-0.5*numpy.power(y/self._sigma_gauss,2.))
				return f1*f2
	
			### Integral range
			if (self._correlation=='i_i' or self._correlation=='f_f'):
				def lim_min_y(x):
					return -numpy.inf+x*0.
			elif (self._correlation=='q_f' or self._correlation=='q_q'):
				def lim_min_y(x):
					return (self._cut_QSO-xi*x/numpy.power(self._sigma_gauss,2.))/sqrt_xi
	
			def lim_max_y(x):
				return numpy.inf+x*0.
					
			return scipy.integrate.dblquad(F1_times_F2, self.lim_min_mean_F1,numpy.inf, lim_min_y, lim_max_y)[0]
	
		return numpy.array([ get_transformation(el) for el in xi0 ])/(mean_F1*mean_F2*2.*numpy.pi*self._sigma_gauss*self._sigma_gauss) - 1.
	def get_bias_vs_cut_QSO(self, cut_QSO_min, cut_QSO_max, cut_QSO_step=0.1):
		
		cut_QSO_before = copy.deepcopy(self._cut_QSO)
		cut_QSO = numpy.arange(cut_QSO_min, cut_QSO_max, cut_QSO_step)
		xi_r = numpy.arange(55., 65., 1.)
			
		bias = numpy.zeros(cut_QSO.size)
		for i in numpy.arange(cut_QSO.size):
			print i
			self._cut_QSO = cut_QSO[i]
			transformed_xi_0 = self.get_xi0_transformed(xi_r)
			bias[i] = numpy.interp(60.,xi_r,transformed_xi_0) / numpy.interp(60.,self._CAMB._xi0[:,0],self._CAMB._xi0[:,1])

		self._cut_QSO = cut_QSO_before

		return numpy.array( zip(cut_QSO,bias) )
	def plot_function(self,function_index):
		"""
		Plot the F1 and F2 function

		"""

		
		if (function_index==0):
			func = self._F1
		elif (function_index==1):
			func = self._F2

		### Plot selection function
		def func2(x):
			return numpy.exp(-0.5*numpy.power(x/self._sigma_gauss,2.))/(numpy.sqrt(2.*numpy.pi)*self._sigma_gauss)
		xxx = numpy.arange(-10.,10.,0.01)
		yyy = func(xxx)
		plt.plot(xxx,yyy, markersize=10,linewidth=2,label=r'$Selection$')
		yyy = func2(xxx)
		plt.plot(xxx,yyy, markersize=10,linewidth=2,label=r'$Gauss$')
		plt.xlabel(r'$\delta$')
		plt.ylabel(r'$Proba$')
		myTools.deal_with_plot(False,False,True)
		plt.show()

		### Plot results
		meanF = self.get_mean_function(function_index)
		def func2(x):
			return func(x)*numpy.exp(-0.5*numpy.power(x/self._sigma_gauss,2.))/(numpy.sqrt(2.*numpy.pi)*self._sigma_gauss)
		
		xxx = numpy.arange(-10.,10.,0.01)
		yyy = func2(xxx)
		def func2(x):
			return numpy.exp(-0.5*numpy.power(x/self._sigma_gauss,2.))/(numpy.sqrt(2.*numpy.pi)*self._sigma_gauss)
		meanF =  scipy.integrate.quad(func2,-numpy.inf,numpy.inf)[0]

		xxx2 = numpy.arange(-10.,10.,0.01)
		yyy2 = func2(xxx)

		plt.plot(xxx2,yyy2, markersize=10,linewidth=2)
		plt.plot(xxx,yyy, markersize=10,linewidth=2)
		plt.xlabel(r'$\delta$')
		plt.ylabel(r'$F(\delta) \cdot Gauss$')
		myTools.deal_with_plot(False,False,False)
		plt.show()

		return
	def get_PDF(self, nbBin_T, min_z, max_z, nbBin_z):

		"""
			Get the PDF for the given Gunn-Peterson model
		"""

		
		def gaussian(x):
			return numpy.exp(-0.5*numpy.power(x/self._sigma_gauss,2.))/(numpy.sqrt(2.*numpy.pi)*self._sigma_gauss)
		def from_T_to_delta(x):
			to_return = numpy.array(x)
			to_return[ (x==0.) ] = numpy.inf + x*0.
			to_return[ (x!=0.) ] = -(1./self._a)*numpy.log(x[ (x!=0.) ])
			to_return[ (to_return!=0.) ] = (1./self._gamma)*numpy.log(to_return[ (to_return!=0.) ])
			to_return[ (to_return==0.) ] = -numpy.inf + x*0.
			return to_return

		extent=[min_z, max_z, 0., 1.]

		xxx = numpy.append( numpy.arange(0.,1.,1./nbBin_T), [1.])
		yyy = numpy.append( numpy.arange(min_z, max_z, (max_z-min_z)/nbBin_z ), [max_z])
		yyy = numpy.array([ yyy[i]+(yyy[i+1]-yyy[i])/2. for i in range(0,yyy.size-1) ])
		zzz = numpy.zeros( shape=(nbBin_T,nbBin_z) )

		for i in numpy.arange(0,nbBin_z):
			self._a     = self.get_parameter_a(yyy[i])
			self._gamma = self.get_parameter_gamma(yyy[i])
			zzz[:,i]    = numpy.array([  scipy.integrate.quad( gaussian, from_T_to_delta(xxx[j+1]), from_T_to_delta(xxx[j]) )[0] for j in numpy.arange(0,nbBin_T) ])

			plt.plot(xxx[:-1], zzz[:,i])
		plt.show()

		xxx  = numpy.array([ xxx[i]+(xxx[i+1]-xxx[i])/2. for i in range(0,xxx.size-1) ])

		plt.imshow(zzz, origin='lower',extent=extent, interpolation='None')
		plt.show()
		plt.imshow(numpy.log10(zzz), origin='lower',extent=extent, interpolation='None')
		plt.show()





"""
### Test redshift evolution
sig8_grf = 0.794961095
lambda_Lya = 1215.67
### for v1575
sigma_LR = 1.96769
sigma_HR = 8.606
factor = 1. #15
### for v1588
#sigma_LR = 1.96769
#sigma_HR = 6.578
#factor = 1. #1.02

sigma = numpy.sqrt( numpy.power(sigma_LR/sig8_grf,2.) + numpy.power(sigma_HR,2.) )

###
dic = {
	'sigma_gauss'     : sigma,
	'z'               : -100000.,
	'a'               : numpy.inf,
	'gamma'           : 1.6,
	'cut_QSO'         : numpy.inf,
	'h_0'             : 0.70,
	'omega_matter_0'  : 0.27,
	'omega_lambda_0'  : 0.73,
	'correlation'     : 'f_f',
	'source_for_CAMB' : 'CAMB_Mocks_me',
	'mock_type'       : 'JMC',
}
a = TRANSFORM_CAMB_BIAS(dic)

#data = a.get_mean_function_vs_redshift2()
#plt.plot( data[:,0], data[:,1]/factor , markersize=10,linewidth=2, label=r'$calcul \, analytique$' )

#a.plot_function(0)
#a.get_PDF(1000, 0., 5., 100)

shift = 3547.5
data = a.get_mean_function_vs_redshift(0, 1.7, 4.) #7234./1215.67-1.)
plt.plot( data[:,0], data[:,1]/factor , markersize=10,linewidth=2, label=r'$Prediction$' )

data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/v1575_hDeltaVsLambdaObs_LYA_JMC.txt')
plt.errorbar(data[:,1][ data[:,2]!=0. ]/1215.67-1., data[:,2][ data[:,2]!=0. ], label=r'$Simulation$', color='red', markersize=10,linewidth=2)

plt.xlabel(r'$z$')
plt.ylabel(r'$\overline{F}$')
myTools.deal_with_plot(False,False,True)
plt.show()



### get variance
factor = 1. #4.
data = a.get_variance_function_vs_redshift(0, 1.7, 4.)
plt.plot( data[:,0], data[:,3]/factor , markersize=10,linewidth=2, label=r'$Prediction$' )

#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575/Box_000/Simu_000/Run/sigma2LSS_LYA_0_0.txt')
#plt.errorbar(data[:,0], data[:,1], yerr=data[:,2] ,fmt='o', markersize=10,linewidth=2, label=r'$Simulation$' )

data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/v1575_h_sigma_LSS_VsRedshift_LYA_JMC.txt')
plt.errorbar(data[:,1], data[:,2], markersize=10,linewidth=2, label=r'$Simulation$' )

plt.xlabel(r'$z$')
plt.ylabel(r'$\sigma_{L.S.S.}^{2}$')
myTools.deal_with_plot(False,False,True)
plt.show()




### Get a evolution
z = numpy.arange(0.,100.,0.1)
aa = a.get_parameter_a(z)

plt.errorbar(z,aa)
myTools.deal_with_plot(False,False,True)
plt.show()
"""



### Test for q_q
"""
dic_CAMB_corr = {
	'z' : 0.,
	'h_0' : 0.70,
	'omega_matter_0' : 0.27,
	'omega_lambda_0' : 0.73,
	'source'         : 'CAMB_Mocks_me'
}
CAMB2 = CAMB.CAMB(dic_CAMB_corr)
z_ref = 2.50
d_growth = cp.fgrowth(z_ref,0.27)
QSO_bias=3.6
mu_grf        = 1.96769
QSO_bias_corr = 0.82
sig8_bias     = 0.8
sig8_grf      = 0.794961095
b_eff         = CAMB2.get_bias_effectif_knowing_bias_QSO(QSO_bias,z_ref)
nu = mu_grf * d_growth * b_eff * QSO_bias_corr * (sig8_bias/sig8_grf)
print '  nu = ', nu
print

z_ref = 0.
dic = {
	'sigma_gauss'     : 1.96769,
	'z'               : 0.,
	'a'               : numpy.inf,
	'gamma'           : numpy.inf,
	'cut_QSO'         : nu*1.96769,
	'h_0'             : 0.70,
	'omega_matter_0'  : 0.27,
	'omega_lambda_0'  : 0.73,
	'correlation'     : 'q_q',
	'source_for_CAMB' : 'CAMB_Mocks_me',
	'mock_type'       : 'JMC',
}
a = TRANSFORM_CAMB_BIAS(dic)

#a.plot_function(0)

r_max = 200.
xi_r = numpy.arange(0.,r_max,1.)
transformed_xi_0 = a.get_xi0_transformed(xi_r,z_ref)
xxx_CAMB = a._CAMB._xi0[:,0]
yyy_CAMB = a._CAMB._xi0[:,1]*d_growth*d_growth

z_ref = 2.5
### Bias
#bias = numpy.interp(60.,xi_r,transformed_xi_0) / numpy.interp(60.,xxx_CAMB,yyy_CAMB)
bias = numpy.mean(transformed_xi_0[ numpy.logical_and(xi_r>180.,xi_r<200.) ]) / numpy.mean(yyy_CAMB[ numpy.logical_and(xxx_CAMB>180.,xxx_CAMB<200.) ])
print bias
print '  bias  = ', numpy.sqrt(bias)
b_eff = CAMB2.get_bias_effectif_knowing_bias_QSO(numpy.sqrt(bias),z_ref)
print '  b_eff = ', b_eff

###
for i in range(0,3):

	coef = numpy.power(xxx_CAMB,i)
	plt.plot(xxx_CAMB,coef*yyy_CAMB)
	coef = numpy.power(xi_r,i)
	plt.errorbar(xi_r, coef*transformed_xi_0/bias, fmt='o')

	plt.xlim([0.,r_max])
	plt.show()

### Look bias evolution
f   = cp.fgrowth(z_ref,a._omega_matter_0)
xi0 = numpy.interp(xi_r,a._CAMB._xi0[:,0],a._CAMB._xi0[:,1])*d_growth*d_growth
yyy = (transformed_xi_0/xi0 - bias)/bias*100.
plt.errorbar(xi_r, yyy, fmt='o')
plt.xlim([0.,r_max])
plt.show()
"""



### Test for f_f
"""
#sig8_grf = 0.794961095
#coef     = 0.951



dic = {
	'sigma_gauss'     : 1.96769,  #*d_growth
	'z'               : 0.,
	'a'               : numpy.inf,
	'gamma'           : 1.6,
	'cut_QSO'         : numpy.inf,
	'h_0'             : 0.70,
	'omega_matter_0'  : 0.27,
	'omega_lambda_0'  : 0.73,
	'correlation'     : 'f_f',
	'source_for_CAMB' : 'CAMB_Mocks_me',
	'mock_type'       : 'JMC',
}
a = TRANSFORM_CAMB_BIAS(dic)

r_max = 200.
xi_r = numpy.arange(1.,r_max,10.)
transformed_xi_0 = a.get_xi0_transformed(xi_r,0.)

d_growth = cp.fgrowth(z_ref,0.27)
z_ref    = 2.2920604868
xxx_CAMB = a._CAMB._xi0[:,0]
yyy_CAMB = a._CAMB._xi0[:,1]*d_growth*d_growth

### Bias
bias = numpy.interp(190.,xi_r,transformed_xi_0) / numpy.interp(190.,xxx_CAMB,yyy_CAMB)
print bias
print numpy.sqrt(bias)

###
for i in range(0,3):

	coef = numpy.power(xxx_CAMB,i)
	plt.plot(xxx_CAMB,coef*yyy_CAMB)
	coef = numpy.power(xi_r,i)
	plt.errorbar(xi_r, coef*transformed_xi_0/bias, fmt='o')

	plt.xlim([0.,200.])
	myTools.deal_with_plot(False,False,True)
	plt.show()

### Look bias evolution
d_growth   = cp.fgrowth(z_ref,a._omega_matter_0)
xi0 = numpy.interp(xi_r,xxx_CAMB,yyy_CAMB)
#yyy = (transformed_xi_0/xi0 - bias)/bias*100.
yyy = transformed_xi_0/xi0
plt.errorbar(xi_r, yyy, fmt='o')
plt.xlim([0.,200.])
myTools.deal_with_plot(False,False,True)
plt.show()
"""


### Test for q_f
"""
z_ref = 2.50
d_growth = cp.fgrowth(z_ref,0.27)
dic = {
	'sigma_gauss'     : 1.96769,
	'z'               : z_ref,
	'a'               : numpy.inf,
	'gamma'           : 1.6,
	'cut_QSO'         : 2.3793413839/d_growth,
	'h_0'             : 0.70,
	'omega_matter_0'  : 0.27,
	'omega_lambda_0'  : 0.73,
	'correlation'     : 'q_f',
	'source_for_CAMB' : 'CAMB_Mocks_me',
	'mock_type'       : 'JMC',
}

a = TRANSFORM_CAMB_BIAS(dic)
a.plot_function(0)
a.plot_function(1)


r_max = 200.
xi_r = numpy.arange(1.,r_max,1.)
transformed_xi_0 = a.get_xi0_transformed(xi_r,z_ref)
xxx_CAMB = a._CAMB._xi0[:,0]
yyy_CAMB = a._CAMB._xi0[:,1]*d_growth*d_growth

### Bias
bias = numpy.interp(60.,xi_r,transformed_xi_0) / numpy.interp(60.,xxx_CAMB,yyy_CAMB)
print bias

###
for i in range(0,3):

	coef = numpy.power(xxx_CAMB,i)
	plt.plot(xxx_CAMB,coef*yyy_CAMB*bias)
	coef = numpy.power(xi_r,i)
	plt.errorbar(xi_r, coef*transformed_xi_0, fmt='o')

	plt.xlim([0.,200.])
	myTools.deal_with_plot(False,False,True)
	plt.show()

### Look bias evolution
d_growth   = cp.fgrowth(z_ref,a._omega_matter_0)
xi0 = numpy.interp(xi_r,xxx_CAMB,yyy_CAMB)
#yyy = (transformed_xi_0/xi0 - bias)/bias*100.
yyy = transformed_xi_0/xi0
plt.errorbar(xi_r, yyy, fmt='o')
plt.xlim([0.,200.])
myTools.deal_with_plot(False,False,True)
plt.show()
"""


"""
const double alpha = .00832;  // F = exp (-alpha*tau),  
const double beta = 0.8;  // tau = [exp(beta drho/rho)]^2, physically beta should be 0.8

beta = 1.6;
alpha = pow(10, -7.04712 + 2.88486*z -0.300941*z*z);

"""




