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
import scipy
from scipy import interpolate
import cosmolopy.perturbation as cp
import copy

### Perso lib
import const
import myTools

raw_dic = {
	'z' : 2.25,
	'h_0' : 0.70,
	'omega_matter_0' : 0.27,
	'omega_lambda_0' : 0.73,
	'source'         : 'CAMB'
}

class CAMB:
	"""

	### BAOFIT (Kirkby et al. 2013)
	http://arxiv.org/pdf/1301.3456v1.pdf
	### Hamilton 1992
	http://cdsads.u-strasbg.fr/cgi-bin/nph-iarticle_query?1992ApJ...385L...5H&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf

	"""
	
	def __init__(self,dic=None):

		if (dic is None):
			dic = copy.deepcopy(raw_dic)
		self._z              = dic['z']
		self._h_0            = dic['h_0']
		self._omega_matter_0 = dic['omega_matter_0']
		self._omega_lambda_0 = dic['omega_lambda_0']
		self._source         = dic['source']

		if (self._source=='CAMB'):
			self._xi0 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.0.dat')
			self._xi2 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.2.dat')
			self._xi4 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.4.dat')
		elif (self._source=='CAMB_Mocks'):
			self._xi0 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocks.0.dat')
			self._xi2 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocks.2.dat')
			self._xi4 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocks.4.dat')
		elif (self._source=='CAMB_Mocks_me'):
			path_to_pk = '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Mock_JMLG/Produce_CAMB/DR9LyaMocks_matterpower.dat'
			self._xi0,self._xi2,self._xi4 = self.get_xi_0_2_4_from_pk(path_to_pk,[0,1],True)
		elif (self._source=='CHRISTOPHE'):
			path_to_pk = '/home/gpfs/manip/mnt0607/bao/cmv/Helion/ginit3d_67_0p0_7859_ntpk.txt'
			self._xi0,self._xi2,self._xi4 = self.get_xi_0_2_4_from_pk(path_to_pk,[0,3],False)

		return
	### Window functions
	def window_spherical(self,k,R):
		"""
			http://arxiv.org/pdf/astro-ph/9710252v1.pdf
		"""

		#RR         = numpy.power(3./(4*numpy.pi),1./3.)*R
		RR         = R
		a          = numpy.array(k)
		kk         = k[(k!=0.)]
		a[(k==0.)] = 1.
		a[(k!=0.)] = 3.*( numpy.sin(kk*RR)-kk*RR*numpy.cos(kk*RR) )/numpy.power(kk*RR,3.)

		return a
	def window_gauss(self,k,R):
		"""
			http://arxiv.org/pdf/astro-ph/9710252v1.pdf
		"""

		#RR         = 1./numpy.sqrt(2.*numpy.pi)*R
		return numpy.exp( -0.5*numpy.power(k*RR,2.) )
	def get_xi_0_2_4_from_pk(self,path_to_pk,index=[0,1],comoving=True):

		f = cp.fgrowth(self._z,self._omega_matter_0)
		print '  The growth factor is : f = ', f

		data = numpy.loadtxt(path_to_pk)
		k  = numpy.zeros( data[:,1].size+1 )
		pk = numpy.zeros( data[:,1].size+1 )

		if (comoving) :
			k = data[:,index[0]]
			pk = data[:,index[1]]*f*f
		else:
			k = data[:,index[0]]/self._h_0
			pk = data[:,index[1]]*numpy.power(self._h_0,3.)*f*f

		pk *= numpy.power(self.window_spherical(k,4.5*0.7),2.)

		r2,cric2 = self.xi_from_pk(k,pk)
		index_of_bin_greeter_than_half_max = numpy.arange(r2.size)[ r2>numpy.max(r2)/2. ][0]
		size = r2[1:index_of_bin_greeter_than_half_max].size
		print '  index of bin greeter than half max  =  ', index_of_bin_greeter_than_half_max
		print '  size = ', size
		xi0 = numpy.zeros( shape=(size,2) )
		xi0[:,0] = r2[1:index_of_bin_greeter_than_half_max]
		xi0[:,1] = cric2[1:index_of_bin_greeter_than_half_max]

		xi2, xi4 = self.get_xi2_and_xi4_from_xi0(xi0)

		return xi0, xi2, xi4
	def get_xi2_and_xi4_from_xi0(self,xi0):

		size = xi0[:,0].size

		step = xi0[1,0]-xi0[0,0]
		diff = 3.*numpy.power(xi0[:,0],-3.)*numpy.cumsum( numpy.power(xi0[:,0],2.)*xi0[:,1])*step
		xi2 = numpy.zeros( shape=(size,2) )
		xi2[:,0] = xi0[:,0]
		xi2[:,1] = xi0[:,1]-diff
		xi2[:,1] -= (xi2[0,1]-xi0[0,1])*numpy.power(xi2[0,0]/xi2[:,0],3.)

		step = xi0[1,0]-xi0[0,0]
		diff = 5.*numpy.power(xi0[:,0],-5.)*numpy.cumsum( numpy.power(xi0[:,0],4.)*(xi0[:,1]+xi2[:,1]) )*step
		xi4 = numpy.zeros( shape=(size,2) )
		xi4[:,0] = xi0[:,0]
		xi4[:,1] = xi0[:,1] - diff
		xi4[:,1] -= (xi4[0,1]-xi0[0,1])*numpy.power(xi4[0,0]/xi4[:,0],5.)

		return xi2, xi4
	def get_beta_knowing_bias_QSO(self,bias,z_eff):

		O_m  = self._omega_matter_0 * numpy.power(1+z_eff,3) / (self._omega_matter_0 * numpy.power(1+z_eff,3) + self._omega_lambda_0)
		f = numpy.power(O_m,0.55)
		print '  at z = ', z_eff
		print '  f = ', f
		print '  beta = ', f/bias
		beta = f/bias

		return beta
	def get_bias_knowing_beta_QSO(self,beta,z_eff):

		O_m  = self._omega_matter_0 * numpy.power(1+z_eff,3) / (self._omega_matter_0 * numpy.power(1+z_eff,3) + self._omega_lambda_0)
		f = numpy.power(O_m,0.55)
		print '  at z = ', z_eff, '  f = ', f
		bias = f/beta

		return bias
	def get_bias_effectif_knowing_bias_QSO(self,bias,z_eff):

		beta     = self.get_beta_knowing_bias_QSO(bias,z_eff)
		bias_eff = bias*numpy.sqrt(1.+(2./3.)*beta + (1./5.)*beta*beta )

		return bias_eff
	def get_bias_effectif_knowing_beta_QSO(self,beta,z_eff):

		bias     = self.get_bias_knowing_beta_QSO(beta,z_eff)
		bias_eff = bias*numpy.sqrt(1.+(2./3.)*beta + (1./5.)*beta*beta )

		return bias_eff
	def get_multipol_coef_knowing_bias_QSO(self,bias,z_eff):

		beta     = self.get_beta_knowing_bias_QSO(bias,z_eff)

		c0 = bias*bias*( 1.+(2./3.)*beta + (1./5.)*beta*beta )
		c2 = bias*bias*( (4./3.)*beta + (4./7.)*beta*beta )
		c4 = bias*bias*( (8./35.)*beta*beta )

		return c0, c2, c4
	def get_multipol_coef_knowing_beta_QSO(self,bias,z_eff):

		bias = self.get_bias_knowing_beta_QSO(beta,z_eff)

		c0 = bias*bias*( 1.+(2./3.)*beta + (1./5.)*beta*beta )
		c2 = bias*bias*( (4./3.)*beta + (4./7.)*beta*beta )
		c4 = bias*bias*( (8./35.)*beta*beta )

		return c0, c2, c4
	def get_multipol_coef(self,bias,beta):

		c0 = bias*bias*( 1.+(2./3.)*beta + (1./5.)*beta*beta )
		c2 = bias*bias*( (4./3.)*beta + (4./7.)*beta*beta )
		c4 = bias*bias*( (8./35.)*beta*beta )

		return c0, c2, c4
	def xi_from_pk(self,k,pk):
		"""
		------------------------------------------------------------
			from P(k) to xi(r) for uneven spaced k points
			From Etienne CamLib.py
		------------------------------------------------------------
		"""

		pkInter=scipy.interpolate.InterpolatedUnivariateSpline(k,pk) #,kind='cubic')
		nk=1000000
		kmax=numpy.max(k)
		kmin=numpy.min(k)
		kIn=numpy.linspace(kmin,kmax,nk)
		pkIn=pkInter(kIn)
		kIn[0]=0.
		pkIn[0]=0.
		r=2.*numpy.pi/kmax*numpy.arange(nk)
		pkk=kIn*pkIn
		cric = numpy.append([0.], -numpy.imag(numpy.fft.fft(pkk)/nk)[1:]/r[1:]/2./numpy.pi**2*kmax )

		return r,cric
	def plot_1d(self,x_power=0):

		xxx = self._xi0[:,0]
		yyy = self._xi0[:,1]
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,fmt='o')

		xxx = self._xi2[:,0]
		yyy = self._xi2[:,1]
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,fmt='o')

		xxx = self._xi4[:,0]
		yyy = self._xi4[:,1]
		coef = numpy.power(xxx,x_power)

		plt.errorbar(xxx,coef*yyy,fmt='o')
		if (x_power==0):
			plt.ylabel(r'$ \xi (|s|)$', fontsize=40)
		if (x_power==1):
			plt.ylabel(r'$|s|.\xi (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (x_power==2):
			plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
		plt.xlim([ numpy.amin(xxx)-10., numpy.amax(self._xi0[:,0])+10. ])
		myTools.deal_with_plot(False,False,False)
		plt.show()

		return
	def test_cosmolopy(self):


		import cosmolopy.perturbation as cp
		cosmo = {
			'omega_M_0' : 0.27,
			'omega_lambda_0' : 1.-0.27,
			'omega_b_0' : 0.0463265,
			'omega_n_0' : 0.0,
			'N_nu' : 0,
			'h' : 0.7,
			'n' : 1.0,
			'sigma_8' : 0.794961,
			'baryonic_effects' : True
		}

		z = 2.4
		za = numpy.arange(0.,10.,0.01)
		k = numpy.arange(0.,100.,0.0001)

		f = cp.fgrowth(za, cosmo['omega_M_0'])
		plt.plot(za,f)
		plt.show()

		pk = cp.power_spectrum(k, z, **cosmo)

		plt.plot(k,k*k*pk)
		plt.show()

		print cp.transfer_function_EH(k, **cosmo)
		plt.plot(k,cp.transfer_function_EH(k, **cosmo)[0] )
		plt.show()

		r2,cric2 = self.xi_from_pk(k,pk)
		r2 *= 0.7

		plt.plot(r2,r2*r2*cric2)
		plt.show()

		return
	def save_in_one_file(self, path_to_save):

		numpy.savetxt(path_to_save,zip(self._xi0[:,0], self._xi0[:,1], self._xi2[:,1], self._xi4[:,1]),fmt='%1.20e %1.20e %1.20e %1.20e')
		

		return

"""
'''
camb = CAMB('CAMB_Mocks_me')
camb.plot_1d(0)
camb.plot_1d(1)
camb.plot_1d(2)
'''

dic = {
	'z' : 0.,
	'h_0' : 0.70,
	'omega_matter_0' : 0.27,
	'omega_lambda_0' : 0.73,
	'source'         : 'CAMB_Mocks_me'
}
camb  = [ CAMB(dic) ]

for x_power in [0,1,2]:
	for el in camb:
		
		xxx = el._xi0[:,0]
		yyy = el._xi0[:,1]
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,fmt='o')

		
		xxx = el._xi2[:,0]
		yyy = el._xi2[:,1]
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,fmt='o')

		xxx = el._xi4[:,0]
		yyy = el._xi4[:,1]
		coef = numpy.power(xxx,x_power)
		plt.errorbar(xxx,coef*yyy,fmt='o')
		

	#plt.xlim([ numpy.amin(xxx)-10., 200.+10. ])
	myTools.deal_with_plot(False,False,False)
	plt.show()

camb[0].save_in_one_file('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/CAMB/CAMB_0/camb.txt')
"""











