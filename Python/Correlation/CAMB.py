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


### Perso lib
import const
import myTools

class CAMB:

	##http://arxiv.org/pdf/1301.3456v1.pdf
	
	def __init__(self,source='CAMB'):

		if (source=='CAMB'):
			self._xi0 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.0.dat')
			self._xi2 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.2.dat')
			self._xi4 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocksLCDM.4.dat')
		if (source=='CAMB_Mocks'):
			self._xi0 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocks.0.dat')
			self._xi2 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocks.2.dat')
			self._xi4 = numpy.loadtxt( const.path_to_BAOFIT_model__ + 'DR9LyaMocks.4.dat')
		if (source=='CAMB_Mocks_me'):
			
			path_to_pk = '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Mock_JMLG/Produce_CAMB/DR9LyaMocks_matterpower.dat'
			self._xi0,self._xi2,self._xi4 = self.get_xi_0_2_4_from_pk(path_to_pk,[0,1],True)

		elif (source=='CHRISTOPHE'):

			path_to_pk = '/home/gpfs/manip/mnt0607/bao/cmv/Helion/ginit3d_67_0p0_7859_ntpk.txt'
			self._xi0,self._xi2,self._xi4 = self.get_xi_0_2_4_from_pk(path_to_pk,[0,3],False)

		return
	def get_xi_0_2_4_from_pk(self,path_to_pk,index=[0,1],comoving=True):

		h         = 0.70
		omega_m_0 = 0.27
		f         = cp.fgrowth(2.25,omega_m_0)

		data = numpy.loadtxt(path_to_pk)
		data_me = numpy.zeros( shape=(data[:,1].size+1,2) )

		if (comoving) :
			data_me[1:,0] = data[:,index[0]]
			data_me[1:,1] = data[:,index[1]]*f*f
		else:
			data_me[1:,0] = data[:,index[0]]/h
			data_me[1:,1] = data[:,index[1]]*numpy.power(h,3.)*f*f

		r2,cric2 = self.xi_from_pk(data_me[:,0],data_me[:,1])
		index_of_bin_greeter_than_half_max = numpy.arange(r2.size)[ r2>numpy.max(r2)/2. ][0]
		size = r2[1:index_of_bin_greeter_than_half_max].size
		print '  index of bin greeter than half max  =  ', index_of_bin_greeter_than_half_max
		print '  size = ', size
		xi0 = numpy.zeros( shape=(size,2) )
		xi0[:,0] = r2[1:index_of_bin_greeter_than_half_max]
		xi0[:,1] = cric2[1:index_of_bin_greeter_than_half_max]

		print '  starting integration with r = ', xi0[0,0]
		step = xi0[1,0]-xi0[0,0]
		print '  step = ', step
		diff = 3.*numpy.power(xi0[:,0],-3.)*numpy.cumsum( numpy.power(xi0[:,0],2.)*xi0[:,1])*step
		xi2 = numpy.zeros( shape=(size,2) )
		xi2[:,0] = xi0[:,0]
		xi2[:,1] = xi0[:,1]-diff

		#plt.plot(xi2[:,0], -xi2[:,1] )
		#plt.plot(xi2[:,0], -numpy.power(xi0[0,0]/xi2[:,0],3.)*(xi2[:,1]-xi0[:,1]) )
		#plt.show()
		#xi2[:,1] += numpy.power(xi2[0,0]/xi2[:,0],3.)*(xi2[:,1]-xi0[:,1])

		step = xi0[1,0]-xi0[0,0]
		diff = 5.*numpy.power(xi0[:,0],-5.)*numpy.cumsum( numpy.power(xi0[:,0],4.)*(xi0[:,1]+xi2[:,1]) )*step
		xi4 = numpy.zeros( shape=(size,2) )
		xi4[:,0] = xi0[:,0]
		xi4[:,1] = xi0[:,1] - diff

		#plt.plot(xi4[:,0], abs( xi4[:,1] ) )
		#plt.plot(xi4[:,0], abs( numpy.power(xi0[0,0]/xi4[:,0],5.)*(xi4[:,1]-xi0[:,1]) ) )
		#plt.show()
		#xi4[:,1] += numpy.power(xi4[0,0]/xi4[:,0],5.)*(xi4[:,1]-xi0[:,1])

		return xi0, xi2, xi4
	def xi_from_pk(self,k,pk):
		#------------------------------------------------------------
		#	from P(k) to xi(r) for uneven spaced k points
		#	From Etienne CamLib.py
		#------------------------------------------------------------
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
	'''
	def get_xi2(self):

		### http://cdsads.u-strasbg.fr/cgi-bin/nph-iarticle_query?1992ApJ...385L...5H&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf

		rr  = self._xi0[:,0]
		xi = self._xi0[:,1]
		interp0 = interpolate.interp1d(rr,xi,bounds_error=False,fill_value=1)
		interp  = interpolate.interp1d(rr,rr*rr*xi,bounds_error=False,fill_value=1)

		def xi_bar(r):
			y = scipy.integrate.quad(interp,0.,r)[0]
			return y

		self._xi2[:,0] = self._xi0[:,0]

		a = 0.
		for i in numpy.arange(1,100):
			self._xi2[i,1] = scipy.integrate.quad(interp,0.,rr[i+1])[0]

		self._xi2[1:,1] = self._xi0[1:,1] - 3.*numpy.power(self._xi0[1:,0],-3.)*self._xi2[1:,1]


		plt.plot(rr,self._xi2[:,1],marker='o')
		plt.plot(rr,self._xi0[:,1],marker='o')
		plt.show()
		plt.plot(rr,rr*self._xi2[:,1],marker='o')
		plt.plot(rr,rr*self._xi0[:,1],marker='o')
		plt.show()
		plt.plot(rr,rr*rr*self._xi2[:,1],marker='o')
		plt.plot(rr,rr*rr*self._xi0[:,1],marker='o')
		plt.show()

		return
	'''
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

'''
camb = CAMB('CAMB_Mocks_me')
camb.plot_1d(0)
camb.plot_1d(1)
camb.plot_1d(2)
'''

'''
camb = [ CAMB('CAMB_Mocks_me') ]

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
		

	plt.xlim([ numpy.amin(xxx)-10., 200.+10. ])
	myTools.deal_with_plot(False,False,False)
	plt.show()
'''

#camb[0].save_in_one_file('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/CAMB/CAMB_2_25/camb.txt')












