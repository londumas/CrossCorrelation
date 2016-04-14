# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_1D.py
#


### Python lib
import cosmolopy
import numpy
import matplotlib.pyplot as plt
import cosmolopy.perturbation as cp
from scipy import interpolate
import scipy

import myTools
import const



xxx = numpy.arange(-100.,100.,0.01)
yyy = numpy.sinc(xxx)
plt.plot( xxx, yyy )
plt.show()

def xi_from_pk(k,pk):
	import scipy
	#------------------------------------------------------------
	#	from P(k) to xi(r) for uneven spaced k points
	#	From Etienne CamLib.py
	#------------------------------------------------------------
	pkInter=scipy.interpolate.InterpolatedUnivariateSpline(k,pk) #,kind='cubic')
	#nk=32768
	nk=1000000
	kmax=numpy.max(k)
	kmin=numpy.min(k)
	kIn=numpy.linspace(kmin,kmax,nk)
	pkIn=pkInter(kIn)
	kIn[0]=0.
	pkIn[0]=0.
	r=2.*numpy.pi/kmax*numpy.arange(nk)
	pkk=kIn*pkIn
	cric=-numpy.imag(numpy.fft.fft(pkk)/nk)/r/2./numpy.pi**2*kmax
	cric[0]=0

	return r,cric

# CAMB de juillet 2013
# /home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/models/DR9LyaMocks_params.ini

# Pour Christophe Magneville :
# change :
# 	transfer_redshift(1) = 0.
# 
#
#
# Run :
#	ombh2 = 0.0227
#	omch2 = 0.1096
#transfer_high_precision = T
#transfer_kmax = 1000
#transfer_k_per_logint = 1000
#transfer_num_redshifts = 1
#transfer_interp_matterpower = F
#transfer_power_var = 7
#transfer_redshift(1) = 0.
#transfer_filename(1) = transfer_out.dat
#transfer_matterpower(1) = matterpower.dat
#	/home/gpfs/manip/mnt0607/bao/hdumasde/Program/camb/camb DR9LyaMocks_params.ini

h = 0.70
omega_m_0 = 0.27

f = cp.fgrowth(2.25,omega_m_0)
print f

camb_baofit = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/models/DR9LyaMocks_matterpower.dat')
k_baofit  = numpy.append( [0.], camb_baofit[:,0])
pk_baofit = numpy.append( [0.], camb_baofit[:,1])
interp = interpolate.interp1d(k_baofit,pk_baofit,bounds_error=False,fill_value=1)

camb_baofit2 = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/models/DR9LyaMocksSB_matterpower.dat')
k_baofit2  = numpy.append( [0.], camb_baofit2[:,0])
pk_baofit2 = numpy.append( [0.], camb_baofit2[:,1])

camb = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Mock_JMLG/Produce_CAMB/DR9LyaMocks_matterpower.dat')
k  = numpy.append( [0.], camb[:,0])
pk = numpy.append( [0.], camb[:,1]*f*f)

camb_CMV = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/cmv/Helion/ginit3d_67_0p0_7859_ntpk.txt')
k_CMV  = numpy.append( [0.], camb_CMV[:,0]/h)
pk_CMV = numpy.append( [0.], camb_CMV[:,3]*numpy.power(h,3.)*f*f)


'''
for i in [0,1,2]:

	coef = numpy.power(k_baofit,i)
	plt.errorbar(k_baofit,coef*pk_baofit, label='with Wiggles',fmt='o')

	coef = numpy.power(k_baofit2,i)
	plt.errorbar(k_baofit2,coef*pk_baofit2, label='no Wiggles',fmt='o')

	coef = numpy.power(k,i)
	plt.errorbar(k,coef*pk, label='mine',fmt='o')

	coef = numpy.power(k_CMV,i)
	plt.errorbar(k_CMV,coef*pk_CMV, label='CMV',fmt='o')


	myTools.deal_with_plot(True,True,True)
	plt.legend(fontsize=30, numpoints=1,ncol=2, loc=3)
	plt.xlabel(r'$k$')
	plt.ylabel(r'$P(k)$')
	plt.show()

for i in [0,1,2]:

	coef = numpy.power(k_baofit,i)
	plt.plot(k_baofit,coef*(pk_baofit-pk_baofit2), label='with Wiggles - no wiggles')

	coef = numpy.power(k,i)
	plt.plot(k,coef*(pk-interp(k)), label='with Wiggles - mine')

	myTools.deal_with_plot(True,True,True)
	plt.legend(fontsize=30, numpoints=1,ncol=2, loc=3)
	plt.xlabel(r'$k$')
	plt.ylabel(r'$P(k)$')
	plt.show()
'''

###
r,cric   = xi_from_pk(k_baofit,pk_baofit)
r2,cric2 = xi_from_pk(k_baofit2,pk_baofit2)
r_mine,cric_mine   = xi_from_pk(k,pk)
r_CMV,cric_CMV = xi_from_pk(k_CMV,pk_CMV)

data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/CAMB/CAMB_2_25/camb.txt')
r3 = data[:,0]
cric3 = data[:,1]

for i in [0,1,2]:

	cut = (r<200.)
	coef = numpy.power(r[cut],i)
	plt.plot(r[cut],coef*cric[cut], label='with Wiggles')

	cut = (r2<200.)
	coef = numpy.power(r2[cut],i)
	plt.plot(r2[cut],coef*cric2[cut], label='no Wiggles')

	cut = (r_mine<200.)
	coef = numpy.power(r_mine[cut],i)
	plt.plot(r_mine[cut],coef*cric_mine[cut], label='Mine')

	cut = (r_CMV<200.)
	coef = numpy.power(r_CMV[cut],i)
	plt.plot(r_CMV[cut],coef*cric_CMV[cut], label='CMV')

	cut = (r3<200.)
	coef = numpy.power(r3[cut],i)
	plt.plot(r3[cut],coef*cric3[cut], label='file')

	myTools.deal_with_plot(False,False,True)
	plt.legend(fontsize=30, numpoints=1,ncol=1, loc=3)
	plt.show()











