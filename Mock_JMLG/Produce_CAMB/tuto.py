# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_1D.py
#
"""
		f1 = cp.fgrowth(self.zref,0.27)
		f2 = cp.fgrowth(z,0.27)
		ff = (f2/f1)**2.
*ff
"""

### Python lib
import cosmolopy
import numpy
import matplotlib.pyplot as plt
import cosmolopy.perturbation as cp
from scipy import interpolate
import scipy

import myTools
import const

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
#	/home/gpfs/manip/mnt0607/bao/hdumasde/Program/camb/camb DR9LyaMocks_params.ini

h = 0.70
omega_m_0 = 0.27

f = cp.fgrowth(2.25,omega_m_0)
print f

'''
camb_CMV = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/cmv/Helion/cmvtstpk_Helion.data')
k_CMV  = camb_CMV[:,1]/h
pk_CMV = camb_CMV[:,4]*numpy.power(h,3.)*f*f

camb_baofit = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/models/DR9LyaMocks_matterpower.dat')
k_baofit  = camb_baofit[:,0]
pk_baofit = camb_baofit[:,1]
print '  nb of points baofit : ', k_baofit.size


camb_init = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Mock_JMLG/Produce_CAMB/DR9LyaMocks_matterpower.dat_init')
k_init  = camb_init[:,0]
pk_init = camb_init[:,1]
print '  nb of points init : ', k_init.size
'''

camb_CMV = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/cmv/Helion/cmvtstpk_Helion.data')
k_CMV  = numpy.append( [0.], camb_CMV[:,1]/h)
pk_CMV = numpy.append( [0.], camb_CMV[:,4]*numpy.power(h,3.)*f*f)

camb_CMV3 = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/cmv/Helion/cmvtstpk_Helion.data')
k_CMV3  = numpy.append( [0.], camb_CMV[:,1]/h)
pk_CMV3 = numpy.append( [0.], camb_CMV[:,7]*numpy.power(h,3.)*f*f)




camb_baofit = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Program/baofit/models/DR9LyaMocks_matterpower.dat')
#k_baofit  = numpy.append( [0.], camb_baofit[:,0])
#pk_baofit = numpy.append( [0.], camb_baofit[:,1])
k_baofit  = camb_baofit[:,0]
pk_baofit = camb_baofit[:,1]
print '  nb of points baofit : ', k_baofit.size

### Save CHristophe, no Oscillation
k_new = k_baofit[ k_baofit>k_CMV3[0] ]
k_new = k_new[ k_new<k_CMV3[-1] ]
interp = interpolate.interp1d(k_CMV3,pk_CMV3,bounds_error=False,fill_value=1)
numpy.savetxt('christophe_Eisentien_Hu_Simu_SB_matterpower.dat',zip( k_new,interp(k_new) ))
k_CMV3  = k_new
pk_CMV3 = interp(k_new)
### Save Christophe, with Oscillation
interp = interpolate.interp1d(k_CMV,pk_CMV,bounds_error=False,fill_value=1)
numpy.savetxt('christophe_Eisentien_Hu_Simu_matterpower.dat',zip( k_new,interp(k_new) ))
k_CMV  = k_new
pk_CMV = interp(k_new)

camb = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/cmv/Helion/DR9LyaMocks_matterpower.dat')
k  = numpy.append( [0.], camb[:,0])
pk = numpy.append( [0.], camb[:,1]*f*f)
print '  nb of points      : ', k.size
print k[0]
print k[-1]

camb_CMV_2 = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/cmv/Helion/ginit3d_67_0p0_7859_ntpk.txt')
k_CMV_2  = numpy.append( [0.], camb_CMV_2[:,0]/h)
pk_CMV_2 = numpy.append( [0.], camb_CMV_2[:,3]*numpy.power(h,3.)*f*f)
print '  nb of points      : ', k_CMV_2.size
print k_CMV_2[0]-k_CMV_2[1]

camb_baofit2 = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Tests/helion_matterpower.dat')
k_baofit2  = numpy.append( [0.], camb_baofit2[:,0])
pk_baofit2 = numpy.append( [0.], camb_baofit2[:,1])

'''
for i in numpy.arange(0,k.size):
	if (k[i]>k_CMV_2[0]):
		print i-1
		a = i-1
		break
'''
#interp = interpolate.interp1d(k,pk,bounds_error=False,fill_value=1)
#k_test  = numpy.append( [7.8e-4], k_CMV_2)
#pk_test = numpy.append( interp([7.8e-4]), pk_CMV_2)

#interp = interpolate.interp1d(k,pk,bounds_error=False,fill_value=1)
#k_test2  = numpy.append( [0.], k_CMV_2)
#pk_test2 = numpy.append( [0.], pk_CMV_2)

for i in [0,1,2]:


	coef = numpy.power(k,i)
	plt.plot(k,coef*pk, label='CAMB',marker='o')

	coef = numpy.power(k_baofit,i)
	plt.plot(k_baofit,coef*pk_baofit, label='BAOFIT',marker='o')

	coef = numpy.power(k_baofit2,i)
	plt.plot(k_baofit2,coef*pk_baofit2, label='BAOFIT 2',marker='o')

	coef = numpy.power(k_CMV_2,i)
	plt.plot(k_CMV_2,coef*pk_CMV_2,marker='o',label='Christophe')

	coef = numpy.power(k_CMV,i)
	plt.plot(k_CMV,coef*pk_CMV,marker='o',label='Christophe before')

	coef = numpy.power(k_CMV3,i)
	plt.plot(k_CMV3,coef*pk_CMV3,marker='o',label='Christophe before3')


	myTools.deal_with_plot(True,True,True)
	plt.legend(fontsize=30, numpoints=1,ncol=2, loc=3)
	plt.show()


###
r_CMV,cric_CMV = xi_from_pk(k_CMV,pk_CMV)
r_CMV3,cric_CMV3 = xi_from_pk(k_CMV3,pk_CMV3)
r_CMV_2,cric_CMV_2 = xi_from_pk(k_CMV_2,pk_CMV_2)
r,cric = xi_from_pk(k,pk)
r_baofit,cric_baofit = xi_from_pk(k_baofit,pk_baofit)
r_baofit2,cric_baofit2 = xi_from_pk(k_baofit2,pk_baofit2)

for i in [0,1,2]:

	coef = numpy.power(r,i)
	plt.plot(r,coef*cric, label='CAMB')

	coef = numpy.power(r_CMV,i)
	plt.plot(r_CMV,coef*cric_CMV, label='Christophe before')
	coef = numpy.power(r_CMV3,i)
	plt.plot(r_CMV3,coef*cric_CMV3, label='Christophe before3')

	coef = numpy.power(r_CMV_2,i)
	plt.plot(r_CMV_2,coef*cric_CMV_2, label='Christophe')

	coef = numpy.power(r_baofit,i)
	plt.plot(r_baofit,coef*cric_baofit, label='BAOFIT')

	coef = numpy.power(r_baofit2,i)
	plt.plot(r_baofit2,coef*cric_baofit2, label='BAOFIT 2')

	myTools.deal_with_plot(False,False,True)
	plt.xlim( [0.,200.] )
	cut = (r_CMV_2<200.)
	plt.ylim( [ -5.,10. ] )
	plt.legend(fontsize=30, numpoints=1,ncol=1, loc=3)
	plt.show()



























