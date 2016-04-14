# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

### My tools
import myTools


import astropy.io.fits as pyfits
import math
import numpy
import os
import decimal ## To set the precision of the double values
import cosmolopy.distance as cosmology
import matplotlib.pyplot as plt
from iminuit import Minuit
import scipy
import sys

import myTools
import CAMB


def get_parameter_a(z):
	return numpy.power(10., -7.04712 + 2.88486*z -0.300941*z*z)

correlation_type = int(sys.argv[1])
### Constants
size_gauss = 10000000

a       = 0.077
gamma   = 2.16
cut_QSO = 1.

"""
z = 2.25
a       = get_parameter_a(z)
print '  a = ', a
gamma   = 1.6
cut_QSO = 2.3793
"""


### Identite
if (correlation_type==0):
	prefix = 'identite'
	def F1(x):
		return (1.+x)
	def F2(x):
		return (1.+x)
	def lim_min_y(x):
		return -numpy.inf+x*0.
	def lim_max_y(x):
		return numpy.inf+x*0.

### Auto-correlation
if (correlation_type==1):
	prefix = 'autoCorrelation'
	def F1(x):
		return numpy.exp(-a*numpy.exp(gamma*x) )
	def F2(x):
		return numpy.exp(-a*numpy.exp(gamma*x) )
	def lim_min_y(x):
		return -numpy.inf+x*0.
	def lim_max_y(x):
		return numpy.inf+x*0.

### Cross-correlation
if (correlation_type==2):
	prefix = 'CrossCorrelation'
	def F1(x):
		x = numpy.array(x)
		x[(x<=cut_QSO)] = 0.
		return x
	def F2(x):
		return numpy.exp(-a*numpy.exp(gamma*x) )
	def lim_min_y(x):
		return cut_QSO+x*0.
	def lim_max_y(x):
		return numpy.inf+x*0.

### Auto-correlation QSO
if (correlation_type==3):
	prefix = 'autoCorrelationQSO'
	def F1(x):
		x = numpy.array(x)
		x[(x<=cut_QSO)] = 0.
		return x
	def F2(x):
		x = numpy.array(x)
		x[(x<=cut_QSO)] = 0.
		return x
	def lim_min_y(x):
		return cut_QSO+x*0.
	def lim_max_y(x):
		return numpy.inf+x*0.

def F1_for_Mean(x):
	return F1(x)*numpy.exp(-0.5*x*x)/numpy.sqrt(2.*numpy.pi)
def F2_for_Mean(x):
	return F2(x)*numpy.exp(-0.5*x*x)/numpy.sqrt(2.*numpy.pi)
'''
def F1_times_F2(x,y,z):
	t  = z*x+numpy.sqrt(1.-z*z)*y
	f1 = F1(x)*numpy.exp(-0.5*x*x)
	f2 = F2(t)*numpy.exp(-0.5*y*y)
	return f1*f2/(2.*numpy.pi)
'''

def correlation_with_random(xi_gauss):
	print '  inside corr'
	delta_gauss_1 = numpy.random.normal(size=(size_gauss,xi_gauss.size))
	delta_gauss_2 = numpy.random.normal(size=(size_gauss,xi_gauss.size))
	print '  inside corr 2'
	integral_1 = F1(delta_gauss_1)
	integral_2 = F2(xi_gauss*delta_gauss_1 + numpy.sqrt(1.-numpy.power(xi_gauss,2.))*delta_gauss_2)
	print '  inside corr 3'

	return numpy.mean( integral_1*integral_2, axis=0  )
def produce_correlation_with_random():

	### Get mean F1
	delta_gauss = numpy.random.normal(size=size_gauss)
	plt.plot(delta_gauss,F1(delta_gauss))
	plt.show()
	mean_F1 = numpy.mean( F1(delta_gauss) )
	print '  < F1 >  = ', mean_F1
	### Get mean F2
	delta_gauss = numpy.random.normal(size=size_gauss)
	mean_F2 = numpy.mean( F2(delta_gauss) )
	print '  < F2 >  = ', mean_F2

	transformed_xi_0 = correlation_with_random(xi_0)
	transformed_xi_0 = transformed_xi_0/(mean_F1*mean_F2) - 1.

	numpy.savetxt('transformed_function_'+prefix+'.txt', zip(xi_r,transformed_xi_0))

	return transformed_xi_0
def produce_correlation():

	x = numpy.arange(-10.,10.,0.01)
	y = F1_for_Mean(x)
	plt.plot(x,y,marker='o')
	y = F2_for_Mean(x)
	plt.plot(x,y,marker='o')
	plt.show()

	### Get mean F1
	mean_F1 = scipy.integrate.quad(F1_for_Mean,-numpy.inf,numpy.inf)[0]
	print '  < F1 >  = ', mean_F1
	### Get mean F2
	mean_F2 = scipy.integrate.quad(F2_for_Mean,-numpy.inf,numpy.inf)[0]
	print '  < F2 >  = ', mean_F2

	def get_transformation(z):

		def lim_min_y(x):
			return cut_QSO+x*0.
		def lim_max_y(x):
			return numpy.inf+x*0.

		def F1_times_F2(x,y):
			t  = z*x+numpy.sqrt(1.-z*z)*y
			f1 = F1(x)*numpy.exp(-0.5*x*x)
			f2 = F2(t)*numpy.exp(-0.5*y*y)
			return f1*f2/(2.*numpy.pi)
		

		a = scipy.integrate.dblquad(F1_times_F2, -numpy.inf,numpy.inf, lim_min_y, lim_max_y)[0]/(mean_F1*mean_F2) - 1.

		def lim_min_y(x):
			return -numpy.inf+x*0.
		def lim_max_y(x):
			return numpy.inf+x*0.
		print scipy.integrate.dblquad(F1_times_F2, -numpy.inf,numpy.inf, lim_min_y, lim_max_y)[0]/(mean_F1*mean_F2) - 1.
		print z, a

		return a

	transformed_xi_0 = [ get_transformation(el) for el in xi_0 ]
	numpy.savetxt('transformed_function_'+prefix+'.txt', zip(xi_r,transformed_xi_0))

	return transformed_xi_0

def read_correlation():

	transformed_xi_0 = numpy.loadtxt('transformed_function_'+prefix+'.txt')[:,1]
	return transformed_xi_0



























camb = CAMB.CAMB('CAMB_Mocks_me')
xi_r = numpy.arange(1,200.,10.)
xi_0 = numpy.interp(xi_r,camb._xi0[:,0],camb._xi0[:,1])

#transformed_xi_0 = produce_correlation_with_random()
transformed_xi_0 = produce_correlation()
#transformed_xi_0 = read_correlation()

bias = 1. #-0.2

for coef_index in [0,1,2]:
	###
	coef = numpy.power(camb._xi0[:,0],coef_index)
	plt.plot(camb._xi0[:,0],coef*camb._xi0[:,1], label=r'$Input$')
	###
	coef = numpy.power(xi_r,coef_index)
	plt.errorbar(xi_r,coef*xi_0, label=r'$Coarse$',fmt='o')
	###
	coef = numpy.power(xi_r,coef_index)
	plt.errorbar(xi_r,coef*transformed_xi_0/bias, label=r'$Transformed$',fmt='o')
	###
	myTools.deal_with_plot(False,False,True)
	plt.xlim([0,200.])
	plt.show()



plt.errorbar(xi_r,transformed_xi_0/xi_0,fmt='o')
myTools.deal_with_plot(False,False,True)
plt.show()

###
plt.errorbar(xi_0,transformed_xi_0, label=r'$input \, vs. \, output$',fmt='o')

xxx = xi_0
yyy = transformed_xi_0
def chi2(a0):
	return numpy.sum(numpy.power( yyy-(a0*xxx) ,2.))
m = Minuit(chi2,a0=1.,error_a0=1.,print_level=-1, errordef=0.01) 	
m.migrad()
a0 = m.values['a0']
print '  a0 = ', a0

plt.plot( xxx,a0*xxx )

###
myTools.deal_with_plot(False,False,True)
plt.show()















