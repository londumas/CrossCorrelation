# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import numpy
import cosmolopy.distance as cosmology
from const_delta import *
import matplotlib.pyplot as plt
import cosmolopy.constants as cc
from scipy import interpolate


### Cosmology
cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}



z = numpy.arange(0.,6.,0.001)
d = cosmology.comoving_distance(z, **cosmo)*h__
convert_d_to_z = interpolate.interp1d(d,z,bounds_error=False,fill_value=0)



def Get_speed(v,zTrue):
	return (1.+zTrue)*v+zTrue
def Get_speed2(v,zTrue):
	vCosmo = ((1.+zTrue)**2.-1.)/((1.+zTrue)**2.+1.)
	tmpv = (vCosmo+v)/( 1.+ vCosmo*v )
	return numpy.sqrt( (1.+tmpv)/(1.-tmpv) ) - 1.


z = 2.5
v = numpy.arange(0.,1.,0.00001)
zObs  = numpy.array( [ Get_speed(i,z) for i in v ] )
zObs2 = numpy.array( [ Get_speed2(i,z) for i in v ] )

plt.plot(v,zObs)
plt.plot(v,zObs2)
plt.grid()
plt.show()


