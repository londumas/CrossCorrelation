# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import numpy
import cosmolopy.distance as cosmology
from const_delta import *

def deltaLYA(line):

	verbose = False

	### Cosmology
	cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}

	z_pixel_1 = (4000./1215.67)-1.
	z_pixel_2 = (4000./line)-1.

	if (verbose): print
	if (verbose): print '  z_pixel_1 = ', z_pixel_1
	if (verbose): print '  z_pixel_2 = ', z_pixel_2

	d_pixel_1 = cosmology.comoving_distance(z_pixel_1, **cosmo)*h__
	d_pixel_2 = cosmology.comoving_distance(z_pixel_2, **cosmo)*h__

	if (verbose): print
	if (verbose): print '  d_pixel_1 = ', d_pixel_1
	if (verbose): print '  d_pixel_2 = ', d_pixel_2
	if (verbose): print
	if (verbose): print d_pixel_1-d_pixel_2

	return d_pixel_1-d_pixel_2


LYA_lines_names   = numpy.array(['MGII_a' ,'MGII_b' ,'FeII_a' ,'FeII_b' ,'AlII' ,'CIV_a' ,'CIV_b' ,'SiII' , 'SiIV' ,'CII' ,'OI_b' ,'SiII_b' ,'SiII_a' ,'NV'  ,'LYa'   ,'SiIII'])
LYA_lines         = numpy.array([2804.    ,2796.    ,2383.    ,2344.    ,1671.  ,1551.   ,1548.   ,1527.  , 1394.  ,1335. ,1302.  ,1304.    ,1260.    ,1243. ,1215.67 ,1207.])

LYA_lines         = numpy.array([1193.,1260., 1190.,1207., 1216.24])
for i in range(0,LYA_lines.size):
	print LYA_lines[i], deltaLYA(LYA_lines[i])




