# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import myTools
from const_delta import *

import sys
import numpy
import matplotlib.pyplot as plt
import pandas
from iminuit import Minuit


nbBin = 100
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_coAdd_2016_05_26/xi_delta_QSO_theta_LYA_QSO.txt'
data = numpy.loadtxt(path)

xi = numpy.zeros( shape=(nbBin,3) )
xi[:,0] = data[:,2]/data[:,3]
xi[:,1] = data[:,0]/data[:,3]
xi[:,2] = numpy.sqrt( (data[:,1]/data[:,3] - xi[:,1]*xi[:,1])/data[:,4])
cut = (data[:,4] == 0.)
xi[:,0][cut] = 0.
xi[:,1][cut] = 0.
xi[:,2][cut] = 0.

plt.errorbar( xi[:,0],xi[:,1],yerr=xi[:,2],fmt='o' )
plt.xlabel(r'$\theta \, [degree]$', fontsize=40)
plt.ylabel(r'$\xi(\theta)$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.show()

