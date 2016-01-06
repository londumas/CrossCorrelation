
# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >



import subprocess
import sys
import os

import astropy.io.fits as pyfits
import numpy
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#import decimal
#import profile
from iminuit import Minuit
import root_numpy
import ROOT
#import array
#import re

#import warnings
#warnings.filterwarnings("error")

### My tools
import myTools
from annalys_Xi import *

path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_noFlat2/correlation_notDebug_allSky_DR12.root"
dr12_xi1D, dr12_xi2D, dr12_xi2D_2, dr12_xi2D_2_distX, dr12_xi2D_2_distY = get_Data(path,False)
path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_allFlat2/correlation_notDebug_allSky_DR12.root"
eBOSS_xi1D, eBOSS_xi2D, eBOSS_xi2D_2, eBOSS_xi2D_2_distX, eBOSS_xi2D_2_distY = get_Data(path,False)



path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_Flat/correlation_notDebug_allSky_DR12.root"
dr122_xi1D, dr122_xi2D, dr122_xi2D_2, dr122_xi2D_2_distX, dr122_xi2D_2_distY = get_Data(path,False)
path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodNicolas_Upgrade2/correlation_notDebug_allSky_DR12.root"
eBOSS2_xi1D, eBOSS2_xi2D, eBOSS2_xi2D_2, eBOSS2_xi2D_2_distX, eBOSS2_xi2D_2_distY = get_Data(path,False)

path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m1_FitsFile_DR12_methodHelion_Upgrade/correlation_notDebug_allSky_DR12.root"
eBOSS3_xi1D, eBOSS3_xi2D, eBOSS3_xi2D_2, eBOSS3_xi2D_2_distX, eBOSS3_xi2D_2_distY = get_Data(path,False)

fig = plt.figure()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.)

### Xi1D
coef = -1.
ax = fig.add_subplot(121)
plt.errorbar(dr12_xi1D[0,:],coef*dr12_xi1D[1,:],yerr=coef*dr12_xi1D[2,:],marker="o",label='DR12Q')
plt.errorbar(eBOSS_xi1D[0,:],coef*eBOSS_xi1D[1,:],yerr=coef*eBOSS_xi1D[2,:],marker="o",label='eBOSS')
plt.errorbar(dr122_xi1D[0,:],coef*dr122_xi1D[1,:],yerr=coef*dr122_xi1D[2,:],marker="o",label='DR12Q2')
plt.errorbar(eBOSS2_xi1D[0,:],coef*eBOSS2_xi1D[1,:],yerr=coef*eBOSS2_xi1D[2,:],marker="o",label='eBOSS2')
plt.errorbar(eBOSS3_xi1D[0,:],coef*eBOSS3_xi1D[1,:],yerr=coef*eBOSS3_xi1D[2,:],marker="o",label='eBOSS3')
plt.title(r'', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
plt.ylabel(r'$\xi(|s|)$', fontsize=40)
myTools.deal_with_plot(False,False,True)

### Xi1D*|s|^{2}
coef = -1.*dr12_xi1D[0,:]*dr12_xi1D[0,:]
ax = fig.add_subplot(122)
plt.errorbar(dr12_xi1D[0,:],coef*dr12_xi1D[1,:],yerr=coef*dr12_xi1D[2,:],marker="o",label='DR12Q')
coef = -1.*eBOSS_xi1D[0,:]*eBOSS_xi1D[0,:]
plt.errorbar(eBOSS_xi1D[0,:],coef*eBOSS_xi1D[1,:],yerr=coef*eBOSS_xi1D[2,:],marker="o",label='eBOSS')

coef = -1.*dr122_xi1D[0,:]*dr122_xi1D[0,:]
plt.errorbar(dr122_xi1D[0,:],coef*dr122_xi1D[1,:],yerr=coef*dr122_xi1D[2,:],marker="o",label='DR12Q2')
coef = -1.*eBOSS2_xi1D[0,:]*eBOSS2_xi1D[0,:]
plt.errorbar(eBOSS2_xi1D[0,:],coef*eBOSS2_xi1D[1,:],yerr=coef*eBOSS2_xi1D[2,:],marker="o",label='eBOSS2')
coef = -1.*eBOSS3_xi1D[0,:]*eBOSS3_xi1D[0,:]
plt.errorbar(eBOSS3_xi1D[0,:],coef*eBOSS3_xi1D[1,:],yerr=coef*eBOSS3_xi1D[2,:],marker="o",label='eBOSS3')
plt.title(r'', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi(|s|)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()
