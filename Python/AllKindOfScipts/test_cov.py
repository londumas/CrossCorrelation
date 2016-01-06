# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import myTools
from const_delta import *

import subprocess
import time
import sys
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit

nbRegion__ = 0

### 1D
min1D__   = 0.
max1D__   = 200. #200.
nbBin1D__ = 50
binSize__ = max1D__/nbBin1D__

### 2D
minX2D__   = min1D__
maxX2D__   = max1D__
nbBinX2D__ = nbBin1D__

minY2D__   = -max1D__
maxY2D__   = max1D__
nbBinY2D__ = 2*nbBin1D__

nbBin2D__  = nbBinX2D__*nbBinY2D__

### Mu
nbBinM__ = 25;


path1__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/'


i = 0
j = 0	
cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/subSampling_LYA_QSO_cov_2D.npy')
### If correlation matrix from another matrix
cor = myTools.getCorrelationMatrix(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_meanSubSampling.npy'))
cov = myTools.getCovarianceMatrix(cor,numpy.diag(cov))

#myTools.plot2D(cov)
#myTools.plot2D(cor)

print cor.flatten()[ cor.flatten()<0.99 ].size
plt.hist(cor.flatten()[ cor.flatten()<0.99 ], bins=10000, histtype='step', label='<corr each simu> ~ 8000 realisation')

cor = myTools.getCorrelationMatrix(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results_RandomPosInCell/subSampling_LYA_QSO_cov_2D.npy'))
plt.hist(cor.flatten()[ cor.flatten()<0.99 ], bins=10000, histtype='step', label='Simu 0 0 = 80 realisation')

cor = myTools.getCorrelationMatrix(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D.npy'))
plt.hist(cor.flatten()[ cor.flatten()<0.99 ], bins=10000, histtype='step', label='100 realisation')

cor = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSampling.npy')
plt.hist(cor.flatten()[ cor.flatten()<0.99 ], bins=10000, histtype='step', label='8000 real')

cor = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy')
plt.hist(cor.flatten()[ cor.flatten()<0.99 ], bins=10000, histtype='step', label='from fit')

plt.xlabel(r'$Cor(s_{1},s_{2})$', fontsize=40)
plt.ylabel(r'$\#$', fontsize=40)
plt.yscale('log')
myTools.deal_with_plot(False, True, True)
plt.show()
















