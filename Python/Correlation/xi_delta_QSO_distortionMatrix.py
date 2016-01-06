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
max1D__   = 200.
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
nbBinM__ = 100;


path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/'


path = path__ + 'xi_delta_QSO_distortionMatrix_2D_LYA_QSO.txt'
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_delta_QSO_distortionMatrix_2D_LYA_QSO.txt'
print path
data = numpy.loadtxt(path)

print data.size
print data[0].size
print data

data[ data==0. ] = numpy.float('nan')
myTools.plot2D(data)

print '  is diag <0.   : ', numpy.diag(data)[ numpy.diag(data)<0. ].size
print '  is diag ==0.  : ', numpy.diag(data)[ numpy.diag(data)==0. ].size


cor = myTools.getCorrelationMatrix(data)
print cor
print cor[ cor==1. ].size
print cor[ cor>1. ].size
print cor[ cor==0. ].size

cor[ cor==0. ] = numpy.float('nan')
#cor[ cor==1. ] = numpy.float('nan')
#myTools.plot2D(cor)

plt.hist( cor[ numpy.logical_and( cor!=0., cor!=1.) ], bins=100  )
plt.show()

cor2 = numpy.array(cor)
### Test if symetric
size = cor[:,0].size
for i in range(0,size):
	for j in range(0,size):
		cor2[i,j] -= cor[i,j]

cor2[ cor2==0. ] = numpy.float('nan')
cor2[ cor2==1. ] = numpy.float('nan')
myTools.plot2D(cor2)
		

plt.hist( cor2[ numpy.logical_and( cor!=0., cor!=1.) ], bins=100  )
plt.show()
















