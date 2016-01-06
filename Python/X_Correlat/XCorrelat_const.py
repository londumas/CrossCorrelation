# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

nbRegion = 80

### 1D
min1D   = 0.
max1D   = 160.
nbBin1D = 16

### 2D
minX2D   = min1D
maxX2D   = max1D
nbBinX2D = nbBin1D

minY2D   = -max1D
maxY2D   = max1D
nbBinY2D = 2*nbBin1D

nbBin2D  = nbBinX2D*nbBinY2D
