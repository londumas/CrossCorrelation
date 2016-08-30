# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#


import scipy
from scipy import interpolate
import numpy
import sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import myTools
import const_delta


pathToData = 'test_1D.fits'
data = pyfits.open(pathToData, memmap=True )

delta = data[0].data.flatten()

plt.plot(delta)
plt.show()
