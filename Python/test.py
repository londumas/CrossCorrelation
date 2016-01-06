import subprocess
import sys
import os

import numpy
import astropy.io.fits as pyfits
from iminuit import Minuit
import matplotlib.pyplot as plt
import re
#import array
#import matplotlib.patches as mpatches
#import decimal
#import profile
import warnings
warnings.filterwarnings("error")

### My tools
import myTools
from myTools import Get_TProfile
### The Constants
from delta_const import *

### "HOME" or "ICLUST"
location__ = 'ICLUST'
### "DR12" or 'Guy' or 'Margala' or "MOCK" or 'eBOSS'
pipeline__ = 'eBOSS'
### False or True
reObs__ = False

##/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_8_0/7153/spec-7153-56904-0004.fits
cat = pyfits.open('/home/gpfs/manip/mnt/bao/Spectra/SpectraV5_8_0/spPlate-7153-56904.fits', memmap=True)

fiber = 3

startValue = 3.55190000000
stepValue = 0.000100000000000
print cat[1].header

print cat[0].data[fiber]
print cat[1].data[fiber]
print cat[2].data[fiber]
print cat[3].data[fiber]
print cat[4].data[fiber]
print cat[5].data[fiber]
print cat[6].data[fiber]


nbPixel = len(cat[1].data[0])
logLambda = numpy.power(10.,numpy.arange(0,nbPixel)*stepValue+startValue)
print len(logLambda)

plt.plot(logLambda,cat[1].data[4] )
plt.plot(logLambda,cat[1].data[3] )
plt.show()
plt.plot(logLambda,cat[2].data[4] )
plt.plot(logLambda,cat[2].data[3] )
plt.show()
plt.plot(logLambda,cat[3].data[4] )
plt.plot(logLambda,cat[3].data[3] )
plt.show()
plt.plot(logLambda,cat[4].data[4] )
plt.plot(logLambda,cat[4].data[3] )
plt.show()




