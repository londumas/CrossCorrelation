# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import astropy.io.fits as pyfits
import numpy
import matplotlib.pyplot as plt


### Values to see histograms
zForHist = []
lambdaForHist = []
lambdaRFForHist = []
fluxForHist = []

### File Jean-Marc, there are 47 spectra
cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_new_generation/TESTS/spectra-780-0.fits')
#'/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/spectra-expander.fits
### Print header, data, ...
print cat[0].header
print cat[0].data
print cat[1].header

for i in range(1,47): ##47

	### Get this hdu
	el = cat[i]

	### Get redshift
	z = el.header['ZQSO']

	### Get data
	el = el.data

	### If you want in lambda_RF
	#el['lambda'] /= (1.+z)

	### plot data
	plt.errorbar( el['lambda'], el['flux'],marker='o')
	plt.ylim([ -0.1, 1.1])

	### Keep values
	zForHist        = numpy.append(z,zForHist)
	lambdaForHist   = numpy.append(el['lambda'],lambdaForHist)
	lambdaRFForHist = numpy.append(el['lambda']/(1.+z),lambdaRFForHist)
	fluxForHist     = numpy.append(el['flux'],fluxForHist)

plt.show()

plt.hist(zForHist)
plt.xlabel(r'$z$', fontsize=40)
plt.ylabel(r'$\#$', fontsize=40)
plt.show()
plt.hist(lambdaForHist)
plt.xlabel(r'$\lambda$', fontsize=40)
plt.ylabel(r'$\#$', fontsize=40)
plt.show()
plt.hist(fluxForHist,bins=50)
plt.xlabel(r'$flux$', fontsize=40)
plt.ylabel(r'$\#$', fontsize=40)
plt.show()


### Only in helion's iclust
import myTools
from myTools import Get_TProfile
from const_delta import *
weight    = numpy.ones(lambdaForHist.size)
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaForHist,fluxForHist, 150,weight)	
plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
myTools.deal_with_plot(False, False, False)
plt.show()

weight    = numpy.ones(lambdaForHist.size)
xxx, yyy, eyyy, nyyy = Get_TProfile(lambdaRFForHist,fluxForHist, 150,weight)	
plt.errorbar(xxx, yyy, yerr=eyyy, marker="o")
myTools.deal_with_plot(False, False, False)
plt.show()

















