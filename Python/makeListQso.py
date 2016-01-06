# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

### My tools
import myTools


import astropy.io.fits as pyfits
import math
import numpy
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import decimal ## To set the precision of the double values

## Cosmology
import cosmolopy.distance as cosmology


bit16__ = 65536
lenSpectrum = 4619

### Lambda_{R.F.}
########################################################################

### Observed wavelenght range
lambdaObsMin__ = 3600.
lambdaObsMax__ = 7235.
log10lambdaObsMin__ = numpy.log10(lambdaObsMin__)
log10lambdaObsMax__ = numpy.log10(lambdaObsMax__)

### Sky lines
#############
skyLines__ = [ (3615,3619),(3932,3937),(3966,3972),(4042,4050),
	(4357,4362),(5458,5467),(5573,5585),(5682,5695),
	(5885,5902),(6235,6241),(6256,6263),(6296,6311),
	(6320,6334),(6362,6369),(6498,6502),(6554,6557),
	(6825,6840),(6862,6870),(6922,6928),(6948,6954),
	(6977,6982) ]
skyLines__ = numpy.log10(skyLines__)

### Nuber of bin for mean transmission flux
nbBinObs__ = int((lambdaObsMax__-lambdaObsMin__)/1.)

########################################################################
### Lambda_{R.F.}
########################################################################

### Rest Frame wavelenght range for template
lambdaRFTemplateMin__ = 1037.
lambdaRFTemplateMax__ = 1203.
log10lambdaRFTemplateMin__ = numpy.log10(lambdaRFTemplateMin__)
log10lambdaRFTemplateMax__ = numpy.log10(lambdaRFTemplateMax__)

### Rest Frame wavelenght range for data
lambdaRFMin__  = 1040.
lambdaRFMax__  = 1200.
lambdaRFMean__ = (lambdaRFMax__+lambdaRFMin__)/2.

### Normalisation region
lambdaRFNormaMin__ = 1275.
lambdaRFNormaMax__ = 1285.
log10lambdaRFNormaMin__ = numpy.log10(lambdaRFNormaMin__)
log10lambdaRFNormaMax__ = numpy.log10(lambdaRFNormaMax__)

### Lyman-alpha line
lambdaRFLya__ = 1215.67

### Min and max number of pixels in the forest
nbBinRFMin__ = 50
nbBinRFMax__ = 700

### Nuber of bin for template
nbBinRFTemplate__ = int((lambdaRFTemplateMax__-lambdaRFTemplateMin__)/1.)


########################################################################
### Weight
########################################################################

### Number of pixel for the weight
nbBinWeight__ = 100

### Min of delta_Ivar_pipeline
minDeltaIvar__ = 0.001

### Z_{0} is used to normalise
z0__        = 2.25
onePlusZ0__ = 1.+z0__
gama__      = 3.8
halfGama__  = gama__/2.

########################################################################

### Minimal redshift to get first pixel
minRedshift__ = lambdaObsMin__/lambdaRFMax__ - 1.
maxRedshift__ = lambdaObsMax__/lambdaRFNormaMin__ - 1.




def find_range_QSO_redshift():
	'''
	'''
	
	### range of the correlation
	maxX__ = 160.

	### Cosmology
	omegaM0__      = 0.27
	omegaLambda0__ = 0.73
	h__            = 0.71
	cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}

	### Distance max of correlation
	dist_max = maxX__*numpy.sqrt(2.)/h__
	print '  dist_max = ', dist_max

	z_pixel_min = (lambdaObsMin__/lambdaRFLya__)-1.
	z_pixel_max = (lambdaObsMax__/lambdaRFLya__)-1.

	print
	print '  z_pixel_min = ', z_pixel_min
	print '  z_pixel_max = ', z_pixel_max

	d_pixel_min = cosmology.comoving_distance(z_pixel_min, **cosmo)
	d_pixel_max = cosmology.comoving_distance(z_pixel_max, **cosmo)

	print
	print '  d_pixel_min = ', d_pixel_min
	print '  d_pixel_max = ', d_pixel_max

	d_QSO_min = d_pixel_min-dist_max
	d_QSO_max = d_pixel_max+dist_max

	print
	print '  d_QSO_min = ', d_QSO_min
	print '  d_QSO_max = ', d_QSO_max

	distfunc, redfunc = cosmology.quick_distance_function(cosmology.comoving_distance, zmax=20.0, zmin=0.0, zstep=0.001, return_inverse=True, **cosmo)
	d1 = distfunc(z_pixel_min)
	d2 = distfunc(z_pixel_max)
	print '\n',d1, d2

	z_QSO_min = redfunc(d_QSO_min)
	z_QSO_max = redfunc(d_QSO_max)

	print
	print '  z_QSO_min = ', z_QSO_min
	print '  z_QSO_max = ', z_QSO_max
	
	return z_QSO_min,z_QSO_max

def plot_all_cat():
	'''

	'''

	redMin = 1.7
	redMax = 6.

	### All known qsos
	print "  --- All known ---"
	redShiftKey = 'ZEM'
	cat = pyfits.open("/home/gpfs/manip/mnt/bao/hdumasde/Catalogue/knownquasarstar.060910.fits")[1].data
	print " Size of cat = " + str(len(cat))
	
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(redMax)

	dr7_ra  = cat['RA']
	dr7_dec = cat['DEC']
	dr7_z   = cat['ZEM']
	
	### DR12Q
	print "  --- DR12 ---"
	redShiftKey = 'Z_VI'
	cat = pyfits.open("/home/gpfs/manip/mnt/bao/hdumasde/Catalogue/DR12Q_v2_10.fits")[1].data
	print " Size of cat = " + str(len(cat))
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(redMax)
	
	dr12_ra  = cat['RA']
	dr12_dec = cat['DEC']
	dr12_z   = cat['Z_VI']
	
	### eBOSS
	print "  --- eBOSS ---"
	redShiftKey = 'Z'
	cat = pyfits.open("/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_8_0.fits")[1].data
	print " Size of cat = " + str(len(cat))
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(redMax)
	print len(cat)
	cat = cat[ cat["CLASS"] == "QSO" ]
	print len(cat)
	cat = cat[ cat["ZWARNING"] == 0 ]
	print len(cat)
	cat = cat[ cat["Z_ERR"] > 0 ]
	print len(cat)
	print
	
	eBOSS_ra  = cat['RA']
	eBOSS_dec = cat['DEC']
	eBOSS_z   = cat['Z']
	
	## Map
	plt.plot(dr7_ra,dr7_dec,     linestyle="", marker="o", label=r'$Known \, before \, DR12$')
	plt.plot(dr12_ra,dr12_dec,   linestyle="", marker="o", label=r'$DR12Q$')
	plt.plot(eBOSS_ra,eBOSS_dec, linestyle="", marker="o", label=r'$eBOSS$')
	plt.xlabel(r'$R. A. \, [^\circ]$', fontsize=40)
	plt.ylabel(r'$Dec. \, [^\circ]$', fontsize=40)
	plt.xlim([0.,360.])
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	## Distribution redshift
	plt.hist(dr7_z,                                            bins=100, histtype='step', label=r'$Known \, before \, DR12$', linewidth=2.)
	plt.hist(dr12_z,                                           bins=100, histtype='step', label=r'$DR12Q$', linewidth=2.)
	plt.hist(eBOSS_z,                                          bins=100, histtype='step', label=r'$eBOSS$', linewidth=2.)
	plt.hist(numpy.append(numpy.append(dr7_z,dr12_z),eBOSS_z), bins=100, histtype='step', label=r'$All$', linewidth=2.)
	plt.xlabel(r'$z_{QSO}$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
def plot_all_cat_Forest():
	'''

	'''
	
	### DR12Q
	path = "/home/gpfs/manip/mnt/bao/hdumasde/Lists/DR12Q_v2_10.fits"		
	cat = pyfits.open(path, memmap=True )
	cat = cat[1].data
	print "  The size of the catalogue is           : " + str(len(cat))
	cat = cat[ cat["Z_VI"]>minRedshift__]
	print "  We keep Z_VI > " + str(minRedshift__) + "  , the size is      : " + str(len(cat))
	#cat = cat[ cat["Z_VI"]<=maxRedshift__]
	cat = cat[ cat["Z_VI"]<=4.7]
	print "  We keep Z_VI <= " + str(maxRedshift__) + "  , the size is      : " + str(len(cat))
	cat = cat[ cat["BAL_FLAG_VI"]==0]
	print "  We keep BAL_FLAG_VI == 0 , the size is : " + str(len(cat))

	dr12_ra  = cat['RA']
	dr12_dec = cat['DEC']
	dr12_z   = cat['Z_VI']

	### eBOSS
	print "  --- eBOSS ---"
	redShiftKey = 'Z'
	cat = pyfits.open("/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_8_0.fits")[1].data
	print " Size of cat = " + str(len(cat))
	cat = cat[ cat[redShiftKey] > minRedshift__ ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(minRedshift__)
	cat = cat[ cat[redShiftKey] < 4.7 ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(maxRedshift__)
	print len(cat)
	cat = cat[ cat["CLASS"] == "QSO" ]
	print len(cat)
	cat = cat[ cat["ZWARNING"] == 0 ]
	print len(cat)
	cat = cat[ cat["Z_ERR"] > 0 ]
	print len(cat)
	print
	
	eBOSS_ra  = cat['RA']
	eBOSS_dec = cat['DEC']
	eBOSS_z   = cat['Z']
	
	## Map
	plt.plot(dr12_ra,dr12_dec,   linestyle="", marker="o", label=r'$DR12Q$')
	plt.plot(eBOSS_ra,eBOSS_dec, linestyle="", marker="o", label=r'$eBOSS$')
	plt.xlabel(r'$R. A. \, [^\circ]$', fontsize=40)
	plt.ylabel(r'$Dec. \, [^\circ]$', fontsize=40)
	plt.xlim([0.,360.])
	myTools.deal_with_plot(False,False,True)
	plt.show()
	## Distribution redshift
	plt.hist(dr12_z,                       bins=100, histtype='step', label=r'$DR12Q$', linewidth=2.)
	plt.hist(eBOSS_z,                      bins=100, histtype='step', label=r'$eBOSS$', linewidth=2.)
	plt.hist(numpy.append(dr12_z,eBOSS_z), bins=100, histtype='step', label=r'$All$', linewidth=2.)
	plt.xlabel(r'$z_{forest}$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	plt.xlim([1.6,5.])
	plt.title(r'$Forests \, Redshift \, Distribution$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	cat = pyfits.open('/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile/FitsFile_DR12/allDR12.fits')[1].data
	ar_cut          = (cat['DELTA_IVAR']>minDeltaIvar__)
	ar_lambdaObs    = cat['LAMBDA_OBS'][ ar_cut ]
	ar_weight       = numpy.power(ar_lambdaObs/lambdaRFLya__, halfGama__)/(0.10+1./(1.*cat['DELTA_IVAR'][ar_cut]))
	z               = ar_lambdaObs/lambdaRFLya__ -1.

	plt.hist(z, bins=100, histtype='step', linewidth=2., normed=True, label=r'$DR12Q$')
	plt.hist(z, bins=100, histtype='step', linewidth=2., normed=True, weights=ar_weight,label=r'$DR12Q, weighted$')
	plt.xlabel(r'$z_{pixel}$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	plt.xlim([1.6,5.])
	plt.ylim([0.,1.7])
	plt.title(r'$Pixel \, Redshift \, Distribution \, (normalized)$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	

	return
def getQsoCatalogueDR7():
	'''
		Get the catalogue of qsos for DR7
	'''
	
	redMin = 0.
	redMax = 20.
	
	redShiftKey = "Z"
	
	#cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/dr7qso.fit")[1].data
	#print pyfits.open("/home/helion/Documents/Thèse/Data/Fits/dr7qso.fit")[1].header
	
	cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/knownquasarstar.060910.fits")[1].data
	#print pyfits.open("/home/helion/Documents/Thèse/Data/Fits/knownquasarstar.060910.fits")[1].header
	redShiftKey = "ZEM"
	
	print "  --- DR7 ---"
	print " Size of cat = " + str(len(cat))
	
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(redMax)
	print " Size of cat after cut in Z = " + str(len(cat[ cat[redShiftKey] >= 2. ])) + " qso with z >= " + str(2.)
	
	file_name = 'QSO_DR7.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat["RA"],cat["DEC"],cat[redShiftKey]))
	file.close()
	
	print min(cat[redShiftKey])
	print max(cat[redShiftKey])
	
	## Map
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.plot(cat["RA"], cat["DEC"], linestyle="", marker="o")
	plt.xlabel("Right Ascension (degree)")
	plt.ylabel("Declination (degree)")
	plt.title("SDSS DR7")
	plt.show()
	## Distribution redshift
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.hist(cat[redShiftKey], bins=50)
	plt.xlabel("Z")
	plt.ylabel("#")
	plt.title("SDSS DR7")
	plt.show()
	
	return
def getQsoCatalogueDR12():
	'''
		Get the catalogue of qsos
	'''
	redMin = 0.
	redShiftKey = "Z_VI"
	
	cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/DR12Q_v2_10.fits")[1].data
	
	print "  --- DR12 ---"
	print " Size of cat = " + str(len(cat))
	
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	print " Size of cat after cut in Z = " + str(len(cat[ cat[redShiftKey] >= 2. ])) + " qso with z >= " + str(2.)
	
	
	
	file_name = 'QSO_DR12.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat["RA"],cat["DEC"],cat[redShiftKey]))
	file.close()
	
	
	print min(cat[redShiftKey])
	print max(cat[redShiftKey])
	
	## Map
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.plot(cat["RA"], cat["DEC"], linestyle="", marker="o")
	plt.xlabel("Right Ascension (degree)")
	plt.ylabel("Declination (degree)")
	plt.title("BOSS DR12")
	plt.show()
	## Distribution redshift
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.hist(cat[redShiftKey], bins=50)
	plt.xlabel("Z")
	plt.ylabel("#")
	plt.title("BOSS DR12")
	plt.show()
	
	return
def getQsoCatalogueEBOSS():
	'''
		Get the catalogue of qsos
	'''
	redMin = 0.
	redShiftKey = "Z"
	
	cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/spAll-v5_8_0.fits")[1].data
	
	##print pyfits.open("/home/helion/Documents/Thèse/Data/Fits/spAll-v5_7_9.fits")[1].header
	
	print "  --- eBOSS ---"
	print " Size of cat = " + str(len(cat))
	
	cat = cat[ cat["CLASS"] == "QSO" ]
	print " Size of cat only qso = " + str(len(cat))
	
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	print " Size of cat after cut in Z = " + str(len(cat[ cat[redShiftKey] >= 2. ])) + " qso with z >= " + str(2.)
	
	#cat = cat[ cat["ZWARNING"] == 0 ]
	#print " Size of cat after no ZWARNING = " + str(len(cat))
	
	file_name = 'QSO_eBOSS.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat["RA"],cat["DEC"],cat[redShiftKey]))
	file.close()
	
	print min(cat[redShiftKey])
	print max(cat[redShiftKey])
	
	## Map
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.plot(cat["RA"], cat["DEC"], linestyle="", marker="o")
	plt.xlabel("Right Ascension (degree)")
	plt.ylabel("Declination (degree)")
	plt.title("eBOSS")
	plt.show()
	## Distribution redshift
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.hist(cat[redShiftKey], bins=50)
	plt.xlabel("Z")
	plt.ylabel("#")
	plt.title("eBOSS")
	plt.show()
	
	return
def getQsoCatalogue2EBOSS():
	'''
		Get the catalogue of qsos
	'''
	redMin = 1.8
	
	cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/Superset_EDR_v1_2.fits")[1].data
	
	print pyfits.open("/home/helion/Documents/Thèse/Data/Fits/Superset_EDR_v1_2.fits")[1].header
	
	## if Superset_EDR_v1_2.fit
	
	print "  --- eBOSS ---"
	print " Size of cat = " + str(len(cat))
	
	z_VI         = cat["Z_VI"]
	class_VI     = cat["CLASS_PERSON"]
	z_conf       = cat["Z_CONF_PERSON"]
	z_auto       = cat["Z_AUTO"]
	class_auto   = cat["CLASS_AUTO"]
	redshift     = [0. for elem in z_VI]
	
	for i in range(0, len(z_VI)):
		if (z_VI[i]>0.0 and (class_VI[i]==3 or class_VI[i]==30) and z_conf[i]==3):
			redshift[i] = z_VI[i]
		elif (z_auto[i]>0.0 and class_auto[i]==3):
			redshift[i] = z_auto[i]
			
	ra = []
	de = []
	zz = []
	for i in range(0, len(redshift)):
		if (redshift[i] >= redMin):
			ra += [cat["RA"][i]]
			de += [cat["DEC"][i]]
			zz += [redshift[i]]
	
	file_name = 'QSO_eBOSS_test.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(ra,de,zz))
	file.close()
	
	print min(zz)
	print max(zz)
	
	## if spAll-v5_7_9.fits
	'''
	print "  --- eBOSS ---"
	print " Size of cat = " + str(len(cat))
	
	##cat = cat[ cat["CLASS"] == "QSO" ]
	##print " Size of cat only qso = " + str(len(cat))
	
	z_VI         = cat["Z_VI"]
	class_VI     = cat["CLASS_PERSON"] ##CLASS_PERSON ##class_VI
	z_conf       = cat["Z_CONF_PERSON"]
	z_auto       = cat["Z_AUTO"]
	class_auto   = cat["CLASS_AUTO"]
	redshift     = [0. for elem in z_VI]
	
	print z_auto
	
	for i in range(0, len(z_VI)):
		if (z_VI[i]>0.0 and (class_VI[i]==3 or class_VI[i]==30) and z_conf[i]==3):
			redshift[i] = z_VI[i]
		elif (z_auto[i]>0.0 and class_auto[i]==3):
			redshift[i] = z_auto[i]
	
	#cat = cat[ cat["Z"] >= redMin ]
	#print " Size of cat after cut in Z = " + str(len(cat))
	
	##cat = cat[ cat["ZWARNING"] == 0 ]
	##print " Size of cat after no ZWARNING = " + str(len(cat))
	
	file_name = 'QSO_eBOSS_test.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat["RA"],cat["DEC"],cat["Z"]))
	file.close()
	
	print min(cat["Z"])
	print max(cat["Z"])
	'''
	
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.plot(cat["RA"], cat["DEC"], linestyle="", marker="o")
	plt.xlabel("Right Ascension (degree)")
	plt.ylabel("Declination (degree)")
	plt.title("eBOSS")
	plt.show()
	
	return

def plotAllCatalog():
	'''
	'''
	
	array = numpy.loadtxt("/home/helion/Documents/Thèse/Data/List/QSO_DR12_DR7_eBOSS_Helion.txt")
	
	zz = array[:,2]
	print zz
	print min(zz)
	print max(zz)
	
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.plot(array[:,0], array[:,1], linestyle="", marker="o")
	plt.xlabel("Right Ascension (degree)")
	plt.ylabel("Declination (degree)")
	plt.title("Final catalog")
	plt.show()
	
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.hist(zz, bins=50)
	plt.grid()
	plt.xlabel("Z")
	plt.ylabel("#")
	plt.title("Distribution of Quasars's redshift")
	plt.show()
	
	return

def saveFlatFileEBOSS():
	'''
		Save information in an ascii file so Christophe can run his jobs
	'''
	
	#cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/spAll-v5_8_0.fits")[1].data
	cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/Superset_EDR_v1_2.fits")[1].data
	
	print len(cat)
	cat = cat[ cat["Z_VI"] >= 1.8 ]
	print len(cat)
	##cat = cat[ cat["CLASS"] == "QSO" ]
	##cat = cat[ cat["ZWARNING"] == 0 ]
	
	iBAL    = [0 for i in range(0, len(cat)) ]
	iDLA    = [0 for i in range(0, len(cat)) ]
	
	magArr  = cat["PSFMAG"]
	umag    = [elem for elem in magArr[:,0 ] ]
	gmag    = [elem for elem in magArr[:,1 ] ]
	rmag    = [elem for elem in magArr[:,2 ] ]
	imag    = [elem for elem in magArr[:,3 ] ]
	zmag    = [elem for elem in magArr[:,4 ] ]
	
	file_name = "list_forest_eBOSS_Superset.txt"
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat["PLATE"], cat["MJD"], cat["FIBERID"], cat["RA"], cat["DEC"], cat["Z_VI"], umag, gmag, rmag, imag, zmag, iBAL, iDLA),fmt='%d %d %d %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %d %d')
	file.close()
	
	return
def saveFlatFileEBOSS_Superset():
	'''
		Save information in an ascii file so Christophe can run his jobs
	'''
	
	cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/Superset_EDR_v1_2.fits")[1].data
	print pyfits.open("/home/helion/Documents/Thèse/Data/Fits/Superset_EDR_v1_2.fits")[1].header
	print "  --- eBOSS ---"
	print " Size of cat = " + str(len(cat))
	
	z_VI         = cat["Z_VI"]
	class_VI     = cat["CLASS_PERSON"]
	z_conf       = cat["Z_CONF_PERSON"]
	z_auto       = cat["Z_AUTO"]
	class_auto   = cat["CLASS_AUTO"]
	redshift     = [0. for elem in z_VI]
	
	for i in range(0, len(z_VI)):
		if (z_VI[i]>0.0 and (class_VI[i]==3 or class_VI[i]==30) and z_conf[i]==3):
			redshift[i] = z_VI[i]
		elif (z_auto[i]>0.0 and class_auto[i]==3):
			redshift[i] = z_auto[i]
	
	iBAL    = cat["BAL_FLAG_VI"]
	iDLA    = cat["DLA_FLAG_VI"]
	
	magArr  = cat["PSFMAG"]
	umag    = [elem for elem in magArr[:,0 ] ]
	gmag    = [elem for elem in magArr[:,1 ] ]
	rmag    = [elem for elem in magArr[:,2 ] ]
	imag    = [elem for elem in magArr[:,3 ] ]
	zmag    = [elem for elem in magArr[:,4 ] ]
	
	file_name = "list_forest_eBOSS_all.txt"
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat["PLATE"], cat["MJD"], cat["FIBERID"], cat["RA"], cat["DEC"], redshift, umag, gmag, rmag, imag, zmag, iBAL, iDLA),fmt='%d %d %d %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %d %d')
	file.close()
	
	return
def saveFlatFile():
	'''
		Save information in an ascii file so Christophe can run his jobs
	'''
	
	fitsList = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/DR12Q_v2_10.fits")[1].data
	
	plateId = [str(elem) for elem in fitsList["PLATE"] ]
	mjd     = [str(elem) for elem in fitsList["MJD"] ]
	fiberId = [str(elem) for elem in fitsList["FIBERID"] ]
	ra      = [str(elem) for elem in fitsList["RA"] ]
	dec     = [str(elem) for elem in fitsList["DEC"] ]
	z       = [str(elem) for elem in fitsList["Z_VI"] ]
	iBAL    = [str(elem) for elem in fitsList["BAL_FLAG_VI"] ]
	iDLA    = str(0)
	
	magArr  = fitsList["PSFMAG"]
	umag    = [str(elem) for elem in magArr[:,0 ] ]
	gmag    = [str(elem) for elem in magArr[:,1 ] ]
	rmag    = [str(elem) for elem in magArr[:,2 ] ]
	imag    = [str(elem) for elem in magArr[:,3 ] ]
	zmag    = [str(elem) for elem in magArr[:,4 ] ]
	
	plateReObs = fitsList["PLATE_DUPLICATE"]
	mjdReObs   = fitsList["MJD_DUPLICATE"]
	fiberReObs = fitsList["FIBERID_DUPLICATE"]
	
	print len(ra)
	print len(fitsList)
	
	del fitsList

	string = ""

	for k in range(0, len(ra) ):
		
		if (int(iBAL[k]) > 0):
			print iBAL[k]
			print k
			print 
			break
		
		plateReObs_k = plateReObs[k][ plateReObs[k] > -1 ]
		mjdReObs_k   = mjdReObs[k][ mjdReObs[k] > -1 ]
		fiberReObs_k = fiberReObs[k][ fiberReObs[k] > -1 ]
		NbObs_k      = len(plateReObs_k)+1
		
		string += plateId[k] + "  " + mjd[k] + "  " + fiberId[k] + "  " + ra[k] + "  " + dec[k] + "  " + z[k] 
		string += "  " + umag[k] + "  " + gmag[k] + "  " + rmag[k] + "  " + imag[k] + "  " + zmag[k]
		string += "  " + iBAL[k] + "  " + iDLA + "  " + str(NbObs_k)
		
		if (NbObs_k > 1):
			for j in range(0, NbObs_k-1):
				plateReObs_j = str(plateReObs_k[j])
				mjdReObs_j = str(mjdReObs_k[j])
				fiberReObs_j = str(fiberReObs_k[j])
					
				string += "  " + plateReObs_j + "  " + mjdReObs_j + "  " + fiberReObs_j
		
		string += "\n"
	
	'''
	fileOpen = open("list.txt", "w")
	fileOpen.write(string)
	fileOpen.close()
	'''
	'''
	ofstream fileList("NewList.txt");
 fileList<<plate<<"  "<<mjd<<"  "<<fiber<<"  "<<setprecision(9)<<ra<<" "<<dec<<setprecision(6)
	    <<"  "<<z<<"  "<<umag<<"  "<<gmag<<"  "<<rmag<<"  "<<imag<<"  "<<zmag
	    <<"  "<<iBAL<<"  "<<iDLA<<"  "<<NbObs;
    if(NbObs>1){
      for(int iqso=0; iqso<NbObs-1; iqso++) 
	fileList<<"  "<<plateReObs[iqso]<<"  "<<mjdReObs[iqso]<<"  "<<fiberReObs[iqso];
    }
    fileList<<endl;
	'''
	
	
	return

def savePlateList_eBOSS():
	'''
	'''
	
	cat = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/platelist-v5_8_0.fits")[1].data
	
	## Data release 
	plate = [str(elem) for elem in cat["PLATE"] ]
	mjd   = [str(elem) for elem in cat["MJD"] ]
	chunk = [str(elem) for elem in cat["CHUNK"] ]
	dr    = [str(0) for i in range(0, len(cat)) ]

	file_name = "PlateList_eBOSS_v5_8_0.txt"
	file = open(file_name,'w')
	numpy.savetxt(file,zip(plate, mjd, cat["PLATEQUALITY"], chunk, dr),fmt='%s %s %s %s %s')
	file.close()
	
	
	return

def getPlateList():
	'''
	
	'''
	
	cat  = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/platelist.fits")[1].data
	cat2 = pyfits.open("/home/helion/Documents/Thèse/Data/Fits/platelist-v5_7_9.fits")[1].data
	##print pyfits.open("/home/helion/Documents/Thèse/Data/Fits/platelist-v5_7_9.fits")[1].header
	
	cat = cat[ cat["MJD"] <= 56658]
	
	plate = cat["PLATE"]
	#plate = numpy.append(plate,cat2["PLATE"])
	
	file_name = 'plate_before56658.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(plate),"%d")
	file.close()
	
	return

#plot_all_cat()
plot_all_cat_Forest()



