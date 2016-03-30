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
import decimal ## To set the precision of the double values
import cosmolopy.distance as cosmology
import matplotlib.pyplot as plt

import myTools
from myTools import Get_TProfile
from const_delta import *

def main():
	'''
	'''

	zKEY = 'Z'
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_ALL_TESTS.fits')[1].data
	cat = cat[ cat['Z']>1.7 ]
	cat = cat[ cat['Z']<7. ]
	cat2 = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR7_DR12_EBOSS.fits')[1].data
	cat2 = cat2[ cat2['Z']>1.7 ]
	cat2 = cat2[ cat2['Z']<7. ]
	cat3 = pyfits.open('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits',memmap=True)[1].data[:1000]
	cat3 = cat3[ cat3['Z_VI']>1.7 ]
	cat3 = cat3[ cat3['Z_VI']<7. ]

	print cat.size
	print cat2.size
	print cat3.size
	
	## Map
	plt.plot(cat['RA'],  cat['DEC'], linestyle="", marker="o",label='DR7+DR12: QSO')
	plt.plot(cat3['RA'],  cat3['DEC'], linestyle="", marker="o",label='DR7+DR12: Forest',color='red')
	plt.xlabel(r'$R.A. (\degree)$')
	plt.ylabel(r'$Dec. (\degree)$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	## Distribution redshift
	plt.hist(cat[zKEY], bins=50,label='DR7+DR12: QSO',histtype='step')
	plt.hist(cat2[zKEY], bins=50,label='DR7+DR12+eBOSS: QSO',histtype='step')
	plt.hist(cat3['Z_VI'], bins=50,label='DR7+DR12: Forest',histtype='step')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

	'''
	col_ra              = pyfits.Column(name='RA',  format='D', array=cat['RA'], unit='deg')
	col_de              = pyfits.Column(name='DEC', format='D', array=cat['DEC'], unit='deg')
	col_zz              = pyfits.Column(name='Z',   format='D', array=cat['Z_BEST'] )
	
	tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
	tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_3DHST.fits', clobber=True)
	'''
	return
def getQsoCatalogueDR7():
	'''
		Get the catalogue of qsos for DR7
	'''
	
	redMin = 0.
	redMax = 20.
	
	redShiftKey = "Z"
	
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/knownquasarstar.060910.fits')[1].data
	redShiftKey = "ZEM"
	
	return numpy.asarray(zip(cat["RA"],cat["DEC"],cat["ZEM"],cat['GMAG']))

	print "  --- DR7 ---"
	print " Size of cat = " + str(len(cat))
	
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(redMax)
	print " Size of cat after cut in Z = " + str(len(cat[ cat[redShiftKey] >= 2. ])) + " qso with z >= " + str(2.)
	
	file_name = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR7.txt'
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
	
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits')[1].data

	return numpy.asarray(zip(cat["RA"],cat["DEC"],cat["Z_VI"],cat['PSFMAG'][:,1]))
	print "  --- DR12 ---"
	print " Size of cat = " + str(len(cat))
	
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	print " Size of cat after cut in Z = " + str(len(cat[ cat[redShiftKey] >= 2. ])) + " qso with z >= " + str(2.)
	
	'''
	file_name = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR12.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat["RA"],cat["DEC"],cat[redShiftKey]))
	file.close()
	'''
	
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

	path = '/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_7_9.fits'
	redShiftKey = 'Z'
	cat = pyfits.open(path, memmap=True )[1].data

	'''
	print "  The size of the catalogue is           : ", cat.size
	## cat = cat[ cat[redShiftKey]>minRedshift__]
	## print "  We keep Z > " + str(minRedshift__) + "  , the size is      : ", cat.size
	## cat = cat[ cat[redShiftKey]<=maxRedshift__]
	## print "  We keep Z <= " + str(maxRedshift__) + "  , the size is      : ", cat.size
	cat = cat[ cat['CLASS'] == 'QSO' ]
	print '  We keep CLASS == QSO, the size is      : ', cat.size
	cat = cat[ cat['ZWARNING'] == 0 ]
	print '  We keep ZWARNING == 0, the size is      : ', cat.size
	cat = cat[ cat['Z_ERR'] > 0 ]
	print '  We keep Z_ERR > 0, the size is      : ', cat.size
	'''

	### BOSS_TARGET1
	selection_BOSS_TARGET1 = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
	bits_BOSS_TARGET1 = 0
	for el in selection_BOSS_TARGET1:
		bits_BOSS_TARGET1 += 2**el
	### ANCILLARY_TARGET1
	selection_ANCILLARY_TARGET1 = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,50,51,52,53,54,55,58,59]
	bits_ANCILLARY_TARGET1 = 0
	for el in selection_ANCILLARY_TARGET1:
		bits_ANCILLARY_TARGET1 += 2**el
	### ANCILLARY_TARGET2
	selection_ANCILLARY_TARGET2 = [0,1,2,3,4,5,7,8,9,10,13,14,15,24,25,26,27,31,32,33,53,54,55,56]
	bits_ANCILLARY_TARGET2 = 0
	for el in selection_ANCILLARY_TARGET2:
		bits_ANCILLARY_TARGET2 += 2**el
	### EBOSS_TARGET0
	selection_EBOSS_TARGET0 = [10,11,12,13,14,15,16,17,18,20,22,30,31,33,34,35,40]
	bits_EBOSS_TARGET0 = 0
	for el in selection_EBOSS_TARGET0:
		bits_EBOSS_TARGET0 += 2**el
	### EBOSS_TARGET1
	selection_EBOSS_TARGET1 = [9,10,11,12,13,14,15,16,17,18,30,31]
	bits_EBOSS_TARGET1 = 0
	for el in selection_EBOSS_TARGET1:
		bits_EBOSS_TARGET1 += 2**el
	### EBOSS_TARGET2
	selection_EBOSS_TARGET2 = [0,2,4,20,21,23,24,25,26,27,31]
	bits_EBOSS_TARGET2 = 0
	for el in selection_EBOSS_TARGET2:
		bits_EBOSS_TARGET2 += 2**el








	print cat.size
	select = ((cat['BOSS_TARGET1'] & bits_BOSS_TARGET1 )>0 )            \
		| ((cat['ANCILLARY_TARGET1'] & bits_ANCILLARY_TARGET1 )>0)  \
		| ((cat['ANCILLARY_TARGET2'] & bits_ANCILLARY_TARGET2 )>0 ) \
		| ((cat['EBOSS_TARGET0'] & bits_EBOSS_TARGET0 )>0 )         \
		| ((cat['EBOSS_TARGET1'] & bits_EBOSS_TARGET1 )>0 )         \
		| ((cat['EBOSS_TARGET2'] & bits_EBOSS_TARGET2 )>0 )
	cat = cat[select]
	print cat.size

	## cat = cat[ cat['ZWARNING'] == 0 ]
	print '  We keep ZWARNING == 0, the size is      : ', cat[ cat['ZWARNING'] != 0 ].size
	## cat = cat[ cat['Z_ERR'] > 0 ]
	print '  We keep Z_ERR > 0, the size is      : ', cat[ cat['Z_ERR'] <= 0 ].size
	print ' Z < 0  ', cat[ (cat[redShiftKey]<0.) ].size
	print ' Z < 0.1  ', cat[ (cat[redShiftKey]<0.1) ].size
	print ' Z > 1.8  ', cat[ (cat[redShiftKey]>1.8) ].size


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


	col_ra              = pyfits.Column(name='RA',  format='D', array=cat['RA'], unit='deg')
	col_de              = pyfits.Column(name='DEC', format='D', array=cat['DEC'], unit='deg')
	col_zz              = pyfits.Column(name='Z',   format='D', array=cat[redShiftKey] )
	
	tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
	tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_EBOSS.fits', clobber=True)

	
	return
def plot_all_cat():
	'''

	'''

	redMin = 1.7
	redMax = 6.
	
	### All known qsos
	print "  --- All known ---"
	redShiftKey = 'ZEM'
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/knownquasarstar.060910.fits')[1].data
	print " Size of cat = " + str(cat.size)
	
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(cat.size) + " qso with z > " + str(redMin)
	
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(cat.size) + " qso with z < " + str(redMax)
	cat = cat[ cat['DEC'] > -20. ]
	print cat.size
	cat = cat[ cat['DEC'] < 80. ]
	print cat.size

	dr7_ra  = cat['RA']
	dr7_dec = cat['DEC']
	dr7_z   = cat['ZEM']
	
	'''
	file_name = 'QSO_DR7.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat['RA'],cat['DEC'],cat[redShiftKey]))
	file.close()
	'''
	
	### DR12Q
	print "  --- DR12 ---"
	redShiftKey = 'Z_VI'
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits')[1].data
	print " Size of cat = " + str(len(cat))
	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(redMax)
	
	dr12_ra  = cat['RA']
	dr12_dec = cat['DEC']
	dr12_z   = cat['Z_VI']
	
	'''
	file_name = 'QSO_DR12.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat['RA'],cat['DEC'],cat[redShiftKey]))
	file.close()
	'''
	
	plate_list_DR12 = []
	for i in cat['PLATE']:
		if i not in plate_list_DR12:
			plate_list_DR12.append(i)
	plate_list_DR12 = numpy.sort(plate_list_DR12)
	print '  Size = ', plate_list_DR12.size
	
	print '  density DR12 = ', cat.size/(plate_list_DR12.size*5.)
	
	### eBOSS
	print "  --- eBOSS ---"
	redShiftKey = 'Z'
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_7_9.fits')[1].data

	cat = cat[ cat[redShiftKey] > redMin ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(redMin)
	cat = cat[ cat[redShiftKey] < redMax ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(redMax)

	### BOSS_TARGET1
	selection_BOSS_TARGET1 = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
	bits_BOSS_TARGET1 = 0
	for el in selection_BOSS_TARGET1:
		bits_BOSS_TARGET1 += 2**el
	### ANCILLARY_TARGET1
	selection_ANCILLARY_TARGET1 = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,50,51,52,53,54,55,58,59]
	bits_ANCILLARY_TARGET1 = 0
	for el in selection_ANCILLARY_TARGET1:
		bits_ANCILLARY_TARGET1 += 2**el
	### ANCILLARY_TARGET2
	selection_ANCILLARY_TARGET2 = [0,1,2,3,4,5,7,8,9,10,13,14,15,24,25,26,27,31,32,33,53,54,55,56]
	bits_ANCILLARY_TARGET2 = 0
	for el in selection_ANCILLARY_TARGET2:
		bits_ANCILLARY_TARGET2 += 2**el
	### EBOSS_TARGET0
	selection_EBOSS_TARGET0 = [10,11,12,13,14,15,16,17,18,20,22,30,31,33,34,35,40]
	bits_EBOSS_TARGET0 = 0
	for el in selection_EBOSS_TARGET0:
		bits_EBOSS_TARGET0 += 2**el
	### EBOSS_TARGET1
	selection_EBOSS_TARGET1 = [9,10,11,12,13,14,15,16,17,18,30,31]
	bits_EBOSS_TARGET1 = 0
	for el in selection_EBOSS_TARGET1:
		bits_EBOSS_TARGET1 += 2**el
	### EBOSS_TARGET2
	selection_EBOSS_TARGET2 = [0,2,4,20,21,23,24,25,26,27,31]
	bits_EBOSS_TARGET2 = 0
	for el in selection_EBOSS_TARGET2:
		bits_EBOSS_TARGET2 += 2**el
	print cat.size
	select = ((cat['BOSS_TARGET1'] & bits_BOSS_TARGET1 )>0 )            \
		| ((cat['ANCILLARY_TARGET1'] & bits_ANCILLARY_TARGET1 )>0)  \
		| ((cat['ANCILLARY_TARGET2'] & bits_ANCILLARY_TARGET2 )>0 ) \
		| ((cat['EBOSS_TARGET0'] & bits_EBOSS_TARGET0 )>0 )         \
		| ((cat['EBOSS_TARGET1'] & bits_EBOSS_TARGET1 )>0 )         \
		| ((cat['EBOSS_TARGET2'] & bits_EBOSS_TARGET2 )>0 )
	cat = cat[select]
	print cat.size
	
	plate_list_eBOSS = []
	for i in cat['PLATE']:
		if i not in plate_list_eBOSS:
			plate_list_eBOSS.append(i)
	plate_list_eBOSS = numpy.sort(plate_list_eBOSS)
	print '  Size = ', plate_list_eBOSS.size
	
	print '  density eBOSS = ', cat.size/(plate_list_eBOSS.size*5.)
	
	eBOSS_ra  = cat['RA']
	eBOSS_dec = cat['DEC']
	eBOSS_z   = cat['Z']
	
	
	file_name = 'QSO_eBOSS.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat['RA'],cat['DEC'],cat[redShiftKey]))
	file.close()
	'''
	
	
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
	#plt.hist(dr7_z,                                            bins=100, histtype='step', label=r'$Known \, before \, DR12$', linewidth=2.)
	plt.hist(dr12_z,                                           bins=100, histtype='step', label=r'$DR12Q$', linewidth=2.,weights=numpy.ones(dr12_z.size)/plate_list_DR12.size)
	plt.hist(eBOSS_z,                                          bins=100, histtype='step', label=r'$eBOSS$', linewidth=2.,weights=numpy.ones(eBOSS_z.size)/plate_list_eBOSS.size)
	
	plt.xlabel(r'$z_{QSO}$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	'''

	dr12_z_h = myTools.GetHisto(dr12_z,numpy.arange(1.6,6.1,0.1))
	dr12_z_h[:,1] /= plate_list_DR12.size*5.
	print numpy.sum( dr12_z_h[:,1] )
	print '  area = ', plate_list_DR12.size*5.
	plt.plot(dr12_z_h[:,0],dr12_z_h[:,1], label=r'$DR12$', drawstyle='steps')

	eBOSS_z_h = myTools.GetHisto(eBOSS_z,numpy.arange(1.6,6.1,0.1))
	eBOSS_z_h[:,1] /= plate_list_eBOSS.size*5.
	print numpy.sum( eBOSS_z_h[:,1] )
	print '  area = ', plate_list_eBOSS.size*5.
	plt.plot(eBOSS_z_h[:,0],eBOSS_z_h[:,1], label=r'$eBOSS$', drawstyle='steps',color='red')

	plt.xlabel(r'$z_{QSO}$', fontsize=40)
	plt.ylabel(r'$\#.deg^{-2}$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()

	
def plot_all_cat_Forest():
	'''

	'''
	
	### DR12Q
	#path = "/home/gpfs/manip/mnt/bao/hdumasde/Lists/DR12Q_v2_10.fits"
	path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits"		
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
	
	file_name = 'QSO_DR12.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat['RA'],cat['DEC'],cat['Z_VI']))
	file.close()
	
	
	plate_list_DR12 = []
	for i in cat['PLATE']:
		if i not in plate_list_DR12:
			plate_list_DR12.append(i)
	plate_list_DR12 = numpy.sort(plate_list_DR12)
	print '  Size = ', plate_list_DR12.size
	
	print '  density DR12 = ', cat.size/(plate_list_DR12.size*5.)
	
	'''
	### eBOSS
	print "  --- eBOSS ---"
	redShiftKey = 'Z'
	#path = "/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_8_0.fits"
	path = "/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_7_9.fits"
	cat = pyfits.open(path)[1].data
	
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
	
	file_name = 'QSO_eBOSS.txt'
	file = open(file_name,'w')
	numpy.savetxt(file,zip(cat['RA'],cat['DEC'],cat['Z']))
	file.close()
	
	
	plate_list_eBOSS = []
	for i in cat['PLATE']:
		if i not in plate_list_eBOSS:
			plate_list_eBOSS.append(i)
	plate_list_eBOSS = numpy.sort(plate_list_eBOSS)
	print '  Size = ', plate_list_eBOSS.size
	
	print '  density eBOSS = ', cat.size/(plate_list_eBOSS.size*5.)
	
	
	eBOSS_ra  = cat['RA']
	eBOSS_dec = cat['DEC']
	eBOSS_z   = cat['Z']
	'''

	print "  --- eBOSS ---"
	redShiftKey = 'Z'
	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_7_9.fits')[1].data

	cat = cat[ cat[redShiftKey] > minRedshift__ ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z > " + str(minRedshift__)
	cat = cat[ cat[redShiftKey] < maxRedshift__ ]
	print " Size of cat after cut in Z = " + str(len(cat)) + " qso with z < " + str(maxRedshift__)

	### BOSS_TARGET1
	selection_BOSS_TARGET1 = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
	bits_BOSS_TARGET1 = 0
	for el in selection_BOSS_TARGET1:
		bits_BOSS_TARGET1 += 2**el
	### ANCILLARY_TARGET1
	selection_ANCILLARY_TARGET1 = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,50,51,52,53,54,55,58,59]
	bits_ANCILLARY_TARGET1 = 0
	for el in selection_ANCILLARY_TARGET1:
		bits_ANCILLARY_TARGET1 += 2**el
	### ANCILLARY_TARGET2
	selection_ANCILLARY_TARGET2 = [0,1,2,3,4,5,7,8,9,10,13,14,15,24,25,26,27,31,32,33,53,54,55,56]
	bits_ANCILLARY_TARGET2 = 0
	for el in selection_ANCILLARY_TARGET2:
		bits_ANCILLARY_TARGET2 += 2**el
	### EBOSS_TARGET0
	selection_EBOSS_TARGET0 = [10,11,12,13,14,15,16,17,18,20,22,30,31,33,34,35,40]
	bits_EBOSS_TARGET0 = 0
	for el in selection_EBOSS_TARGET0:
		bits_EBOSS_TARGET0 += 2**el
	### EBOSS_TARGET1
	selection_EBOSS_TARGET1 = [9,10,11,12,13,14,15,16,17,18,30,31]
	bits_EBOSS_TARGET1 = 0
	for el in selection_EBOSS_TARGET1:
		bits_EBOSS_TARGET1 += 2**el
	### EBOSS_TARGET2
	selection_EBOSS_TARGET2 = [0,2,4,20,21,23,24,25,26,27,31]
	bits_EBOSS_TARGET2 = 0
	for el in selection_EBOSS_TARGET2:
		bits_EBOSS_TARGET2 += 2**el
	print cat.size
	select = ((cat['BOSS_TARGET1'] & bits_BOSS_TARGET1 )>0 )            \
		| ((cat['ANCILLARY_TARGET1'] & bits_ANCILLARY_TARGET1 )>0)  \
		| ((cat['ANCILLARY_TARGET2'] & bits_ANCILLARY_TARGET2 )>0 ) \
		| ((cat['EBOSS_TARGET0'] & bits_EBOSS_TARGET0 )>0 )         \
		| ((cat['EBOSS_TARGET1'] & bits_EBOSS_TARGET1 )>0 )         \
		| ((cat['EBOSS_TARGET2'] & bits_EBOSS_TARGET2 )>0 )
	cat = cat[select]
	print cat.size
	
	plate_list_eBOSS = []
	for i in cat['PLATE']:
		if i not in plate_list_eBOSS:
			plate_list_eBOSS.append(i)
	plate_list_eBOSS = numpy.sort(plate_list_eBOSS)
	print '  Size = ', plate_list_eBOSS.size
	
	print '  density eBOSS = ', cat.size/(plate_list_eBOSS.size*5.)
	
	eBOSS_ra  = cat['RA']
	eBOSS_dec = cat['DEC']
	eBOSS_z   = cat['Z']

	#allQSO = numpy.loadtxt("/home/helion/Documents/Thèse/Data/List/QSO_DR12_DR7_Helion_LYA.txt")
	
	## Map
	plt.plot(dr12_ra,dr12_dec,   linestyle="", marker="o", label=r'$DR12Q$')
	plt.plot(eBOSS_ra,eBOSS_dec, linestyle="", marker="o", label=r'$eBOSS$')
	plt.xlabel(r'$R. A. \, [^\circ]$', fontsize=40)
	plt.ylabel(r'$Dec. \, [^\circ]$', fontsize=40)
	plt.xlim([0.,360.])
	myTools.deal_with_plot(False,False,True)
	plt.show()
	## Distribution redshift
	plt.hist(dr12_z,                       bins=100, histtype='step', label=r'$DR12Q$', linewidth=2.,weights=numpy.ones(dr12_z.size)/plate_list_DR12.size)
	plt.hist(eBOSS_z,                      bins=100, histtype='step', label=r'$eBOSS$', linewidth=2.,weights=numpy.ones(eBOSS_z.size)/plate_list_eBOSS.size)
	#plt.hist(allQSO[:,2], bins=100, histtype='step', label=r'$All$', linewidth=2.,weights=numpy.ones(allQSO[:,2].size)/plate_list_DR12.size)
	plt.xlabel(r'$z_{forest}$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	plt.xlim([1.6,5.])
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	'''
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
	'''


	dr12_z_h = myTools.GetHisto(dr12_z,numpy.arange(minRedshift__-0.1,6.1,0.1))
	dr12_z_h[:,1] /= plate_list_DR12.size*5.
	print dr12_z_h
	print numpy.sum( dr12_z_h[:,1] )
	print '  area = ', plate_list_DR12.size*5.
	plt.plot(dr12_z_h[:,0],dr12_z_h[:,1], label=r'$DR12$', drawstyle='steps')

	eBOSS_z_h = myTools.GetHisto(eBOSS_z,numpy.arange(minRedshift__-0.1,6.1,0.1))
	eBOSS_z_h[:,1] /= plate_list_eBOSS.size*5.
	print eBOSS_z_h
	print numpy.sum( eBOSS_z_h[:,1] )
	print '  area = ', plate_list_eBOSS.size*5.
	plt.plot(eBOSS_z_h[:,0],eBOSS_z_h[:,1], label=r'$eBOSS$', drawstyle='steps',color='red')

	plt.xlabel(r'$z_{forest}$', fontsize=40)
	plt.ylabel(r'$\#.deg^{-2}$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()




	return
def bias_evolution():
	'''
	Loi Croom et al.

	'''

	### Data points
	#data = [2.39685630578, 3.246e+00, 0., 1.863e-01]
	data = [2.39686216554, 3.407847, 0., 0.182078]
	### Data Pierre
	dataPierre = numpy.asarray( [ [1.0573106, 1.978919009022121, 0.1326168218002746], [1.3519373, 2.3978591946019843, 0.11460107023841509], [1.651769, 2.613084388995206, 0.12670393048129994], [1.9889176, 2.9103797886104594, 0.15839821556293388] ] )
	### Simu points
	simu = numpy.zeros(shape=(4,100))
	simu[0,:] = 2.40673056001
	'''
	param = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Annalyse_BAOFIT/param.npy')
	print param[:,0,0]
	simu[2,:] = param[5,0,:]
	simu[3,:] = param[5,1,:]
	meanSimu = numpy.asarray([ numpy.average(simu[0,:]),  numpy.average(simu[2,:],weights=numpy.power(simu[3,:],-2.)),  numpy.power(numpy.sum(numpy.power(simu[3,:],-2.)),-0.5) ])
	print numpy.average(simu[2,:],weights=numpy.power(simu[3,:],-2.))
	print 3.246e+00
	'''
	### Croom et al. Law
	def croom(x):
		return 0.53 + 0.289*(1.+x)**2
	xxx = numpy.arange(0.,10.,0.01)
	yyy = numpy.apply_along_axis(croom, -1, xxx)
	print '  Croom et al. interpolation : z=2.4 ', croom(2.4)

	plt.plot(xxx,yyy, linestyle='dashed', label=r'$Croom \, et \, al.$')
	plt.errorbar(simu[0,:],simu[2,:],yerr=simu[3,:],fmt='o',color='green',label=r'$Simu$')
	#plt.errorbar(meanSimu[0],meanSimu[1],yerr=meanSimu[2], fmt='o',color='blue',label=r'$<Simu>$')
	plt.errorbar(data[0],data[1],yerr=data[3],fmt='o',color='red', label=r'$Data$')
	plt.errorbar(dataPierre[:,0],dataPierre[:,1],yerr=dataPierre[:,2],fmt='o',color='black',label=r'$data \, Pierre$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

plot_all_cat_Forest()
#plot_all_cat()


bias_evolution()



























