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
def get_QSO_cat_DR14():

        path = '/home/gpfs/manip/mnt0607/bao/Spectra/DR14Q_v1_0.fits'
        cat = pyfits.open(path, memmap=True)[1].data

        print cat.size
	cat = cat[(cat['RA']!=0.)]
	cat = cat[(cat['DEC']!=0.)]
	cat = cat[(cat['Z']>0.)]

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
        plt.hist(cat['Z'], bins=50)
        plt.xlabel("Z")
        plt.ylabel("#")
        plt.title("BOSS DR14")
        plt.show()

	### Save        
        col_ra              = pyfits.Column(name='RA',  format='D', array=cat['RA'], unit='deg')
        col_de              = pyfits.Column(name='DEC', format='D', array=cat['DEC'], unit='deg')
        col_zz              = pyfits.Column(name='Z',   format='D', array=cat['Z'])
        tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
        tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR14_v1_0.fits', clobber=True)

        return
def getQsoCatalogueEBOSS():
	'''
		Get the catalogue of qsos
	'''

	path = '/home/gpfs/manip/mnt0607/bao/Spectra/spAll-v5_9_1.fits'
	redShiftKey = 'Z'
	cat = pyfits.open(path, memmap=True )[1].data

	print "  The size of the catalogue is           : ", cat.size
	cat = cat[ cat['CLASS'] == 'QSO' ]
	print '  We keep CLASS == QSO, the size is      : ', cat.size
	cat = cat[ cat['ZWARNING'] == 0 ]
	print '  We keep ZWARNING == 0, the size is      : ', cat.size
	cat = cat[ cat['Z_ERR'] > 0 ]
	print '  We keep Z_ERR > 0, the size is      : ', cat.size

	print  cat[ cat['Z'] > 1.7 ].size


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
	tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_EBOSS_updated_2016_05_24.fits', clobber=True)
	
	
	return
def getQsoCatalogueAllQSO():

	cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR7_DR12_EBOSS_2016_01_08.fits', memmap=True )[1].data

	print cat.size
	cat = cat[ cat['RA']!=0. ]
	print cat.size

	col_ra              = pyfits.Column(name='RA',  format='D', array=cat['RA'], unit='deg')
	col_de              = pyfits.Column(name='DEC', format='D', array=cat['DEC'], unit='deg')
	col_zz              = pyfits.Column(name='Z',   format='D', array=cat['Z'] )

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
	plt.hist(cat['Z'][ cat['Z']>1.7 ], bins=50)
	plt.xlabel("Z")
	plt.ylabel("#")
	plt.title("BOSS DR12")
	plt.show()
	
	tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
	tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR7_DR12_EBOSS_2016_01_08.fits', clobber=True)
def getQsoCatalogueAllObjects():
	'''

	'''

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/'
	listPAth = [path+'QSO_DR7_DR12_EBOSS_2016_01_08.fits',
			path+'DLA_all.fits',
			path+'all_Britt.fits',
			path+'VIPERS.fits',
			path+'QSO_3DHST.fits',
			path+'LOWZ_all.fits',
			path+'CMASS_all.fits']
	name = ['QSO','DLA','Britt','VIPERS','3DHST','LOWZ','CMASS']

	## Distribution redshift
	for i in numpy.arange(len(listPAth)):
		cat = pyfits.open(listPAth[i], memmap=True )[1].data
		cat = cat[ (cat['Z']>0.1) ]
		cat = cat[ (cat['Z']<7.) ]
		if (cat.size==0): continue
		plt.hist(cat['Z'], bins=100,histtype='step',label=name[i])

	plt.xlabel("Z")
	plt.ylabel("#")
	myTools.deal_with_plot(False,False,True)
	plt.show()

	### Merge everyThing
	cat = pyfits.open(listPAth[0], memmap=True )[1].data
	ra = cat['RA']
	de = cat['DEC']
	zz = cat['Z']
	for i in numpy.arange(1,len(listPAth)):
		cat = pyfits.open(listPAth[i], memmap=True )[1].data
		ra = numpy.append(ra, cat['RA'])
		de = numpy.append(de, cat['DEC'])
		zz = numpy.append(zz, cat['Z'])

	## Map
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))	
	plt.grid()
	plt.plot(ra, de, linestyle="", marker="o")
	plt.xlabel("Right Ascension (degree)")
	plt.ylabel("Declination (degree)")
	plt.show()
	## Distribution redshift
	plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
	plt.grid()
	plt.hist(zz, bins=200)
	plt.xlabel("Z")
	plt.ylabel("#")
	plt.show()

	### Save	
	col_ra              = pyfits.Column(name='RA',  format='D', array=ra, unit='deg')
	col_de              = pyfits.Column(name='DEC', format='D', array=de, unit='deg')
	col_zz              = pyfits.Column(name='Z',   format='D', array=zz)
	tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
	tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/ALL_EVERY_OBJECTS_2016_01_08.fits', clobber=True)
def create_fits_qso(cat_path, sizeMax):

	### Create a list of QSOs
	col_ra              = pyfits.Column(name='RA',  format='D', array=numpy.zeros(sizeMax), unit='deg')
	col_de              = pyfits.Column(name='DEC',  format='D', array=numpy.zeros(sizeMax), unit='deg')
	col_zz              = pyfits.Column(name='Z',  format='D', array=numpy.zeros(sizeMax), unit='0 (redshift)')

	tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
	tbhdu.writeto(cat_path, clobber=True)

	return
def remove_empty_cell_in_cat(cat_path):

	cat = pyfits.open(cat_path, memmap=True)[1].data
	print '  nb of QSOs before = ', cat.size
	cut = numpy.invert(numpy.logical_and( cat['RA']==0., numpy.logical_and( cat['DEC']==0.,cat['Z'] == 0.)))
	cat = cat[ (cut) ]
	print '  nb of QSOs after  = ', cat.size
	nbQSO = cat.size

	plt.errorbar(cat['RA'],cat['DEC'],fmt='o')
	plt.show()
	plt.hist(cat['Z'][ (cat['Z']>=0.01) ],bins=1000)
	plt.show()
	
	pyfits.writeto(cat_path, cat, clobber=True)
	
	return
def plot_cat():

	path = '/home/gpfs/manip/mnt0607/bao/Spectra/DR14Q_v1_0.fits'
	cat = pyfits.open(path, memmap=True)[1].data

	print cat.size
	print
	print cat[ (cat['RA']<=0.) ].size
	print cat[ (cat['RA']>=360.) ].size
	print cat[ (cat['DEC']<=-90.) ].size
	print cat[ (cat['DEC']>=90.) ].size
	print cat[ (cat['RA']==0.) ].size
	print cat[ (cat['DEC']==0.) ].size
	print cat[ (cat['Z']<=0.) ].size
	print
	cat = cat[ cat['RA']==cat['DEC'] ]
	print cat.size
	print cat['Z']
	print cat['MJD']
	print cat['PLATE']

	### Geometry
	plt.errorbar(cat['RA'],cat['DEC'],fmt='o')
        plt.show()
	### redshift
	plt.hist(cat['Z'],bins=100)
	plt.show()


	return
def compare_DR12_DR14():


	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits'
	cat_DR12 = pyfits.open(path, memmap=True)[1].data
	path = '/home/gpfs/manip/mnt0607/bao/Spectra/DR14Q_v1_0.fits'
	cat_DR14 = pyfits.open(path, memmap=True)[1].data
	print cat_DR12.size
	print cat_DR14.size

	#cat_DR12 = cat_DR12[ cat_DR12['BOSS_TARGET1']==1 ]
	#cat_DR14 = cat_DR14[ cat_DR14['BOSS_TARGET1']==1 ]
	cat_DR12 = cat_DR12[ cat_DR12['EBOSS_TARGET0']==0 ]
	cat_DR14 = cat_DR14[ cat_DR14['EBOSS_TARGET0']==0 ]
	print cat_DR12.size
	print cat_DR14.size

	
	#print cat_DR12['PLATE']
	PLATE = 6173
	cat_DR12 = cat_DR12[ cat_DR12['PLATE']==PLATE ]
	cat_DR14 = cat_DR14[ cat_DR14['PLATE']==PLATE ]
	print cat_DR12.size
	print cat_DR14.size
	cat_DR12 = cat_DR12[1]
	cat_DR14 = cat_DR14[0]
	print
	print cat_DR12
	print
	print cat_DR14
	print
	print cat_DR14['RA']
	print cat_DR12['RA']
	print
	print cat_DR14['DEC']
	print cat_DR12['DEC']
	print
	

	### Geometry
	plt.errorbar(cat_DR14['RA'],cat_DR14['DEC'],fmt='o')
	plt.errorbar(cat_DR12['RA'],cat_DR12['DEC'],fmt='o')
        plt.show()
	### redshift
	plt.hist(cat_DR14['Z'],bins=100,histtype='step')
	plt.hist(cat_DR12['Z_VI'],bins=100,histtype='step')
	plt.show()



	return
def get_DLA():

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DLA_all.fits'
	cat = pyfits.open(path, memmap=True)[1].data

	print cat.size
	print

	### Geometry
	plt.errorbar(cat['RA'],cat['DEC'],fmt='o')
        plt.show()
	### redshift
	plt.hist(cat['Z'],bins=100)
	plt.show()

	return


get_DLA()
#get_QSO_cat_DR14()
#plot_cat()
#compare_DR12_DR14()

#a = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/PDF/DF_3e-4.fits"
#cat = pyfits.open(a, memmap=True)[1].data

#print cat['DABS'][0]
#myTools.plot2D(numpy.log10(cat['DABS'][0]))

"""
### Create a catalogue:
cat_path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR7_DR12_V5_10_0.fits'
sizeMax  = 5000000
#create_fits_qso(cat_path, sizeMax)
remove_empty_cell_in_cat(cat_path)

#getQsoCatalogueAllObjects()

#getQsoCatalogueEBOSS()
#main()
"""

'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/'
cat = pyfits.open(path + 'QSO_DR7_DR12_EBOSS.fits', memmap=True )[1].data
ra = cat['RA']
de = cat['DEC']
zz = cat['Z']
cat = pyfits.open(path + 'DLA_all.fits', memmap=True )[1].data
ra = numpy.append( ra, cat['RA'])
de = numpy.append( de, cat['DEC'])
zz = numpy.append( zz, cat['Z'])
cat = pyfits.open(path + 'all_Britt.fits', memmap=True )[1].data
ra = numpy.append( ra, cat['RA'])
de = numpy.append( de, cat['DEC'])
zz = numpy.append( zz, cat['Z'])



print ra.size

## Map
plt.plot(ra,  de, linestyle="", marker="o")
plt.xlabel(r'$R.A. (\degree)$')
plt.ylabel(r'$Dec. (\degree)$')
myTools.deal_with_plot(False,False,True)
plt.show()
## Distribution redshift
plt.hist(zz[ numpy.logical_and( zz>1., zz<7. ) ], bins=50)
plt.xlabel(r'$z$')
plt.ylabel(r'$\#$')
myTools.deal_with_plot(False,False,True)
plt.show()


col_ra              = pyfits.Column(name='RA',  format='D', array=ra, unit='deg')
col_de              = pyfits.Column(name='DEC', format='D', array=de, unit='deg')
col_zz              = pyfits.Column(name='Z',   format='D', array=zz )
	
tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
tbhdu.writeto('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/ALL_OBJECTS.fits', clobber=True)
'''

'''
cat3 = pyfits.open('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits',memmap=True)[1].data
cat3 = cat3[ cat3['Z_VI']>1.7 ]
cat3 = cat3[ cat3['Z_VI']<7. ]
RA  = cat3['RA']
Dec = cat3['DEC']
plt.errorbar(RA,Dec, fmt='o')
plt.show()
'''

"""
import numpy as np
import ephem

org=0
projection='mollweide'

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection=projection, axisbg ='white')


cat = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_ALL_TESTS.fits')[1].data
cat = cat[ cat['Z']>1.7 ]
cat = cat[ cat['Z']<7. ]
RA  = cat['RA']
Dec = cat['DEC']
RA[ RA>180. ] -= 360.
ax.scatter(np.radians(RA),np.radians(Dec), label=r'$DR7+DR12: \, QSO$')  # convert degrees to radians


cat3 = pyfits.open('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits',memmap=True)[1].data
cat3 = cat3[ cat3['Z_VI']>2. ]
cat3 = cat3[ cat3['Z_VI']<7. ]
RA  = cat3['RA']
Dec = cat3['DEC']
RA[ RA>180. ] -= 360.
ax.scatter(np.radians(RA),np.radians(Dec), label=r'$DR12: \, Forest$', color='red')  # convert degrees to radians


cat4 = pyfits.open('/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits',memmap=True)[1].data
cat4 = cat4[ cat4['Z_VI']>1.7 ]
cat4 = cat4[ cat4['Z_VI']<7. ]
RA  = cat4['RA']
Dec = cat4['DEC']
RA[ RA>180. ] -= 360.
ax.scatter(np.radians(RA),np.radians(Dec), label=r'$eBOSS \, forest$', color='green')  # convert degrees to radian


lon_array = np.arange(0,360)
lat = 0.
eq_array = np.zeros((360,2))
for lon in lon_array:
    ga = ephem.Galactic(np.radians(lon), np.radians(lat))
    eq = ephem.Equatorial(ga)
    eq_array[lon] = np.degrees(eq.get())

RA  = eq_array[:,0]
Dec = eq_array[:,1]
RA[ RA>180 ] -= 360.
coord = numpy.zeros(shape=(RA.size,2))
sort = RA.argsort()
RA  = RA[sort]
Dec = Dec[sort]
ax.plot(np.radians(RA),np.radians(Dec), color='green', linestyle='dashed', linewidth=5)  # convert degrees to radians

plt.xlabel(r'$R.A. (\degree)$')
plt.ylabel(r'$Dec. (\degree)$')

myTools.deal_with_plot(False,False,True)

plt.legend(bbox_to_anchor=(0., 1., 1., 0.), fontsize=30, frameon=False, scatterpoints=1,ncol=1, loc=4)
plt.show()




"""











