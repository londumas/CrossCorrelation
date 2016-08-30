# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#
#  /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/correlation_3D.py
#

### Python lib
import subprocess
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit
import sys
import copy
import time

### Perso lib
import const_delta
import myTools
import correlation_3D
import correlation_3D_Q
import annalyse_BAOFIT
import annalyse_many_BAOFIT
import CAMB

import matplotlib.pyplot as plt




def get_number_in_files_not_found():

	'''
	grep 'File not found:' echo.o* > not_found.txt
	'''

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/CIV/DR14/DR14_primery/Log/not_found.txt'
	rawPat = '/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_10_0/'
	lenPath = len(rawPat)

	plate = []
	mjd   = []
	fiber = []
	ra = []
	de = []
	zz = []

	with open(path, 'r') as f:
		for line in f:
			line = line.split()
			ra    += [ line[-3] ]
			de    += [ line[-2] ]
			zz    += [ line[-1] ]
			line   = line[4][lenPath:].split('/')[1].split('-')
			plate += [ line[1] ]
			mjd   += [ line[2] ]
			fiber += [ line[3].split('.')[0] ]

	plate = numpy.array(plate).astype('int')
	mjd   = numpy.array(mjd).astype('int')
	fiber = numpy.array(fiber).astype('int')
	ra = numpy.array(ra).astype('float')
	de = numpy.array(de).astype('float')
	zz = numpy.array(zz).astype('float')
	print '  size = ', plate.size
	
	
	cut = mjd>55000
	plate = plate[cut]
	mjd = mjd[cut]
	fiber = fiber[cut]
	print plate.size
	
	
	for i in range(plate.size):

		'''
		### Make folder
		path_to_make = '/home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_10_0/'+str(plate[i])+'/'
		print path_to_make
		subprocess.call( 'mkdir \"'+path_to_make+'\"' , shell=True)
		'''

		
		fileName = str(plate[i]) + '/spec-'+str(plate[i])+'-'+str(mjd[i])+'-'+str(fiber[i]).zfill(4)+'.fits'
		string_to_download = 'wget -nv -N --user=sdss --password=2.5-meters https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_0/spectra/'+fileName
		string_to_download += ' -O /home/gpfs/manip/mnt0607/bao/Spectra/SpectraV5_10_0/' + fileName
		#string_to_download = 'rsync -aLrz --password-file ../password-eboss.txt rsync://sdss@dtn01.sdss.utah.edu/sas/ebosswork/eboss/spectro/redux/v5_10_0/spectra/'+fileName + ' '+ str(plate[i]) + '/'
		
		string_to_download = 'ls ' + fileName
		#tmp_command = "clubatch \"time ; hostname ; " + string_to_download + "\""
		tmp_command = string_to_download
		print tmp_command

		subprocess.call(tmp_command, shell=True)
		#time.sleep(0.001)
		#myTools.isReadyForNewJobs(200, 1000,'time')
		

	
	###
	plt.hist(plate,bins=100)
	plt.show()
	###
	plt.hist(mjd,bins=100)
	plt.show()
	###
	plt.hist(fiber,bins=100)
	plt.show()
	###
	plt.hist(zz,bins=100)
	plt.show()
	###
	plt.errorbar(ra,de,fmt='o')
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.show()

	return
def get_number_in_files_pixels():

	'''
	grep 'less or equal to 50:' echo.o* > less.txt
	'''

	path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/DR14/DR14_primery/Log/less.txt'

	step     = []
	nb_pixel = []
	z_forest = []
	with open(path, 'r') as f:
		for line in f:
			line = line.split()
			step     += [ line[1][0] ]
			nb_pixel += [ line[7] ]
			z_forest += [ line[-1] ]

	step     = numpy.array(step).astype('int')
	nb_pixel = numpy.array(nb_pixel).astype('int')
	z_forest = numpy.array(z_forest).astype('float')
	print '  total cut ', nb_pixel.size
	z_forest = z_forest[ (nb_pixel!=0) ]
	nb_pixel = nb_pixel[ (nb_pixel!=0) ]
	print '  cut with nb pixel >0 ', nb_pixel.size

	###
	plt.hist(step,bins=2, range=(-0.5,1.5))
	plt.xlabel('step')
	plt.show()
	###
	plt.hist(nb_pixel,bins=50)
	plt.xlabel('nb pixel')
	plt.show()
	###
	plt.hist(z_forest,bins=100)
	plt.xlabel('z forest')
	plt.show()
	###
	plt.errorbar(z_forest,nb_pixel,fmt='o')
	plt.xlabel('z forest')
	plt.ylabel('nb pixel')
	plt.show()

	return
def print_low():

	

	return

get_number_in_files_not_found()
#get_number_in_files_pixels()
















