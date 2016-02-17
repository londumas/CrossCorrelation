# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# created by Hélion du Mas des Bourboux
# <Helion.du-Mas-des-Bourboux@cea.fr>
# module add scipy/0.13.1_Python_2.7.6-gnu

import subprocess
import time
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy
import myTools
from const_delta import *

import subprocess
import sys
import os
import time
import numpy
import scipy
import astropy.io.fits as pyfits
from iminuit import Minuit
import matplotlib.pyplot as plt
import re
import cosmolopy.distance as cosmology



### Constants
sizeMax = 300000
sizeMaxForest = 300000
nbQSO__ = 238929
nbFor__ = 170000
nbPixel = 645   ###2148  ###645
ratioForestToQSO__    = 1.*nbFor__/nbQSO__;
pathToFolder = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/' ##noMockExpander



### --------------------------------------------------
###
### idex_pass: 
###		- 0: create files
###		- 1: fill files
###		- 2: remove empty fields
###
index_function = -1
index_pass = -1
index_function = int(sys.argv[1])
index_pass = int(sys.argv[2])
print index_function
print index_pass

def main():

	print '  index_pass = ', index_pass
	
	tmp_command = "echo \" \n ------ Start ------ \n \" " 
	subprocess.call(tmp_command, shell=True)

	## Create the folders
	#####################
	subprocess.call('mkdir ' +pathToFolder+'Results/', shell=True)

	for i in range(0,1):		

		path = pathToFolder + 'Box_00' + str(i) + '/'

		## Create the folders
		#####################
		subprocess.call('mkdir ' +path, shell=True)

		for j in range(0,1):

			print i, j
			path = pathToFolder+ 'Box_00' + str(i) + '/Simu_00'+str(j) + '/'

			if (index_pass==0):
				subprocess.call('mkdir ' + path, shell=True)
				subprocess.call('mkdir ' + path + 'Raw', shell=True)
				subprocess.call('mkdir ' + path + 'Data', shell=True)
				subprocess.call('mkdir ' + path + 'Run', shell=True)
				subprocess.call('mkdir ' + path + 'Results', shell=True)
				subprocess.call('mkdir ' + path + 'Results_NicolasDistortion', shell=True)
				subprocess.call('mkdir ' + path + '/Results/BaoFit_q_f__LYA__QSO', shell=True)
				subprocess.call('mkdir ' + path + '/Results_NicolasDistortion/BaoFit_q_f__LYA__QSO', shell=True)

				if (i==0 and j==0):
					tbhduQSO    = create_fits_qso(sizeMax)
					tbhduQSO.writeto(path + 'Data/QSO_withRSD.fits', clobber=True)
					tbhduForest = create_fits_forest(sizeMaxForest)
					tbhduForest.writeto(path + 'Data/delta.fits', clobber=True)
				else:
					command = 'clubatch cp ' + pathToFolder + 'Box_000/Simu_000/Data/delta.fits ' + path + 'Data/delta.fits'
					subprocess.call(command, shell=True)
					command = 'clubatch cp ' + pathToFolder + 'Box_000/Simu_000/Data/QSO_withRSD.fits ' + path + 'Data/QSO_withRSD.fits'
					subprocess.call(command, shell=True)
					myTools.isReadyForNewJobs(10, 430,'cp')
					time.sleep(20)

			elif (index_pass==1):

				'''
				command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Mock_JMLG/ReadFits/bin/main.exe " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)
				myTools.isReadyForNewJobs(10, 430,'echo')
				time.sleep(20)
				'''

			elif (index_pass==2):
			
				cat = pyfits.open(path + 'Data/QSO_withRSD.fits', memmap=True)[1].data
				print cat.size
				cat = cat[ (cat['Z'] != 0.) ]
				print cat.size
				nbQSO = cat.size
				pyfits.writeto(path + 'Data/QSO_withRSD.fits', cat, clobber=True)
				print nbQSO
			
			
				### Remove useless lines in Forest
				cat = pyfits.open(path + 'Data/delta.fits', memmap=True)[1].data[:nbQSO]
				print cat.size
				'''
				tmp_idx  = numpy.arange(cat.size)
				tmp_bool = (cat['Z_VI'] != 0.)
				tmp_idx = tmp_idx[ (tmp_bool) ]
				last = tmp_idx[-1]
				print last
				cat = cat[:last+1]
				print cat.size
				nbFor = cat.size
				print nbFor-int(nbQSO*ratioForestToQSO__)
				rand = numpy.random.choice(nbFor, nbFor-int(nbQSO*ratioForestToQSO__), replace=False)
				cat['Z_VI'][rand] = 0.

				'''
				cat = cat[ (cat['Z_VI'] != 0.) ]
				print cat.size
				pyfits.writeto(path + 'Data/delta.fits', cat, clobber=True)
			

			'''
			### Get the data to 'good' files
			#####################
			#command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Mock_JMLG/launchjob.sh " + str(i) + ' ' + str(j) + "\""
			command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Mock_JMLG/ReadFits/bin/main.exe " + str(i) + ' ' + str(j) + "\""
			print command
			subprocess.call(command, shell=True)
			myTools.isReadyForNewJobs(10, 430)
			time.sleep(30)
			'''





	print
	print " ------ End ------"
	print


def sendCalculDelta():

	tmp_command = "echo \" \n ------ Start ------ \n \" " 
	subprocess.call(tmp_command, shell=True)

	for i in range(0,1):

		path = pathToFolder + 'Box_00' + str(i) + '/'

		for j in range(0,1):

			tmp_command = "echo " + str(i) + " " + str(j)
			subprocess.call(tmp_command, shell=True)

			path = pathToFolder + 'Box_00' + str(i) + '/Simu_00'+str(j) + '/'

			
			if (index_pass==0):
				### Get the data to 'good' files
				#####################
				#command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Get_delta/bin/main.exe " + str(i) + ' ' + str(j) + " 0 0 2 0\""
				command = "time /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Get_delta/bin/main.exe " + str(i) + ' ' + str(j) + " 0 0 2 0"
				print command
				subprocess.call(command, shell=True)
				myTools.isReadyForNewJobs(20, 430)
				time.sleep(10)
			

			if (index_pass==1):
				### Do the fits
				nbSpectra = 250000
				step  = 2000
				first = 0
				last  = step
				while (first <= nbSpectra):

					tmp_command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Get_delta/bin/main.exe " + str(i) + " " + str(j) + " " + str(first) + " " + str(last) +" 2 1 \""
					subprocess.call(tmp_command, shell=True)
	
					tmp_command = "echo " + tmp_command
					subprocess.call(tmp_command, shell=True)
					tmp_command = "echo " + time.ctime()
					subprocess.call(tmp_command, shell=True)
				
					first = first + step
					last  = last  + step

					myTools.isReadyForNewJobs(200, 430)
					time.sleep(0.1)
			

			if (index_pass==2):
				folder = path + '/Run/'
				scheme = 'alphaAndBeta_LYA_'
				lenScheme = len(scheme)
			
				tmp_command = "rm " + folder + scheme + 'all.txt'
				subprocess.call(tmp_command, shell=True)
			
				### Get the list of files
				all_t = os.listdir(folder)	
				tmp_all_t = []
				for el in all_t:
					if (el[:lenScheme]==scheme):
									tmp_all_t.append(el)
				all_t = tmp_all_t
			
				### Sort the list of files
				convert      = lambda text: int(text) if text.isdigit() else text
				alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
				all_t = sorted(all_t, key = alphanum_key)
			
				for el in all_t:
					print el
			
				data  = numpy.loadtxt(folder+all_t[0])
				idx   = data[:,0].astype(int)
				alpha = data[:,1]
				beta  = data[:,2]
				chi   = data[:,3]
				alphaErr = data[:,4]
				betaErr  = data[:,5]
				flag     = data[:,6].astype(int)
				nbPixel  = data[:,7].astype(int)
			
				for el in all_t[1:]:
					fo = open(folder+el, "r")
					print folder+el
					if os.fstat(fo.fileno()).st_size:
						data   = numpy.loadtxt(folder+el)
						if (data.size/4 != 1):
							idx    = numpy.append( idx,   data[:,0].astype(int) )
							alpha  = numpy.append( alpha, data[:,1]  )
							beta   = numpy.append( beta,  data[:,2]  )
							chi    = numpy.append( chi,   data[:,3]  )
							alphaErr  = numpy.append( alphaErr, data[:,4]  )
							betaErr   = numpy.append( betaErr,  data[:,5]  )
							flag      = numpy.append( flag,   data[:,6].astype(int)  )
							nbPixel   = numpy.append( nbPixel,   data[:,7].astype(int)  )
						else:
							idx    = numpy.append( idx,   data[0].astype(int) )
							alpha  = numpy.append( alpha, data[1]  )
							beta   = numpy.append( beta,  data[2]  )
							chi    = numpy.append( chi,   data[3]  )
							alphaErr  = numpy.append( alphaErr, data[4]  )
							betaErr   = numpy.append( betaErr,  data[5]  )
							flag      = numpy.append( flag,   data[6].astype(int)  )
							nbPixel   = numpy.append( nbPixel,   data[7].astype(int)  )
			
				### change 'nan' values into negative ones for the errors
				alphaErr[ numpy.isinf(alphaErr) ] = -1.
				betaErr[ numpy.isinf(betaErr) ] = -1.
				alphaErr[ numpy.isnan(alphaErr) ] = -1.
				betaErr[ numpy.isnan(betaErr) ] = -1.
			
			
				print
				print
				print '  nb spectra         : ', alpha.size
				print '  alpha              : ', alpha
				print '  beta               : ', beta
				print
				print
				##
				print '  nb pixel <=0              : ', chi[ (nbPixel<=0.) ].size
				print '  idx nbPixel <=0           : ', idx[ (nbPixel<=0.) ]
				print '  nbPixel of (nbPixel <=0)  : ', nbPixel[ (nbPixel<=0.) ]
				print
				##
				print '  nb of flag!=0      : ', chi[ (flag!=0) ].size
				print '  idx of flag!=0     : ', idx[ (flag!=0) ]
				##
				print '  nb of chi2=inf     : ', chi[ (numpy.isinf(chi)) ].size
				print '  idx of chi2=inf    : ', idx[ (numpy.isinf(chi)) ]
				##
				print '  nb of chi2==0      : ', chi[ (chi==0.) ].size
				print '  idx of chi2==0     : ', idx[ (chi==0.) ]
				##
				print '  nb of alpha==alphaStart__     : ', chi[ (alpha==alphaStart__) ].size
				print '  idx of alpha==alphaStart__    : ', idx[ (alpha==alphaStart__) ]
				##
				print '  nb of beta==0      : ', chi[ (beta==0.) ].size
				print '  idx of beta==0     : ', idx[ (beta==0.) ]
				##
				print '  nb of alphaErr<=0.      : ', chi[ (alphaErr<=0.) ].size
				print '  idx of alphaErr<=0.     : ', idx[ (alphaErr<=0.) ]
				##
				print '  nb of betaErr<=0.      : ', chi[ (betaErr<=0.) ].size
				print '  idx of betaErr<=0.     : ', idx[ (betaErr<=0.) ]
				##
				print '  nb of alphaErr>=alpha      : ', chi[ (alphaErr>=numpy.abs(alpha)) ].size
				print '  idx of alphaErr>=alpha     : ', idx[ (alphaErr>=numpy.abs(alpha)) ]
				##
				print '  nb of betaErr>=beta      : ', chi[ (betaErr>=numpy.abs(beta)) ].size
				print '  idx of betaErr>=beta     : ', idx[ (betaErr>=numpy.abs(beta)) ]
			
				saveAlpha = numpy.array(alpha)
			
				alpha[ (alpha==1.) ] = alphaStart__
				beta[ (alpha==1.) ]  = 0.
				alpha[ (nbPixel<=0.) ] = -600.
				beta[ (nbPixel<=0.) ]  = -600.
				alpha[ (flag!=0) ] = -300.
				beta[ (flag!=0) ]  = -300.
				alpha[ (alphaErr<=0.) ] = -400.
				beta[ (alphaErr<=0.) ]  = -400.
				alpha[ (betaErr<=0.) ] = -500.
				beta[ (betaErr<=0.) ]  = -500.
				alpha[ (numpy.isinf(chi)) ] = -100.
				beta[ (numpy.isinf(chi)) ]  = -100.
				alpha[ (numpy.isnan(chi)) ] = -100.
				beta[ (numpy.isnan(chi)) ]  = -100.
				alpha[ (chi==0.) ] = -200.
				beta[ (chi==0.) ]  = -200.
				
				path     = path + 'Data/delta.fits'
				file_cat = pyfits.open(path,mode='update')
				cat      = file_cat[1].data
				print
				print
				print path
				print '  nb spectra         : ', cat.size
				print '  alpha diff         : ', cat['ALPHA_2']-alpha
				print '  beta diff          : ', cat['BETA_2']-beta
				
				cat['ALPHA_2'] = alpha
				cat['BETA_2']  = beta
			
				### Flag for spectra with alpha_error > alpha
				cat['BETA_1'][ (alphaErr>=numpy.abs(saveAlpha)) ] = -600.
				### Keep the chi^{2} if not reObs
				cat['ALPHA_1'][ (nbPixel>0.) ]  = chi[ (nbPixel>0.) ]/nbPixel[ (nbPixel>0.) ]
				cat['ALPHA_1'][ (nbPixel<=0.) ] = 0.
			
				file_cat.close()
			

			if (index_pass==3):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"time  ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Correlation/xi_delta_QSO.py a a a a 0 " + str(i) + ' ' + str(j) + "\""
				print command
				time.sleep(1.)
				subprocess.call(command, shell=True)
			if (index_pass==4):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"time  ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Correlation/xi_Q_Q.py a a a a " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)

		time.sleep(10.)

def create_fits_qso(sizeMax=0):

	### Create a list of QSOs
	col_xx              = pyfits.Column(name='X',  format='D', array=numpy.zeros(sizeMax), unit='Mpc/h')
	col_yy              = pyfits.Column(name='Y',  format='D', array=numpy.zeros(sizeMax), unit='Mpc/h')
	col_zz              = pyfits.Column(name='Z',  format='D', array=numpy.zeros(sizeMax), unit='0 (redshift)')
	tbhduQSO = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz])

	return tbhduQSO

def create_fits_forest(sizeMax=0):

	### Create a list of forests
	tmp_nbBinForest     = str(nbPixel)+'D'

	plate                = pyfits.Column(name='PLATE',           format='J', array=numpy.zeros(sizeMax) )
	mjd                  = pyfits.Column(name='MJD',             format='J', array=numpy.zeros(sizeMax) )
	fiber                = pyfits.Column(name='FIBERID',         format='J', array=numpy.zeros(sizeMax) )
	
	ra                   = pyfits.Column(name='RA',              format='D', array=numpy.zeros(sizeMax))
	de                   = pyfits.Column(name='DEC',             format='D', array=numpy.zeros(sizeMax))
	zz                   = pyfits.Column(name='Z_VI',            format='D', array=numpy.zeros(sizeMax))
	nb                   = pyfits.Column(name='NB_PIXEL',        format='I', array=numpy.zeros(sizeMax))
	meanLambdaRF         = pyfits.Column(name='MEAN_FOREST_LAMBDA_RF', format='D', array=numpy.zeros(sizeMax), unit='angstrom')
	alpha1               = pyfits.Column(name='ALPHA_1',         format='D', array=numpy.ones(sizeMax) )
	beta1                = pyfits.Column(name='BETA_1',          format='D', array=numpy.zeros(sizeMax) )
	alpha2               = pyfits.Column(name='ALPHA_2',         format='D', array=alphaStart__*numpy.ones(sizeMax) )
	beta2                = pyfits.Column(name='BETA_2',          format='D', array=numpy.zeros(sizeMax) )

	tmp_nbBinForest     = str(nbPixel)+'D'
	lambdaForest        = pyfits.Column(name='LAMBDA_OBS',       format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)), unit='angstrom' )
	lambdaRFForest      = pyfits.Column(name='LAMBDA_RF',        format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)),  unit='angstrom' )
	normFluxForest      = pyfits.Column(name='NORM_FLUX',        format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	normFluxIvarForest  = pyfits.Column(name='NORM_FLUX_IVAR',   format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	fluxDLA             = pyfits.Column(name='FLUX_DLA',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)))
	deltaForest         = pyfits.Column(name='DELTA',            format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	deltaIvarForest     = pyfits.Column(name='DELTA_IVAR',       format=tmp_nbBinForest, array=numpy.zeros((sizeMax,nbPixel)))
	deltaWeight         = pyfits.Column(name='DELTA_WEIGHT',     format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)))
	template            = pyfits.Column(name='TEMPLATE',         format=tmp_nbBinForest, array=numpy.ones((sizeMax,nbPixel)))	
	
	tbhduForest = pyfits.BinTableHDU.from_columns([plate, mjd, fiber, ra, de, zz, nb, meanLambdaRF, alpha1, beta1, alpha2, beta2, lambdaForest, lambdaRFForest, normFluxForest, normFluxIvarForest, fluxDLA, deltaForest, deltaIvarForest, deltaWeight, template])

	return tbhduForest






## Execute the code
if (index_function==0):
	main()
elif (index_function==1):
	sendCalculDelta()






