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




pathToFolder = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/'
#pathToFolder = '/home/gpfs/manip/mnt0607/bao/Mock_JMLG/v1528/'

def main():


	'''
	### Constants
	sizeMax = 300000
	nbQSO__ = 238929
	nbFor__ = 170898
	ratioForestToQSO__    = 1.*nbFor__/nbQSO__;
	
	
	### Create a list of QSOs
	col_xx              = pyfits.Column(name='X',  format='D', array=numpy.zeros(sizeMax), unit='Mpc/h')
	col_yy              = pyfits.Column(name='Y',  format='D', array=numpy.zeros(sizeMax), unit='Mpc/h')
	col_zz              = pyfits.Column(name='Z',  format='D', array=numpy.zeros(sizeMax), unit='0 (redshift)')
	tbhduQSO = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz])

	
	### Create a list of forests
	sizeMax = 250000
	nbPixel = 645
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
	alpha2               = pyfits.Column(name='ALPHA_2',         format='D', array=numpy.ones(sizeMax) )
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
	'''


	tmp_command = "echo \" \n ------ Start ------ \n \" " 
	subprocess.call(tmp_command, shell=True)

	## Create the folders
	#####################
	#subprocess.call('mkdir ' +pathToFolder+'Results/', shell=True)

	for i in range(0,10):		

		path = pathToFolder + 'Box_00' + str(i) + '/'

		## Create the folders
		#####################
		#subprocess.call('mkdir ' +path, shell=True)

		for j in range(0,10):

			print i, j

			path = pathToFolder+ 'Box_00' + str(i) + '/Simu_00'+str(j) + '/'

			
			### Create the folders and FITS file
			#####################
			#subprocess.call('mkdir ' + path, shell=True)
			#subprocess.call('mkdir ' + path + 'Data', shell=True)
			#subprocess.call('mkdir ' + path + 'Results', shell=True)
			#tbhduQSO.writeto(path + 'Data/QSO_withRSD.fits', clobber=True)
			#tbhduQSO.writeto(path + 'Data/QSO_noRSD.fits', clobber=True)
			#tbhduForest.writeto(path + 'Data/delta.fits', clobber=True)
			#tbhduForest.writeto(path + 'Data/delta_noNoise_noCont.fits', clobber=True)
			
			subprocess.call('mkdir ' + path + '/Results/BaoFit_q_f', shell=True)
			#subprocess.call('mkdir ' + path + 'BaoFit_q_f_covFromFit_fixedBAO_FixAlphaParal', shell=True)
			#subprocess.call('mkdir ' + path + 'Results_RandomPosInCell', shell=True)

			'''
			command = 'clubatch cp ' + pathToFolder + 'Box_000/Simu_000/Data/delta_noNoise_noCont.fits ' + path + 'Data/delta_noNoise_noCont.fits'
			subprocess.call(command, shell=True)
			command = 'clubatch cp ' + pathToFolder + 'Box_000/Simu_000/Data/QSO_withRSD.fits ' + path + 'Data/QSO_withRSD.fits'
			subprocess.call(command, shell=True)
			myTools.isReadyForNewJobs(10, 430,'cp')
			time.sleep(20)
			'''

			'''
			command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Mock_JMLG/ReadFits/bin/main.exe " + str(i) + ' ' + str(j) + "\""
			print command
			subprocess.call(command, shell=True)
			myTools.isReadyForNewJobs(10, 430,'echo')
			time.sleep(20)
			'''
			

			'''
			nbQSO = 0
			### Remove useless lines in QSO
			#####################
			for el in ['QSO_withRSD.fits']: # ,'QSO_noRSD.fits']:
				cat = pyfits.open(path + 'Data/'+el, memmap=True)[1].data
				print cat.size
				cat = cat[ (cat['Z'] != 0.) ]
				print cat.size
				nbQSO = cat.size
				pyfits.writeto(path + 'Data/'+el, cat, clobber=True)
			
			print nbQSO
			'''
			'''
			### Remove useless lines in Forest
			#####################
			for el in ['delta_noNoise_noCont.fits']: #,'delta.fits','delta_noNoise_noCont.fits']:
				cat = pyfits.open(path + 'Data/'+el, memmap=True)[1].data[:nbQSO]
				print cat.size
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
				cat = cat[ (cat['Z_VI'] != 0.) ]
				print cat.size
				pyfits.writeto(path + 'Data/'+el, cat, clobber=True)
			'''

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

	stepIdx = 3


	tmp_command = "echo \" \n ------ Start ------ \n \" " 
	subprocess.call(tmp_command, shell=True)

	for i in range(0,10):

		path = pathToFolder + 'Box_00' + str(i) + '/'

		for j in range(0,10):

			tmp_command = "echo " + str(i) + " " + str(j)
			subprocess.call(tmp_command, shell=True)

			path = pathToFolder + 'Box_00' + str(i) + '/Simu_00'+str(j) + '/'

			if (stepIdx==0):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Get_delta/bin/main.exe " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)
				myTools.isReadyForNewJobs(20, 430)
				time.sleep(10)
			

			if (stepIdx==1):
				### Do the fits
				nbSpectra = 200000
				step  = 2000
				first = 0
				last  = step
				while (first <= nbSpectra):

					tmp_command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Get_delta/bin/main.exe " + str(i) + " " + str(j) + " " + str(first) + " " + str(last) +"\""
					subprocess.call(tmp_command, shell=True)
	
					tmp_command = "echo " + tmp_command
					subprocess.call(tmp_command, shell=True)
					tmp_command = "echo " + time.ctime()
					subprocess.call(tmp_command, shell=True)
				
					first = first + step
					last  = last  + step

					myTools.isReadyForNewJobs(200, 430)
					time.sleep(0.1)
			

			if (stepIdx==2):
				### Put files of alpha and beta together
		
				folder = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/'
				scheme = 'alphaAndBeta_'
				lenScheme = len(scheme)
				endScheme = '_'+str(i)+'_'+str(j)+'.txt'
				lenEndScheme = -len(endScheme)
		
				tmp_command = "rm " + folder + scheme + 'all.txt'
				subprocess.call(tmp_command, shell=True)
		
				### Get the list of files
				all_t = os.listdir(folder)	
				tmp_all_t = []
				for el in all_t:
					if (el[:lenScheme]==scheme and el[lenEndScheme:]==endScheme):
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
						else:
							idx    = numpy.append( idx,   data[0].astype(int) )
							alpha  = numpy.append( alpha, data[1]  )
							beta   = numpy.append( beta,  data[2]  )
							chi    = numpy.append( chi,   data[3]  )
		
				print idx
				print alpha
				print beta
				print idx.size
		
				numpy.savetxt(folder+scheme+'all.txt', zip(idx,alpha,beta,chi))
				'''
				plt.hist(alpha,bins=100,label='alpha')
				plt.show()
				plt.hist(beta,bins=100,label='beta')
				plt.show()
				plt.hist(chi[ (numpy.isfinite(chi)) ],bins=100,label='beta')
				plt.show()
				'''
				
				### Put the new alpha and beta in fits file
				path     = path + 'Data/delta.fits'
				print path
				file_cat = pyfits.open(path,mode='update')
				cat      = file_cat[1].data
		
				print cat['ALPHA_2']-alpha
				
				cat['ALPHA_2'] = alpha
				cat['BETA_2']  = beta
		
				print cat['ALPHA_2']
		
				file_cat.close()
			

			if (stepIdx==3):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"time  ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Correlation/xi_delta_QSO.py a a a a 0 " + str(i) + ' ' + str(j) + "\""
				print command
				time.sleep(1.)
				subprocess.call(command, shell=True)
			if (stepIdx==4):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"time  ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Correlation/xi_Q_Q.py a a a a " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)

		time.sleep(10.)



## Execute the code
main()
#sendCalculDelta()







