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




pathToFolder = '/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/'

def main():


	tmp_command = "echo \" \n ------ Start ------ \n \" " 
	subprocess.call(tmp_command, shell=True)

	for i in range(0,10):		

		path = pathToFolder + 'M3_0_' + str(i) + '/'

		## Create the folders
		#####################
		#subprocess.call('mkdir ' +path, shell=True)

		for j in range(0,10):

			print i, j

			path = pathToFolder+ 'M3_0_' + str(i) + '/00'+str(j) + '/'

			### Create the folders and FITS file
			#####################
			#subprocess.call('mkdir ' + path, shell=True)
			#subprocess.call('mkdir ' + path + 'Data', shell=True)
			#subprocess.call('mkdir ' + path + 'Results', shell=True)
			


def sendCalculDelta():

	stepIdx = 3


	tmp_command = "echo \" \n ------ Start ------ \n \" " 
	subprocess.call(tmp_command, shell=True)

	for i in range(0,10):		

		path = pathToFolder + 'M3_0_' + str(i) + '/'

		for j in range(0,10):

			tmp_command = "echo " + str(i) + " " + str(j)
			subprocess.call(tmp_command, shell=True)

			path = pathToFolder+ 'M3_0_' + str(i) + '/00'+str(j) + '/'

			if (stepIdx==0):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"echo ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Get_delta/bin/main.exe " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)
				myTools.isReadyForNewJobs(10, 430)
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

					myTools.isReadyForNewJobs(150, 430)
					time.sleep(1.)
			

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
		
				print cat['ALPHA_2']
				
				cat['ALPHA_2'] = alpha
				cat['BETA_2']  = beta
		
				print cat['ALPHA_2']
		
				file_cat.close()
			
			if (stepIdx==3):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"time  ; hostname ; /home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Correlation/bin/main.exe 5 0 0 0 " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)
				time.sleep(10)
			if (stepIdx==4):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"time  ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Correlation/xi_delta_QSO.py a a a a 0 " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)
			if (stepIdx==5):
				### Get the data to 'good' files
				#####################
				command = "clubatch \"time  ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Correlation/xi_Q_Q.py a a a a " + str(i) + ' ' + str(j) + "\""
				print command
				subprocess.call(command, shell=True)

		time.sleep(120)


## Execute the code
#main()
sendCalculDelta()







