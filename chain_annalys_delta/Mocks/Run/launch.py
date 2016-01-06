# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# created by Hélion du Mas des Bourboux
# <Helion.du-Mas-des-Bourboux@cea.fr>
# module add scipy/0.13.1_Python_2.7.6-gnu

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

commandProd     = "python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Mocks/test.py "
def main():

	for i in range(0,10):

		tmp_command = "echo \" \n ------ Start ------ \n \" " 
		subprocess.call(tmp_command, shell=True)
	
		for j in range(0,10):
			
			nameMock = str(i) + "_" + str(j)
			
			tmp_command = "echo \" \n ------ Mock: " + str(i) + "  " + str(j) + " ------\n \" " 
			subprocess.call(tmp_command, shell=True)	

			tmp_command = "clubatch \"echo ; hostname ; " + commandProd + str(i) + " " + str(j) + "\""
			subprocess.call(tmp_command, shell=True)
	
			tmp_command = "echo " + tmp_command
			subprocess.call(tmp_command, shell=True)
			tmp_command = "echo " + time.ctime()
			subprocess.call(tmp_command, shell=True)
					
			myTools.isReadyForNewJobs(150, 430)
			time.sleep(5)

main()
