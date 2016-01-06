# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# created by Hélion du Mas des Bourboux
# <Helion.du-Mas-des-Bourboux@cea.fr>
#
# Allows to send many jobs
#
# module add scipy/0.13.1_Python_2.7.6-gnu

import subprocess
import time
import myTools

#commandProd     = "/home/gpfs/manip/mnt/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/bin/ReadSpectra"
commandProd     = "python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/creatFitsFile.py"

## Parameters
#step = 7000
step = 5000

## nbSpectra
#############

#### DR12
#################
#nbSpectra = 298000

## Mock DR11
#nbSpectra = 245000

#### eBOSS
#################
#nbSpectra = 58000
#nbSpectra = 41150
#nbSpectra = 10000

#### DR12 + eBOSS
#################
nbSpectra = 300000

chunckNb = 0
simulNb  = 9

##
#nbSpectra = 200000

def main():

	print
	print "------ Start ------"
	print

	for i in range(0,2):
		for j in range(0,10):

			'''
			tmp_command = "clubatch \"echo ; hostname ; " + commandProd + ' ' + str(0) + ' ' + str(0) + ' ' + str(i) + ' ' + str(j) + "\""
                        subprocess.call(tmp_command, shell=True)

			tmp_command = "echo " + tmp_command
                        subprocess.call(tmp_command, shell=True)
                        tmp_command = "echo " + time.ctime()
                        subprocess.call(tmp_command, shell=True)
			time.sleep(60)
			'''
			
			## Launch the production chain
			#####################
			last  = step-1
			first = 0
			while (first <= nbSpectra):
		
				myTools.isReadyForNewJobs(150, 430)
				time.sleep(0.5)
		
				tmp_command = "clubatch \"echo ; hostname ; " + commandProd + ' ' + str(first) + ' ' + str(last) + ' ' + str(i) + ' ' + str(j) + "\""
				subprocess.call(tmp_command, shell=True)
		
				tmp_command = "echo " + tmp_command
				subprocess.call(tmp_command, shell=True)
				tmp_command = "echo " + time.ctime()
				subprocess.call(tmp_command, shell=True)
		
				first = first + step
				last  = last  + step
			
	print
	print " ------ End ------"
	print

def sendDR12():

	## Launch the production chain
	#####################
	last  = step-1
	first = 0
	while (first <= nbSpectra):
		
		myTools.isReadyForNewJobs(150, 430)
		time.sleep(0.5)
		
		tmp_command = "clubatch \"echo ; hostname ; " + commandProd + ' ' + str(first) + ' ' + str(last) + ' ' + str(0) + ' ' + str(0) + "\""
		subprocess.call(tmp_command, shell=True)
	
		tmp_command = "echo " + tmp_command
		subprocess.call(tmp_command, shell=True)
		tmp_command = "echo " + time.ctime()
		subprocess.call(tmp_command, shell=True)
		
		first = first + step
		last  = last  + step

	print
	print " ------ End ------"
	print




## Execute the code
#main()
sendDR12()







