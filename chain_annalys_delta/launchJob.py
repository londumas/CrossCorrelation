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
#commandProd     = "python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/creatFitsFile.py"
#commandProd     = "python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/creat_Fits_File.py"
#commandProd     = "python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/All_Spectrum/get_all_spectrum.py"
#commandProd     = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/chain_annalys_delta/Get_delta/bin/main.exe"
#commandProd     = "python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/All_Spectrum/Iterative_spectrum/iterative_spectrum.py"
commandProd     = "python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/creat_Fits_File.py "

commandQstatMe  = "qstat | grep echo | wc -l"
commandQstatAll = "qstat -u '*' | wc -l"

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
nbSpectra = 600000

##
#nbSpectra = 200000

def main():

	print
	print "------ Start ------"
	print

	## Launch the production chain
	#####################
	last  = step-1
	first = 0
	while (first <= nbSpectra):
	
		#myTools.isReadyForNewJobs(150, 430)
		time.sleep(0.5)
		
		tmp_command = "clubatch \"echo ; hostname ; " + commandProd + ' ' + str(first) + " " + str(last) + " 0 0 \""
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
main()








