# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# created by Hélion du Mas des Bourboux
# <Helion.du-Mas-des-Bourboux@cea.fr>
# module add scipy/0.13.1_Python_2.7.6-gnu

import subprocess
import time
import myTools

commandProd     = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Correlation/bin/main.exe"

def main():

	
	### For mocks

	### mock catalogue
	for i in range(0,1):
		### mock catalogue
		for j in range(0,1):

			### for region
			for k in range(0,80):

			### for metals 
			#for k in [21,24,25,26,27]:
				tmp_command = "echo \" \n ------ Start ------ \n \" " 
				subprocess.call(tmp_command, shell=True)

				tmp_command = "clubatch \"time ; hostname ; " + commandProd + ' 21 0 ' + str(k) + ' 0 ' + str(i) + " " + str(j) + "\""
				subprocess.call(tmp_command, shell=True)

				tmp_command = "echo " + tmp_command
				subprocess.call(tmp_command, shell=True)
				tmp_command = "echo " + time.ctime()
				subprocess.call(tmp_command, shell=True)

				#time.sleep(120)
                	        #myTools.isReadyForNewJobs(60, 1000,'time')
		#time.sleep(120)

	
	'''
	for i in range(0,29):
		for j in range(0,29):
			tmp_command = "echo \" \n ------ Start ------ \n \" " 
			subprocess.call(tmp_command, shell=True)

			tmp_command = "clubatch \"time ; hostname ; " + commandProd + " 20 " + str(i) + ' ' + str(j) + " 0 0 0 \""
			subprocess.call(tmp_command, shell=True)
	
			tmp_command = "echo " + tmp_command
			subprocess.call(tmp_command, shell=True)
			tmp_command = "echo " + time.ctime()
			subprocess.call(tmp_command, shell=True)

			time.sleep(1)
			myTools.isReadyForNewJobs(200, 430,'time')
	'''
	
	'''
	### for metals templates
	for k in range(0,29):
		tmp_command = "echo \" \n ------ Start ------ \n \" " 
		subprocess.call(tmp_command, shell=True)

		tmp_command = "clubatch \"time ; hostname ; " + commandProd + ' 15 0 ' + str(k) + " 0 0 0 \""
		subprocess.call(tmp_command, shell=True)

		tmp_command = "echo " + tmp_command
		subprocess.call(tmp_command, shell=True)
		tmp_command = "echo " + time.ctime()
		subprocess.call(tmp_command, shell=True)

		time.sleep(10)
		myTools.isReadyForNewJobs(100, 430,'time')
	'''

	'''	
	for i in range(0,10):
		for j in range(0,10):
			### Get the data to 'good' files
			#####################
			command = "clubatch \" time ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/main_Q.py " + str(i) + ' ' + str(j) + "\""
			print command
			subprocess.call(command, shell=True)
			time.sleep(0.1)
		time.sleep(10)
	'''

	'''
	for i in range(0,10000):
		### Get the data to 'good' files
		#####################
		command = "clubatch \"time  ; hostname ; python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/Correlation/xi_delta_QSO.py a a a a " + str(i) + "\""
		print command
		subprocess.call(command, shell=True)
		isReadyForNewJobs(100, 430,'time')
		time.sleep(10)
	'''

def isReadyForNewJobs(max_nbJobsMe, max_nbJobsAll,word='echo'):
	'''
		This function ends when there is room in the cluster for a job.
		

	'''

	commandQstatMe  = "qstat | grep "+ word +" | wc -l"
	commandQstatAll = "qstat -u '*' | wc -l"

	## Duration to wait before checking again
	waitTime = 1
	## Number of lines that aren't jobs
	noJobsLine = 2
	
	notDoProd = True
	while (notDoProd):
		
		## Counts the jobs from me
		nbJobsMe = int(subprocess.check_output(commandQstatMe, shell=True))
		
		## Counts jobs from everybody
		nbJobsAll = int(subprocess.check_output(commandQstatAll, shell=True))
		if (nbJobsAll>0):
			nbJobsAll = nbJobsAll-noJobsLine

		## Look if there are too many jobs
		if (nbJobsMe < max_nbJobsMe and nbJobsAll < max_nbJobsAll):
			notDoProd = False
		else:
			time.sleep(waitTime)

	## Return if there is room for a new job
	return




	print
	print " ------ End ------"
	print


## Execute the code
main()
