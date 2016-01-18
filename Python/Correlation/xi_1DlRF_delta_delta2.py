# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import sys
import numpy
import matplotlib.pyplot as plt
import myTools
from const_delta import *

forest1 = sys.argv[1]
forest2 = sys.argv[2]
qso1    = sys.argv[3]
qso2    = sys.argv[4]

path1 = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'

def plot():

	if (forest1=='LYA'):
		lines1 = LYA_lines
		names1 = LYA_lines_names
	if (forest1=='CIV'):
		lines1 = CIV_lines
		names1 = CIV_lines_names
	if (forest1=='MGII'):
		lines1 = MGII_lines
		names1 = MGII_lines_names
	if (forest1=='LYB'):
		lines1 = LYB_lines
		names1 = LYB_lines_names
	if (forest1=='SIIV'):
		lines1 = SIIV_lines
		names1 = SIIV_lines_names

	if (forest2=='LYA'):
		lines2 = LYA_lines
		names2 = LYA_lines_names
	if (forest2=='CIV'):
		lines2 = CIV_lines
		names2 = CIV_lines_names
	if (forest2=='MGII'):
		lines2 = MGII_lines
		names2 = MGII_lines_names
	if (forest2=='LYB'):
		lines2 = LYB_lines
		names2 = LYB_lines_names
	if (forest2=='SIIV'):
		lines2 = SIIV_lines
		names2 = SIIV_lines_names

	path = path1 +'xi_1DlRF_delta_delta2_'+forest1+'_'+forest2+'.txt'
	print path
	data = numpy.loadtxt(path)
	xxx  = data[:,0]
	yyy  = data[:,1]
	yer  = data[:,2]
	
	plt.errorbar(xxx, yyy, yerr=yer, marker='o')
	
	### Show lines in the correlation
	xMin    = numpy.amin(xxx)
	xMax    = numpy.amax(xxx)
	yMin    = numpy.amin(yyy)
	yMax    = numpy.amax(yyy)
	for i in range(0,len(lines1)):
		for j in range(0,len(lines2)):
			#if (names1[i][:4]!='SiIV' and names1[j][:4]!='SiIV'): continue
			line = lines1[i]/lines2[j]
			if (line<xMin or line>xMax): continue
			xLi = [line,line]
			yLi = [yMin,yMax]
			name = names1[i]+' - '+names2[j]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, name, rotation='vertical', fontsize=20)
	

	plt.title(r'$1D: \, \delta_{'+forest1+'} \, - \, \delta_{'+forest2+'} $', fontsize=40)
	plt.xlabel(r'$\lambda_{1}/\lambda_{2}$', fontsize=40)
	plt.ylabel(r'$\xi(\lambda_{1}/\lambda_{2})$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	#plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()
	
	return

plot()
