# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import sys
import numpy
import matplotlib.pyplot as plt
import myTools
from const_delta import *

forest__ = sys.argv[1]
#path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_000/Simu_000/Results/'
path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_test_PDFMocksJMC/'


#lineSeen_CIV = [1.00166645,1.00522684,1.01407602,1.01576593,1.02024425,1.03718845,1.03891687,1.05354069,1.07544891,1.07738686,1.07918228]
if (forest__=='LYA'):
	lines = LYA_lines
	names = LYA_lines_names
if (forest__=='CIV'):
	lines = CIV_lines
	names = CIV_lines_names
if (forest__=='MGII'):
	lines = MGII_lines
	names = MGII_lines_names
if (forest__=='LYB'):
	lines = LYB_lines
	names = LYB_lines_names
if (forest__=='SIIV'):
	lines = SIIV_lines
	names = SIIV_lines_names

def plot():

	path = path__ +'xi_1DlRFDevide_delta_delta_'+forest__+'.txt'
	print path
	data = numpy.loadtxt(path)

	### remove empty pixels
	xxx = data[:,2][ (data[:,5]!=0.) ]/data[:,4][ (data[:,5]!=0.) ]
	yyy = data[:,0][ (data[:,5]!=0.) ]/data[:,4][ (data[:,5]!=0.) ]
	yer = numpy.sqrt( (data[:,1][ (data[:,5]!=0.) ]/data[:,4][ (data[:,5]!=0.) ] -yyy**2.)/data[:,5][ (data[:,5]!=0.) ]  )	
	
	plt.errorbar(xxx, yyy, yerr=yer, marker='o')

	### Show lines in the correlation
	xMax    = numpy.amax(xxx)
	yMin    = numpy.amin(yyy)
	yMax    = numpy.amax(yyy)
	nbLines = lines.size
	for i in range(0,nbLines):
		for j in range(0,i):
			#if (names[i][:3]!=forest__ and names[j][:3]!=forest__): continue
			line = lines[j]/lines[i]

			'''
			isPresent = False
			for k in lineSeen_CIV:
				if (abs(line-k)<=0.00005): isPresent = True
			if (not isPresent): continue
			'''

			if (line>xMax): line = 1./line
			#print ' ||  ', names[i] ,' - ', names[j], ' || ', line, ' || ', lines[i], ' || ', lines[j], ' || '

			xLi = [line,line]
			yLi = [yMin,yMax]
			name = names[i]+' - '+names[j]
			plt.plot(xLi,yLi,color='green')
			plt.text(line, yMax, name, rotation='vertical', fontsize=20)


	plt.title(r'$1D: \, \delta_{'+forest__+'} \, - \, \delta_{'+forest__+'} $', fontsize=40)
	plt.xlabel(r'$\lambda_{1}/\lambda_{2}$', fontsize=40)
	plt.ylabel(r'$\xi(\lambda_{1}/\lambda_{2})$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.xlim([ 0.99*numpy.min(xxx), 1.01*numpy.max(xxx) ])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
	
	
	return xxx, yyy, yer
def loadData(path1D):


	print path1D
	data = numpy.loadtxt(path1D)

	cut = (data[:,5]!=0.)
	size = cut[cut].size
	xi1D = numpy.zeros(shape=(size,3))

	xi1D[:,0] = data[:,2][cut]/data[:,4][cut]
	xi1D[:,1] = data[:,0][cut]/data[:,4][cut]
	xi1D[:,2] = numpy.sqrt( (data[:,1][cut]/data[:,4][cut] -xi1D[:,1]**2.)/data[:,5][cut]  )

	return xi1D
def saveListRealMocks(ni,nj):

	pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_1DlRF_delta_delta_LYA_result'

	listSimu = numpy.zeros( shape=(161,ni*nj) )

	for i in numpy.arange(ni):
		for j in numpy.arange(nj):
			path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_1DlRF_delta_delta_LYA.txt'
			xi1D = loadData(path)
			listSimu[:,i*10+j] = xi1D[:,1]

	numpy.save(pathToSave,listSimu)
	numpy.save(pathToSave+'_cov',numpy.cov(listSimu))

	return


path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests6/'
xxx, yyy, yer = plot()
path__ = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
forest__ = 'CIV'
xxx3, yyy3, yer3 = plot()

'''
xxx = xxx[50:]
yyy = yyy[50:]
yer = yer[50:]
xxx = xxx[:200]
yyy = yyy[:200]
yer = yer[:200]
xxx2 = xxx2[50:]
yyy2 = yyy2[50:]
yer2 = yer2[50:]
xxx3 = xxx3[50:]
yyy3 = yyy3[50:]
yer3 = yer3[50:]
'''
plt.errorbar(xxx, yyy, yerr=yer, marker='o',label='eBOSS')
#plt.errorbar(xxx2,yyy2, yerr=yer2, marker='o',label='SIIV-SIIV')
plt.errorbar(xxx3,yyy3, yerr=yer3, marker='o',label='DR12',color='red')

'''
### Show lines in the correlation
xMax    = numpy.amax(xxx)
yMin    = numpy.amin(yyy)
yMax    = 0.0015 #numpy.amax(yyy)
nbLines = lines.size
for i in range(0,nbLines):
	for j in range(0,i):
		if (names[i]=='NV_a' or names[j]=='NV_a'): continue
		if (names[i]=='NV_b' or names[j]=='NV_b'): continue
		line = lines[j]/lines[i]
		if (line>xMax): line = 1./line
		print ' ||  ', names[i] ,' - ', names[j], ' || ', line, ' || ', lines[i], ' || ', lines[j], ' || '

		xLi = [line,line]
		yLi = [yMin,yMax]
		name = names[i]+' - '+names[j]
		plt.plot(xLi,yLi,color='green')
		plt.text(line, yMax, name, rotation='vertical', fontsize=20)
'''

plt.xlabel(r'$\lambda_{1}/\lambda_{2}$', fontsize=40)
plt.ylabel(r'$\xi(\lambda_{1}/\lambda_{2})$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ 0.99*numpy.min(xxx), 1.01*numpy.max(xxx) ])
plt.show()









