# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

### My tools
import myTools


import astropy.io.fits as pyfits
import math
import numpy
import os
import decimal ## To set the precision of the double values
import cosmolopy.distance as cosmology
import matplotlib.pyplot as plt
from scipy import interpolate


import myTools
from myTools import Get_TProfile
from const_delta import *



pathQSO   = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_ALL_TESTS.fits'
path = '/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/histos/'
pathData  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits'
area_data = 9078.1
nbQSO = 232399
nbFor = 166968

pathMocks = '/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/mock.fits'


pathSimu    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_000/Simu_000/Data/'
rawPathSimu = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/'
zKey = 'Z'
alphaKey = 'ALPHA'
area_simu = 9078.1
chunckNb = 10
simulNb  = 10
show_mean = False

num_plots = 10
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in numpy.linspace(0, 0.9, num_plots)])

def comparePlot():

	'''
	## Map data
	catSimu = pyfits.open(pathQSO,memmap=True)[1].data
	#catMock = pyfits.open(pathMocks,memmap=True)[1].data
	plt.plot(catSimu['RA'],  catSimu['DEC'], linestyle="", marker="o",label=r'$Data \, QSO$')
	catSimu = pyfits.open(pathData,memmap=True)[1].data
	plt.plot(catSimu['RA'],  catSimu['DEC'], linestyle="", marker="o",label=r'$Data \, Forest$')
	plt.plot(catMock['RA'],  catMock['DEC'], linestyle="", marker="o",label=r'$Mock \, Forest$')
	plt.xlabel(r'$R.A.$')
	plt.ylabel(r'$Dec.$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	del catSimu
	del catMock
	
	## Map simu
	catSimu = pyfits.open(pathSimu+'QSO_withRSD.fits',memmap=True)[1].data
	print int((numpy.amax(catSimu['X'])-numpy.amin(catSimu['X']))/(4.5*0.7)) + 1
	print int((numpy.amax(catSimu['Y'])-numpy.amin(catSimu['Y']))/(4.5*0.7)) + 1
	plt.plot(catSimu['X'],  catSimu['Y'], linestyle="", marker="o",label=r'$Simu \, QSO$')
	catSimu = pyfits.open('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Data/delta.fits',memmap=True)[1].data
	plt.plot(catSimu['RA'],  catSimu['DEC'], linestyle="", marker="o",label=r'$Simu \, Forest$')
	plt.xlabel(r'$X$')
	plt.ylabel(r'$Y$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	del catSimu
	'''
	
	'''
	## Distribution redshift QSO
	catData = pyfits.open(pathQSO,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'QSO_withRSD.fits',memmap=True)[1].data
	print '  catData = ', catData.size
	print '  catSimu = ', catSimu.size
	print numpy.min(catSimu['Z']), numpy.max(catSimu['Z'])
	print ' QSO between 1.7 and 3.7  = ', catData['Z'][ numpy.logical_and(catData['Z']>1.7,catData['Z']<3.7) ].size
	hist = myTools.GetHisto(catData['Z'],numpy.arange(1.6,5.8,0.1))
	hist[:,1] /= area_data
	plt.plot(hist[:,0],hist[:,1],drawstyle='steps',label=r'$Data$',color='blue',linewidth=2)
	hist = myTools.GetHisto(catSimu['Z'],numpy.arange(1.6,5.8,0.1))
	hist[:,1] /= area_simu
	plt.plot(hist[:,0],hist[:,1],drawstyle='steps',label=r'$Simulation$',color='red',linewidth=2)
	plt.xlabel(r'$z_{QSO}$')
	plt.ylabel(r'$\# \cdot deg^{-2}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()

        print numpy.amin(catSimu['X']), numpy.amax(catSimu['X'])
        print numpy.amin(catSimu['Y']), numpy.amax(catSimu['Y'])
        print
        print (numpy.amax(catSimu['X'])+numpy.amin(catSimu['X']))/(4.5*0.7)
        print (numpy.amax(catSimu['Y'])+numpy.amin(catSimu['Y']))/(4.5*0.7)
        print
        print numpy.amin(catSimu['Z']), numpy.amax(catSimu['Z'])

	del catData
	del catSimu
	
	## Get number of QSO
	nb = numpy.array([])
	vel = numpy.array([])
	vel_all = numpy.array([])
	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_withRSD.fits'
			catSimu = pyfits.open(tmp_path,memmap=True)[1].data

			nb = numpy.append( nb,[catSimu.size] )
			z_withRSD = catSimu['Z']

			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_noRSD.fits'
                        catSimu = pyfits.open(tmp_path,memmap=True)[1].data
			z_noRSD = catSimu['Z']
			vel = numpy.append( vel,[ numpy.mean( (z_withRSD-z_noRSD)/(1.+z_noRSD) )  ] )
			vel_all = numpy.append( vel_all, (z_withRSD-z_noRSD)/(1.+z_noRSD) )

			print i, j, catSimu.size
	
			del catSimu
	
	print '  mean number of QSO = ', numpy.mean(nb)
	print '  compare to data it is = ', (numpy.mean(nb)-nbQSO)/nbQSO*100.
	plt.plot(numpy.arange(chunckNb*simulNb), nb,linewidth=2, label='Simulation',color='blue',marker='o')
	plt.plot( [0.,chunckNb*simulNb], [nbQSO,nbQSO],color='red' ,linewidth=2, label='Data')
	plt.xlabel(r'$Mock \, index$')
	plt.ylabel(r'$\# \, QSO$')
	myTools.deal_with_plot(False,False,True)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()

	plt.hist(nb,bins=10,linewidth=2, label='Simulation',color='blue')
	plt.plot( [nbQSO,nbQSO], [0.,30.],color='red' ,linewidth=2, label='Data')
	plt.xlabel(r'$\# \, QSO$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	plt.show()
	

	print vel
	print '  mean velocity = ', numpy.mean(vel)
        plt.plot(numpy.arange(chunckNb*simulNb), vel,linewidth=2, label='Simulation',color='blue',marker='o')
        plt.xlabel(r'$Mock \, index$')
        plt.ylabel(r'$v_{QSO}/c$')
        myTools.deal_with_plot(False,False,True)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.show()

        plt.hist(vel,bins=10,linewidth=2, label='Simulation',color='blue')
        plt.xlabel(r'$v_{QSO}/c$')
        plt.ylabel(r'$\#$')
        myTools.deal_with_plot(False,False,True)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.show()

        plt.hist(vel_all,bins=10000,linewidth=2, label='Simulation',color='blue')
        plt.xlabel(r'$v_{QSO}/c$')
        plt.ylabel(r'$\#$')
        myTools.deal_with_plot(False,False,True)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.show()

	
	## Get number of Forest
	nb = numpy.array([])
	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/delta.fits'
			catSimu = pyfits.open(tmp_path,memmap=True)[1].data

			nb = numpy.append( nb,[catSimu.size] )
			print i, j, catSimu.size

			del catSimu

	print '  mean number of forest = ', numpy.mean(nb)
	plt.plot(numpy.arange(chunckNb*simulNb), nb,linewidth=2, label='Simulation',color='blue',marker='o')
	plt.plot( [0.,chunckNb*simulNb], [nbFor,nbFor],color='red' ,linewidth=2, label='Data')
	plt.xlabel(r'$Mock \, index$')
	plt.ylabel(r'$\# \, Forest$')
	myTools.deal_with_plot(False,False,True)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()

	plt.hist(nb,bins=10,linewidth=2, label='Simulation',color='blue')
	plt.plot( [nbFor,nbFor], [0.,30.],color='red' ,linewidth=2, label='Data')
	plt.xlabel(r'$\# \, Forest$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	plt.show()

	
	## Distribution redshift Forest
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	#catMock = pyfits.open(pathMocks,memmap=True)[1].data
	print '  catData = ', catData.size
	print '  catSimu = ', catSimu.size
	#print '  catMock = ', catMock.size
	print ' Forest between 1.7 and 3.7  = ', catData['Z'][ numpy.logical_and(catData['Z']>1.7,catData['Z']<3.7) ].size
	hist = myTools.GetHisto(catData['Z'],numpy.arange(1.9,5.8,0.1))
	hist[:,1] /= area_data
	plt.plot(hist[:,0],hist[:,1],drawstyle='steps',label=r'$Data$',color='blue',linewidth=2)
	hist = myTools.GetHisto(catSimu[zKey],numpy.arange(1.9,5.8,0.1))
	hist[:,1] /= area_simu
	plt.plot(hist[:,0],hist[:,1],drawstyle='steps',label=r'$Simulation$',color='red',linewidth=2)
	#plt.hist(catMock['Z_VI'], bins=numpy.arange(1.7,4.,0.1),histtype='step',label=r'$Mock$',color='green',linewidth=2)
	plt.xlabel(r'$z_{Forest}$')
	plt.ylabel(r'$\# \cdot deg^{-2}$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	del catData
	del catSimu
	#del catMock
	
	
	### Distribution ALPHA
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	print ' data <alpha> = ', numpy.mean(catData['ALPHA'])
	print ' simu <alpha> = ', numpy.mean(catSimu['ALPHA'])
	data = myTools.GetHisto(catData['ALPHA'], numpy.arange(-50.,50.,0.01))
	plt.plot(data[:,0],data[:,1],label=r'$Data$',color='blue', markersize=8,linewidth=2)
	data = myTools.GetHisto(catSimu[alphaKey], numpy.arange(-50.,50.,0.01))
	plt.plot(data[:,0],data[:,1],label=r'$Simulation$',color='red', markersize=8,linewidth=2)
	plt.xlabel(r'$\alpha$')
	plt.ylabel(r'$nb \, / \, Nb_{tot}$')
	myTools.deal_with_plot(False,True,True)
	plt.show()
	del catData
	del catSimu
	
	### Distribution Beta
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	print ' data <beta> = ', numpy.mean(catData['BETA'])
	print ' simu <beta> = ', numpy.mean(catSimu['BETA'])
	data = myTools.GetHisto(catData['BETA'], numpy.arange(-0.4,0.4,0.001))
	#data[:,1] /= catData.size
	plt.plot(data[:,0],data[:,1],label=r'$Data$',color='blue', markersize=8,linewidth=2)
	data = myTools.GetHisto(catSimu['BETA'], numpy.arange(-0.4,0.4,0.001))
	#data[:,1] /= catSimu.size
	plt.plot(data[:,0],data[:,1],label=r'$Simulation$',color='red', markersize=8,linewidth=2)
	#data = myTools.GetHisto(catMock['BETA_2'], numpy.arange(-0.4,0.4,0.001))
	#data[:,1] /= catMock.size
	#plt.plot(data[:,0],data[:,1],label=r'$Mock$',color='green', markersize=8,linewidth=2)
	plt.xlabel(r'$\beta$')
	plt.ylabel(r'$nb \, / \, Nb_{tot}$')
	myTools.deal_with_plot(False,True,True)
	plt.show()
	del catData
	del catSimu
	#del catMock	

	### Distribution MEAN_FOREST_LAMBDA_RF
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	#catMock = pyfits.open(pathMocks,memmap=True)[1].data
	print '  catData = ', catData.size
	print '  catSimu = ', catSimu.size
	#print '  catMock = ', catMock.size
	plt.hist(catData['MEAN_FOREST_LAMBDA_RF'], bins=numpy.arange(1040.,1200.,1.),histtype='step',label=r'$Data$',color='blue',linewidth=2)
	plt.hist(catSimu['MEAN_FOREST_LAMBDA_RF'], bins=numpy.arange(1040.,1200.,1.),histtype='step',label=r'$Simulation$',color='red',linewidth=2)
	#plt.hist(catMock['MEAN_FOREST_LAMBDA_RF'], bins=numpy.arange(1040.,1200.,1.),histtype='step',label=r'$Mock$',color='green',linewidth=2)
	plt.xlabel(r'$\overline{\lambda_{R.F.}}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	del catData
	del catSimu
	#del catMock
	


	## Distribution redshift pixels
	catData = pyfits.open(pathData,memmap=True)[1].data[:3000]
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data[:3000]

	print '  catData = ', catData.size
	print '  catSimu = ', catSimu.size

	cut = (catData['NORM_FLUX_IVAR']>0.)
	zi  = catData['LAMBDA_OBS'][ cut ]/lambdaRFLine__ -1.
	we  = catData['DELTA_WEIGHT'][ cut ]
	hist = myTools.GetHisto(zi,numpy.arange(1.8,5.8,0.1),we)
	#hist[:,1] /= numpy.sum(hist[:,1])
	plt.plot(hist[:,0],hist[:,1],drawstyle='steps',label=r'$Data$',color='blue',linewidth=2)

	cut = (catSimu['NORM_FLUX_IVAR']>0.)
	zi  = catSimu['LAMBDA_OBS'][ cut ]/lambdaRFLine__ -1.
	we  = catSimu['DELTA_WEIGHT'][ cut ]
	hist = myTools.GetHisto(zi,numpy.arange(1.8,5.8,0.1),we)
	#hist[:,1] /= numpy.sum(hist[:,1])
	plt.plot(hist[:,0],hist[:,1],drawstyle='steps',label=r'$Simulation$',color='red',linewidth=2)

	plt.xlabel(r'$z_{Forest}$')
	#plt.ylabel(r'$Nb \, / \, integral$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
	del catData
	del catSimu


	distribSoverN()
	distribSoverN_delta()
	distribDelta()
	'''



	



	
	### Template
	data = numpy.loadtxt(path + 'template_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	meanData = []
	saveData = 1.
	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			#print i, j
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'template_LYA_'+str(i)+'_'+str(j)+'.txt')
				if (i==1 and j==1): continue
			except:
				#print ' ERROR: ', i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simulation$',color='red')
			if (i==0 and j==0):
				meanData = data[:,1]
				saveData = numpy.array(data[:,1])
			else: meanData += data[:,1]
	if (show_mean): meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$Simulation$',color='red')

	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$f(\lambda_{R.F.})$', fontsize=40) ##/f(1150.)
	myTools.deal_with_plot(False,False,True)
	plt.plot([lambdaRFMin__,lambdaRFMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.plot([lambdaRFMax__,lambdaRFMax__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.xlim( [lambdaRFMin__-10.,lambdaRFMax__+10.] )
	plt.show()
	
	
	### delta vs. lambda_RF
	data = numpy.loadtxt(path + 'deltaVSLambdaRF_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			#print i, j
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'deltaVSLambdaRF_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print ' ERROR: ', i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simulation$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	if (show_mean): meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$Mean \, Simulation$',color='red')

	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$\delta$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.plot([lambdaRFMin__,lambdaRFMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.plot([lambdaRFMax__,lambdaRFMax__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.xlim( [lambdaRFMin__-10.,lambdaRFMax__+10.] )
	plt.show()

	### delta+1 vs. lambda_RF
	data = numpy.loadtxt(path + 'hDeltaVsLambdaRF_LYA.txt')
	plt.errorbar(data[:,0]+1037.5, data[:,1], label=r'$Data$',color='blue', markersize=8,linewidth=4)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	meanData = []
	saveData = 1.
	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			#print i, j
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'hDeltaVsLambdaRF_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print ' ERROR: ', i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0]+1037.5, data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simu \, '+str(i)+' ' + str(j) + '$')
			if (i==0 and j==0):
				meanData = data[:,1]
				saveData = numpy.array(data[:,1])
			else: meanData += data[:,1]
	if (show_mean): meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0]+1037.5, meanData, markersize=8,linewidth=4, label=r'$Mean \, Simulation$',color='red')

	data = numpy.loadtxt(path + 'hDeltaVsLambdaRF_LYA.txt')
	plt.errorbar(data[:,0]+1037.5, data[:,1],color='blue', markersize=8,linewidth=4)

	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$\overline{C}(\lambda_{R.F.})$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.plot([lambdaRFMin__,lambdaRFMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.plot([lambdaRFMax__,lambdaRFMax__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.xlim( [lambdaRFMin__-10.,lambdaRFMax__+10.] )
	plt.show()
	
	
	### delta+1 vs. lambda_Obs

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twiny()

	data = numpy.loadtxt(path+'hDeltaVsLambdaObs_LYA.txt')
	ax1.errorbar(data[100:-100,0]+3500.5, data[100:-100,1], label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/PDF/hDeltaVsLambdaObs_LYA_JMC.txt')
	#print data[:,1]
	ax1.errorbar(data[:,1][ numpy.logical_and(data[:,2]!=0.,data[:,1]<=5589.5) ], data[:,2][ numpy.logical_and(data[:,2]!=0.,data[:,1]<=5589.5) ], label=r'$Mean \, Simulation \, Input$', color='orange',linewidth=2)

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			#print i, j
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'hDeltaVsLambdaObs_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print ' ERROR: ', i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): ax1.errorbar(data[:,0]+3500.5, data[:,1], markersize=8,linewidth=2) #, label=r'$Simulation$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	if (show_mean): meanData /= chunckNb*simulNb
	#print data[100:2090,0]+3500.5
	if (show_mean): ax1.errorbar(data[100:2090,0]+3500.5, meanData[100:2090], markersize=8,linewidth=2, label=r'$Mean \, Simulation \, Recovered$',color='red')  ##[100:2089]

	#plt.plot([lambdaObsMin__,lambdaObsMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	#plt.plot([7235.,7235.],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)

	ax1.set_xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	ax1.set_ylabel(r'$\overline{F}(\lambda_{Obs.})$', fontsize=40)
	myTools.deal_with_plot(False,False,False)

	
	new_tick_locations = numpy.append( numpy.arange(3600., 7235., 800.), [7235.] )
	def tick_function(X):
		V = X/1215.67-1.
		return ["%.2f" % z for z in V]
	def tick_function1(X):
		return ["%.0f" % z for z in X]
	ax1.set_xticks(new_tick_locations)
	ax1.set_xticklabels(tick_function1(new_tick_locations))
	ax2.set_xlim(ax1.get_xlim())
	ax2.set_xticks(new_tick_locations)
	ax2.set_xticklabels(tick_function(new_tick_locations))
	ax2.set_xlabel(r'$z_{pixel}$')

	ax1.grid(True, which='both')
	ax1.linewidth = 40

	myTools.deal_with_plot(False,False,False)
	ax1.legend(fontsize=40, numpoints=1,ncol=1)
	plt.subplots_adjust(left=0.1, right=0.9, top=0.88, bottom=0.12)
	
	plt.show()

	
	
	### delta vs. lambda_Obs
	data = numpy.loadtxt(path + 'deltaVSLambdaObs_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			#print i, j
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'deltaVSLambdaObs_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print ' ERROR: ', i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], markersize=8,linewidth=2) #, label=r'$Simulation$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	if (show_mean): meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, fmt='o', markersize=8,linewidth=2, label=r'$Simulation$',color='red')

	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$\delta$', fontsize=40)
	plt.plot([lambdaObsMin__,lambdaObsMin__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='green', markersize=8,linewidth=2)
	plt.plot([7235.,7235.],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='green', markersize=8,linewidth=2)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	### eta
	data = numpy.loadtxt(path + 'eta_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)

	#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_000/Simu_000/Run/eta_LYA_0_0.txt')
	#plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Simulation \, before$',color='orange', markersize=8,linewidth=2)

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			#print i, j
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'eta_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print ' ERROR: ', i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simulation$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	if (show_mean): meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$Simulation \, now$',color='red')

	plt.xlabel(r'$z_{pixel}$', fontsize=40)
	plt.ylabel(r'$\eta(z_{pixel})$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim( [1.9,4.7] )
	plt.show()
	
	
	### sigma
	data = numpy.loadtxt(path + 'sigma2LSS_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)

	#data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_000/Simu_000/Run/sigma2LSS_LYA_0_0.txt')
	#plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Simulation \, before$',color='orange', markersize=8,linewidth=2)

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			#print i, j
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'sigma2LSS_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print ' ERROR: ', i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simulation$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	if (show_mean): meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$Simulation \, now$',color='red')

	plt.xlabel(r'$z_{pixel}$', fontsize=40)
	plt.ylabel(r'$\sigma_{L.S.S.}^{2}(z_{pixel})$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim( [1.9,4.7] )
	plt.show()


def distribDelta():

	print
	print "------ Start ------"
	print

	name = [pathData,
		pathSimu+'delta.fits'
		]
	label = ['data','simu']

	saveHist = []
	saveVar  = []

	for i in numpy.arange(len(name)):
		print '\n\n'
		print label[i]
		### Data
		cat = pyfits.open(name[i],memmap=True)[1].data[:1000]
		print cat.size
	
		#cut = numpy.logical_and(numpy.logical_and( numpy.logical_and((cat['DELTA_IVAR'] > 0.), (cat['DELTA_WEIGHT']>0.)), numpy.logical_and((cat['NORM_FLUX_IVAR'] > 0.), (cat['FLUX_DLA']>=0.8)) ), numpy.logical_and( (cat['LAMBDA_RF']>=1040.), (cat['LAMBDA_RF']<1200.)))
		cut = numpy.logical_and( cat['DELTA_WEIGHT']>0., numpy.logical_and(cat['NORM_FLUX_IVAR'] > 0., cat['FLUX_DLA']>=0.8) )
		delta   = cat['DELTA'][cut].flatten()
		weights = cat['DELTA_WEIGHT'][cut].flatten()


		m = numpy.average(delta,weights=weights)
		v = numpy.average((delta-m)**2, weights=weights)

		print delta.size
		print m, v

		
		### Plot
		cut = numpy.abs(delta)<10.
		delta   = delta[cut]
		weights = weights[cut]

		print '  with cut'
		m = numpy.average(delta,weights=weights)
		v = numpy.average((delta-m)**2, weights=weights)
		print delta.size
		print m, v

		delta -= m
		saveVar += [v]
		

		hist, axisX = numpy.histogram(delta,bins=numpy.arange(-10.,10.,0.1),weights=weights)
	        xxx  = numpy.array([ axisX[j]+(axisX[j+1]-axisX[j])/2. for j in range(0,axisX.size-1) ])
		hist = numpy.asarray(zip(xxx,hist))
		hist[:,1] /= numpy.sum(hist[:,1])
		#hist[:,1] /= numpy.amax(hist[:,1])
		saveHist += [hist]

	print '\n\n'


        for i in numpy.arange(len(name)):
		plt.errorbar(saveHist[i][:,0],saveHist[i][:,1],fmt='o',label=label[i])
		print numpy.sum(saveHist[i][:,1][ numpy.abs(saveHist[i][:,0])<2.*saveVar[i] ]  )

	plt.grid(True, which='both')
	plt.xlabel(r'$\delta-<\delta>$', fontsize=40)
	plt.ylabel(r'$Nb/integral$', fontsize=40)
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	plt.yscale('log')
	plt.show()
	for i in numpy.arange(len(name)):
		plt.errorbar(saveHist[i][:,0],saveHist[i][:,1],label=label[i])
	plt.grid(True, which='both')
	plt.xlabel(r'$\delta-<\delta>$', fontsize=40)
	plt.ylabel(r'$Nb/integral$', fontsize=40)
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	plt.show()

	for i in numpy.arange(len(name)):
                plt.errorbar(saveHist[i][:,0],(saveHist[i][:,1]-saveHist[0][:,1])/saveHist[0][:,1],fmt='o',label=label[i])
	plt.grid(True, which='both')
        plt.xlabel(r'$\delta-<\delta>$', fontsize=40)
        plt.ylabel(r'$(nb_{i}-nb_{data})/nb_{data}$', fontsize=40)
        plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
        plt.show()

	return

def distribSoverN():


	print
	print "------ Start ------"
	print

	name = [pathData,
		pathSimu+'delta.fits',
		pathMocks,
		]
	label = ['Data','Simulation','Mock \, auto-correlation']

	saveHist = []

	for i in numpy.arange(len(name)):
		cat = pyfits.open(name[i],memmap=True)[1].data[:3000]
	
		#cut = numpy.logical_and(numpy.logical_and( cat['DELTA_WEIGHT']>0., numpy.logical_and((cat['NORM_FLUX_IVAR'] > 0.), (cat['FLUX_DLA']>=0.8)) ), numpy.logical_and( (cat['LAMBDA_RF']>=1040.), (cat['LAMBDA_RF']<1200.)))
		cut = numpy.logical_and( cat['DELTA_WEIGHT']>0., numpy.logical_and(cat['NORM_FLUX_IVAR'] > 0., cat['FLUX_DLA']>=0.8) )
		yyy = cat['NORM_FLUX'][cut].flatten()/numpy.power(cat['NORM_FLUX_IVAR'][cut].flatten(),-0.5)

		print '  mean S/N = ', numpy.mean(yyy)	

		hist, axisX = numpy.histogram(yyy,bins=numpy.arange(-10.,100.,0.1))
	        xxx  = numpy.array([ axisX[j]+(axisX[j+1]-axisX[j])/2. for j in range(0,axisX.size-1) ])
		hist = numpy.asarray(zip(xxx,hist))
		hist[:,1] /= numpy.sum(hist[:,1])
		saveHist += [hist]

	color = ['blue','red','green','black']
        for i in numpy.arange(len(name)):
		plt.plot(saveHist[i][:,0],saveHist[i][:,1],label=r'$'+label[i]+'$',color=color[i])

	plt.grid(True, which='both')
	plt.xlabel(r'$S/N \, = \, Flux/err$', fontsize=40)
	plt.ylabel(r'$Nb \, / \, integral$', fontsize=40)
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	plt.yscale('log')
	plt.show()

	return
def distribSoverN_delta():


	print
	print "------ Start ------"
	print

	name = [pathData,
		pathSimu+'delta.fits',
		pathMocks,
		]
	label = ['data','simu','Mock \, pipeline']

	saveHist = []

	for i in numpy.arange(len(name)):
		cat = pyfits.open(name[i],memmap=True)[1].data[:1000]
	
		#cut = numpy.logical_and(numpy.logical_and( numpy.logical_and((cat['DELTA_IVAR'] > 0.), (cat['DELTA_WEIGHT']>0.)), numpy.logical_and((cat['NORM_FLUX_IVAR'] > 0.), (cat['FLUX_DLA']>=0.8)) ), numpy.logical_and( (cat['LAMBDA_RF']>=1040.), (cat['LAMBDA_RF']<1200.)))
		cut = numpy.logical_and( cat['DELTA_WEIGHT']>0., numpy.logical_and(cat['NORM_FLUX_IVAR'] > 0., cat['FLUX_DLA']>=0.8) )
		template2 = cat['NORM_FLUX'][cut]/(cat['DELTA'][cut]+1.)
		delta_ivar = cat['NORM_FLUX_IVAR'][cut]*template2*template2
		yyy = cat['DELTA'][cut].flatten()/numpy.power(delta_ivar.flatten(),-0.5)		

		hist, axisX = numpy.histogram(yyy,bins=numpy.arange(-10.,100.,0.1))
	        xxx  = numpy.array([ axisX[j]+(axisX[j+1]-axisX[j])/2. for j in range(0,axisX.size-1) ])
		hist = numpy.asarray(zip(xxx,hist))
		hist[:,1] /= numpy.sum(hist[:,1])
		saveHist += [hist]

	color = ['blue','red','green','black']
        for i in numpy.arange(len(name)):
		plt.plot(saveHist[i][:,0],saveHist[i][:,1],label=r'$'+label[i]+'$',color=color[i])

	plt.grid(True, which='both')
	plt.xlabel(r'$S/N \, = \, \delta/\sigma_{pipeline} $', fontsize=40)
	plt.ylabel(r'$Nb \, / \, integral$', fontsize=40)
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	plt.yscale('log')
	plt.show()

	return

def compare_each_simu():

	pathSimu    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_000/Simu_000/Data/'
	rawPathSimu = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/'
	chunckNb = 10
	simulNb  = 10

	## Distribution redshift QSO
	nb_qso = numpy.zeros(chunckNb*simulNb)
	for i in range(0,chunckNb):
		for j in range(0,simulNb):
			catSimu = pyfits.open(rawPathSimu+'Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_withRSD.fits',memmap=True)[1].data
			nb_qso[i*10+j] = catSimu.size
	print 	numpy.mean(nb_qso)

	###
	plt.errorbar(numpy.arange(chunckNb*simulNb),nb_qso,fmt='o')
	plt.plot(numpy.arange(chunckNb*simulNb), numpy.ones(chunckNb*simulNb)*numpy.mean(nb_qso),color='red',label='Mean')
	plt.xlabel(r'$Mock \, index$', fontsize=40)
	plt.ylabel(r'$\# \, nb \, QSO$', fontsize=40)
	plt.xlim( [-1,chunckNb*simulNb] )
	myTools.deal_with_plot(False,False,True)
	plt.show()
	###
	plt.hist(nb_qso)
	plt.xlabel(r'$nb \, QSO$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()




	return


comparePlot()
#compare_each_simu()






















