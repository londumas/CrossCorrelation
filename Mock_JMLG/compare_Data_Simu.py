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

pathMocks = '/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/mock.fits'


pathSimu    = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/Box_000/Simu_000/Data/'
rawPathSimu = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v_second_generation/'
chunckNb = 4
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
	print int((numpy.amax(catSimu['X'])-numpy.amin(catSimu['X']))/(4.5*0.71)) + 1
	print int((numpy.amax(catSimu['Y'])-numpy.amin(catSimu['Y']))/(4.5*0.71)) + 1
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
	plt.hist(catData['Z'], bins=numpy.arange(1.7,4.,0.1),histtype='step',label=r'$Data$',color='blue',linewidth=2)
	plt.hist(catSimu['Z'], bins=numpy.arange(1.7,4.,0.1),histtype='step',label=r'$Simu$',color='red',linewidth=2)
	plt.xlabel(r'$z_{QSO}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	del catData
	del catSimu
	
	## Distribution redshift Forest
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	#catMock = pyfits.open(pathMocks,memmap=True)[1].data
	print '  catData = ', catData.size
	print '  catSimu = ', catSimu.size
	#print '  catMock = ', catMock.size
	plt.hist(catData['Z_VI'], bins=numpy.arange(1.7,4.,0.1),histtype='step',label=r'$Data$',color='blue',linewidth=2)
	plt.hist(catSimu['Z_VI'], bins=numpy.arange(1.7,4.,0.1),histtype='step',label=r'$Simu$',color='red',linewidth=2)
	#plt.hist(catMock['Z_VI'], bins=numpy.arange(1.7,4.,0.1),histtype='step',label=r'$Mock$',color='green',linewidth=2)
	plt.xlabel(r'$z_{Forest}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	del catData
	del catSimu
	#del catMock

	### Distribution ALPHA
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	#catMock = pyfits.open(pathMocks,memmap=True)[1].data
	#print '  catData = ', catData.size, catData[ (catData['ALPHA_2']==1.) ].size
	#print '  catSimu = ', catSimu.size, catSimu[ (catSimu['ALPHA_2']==1.) ].size
	#print '  catMock = ', catMock.size, catMock[ (catMock['ALPHA_2']==1.) ].size
	data = myTools.GetHisto(catData['ALPHA_2'], numpy.arange(-50.,50.,0.01))
	#data[:,1] /= catData.size
	plt.plot(data[:,0],data[:,1],label=r'$Data$',color='blue', markersize=8,linewidth=2)
	data = myTools.GetHisto(catSimu['ALPHA_2'], numpy.arange(-50.,50.,0.01))
	#data[:,1] /= catSimu.size
	plt.plot(data[:,0],data[:,1],label=r'$Simu$',color='red', markersize=8,linewidth=2)
	#data = myTools.GetHisto(catMock['ALPHA_2'], numpy.arange(-50.,50.,0.01))
	#data[:,1] /= catMock.size
	#plt.plot(data[:,0],data[:,1],label=r'$Mock$',color='green', markersize=8,linewidth=2)
	plt.xlabel(r'$\alpha$')
	plt.ylabel(r'$nb \, / \, Nb_{tot}$')
	myTools.deal_with_plot(False,True,True)
	plt.show()
	del catData
	del catSimu
	#del catMock
	
	### Distribution Beta
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	#catMock = pyfits.open(pathMocks,memmap=True)[1].data
	#print '  catData = ', catData.size, catData[ (catData['BETA_2']==0.) ].size
	#print '  catSimu = ', catSimu.size, catSimu[ (catSimu['BETA_2']==0.) ].size
	#print '  catMock = ', catMock.size, catMock[ (catMock['BETA_2']==0.) ].size
	data = myTools.GetHisto(catData['BETA_2'], numpy.arange(-0.4,0.4,0.001))
	#data[:,1] /= catData.size
	plt.plot(data[:,0],data[:,1],label=r'$Data$',color='blue', markersize=8,linewidth=2)
	data = myTools.GetHisto(catSimu['BETA_2'], numpy.arange(-0.4,0.4,0.001))
	#data[:,1] /= catSimu.size
	plt.plot(data[:,0],data[:,1],label=r'$Simu$',color='red', markersize=8,linewidth=2)
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
	

	### Distribution NB_PIXEL
	catData = pyfits.open(pathData,memmap=True)[1].data
	catSimu = pyfits.open(pathSimu+'delta.fits',memmap=True)[1].data
	#catMock = pyfits.open(pathMocks,memmap=True)[1].data
	print '  catData = ', catData.size
	print '  catSimu = ', catSimu.size
	#print '  catMock = ', catMock.size
	plt.hist(catData['NB_PIXEL'], bins=numpy.arange(0.,700.,4.),histtype='step',label=r'$Data$',color='blue',linewidth=2)
	plt.hist(catSimu['NB_PIXEL'], bins=numpy.arange(0.,700.,4.),histtype='step',label=r'$Simu$',color='red',linewidth=2)
	#plt.hist(catMock['NB_PIXEL'], bins=numpy.arange(0.,700.,4.),histtype='step',label=r'$Mock$',color='green',linewidth=2)
	plt.xlabel(r'$number \, pixel$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
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
	plt.hist(catSimu['MEAN_FOREST_LAMBDA_RF'], bins=numpy.arange(1040.,1200.,1.),histtype='step',label=r'$Simu$',color='red',linewidth=2)
	#plt.hist(catMock['MEAN_FOREST_LAMBDA_RF'], bins=numpy.arange(1040.,1200.,1.),histtype='step',label=r'$Mock$',color='green',linewidth=2)
	plt.xlabel(r'$< \lambda_{R.F.} >$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	del catData
	del catSimu
	#del catMock
	
	
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
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'template_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simu$',color='red')
			if (i==0 and j==0):
				meanData = data[:,1]
				saveData = numpy.array(data[:,1])
			else: meanData += data[:,1]
	meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$mean \, simu$',color='red')

	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$f(\lambda_{R.F.})$', fontsize=40) ##/f(1150.)
	myTools.deal_with_plot(False,False,True)
	plt.plot([lambdaRFMin__,lambdaRFMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.plot([lambdaRFMax__,lambdaRFMax__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.show()
	
	
	### delta vs. lambda_RF
	data = numpy.loadtxt(path + 'deltaVSLambdaRF_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'deltaVSLambdaRF_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simu$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$mean \, simu$',color='red')

	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$\delta$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.plot([lambdaRFMin__,lambdaRFMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.plot([lambdaRFMax__,lambdaRFMax__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.show()

	### delta+1 vs. lambda_RF
	data = numpy.loadtxt(path + 'hDeltaVsLambdaRF_LYA.txt')
	plt.errorbar(data[:,0]+1037.5, data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	meanData = []
	saveData = 1.
	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'hDeltaVsLambdaRF_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0]+1037.5, data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simu \, '+str(i)+' ' + str(j) + '$')
			if (i==0 and j==0):
				meanData = data[:,1]
				saveData = numpy.array(data[:,1])
			else: meanData += data[:,1]
	meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0]+1037.5, meanData, marker='o', markersize=8,linewidth=2, label=r'$mean \, simu$',color='red')

	plt.xlabel(r'$\lambda_{R.F.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$\delta+1$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.plot([lambdaRFMin__,lambdaRFMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.plot([lambdaRFMax__,lambdaRFMax__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	plt.show()
	
	
	### delta+1 vs. lambda_Obs
	data = numpy.loadtxt(path+'hDeltaVsLambdaObs_LYA.txt')
	plt.errorbar(data[100:-100,0]+3500.5, data[100:-100,1], label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'hDeltaVsLambdaObs_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0]+3500.5, data[:,1], markersize=8,linewidth=2) #, label=r'$Simu$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[100:2089,0]+3500.5, meanData[100:2089], markersize=8,linewidth=2, label=r'$mean \, simu$',color='red')

	#plt.plot([lambdaObsMin__,lambdaObsMin__],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)
	#plt.plot([7235.,7235.],[min(numpy.min(data[:,1]),minA),max(numpy.max(data[:,1]),maxA)],color='green', markersize=8,linewidth=2)

	data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/hDeltaVsLambdaObs_LYA_JMC.txt')
	plt.errorbar(data[:,1][ data[:,2]!=0. ], data[:,2][ data[:,2]!=0. ], label=r'$Simu \, input$', color='orange',linewidth=2)
	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$f(\lambda_{Obs.})$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	### delta vs. lambda_Obs
	data = numpy.loadtxt(path + 'deltaVSLambdaObs_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], fmt='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)
	minA = numpy.min(data[:,1])
	maxA = numpy.max(data[:,1])

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'deltaVSLambdaObs_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], markersize=8,linewidth=2) #, label=r'$Simu$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, fmt='o', markersize=8,linewidth=2, label=r'$mean \, simu$',color='red')

	plt.xlabel(r'$\lambda_{Obs.} \, [\AA]$', fontsize=40)
	plt.ylabel(r'$\delta$', fontsize=40)
	plt.plot([lambdaObsMin__,lambdaObsMin__],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='green', markersize=8,linewidth=2)
	plt.plot([7235.,7235.],[numpy.min(data[:,1]),numpy.max(data[:,1])],color='green', markersize=8,linewidth=2)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	### eta
	data = numpy.loadtxt(path + 'eta_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'eta_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simu$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$mean \, simu$',color='red')

	plt.xlabel(r'$z_{pixel}$', fontsize=40)
	plt.ylabel(r'$\eta(z_{pixel})$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	
	### sigma
	data = numpy.loadtxt(path + 'sigma2LSS_LYA.txt')
	plt.errorbar(data[:,0], data[:,1], marker='o', label=r'$Data$',color='blue', markersize=8,linewidth=2)

	for i in range (0,chunckNb):
		for j in range(0,simulNb):
			tmp_path = rawPathSimu + 'Box_00'+str(i)+'/Simu_00'+str(j)+'/Run/'
			try:
				data = numpy.loadtxt(tmp_path + 'sigma2LSS_LYA_'+str(i)+'_'+str(j)+'.txt')
			except:
				print i, j
				continue
			if (len(data)==0):
				print i, j
				continue
			if (not show_mean): plt.errorbar(data[:,0], data[:,1], marker='o', markersize=8,linewidth=2) #, label=r'$Simu$',color='red')
			if (i==0 and j==0): meanData = data[:,1]
			else: meanData += data[:,1]
	meanData /= chunckNb*simulNb
	if (show_mean): plt.errorbar(data[:,0], meanData, marker='o', markersize=8,linewidth=2, label=r'$mean \, simu$',color='red')

	plt.xlabel(r'$z_{pixel}$', fontsize=40)
	plt.ylabel(r'$\sigma_{L.S.S.}^{2}(z_{pixel})$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
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
	
		cut = numpy.logical_and(numpy.logical_and( numpy.logical_and((cat['DELTA_IVAR'] > 0.), (cat['DELTA_WEIGHT']>0.)), numpy.logical_and((cat['NORM_FLUX_IVAR'] > 0.), (cat['FLUX_DLA']>=0.8)) ), numpy.logical_and( (cat['LAMBDA_RF']>=1040.), (cat['LAMBDA_RF']<1200.)))
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
	label = ['data','simu','Mock \, pipeline']

	saveHist = []

	for i in numpy.arange(len(name)):
		cat = pyfits.open(name[i],memmap=True)[1].data[:1000]
	
		cut = numpy.logical_and(numpy.logical_and( numpy.logical_and((cat['DELTA_IVAR'] > 0.), (cat['DELTA_WEIGHT']>0.)), numpy.logical_and((cat['NORM_FLUX_IVAR'] > 0.), (cat['FLUX_DLA']>=0.8)) ), numpy.logical_and( (cat['LAMBDA_RF']>=1040.), (cat['LAMBDA_RF']<1200.)))
		yyy = cat['NORM_FLUX'][cut].flatten()/numpy.power(cat['NORM_FLUX_IVAR'][cut].flatten(),-0.5)		

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
		cat = pyfits.open(name[i],memmap=True)[1].data[:10000]
	
		cut = numpy.logical_and(numpy.logical_and( numpy.logical_and((cat['DELTA_IVAR'] > 0.), (cat['DELTA_WEIGHT']>0.)), numpy.logical_and((cat['NORM_FLUX_IVAR'] > 0.), (cat['FLUX_DLA']>=0.8)) ), numpy.logical_and( (cat['LAMBDA_RF']>=1040.), (cat['LAMBDA_RF']<1200.)))
		yyy = cat['DELTA'][cut].flatten()/numpy.power(cat['DELTA_IVAR'][cut].flatten(),-0.5)		

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






















