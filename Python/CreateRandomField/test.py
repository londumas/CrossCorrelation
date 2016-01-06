# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >
#


import scipy
from scipy import interpolate
import numpy
import sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import myTools
import const_delta


def seePowerSpectrum():
	'''

	'''

	### Constants
	pathToPk = '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/CreateRandomField/DR9LyaMocks_matterpower.dat'

	### Get the power spectra
	pk = scipy.loadtxt(pathToPk)
	pkInterpolate = interpolate.interp1d(pk[:,0],pk[:,1],bounds_error=False,fill_value=0)

	for i in numpy.arange(3):
		coef = numpy.power(pk[:,0],i)
		yLabel = 'P(k)'
		if (i==1):
			yLabel = 'k.' + yLabel + '\, [h/Mpc]'
		elif (i==2):
			yLabel = 'k^{2}.' + yLabel + '\, [(h/Mpc)^{2}]'

		plt.plot(pk[:,0],coef*pk[:,1],marker='o')
		plt.xlabel(r'$k \, [h/Mpc]$', fontsize=40)
		plt.ylabel(r'$'+yLabel+'$', fontsize=40)
		myTools.deal_with_plot(False,False,False)
		plt.show()

	### Get the FT
	


def getListQSO():

	
	### Constants
	nameIn  = '4_5__withBiais_3_5'
	nameOut = 'box_512_512_512__size_4_5__sigma_1_5__replace_no__nbQSO_200000__properSTD_yes__biaisSet__3_5'
	std_4_5 = 1. #1.18045629669e-05
	nbQSO = 200000
	seed = 42
	
	pathToSave = 'QSO__'+nameOut+'.fits'
	pathToData = '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/CreateRandomField/delta_'+nameIn+'.fits'
	
	### Set the seed
	scipy.random.seed(seed)
	
	data = pyfits.open(pathToData, memmap=True )
	
	
	### Number pixel in box
	sizeCell = data[0].header['SIZE']
	sx = data[0].header['NX']
	sy = data[0].header['NY']
	sz = data[0].header['NZ']
	print sizeCell
	print '  Box size X: ',sx
	print '  Box size Y: ',sy
	print '  Box size Z: ',sz
	print '  Box number of cells = ', sx*sy*sz
	
	### Get data
	data = data[0].data
	print data.size
	
	mean = numpy.mean(data)
	print mean
	std = numpy.std(data)
	print std
	print std, std_4_5, std/std_4_5
	nu = 1.5 #*std/std_4_5
	print nu
	minDelta = nu*std
	print minDelta
	
	### Get index of bins
	ik = scipy.arange(sx*sy*sz)
	ik = scipy.reshape(ik,[sx,sy,sz])
	
	### Apply the biais
	ik = ik[ (data>minDelta) ]
	print ik.size
	print '  it represents (in %) ', 100.*ik.size/data.size
	
	### Get the index
	indexQSO_z = ik%sz
	indexQSO_y = ((ik-ik%sz)/sz)%sy
	indexQSO_x = (ik-ik%sz-(((ik-ik%sz)/sz)%sy)*sz)/sz/sy

	### Select randomly nbQSO
	if (ik.size>nbQSO):
		rand = numpy.random.choice(ik.size, nbQSO, replace=False)
		indexQSO_x = indexQSO_x[rand]
		indexQSO_y = indexQSO_y[rand]
		indexQSO_z = indexQSO_z[rand]
	
	### Look at the values of delta around a QSO
	deltaQSO = data[ (data>minDelta) ]
	deltaNeighboor = []
	for idx in numpy.arange(indexQSO_x.size):
		xxx = indexQSO_x[idx]
		yyy = indexQSO_y[idx]
		zzz = indexQSO_z[idx]
		deltaNeighboor += [ data[numpy.fmod(sx+xxx-1, sx)][yyy][zzz] ]
		deltaNeighboor += [ data[numpy.fmod(sx+xxx+1, sx)][yyy][zzz] ]
		deltaNeighboor += [ data[xxx][numpy.fmod(sy+yyy-1,sy)][zzz] ]
		deltaNeighboor += [ data[xxx][numpy.fmod(sy+yyy+1,sy)][zzz] ]
		deltaNeighboor += [ data[xxx][yyy][numpy.fmod(sz+zzz-1,sz)] ]
		deltaNeighboor += [ data[xxx][yyy][numpy.fmod(sz+zzz+1,sz)] ]

	deltaNeighboor = numpy.asarray(deltaNeighboor)

	### plot distrib neighboor and QSOs
	maxxx = numpy.amax(deltaQSO)
	bins = numpy.arange(-maxxx*1.1,maxxx*1.1,2.*maxxx/200)
	plt.hist(deltaQSO,bins=bins,histtype='step',label='QSO',log=True)
	plt.hist(deltaNeighboor,bins=bins,histtype='step',label='neighbor',log=True)
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2)
	plt.xlabel(r'$\delta$', fontsize=40)
	plt.ylabel(r'$\#$', fontsize=40)
	plt.show()
	
	### Show the gaussian
	plt.hist(data.flatten(),bins=500)
	plt.hist(data[ (data>minDelta) ].flatten(),bins=500)
	a = [data[indexQSO_x[i]][indexQSO_y[i]][indexQSO_z[i]] for i in numpy.arange(indexQSO_x.size) ]
	plt.hist(a,bins=500)
	plt.show()
	
	
	### Pass from index to distance
	xxx = indexQSO_x*sizeCell + sizeCell/2.
	yyy = indexQSO_y*sizeCell + sizeCell/2.
	zzz = indexQSO_z*sizeCell + sizeCell/2.
	
	
	### For a random position in the cell
	a = -sizeCell/2.
	b = sizeCell/2.
	ranPosXXX = (b - a) * numpy.random.random_sample() + a
	ranPosYYY = (b - a) * numpy.random.random_sample() + a
	ranPosZZZ = (b - a) * numpy.random.random_sample() + a
	xxx += ranPosXXX
	yyy += ranPosYYY
	zzz += ranPosZZZ
	
	
	print xxx.size
	### Save
	col_xx              = pyfits.Column(name='X', format='D', array=xxx, unit='Mpc.h^-1')
	col_yy              = pyfits.Column(name='Y', format='D', array=yyy, unit='Mpc.h^-1')
	col_zz              = pyfits.Column(name='Z', format='D', array=zzz, unit='Mpc.h^-1' )
		
	tbhdu = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz])
	tbhdu.header['NX'] = (sx,"number of pixels along the X direction")
	tbhdu.header['NY'] = (sy,"number of pixels along the Y direction")
	tbhdu.header['NZ'] = (sz,"number of pixels along the Z direction")
	tbhdu.header['SIZE'] = (sizeCell,"Size of a pixel [Mpc.h^-1]")
	
	tbhdu.writeto(pathToSave, clobber=True)
	
def haveAlookQSO():
	
	### Constants
	name = '4_5_big'

	pathToFile = 'QSO_'+name+'.fits'
	cat  = pyfits.open(pathToFile, memmap=True)[1].data

	print cat.size
	
	print numpy.amin(cat['X']), numpy.amax(cat['X'])
	print numpy.amin(cat['Y']), numpy.amax(cat['Y'])
	print numpy.amin(cat['Z']), numpy.amax(cat['Z'])
	print 1280*1280
	print cat.size
	print 1.*cat.size/(1280.*1280.)

	print

	### Map
	plt.plot(cat['X'],cat['Y'], linestyle="", marker="o")
	plt.xlabel(r'$X \, [h^{-1}.Mpc]$')
	plt.ylabel(r'$Y \, [h^{-1}.Mpc]$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution X
	plt.hist(cat['X'], bins=50)
	plt.xlabel(r'$X$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution Y
	plt.hist(cat['Y'], bins=50)
	plt.xlabel(r'$Y$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	### Distribution redshift
	plt.hist(cat['Z'], bins=50)
	plt.xlabel(r'$z_{QSO}$')
	plt.ylabel(r'$\#$')
	myTools.deal_with_plot(False,False,False)
	plt.show()
	
	return

def pickListQSO():
	
	### Constants
	nbQSO = 200000
	std_4_5 = 1.18045629669e-05  ### variance for a 4.5 Mpc box
	seed = 42
	pathToSave = 'QSO.fits' ### path to save
	pathToData = ''  ### path to the univers
	
	### Set the seed
	scipy.random.seed(seed)
	data = pyfits.open(pathToData, memmap=True )
	
	
	### Number pixel in box
	sizeCell = data[0].header['SIZE']
	sx = data[0].header['NX']
	sy = data[0].header['NY']
	sz = data[0].header['NZ']
	print sizeCell
	print '  Box size X: ',sx
	print '  Box size Y: ',sy
	print '  Box size Z: ',sz
	print '  Box number of cells = ', sx*sy*sz
	
	### Get data
	data = data[0].data
	print data.size
	
	mean = numpy.mean(data)
	print mean
	std = numpy.std(data)
	print std
	print std, std_4_5, std/std_4_5
	nu = 2.5*std/std_4_5
	print nu
	minDelta = nu*std
	print minDelta
	
	### Get index of bins
	ik = scipy.arange(sx*sy*sz)
	ik = scipy.reshape(ik,[sx,sy,sz])
	
	### Apply the biais
	ik = ik[ (data>minDelta) ]
	print ik.size
	print '  it represents (in %) ', 100.*ik.size/data.size
	
	### Get the index
	indexQSO_z = ik%sz
	indexQSO_y = ((ik-ik%sz)/sz)%sy
	indexQSO_x = (ik-ik%sz-(((ik-ik%sz)/sz)%sy)*sz)/sz/sy

	### Select randomly nbQSO in the available
	if (ik.size>nbQSO):
		rand = numpy.random.choice(ik.size, nbQSO, replace=False)
		indexQSO_x = indexQSO_x[rand]
		indexQSO_y = indexQSO_y[rand]
		indexQSO_z = indexQSO_z[rand]
	
	### Pass from index to distance
	xxx = indexQSO_x*sizeCell + sizeCell/2.
	yyy = indexQSO_y*sizeCell + sizeCell/2.
	zzz = indexQSO_z*sizeCell + sizeCell/2.
	
	print xxx.size

	col_xx              = pyfits.Column(name='X', format='D', array=xxx, unit='Mpc.h^-1')
	col_yy              = pyfits.Column(name='Y', format='D', array=yyy, unit='Mpc.h^-1')
	col_zz              = pyfits.Column(name='Z', format='D', array=zzz, unit='Mpc.h^-1' )
		
	tbhdu = pyfits.BinTableHDU.from_columns([col_xx, col_yy, col_zz])
	tbhdu.header['NX'] = (sx,"number of pixels along the X direction")
	tbhdu.header['NY'] = (sy,"number of pixels along the Y direction")
	tbhdu.header['NZ'] = (sz,"number of pixels along the Z direction")
	tbhdu.header['SIZE'] = (sizeCell,"Size of a pixel [Mpc.h^-1]")
	
	tbhdu.writeto(pathToSave, clobber=True)
	
	
#seePowerSpectrum()
getListQSO()
#haveAlookQSO()
	

'''
pathToData = '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/CreateRandomField/delta_4_5.fits'
data = pyfits.open(pathToData, memmap=True )[0].data[:100]
plt.hist(data.flatten(),bins=500,histtype='step')
pathToData = '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/CreateRandomField/delta_4_5__withBiais_3_5.fits'
data = pyfits.open(pathToData, memmap=True )[0].data[:100]
plt.hist(data.flatten(),bins=500,histtype='step')
plt.show()
'''





	
	
