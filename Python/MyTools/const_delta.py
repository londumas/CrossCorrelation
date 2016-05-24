# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import numpy
import cosmolopy.distance as cosmology

#forest__ = 'LYB'
#forest__ = 'LYB_LYA'
forest__ = 'LYA'
#forest__ = 'SIIV'
#forest__ = 'CIV'
#forest__ = 'MGII'

takeNormaZone__ = False





bit16__ = 65536



########################################################################
### Lambda_{Obs.}
########################################################################

### Observed wavelenght range
lambdaObsMin__ = 3600.
lambdaObsMax__ = 9350.  ##7235
log10lambdaObsMin__ = numpy.log10(lambdaObsMin__)
log10lambdaObsMax__ = numpy.log10(lambdaObsMax__)

### Array of edges for the meanTransFlux
lambdaObsBinEdges__ = numpy.arange(lambdaObsMin__, lambdaObsMax__+1., 1.)

### Bin Center for the meanTransFlux
lambdaObsBinCenter__  = numpy.array([ lambdaObsBinEdges__[i]+(lambdaObsBinEdges__[i+1]-lambdaObsBinEdges__[i])/2. for i in range(0,lambdaObsBinEdges__.size-1) ])

### Sky lines
#############
skyLines__ = [ (3615,3619),(3932,3937),(3966,3972),(4042,4050),
	(4357,4362),(5458,5467),(5573,5585),(5682,5695),
	(5885,5902),(6235,6241),(6256,6263),(6296,6311),
	(6320,6334),(6362,6369),(6498,6502),(6554,6557),
	(6825,6840),(6862,6870),(6922,6928),(6948,6954),
	(6977,6982),(9336,9344),(9369,9383),(9435,9449),
	(9460,9465),(9695,9705) ]
skyLinesNames__ = numpy.array(['?','Ca absorbtion','Ca absorbtion','Hg line','Hg line','Hg line','0I line Kirkby','?','?','?','?','O lines Kirkby','?','?','?','? Kirkby','?','?','?','?' ])
skyLines__ = numpy.log10(skyLines__)


########################################################################
### Lambda_{R.F.}
########################################################################

### Correlation lines
MGII_lines_names  = numpy.array(['MgI',   'MgII_a',  'MgII_b',  'FeII_a',      'FeII_b',      'MnII',   'FeII_c',      'FeII_d',      'FeII_e'])
MGII_lines        = numpy.array([2852.96, 2803.5324, 2796.3511, 2600.1724835,  2586.6495659,  2576.877, 2382.7641781,  2374.4603294,  2344.2129601])
lambda_RF_line_MGII = 2803.5324
###
CIV_lines_names   = numpy.append( MGII_lines_names, numpy.array(['AlIII(1863)' ,  'AlIII(1855)' , 'AlII(1671)' ,   'FeII(1609)',  'CIV(1551)' ,   'CIV(1548)',   'SiII(1527)']))
CIV_lines         = numpy.append( MGII_lines, numpy.array([      1862.79113,  1854.71829, 1670.7886, 1608.4511, 1550.77845, 1548.2049, 1526.70698 ]))
lambda_RF_line_CIV = 1548.2049
###
SIIV_lines_names  = numpy.append( CIV_lines_names, numpy.array(['SiIV(1403)',   'SiIV(1394)',   'CII(1335)' ,'SiII(1304)',    'OI(1302)',     'SiII(1260)',   'NV(1243)',   'NV(1239)']))
SIIV_lines        = numpy.append( CIV_lines, numpy.array([      1402.77291, 1393.76018, 1335., 1304.3702,   1302.1685,  1260.4221,  1242.804, 1238.821]))
lambda_RF_line_SIIV = 1393.76018
###
LYA_lines_names   = numpy.append( SIIV_lines_names, numpy.array(['LYA'   ,'SiIII(1207)',  'SiII(1193)',  'SiII(1190)']))
LYA_lines         = numpy.append( SIIV_lines, numpy.array([      1215.67 ,1206.500, 1193.2897, 1190.4158]))
lambda_RF_line_LYA = 1215.67
###
LYB_lines_names   = numpy.append( LYA_lines_names, numpy.array(['LYB','LY3','LY4','LY5']))
LYB_lines         = numpy.append( LYA_lines, numpy.array([1025.72, 972.537, 949.7431, 937.8035]))
lambda_RF_line_LYB = 1025.72


if (forest__ == 'LYB'):
	lambdaRFLine__     = 1025.72
	lambdaRFMin__      = 800.
	lambdaRFMax__      = 1020.
	lambdaRFNormaMin__ = 1050.
	lambdaRFNormaMax__ = 1060.
	nbBinRFMax__       = 1085
	correlationLines     = LYB_lines
	correlationLinesName = LYB_lines_names
	alphaStart__         = 1.
if (forest__ == 'LYB_LYA'):
        lambdaRFLine__     = 1215.67
        lambdaRFMin__      = 800.
        lambdaRFMax__      = 1200.
        lambdaRFNormaMin__ = 1275.
        lambdaRFNormaMax__ = 1295.
        nbBinRFMax__       = 1085
        correlationLines     = LYB_lines
        correlationLinesName = LYB_lines_names
        alphaStart__         = 1.
elif (forest__ == 'LYA'):
	lambdaRFLine__     = 1215.67
	lambdaRFMin__      = 1040.
	lambdaRFMax__      = 1200.
	lambdaRFNormaMin__ = 1275.
	lambdaRFNormaMax__ = 1295.
	nbBinRFMax__       = 647
	correlationLines     = LYA_lines
	correlationLinesName = LYA_lines_names
	alphaStart__         = 1.3
elif (forest__ == 'SIIV'):
	lambdaRFLine__     = 1393.76018
	lambdaRFMin__      = 1286.
	lambdaRFMax__      = 1380.
	lambdaRFNormaMin__ = 1415.
	lambdaRFNormaMax__ = 1425.
	nbBinRFMax__       = 326
	correlationLines     = SIIV_lines
	correlationLinesName = SIIV_lines_names
	alphaStart__         = 1.
elif (forest__ == 'CIV'):
	lambdaRFLine__     = 1548.2049
	lambdaRFMin__      = 1410.
	lambdaRFMax__      = 1530.
	lambdaRFNormaMin__ = 1600.
	lambdaRFNormaMax__ = 1630.
	nbBinRFMax__       = 373
	correlationLines     = CIV_lines
	correlationLinesName = CIV_lines_names
	alphaStart__         = 1.
elif (forest__ == 'MGII'):
	lambdaRFLine__     = 2803.5324
	lambdaRFMin__      = 1570.
	lambdaRFMax__      = 2790.
	lambdaRFNormaMin__ = 2860.
	lambdaRFNormaMax__ = 2880.
	nbBinRFMax__       = 2511
	correlationLines     = MGII_lines
	correlationLinesName = MGII_lines_names
	alphaStart__         = 1.



### Rest Frame wavelenght range for data

lambdaRFMean__ = (lambdaRFMax__+lambdaRFMin__)/2.

### Rest Frame wavelenght range for template
lambdaRFTemplateMin__ = lambdaRFMin__-3.
lambdaRFTemplateMax__ = lambdaRFMax__+3.
log10lambdaRFTemplateMin__ = numpy.log10(lambdaRFTemplateMin__)
log10lambdaRFTemplateMax__ = numpy.log10(lambdaRFTemplateMax__)

### Number of bin for template
#nbBinRFTemplate__ = int((lambdaRFTemplateMax__-lambdaRFTemplateMin__)/1.)

### Array of edges for the template
lambdaRFTemplateBinEdges__ = numpy.arange(lambdaRFTemplateMin__, lambdaRFTemplateMax__+1., 1.)

### Bin Center for the template
lambdaRFTemplateBinCenter__  = numpy.array([ lambdaRFTemplateBinEdges__[i]+(lambdaRFTemplateBinEdges__[i+1]-lambdaRFTemplateBinEdges__[i])/2. for i in range(0,lambdaRFTemplateBinEdges__.size-1) ])

### Normalisation region
log10lambdaRFNormaMin__ = numpy.log10(lambdaRFNormaMin__)
log10lambdaRFNormaMax__ = numpy.log10(lambdaRFNormaMax__)

### Min and max number of pixels in the forest
nbBinRFMin__ = 50


strongLines_CIV = [1.00166645,1.00522684,1.01407602,1.01576593,1.02024425,1.03718845,1.03891687,1.05354069,1.07544891,1.07738686,1.07918228]




########################################################################
### Redshift
########################################################################

### Minimal redshift to get first pixel
minRedshift__ = lambdaObsMin__/lambdaRFLine__ - 1.
maxRedshift__ = lambdaObsMax__/lambdaRFMin__ - 1.

redshiftBinEdges__ = numpy.array( [minRedshift__-0.01, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50,
	2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5,
	3.6, 3.7, 3.9, 4.1, 4.3, 4.5, maxRedshift__+0.01])

### Bin Center
redshiftBinCenter__  = numpy.array([ redshiftBinEdges__[i]+(redshiftBinEdges__[i+1]-redshiftBinEdges__[i])/2. for i in range(0,redshiftBinEdges__.size-1) ])





########################################################################
### Weight
########################################################################

### Number of pixel for the weight
nbBinWeight__ = 100

### Min of delta_Ivar_pipeline
minDeltaIvar__ = 1.e-3
maxDeltaIvar__ = 1.e6

### Array of bins
deltaIvarBinEdges__ = numpy.arange(numpy.log10(minDeltaIvar__),numpy.log10(maxDeltaIvar__)+0.1,0.1)

### Bin Center
deltaIvarBinCenter__  = numpy.array([ deltaIvarBinEdges__[i]+(deltaIvarBinEdges__[i+1]-deltaIvarBinEdges__[i])/2. for i in range(0,deltaIvarBinEdges__.size-1) ])

### Z_{0} is used to normalise
z0__        = 2.25
onePlusZ0__ = 1.+z0__
gama__      = 3.8
halfGama__  = gama__/2.

### Bin in redshift
### Array of bins
tmp_min = lambdaObsMin__/lambdaRFLine__-1.-0.01
tmp_max = lambdaObsMax__/lambdaRFLine__-1.+0.01

redshiftWeightBinEdges__  = numpy.arange(tmp_min,tmp_max+(tmp_max-tmp_min)/14., (tmp_max-tmp_min)/14.)

### Bin Center
redshiftWeightBinCenter__ = numpy.array([ redshiftWeightBinEdges__[i]+(redshiftWeightBinEdges__[i+1]-redshiftWeightBinEdges__[i])/2. for i in range(0,redshiftWeightBinEdges__.size-1) ])



########################################################################
### Cosmology
########################################################################

omegaM0__      = 0.27
omegaLambda0__ = 0.73
h__            = 0.7

########################################################################
### Targets
########################################################################

bitsQSO = numpy.genfromtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/ConfigFile/QTSflag_sdss4.dat',dtype='str')
dtype_out = [('KEY','a50'),('NAME','a50'),('BIT',numpy.int16)]
bitsQSO__ = numpy.ndarray(shape=bitsQSO[:,0].size,dtype=dtype_out)
bitsQSO__['KEY']  = bitsQSO[:,0]
bitsQSO__['NAME'] = bitsQSO[:,1]
bitsQSO__['BIT']  = bitsQSO[:,2].astype(int)

########################################################################
### Quasar
########################################################################

def find_range_QSO_redshift():
	'''
	'''

	verbose = False
	
	### range of the correlation
	maxX__ = 200.

	### Cosmology
	cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}

	### Distance max of correlation
	dist_max = maxX__*numpy.sqrt(2.)/h__
	if (verbose): print '  dist_max = ', dist_max

	z_pixel_min = (lambdaObsMin__/lambdaRFLine__)-1.
	z_pixel_max = (lambdaObsMax__/lambdaRFLine__)-1.

	if (verbose): print
	if (verbose): print '  z_pixel_min = ', z_pixel_min
	if (verbose): print '  z_pixel_max = ', z_pixel_max

	d_pixel_min = cosmology.comoving_distance(z_pixel_min, **cosmo)
	d_pixel_max = cosmology.comoving_distance(z_pixel_max, **cosmo)

	if (verbose): print
	if (verbose): print '  d_pixel_min = ', d_pixel_min
	if (verbose): print '  d_pixel_max = ', d_pixel_max

	d_QSO_min = d_pixel_min-dist_max
	d_QSO_max = d_pixel_max+dist_max

	if (verbose): print
	if (verbose): print '  d_QSO_min = ', d_QSO_min
	if (verbose): print '  d_QSO_max = ', d_QSO_max

	distfunc, redfunc = cosmology.quick_distance_function(cosmology.comoving_distance, zmax=20.0, zmin=0.0, zstep=0.001, return_inverse=True, **cosmo)
	d1 = distfunc(z_pixel_min)
	d2 = distfunc(z_pixel_max)
	if (verbose): print '\n',d1, d2

	z_QSO_min = redfunc(d_QSO_min)
	z_QSO_max = redfunc(d_QSO_max)

	if (verbose): print
	if (verbose): print '  z_QSO_min = ', z_QSO_min
	if (verbose): print '  z_QSO_max = ', z_QSO_max
	
	return z_QSO_min,z_QSO_max

'''
### Min and max redshift of the quasar
minRedshiftQSO__, maxRedshiftQSO__ = find_range_QSO_redshift()
print minRedshiftQSO__
print maxRedshiftQSO__
'''

########################################################################
### Pixels
########################################################################

def find_dist_pixels():
	'''
	'''

	verbose = True
	
	### Cosmology
	cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}

	z_pixel_min = (lambdaObsMin__/lambdaRFLine__)-1.
	z_pixel_max = (lambdaObsMax__/lambdaRFLine__)-1.

	#if (verbose): print
	#if (verbose): print '  z_pixel_min = ', z_pixel_min
	#if (verbose): print '  z_pixel_max = ', z_pixel_max

	d_pixel_min = cosmology.comoving_distance(z_pixel_min, **cosmo)*h__
	d_pixel_max = cosmology.comoving_distance(z_pixel_max, **cosmo)*h__

	if (verbose): print
	if (verbose): print '  d_pixel_min = ', d_pixel_min
	if (verbose): print '  d_pixel_max = ', d_pixel_max

	return d_pixel_min, d_pixel_max

########################################################################
### Pixels
########################################################################

def find_dist_correlation_lines(meanZ,l1,l2):
	'''
	'''

	verbose = False
	
	### Cosmology
	cosmo = {'omega_M_0':omegaM0__, 'omega_lambda_0':omegaLambda0__, 'omega_k_0':0., 'h':h__}

	lObs = (2.*meanZ+2.)/(1./l1+1./l2)
	z1 = lObs/l1-1.
	z2 = lObs/l2-1.

	if (verbose): print
	if (verbose): print '  z1 = ', z1
	if (verbose): print '  z2 = ', z2

	d1 = cosmology.comoving_distance(z1, **cosmo)*h__
	d2 = cosmology.comoving_distance(z2, **cosmo)*h__

	if (verbose): print
	if (verbose): print '  d1 = ',d1 
	if (verbose): print '  d2 = ',d2 

	return d1-d2














