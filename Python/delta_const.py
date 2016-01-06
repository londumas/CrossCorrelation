# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import numpy

bit16__ = 65536




### Lambda_{R.F.}
########################################################################

### Observed wavelenght range
lambdaObsMin__ = 3600.
lambdaObsMax__ = 7235.
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
	(6977,6982) ]
skyLines__ = numpy.log10(skyLines__)







########################################################################
### Lambda_{R.F.}
########################################################################

### Rest Frame wavelenght range for data
lambdaRFMin__  = 1040.
lambdaRFMax__  = 1200.
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
lambdaRFNormaMin__ = 1275.
lambdaRFNormaMax__ = 1285.
log10lambdaRFNormaMin__ = numpy.log10(lambdaRFNormaMin__)
log10lambdaRFNormaMax__ = numpy.log10(lambdaRFNormaMax__)

### Lyman-alpha line
lambdaRFLine__ = 1215.67
### C-IV line
lambdaRFCIV__ = 1550.


### Min and max number of pixels in the forest
nbBinRFMin__ = 50
nbBinRFMax__ = 647










########################################################################
### Redshift
########################################################################

### Minimal redshift to get first pixel
minRedshift__ = lambdaObsMin__/lambdaRFMax__ - 1.
maxRedshift__ = lambdaObsMax__/lambdaRFNormaMin__ - 1.

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
h__            = 0.71


