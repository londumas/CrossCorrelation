#!/usr/bin/env python

import scipy
from scipy import interpolate
import numpy
import sys
import astropy.io.fits as pyfits
import psutil
import matplotlib.pyplot as plt

### Memory
memory = []
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]

pathToSave = 'delta_4_5__withBiais_3_5.fits'
seed = 42

### Number pixel in box
sx = 512
sy = 512
sz = 512

### Size cell [Mpc.h^-1]
sizeCell = 4.5*0.71


### Path to p(k)
#pathToPk = '/home/hdumasde/Documents/Program/baofit/models/DR9LyaMocks_matterpower.dat'
pathToPk = '/home/gpfs/manip/mnt0607/bao/hdumasde/Code/Python/CreateRandomField/DR9LyaMocks_matterpower.dat'

### Biais
bi = 3.5
### beta
be = 0.

















### Size box
lx = sizeCell*sx
ly = sizeCell*sy
lz = sizeCell*sz

print '  Box size X: ',sx,lx
print '  Box size Y: ',sy,ly
print '  Box size Z: ',sz,lz
print '  Box number of cells = ', sx*sy*sz
print '  sizeCell = ', sizeCell
print '  sizeCell/2. = ', sizeCell/2.


kx_max = 2*scipy.pi/lx*sx
ky_max = 2*scipy.pi/ly*sy
kz_max = 2*scipy.pi/lz*sz

dkx = 2*scipy.pi/lx
dky = 2*scipy.pi/ly
dkz = 2*scipy.pi/lz

dx = lx/sx
dy = ly/sy
dz = lz/sz


print '  k_x:  ', kx_max, dkx, dx
print '  k_y:  ', ky_max, dky, dy
print '  k_z:  ', kz_max, dkz, dz

scipy.random.seed(seed)

















### ik
print '  ik'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]

ik = scipy.arange(sx*sy*sz)

### kx
print '  kx'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]

kx                   = dkx*scipy.arange(sx)
memory += [ psutil.virtual_memory()[2] ]
kx[ (kx>kx_max/2) ] -= kx_max
memory += [ psutil.virtual_memory()[2] ]
k                    = kx[ ik%sx ]**2
memory += [ psutil.virtual_memory()[2] ]
del kx
memory += [ psutil.virtual_memory()[2] ]


### ky
print '  ky'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]

ky                   = dky*scipy.arange(sy)
memory += [ psutil.virtual_memory()[2] ]
ky[ (ky>ky_max/2) ] -= ky_max
memory += [ psutil.virtual_memory()[2] ]
k                   += ky[ ((ik-ik%sx)/sx)%sy ]**2
memory += [ psutil.virtual_memory()[2] ]
del ky
memory += [ psutil.virtual_memory()[2] ]


### kz
print '  kz'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]

kz                   = dkz*scipy.arange(sz)
memory += [ psutil.virtual_memory()[2] ]
kz[ (kz>kz_max/2) ] -= kz_max
memory += [ psutil.virtual_memory()[2] ]
k                   += kz[ (ik-ik%sx-(((ik-ik%sx)/sx)%sy)*sx)/sx/sy  ]**2
memory += [ psutil.virtual_memory()[2] ]
del kz, ik
memory += [ psutil.virtual_memory()[2] ]


### k
print '  k'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]
k = numpy.sqrt(k)
memory += [ psutil.virtual_memory()[2] ]
k = scipy.reshape(k,[sx,sy,sz])

### Get the P(k)
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]
pk = scipy.loadtxt(pathToPk)
memory += [ psutil.virtual_memory()[2] ]
pk = interpolate.interp1d(pk[:,0],pk[:,1],bounds_error=False,fill_value=0)


print '  Apply P(k) from BAOFIT'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]
k = pk(k)

k = numpy.sqrt(bi**2.*k)
memory += [ psutil.virtual_memory()[2] ]
del pk
memory += [ psutil.virtual_memory()[2] ]


print '  Initial delta field'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]
k = k*(scipy.random.randn(sx,sy,sz) + scipy.random.randn(sx,sy,sz)*1j)


print '  FFT'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]
k  = scipy.real(dkx*dky*dkz/(2*scipy.pi)**3*numpy.fft.fftn(k))



### Save
print '  Save'
print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]
k    = scipy.rec.array(k,dtype=scipy.float64)
k    = pyfits.PrimaryHDU(data=k)

k.header['NX'] = (sx,"number of pixels along the X direction")
k.header['NY'] = (sy,"number of pixels along the Y direction")
k.header['NZ'] = (sz,"number of pixels along the Z direction")
k.header['LX'] = (lx,"length along the X direction [Mpc.h^-1]")
k.header['LY'] = (ly,"length along the Y direction [Mpc.h^-1]")
k.header['LZ'] = (lz,"length along the Z direction [Mpc.h^-1]")
k.header['SIZE'] = (sizeCell,"Size of a pixel [Mpc.h^-1]")

print '  writing to ', pathToSave
k.writeto(pathToSave,clobber=True)

print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]

del k


print '  % = ', psutil.virtual_memory()[2], '\n'
memory += [ psutil.virtual_memory()[2] ]
print memory

