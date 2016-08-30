import numpy
import myTools
import matplotlib.pyplot as plt

path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26_PlankCosmo/xi_1DlObs_2D_delta_delta_LYA.txt'

for i in range(10):
	for j in range(10):
		path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_1DlObs_2D_delta_delta_LYA.txt'
		if (i==0 and j==0): corr = numpy.loadtxt(path)
		else: corr += numpy.loadtxt(path)
corr /= 100.

###
x = numpy.arange(corr[:,0].size)+3600.
y = numpy.diag(corr)
cut = (y != 0.)
plt.plot( x[cut], y[cut])
plt.xlabel(r'$\lambda_{Obs.}$', fontsize=40)
plt.ylabel(r'$\xi^{1D}(0)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()
###
for i in [0,2]:
	y = numpy.diag(corr,k=i)
	x = numpy.arange(y.size)+3600.
	cut = y!=0.
	plt.plot( x[cut]/1215.67-1., y[cut] )
plt.plot( x[cut]/1215.67-1., numpy.ones(x[cut].size)*0.133177484547 )
plt.show()


###
for i in range(5):
	y = corr[:,i]
	x = numpy.arange(y.size)
	cut = y!=0.
	plt.plot( x[cut], y[cut] )
plt.show()

###
#corr[ (numpy.abs(corr)>0.015)] = numpy.float('nan')
corr[ corr==0. ] = numpy.float('nan')
myTools.plot2D(corr)




















