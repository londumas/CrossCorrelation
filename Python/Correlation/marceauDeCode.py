####
xi1D_ss, xi2D_ss, xiMu_ss, xiWe_ss = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1519/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1519/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
xi1D_dd, xi2D_dd, xiMu_dd, xiWe_dd = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_QSO.txt')


for i in range(0,3):
	coef = numpy.power(xi1D_ss[:,0],i)
	plt.errorbar(xi1D_ss[:,0], coef*(xi1D_ss[:,1]-xi1D_ss[-1,1]), yerr=coef*xi1D_ss[:,2], fmt='o',label='simulation',  markersize=8, elinewidth=6, color='blue')
	plt.errorbar(xi1D_dd[:,0], coef*(xi1D_dd[:,1]-xi1D_dd[-1,1]), yerr=coef*xi1D_dd[:,2], fmt='o',label='data',  markersize=8, elinewidth=6, color='red')
	plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)

	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
	plt.show()



plt.errorbar(xi1D_ss[:,0], (xi1D_ss[:,1])/(xi1D_dd[:,1]), fmt='o',label='data / simulation',  markersize=8, elinewidth=6, color='blue')
plt.ylabel(r'$\xi^{qf} \, \xi^{qf}$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
plt.show()






#############################################################################
a = ['QSO','ALL_BRITT','CIV_BRITT','ALLNOTQSO_BRITT','CIVNOTQSO_BRITT','CIVLOWEW_BRITT', 'CIVHIGHEW_BRITT' ]
a = ['QSO','CIV_BRITT','CIVLOWEW_BRITT', 'CIVHIGHEW_BRITT' ]
###
for el in a:
	xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_'+el+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_'+el+'.txt')
	coef = 1.
	plt.errorbar(xi1D_[:,0], coef*(xi1D_[:,1]-xi1D_[-1,1]), yerr=coef*xi1D_[:,2],label=el, markersize=5, marker='o')

plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([min1D__,max1D__])
plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
plt.show()
###=


xi1Dboot_  = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/shuffleForest_1D.npy')
xi1Dboot2_ = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/bootstrap_1D.npy')
print xi1Dboot_[0].size
xi2Dboot_  = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/shuffleForest_2D.npy')
xi2Dboot2_ = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/bootstrap_2D.npy')
print xi1Dboot_[0].size

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')


cov = numpy.cov(xi1Dboot_)
#cov = numpy.cov(xi1Dboot2_)/xi1Dboot2_[0].size
#cov = plotCovar([numpy.cov(xi1Dboot_),numpy.cov(xi1Dboot2_)/xi1Dboot2_[0].size], ['Schuffle forest','Subsampling'])[1]
tmp_xi2D        = xi1D_[:,1]
tmp_idx         = numpy.arange(0,nbBin1D__)
nbBin1D = nbBin1D__*(nbBin1D__+1.)/2.
tmp_xi2DCov     = numpy.zeros(nbBin1D)
tmp_xi2DCovIdx  = numpy.zeros(nbBin1D)
tmp_xi2DCovIdx2 = numpy.zeros(nbBin1D)
idx = 0


for k1 in range(0,nbBin1D__):
	for k2 in range(0,k1+1):
		#if (k1==k2): tmp_xi2DCov[idx]     = cov[k1][k2]
		#if (k1==k2): tmp_xi2DCov[idx]     = xi1D_[k1,2]**2.
		tmp_xi2DCov[idx]     = cov[k1][k2]
		tmp_xi2DCovIdx[idx]  = k1
		tmp_xi2DCovIdx2[idx] = k2
		idx += 1

numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/BAOFIT/bao1D.data',zip(tmp_idx,tmp_xi2D),fmt='%u %1.20e')
numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/BAOFIT/bao1D.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%u %u %1.20e')

data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_1D_LYA_JMC_QSO_MockJmc.txt')

save2 = data[:,2]
save3 = data[:,3]
save4 = data[:,4]
			
meanRperp  = numpy.zeros(shape=(nbBin1D__,2))
meanRedshift = numpy.zeros(2)
			
for i in range(0,len(save2)):
	iX = i/binSize__

	meanRperp[iX][0]  += save2[i]
	meanRperp[iX][1]  += save4[i]
meanRedshift[0] = numpy.sum(save3)
meanRedshift[1] = numpy.sum(save4)
			
meanRperp[:,0]  /= meanRperp[:,1]
			
### Get the s_perp bin center
stringRperp = ''
for el in meanRperp[:-1,0]:
	stringRperp += str(el) + ','
stringRperp += str(meanRperp[-1,0])
print stringRperp

### Get the redshift center
stringRedshift = str(meanRedshift[0]/meanRedshift[1])
print ' <z> = ', stringRedshift

cov = numpy.cov(xi2Dboot_)
tmp_xi2D        = numpy.zeros(nbBin2D__)
tmp_idx         = numpy.arange(0,nbBin2D__)
nbBin2D = nbBin2D__*(nbBin2D__+1.)/2.
tmp_xi2DCov     = numpy.zeros(nbBin2D)
tmp_xi2DCovIdx  = numpy.zeros(nbBin2D)
tmp_xi2DCovIdx2 = numpy.zeros(nbBin2D)
tmp2_xi2DCov     = numpy.zeros(nbBin2D)
tmp2_xi2DCovIdx  = numpy.zeros(nbBin2D)
tmp2_xi2DCovIdx2 = numpy.zeros(nbBin2D)
idx = 0


for k1 in range(0,nbBin2D__):
	i1       = k1/nbBinY2D__
	j1       = k1%nbBinY2D__
	k11      = j1*nbBinX2D__ + i1
	tmp_xi2D[k11] = xi2D_[i1][j1][1]
		
	for k2 in range(0,k1+1):
		i2       = k2/nbBinY2D__
		j2       = k2%nbBinY2D__
		k22      = j2*nbBinX2D__ + i2
		if (k1==k2): tmp_xi2DCov[idx]     = cov[k1][k2]
		tmp_xi2DCovIdx[idx]  = k11
		tmp_xi2DCovIdx2[idx] = k22
		'''
		tmp2_xi2DCov[idx]     = cov[k1][k2]
		tmp2_xi2DCovIdx[idx]  = k11
		tmp2_xi2DCovIdx2[idx] = k22
		'''
		idx += 1
		
numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/BAOFIT/bao2D.data',zip(tmp_idx,tmp_xi2D),fmt='%u %1.20e')
numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/BAOFIT/bao2D.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%u %u %1.20e')



















#############################################################################
a = ['QSO','ALL_BRITT','CIV_BRITT','ALLNOTQSO_BRITT','CIVNOTQSO_BRITT']
###
for el in a:
	xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_'+el+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_'+el+'.txt')
	coef = 1.
	plt.errorbar(xi1D_[:,0], coef*(xi1D_[:,1]-xi1D_[-1,1]), yerr=coef*xi1D_[:,2], fmt='o',label=el)

plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([min1D__,max1D__])
plt.show()
###
for el in a:
	xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_'+el+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_'+el+'.txt')
	coef = xi1D_[:,0]
	plt.errorbar(xi1D_[:,0], coef*(xi1D_[:,1]-xi1D_[-1,1]), yerr=coef*xi1D_[:,2], fmt='o',label=el)
	plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([min1D__,max1D__])
plt.show()
###
for el in a:
	xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_'+el+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_'+el+'.txt')
	coef = xi1D_[:,0]**2.
	plt.errorbar(xi1D_[:,0], coef*(xi1D_[:,1]-xi1D_[-1,1]), yerr=coef*xi1D_[:,2], fmt='o',label=el)
	plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([min1D__,max1D__])
plt.show()
plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)

'''
xi1Dboot_  = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/shuffleForest_1D.npy')
xi2Dboot_  = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/shuffleForest_2D.npy')
print xi2Dboot_[0].size
xi1Dboot2_ = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/bootstrap_1D.npy')
xi2Dboot2_ = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/bootstrap_2D.npy')
print xi2Dboot2_[0].size

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
mean = numpy.mean(xi1Dboot_,axis=1)

###
plt.errorbar(numpy.arange(0,mean.size), xi1D_[:,1], yerr=xi1D_[:,2], marker='o', label=r'$Simul$')
plt.errorbar(numpy.arange(0,mean.size), mean, yerr=numpy.std(xi1Dboot_,axis=1)/numpy.sqrt(xi1Dboot_[0].size), marker='o', label=r'$schuffle$')
plt.errorbar(numpy.arange(0,mean.size), numpy.mean(xi1Dboot2_,axis=1), marker='o', label=r'$subsampling$')
myTools.deal_with_plot(False,False,True)
plt.show()
mean = numpy.mean(xi2Dboot_,axis=1)
plt.errorbar(numpy.arange(0,mean.size), xi2D_[:,:,1].flatten(), yerr=xi2D_[:,:,2].flatten(), marker='o', label=r'$Simul$')
plt.errorbar(numpy.arange(0,mean.size), mean, yerr=numpy.std(xi2Dboot_,axis=1)/numpy.sqrt(xi2Dboot_[0].size), marker='o', label=r'$schuffle$')
plt.errorbar(numpy.arange(0,mean.size), numpy.mean(xi2Dboot2_,axis=1), marker='o', label=r'$subsampling$')
myTools.deal_with_plot(False,False,True)
plt.show()
###
mean = numpy.diagonal(numpy.cov(xi1Dboot_))
plt.errorbar(numpy.arange(0,mean.size), mean, marker='o', label=r'$schuffle$')
plt.errorbar(numpy.arange(0,mean.size), numpy.diagonal(numpy.cov(xi1Dboot2_))/xi2Dboot2_[0].size, marker='o', label=r'$subsampling$')
plt.errorbar(numpy.arange(0,mean.size), xi1D_[:,2]**2., marker='o', label=r'$Simul$')
myTools.deal_with_plot(False,False,True)
plt.show()
mean = numpy.diagonal(numpy.cov(xi2Dboot_))
plt.errorbar(numpy.arange(0,mean.size), mean, marker='o', label=r'$schuffle$')
plt.errorbar(numpy.arange(0,mean.size), numpy.diagonal(numpy.cov(xi2Dboot2_))/xi2Dboot2_[0].size, marker='o', label=r'$subsampling$')
plt.errorbar(numpy.arange(0,mean.size), xi2D_[:,:,2].flatten()**2., marker='o', label=r'$Simul$')
myTools.deal_with_plot(False,False,True)
plt.show()
###
mean = numpy.diagonal(numpy.cov(xi1Dboot_))
plt.errorbar(numpy.arange(0,mean.size), mean/(xi1D_[:,2]**2.), marker='o', label=r'$schuffle$')
plt.errorbar(numpy.arange(0,mean.size), numpy.diagonal(numpy.cov(xi1Dboot2_))/xi2Dboot2_[0].size/(xi1D_[:,2]**2.), marker='o', label=r'$subsampling$')
myTools.deal_with_plot(False,False,True)
plt.show()
mean = numpy.diagonal(numpy.cov(xi2Dboot_))
plt.errorbar(numpy.arange(0,mean.size), mean/(xi2D_[:,:,2].flatten()**2.), marker='o', label=r'$schuffle$')
plt.errorbar(numpy.arange(0,mean.size), numpy.diagonal(numpy.cov(xi2Dboot2_))/xi2Dboot2_[0].size/(xi2D_[:,:,2].flatten()**2.), marker='o', label=r'$subsampling$')
myTools.deal_with_plot(False,False,True)
plt.show()
###
print numpy.mean(numpy.diagonal(numpy.cov(xi1Dboot_))/(xi1D_[:,2]**2.))
print numpy.mean(numpy.diagonal(numpy.cov(xi1Dboot2_))/xi1Dboot2_[0].size/(xi1D_[:,2]**2.))
print numpy.mean(xi1D_[:,2]**2.)
print numpy.mean(numpy.diagonal(numpy.cov(xi1Dboot_)))
print numpy.mean(numpy.diagonal(numpy.cov(xi1Dboot2_))/xi1Dboot2_[0].size)
print
print numpy.mean(numpy.diagonal(numpy.cov(xi2Dboot_))/(xi2D_[:,:,2].flatten()**2.))
print numpy.mean(numpy.diagonal(numpy.cov(xi2Dboot2_))/xi2Dboot2_[0].size/(xi2D_[:,:,2].flatten()**2.))
print numpy.mean(xi2D_[:,:,2].flatten()**2.)
print numpy.mean(numpy.diagonal(numpy.cov(xi2Dboot_)))
print numpy.mean(numpy.diagonal(numpy.cov(xi2Dboot2_))/xi2Dboot2_[0].size)
###
#plotCovar([numpy.cov(xi1Dboot_),numpy.cov(xi1Dboot2_)/xi2Dboot2_[0].size], ['Schuffle forest','Subsampling'])
#covFromFit = plotCovar([numpy.cov(xi2Dboot_),numpy.cov(xi2Dboot2_)/xi2Dboot2_[0].size], ['Schuffle forest','Subsampling'])
#numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/shuffleForest_subsampling_2D.npy',covFromFit)

cov  = numpy.cov(xi1Dboot_)
icov = numpy.linalg.inv(cov)
'''

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
xi2DCov  = numpy.cov(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/shuffleForest_subsampling_2D.npy')[1]) ##bootstrap
ixi2DCov = numpy.linalg.inv(xi2DCov)

tmp_xi2D        = numpy.zeros(nbBin2D__)
tmp_idx         = numpy.arange(0,nbBin2D__)
nbBin2D = nbBin2D__*(nbBin2D__+1.)/2.
tmp_xi2DCov     = numpy.zeros(nbBin2D)
tmp_xi2DCovIdx  = numpy.zeros(nbBin2D)
tmp_xi2DCovIdx2 = numpy.zeros(nbBin2D)
tmp2_xi2DCov     = numpy.zeros(nbBin2D)
tmp2_xi2DCovIdx  = numpy.zeros(nbBin2D)
tmp2_xi2DCovIdx2 = numpy.zeros(nbBin2D)
idx = 0


for k1 in range(0,nbBin2D__):
	i1       = k1/nbBinY2D__
	j1       = k1%nbBinY2D__
	k11      = j1*nbBinX2D__ + i1
	tmp_xi2D[k11] = xi2D_[i1][j1][1]
		
	for k2 in range(0,k1+1):
		i2       = k2/nbBinY2D__
		j2       = k2%nbBinY2D__
		k22      = j2*nbBinX2D__ + i2
		tmp_xi2DCov[idx]     = xi2DCov[k1][k2]
		tmp_xi2DCovIdx[idx]  = k11
		tmp_xi2DCovIdx2[idx] = k22
		tmp2_xi2DCov[idx]     = ixi2DCov[k1][k2]
		tmp2_xi2DCovIdx[idx]  = k11
		tmp2_xi2DCovIdx2[idx] = k22
		idx += 1
		
numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/BAOFIT/bao2D.data',zip(tmp_idx,tmp_xi2D),fmt='%d %E')
numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/BAOFIT/bao2D.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%d %d %E')
numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/BAOFIT/bao2D.icov',zip(tmp2_xi2DCovIdx,tmp2_xi2DCovIdx2,tmp2_xi2DCov),fmt='%d %d %E')


'''
for ii in range(0,10):
	for jj in range(0,10):
		path = str(ii)+'_'+str(jj)

		xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_DLA_mocks__'+path+'_.txt', path1__ +'xi_delta_QSO_2D_LYA_DLA_mocks__'+path+'_.txt')
		xi2DCov  = numpy.load(path1__ + 'mockDLA_'+path+'_subsampling_fromFit_2D.npy')[0]

		tmp_xi2D        = numpy.zeros(nbBin2D__)
		tmp_idx         = numpy.arange(0,nbBin2D__)
		nbBin2D = nbBin2D__*(nbBin2D__+1.)/2.
		tmp_xi2DCov     = numpy.zeros(nbBin2D)
		tmp_xi2DCovIdx  = numpy.zeros(nbBin2D)
		tmp_xi2DCovIdx2 = numpy.zeros(nbBin2D)
		idx = 0
		
		for k1 in range(0,nbBin2D__):
			i1       = k1/nbBinY2D__
			j1       = k1%nbBinY2D__
			k11      = j1*nbBinX2D__ + i1
			tmp_xi2D[k11] = xi2D_[i1][j1][1]
		
			for k2 in range(0,k1+1):
				i2       = k2/nbBinY2D__
				j2       = k2%nbBinY2D__
				k22      = j2*nbBinX2D__ + i2
				tmp_xi2DCov[idx]     = xi2DCov[k1][k2]
				tmp_xi2DCovIdx[idx]  = k11
				tmp_xi2DCovIdx2[idx] = k22
				idx += 1
		
		numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/DLA_systematiq/bao2D_'+path+'.data',zip(tmp_idx,tmp_xi2D),fmt='%d %E')
		numpy.savetxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Correlation_test2/DLA_systematiq/bao2D_'+path+'.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%d %d %E')
'''	
xi1D_ss, xi2D_ss, xiMu_ss, xiWe_ss = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/2560_QSO/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
xi1D_dd, xi2D_dd, xiMu_dd, xiWe_dd = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_QSO.txt')

coef = 1.
plt.errorbar(xi1D_ss[:,0], coef*(xi1D_ss[:,1]-xi1D_ss[-1,1]), yerr=coef*xi1D_ss[:,2], fmt='o',label='simulation',  markersize=8, elinewidth=6, color='blue')
plt.errorbar(xi1D_dd[:,0], coef*(xi1D_dd[:,1]-xi1D_dd[-1,1]), yerr=coef*xi1D_dd[:,2], fmt='o',label='data',  markersize=8, elinewidth=6, color='red')
plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
plt.show()

coef = xi1D_ss[:,0]
plt.errorbar(xi1D_ss[:,0], coef*(xi1D_ss[:,1]-xi1D_ss[-1,1]), yerr=coef*xi1D_ss[:,2], fmt='o',label='simulation',  markersize=8, elinewidth=6, color='blue' )
plt.errorbar(xi1D_dd[:,0], coef*(xi1D_dd[:,1]-xi1D_dd[-1,1]), yerr=coef*xi1D_dd[:,2], fmt='o',label='data',  markersize=8, elinewidth=6, color='red' )
plt.ylabel(r'|s|.$\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
plt.show()

coef = xi1D_ss[:,0]**2.
plt.errorbar(xi1D_ss[:,0], coef*(xi1D_ss[:,1]-xi1D_ss[-1,1]), yerr=coef*xi1D_ss[:,2], fmt='o',label='simulation',  markersize=8, elinewidth=6, color='blue')
plt.errorbar(xi1D_dd[:,0], coef*(xi1D_dd[:,1]-xi1D_dd[-1,1]), yerr=coef*xi1D_dd[:,2], fmt='o',label='data',  markersize=8, elinewidth=6, color='red')
plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
plt.show()

plt.errorbar(xi1D_ss[:,0], (xi1D_ss[:,1])/(xi1D_dd[:,1]), fmt='o',label='data / simulation',  markersize=8, elinewidth=6, color='blue')
plt.ylabel(r'$\xi^{qf} \, \xi^{qf}$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
plt.show()



plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)




plt.plot(xi1D_[:,0],xi1D_[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],(xi1D2_[:,1]-xi1D2_[-1,1]), marker='o',label=r'$simu2$')
plt.plot(xi1D[:,0],xi1D[:,1], marker='o',label=r'$Data$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plt.plot(xi1D_[:,0],xi1D_[:,0]*xi1D_[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],xi1D2_[:,0]*(xi1D2_[:,1]-xi1D2_[-1,1]), marker='o',label=r'$simu2$')
plt.plot(xi1D[:,0],xi1D[:,0]*xi1D[:,1], marker='o',label=r'$Data$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plt.plot(xi1D_[:,0],xi1D_[:,0]*xi1D_[:,0]*xi1D_[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],xi1D2_[:,0]*xi1D2_[:,0]*(xi1D2_[:,1]-xi1D2_[-1,1]), marker='o',label=r'$simu2$')
plt.plot(xi1D[:,0],xi1D[:,0]*xi1D[:,0]*xi1D[:,1], marker='o',label=r'$Data$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plt.plot(xi1D_[:,0],xi1D_[:,1]/xi1D[:,1], marker='o',label=r'$simu$')
plt.plot(xi1D2_[:,0],(xi1D2_[:,1]-xi1D2_[-1,1])/xi1D[:,1], marker='o',label=r'$simu2$')
myTools.deal_with_plot(False,False,True)
plt.ylabel(r'$\xi^{qf} (|s|) \, / \, \xi^{qf}_{data} (|s|)$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.show()

plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)



xi1Dboot_, xi2Dboot_   = loadBoot(path1__ +'xi_delta_QSO_1D_LYA_QSO_randomQSO_', path1__ +'xi_delta_QSO_2D_LYA_QSO_randomQSO_',2000)
numpy.save(path1__ + 'randomQSO_1D',xi1Dboot_)
numpy.save(path1__ + 'randomQSO_2D',xi2Dboot_)
### For Rand
xi1DList = numpy.load(path1__ + 'randomQSO_1D.npy')
xi2DList = numpy.load(path1__ + 'randomQSO_2D.npy')
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_QSO.txt', path1__ +'xi_delta_QSO_2D_LYA_QSO.txt')
xxx  = xi1D_[:,0]
coef = 1. #(xi1D_[:,0]**2.)
mean = numpy.mean(xi1DList,axis=1)
plt.errorbar(xxx, coef*xi1D_[:,1], yerr=coef*xi1D_[:,2], marker='o', color='blue', label=r'$Data$')
plt.errorbar(xxx, coef*mean, marker='o', color='red', label=r'$mean \, over \, random$')

myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$a$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

### For Rand
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_QSO.txt', path1__ +'xi_delta_QSO_2D_LYA_QSO.txt')
xxx  = numpy.arange(0,xi2D_[:,:,1].size)
coef = 1. #(xi1D_[:,0]**2.)
mean = numpy.mean(xi2DList,axis=1)
plt.errorbar(xxx, coef*xi2D_[:,:,1].flatten(), marker='o', color='blue', label=r'$Data$')
plt.errorbar(xxx, coef*mean, marker='o', color='red', label=r'$mean \, over \, random$')

myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$a$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

### plot 2D
array2D = convert1DTo2D(mean)
fig = plt.figure()
ax = fig.add_subplot(111)

plt.imshow(array2D, origin='lower', interpolation='None')
cbar = plt.colorbar()
plt.grid(True)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()
plt.show()
### plot 2D
array2D = xi2D_[:,:,1]
fig = plt.figure()
ax = fig.add_subplot(111)

plt.imshow(array2D, origin='lower', interpolation='None')
cbar = plt.colorbar()
plt.grid(True)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()
plt.show()


xi1DList    = [ xi1DList ]
xi2DList    = [ xi2DList ]
covNameList = ['hello']
covList1D   = [ numpy.cov(el) for el in xi1DList ]
covList2D   = [ numpy.cov(el) for el in xi2DList ]

covFromFit1D = plotCovar(covList1D, covNameList)
covFromFit2D = plotCovar(covList2D, covNameList)

#### Pilot
'''
for ii in range(0,10):
	for jj in range(0,10):

		xi1DList    = [ numpy.load(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_1D.npy') ]
		xi2DList    = [ numpy.load(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_2D.npy') ]
		covNameList = ['Mock DLA: \, '+str(ii)+' / '+str(jj)]

		size = xi1DList[0][0].size
		print ii, jj, xi1DList[0][0].size, xi2DList[0][0].size, size

		covList1D   = [ numpy.cov(el)/size for el in xi1DList ]
		covList2D   = [ numpy.cov(el)/size for el in xi2DList ]

		covFromFit1D = plotCovar(covList1D, covNameList)
		covFromFit2D = plotCovar(covList2D, covNameList)

		numpy.save(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_fromFit_1D.npy',covFromFit1D)
		numpy.save(path1__ + 'mockDLA_'+str(ii)+'_'+str(jj)+'_subsampling_fromFit_2D.npy',covFromFit2D)
	
'''

















































#a = loadWick()
#plotCovar(a, ['T1','T12','T123','T1234','T12345','T123456'])
#xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt', path1__ +'xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')

### with randomized forest #686
#xi1Dboot_, xi2Dboot_   = loadBoot(path1__ +'xi_delta_QSO_1D_'+forest1__+'_'+qso1__+'_shuffleForest_', path1__ +'xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'_randForest_',340)
### with mocks of delta #100
#xi1Dboot2_, xi2Dboot2_ = loadBoot(path1__ +'xi_delta_QSO_1D_'+forest1__+'_'+qso1__+'_mocks__', path1__ +'xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'_mocks__', 1, 10, 'mock')
### with subsampling #80
#xi1Dboot3_, xi2Dboot3_ = loadBoot(path1__ +'xi_delta_QSO_1D_'+forest1__+'_'+qso1__+'_bootstrap_', path1__ +'xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'_bootstrap_', 80)
### with randomized QSO #221
#xi1Dboot4_, xi2Dboot4_ = loadBoot(path1__ +'xi_delta_QSO_1D_'+forest1__+'_'+qso1__+'_randQSO_', path1__ +'xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'_randQSO_', 221)
### with random QSO #200
#xi1Dboot5_, xi2Dboot5_ = loadBoot(path1__ +'xi_delta_QSO_1D_'+forest1__+'_'+qso1__+'_randomQSO_', path1__ +'xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'_randomQSO_', 200)
### with DLA of mock
#xi1Dboot6_, xi2Dboot6_ = loadBoot(path1__ +'xi_delta_QSO_1D_LYA_DLA_mocks__', path1__ +'xi_delta_QSO_2D_LYA_DLA_mocks__', 10, 10, 'mock')

'''
#xi1Dboot3_ /= numpy.sqrt(xi1Dboot3_[0].size)
#xi2Dboot3_ /= numpy.sqrt(xi2Dboot3_[0].size)

numpy.save(path1__ + 'ShuffleForest_1D',xi1Dboot_)
numpy.save(path1__ + 'ShuffleForest_2D',xi2Dboot_)
numpy.save(path1__ + 'mocks_1D',xi1Dboot2_)
numpy.save(path1__ + 'mocks_2D',xi2Dboot2_)
numpy.save(path1__ + 'subsampling_1D',xi1Dboot3_)
numpy.save(path1__ + 'subsampling_2D',xi2Dboot3_)
numpy.save(path1__ + 'ShuffleQSO_1D',xi1Dboot4_)
numpy.save(path1__ + 'ShuffleQSO_2D',xi2Dboot4_)
numpy.save(path1__ + 'randomQSO_1D',xi1Dboot5_)
numpy.save(path1__ + 'randomQSO_2D',xi2Dboot5_)
numpy.save(path1__ + 'mockDLA_1D',xi1Dboot6_)
numpy.save(path1__ + 'mockDLA_2D',xi2Dboot6_)
'''
'''
xi1DList = [numpy.load(path1__ + 'ShuffleForest_1D.npy'),
		numpy.load(path1__ + 'mocks_1D.npy'),
		numpy.load(path1__ + 'subsampling_1D.npy'),
		numpy.load(path1__ + 'ShuffleQSO_1D.npy'),
		numpy.load(path1__ + 'randomQSO_1D.npy'),
		numpy.load(path1__ + 'mockDLA_1D.npy')  ]
xi2DList = [numpy.load(path1__ + 'ShuffleForest_2D.npy'),
		numpy.load(path1__ + 'mocks_2D.npy'),
		numpy.load(path1__ + 'subsampling_2D.npy'),
		numpy.load(path1__ + 'ShuffleQSO_2D.npy'),
		numpy.load(path1__ + 'randomQSO_2D.npy'),
		numpy.load(path1__ + 'mockDLA_2D.npy')  ]
covNameList = ['Shuffle \, forests', 'Mocks \, QSO','Subsampling','Shuffle \, QSO', 'Random \, QSO', 'Mocks \, DLA']
covList1D[2] /= xi1DList[2][0].size
covList2D[2] /= xi1DList[2][0].size
'''




'''
### with subsampling #80
for i in range(0,10):
	for j in range(0,10):
		xi1Dboot3_, xi2Dboot3_ = loadBoot(path1__ +'xi_delta_QSO_1D_LYA_DLA_bootstrap_', path1__ +'xi_delta_QSO_2D_LYA_DLA_bootstrap_', 80, 1, 'mock', i, j)
		numpy.save(path1__ + 'mockDLA_'+str(i)+'_'+str(j)+'_subsampling_1D',xi1Dboot3_)
		numpy.save(path1__ + 'mockDLA_'+str(i)+'_'+str(j)+'_subsampling_2D',xi2Dboot3_)
		print '  save:  mockDLA_'+str(i)+'_'+str(j)+'_subsampling_*D'
'''



### For Rand
xi1DList = numpy.load(path1__ + 'randomQSO_1D.npy')
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_QSO.txt', path1__ +'xi_delta_QSO_2D_LYA_QSO.txt')
xxx  = xi1D_[:,0]
coef = 1. #(xi1D_[:,0]**2.)
mean = numpy.mean(xi1DList,axis=1)
plt.errorbar(xxx, coef*mean, marker='o', color='red', label=r'$mean \, over \, random$')
plt.errorbar(xxx, coef*xi1D_[:,1], yerr=coef*xi1D_[:,2], marker='o', color='blue', label=r'$Data$')

myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$a$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

coef = 1. #(xi1D_[:,0]**2.)
mean = numpy.mean(xi1DList,axis=1)
plt.errorbar(xxx, coef*(xi1D_[:,1]-mean), yerr=coef*xi1D_[:,2], marker='o', color='blue', label=r'$Data$')

myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$a$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()




for j in range(0,len(xi1DList)):
	xxx = 4.*numpy.arange(0,xi1DList[j][:,0].size)+2.
	mean = numpy.mean(xi1DList[j],axis=1)

	#coef = 1.
	coef = (xxx**2.)

	for i in range(0,xi1DList[j][0].size):
		if ( i==0):
			plt.errorbar(xxx, coef*xi1DList[j][:,i], color='blue', label=r'$Realisation \, N \degree i$', alpha=0.01)
		else:
			plt.errorbar(xxx, coef*xi1DList[j][:,i], color='blue', alpha=0.4)
	plt.errorbar(xxx, coef*mean, yerr=coef*numpy.sqrt(numpy.var(xi1DList[j],axis=1)/xi1DList[j][0].size), marker='o', color='red', label=r'$mean$')

	### plot the overall x-correlation
	xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path1__ +'xi_delta_QSO_Mu_LYA_DLA_mocks__0_0_.txt', path1__ +'xi_delta_QSO_2D_LYA_DLA_mocks__0_0_.txt')
	coef = (xi1D_[:,0]**2.)
	plt.errorbar(xi1D_[:,0], coef*xi1D_[:,1], yerr=coef*xi1D_[:,2], marker='o', color='green', label=r'$all \, sky$')

	myTools.deal_with_plot(False,False,True)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	plt.ylabel(r'$|s|^{2}.\xi (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	plt.title(r'$'+covNameList[j]+'$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.show()



### show the diagonal
for i in range(0,len(xi1DList)):
	diag  = numpy.sqrt(numpy.diagonal(covList1D[i]))
	diag2 = numpy.sqrt(numpy.diagonal(covList1D[2]))
	plt.plot(xi1D_[:,0], diag/diag2, marker='o', label=r'$'+covNameList[i]+'$')
myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$err_{i}/err_{subsampling}$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()
### show the diagonal
for i in range(0,len(xi2DList)):
	diag  = numpy.sqrt(numpy.diagonal(covList2D[i]))
	diag2 = numpy.sqrt(numpy.diagonal(covList2D[2]))
	plt.plot(numpy.arange(0,diag.size), diag/diag2, marker='o', label=r'$'+covNameList[i]+'$')
myTools.deal_with_plot(False,False,True)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$err_{i}/err_{subsampling}$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)





#########################
'''
a = ['QSO_ALL_DR7','QSO_ALL_DR12', 'QSO_ALL_NGC', 'QSO_ALL_SGC', 'QSO_ALL_LowZ', 'QSO_ALL_highZ', 'QSO_ALL_CORE', 'QSO_ALL_LowMg', 'QSO_ALL_HighMg']
import colorsys
aa = []
for el in a:
	aa += [ loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_Mu_LYA_'+el+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/xi_delta_QSO_2D_LYA_'+el+'.txt') ]
for i in range(0,3):
	for j in range(0,len(a)):
		xi1D_, xi2D_, xiMu_, xiWe_ = aa[j]
		coef = numpy.power(xi1D_[:,0],i)
		plt.errorbar(xi1D_[:,0], coef*(xi1D_[:,1]-xi1D_[-1,1]), yerr=coef*xi1D_[:,2], markersize=5, marker='o', label=a[j])

	plt.ylabel(r'$\\xi^{qf} (|s|)$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([min1D__,max1D__])
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2, loc=4)
	plt.show()

nbBoot = 80
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
def send(something,idx):
	aa = a[idx]
	print aa
	shuffleForest1D = numpy.zeros( shape=(50,nbBoot) )
	shuffleForest2D = numpy.zeros( shape=(5000,nbBoot) )

	for j in range(0,nbBoot):

		xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path +'xi_delta_QSO_Mu_LYA_'+aa+'_bootstrap_'+str(j)+'.txt',path +'xi_delta_QSO_2D_LYA_'+aa+'_bootstrap_'+str(j)+'.txt')
	
		shuffleForest1D[:,j] = xi1D_[:,1]
		shuffleForest2D[:,j] = xi2D_[:,:,1].flatten()

	numpy.save(path +'subSampling_'+aa+'_1D',shuffleForest1D)
	numpy.save(path +'subSampling_'+aa+'_2D',shuffleForest2D)
	print '  save, ', aa

import multiprocessing

workers=[multiprocessing.Process(target=send,args=(True,idx)) for idx in numpy.arange(len(a))]
for p in workers:
	p.start()
for p in workers:
	p.join() 

'''
nbBoot = 80
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy/'
#loadBootMap()





















###################### For xi_qso_qso
plt.errorbar(data_xxx,data_yyy, label='dd', fmt='o')
plt.errorbar(rr_xxx,rr_yyy, label='rr', fmt='o')
plt.errorbar(dr_xxx,dr_yyy, label='dr', fmt='o')
myTools.deal_with_plot(False,False,True)
plt.show()
plt.errorbar(data_xxx,data_yyy/(data_xxx**2.), label='dd', fmt='o')
plt.errorbar(rr_xxx,rr_yyy/(rr_xxx**2.), label='rr', fmt='o')
plt.errorbar(dr_xxx,dr_yyy/(dr_xxx**2.), label='dr', fmt='o')
myTools.deal_with_plot(False,False,True)
plt.show()

plt.errorbar(data_xxx,data_xxx-data_xxx, label='dd', fmt='o')
plt.errorbar(data_xxx,rr_xxx-data_xxx, label='rr', fmt='o')
plt.errorbar(data_xxx,dr_xxx-data_xxx, label='dr', fmt='o')
myTools.deal_with_plot(False,False,True)
plt.ylim([-1.,1.])
plt.show()
plt.errorbar(data_xxx,100.*(data_xxx-data_xxx)/data_xxx, label='dd', fmt='o')
plt.errorbar(data_xxx,100.*(rr_xxx-data_xxx)/data_xxx, label='rr', fmt='o')
plt.errorbar(data_xxx,100.*(dr_xxx-data_xxx)/data_xxx, label='dr', fmt='o')
myTools.deal_with_plot(False,False,True)
plt.ylim([-100.,100.])
plt.show()



xxx = numpy.zeros(nbBin1D__)
yyy = numpy.zeros(nbBin1D__)
yer = numpy.zeros(nbBin1D__)
cut = (data_yyy>0.)
xxx[cut] = data_xxx[cut]
yyy[cut] = (data_yyy[cut]-2.*dr_yyy[cut]+rr_yyy[cut])/rr_yyy[cut]
yer[cut] = numpy.sqrt(data_yyy[cut]*coefDD)/(coefDD*rr_yyy[cut])

xxxP = myTools.Rebin(xxx, [25])
yyyP = myTools.Rebin(yyy, [25])
yerP = myTools.Rebin(yer, [25])

#yyy  -= yyy[-1]
#yyyP -= yyyP[-1]

data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/eBOSS/xi_QSO_QSO_1D_QSO.txt')
xxx2 = (data[:,0]+data[:,1])/2.
yyy2 = data[:,2]-data[-1,2]
yer2 = data[:,3]
for rescale in range(0,3):
	if (rescale==0):
		#plt.errorbar(xxx, yyy, yerr=yer, fmt='o', label='Simulation',  markersize=8, elinewidth=6, color='blue')
		plt.errorbar(xxxP, yyyP, yerr=yerP, fmt='o', label='Simulation', color='blue',  markersize=8, elinewidth=6)
		plt.errorbar(xxx2, yyy2, yerr=yer2, fmt='o', label='Data', color='red',  markersize=8, elinewidth=6)
		plt.ylabel(r'$\\xi^{qq} (|s|)$', fontsize=40)
	if (rescale==1):
		#plt.errorbar(xxx, xxx*yyy,yerr=xxx*yer, fmt='o', label='Simulation',  markersize=8, elinewidth=6, color='blue')
		plt.errorbar(xxxP, xxxP*yyyP,yerr=xxxP*yerP, fmt='o', label='Simulation', color='blue',  markersize=8, elinewidth=6)
		plt.errorbar(xxx2, xxx2*yyy2,yerr=xxx2*yer2, fmt='o', label='Data', color='red',  markersize=8, elinewidth=6)
		plt.ylabel(r'$|s|^{1}.\\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		#plt.errorbar(xxx, (xxx**2.)*yyy,yerr=(xxx**2.)*yer, fmt='o', label='Simulation',  markersize=8, elinewidth=6, color='blue')
		plt.errorbar(xxxP, (xxxP**2.)*yyyP,yerr=(xxxP**2.)*yerP, fmt='o', label='Simulation', color='blue',  markersize=8, elinewidth=6)
		plt.errorbar(xxx2, (xxx2**2.)*yyy2,yerr=(xxx2**2.)*yer2, fmt='o', label='Data', color='red',  markersize=8, elinewidth=6)
		plt.ylabel(r'$|s|^{2}.\\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ min1D__, max1D__ ])
	plt.show()

plt.errorbar(xxxP, yyyP/yyy2, yerr=yerP, fmt='o', label='SimulationP', color='blue')
plt.title(r'', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ min1D__, max1D__ ])
plt.show()







######################

pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,0)

result_Multipol = plotMultipol(xiMu_)
testMultipol = numpy.asarray(xi1D_)
testMultipol[:,1] = result_Multipol[:,0,0]
testMultipol[:,2] = result_Multipol[:,0,1]

pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(testMultipol,pathToCamb,0)

result_Multipol = plotMultipol(xiMu_)
testMultipol = numpy.asarray(xi1D_)
testMultipol[:,1] = result_Multipol[:,2,0]
testMultipol[:,2] = result_Multipol[:,2,1]

pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(testMultipol,pathToCamb,2)






##########################
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/xi_delta_QSO_1D_Wick0_LYA_QSO_T1'
var, cov1    = loadWick(path+'_10.txt')
var, cov12   = loadWick(path+'2_10.txt')
var, cov123  = loadWick(path+'23_10.txt')
var, cov1234 = loadWick(path+'234_10.txt')

cor1 = myTools.getCorrelationMatrix(cov1)
cor12 = myTools.getCorrelationMatrix(cov12)
cor123 = myTools.getCorrelationMatrix(cov123)
cor1234 = myTools.getCorrelationMatrix(cov1234)

myTools.plot2D(cor1)
myTools.plot2D(cor12)
myTools.plot2D(cor123)
myTools.plot2D(cor1234)
plotCovar([cov1,cov12,cov123,cov1234], ['1','12','123','1234'])

path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1521/noNoise_noCont_correctLines_moreQSO_withRSD_randomPositionInCell/'
xi1D_ss, xi2D_ss, xiMu_ss, xiWe_ss = loadData(path + 'xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt',path + 'xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1521/noNoise_noCont_correctLines_moreQSO_withRSD/'
xi1D_dd, xi2D_dd, xiMu_dd, xiWe_dd = loadData(path + 'xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt',path + 'xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')


for i in range(0,3):
	coef = numpy.power(xi1D_ss[:,0],i)
	plt.errorbar(xi1D_ss[:,0], coef*(xi1D_ss[:,1]), yerr=coef*xi1D_ss[:,2], fmt='o',label='Simu random')
	plt.errorbar(xi1D_dd[:,0], coef*(xi1D_dd[:,1]), yerr=coef*xi1D_dd[:,2], fmt='o',label='Simu')
	plt.errorbar(xi1D_[:,0], coef*(xi1D_[:,1]), yerr=coef*xi1D_[:,2], fmt='o',label='data')
	plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)

	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
	plt.show()







#####################################################
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/eBOSS/xi_QSO_QSO_1D_QSO.txt')
result_1D_eBOSS = numpy.zeros( shape=(data[:,2].size,3) )
result_1D_eBOSS[:,0] = (data[:,0]+data[:,1])/2.
result_1D_eBOSS[:,1] = data[:,2]
result_1D_eBOSS[:,2] = data[:,3]



rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CreateRandomField/QSO_45/'
nd = 104213
nr = 104213
nbRand = 10
result_1D,result_2D,result_Mu,result_We = getACorrelation(nd,nr,nbRand,rawPath)

plotMultipol(result_Mu)
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
b = fitCamb(result_1D,pathToCamb,0)














rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/noRSD/'
nbRand = 10
result_1D_no,result_2D_no,result_Mu_no,result_We_no = getACorrelation(nd,nr,nbRand,rawPath)


pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
b = fitCamb(result_1D_no,pathToCamb,0)

data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
cut = (data[1:,0] <= 200.)
size = cut[cut].size
result_1D_camb = numpy.zeros( shape=(size,3) )
result_1D_camb[:,0] = data[1:,0][cut]
result_1D_camb[:,1] = b*b*data[1:,1][cut]
result_1D_camb[:,2] = 0.0000000001



for rescale in range(0,3):
	for xi in [[result_1D,'withRSD'],[result_1D_no,'noRSD'],[result_1D_camb,'CAMB - z = 2.4, b = '+str(b)],[result_1D_eBOSS,'eBOSS']]:

		cut = (xi[0][:,2]!=0.)
		xxx = xi[0][:,0][cut]
		yyy = xi[0][:,1][cut]
		yer = xi[0][:,2][cut]

		coef = numpy.power(xxx,rescale)
		if (xi[1]!='CAMB - z = 2.4, b = '+str(b)):
			plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label=xi[1])
		else:
			plt.errorbar(xxx, coef*yyy, label=xi[1])

		if (rescale==0):
			plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
		if (rescale==1):
			plt.ylabel(r'$|s|^{1}.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (rescale==2):
			plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()


for rescale in range(0,3):
	cut = (result_1D_no[:,2]!=0.)
	xxx = result_1D[:,0][cut]
	yyy = result_1D[:,1][cut]/result_1D_no[:,1][cut]
	yer = result_1D[:,2][cut]

	cut = (xxx<=50.)
	print numpy.mean(yyy[cut])

	coef = numpy.power(xxx,rescale)
	plt.errorbar(xxx, coef*yyy, yerr=coef*yer, fmt='o', label='withRSD / noRSD')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|^{1}.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()


### 1D:
plot_Xi_1D(result_1D, 0)
plot_Xi_1D(result_1D, 1)
plot_Xi_1D(result_1D, 2)

### 2D:
plotXi2D(result_2D, 0)
plotXi2D(result_2D, 1)
plotXi2D(result_2D, 2)

### xiMu
plotMu(result_Mu,0)
plotMu(result_Mu,1)
plotMu(result_Mu,2)

### xiWe
plotWe(result_We,0)
plotWe(result_We,1)
plotWe(result_We,2)

### Multipols
plotMultipol(result_Mu)

### 1D:
plot_Xi_1D(result_1D_no, 0)
plot_Xi_1D(result_1D_no, 1)
plot_Xi_1D(result_1D_no, 2)

### 2D:
plotXi2D(result_2D_no, 0)
plotXi2D(result_2D_no, 1)
plotXi2D(result_2D_no, 2)

### xiMu
plotMu(result_Mu_no,0)
plotMu(result_Mu_no,1)
plotMu(result_Mu_no,2)

### xiWe
plotWe(result_We_no,0)
plotWe(result_We_no,1)
plotWe(result_We_no,2)

### Multipols
plotMultipol(result_Mu_no)



##############################################################################














path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Results/'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path+'xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt',path+'xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')


'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/withRSD_RandomPositionInBoxQSO/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/withRSD_RandomPositionInBoxQSO/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
'''
plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)
plotMultipol(xiMu_)
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
b = fitCamb(xi1D_,pathToCamb,0)




'''
nbBoot = 100
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/'
#loadBootMap()


shuffleForest1D = numpy.zeros( shape=(50,nbBoot) )
shuffleForest2D = numpy.zeros( shape=(5000,nbBoot) )

for i in range(0,nbBoot):

	xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path +'xi_delta_QSO_Mu_LYA_QSO_shuffleForest_'+str(i)+'.txt',path +'xi_delta_QSO_2D_LYA_QSO_shuffleForest_'+str(i)+'.txt')
	
	shuffleForest1D[:,i] = xi1D_[:,1]
	shuffleForest2D[:,i] = xi2D_[:,:,1].flatten()

numpy.save(path +'shuffleForest_LYA_QSO_1D',shuffleForest1D)
numpy.save(path +'shuffleForest_LYA_QSO_2D',shuffleForest2D)


a = ['LYA_QSO']

for aaa in a:
	xi1Dboot_ = numpy.load(path +'shuffleForest_'+aaa+'_1D.npy')
	xi2Dboot_ = numpy.load(path +'shuffleForest_'+aaa+'_2D.npy')
	print xi2Dboot_[0].size
	
	xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path +'xi_delta_QSO_Mu_LYA_QSO.txt',path +'xi_delta_QSO_2D_LYA_QSO.txt')
	mean = numpy.mean(xi1Dboot_,axis=1)
	cov1D = numpy.cov(xi1Dboot_)/nbBoot
	cov2D = numpy.cov(xi2Dboot_)/nbBoot
	
	###
	for i in range(0,nbBoot):
		plt.errorbar(numpy.arange(0,mean.size), xi1Dboot_[:,i], marker='o')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	
	###
	plt.errorbar(numpy.arange(0,mean.size), xi1D_[:,1], yerr=xi1D_[:,2], marker='o', label=r'$Simul$')
	plt.errorbar(numpy.arange(0,mean.size), mean, yerr=numpy.std(xi1Dboot_,axis=1)/numpy.sqrt(xi1Dboot_[0].size), marker='o', label=r'$schuffle$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	mean = numpy.mean(xi2Dboot_,axis=1)
	plt.errorbar(numpy.arange(0,mean.size), xi2D_[:,:,1].flatten(), yerr=xi2D_[:,:,2].flatten(), marker='o', label=r'$Simul$')
	plt.errorbar(numpy.arange(0,mean.size), mean, yerr=numpy.std(xi2Dboot_,axis=1)/numpy.sqrt(xi2Dboot_[0].size), marker='o', label=r'$schuffle$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	###
	mean = numpy.diagonal(cov1D)
	plt.errorbar(numpy.arange(0,mean.size), mean, marker='o', label=r'$schuffle$')
	plt.errorbar(numpy.arange(0,mean.size), xi1D_[:,2]**2., marker='o', label=r'$Simul$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	mean = numpy.diagonal(cov2D)
	plt.errorbar(numpy.arange(0,mean.size), mean, marker='o', label=r'$schuffle$')
	plt.errorbar(numpy.arange(0,mean.size), xi2D_[:,:,2].flatten()**2., marker='o', label=r'$Simul$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	###
	mean = numpy.diagonal(cov1D)
	plt.errorbar(numpy.arange(0,mean.size), mean/(xi1D_[:,2]**2.), marker='o', label=r'$schuffle$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	mean = numpy.diagonal(cov2D)
	plt.errorbar(numpy.arange(0,mean.size), mean/(xi2D_[:,:,2].flatten()**2.), marker='o', label=r'$schuffle$')
	myTools.deal_with_plot(False,False,True)
	plt.show()
	

	cov1D = plotCovar([cov1D], ['Schuffle forest'])[0]

	cov2D = plotCovar([cov2D], ['Schuffle forest'])[0]
	
	cov = cov1D
	tmp_xi2D        = xi1D_[:,1]
	tmp_idx         = numpy.arange(0,nbBin1D__)
	nbBin1D = nbBin1D__*(nbBin1D__+1.)/2.
	tmp_xi2DCov     = numpy.zeros(nbBin1D)
	tmp_xi2DCovIdx  = numpy.zeros(nbBin1D)
	tmp_xi2DCovIdx2 = numpy.zeros(nbBin1D)
	idx = 0
	
	
	for k1 in range(0,nbBin1D__):
		for k2 in range(0,k1+1):
			#if (k1==k2): tmp_xi2DCov[idx]     = cov[k1][k2]
			#if (k1==k2): tmp_xi2DCov[idx]     = xi1D_[k1,2]**2.
			tmp_xi2DCov[idx]     = cov[k1][k2]
			tmp_xi2DCovIdx[idx]  = k1
			tmp_xi2DCovIdx2[idx] = k2
			idx += 1
	
	numpy.savetxt(path +'BAOFIT_model_'+aaa+'/bao1D.data',zip(tmp_idx,tmp_xi2D),fmt='%u %1.20e')
	numpy.savetxt(path +'BAOFIT_model_'+aaa+'/bao1D.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%u %u %1.20e')
	
	data = numpy.loadtxt(path +'xi_delta_QSO_1D_LYA_JMC_QSO_MockJmc.txt')

	save2 = data[:,2]
	save3 = data[:,4]
	save4 = data[:,5]
				
	meanRperp  = numpy.zeros(shape=(nbBin1D__,2))
	meanRedshift = numpy.zeros(2)
				
	for i in range(0,len(save2)):
		iX = i/binSize__
	
		meanRperp[iX][0]  += save2[i]
		meanRperp[iX][1]  += save4[i]
	meanRedshift[0] = numpy.sum(save3)
	meanRedshift[1] = numpy.sum(save4)
				
	meanRperp[:,0]  /= meanRperp[:,1]
				
	### Get the s_perp bin center
	stringRperp = ''
	for el in meanRperp[:-1,0]:
		stringRperp += str(el) + ','
	stringRperp += str(meanRperp[-1,0])
	print stringRperp
	
	### Get the redshift center
	stringRedshift = str(meanRedshift[0]/meanRedshift[1])
	print ' <z> = ', stringRedshift
	
	
	
	
	####
	createIni2(path +'xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt',path +'BAOFIT_model_'+aaa+'/bao2D.ini','BAOFIT_'+aaa)
	cov = cov2D
	tmp_xi2D        = numpy.zeros(nbBin2D__)
	tmp_idx         = numpy.arange(0,nbBin2D__)
	nbBin2D = nbBin2D__*(nbBin2D__+1.)/2.
	tmp_xi2DCov     = numpy.zeros(nbBin2D)
	tmp_xi2DCovIdx  = numpy.zeros(nbBin2D)
	tmp_xi2DCovIdx2 = numpy.zeros(nbBin2D)
	tmp2_xi2DCov     = numpy.zeros(nbBin2D)
	tmp2_xi2DCovIdx  = numpy.zeros(nbBin2D)
	tmp2_xi2DCovIdx2 = numpy.zeros(nbBin2D)
	idx = 0
	
	
	for k1 in range(0,nbBin2D__):
		i1       = k1/nbBinY2D__
		j1       = k1%nbBinY2D__
		k11      = j1*nbBinX2D__ + i1
		tmp_xi2D[k11] = xi2D_[i1][j1][1]
			
		for k2 in range(0,k1+1):
			i2       = k2/nbBinY2D__
			j2       = k2%nbBinY2D__
			k22      = j2*nbBinX2D__ + i2
			#if (k1==k2): tmp_xi2DCov[idx]     = xi2D_[i1,j1,2]**2.
			#if (k1==k2): tmp_xi2DCov[idx]     = cov[k1][k2]
			tmp_xi2DCov[idx]     = cov[k1][k2]
			tmp_xi2DCovIdx[idx]  = k11
			tmp_xi2DCovIdx2[idx] = k22
			idx += 1
			
	numpy.savetxt(path +'BAOFIT_model_'+aaa+'/bao2D.data',zip(tmp_idx,tmp_xi2D),fmt='%u %1.20e')
	numpy.savetxt(path +'BAOFIT_model_'+aaa+'/bao2D.cov',zip(tmp_xi2DCovIdx,tmp_xi2DCovIdx2,tmp_xi2DCov),fmt='%u %u %1.20e')
	
'''
	




'''
plotMultipol(xiMu_)
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
b = fitCamb(xi1D_,pathToCamb,0)
plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)

'''

#path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/noRSD_delta/'
#xi1D_ss, xi2D_ss, xiMu_ss, xiWe_ss = loadData(path + 'xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt',path + 'xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/withRSD_delta/'
xi1D_dd, xi2D_dd, xiMu_dd, xiWe_dd = loadData(path + 'xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt',path + 'xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
dd = plotMultipol(xiMu_dd)
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/noRSD_delta/'
xi1D_ss, xi2D_ss, xiMu_ss, xiWe_ss = loadData(path + 'xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt',path + 'xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
ss = plotMultipol(xiMu_ss)
xi1D_rr, xi2D_rr, xiMu_rr, xiWe_rr = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/withRSD_RandomPositionInBoxQSO/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/withRSD_RandomPositionInBoxQSO/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
rr = plotMultipol(xiMu_rr)


path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/'
xi1D_aa, xi2D_aa, xiMu_aa, xiWe_aa = loadData(path + 'xi_delta_QSO_Mu_LYA_QSO.txt',path + 'xi_delta_QSO_2D_LYA_QSO.txt')
aa = plotMultipol(xiMu_aa)





### Get the results
xxx = xi1D_dd[:,0]
color = ['blue','red','red']

### Show the result
for i in range(0,3):

	for j in range(0,5):
		result_xi = dd
		if ( result_xi[:,j,0][ (result_xi[:,j,0]!=0.) ].size == 0): continue

		tmp_xxx = xxx[ (result_xi[:,j,0]!=0.) ]
		tmp_yyy = result_xi[:,j,0][ (result_xi[:,j,0]!=0.) ]
		tmp_yer = result_xi[:,j,1][ (result_xi[:,j,0]!=0.) ]
		coef    = numpy.power(tmp_xxx,i)
		plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$Simu: \xi_{'+str(j)+'}$',color='red')

		result_xi = aa
		tmp_xxx = xxx[ (result_xi[:,j,0]!=0.) ]
		tmp_yyy = result_xi[:,j,0][ (result_xi[:,j,0]!=0.) ]
		tmp_yer = result_xi[:,j,1][ (result_xi[:,j,0]!=0.) ]
		coef    = numpy.power(tmp_xxx,i)
		plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$Data: \xi_{'+str(j)+'}$',color='blue')

		result_xi = ss
		tmp_xxx = xxx[ (result_xi[:,j,0]!=0.) ]
		tmp_yyy = result_xi[:,j,0][ (result_xi[:,j,0]!=0.) ]
		tmp_yer = result_xi[:,j,1][ (result_xi[:,j,0]!=0.) ]
		coef    = numpy.power(tmp_xxx,i)
		plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$Data: \xi_{'+str(j)+'}$',color='green')

		result_xi = rr
		tmp_xxx = xxx[ (result_xi[:,j,0]!=0.) ]
		tmp_yyy = result_xi[:,j,0][ (result_xi[:,j,0]!=0.) ]
		tmp_yer = result_xi[:,j,1][ (result_xi[:,j,0]!=0.) ]
		coef    = numpy.power(tmp_xxx,i)
		plt.errorbar(tmp_xxx, coef*tmp_yyy,  yerr=coef*tmp_yer,  fmt='o', label=r'$Simu rand position: \xi_{'+str(j)+'}$',color='black')

	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	#plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(tmp_xxx)-10., numpy.max(tmp_xxx)+10. ])
	plt.show()


for i in range(0,3):
	coef = numpy.power(xi1D_aa[:,0],i)
	plt.errorbar(xi1D_aa[:,0], coef*(xi1D_aa[:,1]-xi1D_aa[-1,1]), yerr=coef*xi1D_aa[:,2], fmt='o',label='Data',color='blue')
	coef = numpy.power(xi1D_dd[:,0],i)
	plt.errorbar(xi1D_dd[:,0], coef*(xi1D_dd[:,1]-xi1D_dd[-1,1]), yerr=coef*xi1D_dd[:,2], fmt='o',label='Simu',color='red')
	plt.errorbar(xi1D_ss[:,0], coef*(xi1D_ss[:,1]-xi1D_ss[-1,1]), yerr=coef*xi1D_ss[:,2], fmt='o',label='Simu noRSD',color='green')
	plt.errorbar(xi1D_rr[:,0], coef*(xi1D_rr[:,1]-xi1D_rr[-1,1]), yerr=coef*xi1D_rr[:,2], fmt='o',label='Simu rand position',color='black')
	plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)

	#plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xi1D_aa[:,0])-10., numpy.max(xi1D_aa[:,0])+10. ])
	plt.show()






###########################################################

correlations = []
a = numpy.array( ['4_5','4_5_negative'] )
b = numpy.ones(a.size).astype(int)*10
c = numpy.ones(a.size).astype(int)*200000
d = numpy.array( ['4_5','4_5_negative'] )

for i in numpy.arange(a.size):
	rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CreateRandomField/QSO_'+a[i]+'/'
	nd = c[i]
	nr = c[i]
	nbRand = b[i]
	result_1D,result_2D,result_Mu,result_We = getACorrelation(nd,nr,nbRand,rawPath)
	correlations += [ [result_1D,d[i]] ]


data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
cut = (data[1:,0] <= 200.)
size = cut[cut].size
result_1D_camb = numpy.zeros( shape=(size,3) )
result_1D_camb[:,0] = data[1:,0][cut]
result_1D_camb[:,1] = data[1:,1][cut]
result_1D_camb[:,2] = 0.0000000001
cambInterpolate = interpolate.interp1d(result_1D_camb[:,0],result_1D_camb[:,1],bounds_error=False,fill_value=0)
correlations += [ [result_1D_camb,'CAMB'] ]






for rescale in range(0,3):
	for xi in correlations:

		a = numpy.abs( cambInterpolate(xi[0][11,0])/xi[0][11,1] )
		print '  b = ', numpy.sqrt(1./a), a

		cut = (xi[0][:,2]!=0.)
		xxx = xi[0][:,0][cut]
		yyy = xi[0][:,1][cut]
		yer = xi[0][:,2][cut]

		coef = numpy.power(xxx,rescale)
		if ('CAMB'==xi[1]):
			plt.errorbar(xxx, coef*yyy, label=xi[1])
		else:
			yyy *= a
			yer *= a
			plt.errorbar(xxx, coef*yyy, label=xi[1],fmt='o')

		if (rescale==0):
			plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
		if (rescale==1):
			plt.ylabel(r'$|s|^{1}.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (rescale==2):
			plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ -10., 210. ])
	plt.show()







#############################################################################################

'''
i = sys.argv[5]
j = sys.argv[6]

path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_withRSD.fits'
cat = pyfits.open(path)[1].data
nd = cat.size
nr = cat.size
del cat
rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
nbRand = 10
result_1D,result_2D,result_Mu,result_We = getACorrelation(nd,nr,nbRand,rawPath)

pathToSave = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_QSO_QSO_result_'
numpy.save(pathToSave+'1D',result_1D)
numpy.save(pathToSave+'2D',result_2D)
numpy.save(pathToSave+'Mu',result_Mu)
numpy.save(pathToSave+'We',result_We)
'''
i= 0
j=0
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_withRSD.fits'
cat = pyfits.open(path)[1].data
nd = cat.size
nr = cat.size
del cat
rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
nbRand = 10
result_1D,result_2D,result_Mu,result_We = getACorrelation(nd,nr,nbRand,rawPath)

#saveListRealMocks(10,10)
a = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/xi_QSO_QSO_result_1D.npy')
mean = numpy.mean( a, axis=1 )

print a[0,:].size


for j in range(0,3):
	for i in range(0,a[0,:].size):
		xxx = numpy.arange(mean.size)*4.+2.	
		coef = numpy.power(xxx,j)
		plt.errorbar(xxx,coef*a[:,i],fmt='o')
	plt.show()

result_1D_mean = numpy.array(result_1D)
result_1D_mean[:,0] = numpy.arange(mean.size)*4.+2.
result_1D_mean[:,1] = numpy.mean( a, axis=1 )
result_1D_mean[:,2] = numpy.sqrt( numpy.diag( numpy.cov(a)/100. ) )
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
b = fitCamb(result_1D_mean,pathToCamb,0)


for i in range(0,3):
	xxx = numpy.arange(mean.size)*4.+2.
	coef = numpy.power(xxx,i)
	plt.errorbar(xxx,coef*mean,fmt='o')
	plt.show()

cov1D = numpy.cov(a)
myTools.plot2D(cov1D)

nbBin = 50
nbCov = 1
covList = [cov1D]
corList      = [numpy.zeros(shape=(nbBin,nbBin))]*nbCov
diagList     = [numpy.diagonal(covList[i]) for i in range(0,nbCov)]
diagSqrtList = [numpy.sqrt(diagList[i]) for i in range(0,nbCov)]

for i in range(0,nbCov):
	for j in range(0,nbBin):
		for k in range(0,nbBin):
			corList[i][j][k] = covList[i][j][k]/(diagSqrtList[i][j]*diagSqrtList[i][k])
	
	corList[i] = numpy.nan_to_num(corList[i])
myTools.plot2D(corList[0])


'''
### 1D:
plot_Xi_1D(result_1D, 0)
plot_Xi_1D(result_1D, 1)
plot_Xi_1D(result_1D, 2)

### 2D:
plotXi2D(result_2D, 0)
plotXi2D(result_2D, 1)
plotXi2D(result_2D, 2)

### xiMu
plotMu(result_Mu,0)
plotMu(result_Mu,1)
plotMu(result_Mu,2)

### xiWe
plotWe(result_We,0)
plotWe(result_We,1)
plotWe(result_We,2)

### Multipols
plotMultipol(result_Mu)
'''






#############################################################################################




correlations = []
a = numpy.array(	['/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Results/randomPositionOfQSOInCell/'] )
b = numpy.ones(a.size).astype(int)*10
c = numpy.ones(a.size).astype(int)*240597
d = numpy.array( ['3.5'] )

for i in numpy.arange(a.size):
	rawPath = a[i]
	nd = c[i]
	nr = c[i]
	nbRand = b[i]
	result_1D,result_2D,result_Mu,result_We = getACorrelation(nd,nr,nbRand,rawPath)

	
	plot_Xi_1D(result_1D, 0)
	plot_Xi_1D(result_1D, 1)
	plot_Xi_1D(result_1D, 2)

	### 2D:
	plotXi2D(result_2D, 0)
	plotXi2D(result_2D, 1)
	plotXi2D(result_2D, 2)

	### xiMu
	plotMu(result_Mu,0)
	plotMu(result_Mu,1)
	plotMu(result_Mu,2)

	### xiWe
	plotWe(result_We,0)
	plotWe(result_We,1)
	plotWe(result_We,2)

	fitCamb(result_1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
	

	correlations += [ [result_1D,d[i]] ]
	

data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
cut = (data[1:,0] <= 200.)
size = cut[cut].size
result_1D_camb = numpy.zeros( shape=(size,3) )
result_1D_camb[:,0] = data[1:,0][cut]
result_1D_camb[:,1] = data[1:,1][cut]
result_1D_camb[:,2] = 0.0000000001
cambInterpolate = interpolate.interp1d(result_1D_camb[:,0],result_1D_camb[:,1],bounds_error=False,fill_value=0)
correlations += [ [result_1D_camb,'CAMB'] ]






for rescale in range(0,3):
	for xi in correlations:

		idx = 11
		a = numpy.abs( cambInterpolate(xi[0][idx,0])/xi[0][idx,1] )
		print '  b = ', numpy.sqrt(1./a), a

		cut = (xi[0][:,2]!=0.)
		xxx = xi[0][:,0][cut]
		yyy = xi[0][:,1][cut]
		yer = xi[0][:,2][cut]

		coef = numpy.power(xxx,rescale)
		if ('CAMB'==xi[1]):
			plt.errorbar(xxx, coef*yyy, label=xi[1])
		else:
			yyy *= a
			yer *= a
			plt.errorbar(xxx, coef*yyy, label=xi[1],fmt='o')

		if (rescale==0):
			plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
		if (rescale==1):
			plt.ylabel(r'$|s|^{1}.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (rescale==2):
			plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
		
	plt.title(r'', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ -10., 210. ])
	plt.show()



#############################################################################################

#xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit5/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit5/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')
xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_delta_QSO_2D_LYA_QSO.txt')
plotXi(0)
plotXi(1)
plotXi(2)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)
plotMultipol(xiMu_)
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
b = fitCamb(xi1D_,pathToCamb,0)

'''
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_delta_QSO_1D_LYA_QSO.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
saveYYY = yyy
plt.errorbar(xxx,yyy/saveYYY,fmt='o',color='blue',label='Data')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Results/xi_delta_QSO_1D_LYA_JMC_QSO_MockJmc.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy/saveYYY,fmt='o',color='red',label='Simulation')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Results/xi_delta_QSO_1D_LYA_JMC_QSO_MockJmc.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy/saveYYY,fmt='o',color='cyan',label='Simulation corrected')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/Results/xi_delta_QSO_1D_LYA_QSO_mocks_.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy/saveYYY,fmt='o',color='green',label='Pipeline')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit2/xi_delta_QSO_1D_LYA_QSO.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy/saveYYY,fmt='o',color='black',label='Data new')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_1D_LYA_QSO.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy/saveYYY,fmt='o',color='black',label='Data new new')

plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$Var*n_{pairs} \, LYA-QSO$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()
'''




'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/'
name = path
saveListReal(80,name+'xi_delta_QSO_Mu_LYA_QSO_bootstrap_',name+'xi_delta_QSO_2D_LYA_QSO_bootstrap_',name+'subSampling_LYA_QSO_',True)
'''



xi1D, xi2D, xiMu, xiWe = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Results/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Results/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')


cov1D  = []
name   = []
cov1D += [ numpy.cov(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/shuffleForest_LYA_QSO_2D.npy')) ]
name  += ['Data \, shuffle \, Forest']
cov1D += [ numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/subSampling_LYA_QSO_cov_2D.npy') ]
name  += ['Data \, subsampling']

cov1D += [ numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/subSampling_LYA_QSO_cov_2D.npy') ]
name  += ['Data \, subsampling \, corrected' ]

cov1D += [ numpy.cov(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/mocksJMLG_LYA_QSO_2D.npy')) ]
name  += ['Mocks']

tmpcov1D = numpy.zeros( shape=(5000,5000) )
for i in range(0,10):
	for j in range(0,10):
		tmpcov1D += numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/subSampling_LYA_QSO_cov_2D.npy')

tmpcov1D /= 100
cov1D += [tmpcov1D]
name  += ['< Mock \, subsampling >']

for i in numpy.arange(len(cov1D)):
	coef = 1. #numpy.power(xi1D[:,0],2)
	print numpy.diag(cov1D[i]).size
	plt.errorbar(numpy.arange(5000), coef*numpy.diag(cov1D[i]),fmt='o',label=r'$'+name[i]+'$')
plt.ylabel(r'$|s|^{2}.Var(|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
#plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.show()

plotCovar(cov1D,name)







xi1Ddd, xi2Ddd, xiMudd, xiWedd = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_delta_QSO_Mu_LYA_QSO.txt')
xi1D, xi2D, xiMu, xiWe = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Results/xi_delta_QSO_Mu_LYA_JMC_QSO_MockJmc.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Results/xi_delta_QSO_2D_LYA_JMC_QSO_MockJmc.txt')
xi1DnoRand, xi2DnoRand, xiMunoRand, xiWenoRand = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit2/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit2/xi_delta_QSO_2D_LYA_QSO.txt')
xi1DnoRand2, xi2DnoRand2, xiMunoRand2, xiWenoRand2 = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_2D_LYA_QSO.txt')
print xi1DnoRand

'''
list1D = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/shuffleForest_LYA_QSO_1D.npy')
xi1Ddd[:,2] = numpy.sqrt( numpy.diag( numpy.cov(list1D)) )
'''


list1D = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/mocksJMLG_LYA_QSO_1D.npy')
xxx = xi1D[:,0]
yyy = xi1D[:,1]# numpy.mean(list1D,axis=1)
yer = xi1D[:,2]# numpy.sqrt( numpy.diag( numpy.cov(list1D)/100.) )


# yyy -= yyy[-1]
# xi1Ddd[:,1] -= xi1Ddd[-1,1]
# xi1DnoRand[:,1] -= xi1DnoRand[-1,1]

for rescale in [0,1,2]:
	if (rescale==0):
		plt.errorbar(xi1Ddd[:,0], xi1Ddd[:,1], yerr=xi1Ddd[:,2], fmt='o',color='blue',label='Data')
		plt.errorbar(xi1DnoRand[:,0], xi1DnoRand[:,1], yerr=xi1DnoRand[:,2], fmt='o',color='green',label='Data new')
		plt.errorbar(xi1DnoRand2[:,0], xi1DnoRand2[:,1], yerr=xi1DnoRand2[:,2], fmt='o',color='green',label='Data new')
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o',color='red',label='Simu')
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xi1Ddd[:,0], xi1Ddd[:,0]*xi1Ddd[:,1], yerr=xi1Ddd[:,0]*xi1Ddd[:,2], fmt='o',color='blue',label='Data')
		plt.errorbar(xi1DnoRand[:,0], xi1DnoRand[:,0]*xi1DnoRand[:,1], yerr=xi1DnoRand[:,0]*xi1DnoRand[:,2], fmt='o',color='green',label='Data new')
		plt.errorbar(xi1DnoRand2[:,0], xi1DnoRand2[:,0]*xi1DnoRand2[:,1], yerr=xi1DnoRand2[:,0]*xi1DnoRand2[:,2], fmt='o',color='green',label='Data new')
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, fmt='o',color='red',label='Simu')
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xi1Ddd[:,0], xi1Ddd[:,0]*xi1Ddd[:,0]*xi1Ddd[:,1], yerr=xi1Ddd[:,0]*xi1Ddd[:,0]*xi1Ddd[:,2], fmt='o',color='blue',label='Data')
		plt.errorbar(xi1DnoRand[:,0], xi1DnoRand[:,0]*xi1DnoRand[:,0]*xi1DnoRand[:,1], yerr=xi1DnoRand[:,0]*xi1DnoRand[:,0]*xi1DnoRand[:,2], fmt='o',color='green',label='Data new')
		plt.errorbar(xi1DnoRand2[:,0], xi1DnoRand2[:,0]*xi1DnoRand2[:,0]*xi1DnoRand2[:,1], yerr=xi1DnoRand2[:,0]*xi1DnoRand2[:,0]*xi1DnoRand2[:,2], fmt='o',color='green',label='Data new')
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, fmt='o',color='red',label='Simu')
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()

xi1Ddd[:,1] -= yyy

for rescale in [0,1,2]:
	if (rescale==0):
		plt.errorbar(xi1Ddd[:,0], xi1Ddd[:,1], yerr=xi1Ddd[:,2], fmt='o',color='blue',label='Data')
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xi1Ddd[:,0], xi1Ddd[:,0]*xi1Ddd[:,1], yerr=xi1Ddd[:,0]*xi1Ddd[:,2], fmt='o',color='blue',label='Data')
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xi1Ddd[:,0], xi1Ddd[:,0]*xi1Ddd[:,0]*xi1Ddd[:,1], yerr=xi1Ddd[:,0]*xi1Ddd[:,0]*xi1Ddd[:,2], fmt='o',color='blue',label='Data')
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()




list1D = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/mocksJMLG_LYA_QSO_2D.npy')
xxx = xi2D[:,:,0]
yyy = numpy.mean(list1D,axis=1)
yer = numpy.sqrt( numpy.diag( numpy.cov(list1D)) )
yyy = convert1DTo2D(yyy)
yer = convert1DTo2D(yer)

xxx = numpy.transpose(xxx)
yyy = numpy.transpose(yyy)
yer = numpy.transpose(yer)

for rescale in [0,1,2]:

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xticks([ 0.,50.,100.,150.,200.])
	ax.set_yticks([ -200.,-150.,-100.,-50.,0.,50.,100.,150.,200.])

	if (rescale==0):
		plt.imshow(yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__], interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$\xi^{qf}(\, \overrightarrow{s} \,)$',size=40)
	if (rescale==1):
		plt.imshow(xxx*yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|.\xi^{qf}(\, \overrightarrow{s} \,) \, [h^{-1}.Mpc]$',size=40)
	if (rescale==2):
		plt.imshow(xxx**2.*yyy, origin='lower',extent=[minX2D__, maxX2D__, minY2D__, maxY2D__],interpolation='None')
		cbar = plt.colorbar()
		cbar.set_label(r'$|s|^{2}.\xi^{qf}(\, \overrightarrow{s} \,) \, [(h^{-1}.Mpc)^{2}]$',size=40)

	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
	plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
	plt.grid(True)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.update_ticks()
	plt.show()


#############################################################################################
plot()

xi1D = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1536/Box_000/Simu_000/Results/xi_1DlRF_delta_delta_LYA_JMC.txt')
plt.errorbar(xi1D[:,0], xi1D[:,1], yerr=xi1D[:,2], fmt='o',label='simu')
xi1D = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_1DlRF_delta_delta_LYA.txt')
plt.errorbar(xi1D[:,0], xi1D[:,1], yerr=xi1D[:,2], fmt='o',label='data')
xi1D = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit/xi_1DlRF_delta_delta_LYA.txt')
plt.errorbar(xi1D[:,0], xi1D[:,1], yerr=xi1D[:,2], fmt='o',label='data 2')
xi1D = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit2/xi_1DlRF_delta_delta_LYA.txt')
plt.errorbar(xi1D[:,0], xi1D[:,1], yerr=xi1D[:,2], fmt='o',label='data 3')
xi1D = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_1DlRF_delta_delta_LYA.txt')
plt.errorbar(xi1D[:,0], xi1D[:,1], yerr=xi1D[:,2], fmt='o',label='data 4')
myTools.deal_with_plot(False,False,True)
plt.show()




data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_1DlRF_delta_delta_LYA.txt')
plt.errorbar(data[:,2]/data[:,4],data[:,5],fmt='o',color='blue',label='Data')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Results/xi_1DlRF_delta_delta_LYA_JMC.txt')
plt.errorbar(data[:,2]/data[:,4],data[:,5],fmt='o',color='red',label='Simulation')
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$nb \, pairs \, LYA-LYA \, same \, forest$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()

#############################################################################################

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_Mu_'+forest__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_A_delta_delta_2D_'+forest__+'.txt')
plot_Xi_1D(0)
plot_Xi_1D(1)
plot_Xi_1D(2)
plot_Xi_2D(0)
plot_Xi_2D(1)
plot_Xi_2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,0)
plotMultipol(xiMu_)




xi1D_dd, xi2D_dd, xiMu_dd, xiWe_dd = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_A_delta_delta_Mu_LYA.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_A_delta_delta_2D_LYA.txt')
xi1D_ss, xi2D_ss, xiMu_ss, xiWe_ss = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/withRSD_delta/xi_A_delta_delta_Mu_LYA_JMC.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1524/withRSD_delta/xi_A_delta_delta_2D_LYA_JMC.txt')
xi1D_mm, xi2D_mm, xiMu_mm, xiWe_mm = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/Results/xi_A_delta_delta_Mu_LYA.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/Results/xi_A_delta_delta_2D_LYA.txt')
print numpy.mean( xi1D_ss[:4,1]/xi1D_dd[:4,1] )
print numpy.sqrt( numpy.mean( xi1D_ss[:4,1]/xi1D_dd[:4,1]))



for i in range(0,3):
	coef = numpy.power(xi1D_ss[:,0],i)
	plt.errorbar(xi1D_dd[:,0], coef*xi1D_dd[:,1], yerr=coef*xi1D_dd[:,2], fmt='o',label='data',  markersize=8, elinewidth=6, color='blue')
	plt.errorbar(xi1D_ss[:,0], coef*xi1D_ss[:,1], yerr=coef*xi1D_ss[:,2], fmt='o',label='simulation',  markersize=8, elinewidth=6, color='red')
	plt.errorbar(xi1D_mm[:,0], coef*xi1D_mm[:,1], yerr=coef*xi1D_mm[:,2], fmt='o',label='pipeline',  markersize=8, elinewidth=6, color='green')

	if (i==0): plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
	if (i==1): plt.ylabel(r'$|s|.\xi^{ff} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (i==2): plt.ylabel(r'$|s|^{2}.\xi^{ff} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)

	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
	plt.show()


plt.errorbar(xi1D_ss[:,0], (xi1D_ss[:,1])/(xi1D_dd[:,1]), fmt='o',label='data / simulation',  markersize=8, elinewidth=6, color='blue')
plt.ylabel(r'$\xi^{ff} / $\xi^{ff}$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D_ss[:,0])-10., numpy.max(xi1D_ss[:,0])+10. ])
plt.show()




#############################################################
xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit/xi_A_delta_delta_Mu_LYA.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit/xi_A_delta_delta_2D_LYA.txt')
plot_Xi_1D(0)
plot_Xi_1D(1)
plot_Xi_1D(2)
plot_Xi_2D(0)
plot_Xi_2D(1)
plot_Xi_2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,0)
plotMultipol(xiMu_)



data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_A_delta_delta_1D_LYA.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy,fmt='o',color='blue',label='Data')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Results/xi_A_delta_delta_1D_LYA_JMC.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy,fmt='o',color='red',label='Simulation')
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/Results/xi_A_delta_delta_1D_LYA.txt')
xxx = data[:,2]/data[:,5]
yyy = (data[:,1]/data[:,5]-numpy.power(data[:,0]/data[:,5],2.))
print numpy.mean(yyy[ (xxx>50.) ])
plt.errorbar(xxx,yyy,fmt='o',color='green',label='Pipeline')

plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$Var*n_{pairs} \, LYA-LYA$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()









'''
nbBoot = 100
mocksJMLG = numpy.zeros( shape=(50,nbBoot) )

for i in range(0,10):
	for j in range(0,10):
		path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
		xi1D_, xi2D_, xiMu_, xiWe_ = loadData(path+'xi_A_delta_delta_Mu_LYA_JMC.txt',path+'xi_A_delta_delta_2D_LYA_JMC.txt')

		plt.errorbar(xi1D_[:,0], xi1D_[:,1], yerr=xi1D_[:,2], fmt='o')
		mocksJMLG[:,i*10+j] = xi1D_[:,1]

plt.ylabel(r'$\\xi^{ff} (|s|)$', fontsize=40)
plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()


numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/mocksJMLG_LYA_LYA_1D',mocksJMLG)
'''
mocksJMLG = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/mocksJMLG_LYA_LYA_1D.npy')

for i in range(0,3):
	mean = numpy.mean(mocksJMLG,axis=1)
	coef = numpy.power( numpy.arange(mean.size)*4.+2., i )
	plt.errorbar(numpy.arange(mean.size), coef*mean, fmt='o')
	plt.ylabel(r'$\xi^{ff} (|s|)$', fontsize=40)
	plt.title(r'$\delta_{'+forest__+'} \, - \, \delta_{'+forest__+'}$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,False)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()

cov1D = numpy.cov(mocksJMLG)
myTools.plot2D(cov1D)

corList      = numpy.array( cov1D )
diagList     = numpy.diagonal(cov1D)
diagSqrtList = numpy.sqrt(diagList)


for j in range(0,50):
	for k in range(0,50):
		corList[j][k] = cov1D[j][k]/(diagSqrtList[j]*diagSqrtList[k])
myTools.plot2D(corList)

pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,0)
plotMultipol(xiMu_)
plot_Xi_1D(0)
plot_Xi_1D(1)
plot_Xi_1D(2)
plot_Xi_2D(0)
plot_Xi_2D(1)
plot_Xi_2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)













#xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1521/noNoise_noCont_correctLines/xi_A_delta_delta_Mu_LYA_JMC.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1521/noNoise_noCont_correctLines/xi_A_delta_delta_2D_LYA_JMC.txt')

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1521/noNoise_noCont_correctLines/xi_A_delta_delta_Mu_LYA_JMC.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/v1521/noNoise_noCont_correctLines/xi_A_delta_delta_2D_LYA_JMC.txt')


result_Multipol = plotMultipol(xiMu_)
testMultipol = numpy.asarray(xi1D_)
testMultipol[:,1] = result_Multipol[:,0,0]
#testMultipol[:,2] = result_Multipol[:,0,1]

pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,0)

result_Multipol = plotMultipol(xiMu_)
testMultipol = numpy.asarray(xi1D_)
testMultipol[:,1] = result_Multipol[:,2,0]
#testMultipol[:,2] = result_Multipol[:,2,1]

pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,2)


plotMultipol(xiMu_)
plot_Xi_1D(0)
plot_Xi_1D(1)
plot_Xi_1D(2)
plot_Xi_2D(0)
plot_Xi_2D(1)
plot_Xi_2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)









'''
nbBoot = 100
mocksJMLG = numpy.zeros( shape=(161,nbBoot) )

for i in range(0,10):
	for j in range(0,10):
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_1DlRF_delta_delta_LYA_JMC.txt'
		xi1D = loadData(path)
		plt.errorbar(xi1D[:,0], xi1D[:,1], yerr=xi1D[:,2], fmt='o')

		mocksJMLG[:,i*10+j] = xi1D[:,1]

numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/listXi_LYA_LYA_1DlR',mocksJMLG)

plt.title(r'$1D: \, \delta_{'+forest__+'} \, - \, \delta_{'+forest__+'} $', fontsize=40)
plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$\\xi (s)$', fontsize=40)
myTools.deal_with_plot(False,False,False)
plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.show()
'''


mocksJMLG = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/listXi_LYA_LYA_1DlR.npy')
mocksJMLG = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/Results/listXi_LYA_LYA_1DlR.npy')

cov1D = numpy.cov(mocksJMLG)
myTools.plot2D(cov1D)

corList      = numpy.array( cov1D )
diagList     = numpy.diagonal(cov1D)
diagSqrtList = numpy.sqrt(diagList)


for j in range(0,161):
	for k in range(0,161):
		corList[j][k] = cov1D[j][k]/(diagSqrtList[j]*diagSqrtList[k])
myTools.plot2D(corList)





xi1DData = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_1DlRF_delta_delta_LYA.txt')
xi1D = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_000/Simu_000/Results/xi_1DlRF_delta_delta_LYA_JMC.txt') #numpy.mean(mocksJMLG,axis=1)
xi1DMocks = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_0/000/Results/xi_1DlRF_delta_delta_LYA.txt') 





xi1DMocks[:,1] = numpy.mean(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/Results/listXi_LYA_LYA_1DlR.npy'),axis=1)
xi1D[:,1] = numpy.mean(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Results/listXi_LYA_LYA_1DlR.npy'),axis=1)





print xi1D[0,1]/xi1DData[0,1]
print numpy.sqrt( xi1D[0,1]/xi1DData[0,1])

#xi1DData[:,1] /= xi1DData[0,1]
#xi1D[:,1] /= xi1D[0,1]
#xi1DMocks[:,1] /= xi1DMocks[0,1]
plt.errorbar(xi1DData[:,0], xi1DData[:,1], yerr=xi1DData[:,2], fmt='o',color='blue',label='Data')
plt.errorbar(xi1D[:,0], xi1D[:,1], yerr=xi1D[:,2], fmt='o',color='red',label='Simu')
plt.errorbar(xi1DMocks[:,0], xi1DMocks[:,1], yerr=xi1DMocks[:,2], fmt='o',color='green',label='Mocks Colab')
plt.title(r'$1D: \, \delta_{'+forest__+'} \, - \, \delta_{'+forest__+'} $', fontsize=40)
plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$\xi (s)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.show()

plt.errorbar(xi1D[:,0], xi1DData[:,1]-xi1D[:,1], fmt='o')
plt.title(r'$1D: \, \delta_{'+forest__+'} \, - \, \delta_{'+forest__+'} $', fontsize=40)
plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$\xi (s)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.show()
plt.errorbar(xi1D[:,0], (xi1DData[:,1]-xi1D[:,1])/xi1DData[:,1], fmt='o')
plt.title(r'$1D: \, \delta_{'+forest__+'} \, - \, \delta_{'+forest__+'} $', fontsize=40)
plt.xlabel(r'$s \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$\xi (s)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.show()







#####################################################################################################



xi1D, xi2D, xiMu, xiWe = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_2D_LYA_QSO.txt')


cov1D  = []
name   = []
cov1D += [ numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/subSampling_LYA_QSO_cov_1D.npy') ]
name  += ['Data \, subsampling']

cov1D += [ numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_delta_QSO_result_cov_1D.npy') ]
name  += ['Mocks']

'''
tmpcov1D = numpy.zeros( shape=(5000,5000) )
for i in range(0,10):
	for j in range(0,10):
		tmpcov1D += numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1528/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/subSampling_LYA_QSO_cov_2D.npy')

tmpcov1D /= 100
numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_delta_QSO_result_cov_2D_meanSubSampling.npy',tmpcov1D)
'''
tmpcov1D = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/xi_delta_QSO_result_cov_1D_meanSubSampling.npy')
cov1D += [tmpcov1D]
name  += ['< Mock \, subsampling >']


for i in numpy.arange(len(cov1D)):
	coef = numpy.power(xi1D[:,0],1)
	print numpy.diag(cov1D[i]).size
	plt.errorbar(xi1D[:,0], coef*numpy.diag(cov1D[i]),fmt='o',label=r'$'+name[i]+'$')
plt.ylabel(r'$|s|.Var(|s|) \, [h^{-1}.Mpc]$', fontsize=40)
plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.xlim([ numpy.min(xi1D[:,0])-10., numpy.max(xi1D[:,0])+10. ])
plt.show()

plotCovar(cov1D,name)





'''
i = sys.argv[6]
j = sys.argv[7]
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
name = path
saveListReal(80,name+'xi_delta_QSO_Mu_LYA_QSO_subsampling_',name+'xi_delta_QSO_2D_LYA_QSO_subsampling_',name+'subSampling_LYA_QSO_result_',True)
'''

'''
ni = 10
nj = 10
saveListRealMocks(ni,nj)
'''

xi1DD, xi2DD, xiMuDD, xiWeDD = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_2D_LYA_QSO.txt')
result_Multipol_DD = plotMultipol(xiMuDD)

path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
tmpPath =  path + '0/Simu_000/Results/xi_delta_QSO_'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(tmpPath+'Mu_LYA_QSO.txt',tmpPath+'2D_LYA_QSO.txt')
result_Multipol = plotMultipol(xiMu_)

rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/'

xi1D_[:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_1D.npy'),axis=1)
xi1D_[:,2] = numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_1D.npy')))/numpy.sqrt(100.)

print numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_2D.npy'),axis=1).size
xi2D_[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_2D.npy'),axis=1), nbBinX2D__,nbBinY2D__)
xi2D_[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_2D.npy')))/numpy.sqrt(100.), nbBinX2D__,nbBinY2D__)

xiMu_[:,:,2] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_Mu.npy'),axis=1), nbBin1D__,nbBinM__)
xiMu_[:,:,3] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_Mu.npy')))/numpy.sqrt(100.), nbBin1D__,nbBinM__)

result_Multipol[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_Multipol.npy'),axis=2)
result_Multipol[:,:,2] = numpy.var(numpy.load(rawPath+'xi_delta_QSO_result_Multipol.npy'),axis=2)/numpy.sqrt(100.)

xiWe_[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_We.npy'),axis=2)
xiWe_[:,:,2] = numpy.var(numpy.load(rawPath+'xi_delta_QSO_result_We.npy'),axis=2)/numpy.sqrt(100.)


'''
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')



plotXi2D(0)
plotXi2D(1)
plotXi2D(2)

plotMu(0)
plotMu(1)
plotMu(2)

plotWe(0)
plotWe(1)
plotWe(2)
'''
### Multipol
xi1D_[:,0] = result_Multipol[:,0,0]
xi1D_[:,1] = result_Multipol[:,0,1]
xi1D_[:,2] = result_Multipol[:,0,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
### Multipol
xi1D_[:,0] = result_Multipol[:,2,0]
xi1D_[:,1] = result_Multipol[:,2,1]
xi1D_[:,2] = result_Multipol[:,2,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat',2)
### Multipol

for rescale in range(0,3):
	###
	cut = (result_Multipol[:,0,2]!=0.)
	xxx0 = result_Multipol[:,0,0][cut]
	yyy0 = result_Multipol[:,0,1][cut]
	yer0 = result_Multipol[:,0,2][cut]
	###
	cut = (result_Multipol[:,2,2]!=0.)
	xxx1 = result_Multipol[:,2,0][cut]
	yyy1 = result_Multipol[:,2,1][cut]
	yer1 = result_Multipol[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Simu \, \xi_{0}$',color='red')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Simu \, \xi_{2}$',color='red')

	###
	cut = (result_Multipol_DD[:,0,2]!=0.)
	xxx0 = result_Multipol_DD[:,0,0][cut]
	yyy0 = result_Multipol_DD[:,0,1][cut]
	yer0 = result_Multipol_DD[:,0,2][cut]
	###
	cut = (result_Multipol_DD[:,2,2]!=0.)
	xxx1 = result_Multipol_DD[:,2,0][cut]
	yyy1 = result_Multipol_DD[:,2,1][cut]
	yer1 = result_Multipol_DD[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Data \, \xi_{0}$',color='blue')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Data \, \xi_{2}$',color='blue')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()








xxx = xi1D_[:,0]
yyy = xi1D_[:,1]
yer = xi1D_[:,2]



for rescale in [0,1,2]:
	if (rescale==0):
		plt.errorbar(xi1DD[:,0], xi1DD[:,1], yerr=xi1DD[:,2], fmt='o',color='blue',label='Data')
		plt.errorbar(xxx, yyy, yerr=yer, fmt='o',color='red',label='Simu')
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.errorbar(xi1DD[:,0], xi1DD[:,0]*xi1DD[:,1], yerr=xi1DD[:,0]*xi1DD[:,2], fmt='o',color='blue',label='Data')
		plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, fmt='o',color='red',label='Simu')
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.errorbar(xi1DD[:,0], xi1DD[:,0]*xi1DD[:,0]*xi1DD[:,1], yerr=xi1DD[:,0]*xi1DD[:,0]*xi1DD[:,2], fmt='o',color='blue',label='Data')
		plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, fmt='o',color='red',label='Simu')
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.show()









####################################################################################################
'''
saveListRealMocks(1,2)


rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/'
result_1D = numpy.load(rawPath+'xi_QSO_QSO_result_1D.npy')
result_2D = numpy.load(rawPath+'xi_QSO_QSO_result_2D.npy')
result_Mu = numpy.load(rawPath+'xi_QSO_QSO_result_Mu.npy')
result_We = numpy.load(rawPath+'xi_QSO_QSO_result_We.npy')
result_Multipol = numpy.load(rawPath+'xi_QSO_QSO_result_Multipol.npy')


rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/'

listBAO = numpy.load(rawPath+'xi_QSO_QSO_result_BAO.npy')
result_1D[:,1] = numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_1D.npy'),axis=1)
result_1D[:,2] = numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_QSO_QSO_result_cov_1D.npy')))/numpy.sqrt(100.)

result_2D[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_2D.npy'),axis=1), nbBinX2D__,nbBinY2D__)
result_2D[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_QSO_QSO_result_cov_2D.npy')))/numpy.sqrt(100.), nbBinX2D__,nbBinY2D__)

result_Mu[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_Mu.npy'),axis=1), nbBin1D__,nbBinM__)
result_Mu[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_QSO_QSO_result_cov_Mu.npy')))/numpy.sqrt(100.), nbBin1D__,nbBinM__)

result_Multipol[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_Multipol.npy'),axis=2)
result_Multipol[:,:,2] = numpy.var(numpy.load(rawPath+'xi_QSO_QSO_result_Multipol.npy'),axis=2)/numpy.sqrt(100.)

result_We[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_We.npy'),axis=2)
result_We[:,:,2] = numpy.var(numpy.load(rawPath+'xi_QSO_QSO_result_We.npy'),axis=2)/numpy.sqrt(100.)


### beta
plt.hist(listBAO[0,:],bins=20)
myTools.deal_with_plot(False,False,True)
plt.show()
plt.hist(listBAO[1,:],bins=20)
myTools.deal_with_plot(False,False,True)
plt.show()

#print fitCamb(result_1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
arg_fit = {}
arg_fit['a'] = 0.
arg_fit['b'] = 0.
arg_fit['c'] = 0.
arg_fit['A'] = 0.0031#*101.
arg_fit['sigma'] = 10.
arg_fit['mean'] = 105.
arg_fit['fix_A'] = False
arg_fit['fix_sigma'] = False
arg_fit['fix_mean'] = False
arg_fit['background'] = 'inv'
arg_fit['cov'] = 'full'
arg_fit['x_min'] = 60.
arg_fit['x_max'] = 160.
a,b,c,d,e,f,g =  myTools.fit_BAO(result_1D[:,0],result_1D[:,1],numpy.load(rawPath+'xi_QSO_QSO_result_cov_1D.npy')/100.,arg_fit)
print a
print b
plt.errorbar(result_1D[:,0],result_1D[:,0]**2.*result_1D[:,1],yerr=result_1D[:,0]**2.*result_1D[:,2],fmt='o')
plt.plot(c[0],c[0]**2.*c[1])
plt.plot(g[0],g[0]**2.*g[1])
plt.show()


### 2D:
plotXi2D(result_2D, 0)
plotXi2D(result_2D, 1)
plotXi2D(result_2D, 2)

### xiMu
plotMu(result_Mu,0)
plotMu(result_Mu,1)
plotMu(result_Mu,2)

### xiWe
plotWe(result_We,0)
plotWe(result_We,1)
plotWe(result_We,2)

### Multipol
result_1D[:,0] = result_Multipol[:,0,0]
result_1D[:,1] = result_Multipol[:,0,1]
result_1D[:,2] = result_Multipol[:,0,2]
plot_Xi_1D(result_1D, 0)
plot_Xi_1D(result_1D, 1)
plot_Xi_1D(result_1D, 2)
print fitCamb(result_1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
### Multipol
result_1D[:,0] = result_Multipol[:,2,0]
result_1D[:,1] = result_Multipol[:,2,1]
result_1D[:,2] = result_Multipol[:,2,2]
plot_Xi_1D(result_1D, 0)
plot_Xi_1D(result_1D, 1)
plot_Xi_1D(result_1D, 2)
print fitCamb(result_1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat',2)
### Multipol

for rescale in range(0,3):
	###
	cut = (result_Multipol[:,0,2]!=0.)
	xxx0 = result_Multipol[:,0,0][cut]
	yyy0 = result_Multipol[:,0,1][cut]
	yer0 = result_Multipol[:,0,2][cut]
	###
	cut = (result_Multipol[:,2,2]!=0.)
	xxx1 = result_Multipol[:,2,0][cut]
	yyy1 = result_Multipol[:,2,1][cut]
	yer1 = result_Multipol[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$\xi_{0}$')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$\xi_{2}$')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()
'''







################################################################################################


### Numbers
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_1D_LYA_QSO.txt')
data[:,2] /= data[:,5]
plt.errorbar(data[:,2],data[:,6],fmt='o',color='blue',label='Data')

data_SAVE = numpy.array(numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_delta_QSO_1D_LYA_QSO.txt'))
data_SAVE.fill(0.)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_delta_QSO_1D_LYA_QSO.txt'
		data = numpy.loadtxt(path)
		data[:,2] /= data[:,5]
		data_SAVE += data
data_SAVE /= 100.
plt.errorbar(data_SAVE[:,2],data_SAVE[:,6],fmt='o',color='red',label='<Simu>')
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$nb \, pairs \, QSO-LYA$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()



### Errors
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_1D_LYA_QSO.txt')
data[:,2] /= data[:,5]
data[:,3] = data[:,1]/data[:,5] - (data[:,0]/data[:,5])**2.
plt.errorbar(data[:,2],data[:,3],fmt='o',color='blue',label='Data')

data_SAVE = numpy.array(numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_delta_QSO_1D_LYA_QSO.txt'))
data_SAVE.fill(0.)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_delta_QSO_1D_LYA_QSO.txt'
		data = numpy.loadtxt(path)
		data[:,2] /= data[:,5]
		data[:,3] = data[:,1]/data[:,5] - (data[:,0]/data[:,5])**2.
		data_SAVE += data
data_SAVE /= 100.
plt.errorbar(data_SAVE[:,2],data_SAVE[:,3],fmt='o',color='red',label='<Simu>')

plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$Var \, QSO-LYA$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()


################################################################################################
### Numbers
data = numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery/xi_QSO_QSO_1D_QSO.txt')
data[:,1] /= data[:,0]
plt.errorbar(data[:,1],data[:,0],fmt='o',color='blue',label='Data')

data_SAVE = numpy.array(numpy.loadtxt('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_QSO_QSO_1D_QSO_DD.txt'))
data_SAVE.fill(0.)
for i in numpy.arange(10):
	for j in numpy.arange(10):
		
		path  = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/xi_QSO_QSO_1D_QSO_DD.txt'
		data = numpy.loadtxt(path)
		data[:,1] /= data[:,0]
		data_SAVE += data
data_SAVE /= 100.
plt.errorbar(data_SAVE[:,1],data_SAVE[:,0],fmt='o',color='red',label='<Simu>')
plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
plt.ylabel(r'$nb \, pairs \, QSO-QSO$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()



data = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/')



'''
i = sys.argv[5]
j = sys.argv[6]

path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Data/QSO_withRSD.fits'
cat = pyfits.open(path)[1].data
nd = cat.size
nr = cat.size
del cat
rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results/'
nbRand = 10
result_1D,result_2D,result_Mu,result_We,result_Multipol = getACorrelation(nd,nr,nbRand,rawPath)
'''


#saveListRealMocks(10,10)


rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/'
result_1D = numpy.load(rawPath+'xi_QSO_QSO_result_1D.npy')
result_2D = numpy.load(rawPath+'xi_QSO_QSO_result_2D.npy')
result_Mu = numpy.load(rawPath+'xi_QSO_QSO_result_Mu.npy')
result_We = numpy.load(rawPath+'xi_QSO_QSO_result_We.npy')
result_Multipol = numpy.load(rawPath+'xi_QSO_QSO_result_Multipol.npy')


rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/'

listBAO = numpy.load(rawPath+'xi_QSO_QSO_result_BAO.npy')
result_1D[:,1] = numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_1D.npy'),axis=1)
result_1D[:,2] = numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_QSO_QSO_result_cov_1D.npy')))/numpy.sqrt(100.)

result_2D[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_2D.npy'),axis=1), nbBinX2D__,nbBinY2D__)
result_2D[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_QSO_QSO_result_cov_2D.npy')))/numpy.sqrt(100.), nbBinX2D__,nbBinY2D__)

result_Mu[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_Mu.npy'),axis=1), nbBin1D__,nbBinM__)
result_Mu[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_QSO_QSO_result_cov_Mu.npy')))/numpy.sqrt(100.), nbBin1D__,nbBinM__)

result_Multipol[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_Multipol.npy'),axis=2)
result_Multipol[:,:,2] = numpy.var(numpy.load(rawPath+'xi_QSO_QSO_result_Multipol.npy'),axis=2)/numpy.sqrt(100.)

result_We[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_QSO_QSO_result_We.npy'),axis=2)
result_We[:,:,2] = numpy.var(numpy.load(rawPath+'xi_QSO_QSO_result_We.npy'),axis=2)/numpy.sqrt(100.)

'''
### beta
plt.hist(listBAO[0,:],bins=20)
myTools.deal_with_plot(False,False,True)
plt.show()
plt.hist(listBAO[1,:],bins=20)
myTools.deal_with_plot(False,False,True)
plt.show()

#print fitCamb(result_1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
arg_fit = {}
arg_fit['a'] = 0.
arg_fit['b'] = 0.
arg_fit['c'] = 0.
arg_fit['A'] = 0.0031#*101.
arg_fit['sigma'] = 10.
arg_fit['mean'] = 105.
arg_fit['fix_A'] = False
arg_fit['fix_sigma'] = False
arg_fit['fix_mean'] = False
arg_fit['background'] = 'inv'
arg_fit['cov'] = 'full'
arg_fit['x_min'] = 60.
arg_fit['x_max'] = 160.
a,b,c,d,e,f,g =  myTools.fit_BAO(result_1D[:,0],result_1D[:,1],numpy.load(rawPath+'xi_QSO_QSO_result_cov_1D.npy')/100.,arg_fit)
print a
print b
plt.errorbar(result_1D[:,0],result_1D[:,0]**2.*result_1D[:,1],yerr=result_1D[:,0]**2.*result_1D[:,2],fmt='o')
plt.plot(c[0],c[0]**2.*c[1])
plt.plot(g[0],g[0]**2.*g[1])
plt.show()
'''

### 2D:
plotXi2D(result_2D, 0)
plotXi2D(result_2D, 1)
plotXi2D(result_2D, 2)

### xiMu
plotMu(result_Mu,0)
plotMu(result_Mu,1)
plotMu(result_Mu,2)

### xiWe
plotWe(result_We,0)
plotWe(result_We,1)
plotWe(result_We,2)

### Multipol
result_1D[:,0] = result_Multipol[:,0,0]
result_1D[:,1] = result_Multipol[:,0,1]
result_1D[:,2] = result_Multipol[:,0,2]
plot_Xi_1D(result_1D, 0)
plot_Xi_1D(result_1D, 1)
plot_Xi_1D(result_1D, 2)
print fitCamb(result_1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
### Multipol
result_1D[:,0] = result_Multipol[:,2,0]
result_1D[:,1] = result_Multipol[:,2,1]
result_1D[:,2] = result_Multipol[:,2,2]
plot_Xi_1D(result_1D, 0)
plot_Xi_1D(result_1D, 1)
plot_Xi_1D(result_1D, 2)
print fitCamb(result_1D,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat',2)
### Multipol

for rescale in range(0,3):
	###
	cut = (result_Multipol[:,0,2]!=0.)
	xxx0 = result_Multipol[:,0,0][cut]
	yyy0 = result_Multipol[:,0,1][cut]
	yer0 = result_Multipol[:,0,2][cut]
	###
	cut = (result_Multipol[:,2,2]!=0.)
	xxx1 = result_Multipol[:,2,0][cut]
	yyy1 = result_Multipol[:,2,1][cut]
	yer1 = result_Multipol[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$\xi_{0}$')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$\xi_{2}$')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qq} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{qq} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qq} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$'+qso1__+' \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()



#####################################################################################################

'''
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/'
name = path
saveListReal(100,name+'xi_delta_QSO_Mu_LYA_QSO_shuffleQSO_',name+'xi_delta_QSO_2D_LYA_QSO_shuffleQSO_',name+'shuffleQSO_LYA_QSO_',False)


myTools.plot2D( numpy.cov(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/shuffleQSO_LYA_QSO_1D.npy')) )
myTools.plot2D( numpy.corrcoef(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/shuffleQSO_LYA_QSO_1D.npy')) )
plotCovar([numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/shuffleQSO_LYA_QSO_cov_1D.npy')], ['a'])
'''
#saveListRealMocks(10,10)

xi1DD, xi2DD, xiMuDD, xiWeDD = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_2D_LYA_QSO.txt')
result_Multipol_DD = plotMultipol(xiMuDD)

path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
tmpPath =  path + '0/Simu_000/Results_RandomPosInCell/xi_delta_QSO_'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(tmpPath+'Mu_LYA_QSO.txt',tmpPath+'2D_LYA_QSO.txt')
result_Multipol_ = plotMultipol(xiMu_)
replaceValueByMean()




plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')



plotXi2D(0)
plotXi2D(1)
plotXi2D(2)

plotMu(0)
plotMu(1)
plotMu(2)

plotWe(0)
plotWe(1)
plotWe(2)

### Multipol
xi1D_[:,0] = result_Multipol_[:,0,0]
xi1D_[:,1] = result_Multipol_[:,0,1]
xi1D_[:,2] = result_Multipol_[:,0,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
### Multipol
xi1D_[:,0] = result_Multipol_[:,2,0]
xi1D_[:,1] = result_Multipol_[:,2,1]
xi1D_[:,2] = result_Multipol_[:,2,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat',2)
### Multipol

for rescale in range(0,3):
	###
	cut = (result_Multipol_[:,0,2]!=0.)
	xxx0 = result_Multipol_[:,0,0][cut]
	yyy0 = result_Multipol_[:,0,1][cut]
	yer0 = result_Multipol_[:,0,2][cut]
	###
	cut = (result_Multipol_[:,2,2]!=0.)
	xxx1 = result_Multipol_[:,2,0][cut]
	yyy1 = result_Multipol_[:,2,1][cut]
	yer1 = result_Multipol_[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Simu \, \xi_{0}$',color='red')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Simu \, \xi_{2}$',color='red')

	###
	cut = (result_Multipol_DD[:,0,2]!=0.)
	xxx0 = result_Multipol_DD[:,0,0][cut]
	yyy0 = result_Multipol_DD[:,0,1][cut]
	yer0 = result_Multipol_DD[:,0,2][cut]
	###
	cut = (result_Multipol_DD[:,2,2]!=0.)
	xxx1 = result_Multipol_DD[:,2,0][cut]
	yyy1 = result_Multipol_DD[:,2,1][cut]
	yer1 = result_Multipol_DD[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Data \, \xi_{0}$',color='blue')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Data \, \xi_{2}$',color='blue')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()




#####################################################################################################

#prepareForBAOFIT()

### Data
### Data
path = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/'
saveListReal(100,path+'xi_delta_QSO_Mu_LYA_QSO_shuffleForest_',path+'xi_delta_QSO_2D_LYA_QSO_shuffleForest_',path+'shuffleForest_LYA_QSO_',False)
plotCovarDifferentMethod()

'''
list1D = numpy.zeros( shape=(nbBin2D__,80*100) )
m = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results_RandomPosInCell/subSampling_LYA_QSO_2D.npy')
idx = 0
for k in range(0,80):
	list1D[:,idx] = m[:,k]
	idx += 1

for i in range(0,10):
	for j in range(0,10):
		if (i==0 and j==0): continue
		print i, j
		m = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'+str(i)+'/Simu_00'+str(j)+'/Results_RandomPosInCell/subSampling_LYA_QSO_2D.npy')
		for k in range(0,80):
			list1D[:,idx] = m[:,k]
			idx += 1
print idx
print list1D[0,:].size

cov = numpy.cov(list1D)
cor = numpy.corrcoef(list1D)
numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_allSubSampling.npy',cov)
numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSampling.npy',cor)


cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_allSubSampling.npy')
cor = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSampling.npy')

myTools.plot2D(cov)
myTools.plot2D(cor)
covFromFit = myTools.plotCovar([cov],['o'])[0]
corFromFit = myTools.getCorrelationMatrix(covFromFit)

numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_allSubSamplingFromFit.npy',covFromFit)
numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy',corFromFit)


myTools.plot2D(covFromFit)
myTools.plot2D(corFromFit)

myTools.plot2D( covFromFit-cov )
myTools.plot2D( corFromFit-cor )

a = corFromFit-cor
plt.hist( a.flatten(), bins=100 )
plt.show()
'''

#####################################################################################################















'''
listXi1D = []
name = []

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_Mu_LYA_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_2D_LYA_QSO.txt')
listXi1D += [xi1D_]
name += ['LYA-QSO']

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_Mu_CIV_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_2D_CIV_QSO.txt')
listXi1D += [xi1D_]
name += ['CIV-QSO']

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_Mu_SIIV_QSO.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_2D_SIIV_QSO.txt')
listXi1D += [xi1D_]
name += ['SiIV-QSO']

xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/CIV_in_SIIV/xi_delta_QSO_Mu_SIIV_QSO_DR7_DR12_EBOSS.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/CIV_in_SIIV/xi_delta_QSO_2D_SIIV_QSO_DR7_DR12_EBOSS.txt')
listXi1D += [xi1D_]
name += ['CIV-QSO (in SiIV forest)']


xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/CIV_in_LYA/xi_delta_QSO_Mu_LYA_QSO_DR7_DR12_EBOSS.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/CIV_in_LYA/xi_delta_QSO_2D_LYA_QSO_DR7_DR12_EBOSS.txt')
listXi1D += [xi1D_]
name += ['CIV-QSO (in LYA forest)']


for rescale in [0,1,2]:
	for i in range(0, len(name)):

		xxx = listXi1D[i][:,0]
		yyy = listXi1D[i][:,1]
		yer = listXi1D[i][:,2]

		if (i!=0):
			if (rescale==0): print name[i], 1./ (listXi1D[0][1,1]/yyy[1])
			yyy = yyy*listXi1D[0][1,1]/yyy[1]
			yer = yer*listXi1D[0][1,1]/yyy[1]
		else:
			if (rescale==0): print name[i], yyy[1]

		if (rescale==0):
			plt.errorbar(xxx, yyy, yerr=yer, marker='o',label=name[i])
			plt.ylabel(r'$\\xi^{qf} (|s|)$', fontsize=40)
		if (rescale==1):
			plt.errorbar(xxx, xxx*yyy, yerr=xxx*yer, marker='o',label=name[i])
			plt.ylabel(r'$|s|.\\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
		if (rescale==2):
			plt.errorbar(xxx, xxx*xxx*yyy, yerr=xxx*xxx*yer, marker='o',label=name[i])
			plt.ylabel(r'$|s|^{2}.\\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)

	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx)-10., numpy.max(xxx)+10. ])
	plt.legend(fontsize=30, frameon=False, numpoints=1,ncol=2,loc=4)
	plt.show()
'''





























xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_Mu_CIV_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_2D_CIV_'+qso1__+'.txt')
plotXi(0)
plotXi(1)
plotXi(2)

#pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
#fitCamb(xi1D_,pathToCamb,0)

plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)

#prepareForBAOFIT()


'''
cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_allSubSampling.npy')
cor = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSampling.npy')

myTools.plot2D(cov)
myTools.plot2D(cor)
covFromFit = myTools.plotCovar([cov],['o'])[0]
corFromFit = myTools.getCorrelationMatrix(covFromFit)

myTools.plot2D(covFromFit)
myTools.plot2D(corFromFit)

myTools.plot2D( covFromFit-cov )
myTools.plot2D( corFromFit-cor )

a = corFromFit-cor
plt.hist( a.flatten(), bins=100 )
plt.show()

numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cov_2D_allSubSamplingFromFit.npy',covFromFit)
numpy.save('/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results_RandomPosInCell/xi_delta_QSO_result_cor_2D_allSubSamplingFromFit.npy',corFromFit)
'''






model-config = value[beta]=1.1;
model-config = fix[(1+beta)*bias]=-0.336;
model-config = fix[gamma-bias]=0.9; fix[gamma-beta]=0;
model-config = value[bias2]=3.64;
model-config = fix[beta2*bias2]=0.962524;
model-config = value[delta-v]=0;
model-config = fix[BAO amplitude]=1;
model-config = fix[BAO alpha-iso];
model-config = value[BAO alpha-parallel]=1;
model-config = value[BAO alpha-perp]=1;
model-config = fix[gamma-scale]=0;





#####################################################################################################




















xi1D_, xi2D_, xiMu_, xiWe_ = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')

cov = numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/subsampling_LYA_QSO_cov_1D.npy')
print xi1D_[:,2]
print numpy.sqrt(numpy.diag(cov))
xi1D_[:,2] = numpy.sqrt(numpy.diag(cov))

print xiWe_[:,:,2]
print numpy.sqrt(numpy.var(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/subsampling_LYA_ALL_OBJECTS_We.npy'),axis=2)/80.)
xiWe_[:,:,2] = numpy.sqrt(numpy.var(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/subsampling_LYA_ALL_OBJECTS_We.npy'),axis=2)/80.)


plotXi(0)
plotXi(1)
plotXi(2)
pathToCamb = '/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat'
fitCamb(xi1D_,pathToCamb,0)
plotXi2D(0)
plotXi2D(1)
plotXi2D(2)
plotMu(0)
plotMu(1)
plotMu(2)
plotWe(0)
plotWe(1)
plotWe(2)




xi1DD, xi2DD, xiMuDD, xiWeDD = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')
result_Multipol_DD = plotMultipol(xiMuDD)

result_Multipol_DD[:,:,2] = numpy.sqrt(numpy.var(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/subsampling_LYA_QSO_Multipol.npy'),axis=2)/80.)


path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00'
tmpPath =  path + '0/Simu_000/Results/xi_delta_QSO_'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(tmpPath+'Mu_LYA_QSO.txt',tmpPath+'2D_LYA_QSO.txt')
result_Multipol = plotMultipol(xiMu_)


rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Results/'

xi1D_[:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_1D.npy'),axis=1)
xi1D_[:,2] = numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_1D.npy')))/numpy.sqrt(100.)

print numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_2D.npy'),axis=1).size
xi2D_[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_2D.npy'),axis=1), nbBinX2D__,nbBinY2D__)
xi2D_[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_2D.npy')))/numpy.sqrt(100.), nbBinX2D__,nbBinY2D__)

xiMu_[:,:,2] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_Mu.npy'),axis=1), nbBin1D__,nbBinM__)
xiMu_[:,:,3] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_Mu.npy')))/numpy.sqrt(100.), nbBin1D__,nbBinM__)

result_Multipol[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_Multipol.npy'),axis=2)
result_Multipol[:,:,2] = numpy.sqrt(numpy.var(numpy.load(rawPath+'xi_delta_QSO_result_Multipol.npy'),axis=2)/100.)

xiWe_[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_We.npy'),axis=2)
xiWe_[:,:,2] = numpy.sqrt(numpy.var(numpy.load(rawPath+'xi_delta_QSO_result_We.npy'),axis=2)/100.)



plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')



plotXi2D(0)
plotXi2D(1)
plotXi2D(2)

plotMu(0)
plotMu(1)
plotMu(2)

plotWe(0)
plotWe(1)
plotWe(2)

### Multipol
xi1D_[:,0] = result_Multipol[:,0,0]
xi1D_[:,1] = result_Multipol[:,0,1]
xi1D_[:,2] = result_Multipol[:,0,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
### Multipol
xi1D_[:,0] = result_Multipol[:,2,0]
xi1D_[:,1] = result_Multipol[:,2,1]
xi1D_[:,2] = result_Multipol[:,2,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat',2)
### Multipol

for rescale in range(0,3):
	###
	cut = (result_Multipol[:,0,2]!=0.)
	xxx0 = result_Multipol[:,0,0][cut]
	yyy0 = result_Multipol[:,0,1][cut]
	yer0 = result_Multipol[:,0,2][cut]
	###
	cut = (result_Multipol[:,2,2]!=0.)
	xxx1 = result_Multipol[:,2,0][cut]
	yyy1 = result_Multipol[:,2,1][cut]
	yer1 = result_Multipol[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Simu \, \xi_{0}$',color='red')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Simu \, \xi_{2}$',color='red')

	###
	cut = (result_Multipol_DD[:,0,2]!=0.)
	xxx0 = result_Multipol_DD[:,0,0][cut]
	yyy0 = result_Multipol_DD[:,0,1][cut]
	yer0 = result_Multipol_DD[:,0,2][cut]
	###
	cut = (result_Multipol_DD[:,2,2]!=0.)
	xxx1 = result_Multipol_DD[:,2,0][cut]
	yyy1 = result_Multipol_DD[:,2,1][cut]
	yer1 = result_Multipol_DD[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Data \, \xi_{0}$',color='blue')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Data \, \xi_{2}$',color='blue')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()



##########################################################################################################################################




saveListRealMocks(10,10)

xi1DD, xi2DD, xiMuDD, xiWeDD = loadData('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_Mu_'+forest1__+'_'+qso1__+'.txt','/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/xi_delta_QSO_2D_'+forest1__+'_'+qso1__+'.txt')
result_Multipol_DD = plotMultipol(xiMuDD)

result_Multipol_DD[:,:,2] = numpy.sqrt(numpy.var(numpy.load('/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_testNoCutLambdaOBS/subsampling_LYA_QSO_Multipol.npy'),axis=2)/80.)


path       = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_00'
tmpPath =  path + '0/Simu_000/Results/xi_delta_QSO_'
xi1D_, xi2D_, xiMu_, xiWe_ = loadData(tmpPath+'Mu_LYA_QSO.txt',tmpPath+'2D_LYA_QSO.txt')
result_Multipol = plotMultipol(xiMu_)


rawPath = '/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Results/'

xi1D_[:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_1D.npy'),axis=1)
xi1D_[:,2] = numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_1D.npy')))/numpy.sqrt(100.)

print numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_2D.npy'),axis=1).size
xi2D_[:,:,1] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_2D.npy'),axis=1), nbBinX2D__,nbBinY2D__)
xi2D_[:,:,2] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_2D.npy')))/numpy.sqrt(100.), nbBinX2D__,nbBinY2D__)

xiMu_[:,:,2] = myTools.convert1DTo2D( numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_Mu.npy'),axis=1), nbBin1D__,nbBinM__)
xiMu_[:,:,3] = myTools.convert1DTo2D( numpy.sqrt(numpy.diag(numpy.load(rawPath+'xi_delta_QSO_result_cov_Mu.npy')))/numpy.sqrt(100.), nbBin1D__,nbBinM__)

result_Multipol[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_Multipol.npy'),axis=2)
result_Multipol[:,:,2] = numpy.sqrt(numpy.var(numpy.load(rawPath+'xi_delta_QSO_result_Multipol.npy'),axis=2)/100.)

xiWe_[:,:,1] = numpy.mean(numpy.load(rawPath+'xi_delta_QSO_result_We.npy'),axis=2)
xiWe_[:,:,2] = numpy.sqrt(numpy.var(numpy.load(rawPath+'xi_delta_QSO_result_We.npy'),axis=2)/100.)



a = numpy.load(rawPath+'xi_delta_QSO_result_1D.npy')
for i in range(0,a[0,:].size):
	plt.plot(numpy.arange(0,50),a[:,i])
plt.show()
for i in range(0,a[0,:].size):	
	x = numpy.arange(0,50)*4.+2.
	plt.plot(numpy.arange(0,50),x*a[:,i])
plt.show()
for i in range(0,a[0,:].size):
	x = (numpy.arange(0,50)*4.+2.)**2.
	plt.plot(numpy.arange(0,50),x*a[:,i])
plt.show()

plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')



plotXi2D(0)
plotXi2D(1)
plotXi2D(2)

plotMu(0)
plotMu(1)
plotMu(2)

plotWe(0)
plotWe(1)
plotWe(2)

### Multipol
xi1D_[:,0] = result_Multipol[:,0,0]
xi1D_[:,1] = result_Multipol[:,0,1]
xi1D_[:,2] = result_Multipol[:,0,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat')
### Multipol
xi1D_[:,0] = result_Multipol[:,2,0]
xi1D_[:,1] = result_Multipol[:,2,1]
xi1D_[:,2] = result_Multipol[:,2,2]
plotXi(0)
plotXi(1)
plotXi(2)
print fitCamb(xi1D_,'/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/CAMB_2_4/xi-z2.4.dat',2)
### Multipol

for rescale in range(0,3):
	###
	cut = (result_Multipol[:,0,2]!=0.)
	xxx0 = result_Multipol[:,0,0][cut]
	yyy0 = result_Multipol[:,0,1][cut]
	yer0 = result_Multipol[:,0,2][cut]
	###
	cut = (result_Multipol[:,2,2]!=0.)
	xxx1 = result_Multipol[:,2,0][cut]
	yyy1 = result_Multipol[:,2,1][cut]
	yer1 = result_Multipol[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Simu \, \xi_{0}$',color='red')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Simu \, \xi_{2}$',color='red')

	###
	cut = (result_Multipol_DD[:,0,2]!=0.)
	xxx0 = result_Multipol_DD[:,0,0][cut]
	yyy0 = result_Multipol_DD[:,0,1][cut]
	yer0 = result_Multipol_DD[:,0,2][cut]
	###
	cut = (result_Multipol_DD[:,2,2]!=0.)
	xxx1 = result_Multipol_DD[:,2,0][cut]
	yyy1 = result_Multipol_DD[:,2,1][cut]
	yer1 = result_Multipol_DD[:,2,2][cut]

	coef0 = numpy.power(xxx0,rescale)
	coef1 = numpy.power(xxx1,rescale)

	plt.errorbar(xxx0, coef0*yyy0, yerr=coef0*yer0, marker='o', label=r'$Data \, \xi_{0}$',color='blue')
	plt.errorbar(xxx1, coef1*yyy1, yerr=coef1*yer1, marker='o', label=r'$Data \, \xi_{2}$',color='blue')

	if (rescale==0):
		plt.ylabel(r'$\xi^{qf} (|s|)$', fontsize=40)
	if (rescale==1):
		plt.ylabel(r'$|s|.\xi^{qf} (|s|) \, [h^{-1}.Mpc]$', fontsize=40)
	if (rescale==2):
		plt.ylabel(r'$|s|^{2}.\xi^{qf} (|s|) \, [(h^{-1}.Mpc)^{2}]$', fontsize=40)
	
	plt.title(r'$\delta_{'+forest1__+'} \, - \, '+qso1__+'$', fontsize=40)
	plt.xlabel(r'$|s| \, [h^{-1}.Mpc]$', fontsize=40)
	myTools.deal_with_plot(False,False,True)
	plt.xlim([ numpy.min(xxx0)-10., numpy.max(xxx0)+10. ])
	#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.show()























































