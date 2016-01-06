### Var_delta_pow_-1 vs. log10_delta_ivar
		xxx, yyy, eyyy, nyyy = Get_TProfile(log_delta_ivar, ar_delta, deltaIvarBinEdges__,ar_weight)

		yyy   = 1./((eyyy**2.)*nyyy)
		eyyy  = 1./(yyy*numpy.sqrt(nyyy))
		varPipeline.append( ( xxx, yyy, eyyy ) )
		
		ar_cut = numpy.logical_and( (xxx>-1.5), (xxx<2.3) )
		xxx  = xxx[ ar_cut ]
		yyy  = yyy[ ar_cut ]
		eyyy = eyyy[ ar_cut ]

		### Set minuit
		def chi2(TMPeta,TMPsigma2LSS):
			return numpy.sum( numpy.power((yyy-1./(TMPsigma2LSS+1./(TMPeta*numpy.power(10,xxx))))/eyyy ,2.) )
		m = Minuit(chi2, TMPeta=eta[loopIdx], error_TMPeta=0.1,  TMPsigma2LSS=sigma2LSS[loopIdx], error_TMPsigma2LSS=0.1, print_level=-1, errordef=0.01) 	
		m.migrad()
			
		eta[loopIdx+1]       = m.values['TMPeta']
		sigma2LSS[loopIdx+1] = m.values['TMPsigma2LSS']












		yy, xx = numpy.meshgrid(redshiftBinCenter__,lambdaRFTemplateBinCenter__)

		xx = xx[ (template2D[loopIdx][2]>0.) ]
		yy = yy[ (template2D[loopIdx][2]>0.) ]
		zz = template2D[loopIdx][0][ (template2D[loopIdx][2]>0.) ]
		print zip(xx, yy, zz)

		plt.plot(xx, zz, marker="o")
		#plt.plot(cat[0]['LAMBDA_RF'], cat[0]['NORM_FLUX'], marker="o")
		plt.show()
		plt.plot(yy, zz, marker="o")
		plt.show()
		'''
		xx = xx.ravel()
		yy = yy.ravel()
		zz = template2D[loopIdx][0].ravel()
		nb = template2D[loopIdx][2].ravel()

		xx = xx[ (nb>0.) ]
		yy = yy[ (nb>0.) ]
		zz = zz[ (nb>0.) ]

		print xx, yy, zz
		'''
		f = interpolate.interp2d(xx,yy,zz)











		
		#iterp = scipy.interpolate.griddata( (cat[0]['LAMBDA_RF'], 2.), cat[0]['NORM_FLUX'], (lambdaRFTemplateBinEdges__, 2.), method='linear')
		#plt.plot(cat[0]['LAMBDA_RF'], cat[0]['NORM_FLUX'])
		#plt.plot(cat[0]['LAMBDA_RF'], iterp)
		#plt.show()
		from numpy.random import uniform, seed
		def f(x,y):
			return x+y

		seed(1234)

		x = cat['LAMBDA_RF'][ (cat['LAMBDA_RF']!=0.) ]
		y = uniform(2.,6.,416)
		z = uniform(2.,6.,416*415)
		#print z
		#xi, yi = numpy.mgrid[1000.:1200., -10:10.]
		xi, yi = numpy.meshgrid(lambdaRFTemplateBinCenter__,redshiftBinCenter__)
		print xi
		print xi.size
		print yi
		print yi.size
		#print xi
		zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='cubic')
		#print zi
		plt.imshow(zi)
		plt.show()




