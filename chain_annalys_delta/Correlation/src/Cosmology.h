//===================================================================================
//
//         FILE: Cosmology.h
//
//        USAGE: ---
//
//  DESCRIPTION: Cosmology function. Given a cosmology, has a function that returns
//		a conversion table from redshift to distance in h^{-1}.Mpc
//
//      OPTIONS: ---
// REQUIREMENTS: ---
//         BUGS: ---
//        NOTES: ---
//       AUTHOR: Hélion du Mas des Bourboux, Helion.du-Mas-des-Bourboux@cea.fr
//      COMPANY: CEA (France)
//      VERSION: ---
//      CREATED: ---
//     REVISION: ---
//===================================================================================

#include "TH1.h"
#include "TF1.h"

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

class Cosmology
{
	private:

		// Cosmology:
		// Reduced hubble constant 
		double h__;
		// Proportion of dark energy
		double omegaL__;
		// Proportion of matter (all type)
		double omegaM__;
		// Proportion of baryonic matter
		double omegaB__;

		// Estimate d(Chi(z))/dz
		//	in h^{-1}.Mpc
		TF1* dChidz__;

		// histogram to convert redshift into comoving distance in h^{-1}.Mpc
		TH1D* hConvertRedshDist__;

	public:

		// Class constructor:
		//		- double h      : Reduced hubble constant 
		//		- double omegaM : Proportion of matter (all type)
		//		- double omegaB : Proportion of baryonic matter
		Cosmology(double h, double omegaM, double omegaB);

		// Class destructor:
		~Cosmology();
		
		// Return an histogram to convert redshift into comoving distance in h^{-1}.Mpc
		//		- unsigned int nbBinRedsh    : number of redshift bin
		//		- double zExtremaBinConvert0 : minimum of redshift
		//		- double zExtremaBinConvert1 : maximum of redshift
		// return TH1D* : histogram to convert redshift into comoving distance in h^{-1}.Mpc
		TH1D* createHistoConvertRedshDist(unsigned int nbBinRedsh, double zExtremaBinConvert0, double zExtremaBinConvert1);

		// Estimate d(Chi(z))/dz
		//	in h^{-1}.Mpc
		double dChidzfunct(double* zz, double *par);

		// Radial distance Dz*Deltanu/nu0
		//	in h^{-1}.Mpc
		double Dzfunct(double* zz, double *par);

		// Angular distance
		//	in h^{-1}.Mpc
		double DAfunct(double* zz, double *par);

		// Transverse distance
		//	in h^{-1}.Mpc
		double DTfunct(double* zz, double *par);

		// Fill an array with redshift min and max of QSO according to lambda_obs_min and lambda_obs_max of pixels
		//		- double maxDist      : maximum distance of correlation in h^{-1}.Mpc
		//		- double lambdaObsMin : minimum of lambda_obs
		//		- double lambdaObsMax : maximum of lambda_obs
		//		- double lambdaRFLine : lambda_RF of the line
		//		- double* array       : array of size two to fill
		void FindMinMaxRedsift(double maxDist, double lambdaObsMin, double lambdaObsMax, double lambdaRFLine, double* array);

		// Returns the distance minimum of pixels according to lambda_obs_min	
		//		- double lambdaObsMin : minimum of lambda_obs
		//		- double lambdaRFLine : lambda_RF of the line
		double GetMinDistPixels(double lambdaObsMin, double lambdaRFLine);

};

#endif

















