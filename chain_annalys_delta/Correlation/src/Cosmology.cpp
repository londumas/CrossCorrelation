//===================================================================================
//
//         FILE: Cosmology.cpp
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


#include "Cosmology.h"
#include "../../../Root/Library/RootHistoFunctions.h"
#include "../../../Constants/constants.h"

#include <cmath>

// Numerical constants for faster calculation 
const double cOver100_ = 1.*C_C/100.;

Cosmology::Cosmology(double h, double omegaM, double omegaB) {

	// In this cosmology we put omegaL_ = 1.-omegaM_

	h__ = h;
	omegaM__ = omegaM;
	omegaB__ = omegaB;
	omegaL__ = 1.-omegaM__;
}
Cosmology::~Cosmology() {
	delete dChidz__;
}

TH1D* Cosmology::createHistoConvertRedshDist(unsigned int nbBinRedsh, double zExtremaBinConvert0, double zExtremaBinConvert1) {
	
	TH1D* hConvertRedshDist = new TH1D("hConvertRedshDist","",nbBinRedsh,zExtremaBinConvert0,zExtremaBinConvert1);
	R_dealWithPlots_1D(hConvertRedshDist, "z", "r (Mpc.h^{-1})", " Distance versus z ");
	
	dChidz__ = new TF1("dChidz_", this,&Cosmology::dChidzfunct,zExtremaBinConvert0,zExtremaBinConvert1,0,"Cosmology","dChidzfunct");
	TF1* dt = new TF1("dt", this, &Cosmology::DTfunct, zExtremaBinConvert0,zExtremaBinConvert1,0,"Cosmology","DTfunct");
	
	const double ratio = zExtremaBinConvert1/nbBinRedsh;
	for (unsigned int i=0; i<nbBinRedsh; i++) {
		hConvertRedshDist->SetBinContent(i+1,dt->Eval((i+.5)*ratio));
	}
	delete dt;
	
	return hConvertRedshDist;
}

double Cosmology::dChidzfunct(double* zz, double *par) {

	const double z = zz[0];
	return cOver100_/sqrt(pow((1.+z),3.)*omegaM__ + omegaL__);
}
double Cosmology::Dzfunct(double* zz, double *par) {

	const double z = zz[0];
	return cOver100_/sqrt(pow((1.+z),3.)*omegaM__ + omegaL__)*(1.+z)*(1.+z);
}
double Cosmology::DAfunct(double* zz, double *par) {

	const double z = zz[0];
	return dChidz__->Integral(0.,z)/(1.+z);
}
double Cosmology::DTfunct(double* zz, double *par) {

	const double z = zz[0];
	return dChidz__->Integral(0.,z);
}


















