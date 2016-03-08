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
#include "TGraph.h"

// Numerical constants for faster calculation 
const double cOver100_ = 1.*C_C_LIGHT/1000/100.;

Cosmology::Cosmology(double h, double omegaM, double omegaB) {

	// In this cosmology we put omegaL_ = 1.-omegaM_

	h__ = h;
	omegaM__ = omegaM;
	omegaB__ = omegaB;
	omegaL__ = 1.-omegaM__;

	dChidz__ = NULL;
	hConvertRedshDist__ = NULL;
}
Cosmology::~Cosmology() {
	delete dChidz__;
}

TH1D* Cosmology::createHistoConvertRedshDist(unsigned int nbBinRedsh, double zExtremaBinConvert0, double zExtremaBinConvert1) {
	
	hConvertRedshDist__ = new TH1D("hConvertRedshDist__","",nbBinRedsh,zExtremaBinConvert0,zExtremaBinConvert1);
	R_dealWithPlots_1D(hConvertRedshDist__, "z", "r (Mpc.h^{-1})", " Distance versus z ");
	
	dChidz__ = new TF1("dChidz__", this,&Cosmology::dChidzfunct,zExtremaBinConvert0,zExtremaBinConvert1,0,"Cosmology","dChidzfunct");
	TF1* dt = new TF1("dt", this, &Cosmology::DTfunct, zExtremaBinConvert0,zExtremaBinConvert1,0,"Cosmology","DTfunct");
	
	const double ratio = zExtremaBinConvert1/nbBinRedsh;
	for (unsigned int i=0; i<nbBinRedsh; i++) {
		hConvertRedshDist__->SetBinContent(i+1,dt->Eval((i+.5)*ratio));
	}
	delete dt;

	return hConvertRedshDist__;
}

double Cosmology::dChidzfunct(double* zz, double* par) {

	const double z = zz[0];
	return cOver100_/sqrt(pow((1.+z),3.)*omegaM__ + omegaL__);
}
double Cosmology::Dzfunct(double* zz, double* par) {

	const double z = zz[0];
	return cOver100_/sqrt(pow((1.+z),3.)*omegaM__ + omegaL__)*(1.+z)*(1.+z);
}
double Cosmology::DAfunct(double* zz, double* par) {

	const double z = zz[0];
	return dChidz__->Integral(0.,z)/(1.+z);
}
double Cosmology::DTfunct(double* zz, double* par) {

	const double z = zz[0];
	return dChidz__->Integral(0.,z);
}
void Cosmology::FindMinMaxRedsift(double maxDist, double lambdaObsMin, double lambdaObsMax, double lambdaRFLine, double* array) {

	// Maximal distance between pairs
	const double dist_max = maxDist*sqrt(2.);

	// Redsift min and max of pixels
	const double z_pixel_min = (lambdaObsMin/lambdaRFLine)-1.;
	const double z_pixel_max = (lambdaObsMax/lambdaRFLine)-1.;

	if (hConvertRedshDist__==NULL) {
		std::cout << "  Cosmology::FindMinMaxRedsift::  ERROR, hConvertRedshDist__==NULL " << std::endl;
		return;
	}

	// Distance min and max of pixels
	const double d_pixel_min = hConvertRedshDist__->Interpolate(z_pixel_min);
	const double d_pixel_max = hConvertRedshDist__->Interpolate(z_pixel_max);

	// Distance min and max of QSOs
	const double d_QSO_min = d_pixel_min-dist_max;
	const double d_QSO_max = d_pixel_max+dist_max;

	// Histo to convert distance into redshift
	const unsigned int nbBin = hConvertRedshDist__->GetNbinsX();
	double redShift[nbBin];
	double distance[nbBin];
	for (unsigned int i=0; i<nbBin; i++) {
		redShift[i] = hConvertRedshDist__->GetBinCenter(i+1);
		distance[i] = hConvertRedshDist__->GetBinContent(i+1);
	}
	TGraph* gFromDistanceToRedshfit = new TGraph(nbBin, distance, redShift);

	// Redshift min and max of QSOs
	array[0] = std::floor( 10.*gFromDistanceToRedshfit->Eval(d_QSO_min))/10.;
	array[1] = std::ceil( 10.*gFromDistanceToRedshfit->Eval(d_QSO_max))/10.;

	std::cout << "  min_z_qso = " << array[0] << "  , max_z_qso = " << array[1] << std::endl;

	return;
}
double Cosmology::GetMinDistPixels(double lambdaObsMin, double lambdaRFLine) {

	if (hConvertRedshDist__==NULL) {
		std::cout << "  Cosmology::GetMinDistPixels::  ERROR, hConvertRedshDist__==NULL " << std::endl;
		return 0.;
	}

	const double z_pixel_min = (lambdaObsMin/lambdaRFLine)-1.;
	const double d_pixel_min = std::floor( hConvertRedshDist__->Interpolate(z_pixel_min));
	
	std::cout << "\n\n  d_pixel_min = " << d_pixel_min << "\n" << std::endl;
	return d_pixel_min;
}



























