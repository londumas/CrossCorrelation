//===================================================================================
//
//         FILE: rootFunctionsCorrelat.h
//
//        USAGE: #include "/home/helion/Documents/Thèse/Code/Root/Library/rootFunctionsCorrelat.h"
//               #include "/home/usr201/mnt/hdumasde/Code/Root/Library/rootFunctionsCorrelat.h"
//               In a file where you need these functions
//
//  DESCRIPTION: Gathering of functions dealing with root objects
//               for the cross correlations
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

//#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

#ifndef ROOTFUNCTIONSCORRELAT_H
#define ROOTFUNCTIONSCORRELAT_H
	
	void RC_rescale1D(TH1* histo, std::vector<double> meanDist);
	void RC_rescale2D(TH2* histo, std::vector<double> meanRparal, std::vector<double> meanRperp);
	void RC_rescale2DPolar(TH2* histo);
	void RC_rescaleCovariance(TH2* histo);
	
	TGraphErrors* RC_getXCorrelatTGraphFromTH1(TH1D* histo, std::vector<double> meanDist);
	TH1D* RC_getXCorrelat1DFrom2D(TH2* histo);
	
	TH1D* RC_getCrossCorrelatConstantR(TH2* histo, unsigned int idxBinR, double minVar);
	TH1D** RC_getCorrelatMultipol(TH2* histo, std::string variable, double minVar);
	
	void RC_getLegendrePoly(TH2* histo);

#endif


