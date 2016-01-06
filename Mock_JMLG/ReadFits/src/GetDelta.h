//===================================================================================
//
//         FILE: CrossCorrelation.h
//
//        USAGE: ---
//
//  DESCRIPTION: Do the cross correlation between Lya forest and Qso
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

#include "Constants.h"

#include <fstream>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TProfile.h"

#ifndef GETDELTA_H
#define GETDELTA_H

/// Constants
const unsigned int nbBinRFMin__ = 50;
const unsigned int nbPixelsTemplate__ = 645;
const double lambdaRFMin__ = 1040.;
const double lambdaRFMax__ = 1200.;
const double lambdaRFLine__  = 1215.67;
const double lambdaRFTemplateMin__ = lambdaRFMin__-3.;
const double lambdaRFTemplateMax__ = lambdaRFMax__+3.;
const double lambdaRFNormaMin__ = 1275.;
const double lambdaRFNormaMax__ = 1285.;
const double lambdaObsMin__ = 3600.;
const double lambdaObsMax__ = 7235.;
const unsigned int nbBinlambdaObs__ = 3635;
const double c_speedOfLight__ = 2.998e5;

const double minRedshift__ = 1.96;
const double maxRedshift__ = 3.9;
const unsigned int nbBinsRedshift__ = 60;
const double minFlux__ = 0.;
const double maxFlux__ = 1.000001;
const unsigned int nbBinsFlux__ = 100;


/*
///For Christophe
h=0.71,
Olambda=0.73  
W0=-1
Omatter=0.267804  
Obaryon=0.0444356
Orelat=7.9e-05
Otot=0.997883  Ocurv=0.002117

///For Jean-Marc
h=0.71
Omegam = 0.27
OLambda=0.73
Or neglected
c_v = 2.998e08
const double lambdaRFLine__     = 1216.24;
*/

class GetDelta
{
	private:

		double ratioForestToQSO__;

		std::string pathToData__;
		std::string pathToSave__;
		bool noCont__;
		bool noNoise__;
		bool isTest__;
		bool withRSD__;

		/// Histogram to get the PDF
		TProfile* hFluxVsLambdaObs__;
		TH2D* hFluxPDF__;
		TH1D* hRedshift__;
		TH1D* hFlux__;

		void defineHistos(void);
		void saveHistos(void);
		void GetData(void);
		void GetQSO(void);

	public:
		GetDelta(int argc, char** argv);
		~GetDelta();
};


#endif


