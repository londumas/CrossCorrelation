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

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#ifndef GETDELTA_H
#define GETDELTA_H

/// Constants
const double CONVERT_FROM_FLUX_TO_ALPHA = 1.2;
const unsigned int nbPixelsTemplate__ = 647;
const double lambdaRFMin__ = 1040.;
const double lambdaRFMax__ = 1200.;
const double lambdaRFLine__  = 1215.67;
const double lambdaRFNormaMin__ = 1275.;
const double lambdaRFNormaMax__ = 1295.;
const double lambdaObsMin__ = 3600.;
const double lambdaObsMax__ = 7235.;


const double minRedshift__ = 1.96;
const double maxRedshift__ = 3.9;
const unsigned int nbBinsRedshift__ = 60;
const double minFlux__ = 0.;
const double maxFlux__ = 1.000001;
const unsigned int nbBinsFlux__ = 100;


class GetDelta
{
	private:

		double ratioForestToQSO__;

		std::string pathToDataQSO__;
		std::string pathToDataForest__;
		std::string pathToSave__;
		bool noCont__;
		bool noNoise__;
		bool isTest__;
		bool withRSD__;
		bool noMockExpander__;

		/// Histogram to get the PDF
		TProfile* hFluxVsLambdaObs__;
		TH2D* hFluxPDF__;
		TH1D* hRedshift__;
		TH1D* hFlux__;

		void defineHistos(void);
		void saveHistos(void);
		void GetData(void);
		void GetPDF(unsigned int version);
		void GetData_from_Jean_Marc_file(void);
		void GetQSO(void);

	public:
		GetDelta(int argc, char** argv);
		~GetDelta();
};


#endif


