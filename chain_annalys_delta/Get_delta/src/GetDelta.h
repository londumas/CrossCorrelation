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

#include <vector>
#include "TH1D.h"
#include "TGraphErrors.h"

#ifndef GETDELTA_H
#define GETDELTA_H

/*
/// If LYA
const std::string forest__  = "LYA";
const double lambdaRFLine__  = 1215.67;
const double lambdaRFMin__   = 1040.;
const double lambdaRFMax__   = 1200.;
const unsigned int nbBinRFMax__   = 647;
const double alphaStart__ = 1.;
#define C_NBBINZ 12
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5
*/
/// If LYA_JMC
//const unsigned int nbBinRFMax__ = 645;


/// If CIV
const std::string forest__   = "CIV";
const double lambdaRFLine__  = 1550.77845;
const double lambdaRFMin__   = 1410.;
const double lambdaRFMax__   = 1530.;
const unsigned int nbBinRFMax__   = 373;
const double alphaStart__ = 1.;
#define C_NBBINZ 9
#define C_ZEXTREMA0 1.4
#define C_ZEXTREMA1 3.2

/*
/// If MGII
const std::string forest__  = "MGII";
const double lambdaRFLine__  = 2804.;
const double lambdaRFMin__   = 1970.;
const double lambdaRFMax__   = 2760.;
const unsigned int nbBinRFMax__ = 1450;
const double alphaStart__ = 1.;
#define C_NBBINZ 9
#define C_ZEXTREMA0 0.4
#define C_ZEXTREMA1 1.7
*/
/*
/// If LYB
const std::string forest__  = "LYB";
const double lambdaRFLine__  = 1025.72;
const double lambdaRFMin__   = 800.;
const double lambdaRFMax__   = 1020.;
const unsigned int nbBinRFMax__ = 1085;
const double alphaStart__ = 1.;
#define C_NBBINZ 12
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5
*/
/*
/// If SIIV
const std::string forest__  = "SIIV";
const double lambdaRFLine__  = 1402.77291;
const double lambdaRFMin__   = 1286.;
const double lambdaRFMax__   = 1380.;
const unsigned int nbBinRFMax__ = 326;
const double alphaStart__ = 1.;
#define C_NBBINZ 9
#define C_ZEXTREMA0 1.4
#define C_ZEXTREMA1 3.4
*/


const double C_BINSIZEZ = (C_ZEXTREMA1-C_ZEXTREMA0)/C_NBBINZ;

/// Constants for the PDF
const double minRedshift__ = 1.96;
const double maxRedshift__ = 3.9;
const unsigned int nbBinsRedshift__ = 60;
const double minFlux__ = 0.;
const double maxFlux__ = 1.000001;
const unsigned int nbBinsFlux__ = 100;


const double lambdaObsMin__   = 3547.;  //3600.; //
const double lambdaObsMax__   = 10326.;  //7235.; //
 
/// Value allowing to put a roof to the error on a flux
//const double capForError__ = 0.19;

const double z0__        = 2.25;
const double gama__      = 3.8;
const double minPowErrorDelta__ = 1.e-10;
const double maxPowErrorDelta__ = 1.e10;
const double minPowErrorDeltaForFit__ = 1.e-10;
const double maxPowErrorDeltaForFit__ = 200.;
const double etaStart        = 1.;
const double sigma2LSSStart  = 0.1;
const double minAlpha__   = -40.;
const double maxAlpha__   = 40.;
const double betaStart__ = 0.;
const double minBeta__  = -0.3;
const double maxBeta__  =  0.3;

/// Max values of DLA
const unsigned int NbDLA0__ = 15;
const unsigned int nbLoop__ = 10;


class GetDelta
{
	private:

		std::vector< double > v_zz__;
                std::vector< double > v_meanForestLambdaRF__;
                std::vector< double > v_alpha2__;
		std::vector< double > v_beta2__;
		std::vector<unsigned int> v_nbPixel__;
		std::vector<int> v_fromFitsIndexToVectorIndex__;

		std::vector< std::vector< double > > v_LAMBDA_OBS__;
		std::vector< std::vector< double > > v_LAMBDA_RF__;
		std::vector< std::vector< double > > v_NORM_FLUX__;
		std::vector< std::vector< double > > v_NORM_FLUX_IVAR__;
		std::vector< std::vector< double > > v_FLUX_DLA__;
		std::vector< std::vector< double > > v_DELTA__;
		std::vector< std::vector< double > > v_DELTA_IVAR__;
		std::vector< std::vector< double > > v_DELTA_WEIGHT__;
		std::vector< std::vector< double > > v_TEMPLATE__;
		std::vector< std::vector< double > > v_ZZZ__;
		std::vector< std::vector< double > > v_FACTORWEIGHT__;


		TH1D* hTemplate__[nbLoop__];
		TH1D* hDeltaVsLambdaObs__[nbLoop__+1];
		TH1D* hDeltaVsLambdaRF__[nbLoop__+1];
		TGraphErrors* grEta__[nbLoop__];
		TGraphErrors* grSig__[nbLoop__];

		void loadDataForest(std::string fitsnameSpec, unsigned int start, unsigned int end, bool takeNotFittedSpectra=false);
		void getHisto(unsigned int loopIdx);
		void updateDeltaVector(unsigned int loopIdx);
		void fitForests(unsigned int begin, unsigned int end);

		void updateDelta(std::string fitsnameSpec, unsigned int loopIdx, unsigned int start, unsigned int end);
		void putReobsTogether(std::string fitsnameSpec, unsigned int loopIdx);
		void updateDLA(std::string fitsnameSpec, unsigned int start, unsigned int end);
		void updateFlux(std::string fitsnameSpec, unsigned int start, unsigned int end);

		void defineHistos();
		void FitContinuum(float* Flux, float* FluxErr, float* Lambda, float* zHI, int Nflux, float LambdaMean, double redshift, double* param);
		void initFitCont(void);
		double VoigtProfile(float nhi, float lamb, float z_abs);

		void makeCoAdd(unsigned int NbLambda, double* Flux, double* Lambda, double* ErrFlux, unsigned int NbLambdaObs, double* FluxObs, double* LambdaObs, double* ErrFluxObs);
		void sumFlux (double& f1, double& s1, double f2, double s2);


	public:
		GetDelta(int argc, char** argv);
		~GetDelta();
};


#endif


