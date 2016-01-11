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
#include "TGraphErrors.h"


#ifndef GETDELTA_H
#define GETDELTA_H





/// If LYA
const std::string forest__  = "LYA";
const double lambdaRFLine__  = 1215.67;
const double lambdaRFMin__   = 1040.;
const double lambdaRFMax__   = 1200.;
const unsigned int nbBinRFMax__   = 647;
const double alphaStart__ = 1.;
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5
#define C_BINSIZEZ 2.08333333333333343e-01

/// If LYA_JMC
//const unsigned int nbBinRFMax__ = 645;

/*
/// If CIV
const std::string forest__   = "CIV";
const double lambdaRFLine__  = 1550.77845;
const double lambdaRFMin__   = 1410.;
const double lambdaRFMax__   = 1530.;
const unsigned int nbBinRFMax__   = 373;
const double alphaStart__ = 1.;
#define C_ZEXTREMA0 1.4
#define C_ZEXTREMA1 3.2
#define C_BINSIZEZ 0.2
*/

/*
/// If MGII
const std::string forest__  = "MGII";
const double lambdaRFLine__  = 2804.;
const double lambdaRFMin__   = 1970.;
const double lambdaRFMax__   = 2760.;
const unsigned int nbBinRFMax__ = 1450;
const double alphaStart__ = 1.;
#define C_ZEXTREMA0 0.4
#define C_ZEXTREMA1 1.7
#define C_BINSIZEZ 0.15
*/

/*
/// If LYB
const std::string forest__  = "LYB";
const double lambdaRFLine__  = 1025.72;
const double lambdaRFMin__   = 800.;
const double lambdaRFMax__   = 1020.;
const unsigned int nbBinRFMax__ = 1085;
const double alphaStart__ = 1.;
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5
#define C_BINSIZEZ 2.08333333333333343e-01
*/
/*
/// If SIIV
const std::string forest__  = "SIIV";
const double lambdaRFLine__  = 1402.77291;
const double lambdaRFMin__   = 1286.;
const double lambdaRFMax__   = 1380.;
const unsigned int nbBinRFMax__ = 326;
const double alphaStart__ = 1.;
#define C_ZEXTREMA0 1.4
#define C_ZEXTREMA1 3.4
#define C_BINSIZEZ 0.2
*/




/// Constants for the PDF
const double minRedshift__ = 1.96;
const double maxRedshift__ = 3.9;
const unsigned int nbBinsRedshift__ = 60;
const double minFlux__ = 0.;
const double maxFlux__ = 1.000001;
const unsigned int nbBinsFlux__ = 100;


const double lambdaObsMin__   = 3600.; //3547.;  //3600
const double lambdaObsMax__   = 7235.; //10326;  //7235
const unsigned int nbBinRFMin__ = 50;
 
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

/// Sky lines
const unsigned int nbVetoLines__ = 21;
const double vetoLine__[42] = {3615.,3619.,3932.,3937.,3966.,3972.,4042.,4050.,
	4357.,4362.,5458.,5467.,5573.,5585.,5682.,5695.,
	5885.,5902.,6235.,6241.,6256.,6263.,6296.,6311.,
	6320.,6334.,6362.,6369.,6498.,6502.,6554.,6557.,
	6825.,6840.,6862.,6870.,6922.,6928.,6948.,6954.,
	6977.,6982.};


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


