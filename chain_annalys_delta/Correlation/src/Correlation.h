//===================================================================================
//
//         FILE: correlation.h
//
//        USAGE: ---
//
//  DESCRIPTION:
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
#include <cmath>

#ifndef CORRELATION_H
#define CORRELATION_H

class Correlation
{
	private:

		///
		/// in data, angle are in radian
		///

		std::string commandEnd__;
		std::string pathForest__;


		/// Data Forest
		unsigned int nbForest_;
		std::vector<double> v_ra__;
		std::vector<double> v_de__;
		std::vector<double> v_CosDe__;
		std::vector<double> v_SinDe__;
		std::vector<double> v_zz__;
		std::vector<unsigned int> v_nbPixelDelta1__;
		std::vector<unsigned int> v_idx__;
		std::vector< std::vector< double > > v_r__;
		std::vector< std::vector< double > > v_d__;
		std::vector< std::vector< double > > v_w__;
		std::vector< std::vector< double > > v_z__;
		std::vector< std::vector< double > > v_lRF__;
		std::vector< std::vector< double > > v_lObs__;
		std::vector< std::vector< unsigned int > > v_nb__;
		/// Data forest 2
		unsigned int nbForest2__;
		std::vector<double> v_raDelta2__;
		std::vector<double> v_deDelta2__;
		std::vector<double> v_CosDeDelta2__;
		std::vector<double> v_SinDeDelta2__;
		std::vector<double> v_zzDelta2__;
		std::vector<unsigned int> v_nbPixelDelta2__;
		std::vector< std::vector< double > > v_rDelta2__;
		std::vector< std::vector< double > > v_dDelta2__;
		std::vector< std::vector< double > > v_wDelta2__;
		std::vector< std::vector< double > > v_zDelta2__;
		std::vector< std::vector< double > > v_lRFDelta2__;
		std::vector< std::vector< double > > v_lObsDelta2__;
		/// Data QSO
		unsigned int nbQ1__;
		std::vector<double> v_raQ1__;
		std::vector<double> v_deQ1__;
		std::vector<double> v_zzQ1__;
		std::vector<double> v_rQ1__;
		std::vector<double> v_CosDeQ1__;
		std::vector<double> v_SinDeQ1__;
		/// Data Q2
		unsigned int nbQ2__;
		std::vector<double> v_raQ2__;
		std::vector<double> v_deQ2__;
		std::vector<double> v_zzQ2__;
		std::vector<double> v_rQ2__;
		std::vector<double> v_CosDeQ2__;
		std::vector<double> v_SinDeQ2__;

		/// Instance
		void loadDataQ1(void);
		void loadDataQ2(void);
		void loadDataForest(std::string,bool doBootstraps=false, unsigned int bootIdx=0);
		void loadDataDelta2(int dataNeeded=100);

		void removeFalseCorrelations(bool firstPass=true);

	public:
		Correlation(int argc, char** argv);
		~Correlation();

		//
		void xi_1D_delta_delta(void);
		void xi_1DlRF_delta_delta(void);
		void xi_1DlRFDevide_delta_delta(void);
		void xi_1D_delta_delta2(void);
		void xi_1DlRFDevide_delta_delta2(void);
		//
		void xi_A_delta_delta(void);
		void xi_A_delta_delta_lambda(void);
		void xi_A_delta_delta2( unsigned int bootIdx=0 );
		void xi_A_delta_delta2_lambda(void);
		//
		void xi_delta_QSO(         bool doBootstraps=false, unsigned int bootIdx=0);
		void xi_delta_QSO_theta(   bool doBootstraps=false, unsigned int bootIdx=0);
		void xi_delta_QSO_lambda(   bool doBootstraps=false, unsigned int bootIdx=0);
		void xi_delta_QSO_distortionMatrix(void);
		void xi_delta_QSO_distortionMatrix_1D(void);
		void xi_delta_QSO_Metals_Models(double lambdaFrMetal, std::string lambdaFrMetalName);
		void xi_delta_QSO_Wick(unsigned int diagramIdx);
		//
		void xi_QSO_QSO(bool doBootstraps=false, unsigned int bootIdx=0);
		void xi_Q1_Q2(void);
		//
		void xi_A_delta_delta_MockJMc(void);
		void xi_delta_QSO_MockJMc( bool doBootstraps=false, unsigned int bootIdx=0);
		void xi_delta_QSO_lambda_MockJMc( bool doBootstraps=false, unsigned int bootIdx=0);
		void xi_delta_QSO_MockJMc_distortionMatrix( void );
		void xi_delta_QSO_MockJMc_distortionMatrix_1D( void );
		void xi_delta_QSO_Metals_Models_MockJMc(double lambdaFrMetal, std::string lambdaFrMetalName);
		void xi_QSO_QSO_MockJMc(bool doBootstraps=false, unsigned int bootIdx=0);
};

const double lambdaObsMin__   = 3600.; //
const double lambdaObsMax__   = 7235.; //10326.;  //7235.; //
const double maxAlpha__   = 40.;
const double maxBeta__  =  0.3;
const double maxCorrelation__ = 200.;
const double isReobsFlag__ = -100.;
const unsigned int sizeGridX__   = 2560;
const unsigned int sizeGridY__   = 1920;
const unsigned int sizeGridZ__   = 512;
const double sizeCell__         = 4.5*0.71;  //3.195; // in Mpc.h^-1 ( = 4.5*0.71 ) ///0.3195




/*
/// If LYB
const std::string forest__      = "LYB";
const double lambdaRFLine__     = 1025.72;
const double lambdaRFMin__      = 800.;
const double lambdaRFMax__      = 1020.;
const unsigned int nbBinRFMax__ = 1085;
const double alphaStart__       = 1.;
*/

/// If LYA
const std::string forest__      = "LYA";
const double lambdaRFLine__     = 1215.67;
const double lambdaRFMin__      = 1040.;
const double lambdaRFMax__      = 1200.;
const unsigned int nbBinRFMax__ = 647;
const double alphaStart__       = 1.3;
/// If LYA_JMC
//const unsigned int nbBinRFMax__ = 645;

/*
/// If SIIV
const std::string forest__      = "SIIV";
const double lambdaRFLine__     = 1393.76018;
const double lambdaRFMin__      = 1286.;
const double lambdaRFMax__      = 1380.;
const unsigned int nbBinRFMax__ = 326;
const double alphaStart__       = 1.;
*/
/*
/// If CIV
const std::string forest__      = "CIV";
const double lambdaRFLine__     = 1548.2049;
const double lambdaRFMin__      = 1410.;
const double lambdaRFMax__      = 1530.;
const unsigned int nbBinRFMax__ = 373;
const double alphaStart__       = 1.;
*/
/*
/// If MGII
const std::string forest__      = "MGII";
const double lambdaRFLine__  = 2803.5324;
const double lambdaRFMin__   = 1570.;
const double lambdaRFMax__   = 2790.;
const unsigned int nbBinRFMax__ = 2414;
const double alphaStart__       = 1.;
*/







/*
/// For Delta 2 (LYB)
const std::string forest2__  = "LYB";
const double lambdaRFLineDelta2__  = 1025.72;
const double lambdaRFMinDelta2__   = 800.;
const double lambdaRFMaxDelta2__   = 1020.;
const unsigned int nbBinRFMaxDelta2__ = 1080;
const std::string pathDelta2__ = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYB/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits";
const double alphaStart2__       = 1.;
*/
/*
/// For Delta 2 (LYA)
const std::string forest2__      = "LYA";
const double lambdaRFLineDelta2__     = 1215.67;
const double lambdaRFMinDelta2__      = 1040.;
const double lambdaRFMaxDelta2__      = 1200.;
const unsigned int nbBinRFMaxDelta2__ = 647;
const std::string pathDelta2__ = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits";
const double alphaStart2__       = 1.3;
*/
/*
/// For Delta 2 (SIIV)
const std::string forest2__           = "SIIV";
const double lambdaRFLineDelta2__     = 1393.76018;
const double lambdaRFMinDelta2__      = 1286.;
const double lambdaRFMaxDelta2__      = 1380.;
const unsigned int nbBinRFMaxDelta2__ = 326;
const std::string pathDelta2__ = "/home/gpfs/manip/mnt/bao/hdumasde/Data/SIIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits";
const double alphaStart2__       = 1.;
*/

/// For Delta 2 (CIV)
const std::string forest2__  = "CIV";
const double lambdaRFLineDelta2__  = 1548.2049;
const double lambdaRFMinDelta2__   = 1410.;
const double lambdaRFMaxDelta2__   = 1530.;
const unsigned int nbBinRFMaxDelta2__ = 373;
const std::string pathDelta2__ = "/home/gpfs/manip/mnt/bao/hdumasde/Data/CIV/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits";
const double alphaStart2__       = 1.;

/*
/// For Delta 2 (MGII)
const std::string forest2__  = "MGII";
const double lambdaRFLineDelta2__  = 2804;
const double lambdaRFMinDelta2__   = 1570.;
const double lambdaRFMaxDelta2__   = 2790.;
const unsigned int nbBinRFMaxDelta2__ = 2414;
const std::string pathDelta2__ = "/home/gpfs/manip/mnt/bao/hdumasde/Data/MGII/FitsFile_DR12/DR12_primery/DR12_primery.fits";
const double alphaStart2__       = 1.;
*/



/// For Q2
const std::string QSO2__   = "aaa";
const std::string pathQ2__ = "";





#endif


