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

///// LIB
#include <iostream>	// std::cout
#include <fstream>
#include <sstream>	//stringstream
#include <cmath>
#include <assert.h>
#include <vector>
#include <algorithm>	// std::random_shuffle
#include <cstdlib>	// std::rand, std::srand

///// ROOT
#include "TH1D.h"
#include "TFile.h"
#include "fitsio.h"
#include "TString.h"
//#include "TRandom.h"

///// Mine
#include "Correlation.h"
#include "LymanAlphaForestRegion.h"
#include "Cosmology.h"
#include "../../../Root/Library/RootHistoFunctions.h"
#include "../../../Cpp/Library/mathFunctions.h"
#include "../../../Constants/constants.h"
#include "../../../Constants/globalValues.h"

std::string QSO__    = "";
std::string pathQ1__ = "";
std::string alQSOName[21] = {"QSO",
				"DLA",
				"CMASS",
				"CMASSLOWZ",
				"LOWZ",
				"ALL_BRITT",
				"CIV_BRITT",
				"MGII_BRITT",
				"OTHER_BRITT",
				"CIVNOTQSO_BRITT",
				"ALLNOTQSO_BRITT",
				"CIVLOWEW_BRITT",
				"CIVHIGHEW_BRITT",
				"VIPERS",
				"QSO_3DHST",
				"QSO_EBOSS",
				"QSO_DR7_DR12_EBOSS",
				"ALL_OBJECTS",
				"ALL_EVERY_OBJECTS",
				"QSO_DR7_DR12_V5_10_0",
				"QSO_DR14_v1_0"
};
std::string alQSO[21] = {"QSO_ALL_TESTS",
				"DLA_all",
				"CMASS_all",
				"CMASSLOWZ_all",
				"LOWZ_all",
				"all_Britt",
				"CIV_Britt",
				"MGII_Britt",
				"OTHER_Britt",
				"CIVNOTQSO_Britt",
				"ALLNOTQSO_Britt",
				"CIVLOWEW_Britt",
				"CIVHIGHEW_Britt",
				"VIPERS",
				"QSO_3DHST",
				"QSO_EBOSS_updated_2016_05_24",
				"QSO_DR7_DR12_EBOSS_2016_05_24",
				"ALL_OBJECTS",
				"ALL_EVERY_OBJECTS_2016_01_08",
				"QSO_DR7_DR12_V5_10_0",
				"QSO_DR14_v1_0"
};

///// Constants
const unsigned int nbBinlambdaObs__  = int(lambdaObsMax__-lambdaObsMin__);
const double onePlusZ0__             = 1.+z0__;
const double halfGama__              = gama__/2.;
double distMinPixel__                = 0.;
double distMinPixelDelta2__          = 0.;
double distMinQ1__                   = 0.;
double distMinQ2__                   = 0.;
unsigned int idxCommand_[6]          = {0};
const std::string pathToMockJMC__    = "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/";
std::string pathToRaw__              = "";
std::string pathToSave__             = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_Margala__projected_PlanckCosmo_correctedDistance/";
std::string correlation_type__       = "NOTHING";
std::string path_to_rand_QSO_cat__   = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR14_v1_0_RAND/";


///// Type of data
const bool mocks              = false;
const bool mockJMC__          = false;
const bool mockBox__          = false;
//// Attributes of data
 // versin=0 for old version, version=1 for mockExpander, version=2 for files from Jean-Marc
const unsigned int mock_version = 5;
const bool mocks_raw          = false;
const bool meanOver3pixels__  = false;
const bool mocksNoNoiseNoCont = false;
const bool noMetals__         = false;
const bool withRSD__          = false;
const bool randomPositionOfQSOInCellNotBeforeCorrelation__ = false;
unsigned int seed_for_random_position__ = 42;
//// Flags for covariance matrix estimation
const bool shuffleQSO            = false;
const bool shuffleForest         = false;
const bool randomQSO             = false;
const bool randomQSOXYZ          = false;
const bool randomQSO_from_file__ = false;
const bool randomForest          = false;
const bool doBootstraps__        = false;


const bool doVetoLines__          = true;
const bool nicolasEstimator__     = true;
const bool doingCoAddSoRemoving__ = false;

const bool replace_data_by_gaussian_random__ = false;
const bool showRootHistos__      = false;
const bool saveInRootFile__      = false;
const bool cutNotFittedSpectra__ = true;
const bool selectWindowRedshift_bool__ = true;

Correlation::Correlation(int argc, char **argv) {

	std::cout << std::scientific;
	std::cout.precision(std::numeric_limits<double>::digits10);

	// command

	///// 0: correlation, 1: mockCatalogue, 2: bootstraps, 3: Wick, 4: mockChunck, 5: mockSimul
	for (unsigned int i=0; i<6; i++) {
		std::string string = argv[i+1];
		idxCommand_[i] = atoi(string.c_str());
	}

	//// seed for random position
	std::srand(42);
	for (unsigned int i=0; i<1000; i++) {
		const unsigned int seed = rand();
		if (i == (idxCommand_[4]*10 + idxCommand_[5]) ) seed_for_random_position__ = seed;
	}

	
	if (idxCommand_[1]<21) {
		QSO__    = alQSOName[idxCommand_[1]];
		pathQ1__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/";
		pathQ1__ += alQSO[idxCommand_[1]];
		pathQ1__ += ".fits";
	}

	///// Set the number of forest to work on
	nbForest_   = 0;
	nbForest2__ = 0;
	nbQ1__      = 0;
	nbQ2__      = 0;
	if (!mocks && !mockJMC__) {
		pathForest__  = "/home/gpfs/manip/mnt/bao/hdumasde/Data/";
		pathForest__  += forestPath__;
		pathForest__  += "/FitsFile_DR12_Guy_Margala/DR12_primery/DR12_primery.fits";
	}
	else if (mocks) {
		pathForest__   = "/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_";
		pathForest__  += argv[5];
		pathForest__  += "/00";
		pathForest__  += argv[6];
		pathForest__  += "/mock.fits";

		pathToSave__  = "/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_";
		pathToSave__ += argv[5];
		pathToSave__ += "/00";
		pathToSave__ += argv[6];
		pathToSave__ += "/Results/";

		if (QSO__ == "DLA") {
			pathQ1__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/DLA_Catalogue/new_M3_0_";
			pathQ1__ += argv[5];
			pathQ1__ += "__";
			pathQ1__ += argv[6];
			pathQ1__ += ".fits";
		}
	}
	else if (mockJMC__) {

		QSO__     = "QSO";

		pathForest__ = pathToMockJMC__;
		pathForest__ += "Box_00";
		pathForest__ += argv[5];
		pathForest__ += "/Simu_00";
		pathForest__ += argv[6];
		pathForest__ += "/";
		pathToSave__  = pathForest__;
		pathToRaw__   = pathForest__;

		pathQ1__ = pathForest__;
		if (withRSD__) pathQ1__ += "Data/QSO_withRSD.fits";
		else pathQ1__ += "Data/QSO_noRSD.fits";

		if (noMetals__) pathForest__ += "Data_no_metals/delta.fits";
		else            pathForest__ += "Data/delta.fits";

		pathToSave__ += "Results";
		if (mocks_raw) {
			if (mock_version==1) pathToSave__ += "_Raw/";
			if (mock_version==2) pathToSave__ += "_raw_from_JeanMarc/";
			if (mock_version==3) pathToSave__ += "_Flux_Notemplate_WithMetals_WithNoise_WithSpectroResolution/";
			if (mock_version==4) pathToSave__ += "_delta_gaussian/";
			if (mock_version==5) {
				if (randomPositionOfQSOInCellNotBeforeCorrelation__) pathToSave__ += "_only_LR/";
				else {
					if (withRSD__) pathToSave__ += "_only_LR_noRand/";
					else pathToSave__ += "_only_LR_noRand_noRSD/";
				}
			}
		}
		else {
			if (!nicolasEstimator__) pathToSave__ += "_no_projection/";
			else if (noMetals__)     pathToSave__ += "_no_metals/";
			else                     pathToSave__ += "/";
		}
		

		pathToRaw__ += "Raw/mocks-*";

		if (mock_version==0) {
			pathToRaw__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-78";
			pathToRaw__ += argv[5];
			pathToRaw__ += "-";
			pathToRaw__ += argv[6];
			pathToRaw__ += ".fits";
		}
		else if (mock_version==2) {
			pathToRaw__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1575/fits/spectra-785";
			pathToRaw__ += argv[5];
			pathToRaw__ += "-";
			pathToRaw__ += argv[6];
			pathToRaw__ += ".fits";
		}
		else if (mock_version==4) {
			pathToRaw__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1575delta/fits/spectra-785";
			pathToRaw__ += argv[5];
			pathToRaw__ += "-";
			pathToRaw__ += argv[6];
			pathToRaw__ += ".fits";
		}
		else if (mock_version==5) {
			pathToRaw__ = " /home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1652/fits/spectra-785";
			pathToRaw__ += argv[5];
			pathToRaw__ += "-";
			pathToRaw__ += argv[6];
			pathToRaw__ += ".fits";
		}
	}

	///// Command at the end
	commandEnd__  = " ";
	commandEnd__ += forest__;
	commandEnd__ += " ";
	commandEnd__ += forest2__;
	commandEnd__ += " ";
	commandEnd__ += QSO__;
	commandEnd__ += " ";
	commandEnd__ += QSO2__;
	std::cout << "  " << commandEnd__ << std::endl;
	
	const unsigned int command = idxCommand_[0];


	/// Set the type of correlation
	if (command/100==0)      correlation_type__ = "f_f_r";
	else if (command/100==1) correlation_type__ = "f_f2_r";
	else if (command/100==2) correlation_type__ = "f_f";
	else if (command/100==3) correlation_type__ = "f_f2";
	else if (command/100==4) correlation_type__ = "q_f";
	else if (command/100==5) correlation_type__ = "q_q";
	else if (command/100==6) correlation_type__ = "q_q2";

	///
	if      (command == 0) xi_1D_delta_delta();
	else if (command == 1) xi_1DlRF_delta_delta();
	else if (command == 2) xi_1DlRFDevide_delta_delta();
	else if (command == 3) xi_1DlObs_2D_delta_delta();
	else if (command == 4) xi_1D_delta_delta_distortionMatrix();
	else if (command == 5) xi_1DlRF_delta_delta_distortionMatrix();
	else if (command == 6) xi_1DlRFDevide_delta_delta_distortionMatrix();
	///
	else if (command == 100) xi_1D_delta_delta2();
	else if (command == 101) xi_1DlRFDevide_delta_delta2();
	///
	else if (command == 200) xi_A_delta_delta(idxCommand_[2]);
	else if (command == 201) xi_A_delta_delta_lambda();
	else if (command == 202) xi_A_delta_delta_distortionMatrix();
	else if (command == 203) xi_A_delta_delta_Metals_Models(absorber__[idxCommand_[1]],absorberName__[idxCommand_[1]],absorber__[idxCommand_[2]],absorberName__[idxCommand_[2]]);
	else if (command == 206) xi_A_delta_delta_MockJMc(idxCommand_[2]);
	else if (command == 207) xi_A_delta_delta_MockJMc_distortionMatrix();
	else if (command == 208) xi_A_delta_delta_Metals_Models_MockJMc(absorber__[idxCommand_[1]],absorberName__[idxCommand_[1]],absorber__[idxCommand_[2]],absorberName__[idxCommand_[2]]);
	///
	else if (command == 300) xi_A_delta_delta2( idxCommand_[2] );
	else if (command == 301) xi_A_delta_delta2_lambda();
	///
	else if (command == 400) xi_delta_QSO(idxCommand_[2]);
	else if (command == 401) xi_delta_QSO_lambda(idxCommand_[2]);
	else if (command == 402) xi_delta_QSO_lambda_same_LOS(idxCommand_[2]);
	else if (command == 403) xi_delta_QSO_distortionMatrix();
	else if (command == 404) xi_delta_QSO_distortionMatrix_1D();
	else if (command == 405) xi_delta_QSO_Metals_Models(absorber__[idxCommand_[2]],absorberName__[idxCommand_[2]]);
	else if (command == 406) xi_delta_QSO_Wick_T1_with_wi1D_array();
	else if (command == 407) xi_delta_QSO_Wick_T12_with_wi1D_array();
	else if (command == 408) xi_delta_QSO_Wick_T123_with_wi1D_array();
	else if (command == 409) xi_delta_QSO_Wick_T1234_with_wi1D_array();
	else if (command == 410) xi_delta_QSO_Wick_1D(idxCommand_[3]);

	else if (command == 411) xi_delta_QSO_MockJMc(idxCommand_[2]);
	else if (command == 412) xi_delta_QSO_MockJMc_distortionMatrix();
	else if (command == 413) xi_delta_QSO_MockJMc_distortionMatrix_1D();
	else if (command == 414) xi_delta_QSO_Metals_Models_MockJMc(absorber__[idxCommand_[2]],absorberName__[idxCommand_[2]]);
	else if (command == 415) xi_delta_QSO_MockJMc_Wick_T1_with_wi1D_array();
	else if (command == 416) xi_delta_QSO_MockJMc_Wick_T12_with_wi1D_array();
	else if (command == 417) xi_delta_QSO_MockJMc_Wick_T123_with_wi1D_array();
	else if (command == 418) xi_delta_QSO_MockJMc_Wick_T1234_with_wi1D_array();
	///
	else if (command == 500) xi_QSO_QSO(idxCommand_[2]);
	else if (command == 502) xi_QSO_QSO_MockJMc(idxCommand_[2]);
	///
	else if (command == 601) xi_Q1_Q2();


	std::cout << "\n\n\n\n" << std::endl;
}
Correlation::~Correlation(void) {	
}





// ---------------------------------------------------------------------
//
//		1D, same forest - same forest correlation
//
// ---------------------------------------------------------------------
void Correlation::xi_1D_delta_delta(void) {

	std::cout << "\n\n\n\n  ------ xi_1D_delta_delta ------\n" << std::endl;

	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	if (mocksNoNoiseNoCont || mocks_raw) {
		removeFalseCorrelations();
	}
	v_ra__.clear();
	v_de__.clear();
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Constants:
	///// The space between bins is of 1 Mpc.h^-1, the first bin is for deltaR == 0.
	const double max = 2100.;
	const unsigned int nbBins = 2100;
	const double binSize = max/nbBins;
	const double inv_binSize = 1./binSize;

	///// Arrays for data
	double data[nbBins][6];
	for (unsigned int i=0; i<nbBins; i++) {
		for (unsigned int j=0; j<6; j++) {
			data[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;


	for (unsigned int f=0; f<nbForest_; f++) {

		const unsigned int nb = v_nbPixelDelta1__[f];

		for (unsigned int i1=0; i1<nb; i1++) {

			const double r1 = v_r__[f][i1];
			const double d1 = v_d__[f][i1];
			const double w1 = v_w__[f][i1];
			const double z1 = v_z__[f][i1];

			const double w1w1 = w1*w1;
			const double d1d1 = d1*d1;
			const double z1z1 = z1+z1;

			///// Compute the xi1D
			data[0][0] += w1w1*d1d1;
			data[0][1] += w1w1*d1d1*d1d1;
			data[0][3] += w1w1*z1z1;
			data[0][4] += w1w1;
			data[0][5] ++;

			for (unsigned int i2=0; i2<i1; i2++) {

				const double r1r2 = r1-v_r__[f][i2];
				const double w1w2 = w1*v_w__[f][i2];
				const double d1d2 = d1*v_d__[f][i2];
				const double z1z2 = z1+v_z__[f][i2];

				///// Find bin index
				const unsigned int idx = 1 + int(r1r2*inv_binSize);

				///// Compute the xi1D
				data[idx][0] += w1w2*d1d2;
				data[idx][1] += w1w2*d1d2*d1d2;
				data[idx][2] += w1w2*r1r2;
				data[idx][3] += w1w2*z1z2;
				data[idx][4] += w1w2;
				data[idx][5] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_1D_delta_delta_";
	pathToSave += forest__;
	pathToSave += ".txt";
	std::cout << "\n  pathToSave = " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBins; i++) {
		fFile << data[i][0];
		fFile << " " << data[i][1];
		fFile << " " << data[i][2];
		fFile << " " << data[i][3]/2.;
		fFile << " " << data[i][4];
		fFile << " " << data[i][5];
		fFile << std::endl;

		sumOne += data[i][5];
		sumZZZ += data[i][3]/2.;
		sumWeg += data[i][4];
	}
	fFile.close();

	///// Get mean redshift of pairs
	std::cout << "  < z >           = " << sumZZZ/sumWeg << std::endl;
	std::cout << "  Number of pairs = " << sumOne << std::endl;
	
	return;
}
void Correlation::xi_1DlRF_delta_delta(void) {

	std::cout << "\n\n\n\n  ------ xi_1DlRF_delta_delta ------\n" << std::endl;

	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	if (mocksNoNoiseNoCont || mocks_raw) {
		removeFalseCorrelations();
	}
	v_ra__.clear();
	v_de__.clear();
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_r__.clear();
	v_lObs__.clear();
	v_nb__.clear();


	///// Constants:
	const double max           = lambdaRFMax__-lambdaRFMin__;
	const unsigned int nbBins  = nbBinRFMax__-1;
	const double binSize = max/nbBins;
	const double inv_binSize = 1./binSize;

	///// Arrays for data
	double data[nbBins][6];
	for (unsigned int i=0; i<nbBins; i++) {
		for (unsigned int j=0; j<6; j++) {
			data[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;
	
	for (unsigned int f=0; f<nbForest_; f++) {

		const unsigned int nb = v_nbPixelDelta1__[f];

		for (unsigned int i1=0; i1<nb; i1++) {

			const double l1 = v_lRF__[f][i1];
			const double d1 = v_d__[f][i1];
			const double w1 = v_w__[f][i1];
			const double z1 = v_z__[f][i1];

			const double w1w1 = w1*w1;
			const double d1d1 = d1*d1;
			const double z1z1 = z1+z1;

			///// Compute the xi1D
			data[0][0] += w1w1*d1d1;
			data[0][1] += w1w1*d1d1*d1d1;
			data[0][3] += w1w1*z1z1;
			data[0][4] += w1w1;
			data[0][5] ++;

			for (unsigned int i2=0; i2<i1; i2++) {

				const double l1l2  = l1-v_lRF__[f][i2];
				const double w1w2  = w1*v_w__[f][i2];
				const double d1d2  = d1*v_d__[f][i2];
				const double z1z2  = z1+v_z__[f][i2];

				///// Find bin index
				const unsigned int idx = 1 + int(l1l2*inv_binSize);

				///// Compute the xi1D
				data[idx][0] += w1w2*d1d2;
				data[idx][1] += w1w2*d1d2*d1d2;
				data[idx][2] += w1w2*l1l2;
				data[idx][3] += w1w2*z1z2;
				data[idx][4] += w1w2;
				data[idx][5] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_1DlRF_delta_delta_";
	pathToSave += forest__;
	pathToSave += ".txt";
	std::cout << "\n  pathToSave = " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	///// Set the values of data
	for (unsigned int i=0; i<nbBins; i++) {
		fFile << data[i][0];
		fFile << " " << data[i][1];
		fFile << " " << data[i][2];
		fFile << " " << data[i][3]/2.;
		fFile << " " << data[i][4];
		fFile << " " << data[i][5];
		fFile << std::endl;

		sumOne += data[i][5];
		sumZZZ += data[i][3]/2.;
		sumWeg += data[i][4];
	}
	fFile.close();

	std::cout << "  < z >           = " << sumZZZ/sumWeg << std::endl;
	std::cout << "  Number of pairs = " << sumOne << std::endl;

	return;
}
void Correlation::xi_1DlRFDevide_delta_delta(void) {

	std::cout << "\n\n\n\n  ------ xi_1DlRFDevide_delta_delta ------\n" << std::endl;

	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	if (mocksNoNoiseNoCont || mocks_raw) {
		removeFalseCorrelations();
	}
	v_ra__.clear();
	v_de__.clear();
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_r__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Constants:
	//const double min           = lambdaRFMin__/lambdaRFMax__;
	const double max           = lambdaRFMax__/lambdaRFMin__;
	const unsigned int nbBins  = 500;
	const double binSize = (max-1.)/nbBins;
	const double inv_binSize = 1./binSize;

	///// Arrays for data
	double data[nbBins][6];
	for (unsigned int i=0; i<nbBins; i++) {
		for (unsigned int j=0; j<6; j++) {
			data[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;
	
	for (unsigned int f=0; f<nbForest_; f++) {

		const unsigned int nb = v_nbPixelDelta1__[f];

		for (unsigned int i1=0; i1<nb; i1++) {

			const double l1 = v_lRF__[f][i1];
			const double d1 = v_d__[f][i1];
			const double w1 = v_w__[f][i1];
			const double z1 = v_z__[f][i1];

			const double w1w1 = w1*w1;
			const double d1d1 = d1*d1;
			const double z1z1 = z1+z1;

			///// Compute the xi1D
			data[0][0] += w1w1*d1d1;
			data[0][1] += w1w1*d1d1*d1d1;
			data[0][2] += w1w1;
			data[0][3] += w1w1*z1z1;
			data[0][4] += w1w1;
			data[0][5] ++;

			for (unsigned int i2=0; i2<i1; i2++) {

				const double l1l2  = l1/v_lRF__[f][i2];
				const double w1w2  = w1*v_w__[f][i2];
				const double d1d2  = d1*v_d__[f][i2];
				const double z1z2  = z1+v_z__[f][i2];

				///// Find bin index
				const unsigned int idx = 1 + int( (l1l2-1.)*inv_binSize );

				///// Compute the xi1D
				data[idx][0] += w1w2*d1d2;
				data[idx][1] += w1w2*d1d2*d1d2;
				data[idx][2] += w1w2*l1l2;
				data[idx][3] += w1w2*z1z2;
				data[idx][4] += w1w2;
				data[idx][5] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_1DlRFDevide_delta_delta_";
	pathToSave += forest__;
	pathToSave += ".txt";
	std::cout << "\n  pathToSave = " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	///// Set the values of data
	for (unsigned int i=0; i<nbBins; i++) {
		fFile << data[i][0];
		fFile << " " << data[i][1];
		fFile << " " << data[i][2];
		fFile << " " << data[i][3]/2.;
		fFile << " " << data[i][4];
		fFile << " " << data[i][5];
		fFile << std::endl;

		sumOne += data[i][5];
		sumZZZ += data[i][3]/2.;
		sumWeg += data[i][4];
	}
	fFile.close();

	std::cout << "  < z >           = " << sumZZZ/sumWeg << std::endl;
	std::cout << "  Number of pairs = " << sumOne << std::endl;

	return;
}
void Correlation::xi_1DlObs_2D_delta_delta(void) {

	std::cout << "\n\n\n\n  ------ xi_1DlObs_2D_delta_delta ------\n" << std::endl;

	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	if (mocksNoNoiseNoCont || mocks_raw) {
		removeFalseCorrelations();
	}
	v_ra__.clear();
	v_de__.clear();
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_r__.clear();
	v_lRF__.clear();
	v_nb__.clear();


	///// Constants:
	const unsigned int nbBins  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double*** data = new double**[nbBins];
	for (unsigned int i=0; i<nbBins; i++) {
		data[i] = new double*[nbBins];
		for (unsigned int j=0; j<nbBins; j++) {
			data[i][j] = new double[2];
			for (unsigned int k=0; k<2; k++) {
				data[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;
	
	for (unsigned int f=0; f<nbForest_; f++) {

		const unsigned int nb = v_nbPixelDelta1__[f];

		for (unsigned int i1=0; i1<nb; i1++) {

			const double d1 = v_d__[f][i1];
			const double w1 = v_w__[f][i1];
			const unsigned int bin_index_1 = 1+int(v_lObs__[f][i1]-lambdaObsMin__);

			const double w1w1 = w1*w1;
			const double d1d1 = d1*d1;

			///// Compute the xi1D
			data[bin_index_1][bin_index_1][0] += w1w1*d1d1;
			data[bin_index_1][bin_index_1][1] += w1w1;

			for (unsigned int i2=0; i2<i1; i2++) {

				const double w1w2  = w1*v_w__[f][i2];
				const double d1d2  = d1*v_d__[f][i2];
				const unsigned int bin_index_2 = int(v_lObs__[f][i2]-lambdaObsMin__);

				///// Compute the xi1D
				data[bin_index_1][bin_index_2][0] += w1w2*d1d2;
				data[bin_index_1][bin_index_2][1] += w1w2;
			}
		}
	}

	///// Set the values of upper diagonal
	for (unsigned int i=0; i<nbBins; i++) {
		for (unsigned int j=0; j<i; j++) {
			data[j][i][0] = data[i][j][0];
			data[j][i][1] = data[i][j][1];
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_1DlObs_2D_delta_delta_";
	pathToSave += forest__;
	pathToSave += ".txt";
	std::cout << "\n  pathToSave = " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	///// Set the values of data
	for (unsigned int i=0; i<nbBins; i++) {
		for (unsigned int j=0; j<nbBins; j++) {
			if (data[i][j][1]>0.) data[i][j][0] /= data[i][j][1];
			fFile << data[i][j][0] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	return;
}
void Correlation::xi_1D_delta_delta_distortionMatrix(void) {

	std::cout << "\n\n\n\n  ------ xi_1D_delta_delta_distortionMatrix ------\n" << std::endl;

	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	if (mocksNoNoiseNoCont || mocks_raw) {
		removeFalseCorrelations();
	}
	v_ra__.clear();
	v_de__.clear();
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Set usefull vectors
	std::vector<double> v_invSumWeight(nbForest_,0.);
	std::vector<double> v_meanLambda(nbForest_,0.);
	std::vector<double> v_invVarLambda(nbForest_,0.);

	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nbPixel = v_nbPixelDelta1__[f];

		long double sumWeight  = 0.;
		long double meanLambda = 0.;
		long double varLambda  = 0.;

		///// Loops over all pixels of the forest
		for (unsigned int i=0; i<nbPixel; i++) {

			const long double w = v_w__[f][i];
			const long double l = v_lRF__[f][i];

			sumWeight  += w;
			meanLambda += w*l;
			varLambda  += w*l*l;
		}

		if (sumWeight<=0.) std::cout << "   ERROR:: sumWeight <= 0. " << std::endl;
		meanLambda        /= sumWeight;
		varLambda          = varLambda - meanLambda*meanLambda*sumWeight;
		if (varLambda<=0.) std::cout << "   ERROR:: varLambda <= 0. " << std::endl;

		v_invSumWeight[f] = 1./sumWeight;
		v_meanLambda[f]   = meanLambda;
		v_invVarLambda[f] = 1./varLambda;
	}

	///// Constants:
	///// The space between bins is of 1 Mpc.h^-1, the first bin is for deltaR == 0.
	const double max = 2100.;
	const unsigned int nbBins = 2100;
	const double binSize = max/nbBins;
	const double inv_binSize = 1./binSize;

	///// Arrays for distortion matrix
	double weight[nbBins];
	double** dataMatrix = new double*[nbBins];
	for (unsigned int i=0; i<nbBins; i++) {
		weight[i] = 0.;
		dataMatrix[i] = new double[nbBins];
		for (unsigned int j=0; j<nbBins; j++) {
			dataMatrix[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;


	for (unsigned int f=0; f<nbForest_; f++) {

		const unsigned int nb     = v_nbPixelDelta1__[f];
		const double invSumWeight = v_invSumWeight[f];
		const double meanLambda   = v_meanLambda[f];
		const double invVarLambda = v_invVarLambda[f];

		std::vector< std::vector< double > > list_for_pixel_1_0(nb, std::vector<double>(2100,0.));
		std::vector< std::vector< double > > list_for_pixel_2_0(nb, std::vector<double>(2100,0.));
                std::vector< std::vector< double > > list_for_pixel_1_1(nb, std::vector<double>(2100,0.));
                std::vector< std::vector< double > > list_for_pixel_2_1(nb, std::vector<double>(2100,0.));
		std::vector< std::vector< bool > > binTouched_1(nb, std::vector<bool>(2100,false));
		std::vector< std::vector< bool > > binTouched_2(nb, std::vector<bool>(2100,false));
		double value2[2100]   = {0.};
		double value5[2100]   = {0.};
		double value6[2100]   = {0.};
		double value7[2100]   = {0.};
		bool binTouched[2100] = {false};

		///// Index of the bin
		std::vector< unsigned int > binIdx;
		std::vector< unsigned int > fromBinsToPixels1;
		std::vector< unsigned int > fromBinsToPixels2;

		for (unsigned int i1=0; i1<nb; i1++) {

			const double r1 = v_r__[f][i1];
			const double w1 = v_w__[f][i1];
			const double lambda_minus_mean_over_var1 = (v_lRF__[f][i1] -meanLambda)*invVarLambda;

			for (unsigned int i2=0; i2<nb; i2++) {

				///// Find bin index
				const double r1r2 = fabs( r1-v_r__[f][i2] );
				unsigned int idx = 0;
				if (i1!=i2)  idx = 1 + int(r1r2*inv_binSize);

				const double w2   = v_w__[f][i2];
				const double lambda_minus_mean_over_var2 = (v_lRF__[f][i2] -meanLambda)*invVarLambda;

				//// Fill arrays for the distortion
				list_for_pixel_1_0[i1][idx] += w2*invSumWeight;
                                list_for_pixel_2_0[i2][idx] += w1*invSumWeight;
                                list_for_pixel_1_1[i1][idx] += w2*lambda_minus_mean_over_var2;
                                list_for_pixel_2_1[i2][idx] += w1*lambda_minus_mean_over_var1;
				binTouched_1[i1][idx] = true;
				binTouched_2[i2][idx] = true;
				value2[idx] += w1*w2*invSumWeight*invSumWeight;
				value5[idx] += w2*invSumWeight*w1*lambda_minus_mean_over_var1;
				value6[idx] += w1*invSumWeight*w2*lambda_minus_mean_over_var2;
				value7[idx] += w1*w2*lambda_minus_mean_over_var1*lambda_minus_mean_over_var2;
				binTouched[idx] = true;

				if (i2<i1) {
					binIdx.push_back(idx);
					fromBinsToPixels1.push_back(i1);
					fromBinsToPixels2.push_back(i2);
				}

				///// Fill array of weights
				weight[idx] += w1*w2;
			}
		}

		///// Number of pairs of two pixels ( here = nb*nb)
		const unsigned int nbPixelsWithPairs = binIdx.size();

		///// Get a vector of not empty bins
		std::vector< unsigned int> bins_touched;
		for (unsigned int i=0; i<2100; i++) {
			if (binTouched[i]) bins_touched.push_back(i);
		}
		const unsigned int nbBinsTouched = bins_touched.size();

		std::vector< unsigned int > nbBinsTouched_pixel_1(nb);
                std::vector< unsigned int > nbBinsTouched_pixel_2(nb);
		std::vector< std::vector< unsigned int > > bins_Touched_pixel_1(nb);
		std::vector< std::vector< unsigned int > > bins_Touched_pixel_2(nb);
		for (unsigned int i=0; i<nb; i++) {
			std::vector< unsigned int > tmp_bins_Touched_pixel_1;
			std::vector< unsigned int > tmp_bins_Touched_pixel_2;
			for (unsigned int b=0; b<2100; b++) {
				if (binTouched_1[i][b]) tmp_bins_Touched_pixel_1.push_back(b);
				if (binTouched_2[i][b]) tmp_bins_Touched_pixel_2.push_back(b);
			}
			bins_Touched_pixel_1[i] = tmp_bins_Touched_pixel_1;
			bins_Touched_pixel_2[i] = tmp_bins_Touched_pixel_2;
			nbBinsTouched_pixel_1[i] = tmp_bins_Touched_pixel_1.size();
                        nbBinsTouched_pixel_2[i] = tmp_bins_Touched_pixel_2.size();
		}

		///// Loops over all pixels of the forest (Fill the distortion matrix)
		for (unsigned int pa=0; pa<nbPixelsWithPairs; pa++) {

			const unsigned int globalBin1 = binIdx[pa];
			const unsigned int i1 = fromBinsToPixels1[pa];
			const unsigned int i2 = fromBinsToPixels2[pa];

			const double w1w2 = v_w__[f][i1]*v_w__[f][i2];
			const double vali1 = v_lRF__[f][i1] -v_meanLambda[f];
			const double vali2 = v_lRF__[f][i2] -v_meanLambda[f];

			///// Fill the distortion matrix
			dataMatrix[globalBin1][globalBin1] += w1w2;
			
			for (unsigned int b2=0; b2<nbBinsTouched; b2++) {
				const unsigned int globalBin2 = bins_touched[b2];
				dataMatrix[globalBin1][globalBin2] += w1w2*( value2[globalBin2]
					+vali1*value5[globalBin2] +vali2*value6[globalBin2]
					+vali1*vali2*value7[globalBin2] );
			}
			for (unsigned int b2=0; b2<nbBinsTouched_pixel_1[i1]; b2++) {
				const unsigned int globalBin2 = bins_Touched_pixel_1[i1][b2];
				dataMatrix[globalBin1][globalBin2] -= w1w2*( list_for_pixel_1_0[i1][globalBin2] + vali2*list_for_pixel_1_1[i1][globalBin2] );
			}
			for (unsigned int b2=0; b2<nbBinsTouched_pixel_2[i2]; b2++) {
				const unsigned int globalBin2 = bins_Touched_pixel_2[i2][b2];
				dataMatrix[globalBin1][globalBin2] -= w1w2*( list_for_pixel_2_0[i2][globalBin2] + vali1*list_for_pixel_2_1[i2][globalBin2] );
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	///// Save the 2D cross-correlation
	std::ofstream fFile;
        std::string pathToSave = pathToSave__;
	pathToSave += "xi_1D_delta_delta_distortionMatrix_";
	pathToSave += forest__;
	pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	///// Set the values of data
        ///// [0] for value, [1] for error, [2] for bin center
	long double value_total  = 0.;
	long double weight_total = 0.;
        for (unsigned int i=0; i<nbBins; i++) {
                for (unsigned int j=0; j<nbBins; j++) {
                        double value = 0.;
			if (weight[i]!=0.) {
				value = dataMatrix[i][j]/weight[i];
				value_total  += dataMatrix[i][j];
			}
                        fFile << value << " ";
                }
                fFile << std::endl;
		weight_total += weight[i];
        }
        fFile.close();

	std::cout << "  value_total     = " << value_total << std::endl;
	std::cout << "  weight_total    = " << weight_total << std::endl;
	std::cout << "  integral matrix = " << value_total / weight_total << std::endl;

	return;
}
void Correlation::xi_1DlRF_delta_delta_distortionMatrix(void) {

	return;
}
void Correlation::xi_1DlRFDevide_delta_delta_distortionMatrix(void) {

	return;
}
void Correlation::xi_1D_delta_delta2(void) {

	std::cout << "\n\n\n\n  ------ xi_1D_delta_delta2 ------\n" << std::endl;

	///// Forest 1
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	v_CosDe__.clear();
	v_SinDe__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	///// Forest 2
	loadDataDelta2();
	if (v_raDelta2__.size()==0.) return;

	v_CosDeDelta2__.clear();
	v_SinDeDelta2__.clear();
	v_lRFDelta2__.clear();
	v_lObsDelta2__.clear();

	///// Constants:
	///// The space between bins is of 1 Mpc.h^-1, the first bin is for deltaR == 0.
	const double max = 1000.;
	const unsigned int nbBins = 2*int(max);

	///// Arrays for data
	double data[nbBins][6];
	for (unsigned int i=0; i<nbBins; i++) {
		for (unsigned int j=0; j<6; j++) {
			data[i][j] = 0.;
		}
	}

	unsigned int nb_same_forest = 0;
	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {
		
		const double ra1       = v_ra__[f1];
		const double de1       = v_de__[f1];
		const double zz1       = v_zz__[f1];
		const unsigned int nb1 = v_nbPixelDelta1__[f1];
		
		for (unsigned int f2=0; f2<nbForest2__; f2++) {

			if (ra1!=v_raDelta2__[f2] || de1!=v_deDelta2__[f2] || zz1!=v_zzDelta2__[f2]) continue;
			nb_same_forest ++;
			
			const unsigned int nb2 = v_nbPixelDelta2__[f2];

			for (unsigned int i1=0; i1<nb1; i1++) {
		
				const double r1 = v_r__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double w1 = v_w__[f1][i1];
				const double z1 = v_z__[f1][i1];
	
				for (unsigned int i2=0; i2<nb2; i2++) {
					
					const double r1r2 = r1-v_rDelta2__[f2][i2];
					const double w1w2 = w1*v_wDelta2__[f2][i2];
					const double d1d2 = d1*v_dDelta2__[f2][i2];
					const double z1z2 = z1+v_zDelta2__[f2][i2];
	
					///// Find bin index
					const unsigned int idx = int(max+r1r2);

					///// Compute the xi1D
					data[idx][0] += w1w2*d1d2;
					data[idx][1] += w1w2*d1d2*d1d2;
					data[idx][2] += w1w2*r1r2;
					data[idx][3] += w1w2*z1z2;
					data[idx][4] += w1w2;
					data[idx][5] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::cout << "  nb of same forest = " << nb_same_forest << std::endl;
	if (nb_same_forest==0) return;

	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_1D_delta_delta2_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += forest2__;
	pathToSave += ".txt";
	std::cout << "\n  pathToSave = " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	///// Set the values of data
	for (unsigned int i=0; i<nbBins; i++) {
		fFile << data[i][0];
		fFile << " " << data[i][1];
		fFile << " " << data[i][2];
		fFile << " " << data[i][3]/2.;
		fFile << " " << data[i][4];
		fFile << " " << data[i][5];
		fFile << std::endl;

		sumOne += data[i][5];
		sumZZZ += data[i][3]/2.;
		sumWeg += data[i][4];
	}
	fFile.close();

	std::cout << "  < z >           = " << sumZZZ/sumWeg << std::endl;
	std::cout << "  Number of pairs = " << sumOne << std::endl;
	
	return;
}
void Correlation::xi_1DlRFDevide_delta_delta2(void) {

	std::cout << "\n\n\n\n  ------ xi_1DlRF_delta_delta2 ------\n" << std::endl;

	///// Forest 1
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	v_CosDe__.clear();
	v_SinDe__.clear();
	v_idx__.clear();
	v_r__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	///// Forest 2
	loadDataDelta2();
	if (v_raDelta2__.size()==0.) return;

	v_CosDeDelta2__.clear();
	v_SinDeDelta2__.clear();
	v_rDelta2__.clear();
	v_lObsDelta2__.clear();

	///// Constants:
	///// The space between bins is of 1 Mpc.h^-1, the first bin is for deltaR == 0.
	const double min = lambdaRFMin__/lambdaRFMaxDelta2__;
	const double max = lambdaRFMax__/lambdaRFMinDelta2__;
	const unsigned int nbBins = 50000;
	const double binSize = (max-min)/nbBins;

	///// Arrays for data
	double data[nbBins][6];
	for (unsigned int i=0; i<nbBins; i++) {
		for (unsigned int j=0; j<6; j++) {
			data[i][j] = 0.;
		}
	}

	unsigned int nb_same_forest = 0;
	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {
		
		const double ra1       = v_ra__[f1];
		const double de1       = v_de__[f1];
		const double zz1       = v_zz__[f1];
		const unsigned int nb1 = v_nbPixelDelta1__[f1];

		for (unsigned int f2=0; f2<nbForest2__; f2++) {

			if (ra1!=v_raDelta2__[f2] || de1!=v_deDelta2__[f2] || zz1!=v_zzDelta2__[f2]) continue;
			nb_same_forest ++;
			
			const unsigned int nb2 = v_nbPixelDelta2__[f2];

			for (unsigned int i1=0; i1<nb1; i1++) {
		
				const double l1 = v_lRF__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double w1 = v_w__[f1][i1];
				const double z1 = v_z__[f1][i1];
	
				for (unsigned int i2=0; i2<nb2; i2++) {
					
					const double l1l2 = l1/v_lRFDelta2__[f2][i2];
					const double w1w2 = w1*v_wDelta2__[f2][i2];
					const double d1d2 = d1*v_dDelta2__[f2][i2];
					const double z1z2 = z1+v_zDelta2__[f2][i2];
	
					///// Find bin index
					const unsigned int idx = int((l1l2-min)/binSize);
	
					///// Compute the xi1D
					data[idx][0] += w1w2*d1d2;
					data[idx][1] += w1w2*d1d2*d1d2;
					data[idx][2] += w1w2*l1l2;
					data[idx][3] += w1w2*z1z2;
					data[idx][4] += w1w2;
					data[idx][5] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::cout << "  nb of same forest = " << nb_same_forest << std::endl;
	if (nb_same_forest==0) return;

	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_1DlRFDevide_delta_delta2_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += forest2__;
	pathToSave += ".txt";
	std::cout << "\n  pathToSave = " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	///// Set the values of data
	for (unsigned int i=0; i<nbBins; i++) {
		fFile << data[i][0];
		fFile << " " << data[i][1];
		fFile << " " << data[i][2];
		fFile << " " << data[i][3]/2.;
		fFile << " " << data[i][4];
		fFile << " " << data[i][5];
		fFile << std::endl;

		sumOne += data[i][5];
		sumZZZ += data[i][3]/2.;
		sumWeg += data[i][4];
	}
	fFile.close();

	std::cout << "  < z >           = " << sumZZZ/sumWeg << std::endl;
	std::cout << "  Number of pairs = " << sumOne << std::endl;
	
	return;
}





// ---------------------------------------------------------------------
//
//		3D, forest 1 - forest 2 correlation
//
// ---------------------------------------------------------------------
void Correlation::xi_A_delta_delta(unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta ------\n" << std::endl;

        if (doBootstraps__)  std::cout << "  subsampling N° "     << bootIdx << std::endl;
        if ( shuffleForest ) std::cout << "  shuffleForest seed " << bootIdx*10 << std::endl;
        if ( shuffleQSO  )   std::cout << "  shuffleQSO seed    " << bootIdx*10 << std::endl;
        if ( randomQSO )     std::cout << "  randomQSO seed     " << bootIdx*10 << std::endl;

	///// Load forest
	loadDataForest(pathForest__, bootIdx);
	if (nbForest_==0.) return;

	v_zz__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Needed to randomize the forests
	if (shuffleForest) {

		///// Vectors with index of forest
		std::vector<unsigned int> randomIdx(nbForest_);
		for (unsigned int i=0; i<nbForest_; i++) {
			randomIdx[i] = i;
		}

		std::cout << "  Shuffle forest " << bootIdx << std::endl;
		bool doLoop = true;
		unsigned int nbLoop = 0;
		std::srand (bootIdx*10);
		while (doLoop) {
			doLoop = false;
			nbLoop ++;
			std::cout << "  nbLoop = " << nbLoop << std::endl;
			std::random_shuffle ( randomIdx.begin(), randomIdx.end() );
			for (unsigned int i=0; i<nbForest_; i++) {
				if (i==randomIdx[i]) {
					std::cout << nbLoop << " " << i << std::endl;
					doLoop = true;
					break;
				}
			}
		}

		///// Copy data
		std::vector<double> tmp_cosDe(v_CosDe__);
		std::vector<double> tmp_SinDe(v_SinDe__);
		std::vector<double> tmp_ra(v_ra__);
		std::vector<double> tmp_de(v_de__);

		///// Put the new data
		for (unsigned int i=0; i<nbForest_; i++) {
			const unsigned int ii = randomIdx[i];
			v_CosDe__[i] = tmp_cosDe[ii];
			v_SinDe__[i] = tmp_SinDe[ii];
			v_ra__[i]    = tmp_ra[ii];
			v_de__[i]    = tmp_de[ii];
		}
	}

	///// Constants:
	const double  max         = 200.;
	const unsigned int nbBin  = int(max);
	const unsigned int nbBinM = 50;

	const double maxPow2      = max*max;
	const double distMax      = sqrt(2.)*max+1.;
	const double distMaxPow2  = 2.*max*max;

	///// Find the maximal angle of seperation
	const double maxTheta = 2.*asin( distMax/(distMinPixel__*2.) );
	const double cosmaxTheta = cos(maxTheta);
	std::cout << "  maxTheta = " << maxTheta*180./M_PI << " degree" << std::endl;


	///// get an array for nb of pairs for the forest
	double a_nbPairs[nbForest_];
	for (unsigned int i=0; i<nbForest_; i++) {
		a_nbPairs[i] = 0.;
	}

	///// Arrays for data
	double dataMu[nbBin][nbBinM][9];
	double data2D[nbBin][nbBin][9];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<9; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<9; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}
	
	std::cout << "\n  Starting\n" << std::endl;			
	
	for (unsigned int f1=0; f1<nbForest_; f1++) {

		if (doBootstraps__ && v_region_Map__[f1]!=bootIdx) continue;

		const unsigned int nb1      = v_nbPixelDelta1__[f1];
		const double cosDe          = v_CosDe__[f1];
		const double sinDe          = v_SinDe__[f1];
		const double firstPixelPow2 = v_r__[f1][0]*v_r__[f1][0];
		const double ra1            = v_ra__[f1];
		const double de1            = v_de__[f1];

		for (unsigned int f2=0; f2<f1; f2++) {

			///// Usefull for eBOSS data because many re-obs.
			if (ra1==v_ra__[f2] && de1==v_de__[f2]) continue;
			else if (fabs(ra1-v_ra__[f2])<C_AUTOCORRCRIT && fabs(de1-v_de__[f2])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDe__[f2]*cos(ra1-v_ra__[f2]) + sinDe*v_SinDe__[f2];
			if (cosTheta<cosmaxTheta) continue;

			
			const double distMinPow2 = std::min(firstPixelPow2,v_r__[f2][0]*v_r__[f2][0])*(1.-cosTheta*cosTheta);
			if (distMinPow2>distMaxPow2) continue;

			const double sinTheta = sqrt(1.-cosTheta*cosTheta);
			const unsigned int nb2 = v_nbPixelDelta1__[f2];

			for (unsigned int i1=0; i1<nb1; i1++) {

				///// Get the r_perp distance
				const double rPerp = v_r__[f1][i1]*sinTheta;
				if (rPerp>=max) continue;
				const unsigned int rPerpBinIdx = int(rPerp);
				const double tmp_rParral = v_r__[f1][i1]*cosTheta;

				const double w1        = v_w__[f1][i1];
				const double d1        = v_d__[f1][i1];
				const double z1        = v_z__[f1][i1];
				const double res_lRF1  = v_residual_delta_vs_lRF__[f1][i1];
				const double res_lObs1 = v_residual_delta_vs_lObs__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-tmp_rParral);
					if (rParral>=max) continue;
					const double distTotPow2 = rPerp*rPerp + rParral*rParral;

					const double w1w2                 = w1*v_w__[f2][i2];
					const double d1d2                 = d1*v_d__[f2][i2];
					const double w1w2_d1d2            = w1w2*d1d2;
					const double w1w2_d1d2_d1d2       = w1w2_d1d2*d1d2;
					const double w1w2_z1z2            = w1w2*(z1+v_z__[f2][i2]);
					//const double w1w2_res_lRF1_lRF2   = w1w2*res_lRF1*v_d__[f2][i2];
					//const double w1w2_res_lObs1_lObs2 = w1w2*res_lObs1*v_d__[f2][i2];
					const double w1w2_res_lRF1_lRF2   = w1w2*res_lRF1*v_residual_delta_vs_lRF__[f2][i2];
					const double w1w2_res_lObs1_lObs2 = w1w2*res_lObs1*v_residual_delta_vs_lObs__[f2][i2];

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int(mu*50.);

						dataMu[idx][idxM][0] += w1w2_d1d2;
						dataMu[idx][idxM][1] += w1w2_d1d2_d1d2;
						dataMu[idx][idxM][2] += w1w2*distTot;
						dataMu[idx][idxM][3] += w1w2*mu;
						dataMu[idx][idxM][4] += w1w2_z1z2;
						dataMu[idx][idxM][5] += w1w2;
						dataMu[idx][idxM][6] ++;
						dataMu[idx][idxM][7] += w1w2_res_lRF1_lRF2;
						dataMu[idx][idxM][8] += w1w2_res_lObs1_lObs2;
					}
					
					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += w1w2_d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][1] += w1w2_d1d2_d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2_z1z2;
					data2D[rPerpBinIdx][rParralBinIdx][5] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][6] ++;
					data2D[rPerpBinIdx][rParralBinIdx][7] += w1w2_res_lRF1_lRF2;
					data2D[rPerpBinIdx][rParralBinIdx][8] += w1w2_res_lObs1_lObs2;

					///// Get the number of pairs
					a_nbPairs[f1] += w1w2;
					a_nbPairs[f2] += w1w2;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	//// Set the prefix of different forest
	std::string prefix1 = forest__;

	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	//// Set the prefix for different type of runs
	std::string prefix = "_";
	if (doBootstraps__)  prefix += "subsampling";
	if (shuffleForest) prefix += "shuffleForest";
	if (shuffleQSO)    prefix += "shuffleQSO";
	if (randomQSO)     prefix += "randomQSO";
	prefix += "_";
	prefix += strBootIdx;

	std::ofstream fFile;
	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	///// Save the 2D cross-correlation
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_2D_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4]/2.;
			fFile << " " << data2D[i][j][5];
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << " " << data2D[i][j][8];
			fFile << std::endl;

			sumZZZ += data2D[i][j][4]/2.;
			sumWeg += data2D[i][j][5];
			sumOne += data2D[i][j][6];
		}
	}
	fFile.close();


	std::cout << "  < z >                = " << sumZZZ/sumWeg << std::endl;
	std::cout << "  number of pairs      = " << sumOne          << std::endl;
	std::cout << "  < pairs per forest > = " << sumOne/nbForest_ << std::endl;


	///// Mu
	pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_Mu_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << " " << dataMu[i][j][8];
			fFile << std::endl;
		}
	}
	fFile.close();


	if (!doBootstraps__ && !shuffleForest && !shuffleQSO && !randomQSO) {

		std::vector< std::vector< double > > forests;
		std::vector< double > tmp_forests_id;
		std::vector< double > tmp_forests_pa;
		
		for (unsigned int i=0; i<nbForest_; i++) {
			tmp_forests_id.push_back( v_idx__[i] );
			tmp_forests_pa.push_back( a_nbPairs[i] );
		}
		forests.push_back(tmp_forests_id);
		forests.push_back(v_ra__);
		forests.push_back(v_de__);
		forests.push_back(tmp_forests_pa);


		// Save the list of pairs
		pathToSave = pathToSave__;
		pathToSave += "xi_A_delta_delta_list_pairs_";
		pathToSave += prefix1;
		pathToSave += ".txt";
		std::cout << "\n  " << pathToSave << std::endl;

		std::ofstream fFile;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);

		for (unsigned int i=0; i<nbForest_; i++) {
			fFile << v_idx__[i];
			fFile << " " << 0;
			fFile << " " << v_ra__[i];
			fFile << " " << v_de__[i];
			fFile << " " << a_nbPairs[i];
			fFile << std::endl;
		}
		fFile.close();
		

		//// find the index of each forest among the C_NBSUBSAMPLES sub-samples
		LymanForest* lymanForestObject = new LymanForest(forests, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		pathToSave = pathToSave__;
		pathToSave += "xi_A_delta_delta_map_";
		pathToSave += forest__;
		pathToSave += ".txt";
		lymanForestObject->SaveRegionMap(pathToSave);
		delete lymanForestObject;
	}


	return;
}
void Correlation::xi_A_delta_delta_lambda(void) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_lambda ------\n" << std::endl;

	///// Load forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	v_zz__.clear();
	v_r__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_nb__.clear();

	///// Rapport l1/l2
	const double minl1l2 = 1.;
	const double maxl1l2 = 1.11;
	double binSizel1l2      = 1.e-03;
	const unsigned int nbBin = int( (maxl1l2-minl1l2)/binSizel1l2);
	std::cout << "  nbBin = " << nbBin << std::endl;
	binSizel1l2 = (maxl1l2-minl1l2)/nbBin;
	const double inv_binSizel1l2 = 1./binSizel1l2;
	std::cout << "  nbBin = " << nbBin<<std::endl;

	///// Theta
	const double minTheta = 0.;
	const double maxTheta = 0.003;
	const double cosThetaMax = cos(maxTheta);
	double binSizeTheta   = 1.e-04;
	const unsigned int nbBinM = int( (maxTheta-minTheta)/binSizeTheta);
	binSizeTheta = (maxTheta-minTheta)/nbBinM;
	const double inv_binSizeTheta = 1./binSizeTheta;
	std::cout << "  nbBinM = " << nbBinM << std::endl;


	///// Arrays for data
	double dataMu[nbBin][nbBinM][7];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;	
	
	for (unsigned int f1=0; f1<nbForest_; f1++) {
		const unsigned int nb1      = v_nbPixelDelta1__[f1];
		const double cosDe          = v_CosDe__[f1];
		const double sinDe          = v_SinDe__[f1];
		const double ra1            = v_ra__[f1];

		for (unsigned int f2=0; f2<f1; f2++) {

			///// Usefull for eBOSS data because many re-obs.
			if (v_ra__[f1]==v_ra__[f2] && v_CosDe__[f1]==v_CosDe__[f2] && v_SinDe__[f1]==v_SinDe__[f2]) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDe__[f2]*cos(ra1-v_ra__[f2]) + sinDe*v_SinDe__[f2];
			if (cosTheta<cosThetaMax) continue;
			const double theta = acos(cosTheta);
			const unsigned int idxTheta = int( theta*inv_binSizeTheta );

			///// Constants of forest n°2
			const unsigned int nb2 = v_nbPixelDelta1__[f2];

			for (unsigned int i1=0; i1<nb1; i1++) {

				const double w1 = v_w__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double z1 = v_z__[f1][i1];
				const double l1 = v_lObs__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					///// Same observed wavelength
					if (v_lObs__[f2][i2]==l1) continue;
					const double l1Overl2 = std::max(v_lObs__[f2][i2]/l1, l1/v_lObs__[f2][i2]);
					if (l1Overl2>maxl1l2) continue;

					const double w1w2 = w1*v_w__[f2][i2];
					const double d1d2 = d1*v_d__[f2][i2];
					const double z1z2 = z1+v_z__[f2][i2];
					const unsigned int idx  = int( (l1Overl2-minl1l2)*inv_binSizel1l2 );
					
					dataMu[idx][idxTheta][0] += w1w2*d1d2;
					dataMu[idx][idxTheta][1] += w1w2*d1d2*d1d2;
					dataMu[idx][idxTheta][2] += w1w2*l1Overl2;
					dataMu[idx][idxTheta][3] += w1w2*theta;
					dataMu[idx][idxTheta][4] += w1w2*z1z2;
					dataMu[idx][idxTheta][5] += w1w2;
					dataMu[idx][idxTheta][6] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;
	
	double sumOne = 0.;

	///// Mu
	std::ofstream fFile;
	std::string tmp_pathToSave = pathToSave__;
	tmp_pathToSave = pathToSave__;
	tmp_pathToSave += "xi_A_delta_delta_lambda_Mu_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += ".txt";
	std::cout << "\n  " << tmp_pathToSave << std::endl;
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << std::endl;
			sumOne += dataMu[i][j][6];
		}
	}
	fFile.close();

	std::cout << "  number of pairs      = " << sumOne          << std::endl;


	return;
}
void Correlation::xi_A_delta_delta_distortionMatrix(void) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_distortionMatrix ------\n" << std::endl;

	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_d__.clear();
	v_z__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	v_idx__.clear();

	///// Set usefull vectors
	std::vector<double> v_invSumWeight(nbForest_,0.);
	std::vector<double> v_meanLambda(nbForest_,0.);
	std::vector<double> v_invVarLambda(nbForest_,0.);

	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nbPixel = v_nbPixelDelta1__[f];

		long double sumWeight  = 0.;
		long double meanLambda = 0.;
		long double varLambda  = 0.;

		///// Loops over all pixels of the forest
		for (unsigned int i=0; i<nbPixel; i++) {

			const long double w = v_w__[f][i];
			const long double l = v_lRF__[f][i];

			sumWeight  += w;
			meanLambda += w*l;
			varLambda  += w*l*l;
		}

		if (sumWeight==0.) std::cout << "   ERROR:: sumWeight = 0. " << std::endl;
		meanLambda        /= sumWeight;
		varLambda          = varLambda - meanLambda*meanLambda*sumWeight;
		if (varLambda==0.) std::cout << "   ERROR:: varLambda = 0. " << std::endl;

		v_invSumWeight[f] = 1./sumWeight;
		v_meanLambda[f]   = meanLambda;
		v_invVarLambda[f] = 1./varLambda;
	}

	///// Constants:
	const double max         = 200.;
	const double binSize     = 4.;
	const unsigned int nbBin = int(max/binSize);
	const unsigned int nbBin2D = nbBin*nbBin;

	const double distMax      = sqrt(2.)*max+1.;
	const double distMaxPow2  = 2.*max*max;
	const double fromValToIdx  = nbBin/max;

	///// Find the maximal angle of seperation
	const double maxTheta = 2.*asin( distMax/(distMinPixel__*2.) );
	const double cosmaxTheta = cos(maxTheta);
	std::cout << "  maxTheta = " << maxTheta*180./M_PI << " degree" << std::endl;

	///// Arrays for distortion matrix
	double weight2D[nbBin2D];
	double** data2DMatrix = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		weight2D[i] = 0.;
		data2DMatrix[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			data2DMatrix[i][j] = 0.;
		}
	}





	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {

		const unsigned int nb1      = v_nbPixelDelta1__[f1];
		const double cosDe          = v_CosDe__[f1];
		const double sinDe          = v_SinDe__[f1];
		const double firstPixelPow2 = v_r__[f1][0]*v_r__[f1][0];
		const double ra1            = v_ra__[f1];
		const double de1            = v_de__[f1];
		const double invSumWeight1  = v_invSumWeight[f1];
		const double meanLambda1    = v_meanLambda[f1];
		const double invVarLambda1  = v_invVarLambda[f1];

		for (unsigned int f2=0; f2<f1; f2++) {

			///// Usefull for eBOSS data because many re-obs.
			if (ra1==v_ra__[f2] && de1==v_de__[f2]) continue;
			else if (fabs(ra1-v_ra__[f2])<C_AUTOCORRCRIT && fabs(de1-v_de__[f2])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDe__[f2]*cos(ra1-v_ra__[f2]) + sinDe*v_SinDe__[f2];
			if (cosTheta<cosmaxTheta) continue;

			const double distMinPow2 = std::min(firstPixelPow2,v_r__[f2][0]*v_r__[f2][0])*(1.-cosTheta*cosTheta);
			if (distMinPow2>distMaxPow2) continue;

			const double sinTheta = sqrt(1.-cosTheta*cosTheta);
			const unsigned int nb2 = v_nbPixelDelta1__[f2];
			const double invSumWeight2  = v_invSumWeight[f2];
			const double meanLambda2    = v_meanLambda[f2];
			const double invVarLambda2  = v_invVarLambda[f2];

			//// For the distortion-matrix
			std::vector< std::vector< double > > list_for_pixel_1_0(nb1, std::vector<double>(2500,0.));
			std::vector< std::vector< double > > list_for_pixel_2_0(nb2, std::vector<double>(2500,0.));
			std::vector< std::vector< double > > list_for_pixel_1_1(nb1, std::vector<double>(2500,0.));
			std::vector< std::vector< double > > list_for_pixel_2_1(nb2, std::vector<double>(2500,0.));
			std::vector< bool > pixelTouched_1(nb1,false);
			std::vector< bool > pixelTouched_2(nb2,false);
			std::vector< std::vector< bool > > binTouched_1(nb1, std::vector<bool>(2500,false));
			std::vector< std::vector< bool > > binTouched_2(nb2, std::vector<bool>(2500,false));
			double value2[2500]   = {0.};
			double value5[2500]   = {0.};
			double value6[2500]   = {0.};
			double value7[2500]   = {0.};
			bool binTouched[2500] = {false};

			///// Index of the bin
			std::vector< unsigned int > binIdx;
			std::vector< unsigned int > fromBinsToPixels1;
			std::vector< unsigned int > fromBinsToPixels2;

			for (unsigned int i1=0; i1<nb1; i1++) {

				///// Get the r_perp distance
				const double rPerp = v_r__[f1][i1]*sinTheta;
				if (rPerp>=max) continue;
				const unsigned int rPerpBinIdx           = int(rPerp*fromValToIdx);
				const double tmp_rParral                 = v_r__[f1][i1]*cosTheta;
				const double w1                          = v_w__[f1][i1];
				const double lambda_minus_mean_over_var1 = (v_lRF__[f1][i1] -meanLambda1)*invVarLambda1;

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-tmp_rParral);
					if (rParral>=max) continue;
					const unsigned int rParralBinIdx = int(rParral*fromValToIdx);
					const unsigned int globalBin     = rPerpBinIdx*nbBin +rParralBinIdx;

					const double w2                          = v_w__[f2][i2];
					const double lambda_minus_mean_over_var2 = (v_lRF__[f2][i2] -meanLambda2)*invVarLambda2;

					//// Fill arrays for the distortion
					list_for_pixel_1_0[i1][globalBin] += w2*invSumWeight2;
					list_for_pixel_2_0[i2][globalBin] += w1*invSumWeight1;
					list_for_pixel_1_1[i1][globalBin] += w2*lambda_minus_mean_over_var2;
					list_for_pixel_2_1[i2][globalBin] += w1*lambda_minus_mean_over_var1;
					pixelTouched_1[i1] = true;
					pixelTouched_2[i2] = true;
					binTouched_1[i1][globalBin] = true;
					binTouched_2[i2][globalBin] = true;
					value2[globalBin] += w1*w2*invSumWeight1*invSumWeight2;
					value5[globalBin] += w2*invSumWeight2*w1*lambda_minus_mean_over_var1;
					value6[globalBin] += w1*invSumWeight1*w2*lambda_minus_mean_over_var2;
					value7[globalBin] += w1*w2*lambda_minus_mean_over_var1*lambda_minus_mean_over_var2;
					binIdx.push_back(globalBin);
					fromBinsToPixels1.push_back(i1);
					fromBinsToPixels2.push_back(i2);
					binTouched[globalBin] = true;

					///// Fill array of weights
					weight2D[globalBin] += w1*w2;
				}
			}

			///// Number of pairs of two pixels ( here = nb*nb)
			const unsigned int nbPixelsWithPairs = binIdx.size();
	
			///// Get a vector of not empty bins
			std::vector< unsigned int> bins_touched;
			for (unsigned int b=0; b<2500; b++) {
				if (binTouched[b]) bins_touched.push_back(b);
			}
			const unsigned int nbBinsTouched = bins_touched.size();

			//// For forest 1
			std::vector< unsigned int > nbBinsTouched_pixel_1(nb1);
			std::vector< std::vector< unsigned int > > bins_Touched_pixel_1(nb1);
			for (unsigned int i=0; i<nb1; i++) {
				if (!pixelTouched_1[i]) continue;
				std::vector< unsigned int > tmp_bins_Touched_pixel_1;
				for (unsigned int b=0; b<2500; b++) {
					if (binTouched_1[i][b]) tmp_bins_Touched_pixel_1.push_back(b);
				}
				bins_Touched_pixel_1[i] = tmp_bins_Touched_pixel_1;
				nbBinsTouched_pixel_1[i] = tmp_bins_Touched_pixel_1.size();
			}
			//// For forest 2
			std::vector< unsigned int > nbBinsTouched_pixel_2(nb2);
			std::vector< std::vector< unsigned int > > bins_Touched_pixel_2(nb2);
			for (unsigned int i=0; i<nb2; i++) {
				if (!pixelTouched_2[i]) continue;
				std::vector< unsigned int > tmp_bins_Touched_pixel_2;
				for (unsigned int b=0; b<2500; b++) {
					if (binTouched_2[i][b]) tmp_bins_Touched_pixel_2.push_back(b);
				}
				bins_Touched_pixel_2[i] = tmp_bins_Touched_pixel_2;
				nbBinsTouched_pixel_2[i] = tmp_bins_Touched_pixel_2.size();
			}
	
			///// Loops over all pairs
			for (unsigned int p1=0; p1<nbPixelsWithPairs; p1++) {
	
				const unsigned int globalBin1 = binIdx[p1];
				const unsigned int i1 = fromBinsToPixels1[p1];
				const unsigned int i2 = fromBinsToPixels2[p1];
	
				const double w1w2 = v_w__[f1][i1]*v_w__[f2][i2];
				const double valf1i1 = v_lRF__[f1][i1] -v_meanLambda[f1];
				const double valf2i2 = v_lRF__[f2][i2] -v_meanLambda[f2];
	
				///// Fill the distortion matrix
				data2DMatrix[globalBin1][globalBin1] += w1w2;

				/// Loop over all not empty bins
				for (unsigned int b2=0; b2<nbBinsTouched; b2++) {
					const unsigned int globalBin2 = bins_touched[b2];
					data2DMatrix[globalBin1][globalBin2] += w1w2*( value2[globalBin2]
						+valf1i1*value5[globalBin2] +valf2i2*value6[globalBin2]
						+valf1i1*valf2i2*value7[globalBin2] );
				}
				
				for (unsigned int b2=0; b2<nbBinsTouched_pixel_1[i1]; b2++) {
					const unsigned int globalBin2 = bins_Touched_pixel_1[i1][b2];
					data2DMatrix[globalBin1][globalBin2] -= w1w2*( list_for_pixel_1_0[i1][globalBin2] + valf2i2*list_for_pixel_1_1[i1][globalBin2] );
				}
				for (unsigned int b2=0; b2<nbBinsTouched_pixel_2[i2]; b2++) {
					const unsigned int globalBin2 = bins_Touched_pixel_2[i2][b2];
					data2DMatrix[globalBin1][globalBin2] -= w1w2*( list_for_pixel_2_0[i2][globalBin2] + valf1i1*list_for_pixel_2_1[i2][globalBin2] );
				}
			}
		}
std::cout << f1 << std::endl;
	}


        std::cout << "\n  Saving\n" << std::endl;

        std::ofstream fFile;

	///// Save the 2D cross-correlation 
        std::string pathToSave = pathToSave__;
        pathToSave += "xi_A_delta_delta_distortionMatrix_2D_";
        pathToSave += forest__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	///// Set the values of data
        ///// [0] for value, [1] for error, [2] for bin center
	long double value_total  = 0.;
	long double weight_total = 0.;
        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
                        double value = 0.;
			if (weight2D[i]!=0.) {
				value = data2DMatrix[i][j]/weight2D[i];
				value_total  += data2DMatrix[i][j];
			}
                        fFile << value << " ";
                }
                fFile << std::endl;
		weight_total += weight2D[i];
        }
        fFile.close();

	std::cout << "  value_total     = " << value_total << std::endl;
	std::cout << "  weight_total    = " << weight_total << std::endl;
	std::cout << "  integral matrix = " << value_total / weight_total << std::endl;


	return;
}
void Correlation::xi_A_delta_delta_Metals_Models(double lambdaRFMetal1, std::string lambdaRFMetalName1,double lambdaRFMetal2, std::string lambdaRFMetalName2) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_Metals_Models ------\n" << std::endl;

	std::cout << "  line1 = " << lambdaRFMetalName1 << " : " << lambdaRFMetal1 << std::endl;
	std::cout << "  line2 = " << lambdaRFMetalName2 << " : " << lambdaRFMetal2 << std::endl;

	///// Load forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
	
	/// Get an arrey to get te growth factor
	double from_z_to_growth_factor_pow_2[C_NBBINREDSH] = {0.};
	const double g0    = cosmo->get_growth_factor(0.);
	const double growth_factor_step_size = C_ZEXTREMABINCONVERT1/C_NBBINREDSH;
	const double inverse_growth_factor_step_size = 1./growth_factor_step_size;
	for (unsigned int i=0; i<C_NBBINREDSH; i++) {
		const double ZZZ = i*growth_factor_step_size;
		from_z_to_growth_factor_pow_2[i] = (cosmo->get_growth_factor(ZZZ)/g0)*(cosmo->get_growth_factor(ZZZ)/g0);
	}

	//// Get the distance if it was this metal
	std::vector< std::vector< double > > v_r_metal1(v_r__);
	std::vector< std::vector< double > > v_r_metal2(v_r__);
	for (unsigned int i=0; i<v_r__.size(); i++) {
		for (unsigned int j=0; j<v_r__[i].size(); j++) {
			v_r_metal1[i][j] = hConvertRedshDist->Interpolate( v_lObs__[i][j]/lambdaRFMetal1-1. );
			v_r_metal2[i][j] = hConvertRedshDist->Interpolate( v_lObs__[i][j]/lambdaRFMetal2-1. );
		}
	}

	//// Empty useless vectors
	delete cosmo;
	delete hConvertRedshDist;
	v_zz__.clear();
	v_idx__.clear();
	v_d__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	//// Get the monopole, quadrupol, exadecapol
	std::string pathToLoad = PATHTOWORK;
	pathToLoad += C_PATHTOCAMB;
	std::ifstream fileData(pathToLoad.c_str());
	std::vector< double > data_x;
	std::vector< double > data_xi0;
	std::vector< double > data_xi2;
	std::vector< double > data_xi4;
	std::cout << pathToLoad << std::endl;
	while (fileData) {
		double x;
		double xi0;
		double xi2;
		double xi4;
		fileData>>x>>xi0>>xi2>>xi4;
		if (fileData==0) break;
		
		data_x.push_back(x);
		data_xi0.push_back(xi0);
		data_xi2.push_back(xi2);
		data_xi4.push_back(xi4);
	}

	const double inverse_step_size = 1./data_x[0];
	const unsigned int nbBinCAMB = data_x.size();
	//// Set to zero the last pixel
	const double maxDistCAMB = data_x[nbBinCAMB-1];
	data_x.push_back(data_x[nbBinCAMB-1]+1.);
	data_xi0.push_back(0.);
	data_xi2.push_back(0.);
	data_xi4.push_back(0.);

	///// Constants:
	const double max          = 200.;
	const unsigned int nbBin  = int(max);
	const unsigned int nbBinM = 50.;

	const double maxPow2      = max*max;
	const double distMax      = sqrt(2.)*max+1.;
	const double distMaxPow2  = 2.*max*max;

	///// Find the maximal angle of seperation
	const double maxTheta = 2.*asin( distMax/(distMinPixel__*2.) );
	const double cosmaxTheta = cos(maxTheta);
	std::cout << "  maxTheta = " << maxTheta*180./M_PI << " degree" << std::endl;

	///// Arrays for data
	double dataMu[nbBin][nbBinM][8];
	double data2D[nbBin][nbBin][8];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<8; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<8; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}
	
	std::cout << "\n  Starting\n" << std::endl;			
	
	for (unsigned int f1=0; f1<nbForest_; f1++) {
		const unsigned int nb1      = v_nbPixelDelta1__[f1];
		const double cosDe          = v_CosDe__[f1];
		const double sinDe          = v_SinDe__[f1];
		const double firstPixelPow2 = v_r__[f1][0]*v_r__[f1][0];
		const double ra1            = v_ra__[f1];
		const double de1            = v_de__[f1];

		for (unsigned int f2=0; f2<f1; f2++) {

			///// Usefull for eBOSS data because many re-obs.
			if (ra1==v_ra__[f2] && de1==v_de__[f2]) continue;
			else if (fabs(ra1-v_ra__[f2])<C_AUTOCORRCRIT && fabs(de1-v_de__[f2])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDe__[f2]*cos(ra1-v_ra__[f2]) + sinDe*v_SinDe__[f2];
			if (cosTheta<cosmaxTheta) continue;

			
			const double distMinPow2 = std::min(firstPixelPow2,v_r__[f2][0]*v_r__[f2][0])*(1.-cosTheta*cosTheta);
			if (distMinPow2>distMaxPow2) continue;

			const double sinTheta = sqrt(1.-cosTheta*cosTheta);
			const unsigned int nb2 = v_nbPixelDelta1__[f2];

			for (unsigned int i1=0; i1<nb1; i1++) {

				///// Get the r_perp distance
				const double rPerp = v_r__[f1][i1]*sinTheta;
				if (rPerp>=max) continue;
				const unsigned int rPerpBinIdx = int(rPerp);
				const double tmp_rParral = v_r__[f1][i1]*cosTheta;

				const double rm1 = v_r_metal1[f1][i1];
				const double w1  = v_w__[f1][i1];
				const double z1  = v_z__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-tmp_rParral);
					if (rParral>=max) continue;
					const double distTotPow2 = rPerp*rPerp + rParral*rParral;
					const double distTotMetal = sqrt( rPerp*rPerp + (v_r_metal2[f2][i2]-rm1)*(v_r_metal2[f2][i2]-rm1) );

					const double w1w2 = w1*v_w__[f2][i2];
					const double wz   = w1w2*(z1+v_z__[f2][i2]);

					unsigned int idxBinCAMB = nbBinCAMB;
					if (distTotMetal<maxDistCAMB) idxBinCAMB = int(distTotMetal*inverse_step_size);
					const double growth_factor_pow_2 = from_z_to_growth_factor_pow_2[ int(0.5*(z1+v_z__[f2][i2])*inverse_growth_factor_step_size) ];
					const double wxi0 = w1w2*growth_factor_pow_2*data_xi0[idxBinCAMB];
					const double wxi2 = w1w2*growth_factor_pow_2*data_xi2[idxBinCAMB];
					const double wxi4 = w1w2*growth_factor_pow_2*data_xi4[idxBinCAMB];

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int(mu*50.);

						dataMu[idx][idxM][0] += wxi0;
						dataMu[idx][idxM][1] += wxi2;
						dataMu[idx][idxM][2] += wxi4;
						dataMu[idx][idxM][3] += w1w2*distTot;
						dataMu[idx][idxM][4] += w1w2*mu;
						dataMu[idx][idxM][5] += wz;
						dataMu[idx][idxM][6] += w1w2;
						dataMu[idx][idxM][7] ++;
					}

					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += wxi0;
					data2D[rPerpBinIdx][rParralBinIdx][1] += wxi2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += wxi4;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][5] += wz;
					data2D[rPerpBinIdx][rParralBinIdx][6] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][7] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	//// Set the prefix of different forest and QSOs
	std::string prefix1 = forest__;
	prefix1 += "_";
	prefix1 += lambdaRFMetalName1;
	prefix1 += "_";
	prefix1 += lambdaRFMetalName2;

	std::ofstream fFile;
	std::string pathToSave;

	///// Save the 2D cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_Metals_Models_2D_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4];
			fFile << " " << data2D[i][j][5]/2.;
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();


	///// Mu
	pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_Metals_Models_Mu_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4];
			fFile << " " << dataMu[i][j][5]/2.;
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();

	return;
}
void Correlation::xi_A_delta_delta2( unsigned int bootIdx/*=0*/ ) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta2 ------\n" << std::endl;

	///// Forest 1
	loadDataForest(pathForest__, bootIdx);
	if (nbForest_==0.) return;

	if (doBootstraps__) removeFalseCorrelations();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	///// Forest 2
	loadDataDelta2();
	v_lRFDelta2__.clear();
	v_lObsDelta2__.clear();

	///// Constants:
	const double max          = 200.;
	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100;

	const double maxPow2      = max*max;
	const double distMax      = sqrt(2.)*max+1.;
	const double distMaxPow2  = 2.*max*max;

	///// Find the maximal angle of seperation
	const double maxTheta = 2.*asin( distMax/(distMinPixel__+distMinPixelDelta2__) );
	const double cosmaxTheta = cos(maxTheta);
	std::cout << "  maxTheta = " << maxTheta*180./M_PI << " degree" << std::endl;

	///// Arrays for data
	double dataMu[nbBin][nbBinM][9];
	double data2D[nbBinX][nbBinY][9];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<9; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<9; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}

	///// get an array for nb of pairs for the forest 1
	double a_nbPairs[nbForest_];
	for (unsigned int i=0; i<nbForest_; i++) {
		a_nbPairs[i] = 0.;
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {

		double nbPairs = 0.;

		const unsigned int nb1      = v_nbPixelDelta1__[f1];
		const double cosDe          = v_CosDe__[f1];
		const double sinDe          = v_SinDe__[f1];
		const double firstPixelPow2 = v_r__[f1][0]*v_r__[f1][0];
		const double ra1            = v_ra__[f1];
		const double de1            = v_de__[f1];
		const double zz1            = v_zz__[f1];

		for (unsigned int f2=0; f2<nbForest2__; f2++) {
		
			///// Remove the correlation qso-ownForest
			if (ra1==v_raDelta2__[f2] && de1==v_deDelta2__[f2] && zz1==v_zzDelta2__[f2]) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeDelta2__[f2]*cos(ra1-v_raDelta2__[f2]) + sinDe*v_SinDeDelta2__[f2];
			if (cosTheta<cosmaxTheta) continue;

			const double distMinPow2 = std::min(firstPixelPow2,v_rDelta2__[f2][0]*v_rDelta2__[f2][0])*(1.-cosTheta*cosTheta);
			if (distMinPow2>distMaxPow2) continue;

			const double sinTheta = sqrt(1.-cosTheta*cosTheta);
			const unsigned int nb2 = v_nbPixelDelta2__[f2];

			for (unsigned int i1=0; i1<nb1; i1++) {

				///// Get the r_perp distance
				const double rPerp = v_r__[f1][i1]*sinTheta;
				if (rPerp>=max) continue;
				const unsigned int rPerpBinIdx = int(rPerp);
				const double tmp_rParral = v_r__[f1][i1]*cosTheta;

				const double w1 = v_w__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double z1 = v_z__[f1][i1];
				const double res_lRF1  = v_residual_delta_vs_lRF__[f1][i1];
				const double res_lObs1 = v_residual_delta_vs_lObs__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = v_rDelta2__[f2][i2]-tmp_rParral;
					if (fabs(rParral)>=max) continue;
					const double distTotPow2 = rPerp*rPerp + rParral*rParral;

					const double w1w2                 = w1*v_wDelta2__[f2][i2];
					const double d1d2                 = d1*v_dDelta2__[f2][i2];
					const double w1w2_d1d2            = w1w2*d1d2;
					const double w1w2_d1d2_d1d2       = w1w2_d1d2*d1d2;
					const double w1w2_z1z2            = w1w2*(z1+v_zDelta2__[f2][i2]);
					const double w1w2_res_lRF1_lRF2   = w1w2*res_lRF1*v_residual2_delta_vs_lRF__[f2][i2];
					const double w1w2_res_lObs1_lObs2 = w1w2*res_lObs1*v_residual2_delta_vs_lObs__[f2][i2];

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int((1.+mu)*50.);

						dataMu[idx][idxM][0] += w1w2_d1d2;
						dataMu[idx][idxM][1] += w1w2_d1d2_d1d2;
						dataMu[idx][idxM][2] += w1w2*distTot;
						dataMu[idx][idxM][3] += w1w2*mu;
						dataMu[idx][idxM][4] += w1w2_z1z2;
						dataMu[idx][idxM][5] += w1w2;
						dataMu[idx][idxM][6] ++;
						dataMu[idx][idxM][7] += w1w2_res_lRF1_lRF2;
						dataMu[idx][idxM][8] += w1w2_res_lObs1_lObs2;
					}
					
					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(max+rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += w1w2_d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][1] += w1w2_d1d2_d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2_z1z2;
					data2D[rPerpBinIdx][rParralBinIdx][5] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][6] ++;
					data2D[rPerpBinIdx][rParralBinIdx][7] += w1w2_res_lRF1_lRF2;
					data2D[rPerpBinIdx][rParralBinIdx][8] += w1w2_res_lObs1_lObs2;

					///// Get the number of pairs
					nbPairs += w1w2;
				}
			}
		}

		///// Put the number of pairs comming with this forest
		a_nbPairs[f1] = nbPairs;
	}

	std::cout << "\n  Saving\n" << std::endl;
	
	
	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	//// Set the prefix for different type of runs
	std::string prefix = "_";
	if (doBootstraps__)  prefix += "subsampling";
	if (shuffleForest) prefix += "shuffleForest";
	if (shuffleQSO)    prefix += "shuffleQSO";
	if (randomQSO)     prefix += "randomQSO";
	prefix += "_";
	prefix += strBootIdx;
	

	std::string prefix1 = forest__;
	prefix1 += "_";
	prefix1 += forest2__;

	std::ofstream fFile;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;
	long double sumOne = 0.;

	///// Save the 2D cross-correlation
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta2_2D_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	sumZZZ = 0.;
	sumWeg = 0.;
	sumOne = 0.;
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4]/2.;
			fFile << " " << data2D[i][j][5];
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << " " << data2D[i][j][8];
			fFile << std::endl;

			sumZZZ += data2D[i][j][4]/2.;
			sumWeg += data2D[i][j][5];
			sumOne += data2D[i][j][6];
		}
	}
	fFile.close();

	std::cout << "  < z >                = " << sumZZZ/sumWeg    << std::endl;
	std::cout << "  number of pairs      = " << sumOne           << std::endl;
	std::cout << "  < pairs per forest > = " << sumOne/nbForest_ << std::endl;


	///// Mu
	pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta2_Mu_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << " " << dataMu[i][j][8];
			fFile << std::endl;
		}
	}
	fFile.close();


	if (!doBootstraps__ && !shuffleForest && !shuffleQSO && !randomQSO) {

		std::vector< std::vector< double > > forests;
		std::vector< double > tmp_forests_id;
		std::vector< double > tmp_forests_pa;
		
		for (unsigned int i=0; i<nbForest_; i++) {
			tmp_forests_id.push_back( v_idx__[i] );
			tmp_forests_pa.push_back( a_nbPairs[i] );
		}
		forests.push_back(tmp_forests_id);
		forests.push_back(v_ra__);
		forests.push_back(v_de__);
		forests.push_back(tmp_forests_pa);

		// Save the list of pairs
		pathToSave = pathToSave__;
		pathToSave += "xi_A_delta_delta2_list_pairs_";
		pathToSave += prefix1;
		pathToSave += ".txt";
		std::cout << "\n  " << pathToSave << std::endl;

		std::ofstream fFile;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);

		for (unsigned int i=0; i<nbForest_; i++) {
			fFile << v_idx__[i];
			fFile << " " << 0;
			fFile << " " << v_ra__[i];
			fFile << " " << v_de__[i];
			fFile << " " << a_nbPairs[i];
			fFile << std::endl;
		}
		fFile.close();
		

		/// find the index of each forest among the C_NBSUBSAMPLES sub-samples
		LymanForest* lymanForestObject = new LymanForest(forests, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		pathToSave = pathToSave__;
		pathToSave += "xi_A_delta_delta2_map_";
		pathToSave += prefix1;
		pathToSave += ".txt";
		lymanForestObject->SaveRegionMap(pathToSave);
		delete lymanForestObject;
	}


	return;
}
void Correlation::xi_A_delta_delta2_lambda(void) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta2_lambda ------\n" << std::endl;

	///// Load forest 1
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	v_r__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_nb__.clear();
	///// Load forest 2
	loadDataDelta2();
	v_rDelta2__.clear();
	v_lRFDelta2__.clear();

	///// Rapport l1/l2
	const double minl1l2 = 0.89;
	const double maxl1l2 = 1.11;
	double binSizel1l2      = 1.e-03;
	const unsigned int nbBin = int( (maxl1l2-minl1l2)/binSizel1l2);
	std::cout << "  nbBin = " << nbBin << std::endl;
	binSizel1l2 = (maxl1l2-minl1l2)/nbBin;
	const double inv_binSizel1l2 = 1./binSizel1l2;
	std::cout << "  nbBin = " << nbBin<<std::endl;

	///// Theta
	const double minTheta = 0.;
	const double maxTheta = 0.003;
	const double cosThetaMax = cos(maxTheta);
	double binSizeTheta   = 1.e-04;
	const unsigned int nbBinM = int( (maxTheta-minTheta)/binSizeTheta);
	binSizeTheta = (maxTheta-minTheta)/nbBinM;
	const double inv_binSizeTheta = 1./binSizeTheta;
	std::cout << "  nbBinM = " << nbBinM << std::endl;


	///// Arrays for data
	double dataMu[nbBin][nbBinM][7];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;	
	
	for (unsigned int f1=0; f1<nbForest_; f1++) {
		const unsigned int nb1      = v_nbPixelDelta1__[f1];
		const double cosDe          = v_CosDe__[f1];
		const double sinDe          = v_SinDe__[f1];
		const double ra1            = v_ra__[f1];
		const double de1            = v_de__[f1];
		const double zz1            = v_zz__[f1];

		for (unsigned int f2=0; f2<nbForest2__; f2++) {

			///// Remove the correlation qso-ownForest
			if (ra1==v_raDelta2__[f2] && de1==v_deDelta2__[f2] && zz1==v_zzDelta2__[f2]) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeDelta2__[f2]*cos(ra1-v_raDelta2__[f2]) + sinDe*v_SinDeDelta2__[f2];
			if (cosTheta<cosThetaMax) continue;
			const double theta = acos(cosTheta);
			const unsigned int idxTheta = int( theta*inv_binSizeTheta );

			///// Constants of forest n°2
			const unsigned int nb2 = v_nbPixelDelta2__[f2];

			for (unsigned int i1=0; i1<nb1; i1++) {

				const double w1 = v_w__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double z1 = v_z__[f1][i1];
				const double l1 = v_lObs__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					///// Same observed wavelength
					const double l2 = v_lObsDelta2__[f2][i2];
					if (l2==l1) continue;
					const double l1Overl2 = l2/l1;
					if ( l1Overl2<minl1l2 || l1Overl2>maxl1l2 ) continue;

					const double w1w2 = w1*v_wDelta2__[f2][i2];
					const double d1d2 = d1*v_dDelta2__[f2][i2];
					const double z1z2 = z1+v_zDelta2__[f2][i2];
					const unsigned int idx  = int( (l1Overl2-minl1l2)*inv_binSizel1l2 );
					
					dataMu[idx][idxTheta][0] += w1w2*d1d2;
					dataMu[idx][idxTheta][1] += w1w2*d1d2*d1d2;
					dataMu[idx][idxTheta][2] += w1w2*l1Overl2;
					dataMu[idx][idxTheta][3] += w1w2*theta;
					dataMu[idx][idxTheta][4] += w1w2*z1z2;
					dataMu[idx][idxTheta][5] += w1w2;
					dataMu[idx][idxTheta][6] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;
	
	double sumOne = 0.;

	///// Mu
	std::ofstream fFile;
	std::string tmp_pathToSave = pathToSave__;
	tmp_pathToSave = pathToSave__;
	tmp_pathToSave += "xi_A_delta_delta2_lambda_Mu_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += "_";
	tmp_pathToSave += forest2__;
	tmp_pathToSave += ".txt";
	std::cout << "\n  " << tmp_pathToSave << std::endl;
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << std::endl;
			sumOne += dataMu[i][j][6];
		}
	}
	fFile.close();

	std::cout << "  number of pairs      = " << sumOne          << std::endl;


	return;
}



// ---------------------------------------------------------------------
//
//		3D QSO - forest correlation
//
// ---------------------------------------------------------------------
void Correlation::xi_delta_QSO(unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO ------\n" << std::endl;

	if (doBootstraps__)  std::cout << "  subsampling N° "     << bootIdx << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__,bootIdx);
	if (nbForest_==0.) return;

	if (doingCoAddSoRemoving__) removeFalseCorrelations();

	///// Empty useless vectors
	v_zz__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	///// Constants:
	///// The space between bins is of 10 Mpc.h^-1
	const double max          = 200.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100;

	const double maxPow2      = max*max;
	const double dist_diag_max_pow2 = 2*max*max;

	if (shuffleForest || shuffleQSO || randomQSO || randomQSOXYZ) {
		unsigned int seed = 0;
		std::srand(42);
		for (unsigned int i=0; i<10000; i++) {
			const unsigned int tmp_seed = rand();
			if (i == bootIdx) {
				seed = tmp_seed;
				break;
			}
		}
		std::cout << "\n\n" << std::endl;
		std::cout << "  seed = " << seed << std::endl;
		std::srand(seed);
	}

	///// Needed to randomize the forests
	if (shuffleForest) {

		///// Vectors with index of forest
		std::vector<unsigned int> randomIdx(nbForest_);
		for (unsigned int i=0; i<nbForest_; i++) {
			randomIdx[i] = i;
		}

		std::cout << "  Shuffle forest " << bootIdx << std::endl;
		bool doLoop = true;
		unsigned int nbLoop = 0;
		while (doLoop) {
			doLoop = false;
			nbLoop ++;
			std::cout << "  nbLoop = " << nbLoop << std::endl;
			std::random_shuffle ( randomIdx.begin(), randomIdx.end() );
			for (unsigned int i=0; i<nbForest_; i++) {
				if (i==randomIdx[i]) {
					std::cout << nbLoop << " " << i << std::endl;
					doLoop = true;
					break;
				}
			}
		}

		///// Copy data
		std::vector<double> tmp_cosDe(v_CosDe__);
		std::vector<double> tmp_SinDe(v_SinDe__);
		std::vector<double> tmp_ra(v_ra__);
		std::vector<double> tmp_de(v_de__);

		///// Put the new data
		for (unsigned int i=0; i<nbForest_; i++) {
			const unsigned int ii = randomIdx[i];
			v_CosDe__[i] = tmp_cosDe[ii];
			v_SinDe__[i] = tmp_SinDe[ii];
			v_ra__[i]    = tmp_ra[ii];
			v_de__[i]    = tmp_de[ii];
		}
	}
	///// Needed to randomize the QSO
	if (shuffleQSO) {

		///// Vectors with index of QSO
		std::vector<unsigned int> randomIdx(nbQ1__);
		for (unsigned int i=0; i<nbQ1__; i++) {
			randomIdx[i] = i;
		}

		std::cout << "  Shuffle QSO " << bootIdx << std::endl;
		bool doLoop = true;
		unsigned int nbLoop = 0;
		while (doLoop) {
			doLoop = false;
			nbLoop ++;
			std::cout << "  nbLoop = " << nbLoop << std::endl;
			std::random_shuffle ( randomIdx.begin(), randomIdx.end() );
			for (unsigned int i=0; i<nbQ1__; i++) {
				if (i==randomIdx[i]) {
					std::cout << nbLoop << " " << i << std::endl;
					doLoop = true;
					break;
				}
			}
		}

		///// Copy data
		std::vector<double> tmp_zz(v_zzQ1__);
		std::vector<double> tmp_r(v_rrQ1__);
		///// Put the new data
		for (unsigned int i=0; i<nbQ1__; i++) {
			const unsigned int ii = randomIdx[i];
			v_zzQ1__[i] = tmp_zz[ii];
			v_rrQ1__[i] = tmp_r[ii];
		}
	}
	///// Needed to randomize the QSO
	if (randomQSO || randomQSOXYZ) {

		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "\n\n  Random QSO " << bootIdx << std::endl;

		///// Get the edge of the sky
		double raMin = *std::min_element(v_raQ1__.begin(), v_raQ1__.end());
		double raMax = *std::max_element(v_raQ1__.begin(), v_raQ1__.end());
		double deMin = *std::min_element(v_deQ1__.begin(), v_deQ1__.end());
		double deMax = *std::max_element(v_deQ1__.begin(), v_deQ1__.end());
		const double zzMin = *std::min_element(v_zzQ1__.begin(), v_zzQ1__.end());
                const double zzMax = *std::max_element(v_zzQ1__.begin(), v_zzQ1__.end());
		raMin = raMin*0.99;
		raMax = std::min(raMax*1.01,2.*M_PI);
		deMin = std::max(deMin*1.01,-M_PI);
		deMax = std::min(deMax*1.01,M_PI);
		std::cout << "  ra_min = " << raMin << " ra_max = " << raMax << std::endl;
		std::cout << "  de_min = " << deMin << " de_max = " << deMax << std::endl;
		std::cout << "  z_min  = " << zzMin << " z_max  = " << zzMax << std::endl;
		const double coefRA = 1./(raMax-raMin);
		const double coefDE = 1./(deMax-deMin);
		const double coefzz = 1./(zzMax-zzMin);
		for (unsigned int i=0; i<nbQ1__; i++) {
			const double ra = (double)rand()/(double)(RAND_MAX*coefRA) +raMin;
			const double de = (double)rand()/(double)(RAND_MAX*coefDE) +deMin;
			v_raQ1__[i]    = ra;
			v_deQ1__[i]    = de;
			v_CosDeQ1__[i] = cos(de);
			v_SinDeQ1__[i] = sin(de);
			if (randomQSOXYZ) {
				const double zz = (double)rand()/(double)(RAND_MAX*coefzz) +zzMin;
				v_zzQ1__[i] = zz;
				v_rrQ1__[i] = hConvertRedshDist->Interpolate(zz);
			}
		}
		delete hConvertRedshDist;
	}
	if (randomQSO_from_file__) {

		pathQ1__  = path_to_rand_QSO_cat__;
		pathQ1__ += "cat_";
		pathQ1__ += strBootIdx;
		pathQ1__ += ".fits";

		nbQ1__ = 0;
		v_raQ1__.clear();
		v_deQ1__.clear();
		v_zzQ1__.clear();
		v_rrQ1__.clear();
		v_CosDeQ1__.clear();
		v_SinDeQ1__.clear();
		loadDataQ1();
		if (nbQ1__ == 0) return;
	}

	///// get an array for nb of pairs for the forest
	double a_nbPairs[nbForest_];
	for (unsigned int i=0; i<nbForest_; i++) {
		a_nbPairs[i] = 0.;
	}

	///// Arrays for data
	double data2D[nbBinX][nbBinY][9];
	double dataMu[nbBin][nbBinM][9];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<9; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<9; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;


//TH1D* h1 = new TH1D("h1","",70,1.70,5.80);
//double sum[3] = {0.};

	for (unsigned int f=0; f<nbForest_; f++) {

		double nbPairs = 0.;
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			//if (fabs(ra-v_raQ1__[q])<1.e-9 && fabs(dec-v_deQ1__[q])<1.e-9 ) continue;
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;  // && fabs(z-v_zzQ1__[q])<0.001

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= dist_diag_max_pow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;
			const double zQSO = v_zzQ1__[q];

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					else continue;
				}
				
				const double distTotPow2 = distP*distP + distT*distT;
				const double w           = v_w__[f][ii];
				const double d           = v_d__[f][ii];
				const double wd          = w*d;
				const double wdd         = wd*d;
				const double w_res_lRF   = w*v_residual_delta_vs_lRF__[f][ii];
				const double w_res_lObs  = w*v_residual_delta_vs_lObs__[f][ii];

				if (distTotPow2 < maxPow2) {
					const double distTot    = sqrt(distTotPow2);
					const double mu         = distP/distTot;
					const unsigned int idx  = int(distTot);
					const unsigned int idxM = int((mu+1.)*50.);
	
					dataMu[idx][idxM][0] += wd;
					dataMu[idx][idxM][1] += wdd;
					dataMu[idx][idxM][2] += w*distTot;
					dataMu[idx][idxM][3] += w*mu;
					dataMu[idx][idxM][4] += w*(zQSO+v_z__[f][ii]);
					dataMu[idx][idxM][5] += w;
					dataMu[idx][idxM][6] ++;
					dataMu[idx][idxM][7] += w_res_lRF;
					dataMu[idx][idxM][8] += w_res_lObs;
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rPerpBinIdx = int(distT);
				const unsigned int rParralBinIdx = int(distP+max);

				data2D[rPerpBinIdx][rParralBinIdx][0] += wd;
				data2D[rPerpBinIdx][rParralBinIdx][1] += wdd;
				data2D[rPerpBinIdx][rParralBinIdx][2] += w*distT;
				data2D[rPerpBinIdx][rParralBinIdx][3] += w*distP;
				data2D[rPerpBinIdx][rParralBinIdx][4] += w*(zQSO+v_z__[f][ii]);
				data2D[rPerpBinIdx][rParralBinIdx][5] += w;
				data2D[rPerpBinIdx][rParralBinIdx][6] ++;
				data2D[rPerpBinIdx][rParralBinIdx][7] += w_res_lRF;
				data2D[rPerpBinIdx][rParralBinIdx][8] += w_res_lObs;

//h1->Fill( (zQSO+v_z__[f][ii])/2., w );
//sum[0] += w*(zQSO+v_z__[f][ii])/2.;
//sum[1] += w;
//sum[2] ++;		
				///// Get the number of pairs
				nbPairs += w;
			}
		}
	
		///// Put the number of pairs comming with this forest
		a_nbPairs[f] = nbPairs;
	}
/*
	//R_plot1D(h1);
	std::ofstream fFile2;
	fFile2.open("histo_pairs_lya_qso.txt");
	fFile2 << std::scientific;
	fFile2.precision(std::numeric_limits<double>::digits10);
	for (unsigned int i=0; i<70; i++) {
			fFile2 << i << " " << h1->GetBinCenter(i+1) << " " << h1->GetBinContent(i+1) << std::endl;
	}
	fFile2 << 0 << " " << 0 << " " << sum[0]/sum[1] << std::endl;
	fFile2 << 0 << " " << 0 << " " << sum[1] << std::endl;
	fFile2 << 0 << " " << 0 << " " << sum[2] << std::endl;
	fFile2.close();
	//return;
*/
	//// Look if there were some pairs
	bool empty = true;
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			if (data2D[i][j][0]>0.) {
				empty = false;
				break;
			}
		}
		if (!empty) break;
	}
	if (empty) {
		std::cout << "  No pairs" << std::endl;
		return;
	}

	std::cout << "\n  Saving\n" << std::endl;


	//// Set the prefix of different forest and QSOs
	std::string prefix1 = forest__;
	prefix1 += "_";
	prefix1 += QSO__;

	//// Set the prefix for different type of runs
	std::string prefix = "_";
	if (doBootstraps__)           prefix += "subsampling";
	if (shuffleForest)            prefix += "shuffleForest";
	if (shuffleQSO)               prefix += "shuffleQSO";
	if (randomQSO)                prefix += "randomQSO";
	if (randomQSOXYZ)             prefix += "randomQSOXYZ";
	if (randomQSO_from_file__)    prefix += "randomQSOFromFile";
	prefix += "_";
	prefix += strBootIdx;

	std::ofstream fFile;
	std::string pathToSave;
	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;


	///// Save the 2D cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO || randomQSOXYZ || randomQSO_from_file__) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	sumOne = 0.;
	sumZZZ = 0.;
	sumWeg = 0.;
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {

			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4]/2.;
			fFile << " " << data2D[i][j][5];
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << " " << data2D[i][j][8];
			fFile << std::endl;

			sumZZZ += data2D[i][j][4]/2.;
			sumWeg += data2D[i][j][5];
			sumOne += data2D[i][j][6];
		}
	}
	fFile.close();

	std::cout << "  < z >                = " << sumZZZ/sumWeg    << " +- " << 0. << std::endl;
	std::cout << "  number of pairs      = " << sumOne           << std::endl;
	std::cout << "  sum of weight        = " << sumWeg           << std::endl;
	std::cout << "  < pairs per forest > = " << sumOne/nbForest_ << std::endl;
	std::cout << "  < pairs per QSO >    = " << sumOne/nbQ1__    << std::endl;


	///// Save the Mu cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_Mu_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO || randomQSOXYZ || randomQSO_from_file__) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << " " << dataMu[i][j][8];
			fFile << std::endl;
		}
	}
	fFile.close();

	
	if (!doBootstraps__ && !shuffleForest && !shuffleQSO && !randomQSO && !randomQSOXYZ && !randomQSO_from_file__) {

		std::vector< std::vector< double > > forests;
		std::vector< double > tmp_forests_id;
		std::vector< double > tmp_forests_pa;
		
		for (unsigned int i=0; i<nbForest_; i++) {
			tmp_forests_id.push_back( v_idx__[i] );
			tmp_forests_pa.push_back( a_nbPairs[i] );
		}
		forests.push_back(tmp_forests_id);
		forests.push_back(v_ra__);
		forests.push_back(v_de__);
		forests.push_back(tmp_forests_pa);


		// Save the list of pairs
		pathToSave = pathToSave__;
		pathToSave += "xi_delta_QSO_list_pairs_";
		pathToSave += prefix1;
		pathToSave += ".txt";
		std::cout << "\n  " << pathToSave << std::endl;

		std::ofstream fFile;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);

		for (unsigned int i=0; i<nbForest_; i++) {
			fFile << v_idx__[i];
			fFile << " " << 0;
			fFile << " " << v_ra__[i];
			fFile << " " << v_de__[i];
			fFile << " " << a_nbPairs[i];
			fFile << std::endl;
		}
		fFile.close();
		

		//// find the index of each forest among the C_NBSUBSAMPLES sub-samples
		LymanForest* lymanForestObject = new LymanForest(forests, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		pathToSave = pathToSave__;
		pathToSave += "xi_delta_QSO_map_";
		pathToSave += prefix1;
		pathToSave += ".txt";
		lymanForestObject->SaveRegionMap(pathToSave);
		delete lymanForestObject;
	}

	return;
}
void Correlation::xi_delta_QSO_lambda(unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_lambda ------\n" << std::endl;

	if (doBootstraps__)  std::cout << "  subsampling N° "     << bootIdx << std::endl;
	if (shuffleForest) std::cout << "  shuffleForest seed " << bootIdx*10 << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__,bootIdx);
	if (nbForest_==0.) return;

	if (doBootstraps__) removeFalseCorrelations();
	v_zz__.clear();
	v_r__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_nb__.clear();

	///// Set lambda_Obs QSO
	std::vector<double> v_lObsQSO(nbQ1__,0.);
	for (unsigned int q=0; q<nbQ1__; q++) {
		v_lObsQSO[q] = 1./( (v_zzQ1__[q]+1.)*lambdaRFLine__);
	}

	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	///// Rapport l1/l2
	const double minValue = 0.40;
	const double maxValue = 2.30;
	double binSize      = 1.e-03;
	const unsigned int nbBin = int( (maxValue-minValue)/binSize);
	std::cout << "  nbBin = " << nbBin << std::endl;
	binSize = (maxValue-minValue)/nbBin;
	const double inv_binSize = 1./binSize;
	std::cout << "  nbBin = " << nbBin<<std::endl;

	///// Theta
	const double minTheta = 0.;
	const double maxTheta = 0.003;
	const double cosThetaMax = cos(maxTheta);
	double binSizeTheta   = 1.e-04;
	const unsigned int nbBinM = int( (maxTheta-minTheta)/binSizeTheta);
	binSizeTheta = (maxTheta-minTheta)/nbBinM;
	const double inv_binSizeTheta = 1./binSizeTheta;
	std::cout << "  nbBinM = " << nbBinM << std::endl;

	///// Arrays for data
	double dataMu[nbBin][nbBinM][7];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];
		const double firstPixelLambda = v_lObs__[f][0];
		const double lastPixelLambda  = v_lObs__[f][nbPixel-1];
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];
			if ( cosTheta<cosThetaMax ) continue;
			const double theta = acos( cosTheta );
			const unsigned int idxTheta = int( theta*inv_binSizeTheta  );

			///// Transvers distance between the qso and the lya
			const double zQSO = v_zzQ1__[q];
			const double inv_lObsQSO = v_lObsQSO[q];

			// Tests if pixels are in front of QSO
			if (firstPixelLambda*inv_lObsQSO>maxValue) continue;
			else if (lastPixelLambda*inv_lObsQSO<minValue) continue;

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double l1Overl2 = v_lObs__[f][i]*inv_lObsQSO;
				if (l1Overl2<minValue || l1Overl2>maxValue) continue;
				const unsigned int idx  = int( (l1Overl2-minValue)*inv_binSize  );
			
				const double w   = v_w__[f][i];
				const double d   = v_d__[f][i];
				const double wd  = w*d;
				const double wdd = wd*d;
	
				dataMu[idx][idxTheta][0] += wd;
				dataMu[idx][idxTheta][1] += wdd;
				dataMu[idx][idxTheta][2] += w*l1Overl2;
				dataMu[idx][idxTheta][3] += w*theta;
				dataMu[idx][idxTheta][4] += w*(zQSO+v_z__[f][i]);
				dataMu[idx][idxTheta][5] += w;
				dataMu[idx][idxTheta][6] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave;

	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	std::cout << "  < z >                = " << sumZZZ/sumWeg    << " +- " << 0. << std::endl;
	std::cout << "  number of pairs      = " << sumOne          << std::endl;
	std::cout << "  < pairs per forest > = " << sumOne/nbForest_ << std::endl;
	std::cout << "  < pairs per QSO >    = " << sumOne/nbQ1__    << std::endl;

	///// Save the 2D cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_lambda_Mu_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	///// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << std::endl;
		}
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_lambda_same_LOS(unsigned int bootIdx /*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_lambda_same_LOS ------\n" << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__,bootIdx);
	if (nbForest_==0.) return;

	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_r__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_nb__.clear();

	///// Set lambda_Obs QSO
	std::vector<double> v_lObsQSO(nbQ1__,0.);
	for (unsigned int q=0; q<nbQ1__; q++) {
		v_lObsQSO[q] = 1./( (v_zzQ1__[q]+1.)*lambdaRFLine__);
	}

	///// Rapport l1/l2
	const double minValue = 0.40;
	const double maxValue = 2.30;
	double binSize      = 1.e-03;
	const unsigned int nbBin = int( (maxValue-minValue)/binSize);
	std::cout << "  nbBin = " << nbBin << std::endl;
	binSize = (maxValue-minValue)/nbBin;
	const double inv_binSize = 1./binSize;
	std::cout << "  nbBin = " << nbBin<<std::endl;

	///// Arrays for data
	double dataMu[nbBin][6];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int k=0; k<6; k++) {
			dataMu[i][k] = 0.;
		}
	}

	long unsigned int nb_corr = 0;
	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double ra            = v_ra__[f];
		const double de            = v_de__[f];
		const double firstPixelLambda = v_lObs__[f][0];
		const double lastPixelLambda  = v_lObs__[f][nbPixel-1];
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove qso not on the LOS
			if (fabs(ra-v_raQ1__[q])>=C_AUTOCORRCRIT || fabs(de-v_deQ1__[q])>=C_AUTOCORRCRIT ) continue;
			nb_corr ++;

			/////
			const double zQSO = v_zzQ1__[q];
			const double inv_lObsQSO = v_lObsQSO[q];

			// Tests if pixels are in front of QSO
			if (firstPixelLambda*inv_lObsQSO>maxValue) {
				std::cout << firstPixelLambda << " " << inv_lObsQSO << std::endl;
				return;
				continue;
			}
			else if (lastPixelLambda*inv_lObsQSO<minValue) {
				std::cout << lastPixelLambda << " " << inv_lObsQSO << std::endl;
				continue;
			}

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double l1Overl2 = v_lObs__[f][i]*inv_lObsQSO;
				if (l1Overl2<minValue || l1Overl2>maxValue) {
					std::cout << l1Overl2 << std::endl;
					continue;
				}
				const unsigned int idx  = int( (l1Overl2-minValue)*inv_binSize  );
			
				const double w   = v_w__[f][i];
				const double d   = v_d__[f][i];
				const double wd  = w*d;
				const double wdd = wd*d;
	
				dataMu[idx][0] += wd;
				dataMu[idx][1] += wdd;
				dataMu[idx][2] += w*l1Overl2;
				dataMu[idx][3] += w*(zQSO+v_z__[f][i]);
				dataMu[idx][4] += w;
				dataMu[idx][5] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;
	std::cout << "  nb corr = " << nb_corr << std::endl;

	long double sumOne = 0.;
	for (unsigned int i=0; i<nbBin; i++) {
		sumOne += dataMu[i][5];
	}
	if (sumOne==0.) return;

	std::ofstream fFile;
	std::string pathToSave;

	///// Save the 2D cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_lambda_same_LOS_Mu_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	///// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		fFile << dataMu[i][0];
		fFile << " " << dataMu[i][1];
		fFile << " " << dataMu[i][2];
		fFile << " " << dataMu[i][3]/2.;
		fFile << " " << dataMu[i][4];
		fFile << " " << dataMu[i][5];
		fFile << std::endl;
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_distortionMatrix(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_distortionMatrix ------\n" << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_d__.clear();
	v_z__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	v_idx__.clear();

	///// Set usefull vectors
	std::vector<double> v_invSumWeight(nbForest_,0.);
	std::vector< std::vector< double > > v_varLambda(v_lRF__);

	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nbPixel = v_nbPixelDelta1__[f];

		long double sumWeight  = 0.;
		long double meanLambda = 0.;
		long double stdLambda  = 0.;

		///// Loops over all pixels of the forest
		for (unsigned int i=0; i<nbPixel; i++) {

			const long double w = v_w__[f][i];
			const long double l = v_lRF__[f][i];

			sumWeight  += w;
			meanLambda += w*l;
			stdLambda  += w*l*l;
		}

		if (sumWeight==0.) std::cout << "   ERROR:: sumWeight = 0. " << std::endl;
		meanLambda        /= sumWeight;
		stdLambda          = stdLambda/sumWeight - meanLambda*meanLambda;
		if (stdLambda==0.) std::cout << "   ERROR:: stdLambda = 0. " << std::endl;
		stdLambda          = 1./sqrt(stdLambda);
		v_invSumWeight[f]  = 1./sumWeight;

		for (unsigned int i=0; i<nbPixel; i++) {
			v_varLambda[f][i] = (v_lRF__[f][i]-meanLambda)*stdLambda;
		}
	}

	//// Clear usless vector
	v_lRF__.clear();

	///// Constants:
	///// The space between bins is of 10 Mpc.h^-1
	const double max           = 200.;
	const double binSize       = 4.;
	const unsigned int nbBin   = int(max/binSize);
	const unsigned int nbBinX  = nbBin;
	const unsigned int nbBinY  = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double maxPow2       = max*max;
	const double fromValToIdx  = nbBin/max;
	const double dist_diag_max_pow2 = 2*max*max;
	
	///// Arrays for distortion matrix
	double weight2D[nbBin2D];
	double** data2DMatrix = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		weight2D[i] = 0.;
		data2DMatrix[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			data2DMatrix[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];

		const double invSumWeight  = v_invSumWeight[f];

		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= dist_diag_max_pow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;
			const double zQSO = v_zzQ1__[q];

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);



			///// get the weight and mean lambda
			double xValue[5000]  = {0.};
			double xlValue[5000] = {0.};
			bool binTouched[5000] = {false};

			///// Index of the bin
			std::vector< unsigned int > binIdx;
			std::vector< unsigned int > fromBinsToPixels;

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					else continue;
				}

				const double w    = v_w__[f][ii];
				const double val0 = w*invSumWeight;
	
				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rPerpBinIdx = int(distT*fromValToIdx);
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );
				xValue[globalBin]  += val0;
				xlValue[globalBin] += val0*v_varLambda[f][ii];

				///// Fill array of weights
				weight2D[globalBin] += w;

				///// Keep values for the distortion matrix
				binIdx.push_back(globalBin);
				fromBinsToPixels.push_back(ii);

				binTouched[globalBin] = true;
			}

			///// Number of pixels with pairs 
			const unsigned int nbPixelsWithPairs = binIdx.size();

			///// Get a vector of not empty bins
			std::vector< unsigned int> bins_touched;
			for (unsigned int i=0; i<5000; i++) {
				if (binTouched[i]) bins_touched.push_back(i);
			}
			const unsigned int nbBinsTouched = bins_touched.size();

			///// Loops over all pixels of the forest (Fill the distortion matrix)
			for (unsigned int i=0; i<nbPixelsWithPairs; i++) {

				const unsigned int globalBin1 = binIdx[i];
				const unsigned int pixelIdx   = fromBinsToPixels[i];

				const double w    = v_w__[f][pixelIdx];
				const double val1 = v_varLambda[f][pixelIdx];
				
				///// Fill the distortion matrix
				data2DMatrix[globalBin1][globalBin1] += w;

				for (unsigned int j=0; j<nbBinsTouched; j++) {
					const unsigned int globalBin2 = bins_touched[j];
					data2DMatrix[globalBin1][globalBin2] -= w*( xValue[globalBin2] + val1*xlValue[globalBin2] );
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;	

	std::ofstream fFile;

	///// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_distortionMatrix_2D_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			double value = 0.;
			if (weight2D[i]!=0.) value = data2DMatrix[i][j]/weight2D[i];
			fFile << value << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_distortionMatrix_1D(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_distortionMatrix_1D ------\n" << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_d__.clear();
	v_z__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	v_idx__.clear();

	///// Set usefull vectors
	std::vector<double> v_invSumWeight(nbForest_,0.);
	std::vector< std::vector< double > > v_varLambda(v_lRF__);

	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nbPixel = v_nbPixelDelta1__[f];

		double sumWeight  = 0.;
		double meanLambda = 0.;
		double stdLambda  = 0.;

		///// Loops over all pixels of the forest
		for (unsigned int i=0; i<nbPixel; i++) {

			const double w = v_w__[f][i];
			const double l = v_lRF__[f][i];

			sumWeight  += w;
			meanLambda += w*l;
			stdLambda  += w*l*l;
		}

		meanLambda        /= sumWeight;
		stdLambda          = 1./sqrt(stdLambda/sumWeight - meanLambda*meanLambda);
		v_invSumWeight[f]  = 1./sumWeight;

		for (unsigned int i=0; i<nbPixel; i++) {
			v_varLambda[f][i] = (v_lRF__[f][i]-meanLambda)*stdLambda;
		}
	}

	//// Clear usless vector
	v_lRF__.clear();

	///// Constants:
	///// The space between bins is of 10 Mpc.h^-1
	const double max           = 200.;
	const double binSize       = 4.;
	const unsigned int nbBin   = int(max/binSize);
	const double maxPow2       = max*max;
	const double fromValToIdx  = nbBin/max;

	///// Arrays for distortion matrix
	double weight[nbBin];
	double** dataMatrix = new double*[nbBin];
	for (unsigned int i=0; i<nbBin; i++) {
		weight[i] = 0.;
		dataMatrix[i] = new double[nbBin];
		for (unsigned int j=0; j<nbBin; j++) {
			dataMatrix[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];

		const double invSumWeight  = v_invSumWeight[f];

		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= maxPow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;
			const double zQSO = v_zzQ1__[q];

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// get the weight and mean lambda
			double xValue[50]  = {0.};
			double xlValue[50] = {0.};
			bool binTouched[50] = {false};

			///// Index of the bin
			std::vector< unsigned int > binIdx;
			std::vector< unsigned int > fromBinsToPixels;

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					else continue;
				}

				const double distTotPow2 = distT*distT + distP*distP;
				if (distTotPow2 >= maxPow2) continue;

				const double w    = v_w__[f][ii];
				const double val0 = w*invSumWeight;
	
				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int globalBin = int(sqrt(distTotPow2)*fromValToIdx);
				xValue[globalBin]  += val0;
				xlValue[globalBin] += val0*v_varLambda[f][ii];

				///// Fill array of weights
				weight[globalBin] += w;

				///// Keep values for the distortion matrix
				binIdx.push_back(globalBin);
				fromBinsToPixels.push_back(ii);

				binTouched[globalBin] = true;
			}

			///// Number of pixels with pairs 
			const unsigned int nbPixelsWithPairs = binIdx.size();

			///// Get a vector of not empty bins
			std::vector< unsigned int> bins_touched;
			for (unsigned int i=0; i<50; i++) {
				if (binTouched[i]) bins_touched.push_back(i);
			}
			const unsigned int nbBinsTouched = bins_touched.size();			

			///// Loops over all pixels of the forest (Fill the distortion matrix)
			for (unsigned int i=0; i<nbPixelsWithPairs; i++) {

				const unsigned int globalBin1 = binIdx[i];
				const unsigned int pixelIdx   = fromBinsToPixels[i];

				const double w    = v_w__[f][pixelIdx];
				const double val1 = v_varLambda[f][pixelIdx];
				
				///// Fill the distortion matrix
				dataMatrix[globalBin1][globalBin1] += w;

				for (unsigned int j=0; j<nbBinsTouched; j++) {
					const unsigned int globalBin2 = bins_touched[j];
					dataMatrix[globalBin1][globalBin2] -= w*( xValue[globalBin2] + val1*xlValue[globalBin2] );
				}
			}
		}
	}
	
	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;

	///// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_distortionMatrix_1D_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	//// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			double value = 0.;
			if (weight[i]!=0.) value = dataMatrix[i][j]/weight[i];
			fFile << value << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_Metals_Models(double lambdaFrMetal, std::string lambdaFrMetalName) {

	// 21,24,25,26,27
	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Metals_Models ------" << std::endl;
	std::string command = "  ";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	std::cout << "  line = " << lambdaFrMetalName << " : " << lambdaFrMetal << std::endl;

	//// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	//// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	/// Get an arrey to get te growth factor
	double from_z_to_growth_factor_pow_2[C_NBBINREDSH] = {0.};
	const double g0    = cosmo->get_growth_factor(0.);
	const double growth_factor_step_size = C_ZEXTREMABINCONVERT1/C_NBBINREDSH;
	const double inverse_growth_factor_step_size = 1./growth_factor_step_size;
	for (unsigned int i=0; i<C_NBBINREDSH; i++) {
		const double ZZZ = i*growth_factor_step_size;
		from_z_to_growth_factor_pow_2[i] = (cosmo->get_growth_factor(ZZZ)/g0)*(cosmo->get_growth_factor(ZZZ)/g0);
	}
	
	//// Get the distance if it was this metal
	std::vector< std::vector< double > > v_r_metal(v_r__);
	std::vector< std::vector< double > > v_z_metal(v_r__);	
	for (unsigned int i=0; i<v_r_metal.size(); i++) {
		for (unsigned int j=0; j<v_r_metal[i].size(); j++) {
			v_z_metal[i][j] = v_lObs__[i][j]/lambdaFrMetal-1.;
			v_r_metal[i][j] = hConvertRedshDist->Interpolate( v_z_metal[i][j] );
		}
	}

	//// Empty useless vectors
	delete cosmo;
	delete hConvertRedshDist;
	v_zz__.clear();
	v_d__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();


	//// Get the monopole, quadrupol, exadecapol
	std::string pathToLoad = PATHTOWORK;
	pathToLoad += C_PATHTOCAMB;
	std::ifstream fileData(pathToLoad.c_str());
	std::vector< double > data_x;
	std::vector< double > data_xi0;
	std::vector< double > data_xi2;
	std::vector< double > data_xi4;
	std::cout << pathToLoad << std::endl;
	while (fileData) {
		double x;
		double xi0;
		double xi2;
		double xi4;
		fileData>>x>>xi0>>xi2>>xi4;
		if (fileData==0) break;
		
		data_x.push_back(x);
		data_xi0.push_back(xi0);
		data_xi2.push_back(xi2);
		data_xi4.push_back(xi4);
	}
	
	const double inverse_step_size = 1./data_x[0];
	const unsigned int nbBinCAMB = data_x.size();
	//// Set to zero the last pixel
	const double maxDistCAMB = data_x[nbBinCAMB-1];
	data_x.push_back(data_x[nbBinCAMB-1]+1.);
	data_xi0.push_back(0.);
	data_xi2.push_back(0.);
	data_xi4.push_back(0.);

	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max          = 200.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100;

	const double maxPow2      = max*max;
	const double dist_diag_max_pow2 = 2*max*max;

	//// Arrays for data
	double data2D[nbBinX][nbBinY][11];
	double dataMu[nbBin][nbBinM][11];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<8; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<8; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {
	
		//// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= dist_diag_max_pow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;
			const double zQSO = v_zzQ1__[q];

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					else continue;
				}

				const double distP_Metal = v_r_metal[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT_Metal = v_r_metal[f][ii]*sinHalfTheta + distPerpQsoLya;
				const double distTotMetal = sqrt( distT_Metal*distT_Metal + distP_Metal*distP_Metal );

				const double w   = v_w__[f][ii];
				unsigned int idxBinCAMB = nbBinCAMB;
				if (distTotMetal<maxDistCAMB) idxBinCAMB = int(distTotMetal*inverse_step_size);
				const double growth_factor_pow_2 = from_z_to_growth_factor_pow_2[ int(0.5*(zQSO+v_z_metal[f][ii])*inverse_growth_factor_step_size) ];
				const double wxi0 = w*growth_factor_pow_2*data_xi0[idxBinCAMB];
				const double wxi2 = w*growth_factor_pow_2*data_xi2[idxBinCAMB];
				const double wxi4 = w*growth_factor_pow_2*data_xi4[idxBinCAMB];
				const double wz   = w*(zQSO+v_z_metal[f][ii]);

				const double distTotPow2 = distT*distT + distP*distP;
				
				if (distTot < max) {
					const double distTot    = sqrt(distTotPow2);
					const double mu         = distP/distTot;
					const unsigned int idx  = int(distTot);
					const unsigned int idxM = int((mu+1.)*50.);
	
					dataMu[idx][idxM][0] += wxi0;
					dataMu[idx][idxM][1] += wxi2;
					dataMu[idx][idxM][2] += wxi4;
					dataMu[idx][idxM][3] += w*distTot;
					dataMu[idx][idxM][4] += w*mu;
					dataMu[idx][idxM][5] += wz;
					dataMu[idx][idxM][6] += w;
					dataMu[idx][idxM][7] ++;
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rPerpBinIdx = int(distT);
				const unsigned int rParralBinIdx = int(distP+max);
				data2D[rPerpBinIdx][rParralBinIdx][0] += wxi0;
				data2D[rPerpBinIdx][rParralBinIdx][1] += wxi2;
				data2D[rPerpBinIdx][rParralBinIdx][2] += wxi4;
				data2D[rPerpBinIdx][rParralBinIdx][3] += w*distT;
				data2D[rPerpBinIdx][rParralBinIdx][4] += w*distP;
				data2D[rPerpBinIdx][rParralBinIdx][5] += wz;
				data2D[rPerpBinIdx][rParralBinIdx][6] += w;
				data2D[rPerpBinIdx][rParralBinIdx][7] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	//// Set the prefix of different forest and QSOs
	std::string prefix1 = forest__;
	prefix1 += "_";
	prefix1 += QSO__;
	prefix1 += "_";
	prefix1 += lambdaFrMetalName;

	std::ofstream fFile;
	std::string pathToSave;


	///// Save the 2D cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_Metals_Models_2D_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {

			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4];
			fFile << " " << data2D[i][j][5]/2.;
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();

	///// Save the Mu cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_Metals_Models_Mu_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4];
			fFile << " " << dataMu[i][j][5]/2.;
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_Wick_T1_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Wick_T1_with_wi1D_array ------\n" << std::endl;


	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_lRF__.clear();
	v_nb__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;
	const double dist_diag_max_pow2 = 2*max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** cross  = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		cross[i]  = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			cross[i][j]  = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	//// ------------------------------------------------------------------
	//// Load correlation delta-QSO
	pathToFile = pathToSave__;
	pathToFile += "xi_delta_QSO_txt_2D_";
	pathToFile += forest__;
	pathToFile += "_";
	pathToFile += QSO__;
	pathToFile += ".txt";
	std::ifstream fileData2(pathToFile.c_str());

	unsigned int nb_lines_in_file = 0;
	double xiCross[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		long double data0,data1,data2,data3;
		fileData2>>data0>>data1>>data2>>data3;
		xiCross[i] = data2;
		nb_lines_in_file ++;
	}
	fileData2.close();
	if (nb_lines_in_file!=nbBin2D) return;


	std::cout << "  Starting " << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;  // && fabs(z-v_zzQ1__[q])<0.001

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= dist_diag_max_pow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					continue;
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rPerpBinIdx = int(distT*fromValToIdx);
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );

				const double w = v_w__[f][ii];

				nb[globalBin] ++;
				weight[globalBin] += w;

				const unsigned int pixel_lObs_index_1 = 1+int(v_lObs__[f][ii]-lambdaObsMin__);
				cov[globalBin][globalBin]    += w*w*xidd1D[pixel_lObs_index_1][pixel_lObs_index_1];
				cross[globalBin][globalBin]  += w*w*xiCross[globalBin]*xiCross[globalBin];
				nb_all[globalBin][globalBin] ++;
                        	we_all[globalBin][globalBin] += w*w;
			}
		}
	}

	const std::string diagram = "T1_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the covariance with delta-delta
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the covaraince terms with the cross
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_crossTerms_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cross[i][j] /= weight[i]*weight[j];
			fFile << cross[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}
void Correlation::xi_delta_QSO_Wick_T12_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Wick_T12_with_wi1D_array ------\n" << std::endl;


	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_lRF__.clear();
	v_nb__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;
        const double dist_diag_max_pow2 = 2*max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** cross  = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		cross[i]  = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			cross[i][j]  = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	//// ------------------------------------------------------------------
	//// Load correlation delta-QSO
	pathToFile = pathToSave__;
	pathToFile += "xi_delta_QSO_txt_2D_";
	pathToFile += forest__;
	pathToFile += "_";
	pathToFile += QSO__;
	pathToFile += ".txt";
	std::ifstream fileData2(pathToFile.c_str());

	unsigned int nb_lines_in_file = 0;
	double xiCross[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		long double data0,data1,data2,data3;
		fileData2>>data0>>data1>>data2>>data3;
		xiCross[i] = data2;
		nb_lines_in_file ++;
	}
	fileData2.close();
	if (nb_lines_in_file!=nbBin2D) return;


	std::cout << "  Starting " << std::endl;


	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			//// For T1 + T2 + T3 + T4
			std::vector<unsigned int> vect_pair_bin;
			std::vector<unsigned int> vect_pixel_index;

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;  // && fabs(z-v_zzQ1__[q])<0.001

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= dist_diag_max_pow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					continue;
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rPerpBinIdx = int(distT*fromValToIdx);
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );
				const double w = v_w__[f][ii];

				nb[globalBin] ++;
				weight[globalBin] += w;
				
				const unsigned int pixel_lObs_index_1 = 1+int(v_lObs__[f][ii]-lambdaObsMin__);
				cov[globalBin][globalBin]    += w*w*xidd1D[pixel_lObs_index_1][pixel_lObs_index_1];
				cross[globalBin][globalBin]  += w*w*xiCross[globalBin]*xiCross[globalBin];
				nb_all[globalBin][globalBin] ++;
                        	we_all[globalBin][globalBin] += w*w;

				//// For T2
				vect_pair_bin.push_back(globalBin);
				vect_pixel_index.push_back(ii);
			}

			//// Calculate the covariance
			const unsigned int nbPairs = vect_pair_bin.size();
			
			for (unsigned int i=0; i<nbPairs; i++) {

				const unsigned int global_bin_1  = vect_pair_bin[i];
				const unsigned int pixel_index_1 = vect_pixel_index[i];
				const double www1 = v_w__[f][pixel_index_1];
				const unsigned int pixel_lObs_index_1 = int(v_lObs__[f][pixel_index_1]-lambdaObsMin__);

				//// Covariance
				for (unsigned int j=0; j<nbPairs; j++) {
					if (i==j) continue;

					const unsigned int global_bin_2  = vect_pair_bin[j];
					const unsigned int pixel_index_2 = vect_pixel_index[j];
					const double www2 = v_w__[f][pixel_index_2];
					unsigned int pixel_lObs_index_2 = int(v_lObs__[f][pixel_index_2]-lambdaObsMin__);

					if (v_lObs__[f][pixel_index_1]>v_lObs__[f][pixel_index_2]) {
						cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2];
					}
					else cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_2+1][pixel_lObs_index_1];
					cross[global_bin_1][global_bin_2] += www1*www2*xiCross[global_bin_1]*xiCross[global_bin_2];
					nb_all[global_bin_1][global_bin_2] ++;
					we_all[global_bin_1][global_bin_2] += www1*www2;
				}
			}
		}
	}


	const std::string diagram = "T12_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the covariance with delta-delta
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the covaraince terms with the cross
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_crossTerms_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cross[i][j] /= weight[i]*weight[j];
			fFile << cross[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}
void Correlation::xi_delta_QSO_Wick_T123_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Wick_T123_with_wi1D_array ------\n" << std::endl;


	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_lRF__.clear();
	v_nb__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;
        const double dist_diag_max_pow2 = 2*max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** cross  = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		cross[i]  = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			cross[i][j]  = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	//// ------------------------------------------------------------------
	//// Load correlation delta-QSO
	pathToFile = pathToSave__;
	pathToFile += "xi_delta_QSO_txt_2D_";
	pathToFile += forest__;
	pathToFile += "_";
	pathToFile += QSO__;
	pathToFile += ".txt";
	std::ifstream fileData2(pathToFile.c_str());

	unsigned int nb_lines_in_file = 0;
	double xiCross[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		long double data0,data1,data2,data3;
		fileData2>>data0>>data1>>data2>>data3;
		xiCross[i] = data2;
		nb_lines_in_file ++;
	}
	fileData2.close();
	if (nb_lines_in_file!=nbBin2D) return;


	std::cout << "  Starting " << std::endl;


	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];

		//// For T1 + T2 + T3 + T4
		std::vector<unsigned int> vect_pair_bin;
		std::vector<unsigned int> vect_pixel_index;
		std::vector<unsigned int> vect_QSO_index;
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;  // && fabs(z-v_zzQ1__[q])<0.001

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= dist_diag_max_pow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					continue;
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rPerpBinIdx = int(distT*fromValToIdx);
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );
				const double w = v_w__[f][ii];

				nb[globalBin] ++;
				weight[globalBin] += w;

				//// For T3
				vect_pair_bin.push_back(globalBin);
				vect_pixel_index.push_back(ii);
				vect_QSO_index.push_back(q);
			}
		}

		//// Calculate the covariance
		const unsigned int nbPairs = vect_pair_bin.size();
		std::cout << f << " " << nbPairs << std::endl;

		for (unsigned int i=0; i<nbPairs; i++) {

			const unsigned int global_bin_1  = vect_pair_bin[i];
			const unsigned int pixel_index_1 = vect_pixel_index[i];
			const unsigned int QSO_index_1   = vect_QSO_index[i];
			const double www1 = v_w__[f][pixel_index_1];
			const unsigned int pixel_lObs_index_1 = int(v_lObs__[f][pixel_index_1]-lambdaObsMin__);

			//// Covariance
			for (unsigned int j=0; j<nbPairs; j++) {
				
				const unsigned int pixel_index_2 = vect_pixel_index[j];
				const unsigned int QSO_index_2   = vect_QSO_index[j];
				if (pixel_index_1!=pixel_index_2 && QSO_index_1!=QSO_index_2) continue;  // remove T4
				const unsigned int global_bin_2  = vect_pair_bin[j];
				const unsigned int pixel_lObs_index_2 = int(v_lObs__[f][pixel_index_2]-lambdaObsMin__);
				const double www2 = v_w__[f][pixel_index_2];
				
				if (pixel_index_1==pixel_index_2) cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2+1];
				else {
					if (v_lObs__[f][pixel_index_1]>v_lObs__[f][pixel_index_2]) {
						cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2];
					}
					else cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_2+1][pixel_lObs_index_1];
				}
				cross[global_bin_1][global_bin_2] += www1*www2*xiCross[global_bin_1]*xiCross[global_bin_2];
				nb_all[global_bin_1][global_bin_2] ++;
				we_all[global_bin_1][global_bin_2] += www1*www2;
			}
		}

	}


	const std::string diagram = "T123_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the covariance with delta-delta
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the covaraince terms with the cross
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_crossTerms_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cross[i][j] /= weight[i]*weight[j];
			fFile << cross[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}
void Correlation::xi_delta_QSO_Wick_T1234_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Wick_T1234_with_wi1D_array ------\n" << std::endl;


	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_lRF__.clear();
	v_nb__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;
        const double dist_diag_max_pow2 = 2*max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** cross  = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		cross[i]  = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			cross[i][j]  = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	//// ------------------------------------------------------------------
	//// Load correlation delta-QSO
	pathToFile = pathToSave__;
	pathToFile += "xi_delta_QSO_txt_2D_";
	pathToFile += forest__;
	pathToFile += "_";
	pathToFile += QSO__;
	pathToFile += ".txt";
	std::ifstream fileData2(pathToFile.c_str());

	unsigned int nb_lines_in_file = 0;
	double xiCross[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		long double data0,data1,data2,data3;
		fileData2>>data0>>data1>>data2>>data3;
		xiCross[i] = data2;
		nb_lines_in_file ++;
	}
	fileData2.close();
	if (nb_lines_in_file!=nbBin2D) return;


	std::cout << "  Starting " << std::endl;


	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];

		//// For T1 + T2 + T3 + T4
		std::vector<unsigned int> vect_pair_bin;
		std::vector<unsigned int> vect_pixel_index;
			
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;  // && fabs(z-v_zzQ1__[q])<0.001

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			if (v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta) >= dist_diag_max_pow2) continue;

			///// Reject QSO with a perpandicular distance too large
			if ( (firstPixel+v_rrQ1__[q])*(firstPixel+v_rrQ1__[q])*0.5*(1.-cosTheta) >= maxPow2 ) continue;

			///// Parrallel distance between the qso and the lya
			const double cosHalfTheta = sqrt( 0.5*(1.+cosTheta) );
			const double distParalQsoLya = v_rrQ1__[q]*cosHalfTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel*cosHalfTheta - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Perpandicular distance between the qso and the lya
			const double sinHalfTheta = sqrt( 0.5*(1.-cosTheta) );
			const double distPerpQsoLya = v_rrQ1__[q]*sinHalfTheta;

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii]*cosHalfTheta - distParalQsoLya;
				const double distT = v_r__[f][ii]*sinHalfTheta + distPerpQsoLya;
				if (fabs(distP) >= max || distT >= max) {
					if ( fabs(distP) >= max && (infPosBool || supPosBool) ) break;
					else if ( distT >= max && infPosBool) break;
					continue;
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rPerpBinIdx = int(distT*fromValToIdx);
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );
				const double w = v_w__[f][ii];

				nb[globalBin] ++;
				weight[globalBin] += w;

				//// For T3
				vect_pair_bin.push_back(globalBin);
				vect_pixel_index.push_back(ii);
			}
		}

		//// Calculate the covariance
		const unsigned int nbPairs = vect_pair_bin.size();
		std::cout << f << " " << nbPairs << std::endl;

		for (unsigned int i=0; i<nbPairs; i++) {

			const unsigned int global_bin_1  = vect_pair_bin[i];
			const unsigned int pixel_index_1 = vect_pixel_index[i];
			const double www1 = v_w__[f][pixel_index_1];
			const unsigned int pixel_lObs_index_1 = int(v_lObs__[f][pixel_index_1]-lambdaObsMin__);

			//// Covariance
			for (unsigned int j=0; j<nbPairs; j++) {
				
				const unsigned int pixel_index_2 = vect_pixel_index[j];
				const unsigned int global_bin_2  = vect_pair_bin[j];
				const unsigned int pixel_lObs_index_2 = int(v_lObs__[f][pixel_index_2]-lambdaObsMin__);
				const double www2 = v_w__[f][pixel_index_2];
				
				if (pixel_index_1==pixel_index_2) cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2+1];
				else {
					if (v_lObs__[f][pixel_index_1]>v_lObs__[f][pixel_index_2]) {
						cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2];
					}
					else cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_2+1][pixel_lObs_index_1];
				}
				nb_all[global_bin_1][global_bin_2] ++;
				we_all[global_bin_1][global_bin_2] += www1*www2;
			}
		}
	}


	const std::string diagram = "T1234_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the covariance with delta-delta
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the covaraince terms with the cross
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_crossTerms_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cross[i][j] /= weight[i]*weight[j];
			fFile << cross[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}
void Correlation::xi_delta_QSO_Wick_1D(unsigned int diagramIdx) {

	
	//// Convert the diagram index to string
	std::stringstream convert;
	convert << diagramIdx;
	std::string strConvert = convert.str();

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Wick_1D ------\n" << std::endl;


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator_2016_05_26/xi_1D_delta_delta_LYA.txt";
	std::ifstream fileData(pathToFile.c_str());

	std::vector<double> tmp_xidd1D;

	long double save0 = 0.;
	long double save1 = 0.;
	long double save2 = 0.;
	long double save3 = 0.;
	long double save4 = 0.;
	long double save5 = 0.;

	while (fileData) {
		fileData>>save0>>save1>>save2>>save3>>save4>>save5;
		if (fileData==0) break;
		if (save5==0.) continue;
		tmp_xidd1D.push_back(save0/save4);
	}
	fileData.close();

	std::cout << "  size xi_delta_delta_1D = " << tmp_xidd1D.size() << std::endl;
	double xidd1D[tmp_xidd1D.size()];
	for (unsigned int i=0; i<tmp_xidd1D.size(); i++) {
		xidd1D[i] = tmp_xidd1D[i];
	}
	tmp_xidd1D.clear();
	//// ------------------------------------------------------------------
	

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_zz__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;

	//// Array with the covariance
	double nb[nbBin];
	double weight[nbBin];
	double** cov    = new double*[nbBin];
	double** nb_all = new double*[nbBin];
	double** we_all = new double*[nbBin];
	for (unsigned int i=0; i<nbBin; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin];
		nb_all[i] = new double[nbBin];
		we_all[i] = new double[nbBin];
		for (unsigned int j=0; j<nbBin; j++) {
			cov[i][j]    = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {
	
		///// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];

		for (unsigned int q=0; q<nbQ1__; q++) {

			//// For T1 + T2 + T3 + T4
			std::vector<unsigned int> vect_pair_bin;
			std::vector<unsigned int> vect_pixel_index;

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			const double distTransQsoLyaPow2 = v_rrQ1__[q]*v_rrQ1__[q]*(1.-cosTheta*cosTheta);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			///// Parrallel distance between the qso and the lya
			const double distParalQsoLya = v_rrQ1__[q]*cosTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// 'true' if first pixel of forest is further than the QSO
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///// 'true' if last pixel of forest is lower than the QSO
			const bool supPosBool = (distParalQsoLastPixel < 0.);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
	
				const double distP = v_r__[f][ii] - distParalQsoLya;				

				///// Look at the position of the Lya forest regarding the qso
				if (fabs(distP) >= max) {
					if (infPosBool || supPosBool) break;
					else continue;
				}
				const double distTotPow2 = distTransQsoLyaPow2 + distP*distP;
				if (distTotPow2 >= maxPow2) continue;

				const unsigned int globalBin = int(sqrt(distTotPow2)*fromValToIdx);
				const double w = v_w__[f][ii];

				//cov[globalBin][globalBin] += w*w*xidd1D[0];
				cov[globalBin][globalBin] += w*w*v_d__[f][ii]*v_d__[f][ii];
				nb[globalBin] ++;
				weight[globalBin] += w;
				nb_all[globalBin][globalBin] ++;
                        	we_all[globalBin][globalBin] += w*w;

				//// For T1 + T2
				//vect_pair_bin.push_back(globalBin);
				//vect_pixel_index.push_back(ii);
			}

			/*
			//// Calculate the covariance
			const unsigned int nbPairs = vect_pair_bin.size();
			
			for (unsigned int i=0; i<nbPairs; i++) {

				const unsigned int global_bin_1  = vect_pair_bin[i];
				const unsigned int pixel_index_1 = vect_pixel_index[i];
				const double www1 = v_w__[f][pixel_index_1];
				//const double rrr1 = v_w__[f][pixel_index_1];
				const double ddd1 = v_d__[f][pixel_index_1];

				//// Covariance
				for (unsigned int j=0; j<nbPairs; j++) {
					if (i==j) continue;

					const unsigned int global_bin_2  = vect_pair_bin[j];
					const unsigned int pixel_index_2 = vect_pixel_index[j];
					const double www2 = v_w__[f][pixel_index_2];
					//const double rrr2 = v_w__[f][pixel_index_2];
					const double ddd2 = v_d__[f][pixel_index_2];
					//const unsigned int idx_corr_1D = 1 + int(fabs(rrr1-rrr2));
					
					//cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[idx_corr_1D];
					cov[global_bin_1][global_bin_2] += www1*www2*ddd1*ddd2;
					nb_all[global_bin_1][global_bin_2] ++;
					we_all[global_bin_1][global_bin_2] += www1*www2;
				}
			}*/
		}
	}





	std::cout << "\n  Saving\n" << std::endl;	

	std::ofstream fFile;

	///// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_1D_Wick";
	pathToSave += "_T1"; // T12  // T1234  // T124 //withDelta
	pathToSave += "_withDelta_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_1D_Wick";
        pathToSave += "_T1"; // T12  // T1234  // T124
	pathToSave += "_withDelta";
	pathToSave += "_number_"; 
	pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
                for (unsigned int j=0; j<nbBin; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_1D_Wick";
        pathToSave += "_T1"; // T12  // T1234  // T124
	pathToSave += "_withDelta";
	pathToSave += "_weight_"; // T12  // T1234  // T124
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin; i++) {
                for (unsigned int j=0; j<nbBin; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	return;
}



// ---------------------------------------------------------------------
//
//		3D QSO 1 - QSO 2 correlation
//
// ---------------------------------------------------------------------
void Correlation::xi_QSO_QSO(unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_QSO_QSO ------\n" << std::endl;

	loadDataQ1();
	if (nbQ1__ == 0) return;

	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	//// Keep the data
	std::vector<double> dataRa(v_raQ1__);
	std::vector<double> dataCoDe(v_CosDeQ1__);
	std::vector<double> dataSiDe(v_SinDeQ1__);
	std::vector<double> dataZZ(v_zzQ1__);
	std::vector<double> dataRR(v_rrQ1__);

	///// Needed to randomize the QSO
	if (randomQSO) {

		std::cout << "  Random QSO " << bootIdx << std::endl;
		std::srand (bootIdx*10);

		///// Get the edge of the sky
		double raMin = *std::min_element(v_raQ1__.begin(), v_raQ1__.end());
		double raMax = *std::max_element(v_raQ1__.begin(), v_raQ1__.end());
		double deMin = *std::min_element(v_deQ1__.begin(), v_deQ1__.end());
		double deMax = *std::max_element(v_deQ1__.begin(), v_deQ1__.end());
		raMin = raMin*0.99;
		raMax = std::min(raMax*1.01,2.*M_PI);
		deMin = std::max(deMin*1.01,-M_PI);
		deMax = std::min(deMax*1.01,M_PI);
		std::cout << "  " << raMin << " " << raMax << " " << deMin << " " << deMax << std::endl;

		const double coefRA = 1./(raMax-raMin);
		const double coefDE = 1./(deMax-deMin);
		for (unsigned int i=0; i<nbQ1__; i++) {
			const float ra = (float)rand()/(float)(RAND_MAX*coefRA) +raMin;
			const float de = (float)rand()/(float)(RAND_MAX*coefDE) +deMin;
			v_raQ1__[i]    = ra;
			v_deQ1__[i]    = de;
			v_CosDeQ1__[i] = cos(de);
			v_SinDeQ1__[i] = sin(de);
		}
	}
	else if (randomQSO_from_file__) {

		pathQ1__  = path_to_rand_QSO_cat__;
		pathQ1__ += "cat_";
		pathQ1__ += strBootIdx;
		pathQ1__ += ".fits";

		v_raQ1__.clear();
		v_deQ1__.clear();
		v_zzQ1__.clear();
		v_rrQ1__.clear();
		v_CosDeQ1__.clear();
		v_SinDeQ1__.clear();
		loadDataQ1();
		if (nbQ1__ == 0) return;
	}
	else {
		dataRa.clear();
		dataCoDe.clear();
		dataSiDe.clear();
		dataZZ.clear();
		dataRR.clear();
	}



	//// Constants:
	const double max          = 200.;
	const unsigned int nbBin  = int(max);
	const double maxPow2      = max*max;
	const unsigned int nbBinM = 50;

	//// Arrays for data
	double dataMu[nbBin][nbBinM][4];
	double data2D[nbBin][nbBin][4];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<4; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<4; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int q1=0; q1<nbQ1__; q1++) {
		
		const double ra1 = v_raQ1__[q1];
		const double cosDe1 = v_CosDeQ1__[q1];
		const double sinDe1 = v_SinDeQ1__[q1];
		const double zz1 = v_zzQ1__[q1];
		const double rr1 = v_rrQ1__[q1];

		for (unsigned int q2=0; q2<q1; q2++) {

			//// Angle between the two directions of thes QSOs
			const double cosTheta = cosDe1*v_CosDeQ1__[q2]*cos(ra1-v_raQ1__[q2]) + sinDe1*v_SinDeQ1__[q2];

			///// Perpendicular distance
			const double distT = (rr1+v_rrQ1__[q2])*sqrt( 0.5*(1.-cosTheta) );
			if (distT >= max) continue;

			///// Parrallel distance
			const double distP = fabs(rr1-v_rrQ1__[q2])*sqrt( 0.5*(1.+cosTheta) );
			if (distP >= max) continue;

			const unsigned int idxPerp = int(distT);
			const unsigned int idxPara = int(distP);

			//// 1D
			const double distTotPow2 = distT*distT + distP*distP;
			if (distTotPow2<maxPow2) {
				const double distTot    = sqrt(distTotPow2);
				const double mu         = distP/distTot;
				const unsigned int idx  = int(distTot);
				const unsigned int idxM = int(50.*mu);

				dataMu[idx][idxM][0] ++;
				dataMu[idx][idxM][1] += distTot;
				dataMu[idx][idxM][2] += mu;
				dataMu[idx][idxM][3] += zz1+v_zzQ1__[q2];
			}

			//// 2D
			data2D[idxPerp][idxPara][0] ++;
			data2D[idxPerp][idxPara][1] += distT;
			data2D[idxPerp][idxPara][2] += distP;
			data2D[idxPerp][idxPara][3] += zz1+v_zzQ1__[q2];
		}
	}
	
	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave;
	std::string prefix = QSO__;
	if (randomQSO || randomQSO_from_file__) {
		prefix += "_RR_";
		prefix += strBootIdx;
	}
	else prefix += "_DD";


	long double sumOne = 0.;
	long double meanZ = 0.;

	//// 2D
	pathToSave = pathToSave__;
	pathToSave += "xi_QSO_QSO_2D_";
	pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	sumOne = 0.;
	meanZ  = 0.;

	//// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3]/2.;
			fFile << std::endl;

			sumOne += data2D[i][j][0];
			meanZ  += data2D[i][j][3]/2.;
		}
	}
	fFile.close();

	//// Get mean redshift of pairs
	std::cout << "  < z >           = " << meanZ/sumOne << " +- " << 0. << std::endl;
	std::cout << "  number of pairs = " << sumOne          << std::endl;


	// Mu
	pathToSave = pathToSave__;
	pathToSave += "xi_QSO_QSO_Mu_";
	pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	//// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3]/2.;
			fFile << std::endl;
		}
	}
	fFile.close();


	//// Doing the correlation random - data
	if (randomQSO || randomQSO_from_file__) {
	
		//// Reset values to zero
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBinM; j++) {
				for (unsigned int k=0; k<4; k++) {
					dataMu[i][j][k] = 0.;
				}
			}
		}
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBin; j++) {
				for (unsigned int k=0; k<4; k++) {
					data2D[i][j][k] = 0.;
				}
			}
		}

		std::cout << "\n  Starting\n" << std::endl;

		for (unsigned int q1=0; q1<nbQ1__; q1++) {

			const double ra1    = dataRa[q1];
			const double cosDe1 = dataCoDe[q1];
			const double sinDe1 = dataSiDe[q1];
			const double zz1    = dataZZ[q1];
			const double rr1    = dataRR[q1];
	
			for (unsigned int q2=0; q2<nbQ1__; q2++) {
	
				//// Angle between the two directions of thes QSOs
				const double cosTheta = cosDe1*v_CosDeQ1__[q2]*cos(ra1-v_raQ1__[q2]) + sinDe1*v_SinDeQ1__[q2];
	
				///// Perpendicular distance
				const double distT = (rr1+v_rrQ1__[q2])*sqrt( 0.5*(1.-cosTheta) );
				if (distT >= max) continue;
	
				///// Parrallel distance
				const double distP = fabs(rr1-v_rrQ1__[q2])*sqrt( 0.5*(1.+cosTheta) );
				if (distP >= max) continue;
	
				const unsigned int idxPerp = int(distT);
				const unsigned int idxPara = int(distP);
	
				//// 1D
				const double distTotPow2 = distT*distT + distP*distP;
				if (distTotPow2<maxPow2) {
					const double distTot    = sqrt(distTotPow2);
					const double mu         = distP/distTot;
					const unsigned int idx  = int(distTot);
					const unsigned int idxM = int(50.*mu);
	
					dataMu[idx][idxM][0] ++;
					dataMu[idx][idxM][1] += distTot;
					dataMu[idx][idxM][2] += mu;
					dataMu[idx][idxM][3] += zz1+v_zzQ1__[q2];
				}
	
				//// 2D
				data2D[idxPerp][idxPara][0] ++;
				data2D[idxPerp][idxPara][1] += distT;
				data2D[idxPerp][idxPara][2] += distP;
				data2D[idxPerp][idxPara][3] += zz1+v_zzQ1__[q2];
			}
		}
		
		std::cout << "\n  Saving\n" << std::endl;
	
		prefix = QSO__;
		prefix += "_DR_";
		prefix += strBootIdx;
		prefix += ".txt";
	
		//// 2D
		pathToSave = pathToSave__;
		pathToSave += "xi_QSO_QSO_2D_";
		pathToSave += prefix;
		std::cout << "\n  " << pathToSave << std::endl;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);
	
		sumOne = 0.;
		meanZ  = 0.;
	
		//// Set the values of data
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBin; j++) {
				fFile << data2D[i][j][0];
				fFile << " " << data2D[i][j][1];
				fFile << " " << data2D[i][j][2];
				fFile << " " << data2D[i][j][3]/2.;
				fFile << std::endl;
	
				sumOne += data2D[i][j][0];
				meanZ  += data2D[i][j][3]/2.;
			}
		}
		fFile.close();
	
		//// Get mean redshift of pairs
		std::cout << "  < z >           = " << meanZ/sumOne << " +- " << 0. << std::endl;
		std::cout << "  number of pairs = " << sumOne          << std::endl;
	
	
		// Mu
		pathToSave = pathToSave__;
		pathToSave += "xi_QSO_QSO_Mu_";
		pathToSave += prefix;
		std::cout << "\n  " << pathToSave << std::endl;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);
	
		//// Set the values of data
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBinM; j++) {
				fFile << dataMu[i][j][0];
				fFile << " " << dataMu[i][j][1];
				fFile << " " << dataMu[i][j][2];
				fFile << " " << dataMu[i][j][3]/2.;
				fFile << std::endl;
			}
		}
		fFile.close();
	}

	return;
}
void Correlation::xi_Q1_Q2(void) {

	std::cout << "\n\n\n\n  ------ xi_Q1_Q2 ------\n" << std::endl;

	loadDataQ1();
	loadDataQ2();

	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max          = 160.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;

	const double maxPow2 = max*max;
	
	//// Arrays for data
	double meanZZ[2] = {};
	double data1D_0[nbBin];
	double data1D_1[nbBin];
	for (unsigned int i=0; i<nbBin; i++) {
		data1D_0[i] = 0.;
		data1D_1[i] = 0.;
	}

	long double data2D[nbBinX][nbBinY][3];
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			data2D[i][j][0] = 0.;
			data2D[i][j][1] = 0.;
			data2D[i][j][2] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int q1=0; q1<nbQ1__; q1++) {
		
		const double cosDe1 = v_CosDeQ1__[q1];
		const double sinDe1 = v_SinDeQ1__[q1];
		const double ra1    = v_raQ1__[q1];
		const double r1    = v_rrQ1__[q1];
		const double z1    = v_zzQ1__[q1];

		for (unsigned int q2=0; q2<nbQ2__; q2++) {

			//// Angle between the two directions of the qso and the lya
			double cosTheta = cosDe1*v_CosDeQ2__[q2]*cos(ra1-v_raQ2__[q2]) + sinDe1*v_SinDeQ2__[q2];
			if (fabs(1.-cosTheta)<1.e-10) cosTheta=1.;
			
			//// Transvers distance between the qso and the lya
			const double distTransQsoQsoPow2 = r1*r1*(1.-cosTheta*cosTheta);
			if (distTransQsoQsoPow2 >= maxPow2) continue;

			//// Parrallel distance between the qso and the lya
			const double distParalQsoQso = r1*cosTheta-v_rrQ2__[q2];
			if (fabs(distParalQsoQso) >= max) continue;
			const unsigned int rParralBinIdx = int(max+distParalQsoQso);

			//// Transvers distance between the qso and the lya
			const double distTransQsoQso   = sqrt(distTransQsoQsoPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoQso);

			const double distTotPow2 = distTransQsoQsoPow2 + distParalQsoQso*distParalQsoQso;

			if (distTotPow2 < maxPow2) {
				const double distTot = sqrt(distTotPow2);
				const unsigned int idx = int(distTot);
				data1D_0[idx] ++;
				data1D_1[idx] += distTot;

				const double tmp_meanZZ = z1+v_zzQ2__[q2];
				meanZZ[0] += tmp_meanZZ;
				meanZZ[1] += tmp_meanZZ*tmp_meanZZ;
			}
				
			//// Fill the histogramm of xi(r_{perp}, r_{paral}
			data2D[rPerpBinIdx][rParralBinIdx][0] ++;
			data2D[rPerpBinIdx][rParralBinIdx][1] += distTransQsoQso;
			data2D[rPerpBinIdx][rParralBinIdx][2] += distParalQsoQso;

		}
	}
	
	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string tmp_pathToSave = pathToSave__;
	tmp_pathToSave += "xi_Q1_Q2_1D_";
	tmp_pathToSave += QSO__;
	tmp_pathToSave += "_";
	tmp_pathToSave += QSO2__;
	tmp_pathToSave += ".txt";
	std::cout << "\n  " << tmp_pathToSave << std::endl;
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	long double sumOne = 0.;

	//// Set the values of data
	//// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {
		if (data1D_0[i]!=0.) {
			const long double save0 = data1D_0[i];
			const long double save1 = data1D_1[i];

			data1D_1[i] /= data1D_0[i];

			fFile << data1D_1[i];
			fFile << " " << data1D_0[i];
			fFile << " " << save0;
			fFile << " " << save1;
			fFile << " " << std::endl;

			sumOne += save0;
		}
	}
	fFile.close();



	//// Get mean redshift of pairs
	meanZZ[0] = meanZZ[0]/(2.*sumOne);
	meanZZ[1] = sqrt( (meanZZ[1]/(4.*sumOne) - meanZZ[0]*meanZZ[0])/sumOne);
	std::cout << "  < z >           = " << meanZZ[0] << " +- " << meanZZ[1] << std::endl;
	std::cout << "  number of pairs = " << sumOne          << std::endl;


	sumOne = 0.;

	
	//// Save the 2D cross-correlation
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_Q1_Q2_2D_";
	pathToSave += QSO__;
	pathToSave += "_";
	pathToSave += QSO2__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	//// Set the values of data
	//// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			if (data2D[i][j][0]!=0.) {
				const long double save0 = data2D[i][j][0];
				const long double save1 = data2D[i][j][1];
				const long double save2 = data2D[i][j][2];

				data2D[i][j][1] /= data2D[i][j][0];
				data2D[i][j][2] /= data2D[i][j][0];

				fFile << i*nbBinY + j;
				fFile << " " << data2D[i][j][1];
				fFile << " " << data2D[i][j][2];
				fFile << " " << data2D[i][j][0];
				fFile << " " << save0;
				fFile << " " << save1;
				fFile << " " << save2;
				fFile << " " << std::endl;

				sumOne += save0;
			}
		}
	}
	fFile.close();
	

	std::cout << "  number of pairs      = " << sumOne          << std::endl;
	

	return;
}




// ---------------------------------------------------------------------
//
//		Correlation function for Jean-Marc Le Goff
//			(i.e. mocks in euclidean geometry)
//
// ---------------------------------------------------------------------
void Correlation::xi_A_delta_delta_MockJMc(unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_MockJMc ------\n" << std::endl;

	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__, bootIdx);
	if (nbForest_==0.) return;

	if (mocksNoNoiseNoCont || mocks_raw) {
		removeFalseCorrelations();
	}
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Constants:
	const double max          = 200.;
	const unsigned int nbBin  = int(max);
	const unsigned int nbBinM = 50;
	const double maxPow2      = max*max;

	///// get an array for nb of pairs for the forest
	double a_nbPairs[nbForest_];
	for (unsigned int i=0; i<nbForest_; i++) {
		a_nbPairs[i] = 0.;
	}

	///// Arrays for data
	double dataMu[nbBin][nbBinM][9];
	double data2D[nbBin][nbBin][9];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<9; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<9; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}
	

	///// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl; 
		std::srand(seed_for_random_position__);

		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}
	}
	
	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {

		if (doBootstraps__ && v_region_Map__[f1]!=bootIdx) continue;

		const unsigned int nb1  = v_nbPixelDelta1__[f1];
		const double x1         = v_ra__[f1];
		const double y1         = v_de__[f1];
		const double fpix1      = v_r__[f1][0];
		const double lPix1      = v_r__[f1][nb1-1];

		double x11 = x1;
		double y11 = y1;
		///// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x11 = v_raRandForest[f1];
			y11 = v_deRandForest[f1];
		}

		for (unsigned int f2=0; f2<f1; f2++) {

			const double x2 = v_ra__[f2];
			const double y2 = v_de__[f2];

			///// Not in the same line of sight
			if (x1==x2 && y1==y2) continue;

			double x22 = x2;
			double y22 = y2;
			///// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x22 = v_raRandForest[f2];
				y22 = v_deRandForest[f2];
			}

			///// If r_perp is too high
			const double distTransPow2 = (x11-x22)*(x11-x22) + (y11-y22)*(y11-y22);
			if (distTransPow2 >= maxPow2) continue;

			///// If forests too far from one another
			const unsigned int nb2 = v_nbPixelDelta1__[f2];
			const double fpix2 = v_r__[f2][0];
			const double lPix2 = v_r__[f2][nb2-1];
			const double minLastPixel = std::min(lPix1,lPix2);
			if (minLastPixel==lPix1 && fpix2>lPix1 && fpix2-lPix1>=max ) continue;
			if (minLastPixel==lPix2 && fpix1>lPix2 && fpix1-lPix2>=max ) continue;

			///// Transvers distance
			const double rPerp = sqrt(distTransPow2);
			const unsigned int rPerpBinIdx = int( rPerp );

			for (unsigned int i1=0; i1<nb1; i1++) {

				const double r1 = v_r__[f1][i1];
				const double w1 = v_w__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double z1 = v_z__[f1][i1];
				const double res_lRF1  = v_residual_delta_vs_lRF__[f1][i1];
				const double res_lObs1 = v_residual_delta_vs_lObs__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-r1);
					if (rParral>=max) continue;
					const double distTotPow2 = distTransPow2 + rParral*rParral;

					const double w1w2                 = w1*v_w__[f2][i2];
					const double d1d2                 = d1*v_d__[f2][i2];
					const double w1w2_d1d2            = w1w2*d1d2;
					const double w1w2_d1d2_d1d2       = w1w2_d1d2*d1d2;
					const double w1w2_z1z2            = w1w2*(z1+v_z__[f2][i2]);
					const double w1w2_res_lRF1_lRF2   = w1w2*res_lRF1*v_residual_delta_vs_lRF__[f2][i2];
					const double w1w2_res_lObs1_lObs2 = w1w2*res_lObs1*v_residual_delta_vs_lObs__[f2][i2];

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int(mu*50.);

						dataMu[idx][idxM][0] += w1w2_d1d2;
						dataMu[idx][idxM][1] += w1w2_d1d2_d1d2;
						dataMu[idx][idxM][2] += w1w2*distTot;
						dataMu[idx][idxM][3] += w1w2*mu;
						dataMu[idx][idxM][4] += w1w2_z1z2;
						dataMu[idx][idxM][5] += w1w2;
						dataMu[idx][idxM][6] ++;
						dataMu[idx][idxM][7] += w1w2_res_lRF1_lRF2;
						dataMu[idx][idxM][8] += w1w2_res_lObs1_lObs2;
					}
					
					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += w1w2_d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][1] += w1w2_d1d2_d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2_z1z2;
					data2D[rPerpBinIdx][rParralBinIdx][5] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][6] ++;
					data2D[rPerpBinIdx][rParralBinIdx][7] += w1w2_res_lRF1_lRF2;
					data2D[rPerpBinIdx][rParralBinIdx][8] += w1w2_res_lObs1_lObs2;

					///// Get the number of pairs
					a_nbPairs[f1] += w1w2;
					a_nbPairs[f2] += w1w2;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	//// Set the prefix for different type of runs
	std::string prefix = forest__;
	if (doBootstraps__) {

		std::stringstream convert;
		convert << bootIdx;
		const std::string strBootIdx = convert.str();

		prefix += "_subsampling_";
		prefix += strBootIdx;
	}
	prefix += ".txt";

	std::ofstream fFile;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;
	long double sumOne = 0.;

	///// Save the 2D cross-correlation
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_2D_";
	pathToSave += prefix;
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	sumZZZ = 0.;
	sumWeg = 0.;
	sumOne = 0.;
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4]/2.;
			fFile << " " << data2D[i][j][5];
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << " " << data2D[i][j][8];
			fFile << std::endl;

			sumZZZ += data2D[i][j][4]/2.;
			sumWeg += data2D[i][j][5];
			sumOne += data2D[i][j][6];
		}
	}
	fFile.close();


	std::cout << "  < z >                = " << sumZZZ/sumWeg << std::endl;
	std::cout << "  number of pairs      = " << sumOne          << std::endl;
	std::cout << "  < pairs per forest > = " << sumOne/nbForest_ << std::endl;


	///// Mu
	pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_Mu_";
	pathToSave += prefix;
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << " " << dataMu[i][j][8];
			fFile << std::endl;
		}
	}
	fFile.close();

	if (!doBootstraps__ && !shuffleForest && !shuffleQSO && !randomQSO) {

		std::vector< std::vector< double > > forests;
		std::vector< double > tmp_forests_id;
		std::vector< double > tmp_forests_pa;
		
		for (unsigned int i=0; i<nbForest_; i++) {
			tmp_forests_id.push_back( v_idx__[i] );
			tmp_forests_pa.push_back( a_nbPairs[i] );
		}
		forests.push_back(tmp_forests_id);
		forests.push_back(v_ra__);
		forests.push_back(v_de__);
		forests.push_back(tmp_forests_pa);


		// Save the list of pairs
		pathToSave = pathToSave__;
		pathToSave += "xi_A_delta_delta_list_pairs_";
		pathToSave += forest__;
		pathToSave += ".txt";
		std::cout << "\n  " << pathToSave << std::endl;

		std::ofstream fFile;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);

		for (unsigned int i=0; i<nbForest_; i++) {
			fFile << v_idx__[i];
			fFile << " " << 0;
			fFile << " " << v_ra__[i];
			fFile << " " << v_de__[i];
			fFile << " " << a_nbPairs[i];
			fFile << std::endl;
		}
		fFile.close();
		

		//// find the index of each forest among the C_NBSUBSAMPLES sub-samples
		LymanForest* lymanForestObject = new LymanForest(forests, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		pathToSave = pathToSave__;
		pathToSave += "xi_A_delta_delta_map_";
		pathToSave += forest__;
		pathToSave += ".txt";
		lymanForestObject->SaveRegionMap(pathToSave);
		delete lymanForestObject;
	}

	return;
}
void Correlation::xi_A_delta_delta_MockJMc_distortionMatrix(void) {


	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_MockJMc_distortionMatrix ------\n" << std::endl;

	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl; 
		std::srand(seed_for_random_position__);

		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}
	}

	///// Set usefull vectors
	std::vector<double> v_invSumWeight(nbForest_,0.);
	std::vector<double> v_meanLambda(nbForest_,0.);
	std::vector<double> v_invVarLambda(nbForest_,0.);

	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nbPixel = v_nbPixelDelta1__[f];

		long double sumWeight  = 0.;
		long double meanLambda = 0.;
		long double varLambda  = 0.;

		///// Loops over all pixels of the forest
		for (unsigned int i=0; i<nbPixel; i++) {

			const long double w = v_w__[f][i];
			const long double l = v_lRF__[f][i];

			sumWeight  += w;
			meanLambda += w*l;
			varLambda  += w*l*l;
		}

		if (sumWeight==0.) std::cout << "   ERROR:: sumWeight = 0. " << std::endl;
		meanLambda        /= sumWeight;
		varLambda          = varLambda - meanLambda*meanLambda*sumWeight;
		if (varLambda==0.) std::cout << "   ERROR:: varLambda = 0. " << std::endl;

		v_invSumWeight[f] = 1./sumWeight;
		v_meanLambda[f]   = meanLambda;
		v_invVarLambda[f] = 1./varLambda;
	}

	///// Constants:
	const double max         = 200.;
	const double binSize     = 4.;
	const unsigned int nbBin = int(max/binSize);
	const unsigned int nbBin2D = nbBin*nbBin;

	const double maxPow2      = max*max;
	const double fromValToIdx  = nbBin/max;

	///// Arrays for distortion matrix
	double weight2D[nbBin2D];
	double** data2DMatrix = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		weight2D[i] = 0.;
		data2DMatrix[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			data2DMatrix[i][j] = 0.;
		}
	}





	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {

		const unsigned int nb1  = v_nbPixelDelta1__[f1];
		const double x1         = v_ra__[f1];
		const double y1         = v_de__[f1];
		const double fpix1      = v_r__[f1][0];
		const double lPix1      = v_r__[f1][nb1-1];
		const double invSumWeight1  = v_invSumWeight[f1];
		const double meanLambda1    = v_meanLambda[f1];
		const double invVarLambda1  = v_invVarLambda[f1];

		double x11 = x1;
		double y11 = y1;
		///// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x11 = v_raRandForest[f1];
			y11 = v_deRandForest[f1];
		}

		for (unsigned int f2=0; f2<f1; f2++) {

			const double x2 = v_ra__[f2];
			const double y2 = v_de__[f2];

			///// Not in the same line of sight
			if (x1==x2 && y1==y2) continue;

			double x22 = x2;
			double y22 = y2;
			///// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x22 = v_raRandForest[f2];
				y22 = v_deRandForest[f2];
			}

			///// If r_perp is too high
			const double distTransPow2 = (x11-x22)*(x11-x22) + (y11-y22)*(y11-y22);
			if (distTransPow2 >= maxPow2) continue;

			///// If forests too far from one another
			const unsigned int nb2 = v_nbPixelDelta1__[f2];
			const double fpix2 = v_r__[f2][0];
			const double lPix2 = v_r__[f2][nb2-1];
			const double minLastPixel = std::min(lPix1,lPix2);
			if (minLastPixel==lPix1 && fpix2>lPix1 && fpix2-lPix1>=max ) continue;
			if (minLastPixel==lPix2 && fpix1>lPix2 && fpix1-lPix2>=max ) continue;

			///// Transvers distance
			const double rPerp = sqrt(distTransPow2);
			const unsigned int rPerpBinIdx = int( rPerp*fromValToIdx );

			const double invSumWeight2  = v_invSumWeight[f2];
			const double meanLambda2    = v_meanLambda[f2];
			const double invVarLambda2  = v_invVarLambda[f2];

			//// For the distortion-matrix
			std::vector< std::vector< double > > list_for_pixel_1_0(nb1, std::vector<double>(2500,0.));
			std::vector< std::vector< double > > list_for_pixel_2_0(nb2, std::vector<double>(2500,0.));
			std::vector< std::vector< double > > list_for_pixel_1_1(nb1, std::vector<double>(2500,0.));
			std::vector< std::vector< double > > list_for_pixel_2_1(nb2, std::vector<double>(2500,0.));
			std::vector< bool > pixelTouched_1(nb1,false);
			std::vector< bool > pixelTouched_2(nb2,false);
			std::vector< std::vector< bool > > binTouched_1(nb1, std::vector<bool>(2500,false));
			std::vector< std::vector< bool > > binTouched_2(nb2, std::vector<bool>(2500,false));
			double value2[2500]   = {0.};
			double value5[2500]   = {0.};
			double value6[2500]   = {0.};
			double value7[2500]   = {0.};
			bool binTouched[2500] = {false};

			///// Index of the bin
			std::vector< unsigned int > binIdx;
			std::vector< unsigned int > fromBinsToPixels1;
			std::vector< unsigned int > fromBinsToPixels2;

			for (unsigned int i1=0; i1<nb1; i1++) {

				const double r1 = v_r__[f1][i1];
				const double w1 = v_w__[f1][i1];
				const double lambda_minus_mean_over_var1 = (v_lRF__[f1][i1] -meanLambda1)*invVarLambda1;

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-r1);
					if (rParral>=max) continue;
					const unsigned int rParralBinIdx = int(rParral*fromValToIdx);
					const unsigned int globalBin = rPerpBinIdx*nbBin +rParralBinIdx;

					const double w2 = v_w__[f2][i2];
					const double lambda_minus_mean_over_var2 = (v_lRF__[f2][i2] -meanLambda2)*invVarLambda2;

					//// Fill arrays for the distortion
					list_for_pixel_1_0[i1][globalBin] += w2*invSumWeight2;
					list_for_pixel_2_0[i2][globalBin] += w1*invSumWeight1;
					list_for_pixel_1_1[i1][globalBin] += w2*lambda_minus_mean_over_var2;
					list_for_pixel_2_1[i2][globalBin] += w1*lambda_minus_mean_over_var1;
					pixelTouched_1[i1] = true;
					pixelTouched_2[i2] = true;
					binTouched_1[i1][globalBin] = true;
					binTouched_2[i2][globalBin] = true;
					value2[globalBin] += w1*w2*invSumWeight1*invSumWeight2;
					value5[globalBin] += w2*invSumWeight2*w1*lambda_minus_mean_over_var1;
					value6[globalBin] += w1*invSumWeight1*w2*lambda_minus_mean_over_var2;
					value7[globalBin] += w1*w2*lambda_minus_mean_over_var1*lambda_minus_mean_over_var2;
					binIdx.push_back(globalBin);
					fromBinsToPixels1.push_back(i1);
					fromBinsToPixels2.push_back(i2);
					binTouched[globalBin] = true;

					///// Fill array of weights
					weight2D[globalBin] += w1*w2;
				}
			}

			///// Number of pairs of two pixels ( here = nb*nb)
			const unsigned int nbPixelsWithPairs = binIdx.size();
	
			///// Get a vector of not empty bins
			std::vector< unsigned int> bins_touched;
			for (unsigned int b=0; b<2500; b++) {
				if (binTouched[b]) bins_touched.push_back(b);
			}
			const unsigned int nbBinsTouched = bins_touched.size();

			//// For forest 1
			std::vector< unsigned int > nbBinsTouched_pixel_1(nb1);
			std::vector< std::vector< unsigned int > > bins_Touched_pixel_1(nb1);
			for (unsigned int i=0; i<nb1; i++) {
				if (!pixelTouched_1[i]) continue;
				std::vector< unsigned int > tmp_bins_Touched_pixel_1;
				for (unsigned int b=0; b<2500; b++) {
					if (binTouched_1[i][b]) tmp_bins_Touched_pixel_1.push_back(b);
				}
				bins_Touched_pixel_1[i] = tmp_bins_Touched_pixel_1;
				nbBinsTouched_pixel_1[i] = tmp_bins_Touched_pixel_1.size();
			}
			//// For forest 2
			std::vector< unsigned int > nbBinsTouched_pixel_2(nb2);
			std::vector< std::vector< unsigned int > > bins_Touched_pixel_2(nb2);
			for (unsigned int i=0; i<nb2; i++) {
				if (!pixelTouched_2[i]) continue;
				std::vector< unsigned int > tmp_bins_Touched_pixel_2;
				for (unsigned int b=0; b<2500; b++) {
					if (binTouched_2[i][b]) tmp_bins_Touched_pixel_2.push_back(b);
				}
				bins_Touched_pixel_2[i] = tmp_bins_Touched_pixel_2;
				nbBinsTouched_pixel_2[i] = tmp_bins_Touched_pixel_2.size();
			}
	
			///// Loops over all pairs
			for (unsigned int pa=0; pa<nbPixelsWithPairs; pa++) {
	
				const unsigned int globalBin1 = binIdx[pa];
				const unsigned int i1 = fromBinsToPixels1[pa];
				const unsigned int i2 = fromBinsToPixels2[pa];
	
				const double w1w2 = v_w__[f1][i1]*v_w__[f2][i2];
				const double valf1i1 = v_lRF__[f1][i1] -v_meanLambda[f1];
				const double valf2i2 = v_lRF__[f2][i2] -v_meanLambda[f2];
	
				///// Fill the distortion matrix
				data2DMatrix[globalBin1][globalBin1] += w1w2;

				/// Loop over all not empty bins
				for (unsigned int b2=0; b2<nbBinsTouched; b2++) {
					const unsigned int globalBin2 = bins_touched[b2];
					data2DMatrix[globalBin1][globalBin2] += w1w2*( value2[globalBin2]
						+valf1i1*value5[globalBin2] +valf2i2*value6[globalBin2]
						+valf1i1*valf2i2*value7[globalBin2] );
				}
				
				for (unsigned int b2=0; b2<nbBinsTouched_pixel_1[i1]; b2++) {
					const unsigned int globalBin2 = bins_Touched_pixel_1[i1][b2];
					data2DMatrix[globalBin1][globalBin2] -= w1w2*( list_for_pixel_1_0[i1][globalBin2] + valf2i2*list_for_pixel_1_1[i1][globalBin2] );
				}
				for (unsigned int b2=0; b2<nbBinsTouched_pixel_2[i2]; b2++) {
					const unsigned int globalBin2 = bins_Touched_pixel_2[i2][b2];
					data2DMatrix[globalBin1][globalBin2] -= w1w2*( list_for_pixel_2_0[i2][globalBin2] + valf1i1*list_for_pixel_2_1[i2][globalBin2] );
				}
			}
		}
std::cout << f1 << std::endl;
	}


        std::cout << "\n  Saving\n" << std::endl;

        std::ofstream fFile;

	///// Save the 2D cross-correlation 
        std::string pathToSave = pathToSave__;
        pathToSave += "xi_A_delta_delta_distortionMatrix_2D_";
        pathToSave += forest__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	///// Set the values of data
        ///// [0] for value, [1] for error, [2] for bin center
	long double value_total  = 0.;
	long double weight_total = 0.;
        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
                        double value = 0.;
			if (weight2D[i]!=0.) {
				value = data2DMatrix[i][j]/weight2D[i];
				value_total  += data2DMatrix[i][j];
			}
                        fFile << value << " ";
                }
                fFile << std::endl;
		weight_total += weight2D[i];
        }
        fFile.close();

	std::cout << "  value_total     = " << value_total << std::endl;
	std::cout << "  weight_total    = " << weight_total << std::endl;
	std::cout << "  integral matrix = " << value_total / weight_total << std::endl;


	return;

}
void Correlation::xi_A_delta_delta_Metals_Models_MockJMc(double lambdaRFMetal1, std::string lambdaRFMetalName1,double lambdaRFMetal2, std::string lambdaRFMetalName2) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_Metals_Models ------\n" << std::endl;

	std::cout << "  line1 = " << lambdaRFMetalName1 << " : " << lambdaRFMetal1 << std::endl;
	std::cout << "  line2 = " << lambdaRFMetalName2 << " : " << lambdaRFMetal2 << std::endl;

	///// Load forest
	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
	
	/// Get an arrey to get te growth factor
	double from_z_to_growth_factor_pow_2[C_NBBINREDSH] = {0.};
	const double g0    = cosmo->get_growth_factor(0.);
	const double growth_factor_step_size = C_ZEXTREMABINCONVERT1/C_NBBINREDSH;
	const double inverse_growth_factor_step_size = 1./growth_factor_step_size;
	for (unsigned int i=0; i<C_NBBINREDSH; i++) {
		const double ZZZ = i*growth_factor_step_size;
		from_z_to_growth_factor_pow_2[i] = (cosmo->get_growth_factor(ZZZ)/g0)*(cosmo->get_growth_factor(ZZZ)/g0);
	}

	//// Get the distance if it was this metal
	std::vector< std::vector< double > > v_r_metal1(v_r__);
	std::vector< std::vector< double > > v_r_metal2(v_r__);
	for (unsigned int i=0; i<v_r_metal1.size(); i++) {
		for (unsigned int j=0; j<v_r_metal1[i].size(); j++) {
			v_r_metal1[i][j] = hConvertRedshDist->Interpolate( v_lObs__[i][j]/lambdaRFMetal1-1. );
			v_r_metal2[i][j] = hConvertRedshDist->Interpolate( v_lObs__[i][j]/lambdaRFMetal2-1. );
		}
	}

	//// Empty useless vectors
	delete cosmo;
	delete hConvertRedshDist;
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_d__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	//// Get the monopole, quadrupol, exadecapol
	std::string pathToLoad = PATHTOWORK;
	pathToLoad += C_PATHTOCAMB;
	std::ifstream fileData(pathToLoad.c_str());
	std::vector< double > data_x;
	std::vector< double > data_xi0;
	std::vector< double > data_xi2;
	std::vector< double > data_xi4;
	std::cout << pathToLoad << std::endl;
	while (fileData) {
		double x;
		double xi0;
		double xi2;
		double xi4;
		fileData>>x>>xi0>>xi2>>xi4;
		if (fileData==0) break;
		
		data_x.push_back(x);
		data_xi0.push_back(xi0);
		data_xi2.push_back(xi2);
		data_xi4.push_back(xi4);
	}

	const double inverse_step_size = 1./data_x[0];
	const unsigned int nbBinCAMB = data_x.size();
	//// Set to zero the last pixel
	const double maxDistCAMB = data_x[nbBinCAMB-1];
	data_x.push_back(data_x[nbBinCAMB-1]+1.);
	data_xi0.push_back(0.);
	data_xi2.push_back(0.);
	data_xi4.push_back(0.);

	///// Constants:
	const double max          = 200.;
	const unsigned int nbBin  = int(max);
	const unsigned int nbBinM = 50;

	const double maxPow2      = max*max;

	///// Arrays for data
	double dataMu[nbBin][nbBinM][8];
	double data2D[nbBin][nbBin][8];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<8; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<8; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}

	///// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand (seed_for_random_position__);

		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);

		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}
	}

	
	std::cout << "\n  Starting\n" << std::endl;			
	

	for (unsigned int f1=0; f1<nbForest_; f1++) {
		const unsigned int nb1  = v_nbPixelDelta1__[f1];
		const double x1         = v_ra__[f1];
		const double y1         = v_de__[f1];
		const double fpix1      = v_r__[f1][0];
		const double lPix1      = v_r__[f1][nb1-1];

		double x11 = x1;
		double y11 = y1;
		///// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x11 = v_raRandForest[f1];
			y11 = v_deRandForest[f1];
		}

		for (unsigned int f2=0; f2<f1; f2++) {

			const double x2 = v_ra__[f2];
			const double y2 = v_de__[f2];

			///// Not in the same line of sight
			if (x1==x2 && y1==y2) continue;

			double x22 = x2;
			double y22 = y2;
			///// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x22 = v_raRandForest[f2];
				y22 = v_deRandForest[f2];
			}

			///// If r_perp is too high
			const double distTransPow2 = (x11-x22)*(x11-x22) + (y11-y22)*(y11-y22);
			if (distTransPow2 >= maxPow2) continue;

			///// If forests too far from one another
			const unsigned int nb2 = v_nbPixelDelta1__[f2];
			const double fpix2 = v_r__[f2][0];
			const double lPix2 = v_r__[f2][nb2-1];
			const double minLastPixel = std::min(lPix1,lPix2);
			if (minLastPixel==lPix1 && fpix2>lPix1 && fpix2-lPix1>=max ) continue;
			if (minLastPixel==lPix2 && fpix1>lPix2 && fpix1-lPix2>=max ) continue;

			///// Transvers distance
			const double rPerp = sqrt(distTransPow2);
			const unsigned int rPerpBinIdx = int( rPerp );

			for (unsigned int i1=0; i1<nb1; i1++) {

				const double r1  = v_r__[f1][i1];
				const double rm1 = v_r_metal1[f1][i1];
				const double w1  = v_w__[f1][i1];
				const double z1  = v_z__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-r1);
					if (rParral>=max) continue;
					const double distTotPow2 = distTransPow2 + rParral*rParral;
					const double distTotMetal = sqrt( distTransPow2 + (v_r_metal2[f2][i2]-rm1)*(v_r_metal2[f2][i2]-rm1) );

					const double w1w2     = w1*v_w__[f2][i2];
					const double w1w2z1z2 = w1w2*(z1+v_z__[f2][i2]);

					unsigned int idxBinCAMB = nbBinCAMB;
					if (distTotMetal<maxDistCAMB) idxBinCAMB = int(distTotMetal*inverse_step_size);
					const double growth_factor_pow_2 = from_z_to_growth_factor_pow_2[ int(0.5*(z1+v_z__[f2][i2])*inverse_growth_factor_step_size) ];
					const double wxi0 = w1w2*growth_factor_pow_2*data_xi0[idxBinCAMB];
					const double wxi2 = w1w2*growth_factor_pow_2*data_xi2[idxBinCAMB];
					const double wxi4 = w1w2*growth_factor_pow_2*data_xi4[idxBinCAMB];

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int(mu*50.);

						dataMu[idx][idxM][0] += wxi0;
						dataMu[idx][idxM][1] += wxi2;
						dataMu[idx][idxM][2] += wxi4;
						dataMu[idx][idxM][3] += w1w2*distTot;
						dataMu[idx][idxM][4] += w1w2*mu;
						dataMu[idx][idxM][5] += w1w2z1z2;
						dataMu[idx][idxM][6] += w1w2;
						dataMu[idx][idxM][7] ++;
					}
					
					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += wxi0;
					data2D[rPerpBinIdx][rParralBinIdx][1] += wxi2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += wxi4;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][5] += w1w2z1z2;
					data2D[rPerpBinIdx][rParralBinIdx][6] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][7] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	//// Set the prefix of different forest and QSOs
	std::string prefix1 = forest__;
	prefix1 += "_";
	prefix1 += lambdaRFMetalName1;
	prefix1 += "_";
	prefix1 += lambdaRFMetalName2;

	std::ofstream fFile;
	std::string pathToSave;

	///// Save the 2D cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_Metals_Models_2D_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4];
			fFile << " " << data2D[i][j][5]/2.;
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();


	///// Mu
	pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_Metals_Models_Mu_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4];
			fFile << " " << dataMu[i][j][5]/2.;
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_MockJMc(unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc ------\n" << std::endl;

	//// QSO
	loadDataQ1();
	//// Forest
	if (mocks_raw) loadDataForest_Raw(bootIdx);
	else loadDataForest(pathForest__,bootIdx);
	if (nbForest_==0.) return;

	if (mocksNoNoiseNoCont || mocks_raw ) removeFalseCorrelations();
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Schiffle the forest
	if (shuffleForest) {

		//// Vectors with index of forest
		std::vector<unsigned int> randomIdx(nbForest_);
		for (unsigned int i=0; i<nbForest_; i++) {
			randomIdx[i] = i;
		}
		
		std::cout << "  Shuffle forest " << bootIdx << std::endl;
		bool doLoop = true;
		unsigned int nbLoop = 0;
		std::srand (bootIdx*10);
		while (doLoop) {
			doLoop = false;
			nbLoop ++;
			std::cout << "  nbLoop = " << nbLoop << std::endl;
			std::random_shuffle ( randomIdx.begin(), randomIdx.end() );
			for (unsigned int i=0; i<nbForest_; i++) {
				if (i==randomIdx[i]) {
					std::cout << nbLoop << " " << i << std::endl;
					doLoop = true;
					break;
				}
			}
		}

		//// Copy data
		std::vector<double> tmp_ra(v_ra__);
		std::vector<double> tmp_de(v_de__);
		
		//// Put the new data
		for (unsigned int i=0; i<nbForest_; i++) {
			const unsigned int ii = randomIdx[i];
			v_ra__[i]    = tmp_ra[ii];
			v_de__[i]    = tmp_de[ii];
		}
	}


	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max          =  200.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100;

	const double maxPow2      = max*max;
	//// get an array for nb of pairs for the forest
	double a_nbPairs[nbForest_];
	for (unsigned int i=0; i<nbForest_; i++) {
		a_nbPairs[i] = 0.;
	}


	//// Arrays for data
	double dataMu[nbBin][nbBinM][9];
	double data2D[nbBinX][nbBinY][9];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<9; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<9; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}


	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	std::vector<double> v_zzRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		//// Create the conversion table from redshift to distance
		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TGraph* gr_from_dist_to_z = cosmo->create_TGraph_to_convert_distance_to_redshift(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand (seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		v_zzRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_zzRandQSO[i] = gr_from_dist_to_z->Eval(v_rrRandQSO[i]);
		}

		delete gr_from_dist_to_z;
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		double nbPairs = 0.;

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];

		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];
			double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
				z2 = v_zzRandQSO[q];
			}

			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya );
			const double zQSO = z2;
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;

				if (fabs(distP) >= max) continue;

				const double distTotPow2 = distTransQsoLyaPow2 + distP*distP;
				const double w           = v_w__[f][i];
				const double d           = v_d__[f][i];
				const double wd          = w*d;
				const double wdd         = wd*d;
				const double wZpair      = w*(zQSO+v_z__[f][i]);
				const double w_res_lRF   = w*v_residual_delta_vs_lRF__[f][i];
                                const double w_res_lObs  = w*v_residual_delta_vs_lObs__[f][i];

				if (distTotPow2 < maxPow2) {
					const double distTot    = sqrt(distTotPow2);
					const double mu         = distP/distTot;
					const unsigned int idx  = int(distTot);
					const unsigned int idxM = int((1.+mu)*50.);

					dataMu[idx][idxM][0] += wd;
					dataMu[idx][idxM][1] += wdd;
					dataMu[idx][idxM][2] += w*distTot;
					dataMu[idx][idxM][3] += w*mu;
					dataMu[idx][idxM][4] += wZpair;
					dataMu[idx][idxM][5] += w;
					dataMu[idx][idxM][6] ++;
					dataMu[idx][idxM][7] += w_res_lRF;
					dataMu[idx][idxM][8] += w_res_lObs;
				}
	
				//// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rParralBinIdx = int(distP+max);

				data2D[rPerpBinIdx][rParralBinIdx][0] += wd;
				data2D[rPerpBinIdx][rParralBinIdx][1] += wdd;
				data2D[rPerpBinIdx][rParralBinIdx][2] += w*distTransQsoLya;
				data2D[rPerpBinIdx][rParralBinIdx][3] += w*distP;
				data2D[rPerpBinIdx][rParralBinIdx][4] += wZpair;
				data2D[rPerpBinIdx][rParralBinIdx][5] += w;
				data2D[rPerpBinIdx][rParralBinIdx][6] ++;
				data2D[rPerpBinIdx][rParralBinIdx][7] += w_res_lRF;
				data2D[rPerpBinIdx][rParralBinIdx][8] += w_res_lObs;

				//// Get the number of pairs
				nbPairs += w;
			}
		}

		//// Put the number of pairs comming with this forest
		a_nbPairs[f] = nbPairs;
	}

	std::cout << "\n  Saving\n" << std::endl;

	//// Set the prefix of different forest and QSOs
	std::string prefix1 = forest__;
	prefix1 += "_";
	prefix1 += QSO__;

	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	//// Set the prefix for different type of runs
	std::string prefix = "_";
	if (doBootstraps__)  prefix += "subsampling";
	if (shuffleForest) prefix += "shuffleForest";
	if (shuffleQSO)    prefix += "shuffleQSO";
	if (randomQSO)     prefix += "randomQSO";
	prefix += "_";
	prefix += strBootIdx;

	std::ofstream fFile;

	long double sumZZZ = 0.;
	long double sumWeg = 0.;
	long double sumOne = 0.;

	//// Save the 2D cross-correlation
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	sumZZZ = 0.;
	sumWeg = 0.;
	sumOne = 0.;
	
	//// Set the values of data
	//// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4]/2.;
			fFile << " " << data2D[i][j][5];
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << " " << data2D[i][j][8];
			fFile << std::endl;

			sumZZZ += data2D[i][j][4]/2.;
			sumWeg += data2D[i][j][5];
			sumOne += data2D[i][j][6];
		}
	}
	fFile.close();

	std::cout << "  < z >                = " << sumZZZ/sumWeg    << std::endl;
	std::cout << "  number of pairs      = " << sumOne           << std::endl;
	std::cout << "  sum of weight        = " << sumWeg           << std::endl;
	std::cout << "  < pairs per forest > = " << sumOne/nbForest_ << std::endl;
	std::cout << "  < pairs per QSO >    = " << sumOne/nbQ1__    << std::endl;


	//// Save the Mu cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_Mu_";
	pathToSave += prefix1;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	//// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4]/2.;
			fFile << " " << dataMu[i][j][5];
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << " " << dataMu[i][j][8];
			fFile << std::endl;
		}
	}
	fFile.close();

	
	if (!doBootstraps__ && !shuffleForest && !shuffleQSO && !randomQSO) {

		std::vector< std::vector< double > > forests;
		std::vector< double > tmp_forests_id;
		std::vector< double > tmp_forests_pa;
		
		for (unsigned int i=0; i<nbForest_; i++) {
			tmp_forests_id.push_back( v_idx__[i] );
			tmp_forests_pa.push_back( a_nbPairs[i] );
		}
		forests.push_back(tmp_forests_id);
		forests.push_back(v_ra__);
		forests.push_back(v_de__);
		forests.push_back(tmp_forests_pa);

		// Save the list of pairs
		pathToSave = pathToSave__;
		pathToSave += "xi_delta_QSO_list_pairs_";
		pathToSave += prefix1;
		pathToSave += ".txt";
		std::cout << "\n  " << pathToSave << std::endl;

		std::ofstream fFile;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);

		for (unsigned int i=0; i<nbForest_; i++) {
			fFile << v_idx__[i];
			fFile << " " << 0;
			fFile << " " << v_ra__[i];
			fFile << " " << v_de__[i];
			fFile << " " << a_nbPairs[i];
			fFile << std::endl;
		}
		fFile.close();
		

		//// find the index of each forest among the C_NBSUBSAMPLES sub-samples
		LymanForest* lymanForestObject = new LymanForest(forests, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		pathToSave = pathToSave__;
		pathToSave += "xi_delta_QSO_map_";
		pathToSave += prefix1;
		pathToSave += ".txt";
		lymanForestObject->SaveRegionMap(pathToSave);
		delete lymanForestObject;
	}

	return;
}
void Correlation::xi_delta_QSO_MockJMc_distortionMatrix(void) {
	
	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_distortionMatrix ------\n" << std::endl;

	//// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	//// Forest
	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	//// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_d__.clear();
	v_z__.clear();
	v_idx__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Set usefull vectors
	std::vector<double> v_invSumWeight(nbForest_,0.);
	std::vector< std::vector< double > > v_varLambda(v_lRF__);

	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nbPixel = v_nbPixelDelta1__[f];

		double sumWeight  = 0.;
		double meanLambda = 0.;
		double stdLambda  = 0.;

		//// Loops over all pixels of the forest
		for (unsigned int i=0; i<nbPixel; i++) {

			const double w = v_w__[f][i];
			const double l = v_lRF__[f][i];

			sumWeight  += w;
			meanLambda += w*l;
			stdLambda  += w*l*l;
		}

		meanLambda        /= sumWeight;
		stdLambda          = 1./sqrt(stdLambda/sumWeight - meanLambda*meanLambda);
		v_invSumWeight[f]  = 1./sumWeight;

		for (unsigned int i=0; i<nbPixel; i++) {
			v_varLambda[f][i] = (v_lRF__[f][i]-meanLambda)*stdLambda;
		}
	}

	//// Clear usless vector
	v_lRF__.clear();
	
	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand(seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}
	}



	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max           = 200.;
	const double binSize       = 4.;
	const unsigned int nbBin   = int(max/binSize);
	const unsigned int nbBinX  = nbBin;
	const unsigned int nbBinY  = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double maxPow2       = max*max;
	const double fromValToIdx  = nbBin/max;

	//// Arrays for distortion matrix
	double weight[nbBin2D];
	double** dataMatrix = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		weight[i] = 0.;
		dataMatrix[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			dataMatrix[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];
		
		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}
		const double invSumWeight  = v_invSumWeight[f];

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
			}
			
			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya*fromValToIdx );
			
			//// get the weight and mean lambda
			double xValue[5000]  = {0.};
			double xlValue[5000] = {0.};
			bool binTouched[5000] = {false};

			//// Index of the bin
			std::vector< unsigned int > binIdx;
			std::vector< unsigned int > fromBinsToPixels;
	
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;
				if (fabs(distP) >= max) continue;

				const double w    = v_w__[f][i];
				const double val0 = w*invSumWeight;
	
				//// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int globalBin = rPerpBinIdx*nbBinY + int( (distP+max)*fromValToIdx );
				xValue[globalBin]  += val0;
				xlValue[globalBin] += val0*v_varLambda[f][i];

				//// Fill array of weights
				weight[globalBin] += w;

				//// Keep values for the distortion matrix
				binIdx.push_back(globalBin);
				fromBinsToPixels.push_back(i);

				binTouched[globalBin] = true;
			}

			//// Number of pixels with pairs 
			const unsigned int nbPixelsWithPairs = binIdx.size();

			//// Get a vector of not empty bins
			std::vector< unsigned int> bins_touched;
			for (unsigned int i=0; i<5000; i++) {
				if (binTouched[i]) bins_touched.push_back(i);
			}
			const unsigned int nbBinsTouched = bins_touched.size();

			//// Loops over all pixels of the forest (Fill the distortion matrix)
			for (unsigned int i=0; i<nbPixelsWithPairs; i++) {

				const unsigned int globalBin1 = binIdx[i];
				const unsigned int pixelIdx   = fromBinsToPixels[i];

				const double w    = v_w__[f][pixelIdx];
				const double val1 = v_varLambda[f][pixelIdx];
				
				//// Fill the distortion matrix
				dataMatrix[globalBin1][globalBin1] += w;

				for (unsigned int j=0; j<nbBinsTouched; j++) {
					const unsigned int globalBin2 = bins_touched[j];
					dataMatrix[globalBin1][globalBin2] -= w*( xValue[globalBin2] + val1*xlValue[globalBin2] );
				}
			}
		}
	}
	
	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;

	//// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_distortionMatrix_2D_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	//// Set the values of data
	//// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			double value = 0.;
			if (weight[i]!=0.) value = dataMatrix[i][j]/weight[i];
			fFile << value << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_MockJMc_distortionMatrix_1D(void) {
	
	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_distortionMatrix_1D ------\n" << std::endl;

	//// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	//// Forest
	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	//// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_d__.clear();
	v_z__.clear();
	v_idx__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Set usefull vectors
	std::vector<double> v_invSumWeight(nbForest_,0.);
	std::vector< std::vector< double > > v_varLambda(v_lRF__);

	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nbPixel = v_nbPixelDelta1__[f];

		double sumWeight  = 0.;
		double meanLambda = 0.;
		double stdLambda  = 0.;

		//// Loops over all pixels of the forest
		for (unsigned int i=0; i<nbPixel; i++) {

			const double w = v_w__[f][i];
			const double l = v_lRF__[f][i];

			sumWeight  += w;
			meanLambda += w*l;
			stdLambda  += w*l*l;
		}

		meanLambda        /= sumWeight;
		stdLambda          = 1./sqrt(stdLambda/sumWeight - meanLambda*meanLambda);
		v_invSumWeight[f]  = 1./sumWeight;

		for (unsigned int i=0; i<nbPixel; i++) {
			v_varLambda[f][i] = (v_lRF__[f][i]-meanLambda)*stdLambda;
		}
	}

	//// Clear usless vector
	v_lRF__.clear();
	
	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand(seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}
	}



	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max           = 200.;
	const double binSize       = 4.;
	const unsigned int nbBin   = int(max/binSize);
	const double maxPow2       = max*max;
	const double fromValToIdx  = nbBin/max;

	//// Arrays for distortion matrix
	double weight[nbBin];
	double** dataMatrix = new double*[nbBin];
	for (unsigned int i=0; i<nbBin; i++) {
		weight[i] = 0.;
		dataMatrix[i] = new double[nbBin];
		for (unsigned int j=0; j<nbBin; j++) {
			dataMatrix[i][j] = 0.;
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];
		
		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}
		const double invSumWeight  = v_invSumWeight[f];

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
			}

			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;
			
			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;
			
			//// get the weight and mean lambda
			double xValue[50]  = {0.};
			double xlValue[50] = {0.};
			bool binTouched[50] = {false};

			//// Index of the bin
			std::vector< unsigned int > binIdx;
			std::vector< unsigned int > fromBinsToPixels;
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;
				if (fabs(distP) >= max) continue;
				const double distTotPow2 = distTransQsoLyaPow2 + distP*distP;
				if (distTotPow2 >= maxPow2) continue;

				const double w    = v_w__[f][i];
				const double val0 = w*invSumWeight;
	
				//// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int globalBin = int(sqrt(distTotPow2)*fromValToIdx);
				xValue[globalBin]  += val0;
				xlValue[globalBin] += val0*v_varLambda[f][i];

				//// Fill array of weights
				weight[globalBin] += w;

				//// Keep values for the distortion matrix
				binIdx.push_back(globalBin);
				fromBinsToPixels.push_back(i);
				
				binTouched[globalBin] = true;
			}

			//// Number of pixels with pairs 
			const unsigned int nbPixelsWithPairs = binIdx.size();
			
			//// Get a vector of not empty bins
			std::vector< unsigned int> bins_touched;
			for (unsigned int i=0; i<50; i++) {
				if (binTouched[i]) bins_touched.push_back(i);
			}
			const unsigned int nbBinsTouched = bins_touched.size();

			//// Loops over all pixels of the forest (Fill the distortion matrix)
			for (unsigned int i=0; i<nbPixelsWithPairs; i++) {

				const unsigned int globalBin1 = binIdx[i];
				const unsigned int pixelIdx   = fromBinsToPixels[i];

				const double w    = v_w__[f][pixelIdx];
				const double val1 = v_varLambda[f][pixelIdx];
				
				//// Fill the distortion matrix
				dataMatrix[globalBin1][globalBin1] += w;

				for (unsigned int j=0; j<nbBinsTouched; j++) {
					const unsigned int globalBin2 = bins_touched[j];
					dataMatrix[globalBin1][globalBin2] -= w*( xValue[globalBin2] + val1*xlValue[globalBin2] );
				}
			}
		}
	}
	
	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	//// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_distortionMatrix_1D_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	//// Set the values of data
	//// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			double value = 0.;
			if (weight[i]!=0.) value = dataMatrix[i][j]/weight[i];
			fFile << value << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_Metals_Models_MockJMc(double lambdaFrMetal, std::string lambdaFrMetalName) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Metals_Models_MockJMc ------" << std::endl;
	std::string command = "  ";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	std::cout << "  line = " << lambdaFrMetalName << " : " << lambdaFrMetal << std::endl;

	//// QSO
	loadDataQ1();
	//// Forest
	if (mocks_raw) loadDataForest_Raw();
	else loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	/// Get an arrey to get te growth factor
	double from_z_to_growth_factor_pow_2[C_NBBINREDSH] = {0.};
	const double g0    = cosmo->get_growth_factor(0.);
	const double growth_factor_step_size = C_ZEXTREMABINCONVERT1/C_NBBINREDSH;
	const double inverse_growth_factor_step_size = 1./growth_factor_step_size;
	for (unsigned int i=0; i<C_NBBINREDSH; i++) {
		const double ZZZ = i*growth_factor_step_size;
		from_z_to_growth_factor_pow_2[i] = (cosmo->get_growth_factor(ZZZ)/g0)*(cosmo->get_growth_factor(ZZZ)/g0);
	}

	//// Get the distance if it was this metal
	std::vector< std::vector< double > > v_r_metal(v_r__);
	for (unsigned int i=0; i<v_r_metal.size(); i++) {
		for (unsigned int j=0; j<v_r_metal[i].size(); j++) {
			v_r_metal[i][j] = hConvertRedshDist->Interpolate( v_lObs__[i][j]/lambdaFrMetal-1. );
		}
	}

	delete cosmo;
	delete hConvertRedshDist;
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_d__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Get the monopole, quadrupol, exadecapol
	std::string pathToLoad = PATHTOWORK;
	pathToLoad += C_PATHTOCAMB;
	std::ifstream fileData(pathToLoad.c_str());
	std::vector< double > data_x;
	std::vector< double > data_xi0;
	std::vector< double > data_xi2;
	std::vector< double > data_xi4;
	std::cout << pathToLoad << std::endl;
	while (fileData) {
		double x;
		double xi0;
		double xi2;
		double xi4;
		fileData>>x>>xi0>>xi2>>xi4;
		if (fileData==0) break;
		
		data_x.push_back(x);
		data_xi0.push_back(xi0);
		data_xi2.push_back(xi2);
		data_xi4.push_back(xi4);
	}

	const double inverse_step_size = 1./data_x[0];
	const unsigned int nbBinCAMB = data_x.size();
	//// Set to zero the last pixel
	const double maxDistCAMB = data_x[nbBinCAMB-1];
	data_x.push_back(data_x[nbBinCAMB-1]+1.);
	data_xi0.push_back(0.);
	data_xi2.push_back(0.);
	data_xi4.push_back(0.);



	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max          =  200.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100;
	const double maxPow2      = max*max;


	//// Arrays for data
	double dataMu[nbBin][nbBinM][11];
	double data2D[nbBinX][nbBinY][11];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<8; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<8; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}


	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	std::vector<double> v_zzRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		//// Create the conversion table from redshift to distance
		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TGraph* gr_from_dist_to_z = cosmo->create_TGraph_to_convert_distance_to_redshift(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand (seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		v_zzRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_zzRandQSO[i] = gr_from_dist_to_z->Eval(v_rrRandQSO[i]);
		}

		delete gr_from_dist_to_z;
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];

		
		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];
			double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
				z2 = v_zzRandQSO[q];
			}

			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;
			
			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya );
			const double zQSO = z2;
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;

				if (fabs(distP) >= max) continue;

				const double distTotPow2 = distTransQsoLyaPow2 + distP*distP;
				const double distTot    = sqrt(distTotPow2);
				const double distTotMetal = sqrt(distTransQsoLyaPow2 + (v_r_metal[f][i]-rQSO)*(v_r_metal[f][i]-rQSO) );

				const double w   = v_w__[f][i];
				unsigned int idxBinCAMB = nbBinCAMB;
				if (distTotMetal<maxDistCAMB) idxBinCAMB = int(distTotMetal*inverse_step_size);
				const double growth_factor_pow_2 = from_z_to_growth_factor_pow_2[ int(0.5*(zQSO+v_z__[f][i])*inverse_growth_factor_step_size) ];
				const double wxi0 = w*growth_factor_pow_2*data_xi0[idxBinCAMB];
				const double wxi2 = w*growth_factor_pow_2*data_xi2[idxBinCAMB];
				const double wxi4 = w*growth_factor_pow_2*data_xi4[idxBinCAMB];

				if (distTot < max) {
					const double mu         = distP/distTot;
					const unsigned int idx  = int(distTot);
					const unsigned int idxM = int((mu+1.)*50.);
	
					dataMu[idx][idxM][0] += wxi0;
					dataMu[idx][idxM][1] += wxi2;
					dataMu[idx][idxM][2] += wxi4;
					dataMu[idx][idxM][3] += w*distTot;
					dataMu[idx][idxM][4] += w*mu;
					dataMu[idx][idxM][5] += w*(zQSO+v_z__[f][i]);
					dataMu[idx][idxM][6] += w;
					dataMu[idx][idxM][7] ++;
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rParralBinIdx = int(distP+max);
				data2D[rPerpBinIdx][rParralBinIdx][0] += wxi0;
				data2D[rPerpBinIdx][rParralBinIdx][1] += wxi2;
				data2D[rPerpBinIdx][rParralBinIdx][2] += wxi4;
				data2D[rPerpBinIdx][rParralBinIdx][3] += w*distTransQsoLya;
				data2D[rPerpBinIdx][rParralBinIdx][4] += w*distP;
				data2D[rPerpBinIdx][rParralBinIdx][5] += w*(zQSO+v_z__[f][i]);
				data2D[rPerpBinIdx][rParralBinIdx][6] += w;
				data2D[rPerpBinIdx][rParralBinIdx][7] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	//// Set the prefix of different forest and QSOs
	std::string prefix1 = forest__;
	prefix1 += "_";
	prefix1 += QSO__;
	prefix1 += "_";
	prefix1 += lambdaFrMetalName;

	std::ofstream fFile;
	std::string pathToSave;


	///// Save the 2D cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_Metals_Models_2D_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {

			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3];
			fFile << " " << data2D[i][j][4];
			fFile << " " << data2D[i][j][5]/2.;
			fFile << " " << data2D[i][j][6];
			fFile << " " << data2D[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();

	///// Save the Mu cross-correlation
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_Metals_Models_Mu_";
	pathToSave += prefix1;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3];
			fFile << " " << dataMu[i][j][4];
			fFile << " " << dataMu[i][j][5]/2.;
			fFile << " " << dataMu[i][j][6];
			fFile << " " << dataMu[i][j][7];
			fFile << std::endl;
		}
	}
	fFile.close();

	return;
}
void Correlation::xi_delta_QSO_MockJMc_Wick_T1_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_Wick_T1_with_wi1D_array ------\n" << std::endl;


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_region_Map__.clear();
	v_d__.clear();
	v_z__.clear();
	v_lRF__.clear();
	v_residual_delta_vs_lRF__.clear();
	v_residual_delta_vs_lObs__.clear();
	v_nb__.clear();
	v_zzQ1__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	std::vector<double> v_zzRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		//// Create the conversion table from redshift to distance
		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TGraph* gr_from_dist_to_z = cosmo->create_TGraph_to_convert_distance_to_redshift(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand (seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		v_zzRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_zzRandQSO[i] = gr_from_dist_to_z->Eval(v_rrRandQSO[i]);
		}

		delete gr_from_dist_to_z;
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];

		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];
			double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
				z2 = v_zzRandQSO[q];
			}

			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya*fromValToIdx );
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;
				if (fabs(distP) >= max) continue;
				const double w           = v_w__[f][i];
	
				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );

				nb[globalBin] ++;
				weight[globalBin] += w;

				const unsigned int pixel_lObs_index_1 = 1+int(v_lObs__[f][i]-lambdaObsMin__);
				cov[globalBin][globalBin] += w*w*xidd1D[pixel_lObs_index_1][pixel_lObs_index_1];
				nb_all[globalBin][globalBin] ++;
                        	we_all[globalBin][globalBin] += w*w;
			}
		}
	}




	const std::string diagram = "T1_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}
void Correlation::xi_delta_QSO_MockJMc_Wick_T12_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_Wick_T12_with_wi1D_array ------\n" << std::endl;


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_region_Map__.clear();
	v_d__.clear();
	v_z__.clear();
	v_lRF__.clear();
	v_residual_delta_vs_lRF__.clear();
	v_residual_delta_vs_lObs__.clear();
	v_nb__.clear();
	v_zzQ1__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	std::vector<double> v_zzRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		//// Create the conversion table from redshift to distance
		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TGraph* gr_from_dist_to_z = cosmo->create_TGraph_to_convert_distance_to_redshift(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand (seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		v_zzRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_zzRandQSO[i] = gr_from_dist_to_z->Eval(v_rrRandQSO[i]);
		}

		delete gr_from_dist_to_z;
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];

		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}

		for (unsigned int q=0; q<nbQ1__; q++) {

			//// For T1 + T2 + T3 + T4
			std::vector<unsigned int> vect_pair_bin;
			std::vector<unsigned int> vect_pixel_index;

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];
			double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
				z2 = v_zzRandQSO[q];
			}

			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya*fromValToIdx );
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;
				if (fabs(distP) >= max) continue;
				const double w           = v_w__[f][i];
	
				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );

				nb[globalBin] ++;
				weight[globalBin] += w;

				const unsigned int pixel_lObs_index_1 = 1+int(v_lObs__[f][i]-lambdaObsMin__);
				cov[globalBin][globalBin] += w*w*xidd1D[pixel_lObs_index_1][pixel_lObs_index_1];
				nb_all[globalBin][globalBin] ++;
                        	we_all[globalBin][globalBin] += w*w;

				//// For T2
				vect_pair_bin.push_back(globalBin);
				vect_pixel_index.push_back(i);
			}


			//// Calculate the covariance
			const unsigned int nbPairs = vect_pair_bin.size();
			
			for (unsigned int i=0; i<nbPairs; i++) {

				const unsigned int global_bin_1  = vect_pair_bin[i];
				const unsigned int pixel_index_1 = vect_pixel_index[i];
				const double www1 = v_w__[f][pixel_index_1];
				const unsigned int pixel_lObs_index_1 = int(v_lObs__[f][pixel_index_1]-lambdaObsMin__);

				//// Covariance
				for (unsigned int j=0; j<nbPairs; j++) {
					if (i==j) continue;

					const unsigned int global_bin_2  = vect_pair_bin[j];
					const unsigned int pixel_index_2 = vect_pixel_index[j];
					const double www2 = v_w__[f][pixel_index_2];
					unsigned int pixel_lObs_index_2 = int(v_lObs__[f][pixel_index_2]-lambdaObsMin__);

					if (v_lObs__[f][pixel_index_1]>v_lObs__[f][pixel_index_2]) {
						cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2];
					}
					else cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_2+1][pixel_lObs_index_1];
					nb_all[global_bin_1][global_bin_2] ++;
					we_all[global_bin_1][global_bin_2] += www1*www2;
				}
			}

		}
	}




	const std::string diagram = "T12_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}
void Correlation::xi_delta_QSO_MockJMc_Wick_T123_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_Wick_T123_with_wi1D_array ------\n" << std::endl;


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_region_Map__.clear();
	v_d__.clear();
	v_z__.clear();
	v_lRF__.clear();
	v_residual_delta_vs_lRF__.clear();
	v_residual_delta_vs_lObs__.clear();
	v_nb__.clear();
	v_zzQ1__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	std::vector<double> v_zzRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		//// Create the conversion table from redshift to distance
		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TGraph* gr_from_dist_to_z = cosmo->create_TGraph_to_convert_distance_to_redshift(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand (seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		v_zzRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_zzRandQSO[i] = gr_from_dist_to_z->Eval(v_rrRandQSO[i]);
		}

		delete gr_from_dist_to_z;
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];

		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}

		//// For T1 + T2 + T3 + T4
		std::vector<unsigned int> vect_pair_bin;
		std::vector<unsigned int> vect_pixel_index;
		std::vector<unsigned int> vect_QSO_index;

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];
			double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
				z2 = v_zzRandQSO[q];
			}

			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya*fromValToIdx );
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;
				if (fabs(distP) >= max) continue;
				const double w           = v_w__[f][i];
	
				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );

				nb[globalBin] ++;
				weight[globalBin] += w;

				//// For T3
				vect_pair_bin.push_back(globalBin);
				vect_pixel_index.push_back(i);
				vect_QSO_index.push_back(q);
			}
		}

		//// Calculate the covariance
		const unsigned int nbPairs = vect_pair_bin.size();
		std::cout << f << " " << nbPairs << std::endl;

		for (unsigned int i=0; i<nbPairs; i++) {

			const unsigned int global_bin_1  = vect_pair_bin[i];
			const unsigned int pixel_index_1 = vect_pixel_index[i];
			const unsigned int QSO_index_1   = vect_QSO_index[i];
			const double www1 = v_w__[f][pixel_index_1];
			const unsigned int pixel_lObs_index_1 = int(v_lObs__[f][pixel_index_1]-lambdaObsMin__);

			//// Covariance
			for (unsigned int j=0; j<nbPairs; j++) {
				
				const unsigned int pixel_index_2 = vect_pixel_index[j];
				const unsigned int QSO_index_2   = vect_QSO_index[j];
				if (pixel_index_1!=pixel_index_2 && QSO_index_1!=QSO_index_2) continue;  // remove T4
				const unsigned int global_bin_2  = vect_pair_bin[j];
				const unsigned int pixel_lObs_index_2 = int(v_lObs__[f][pixel_index_2]-lambdaObsMin__);
				const double www2 = v_w__[f][pixel_index_2];
				
				if (pixel_index_1==pixel_index_2) cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2+1];
				else {
					if (v_lObs__[f][pixel_index_1]>v_lObs__[f][pixel_index_2]) {
						cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2];
					}
					else cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_2+1][pixel_lObs_index_1];
				}
				nb_all[global_bin_1][global_bin_2] ++;
				we_all[global_bin_1][global_bin_2] += www1*www2;
			}
		}


	}




	const std::string diagram = "T123_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}
void Correlation::xi_delta_QSO_MockJMc_Wick_T1234_with_wi1D_array(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_Wick_T1234_with_wi1D_array ------\n" << std::endl;


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = pathToSave__;
	pathToFile += "xi_1DlObs_2D_delta_delta_";
	pathToFile += forest__;
	pathToFile += ".txt";
	std::ifstream fileData(pathToFile.c_str());

	const unsigned int nbBins_xi1D  = int(lambdaObsMax__-lambdaObsMin__)+1;
	double** xidd1D = new double*[nbBins_xi1D];
	for (unsigned int i=0; i<nbBins_xi1D; i++) {
		xidd1D[i] = new double[nbBins_xi1D];
		for (unsigned int j=0; j<nbBins_xi1D; j++) {
			long double data=0.;
			fileData>>data;
			xidd1D[i][j] = data;
		}
	}
	fileData.close();

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);
	if (nbForest_==0.) return;

	///// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_idx__.clear();
	v_region_Map__.clear();
	v_d__.clear();
	v_z__.clear();
	v_lRF__.clear();
	v_residual_delta_vs_lRF__.clear();
	v_residual_delta_vs_lObs__.clear();
	v_nb__.clear();
	v_zzQ1__.clear();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// Constants:
	const double max                  = 200.;
	const unsigned int nbBin          = 50;
	const unsigned int nbBinX         = nbBin;
	const unsigned int nbBinY         = 2*nbBin;
	static const unsigned int nbBin2D = nbBinX*nbBinY;
	const double fromValToIdx         = nbBin/max;
	const double maxPow2              = max*max;

	//// Array with the covariance
	double nb[nbBin2D];
	double weight[nbBin2D];
	double** cov    = new double*[nbBin2D];
	double** nb_all = new double*[nbBin2D];
	double** we_all = new double*[nbBin2D];
	for (unsigned int i=0; i<nbBin2D; i++) {
		nb[i]     = 0.;
		weight[i] = 0.;
		cov[i]    = new double[nbBin2D];
		nb_all[i] = new double[nbBin2D];
		we_all[i] = new double[nbBin2D];
		for (unsigned int j=0; j<nbBin2D; j++) {
			cov[i][j]    = 0.;
			nb_all[i][j] = 0.;
			we_all[i][j] = 0.;
		}
	}


	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	std::vector<double> v_rrRandQSO;
	std::vector<double> v_zzRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		//// Create the conversion table from redshift to distance
		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TGraph* gr_from_dist_to_z = cosmo->create_TGraph_to_convert_distance_to_redshift(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand (seed_for_random_position__);

		//// Forest
		v_raRandForest.resize(nbForest_,0.);
		v_deRandForest.resize(nbForest_,0.);
		for (unsigned int i=0; i<nbForest_; i++) {
			v_raRandForest[i] = v_ra__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandForest[i] = v_de__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}

		//// QSO
		v_raRandQSO.resize(nbQ1__,0.);
		v_deRandQSO.resize(nbQ1__,0.);
		v_rrRandQSO.resize(nbQ1__,0.);
		v_zzRandQSO.resize(nbQ1__,0.);
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrRandQSO[i] = v_rrQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_zzRandQSO[i] = gr_from_dist_to_z->Eval(v_rrRandQSO[i]);
		}

		delete gr_from_dist_to_z;
	}

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get attributs of forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double x = v_ra__[f];
		const double y = v_de__[f];

		double x1 = x;
		double y1 = y;
		//// If random position in the cell
		if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
			x1 = v_raRandForest[f];
		 	y1 = v_deRandForest[f];
		}

		//// For T1 + T2 + T3 + T4
		std::vector<unsigned int> vect_pair_bin;
		std::vector<unsigned int> vect_pixel_index;

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			double r2 = v_rrQ1__[q];
			double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
				r2 = v_rrRandQSO[q];
				z2 = v_zzRandQSO[q];
			}

			//// distance of QSO
			const double rQSO = r2;
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya*fromValToIdx );
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;
				if (fabs(distP) >= max) continue;
				const double w           = v_w__[f][i];
	
				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int globalBin = rPerpBinIdx*nbBinY+int( (distP+max)*fromValToIdx );

				nb[globalBin] ++;
				weight[globalBin] += w;

				//// For T3
				vect_pair_bin.push_back(globalBin);
				vect_pixel_index.push_back(i);
			}
		}

		//// Calculate the covariance
		const unsigned int nbPairs = vect_pair_bin.size();
		std::cout << f << " " << nbPairs << std::endl;

		for (unsigned int i=0; i<nbPairs; i++) {

			const unsigned int global_bin_1  = vect_pair_bin[i];
			const unsigned int pixel_index_1 = vect_pixel_index[i];
			const double www1 = v_w__[f][pixel_index_1];
			const unsigned int pixel_lObs_index_1 = int(v_lObs__[f][pixel_index_1]-lambdaObsMin__);

			//// Covariance
			for (unsigned int j=0; j<nbPairs; j++) {
				
				const unsigned int pixel_index_2 = vect_pixel_index[j];
				const unsigned int global_bin_2  = vect_pair_bin[j];
				const unsigned int pixel_lObs_index_2 = int(v_lObs__[f][pixel_index_2]-lambdaObsMin__);
				const double www2 = v_w__[f][pixel_index_2];
				
				if (pixel_index_1==pixel_index_2) cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2+1];
				else {
					if (v_lObs__[f][pixel_index_1]>v_lObs__[f][pixel_index_2]) {
						cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_1+1][pixel_lObs_index_2];
					}
					else cov[global_bin_1][global_bin_2] += www1*www2*xidd1D[pixel_lObs_index_2+1][pixel_lObs_index_1];
				}
				nb_all[global_bin_1][global_bin_2] ++;
				we_all[global_bin_1][global_bin_2] += www1*www2;
			}
		}


	}




	const std::string diagram = "T1234_with_wi1D_array";

	std::cout << "\n  Saving\n" << std::endl;
	std::ofstream fFile;

	///// Save the 2D cross-correlation 
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	
	///// Set the values of data
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) cov[i][j] /= weight[i]*weight[j];
			fFile << cov[i][j] << " ";
		}
		fFile << std::endl;
	}
	fFile.close();

	///// Save the number 
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_number_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.)  nb_all[i][j] /= nb[i]*nb[j];
                        fFile << nb_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();


	///// Save the weight
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weight_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

        for (unsigned int i=0; i<nbBin2D; i++) {
                for (unsigned int j=0; j<nbBin2D; j++) {
			if (nb[i]>1. && nb[j]>1.) we_all[i][j] /= weight[i]*weight[j];
                        fFile << we_all[i][j] << " ";
                }
                fFile << std::endl;
        }
        fFile.close();

	
	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_numberCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << nb[i] << " ";
                fFile << std::endl;
        }
        fFile.close();

	///// Save the number for the correlation
	pathToSave = pathToSave__;
        pathToSave += "xi_delta_QSO_2D_Wick";
	pathToSave += "_";
	pathToSave += diagram;
	pathToSave += "_weightCorr_";
        pathToSave += forest__;
        pathToSave += "_";
        pathToSave += QSO__;
        pathToSave += ".txt";
        std::cout << "\n  " << pathToSave << std::endl;
        fFile.open(pathToSave.c_str());
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);

	for (unsigned int i=0; i<nbBin2D; i++) {
                fFile << weight[i] << " ";
                fFile << std::endl;
        }
        fFile.close();


	return;
}


























void Correlation::xi_QSO_QSO_MockJMc(unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_QSO_QSO_MockJMc ------\n" << std::endl;

	loadDataQ1();
	v_CosDeQ1__.clear();
	v_SinDeQ1__.clear();

	//// If doing with random
	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	std::vector<double> dataX(v_raQ1__);
	std::vector<double> dataY(v_deQ1__);
	std::vector<double> dataR(v_rrQ1__);
	std::vector<double> dataZ(v_zzQ1__);

	if (randomQSO) {
		std::cout << "  Random QSO " << bootIdx << std::endl;
		std::srand (bootIdx*10);

		//// number of qso
		unsigned int nbQSO = 0;
		while (nbQSO<nbQ1__) {

			//// Pick a QSO on the grid
			double x = sizeCell__*(int(  1.*sizeGridX__*rand()/RAND_MAX) + 0.5);
			double y = sizeCell__*(int(  1.*sizeGridY__*rand()/RAND_MAX) + 0.5);
			const double r = v_rrQ1__[nbQSO];
			const double z = v_zzQ1__[nbQSO];

			bool cont = false;

			//// Pick a new one if already picked
			for (unsigned int i=0; i<nbQSO; i++) {
				if (v_raQ1__[i]==x && v_deQ1__[i]==y && v_rrQ1__[i]==r) {
					cont = true;
					continue;
				}
			}
			if (cont) continue;

			v_raQ1__[nbQSO] = x;
			v_deQ1__[nbQSO] = y;
			v_rrQ1__[nbQSO] = r;
			v_zzQ1__[nbQSO] = z;
			nbQSO++;
		}

	}
	else {
		dataX.clear();
		dataY.clear();
		dataR.clear();
		dataZ.clear();
	}

	//// Constants:
	const double max          = 200.;
	const unsigned int nbBin  = int(max);
	const double maxPow2      = max*max;
	const unsigned int nbBinM = 50;

	//// Arrays for data
	double dataMu[nbBin][nbBinM][4];
	double data2D[nbBin][nbBin][4];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<4; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<4; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}


	///// Vectors of randomized positions in cell
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		//// Create the conversion table from redshift to distance
		Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
		TGraph* gr_from_dist_to_z = cosmo->create_TGraph_to_convert_distance_to_redshift(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
		delete cosmo;

		std::cout << "  Seed is = " << seed_for_random_position__ << std::endl;
		std::srand(seed_for_random_position__);

		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raQ1__[i] += sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deQ1__[i] += sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_rrQ1__[i] += sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_zzQ1__[i]  = gr_from_dist_to_z->Eval(v_rrQ1__[i]);

			if (randomQSO) {
				dataX[i]    += sizeCell__*(1.*rand()/RAND_MAX-0.5);
                        	dataY[i]    += sizeCell__*(1.*rand()/RAND_MAX-0.5);
				dataR[i]    += sizeCell__*(1.*rand()/RAND_MAX-0.5);
				dataZ[i]     = gr_from_dist_to_z->Eval(dataR[i]);
			}
		}
		delete gr_from_dist_to_z;
	}


	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int q1=0; q1<nbQ1__; q1++) {
		
		const double x1 = v_raQ1__[q1];
		const double y1 = v_deQ1__[q1];
		const double r1 = v_rrQ1__[q1];
		const double z1 = v_zzQ1__[q1];

		for (unsigned int q2=0; q2<q1; q2++) {

			//// s_paral
			const double rParaPow2 = (r1-v_rrQ1__[q2])*(r1-v_rrQ1__[q2]);
			if (rParaPow2>=maxPow2) continue;

			//// s_perp
			const double rPerpPow2 = (x1-v_raQ1__[q2])*(x1-v_raQ1__[q2]) + (y1-v_deQ1__[q2])*(y1-v_deQ1__[q2]);
			if (rPerpPow2>=maxPow2) continue;

			//// |s|
			const double distTotPow2 = rPerpPow2+rParaPow2;
			const double rPerp = sqrt(rPerpPow2);
			const double rPara = sqrt(rParaPow2);
			const unsigned int idxPerp = int( rPerp );
			const unsigned int idxPara = int( rPara );

			//// 1D
			if (distTotPow2<maxPow2) {
				const double distTot    = sqrt(distTotPow2);
				const double mu         = rPara/distTot;
				const unsigned int idx  = int(distTot);

				unsigned int idxM = 49;
				if (rPara != distTot) idxM = int(50.*mu);

				dataMu[idx][idxM][0] ++;
				dataMu[idx][idxM][1] += distTot;
				dataMu[idx][idxM][2] += mu;
				dataMu[idx][idxM][3] += z1+v_zzQ1__[q2];
			}

			//// 2D
			data2D[idxPerp][idxPara][0] ++;
			data2D[idxPerp][idxPara][1] +=rPerp;
			data2D[idxPerp][idxPara][2] +=rPara;
			data2D[idxPerp][idxPara][3] +=z1+v_zzQ1__[q2];
		}
	}
	
	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave;
	long double sumOne = 0.;
	long double meanZ = 0.;

	//// 2D
	pathToSave = pathToSave__;
	pathToSave += "xi_QSO_QSO_2D_";
	pathToSave += QSO__;
	if (randomQSO) {
		pathToSave += "_RR_";
		pathToSave += strBootIdx;
	}
	else {
		pathToSave += "_DD";
	}
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	sumOne = 0.;
	meanZ  = 0.;

	//// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			fFile << data2D[i][j][0];
			fFile << " " << data2D[i][j][1];
			fFile << " " << data2D[i][j][2];
			fFile << " " << data2D[i][j][3]/2.;
			fFile << std::endl;

			sumOne += data2D[i][j][0];
			meanZ  += data2D[i][j][3]/2.;
		}
	}
	fFile.close();

	//// Get mean redshift of pairs
	std::cout << "  < z >           = " << meanZ/sumOne << " +- " << 0. << std::endl;
	std::cout << "  number of pairs = " << sumOne          << std::endl;


	// Mu
	pathToSave = pathToSave__;
	pathToSave += "xi_QSO_QSO_Mu_";
	pathToSave += QSO__;
	if (randomQSO) {
		pathToSave += "_RR_";
		pathToSave += strBootIdx;
	}
	else {
		pathToSave += "_DD";
	}
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	//// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			fFile << dataMu[i][j][0];
			fFile << " " << dataMu[i][j][1];
			fFile << " " << dataMu[i][j][2];
			fFile << " " << dataMu[i][j][3]/2.;
			fFile << std::endl;
		}
	}
	fFile.close();


	//// Doing the correlation random - data
	if (randomQSO) {
	
		//// Reset values to zero
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBinM; j++) {
				for (unsigned int k=0; k<4; k++) {
					dataMu[i][j][k] = 0.;
				}
			}
		}
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBin; j++) {
				for (unsigned int k=0; k<4; k++) {
					data2D[i][j][k] = 0.;
				}
			}
		}

		std::cout << "\n  Starting\n" << std::endl;
	
		//// get a data
		for (unsigned int q1=0; q1<nbQ1__; q1++) {
			const double x1 = dataX[q1];
			const double y1 = dataY[q1];
			const double r1 = dataR[q1];
			const double z1 = dataZ[q1];
			
			//// Get a random
			for (unsigned int q2=0; q2<nbQ1__; q2++) {

				//// s_paral
				const double rParaPow2 = (r1-v_rrQ1__[q2])*(r1-v_rrQ1__[q2]);
                                if (rParaPow2>=maxPow2) continue;

				//// s_perp
				const double rPerpPow2 = (x1-v_raQ1__[q2])*(x1-v_raQ1__[q2]) + (y1-v_deQ1__[q2])*(y1-v_deQ1__[q2]);
				if (rPerpPow2>=maxPow2) continue;
				
				//// |s|
				const double distTotPow2 = rPerpPow2+rParaPow2;
				const double rPerp = sqrt(rPerpPow2);
				const double rPara = sqrt(rParaPow2);
				const unsigned int idxPerp = int( rPerp );
				const unsigned int idxPara = int( rPara );
				
				//// 1D
				if (distTotPow2<maxPow2) {
					const double distTot    = sqrt(distTotPow2);
					const double mu         = rPara/distTot;
					const unsigned int idx  = int(distTot);
					
					unsigned int idxM = 49;
					if (rPara != distTot) idxM = int(50.*mu);
					
					dataMu[idx][idxM][0] ++;
					dataMu[idx][idxM][1] += distTot;
					dataMu[idx][idxM][2] += mu;
					dataMu[idx][idxM][3] += z1+v_zzQ1__[q2];
				}
				
				//// 2D
				data2D[idxPerp][idxPara][0] ++;
				data2D[idxPerp][idxPara][1] +=rPerp;
				data2D[idxPerp][idxPara][2] +=rPara;
				data2D[idxPerp][idxPara][3] +=z1+v_zzQ1__[q2];
			}
		}

		std::cout << "\n  Saving\n" << std::endl;

		//// 2D
		pathToSave = pathToSave__;
		pathToSave += "xi_QSO_QSO_2D_";
		pathToSave += QSO__;
		pathToSave += "_DR_";
		pathToSave += strBootIdx;
		pathToSave += ".txt";
		std::cout << "\n  " << pathToSave << std::endl;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);

		sumOne = 0.;
		meanZ  = 0.;

		//// Set the values of data
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBin; j++) {
				fFile << data2D[i][j][0];
				fFile << " " << data2D[i][j][1];
				fFile << " " << data2D[i][j][2];
				fFile << " " << data2D[i][j][3]/2.;
				fFile << std::endl;

				sumOne += data2D[i][j][0];
				meanZ  += data2D[i][j][3]/2.;
			}
		}
		fFile.close();

		//// Get mean redshift of pairs
		std::cout << "  < z >           = " << meanZ/sumOne << " +- " << 0. << std::endl;
		std::cout << "  number of pairs = " << sumOne          << std::endl;

		// Mu
		pathToSave = pathToSave__;
		pathToSave += "xi_QSO_QSO_Mu_";
		pathToSave += QSO__;
		pathToSave += "_DR_";
		pathToSave += strBootIdx;
		pathToSave += ".txt";
		std::cout << "\n  " << pathToSave << std::endl;
		fFile.open(pathToSave.c_str());
		fFile << std::scientific;
		fFile.precision(std::numeric_limits<double>::digits10);

		//// Set the values of data
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<nbBinM; j++) {
				fFile << dataMu[i][j][0];
				fFile << " " << dataMu[i][j][1];
				fFile << " " << dataMu[i][j][2];
				fFile << " " << dataMu[i][j][3]/2.;
				fFile << std::endl;
			}
		}
		fFile.close();
	}

	return;
}























































// ---------------------------------------------------------------------
//
//		Load data
//
// ---------------------------------------------------------------------
void Correlation::loadDataQ1(void) {

	//// Starting
	std::cout << "\n\n  ------ load Data Q1 ------" << std::endl;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	bool selectWindowRedshift = true;
	double array[2] = {0.};
	// Get min_z, max_z according to pixels positions
	if (correlation_type__=="q_q" || correlation_type__=="q_q2" || mockBox__) {
		selectWindowRedshift = false;
	}
	else cosmo->FindMinMaxRedsift(maxCorrelation__, lambdaObsMin__, lambdaObsMax__, lambdaRFLine__, array);
	const double minRedshiftQSO = array[0];
	const double maxRedshiftQSO = array[1];

	delete cosmo;

	//// Variables for FITS
	const TString TSfitsnameSpec = pathQ1__;
	std::cout << "  pathToFits = " << pathQ1__ << std::endl;
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	if ( nbQ1__ == 0 || nbQ1__ > nrows) nbQ1__ = nrows;

	std::cout << "  number of        Q1    = " << nrows << std::endl;
	std::cout << "  number of loaded Q1    = " << nbQ1__ << std::endl;

	//// Load data
	for (unsigned int i=0; i<nbQ1__; i++) {

		//// Variables for data in FITS
		double ra = 0.;
		double de = 0.;
		double zz = 0.;

		fits_read_col(fitsptrSpec,TDOUBLE, 1,i+1,1,1,NULL,&ra, NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 2,i+1,1,1,NULL,&de, NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 3,i+1,1,1,NULL,&zz, NULL,&sta);

		if (zz<=0.) {
			//std::cout << "  Correlation::loadDataQ1::  ERROR:  zz<=0. , zz = " << zz << std::endl;
			continue;
		}
		if ( !mockBox__ && ra==0. && de==0. ) {
			std::cout << "  Correlation::loadDataQ1::  ERROR:  ra==0. && de==0. , RA = " << ra << " , Dec = " << de << " , idx = " << i << std::endl;
			continue;
		}
		if ( selectWindowRedshift_bool__ && selectWindowRedshift && (zz<minRedshiftQSO || zz>maxRedshiftQSO) ) continue;

		//// If not dealing with Jean-Marc's simulations
		if (!mockJMC__ && !mockBox__) {
			ra *= C_DEGTORAD;
			de *= C_DEGTORAD;
		}

		//// Store the data
		v_raQ1__.push_back(ra);
		v_deQ1__.push_back(de);
		v_CosDeQ1__.push_back(cos(de));
		v_SinDeQ1__.push_back(sin(de));
		v_zzQ1__.push_back(zz);
		if (mockBox__) v_rrQ1__.push_back( zz );
		else v_rrQ1__.push_back( hConvertRedshDist->Interpolate(zz) );
		distMinQ1__ = std::min(distMinQ1__,hConvertRedshDist->Interpolate(zz));
	}

	nbQ1__ = v_raQ1__.size();
	std::cout << "  number of good   Q1    = " << nbQ1__ << std::endl;

	delete hConvertRedshDist;

	return;
}
void Correlation::loadDataQ2(void) {

	//// Starting
	std::cout << "\n\n  ------ load Data Q2 ------" << std::endl;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
	delete cosmo;

	//// Variables for FITS
	const TString TSfitsnameSpec = pathQ2__;
	std::cout << "  pathToFits = " << pathQ2__ << std::endl;
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	if ( nbQ2__ == 0 || nbQ2__ > nrows) nbQ2__ = nrows;

	std::cout << "  number of        Q2    = " << nrows << std::endl;
	std::cout << "  number of loaded Q2    = " << nbQ2__ << std::endl;

	//// Load data
	for (unsigned int i=0; i<nbQ2__; i++) {

		//// Variables for data in FITS
		double ra = 0.;
		double de = 0.;
		double zz = 0.;

		fits_read_col(fitsptrSpec,TDOUBLE, 1,i+1,1,1,NULL,&ra, NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 2,i+1,1,1,NULL,&de, NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 3,i+1,1,1,NULL,&zz, NULL,&sta);

		ra *= C_DEGTORAD;
		de *= C_DEGTORAD;

		//// Store the data
		v_raQ2__.push_back(ra);
		v_deQ2__.push_back(de);
		v_CosDeQ2__.push_back(cos(de));
		v_SinDeQ2__.push_back(sin(de));
		v_zzQ2__.push_back(zz);
		v_rrQ2__.push_back( hConvertRedshDist->Interpolate(zz) );
		distMinQ2__ = std::min(distMinQ2__,hConvertRedshDist->Interpolate(zz));
	}

	nbQ2__ = v_raQ2__.size();
	std::cout << "  number of good   Q2    = " << nbQ2__ << std::endl;

	delete hConvertRedshDist;

	return;
}
void Correlation::loadDataForest(std::string pathToFits,unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n  ------ load Data Forest ------" << std::endl;
	if (replace_data_by_gaussian_random__) std::cout << "  Replacing data by random field" << std::endl;

	//TRandom3* rand = new Trandom3();
	//std::cout << rand->Gaus(0., 1. ); << std::endl;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	// Get the minimal distance of pixels
	distMinPixel__ = cosmo->GetMinDistPixels(lambdaObsMin__, lambdaRFLine__);

	delete cosmo;

	//// Variables for FITS
	const TString TSfitsnameSpec = pathToFits;
	std::cout << "  pathToFits = " << pathToFits << std::endl;
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	if ( nbForest_ == 0 || nbForest_ > nrows) nbForest_ = nrows;

	std::cout << "  number of        forest = " << nrows << std::endl;
	std::cout << "  number of loaded forest = " << nbForest_ << std::endl;

	v_region_Map__.assign(nrows, 100000);
	//// If doing sub-sampling, geting the array to know if the forest is in the given region
	if (doBootstraps__) {

		//// For the cross-correlation
		std::string pathToLoad = pathToSave__;
		if (correlation_type__ == "f_f") {
			pathToLoad += "xi_A_delta_delta_map_";
			pathToLoad += forest__;
			pathToLoad += ".txt";
		}
		else if (correlation_type__ == "f_f2") {
			pathToLoad += "xi_A_delta_delta2_map_";
			pathToLoad += forest__;
			pathToLoad += "_";
			pathToLoad += forest2__;
			pathToLoad += ".txt";
		}
		else if (correlation_type__ == "q_f") {
			pathToLoad += "xi_delta_QSO_map_";
			pathToLoad += forest__;
			pathToLoad += "_";
			pathToLoad += QSO__;
			pathToLoad += ".txt";
		}

		LymanForest* lymanForestObject = new LymanForest(pathToLoad, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		lymanForestObject->GetRegionArray(v_region_Map__);
		lymanForestObject->PrintRegionDetail(bootIdx);
		delete lymanForestObject;
		std::cout << "\n\n" << std::endl;
	}

	long unsigned int nbCutted[6] = {0}; 
	long double meanDelta[3] = {0.};
	long double meanDelta_Nicolas[7] = {0.};

	// Setup the root file where to store the data
	const TString tPathToSave = "../run/test.root";
	TFile* storeFile = new TFile(tPathToSave,"RECREATE","delta");
	storeFile->cd();
	// TH1D
	TH1D* h_delta     = new TH1D("h_delta","",20000,-5000.,5000.);
	TH1D* h_delta_projected     = new TH1D("h_delta_projected","",20000,-5000.,5000.);
	TH1D* h_meandelta = new TH1D("h_meandelta","",1000,-10.,10.);
	TH1D* h_fluxDLA   = new TH1D("h_fluxDLA","",10000,-1.,2.);
	// TProfile
	TProfile* tp_flux_vs_lambdaRF  = new TProfile("tp_flux_vs_lambdaRF","",(int)(20+lambdaRFMax__-lambdaRFMin__),lambdaRFMin__-10.,lambdaRFMax__+10.);
	TProfile* tp_flux_vs_lambdaOBS = new TProfile("tp_flux_vs_lambdaOBS","",(int)(20+lambdaObsMax__-lambdaObsMin__),lambdaObsMin__-10.,lambdaObsMax__+10.);
	TProfile* tp_delta_vs_lambdaRF  = new TProfile("tp_delta_vs_lambdaRF","",(int)(20+lambdaRFMax__-lambdaRFMin__),lambdaRFMin__-10.,lambdaRFMax__+10.);
	TProfile* tp_delta_vs_lambdaOBS = new TProfile("tp_delta_vs_lambdaOBS","",(int)(20+lambdaObsMax__-lambdaObsMin__),lambdaObsMin__-10.,lambdaObsMax__+10.);
	TProfile* tp_delta_vs_z = new TProfile("tp_delta_vs_z","",1000,1.7,6.7);
	TProfile* tp_delta_vs_r = new TProfile("tp_delta_vs_r","",1000,3000.,6000.);
	TProfile* tp_delta_projected_vs_lambdaRF  = new TProfile("tp_delta_projected_vs_lambdaRF","",(int)(20+lambdaRFMax__-lambdaRFMin__),lambdaRFMin__-10.,lambdaRFMax__+10.);
	TProfile* tp_delta_projected_vs_lambdaOBS = new TProfile("tp_delta_projected_vs_lambdaOBS","",(int)(20+lambdaObsMax__-lambdaObsMin__),lambdaObsMin__-10.,lambdaObsMax__+10.);
	TProfile* tp_delta_projected_vs_z = new TProfile("tp_delta_projected_vs_z","",1000,1.7,6.7);
	TProfile* tp_delta_projected_vs_r = new TProfile("tp_delta_projected_vs_r","",1000,3000.,6000.);

	//// Load data
	for (unsigned int i=0; i<nbForest_; i++) {

		if (doBootstraps__ && v_region_Map__[i]!=bootIdx && (correlation_type__=="f_f2" || correlation_type__=="q_f") ) continue;

		//// Variables for data in FITS
		double ra, de, zz, alpha, beta, meanForestLambdaRF;
		int bit;
		double LAMBDA_OBS[nbBinRFMax__];
		double NORM_FLUX[nbBinRFMax__];
		double NORM_FLUX_IVAR[nbBinRFMax__];
		double FLUX_DLA[nbBinRFMax__];
		double DELTA[nbBinRFMax__];
		double DELTA_WEIGHT[nbBinRFMax__];

		fits_read_col(fitsptrSpec,TDOUBLE, 4,i+1,1,1,NULL,&ra,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 5,i+1,1,1,NULL,&de,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 6,i+1,1,1,NULL,&zz,   NULL,&sta);
		fits_read_col(fitsptrSpec,TINT,    7,i+1,1,1,NULL,&bit,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 8,i+1,1,1,NULL,&meanForestLambdaRF,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 9,i+1,1,1,NULL, &alpha,               NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 10,i+1,1,1,NULL,&beta,                NULL,&sta);

		fits_read_col(fitsptrSpec,TDOUBLE, 11,i+1,1,nbBinRFMax__,NULL, &LAMBDA_OBS,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 12,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX_IVAR,  NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 14,i+1,1,nbBinRFMax__,NULL, &FLUX_DLA,        NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMax__,NULL, &DELTA,           NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 17,i+1,1,nbBinRFMax__,NULL, &DELTA_WEIGHT,    NULL,&sta);

		if ( ra==0. && de==0. ) {
			//std::cout << "  Correlation::loadDataForest::  ERROR:  ra==0. && de==0." << std::endl;
			nbCutted[5] ++;
			continue;
		}
		if (cutNotFittedSpectra__ && !mocksNoNoiseNoCont && ( (alpha == alphaStart__ && beta == betaStart__) || (fabs(alpha)>=maxAlpha__-0.5) || (fabs(beta)>=maxBeta__-0.05) ) ) {
			if (alpha == alphaStart__ && beta == betaStart__) nbCutted[0] ++;
			else if (fabs(alpha)>=maxAlpha__-0.5) nbCutted[1] ++;
			else if (fabs(beta)>=maxBeta__-0.05) nbCutted[2] ++;
			continue;
		}

		//// If a reobs
		if (bit==isReobsFlag__) continue;

		const double oneOverOnePlusZ = 1./(1.+zz);
		bool templateHasNegative = false;
		long double tmp_meanDelta[3] = {0.};
		long double meanNeeded[5] = {0.};
		//long double SNR[2] = {0.};
		std::vector< double > v_tmp_r;
		std::vector< double > v_tmp_d;
		std::vector< double > v_tmp_w;
		std::vector< double > v_tmp_z;
		std::vector< double > v_tmp_lRF;
		std::vector< double > v_tmp_lObs;

		for (unsigned int j=0; j<nbBinRFMax__; j++) {

			if (DELTA_WEIGHT[j]<=0. || NORM_FLUX_IVAR[j]<=0. || FLUX_DLA[j]<C_DLACORR || LAMBDA_OBS[j]<lambdaObsMin__ || LAMBDA_OBS[j]>=lambdaObsMax__) continue;

			const double lambdaRF = LAMBDA_OBS[j]*oneOverOnePlusZ;
			if (lambdaRF<lambdaRFMin__ || lambdaRF>=lambdaRFMax__) continue;

			//// Remove veto lines
			bool isLine = false;
			if (doVetoLines__) {
				for (unsigned int k=0; k<nbVetoLines__; k++) {
					 if (LAMBDA_OBS[j]>=vetoLine__[2*k] && LAMBDA_OBS[j]<vetoLine__[2*k+1]) {
						isLine = true;
						break;
					}
				}
			}
			if (isLine) continue;

			const double zi = LAMBDA_OBS[j]/lambdaRFLine__ -1.;

			v_tmp_r.push_back(hConvertRedshDist->Interpolate(zi));
			v_tmp_d.push_back(DELTA[j]);
			v_tmp_w.push_back(DELTA_WEIGHT[j]);
			v_tmp_z.push_back(zi);
			v_tmp_lRF.push_back(lambdaRF);
			v_tmp_lObs.push_back(LAMBDA_OBS[j]);

			//// Get if the template is even one time negative
			if ( alpha+beta*(lambdaRF-meanForestLambdaRF) <= 0.) templateHasNegative = true;

			//// Get Nb Pixel in forest
			tmp_meanDelta[0] += DELTA_WEIGHT[j]*DELTA[j];
			tmp_meanDelta[1] += DELTA_WEIGHT[j];
			tmp_meanDelta[2] ++;

			//// SNR
			//SNR[0] += NORM_FLUX[j]*sqrt(NORM_FLUX_IVAR[j]);
			//SNR[1] ++;

			//// For Nicolas's estimator
			if (nicolasEstimator__) {
				meanNeeded[0] += DELTA_WEIGHT[j]*DELTA[j];
				meanNeeded[1] += DELTA_WEIGHT[j]*lambdaRF;
				meanNeeded[2] += DELTA_WEIGHT[j]*lambdaRF*lambdaRF;
				meanNeeded[3] += DELTA_WEIGHT[j]*lambdaRF*DELTA[j];
				meanNeeded[4] += DELTA_WEIGHT[j];
			}

			tp_flux_vs_lambdaRF->Fill(lambdaRF,NORM_FLUX[j],DELTA_WEIGHT[j]);
			tp_flux_vs_lambdaOBS->Fill(LAMBDA_OBS[j],NORM_FLUX[j],DELTA_WEIGHT[j]);
			h_fluxDLA->Fill(FLUX_DLA[j]);
		}

		const unsigned int tmp_nb = v_tmp_r.size();
		if (!mocksNoNoiseNoCont && cutNotFittedSpectra__ && (tmp_nb<=C_MIN_NB_PIXEL || templateHasNegative) ) {
			if (tmp_nb<=C_MIN_NB_PIXEL) nbCutted[3] ++;
			else if (templateHasNegative) nbCutted[4] ++;
			continue;
		}

		//// Get Nb Pixel in forest
		meanDelta[0] += tmp_meanDelta[0];
		meanDelta[1] += tmp_meanDelta[1];
		meanDelta[2] += tmp_meanDelta[2];

		//// For Nicolas's estimator
		const long double meanDelta   = meanNeeded[0]/meanNeeded[4];
		const long double meanLambda  = meanNeeded[1]/meanNeeded[4];
		const long double numerator   = meanNeeded[3]-meanLambda*meanNeeded[0];
		const long double denominator = meanNeeded[2]-meanLambda*meanLambda*meanNeeded[4];
		const long double coef        = numerator/denominator;
		for (unsigned int j=0; j<tmp_nb; j++) {

			h_delta->Fill(v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_lambdaRF->Fill(v_tmp_lRF[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_lambdaOBS->Fill(v_tmp_lObs[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_z->Fill(v_tmp_z[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_r->Fill(v_tmp_r[j], v_tmp_d[j],v_tmp_w[j]);
			
			if (nicolasEstimator__) v_tmp_d[j] -= meanDelta + (v_tmp_lRF[j]-meanLambda)*coef;

			h_delta_projected->Fill(v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_lambdaRF->Fill(v_tmp_lRF[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_lambdaOBS->Fill(v_tmp_lObs[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_z->Fill(v_tmp_z[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_r->Fill(v_tmp_r[j], v_tmp_d[j],v_tmp_w[j]);

			//// Get Nb Pixel in forest
			meanDelta_Nicolas[0] += v_tmp_w[j]*v_tmp_d[j];
			meanDelta_Nicolas[1] += v_tmp_w[j]*meanDelta;
			meanDelta_Nicolas[2] += v_tmp_w[j]*coef;
			meanDelta_Nicolas[3] += v_tmp_w[j]*(v_tmp_lRF[j]-meanLambda)*coef;
			meanDelta_Nicolas[4] += v_tmp_w[j]*(meanLambda-meanForestLambdaRF);
			meanDelta_Nicolas[5] += v_tmp_w[j];
			meanDelta_Nicolas[6] ++;
		}

		h_meandelta->Fill(meanDelta);
		

		//// If not dealing with Jean-Marc's simulations
		if (!mockJMC__) {
			ra *= C_DEGTORAD;
			de *= C_DEGTORAD;
		}

		v_ra__.push_back(ra);
		v_de__.push_back(de);
		v_CosDe__.push_back(cos(de));
		v_SinDe__.push_back(sin(de));
		v_zz__.push_back(zz);
		v_nbPixelDelta1__.push_back(v_tmp_r.size());
		v_idx__.push_back(i);
		v_r__.push_back(v_tmp_r);
		v_d__.push_back(v_tmp_d);
		v_w__.push_back(v_tmp_w);
		v_z__.push_back(v_tmp_z);
		v_lRF__.push_back(v_tmp_lRF);
		v_lObs__.push_back(v_tmp_lObs);
		v_residual_delta_vs_lRF__.push_back(v_tmp_lRF);
		v_residual_delta_vs_lObs__.push_back(v_tmp_lObs);
		std::vector< unsigned int > v_tmp_nb(tmp_nb,0);
		v_nb__.push_back(v_tmp_nb);
	}

	nbForest_ = v_ra__.size();
	std::cout << "  number of good forest   = " << nbForest_ << std::endl;
	for (unsigned int i=0; i<6; i++) std::cout << "  lost n°" << i << " = " << nbCutted[i] << std::endl;

	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;

	if (nicolasEstimator__) {
		std::cout << "  < delta_Nicolas >              = " << meanDelta_Nicolas[0]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < <delta> >                    = " << meanDelta_Nicolas[1]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < coef >                       = " << meanDelta_Nicolas[2]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < (l_i-l).coef >               = " << meanDelta_Nicolas[3]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < lambda > - < lambda >_before = " << meanDelta_Nicolas[4]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  sum(w_i)_Nicolas               = " << meanDelta_Nicolas[5] << std::endl;
		std::cout << "  nb pixel_Nicolas               = " << (long long unsigned int)meanDelta_Nicolas[6] << std::endl;
	}

	fits_close_file(fitsptrSpec,&sta);

	delete hConvertRedshDist;

	//// Replace the values of delta by the one of the stack of delta vs. lambda_RF 
	for ( unsigned int i=0; i<nbForest_; i++ ) {
		const unsigned int nbDelta = v_d__[i].size();
		for ( unsigned int j=0; j<nbDelta; j++ ) {
			if (nicolasEstimator__) {
				v_residual_delta_vs_lRF__[i][j]  = tp_delta_projected_vs_lambdaRF->Interpolate( v_lRF__[i][j]  );
				v_residual_delta_vs_lObs__[i][j] = tp_delta_projected_vs_lambdaOBS->Interpolate( v_lObs__[i][j]  );
			}
			else {
				v_residual_delta_vs_lRF__[i][j]  = tp_delta_vs_lambdaRF->Interpolate( v_lRF__[i][j]  );
                                v_residual_delta_vs_lObs__[i][j] = tp_delta_vs_lambdaOBS->Interpolate( v_lObs__[i][j]  );
			}
		}
	}

	if (showRootHistos__) {
		// TH1D
		R_plot1D(h_meandelta,"<\\delta>");
		R_plot1D(h_fluxDLA,"< f_{DLA} >");
		R_plot1D(h_delta,"\\delta");
		R_plot1D(h_delta_projected,"\\delta_{proj.}");
		// TProfile
		R_plot1D(tp_flux_vs_lambdaRF,"\\lambda_{R.F.}","flux");
		R_plot1D(tp_flux_vs_lambdaOBS,"\\lambda_{Obs.}","flux");
		R_plot1D(tp_delta_vs_lambdaRF,"\\lambda_{R.F.}","\\delta");
		R_plot1D(tp_delta_vs_lambdaOBS,"\\lambda_{Obs.}","\\delta");
		R_plot1D(tp_delta_vs_z,"z_{pixel}","\\delta");
		R_plot1D(tp_delta_vs_r,"r_{pixel}","\\delta");
		R_plot1D(tp_delta_projected_vs_lambdaRF,"\\lambda_{R.F.}","\\delta_{proj.}");
		R_plot1D(tp_delta_projected_vs_lambdaOBS,"\\lambda_{Obs.}","\\delta_{proj.}");
		R_plot1D(tp_delta_projected_vs_z,"z_{pixel}","\\delta_{proj.}");
		R_plot1D(tp_delta_projected_vs_r,"r_{pixel}","\\delta_{proj.}");
	}
	if (saveInRootFile__) {
		storeFile->Write();
		storeFile->Close();
	}
	delete storeFile;

	if (!(correlation_type__=="f_f") || !doBootstraps__) v_region_Map__.clear();

	return;
}
void Correlation::loadDataDelta2(int dataNeeded/*=100*/) {

	//// Starting
	std::cout << "\n\n  ------ load Data Forest 2 ------" << std::endl;

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	// Get the minimal distance of pixels
	distMinPixelDelta2__ = cosmo->GetMinDistPixels(lambdaObsMin__, lambdaRFLineDelta2__);

	delete cosmo;



	//// Variables for FITS
	const TString TSfitsnameSpec = pathDelta2__;
	std::cout << "  pathToFits = " << pathDelta2__ << std::endl;
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	if ( nbForest2__ == 0 || nbForest2__ > nrows) nbForest2__ = nrows;

	std::cout << "  number of        forest = " << nrows << std::endl;
	std::cout << "  number of loaded forest = " << nbForest2__ << std::endl;

	long unsigned int nbCutted[6] = {0};
	long double meanDelta[3] = {0.};
	long double meanDelta_Nicolas[7] = {0.};

	// Setup the root file where to store the data
	const TString tPathToSave = "../run/test.root";
	TFile* storeFile = new TFile(tPathToSave,"RECREATE","delta");
	storeFile->cd();
	// TH1D
	TH1D* h_delta     = new TH1D("h_delta","",20000,-5000.,5000.);
	TH1D* h_delta_projected     = new TH1D("h_delta_projected","",20000,-5000.,5000.);
	TH1D* h_meandelta = new TH1D("h_meandelta","",1000,-10.,10.);
	TH1D* h_fluxDLA   = new TH1D("h_fluxDLA","",10000,-1.,2.);
	// TProfile
	TProfile* tp_flux_vs_lambdaRF  = new TProfile("tp_flux_vs_lambdaRF","",(int)(20+lambdaRFMax__-lambdaRFMin__),lambdaRFMin__-10.,lambdaRFMax__+10.);
	TProfile* tp_flux_vs_lambdaOBS = new TProfile("tp_flux_vs_lambdaOBS","",(int)(20+lambdaObsMax__-lambdaObsMin__),lambdaObsMin__-10.,lambdaObsMax__+10.);
	TProfile* tp_delta_vs_lambdaRF  = new TProfile("tp_delta_vs_lambdaRF","",(int)(20+lambdaRFMax__-lambdaRFMin__),lambdaRFMin__-10.,lambdaRFMax__+10.);
	TProfile* tp_delta_vs_lambdaOBS = new TProfile("tp_delta_vs_lambdaOBS","",(int)(20+lambdaObsMax__-lambdaObsMin__),lambdaObsMin__-10.,lambdaObsMax__+10.);
	TProfile* tp_delta_vs_z = new TProfile("tp_delta_vs_z","",1000,1.7,6.7);
	TProfile* tp_delta_vs_r = new TProfile("tp_delta_vs_r","",1000,3000.,6000.);
	TProfile* tp_delta_projected_vs_lambdaRF  = new TProfile("tp_delta_projected_vs_lambdaRF","",(int)(20+lambdaRFMax__-lambdaRFMin__),lambdaRFMin__-10.,lambdaRFMax__+10.);
	TProfile* tp_delta_projected_vs_lambdaOBS = new TProfile("tp_delta_projected_vs_lambdaOBS","",(int)(20+lambdaObsMax__-lambdaObsMin__),lambdaObsMin__-10.,lambdaObsMax__+10.);
	TProfile* tp_delta_projected_vs_z = new TProfile("tp_delta_projected_vs_z","",1000,1.7,6.7);
	TProfile* tp_delta_projected_vs_r = new TProfile("tp_delta_projected_vs_r","",1000,3000.,6000.);

	//// Load data
	for (unsigned int i=0; i<nbForest2__; i++) {

		//// Variables for data in FITS
		double ra, de, zz, alpha, beta, meanForestLambdaRF;
		unsigned int bit;

		double LAMBDA_OBS[nbBinRFMaxDelta2__];
		double NORM_FLUX[nbBinRFMaxDelta2__];
		double NORM_FLUX_IVAR[nbBinRFMaxDelta2__];
		double FLUX_DLA[nbBinRFMaxDelta2__];
		double DELTA[nbBinRFMaxDelta2__];
		double DELTA_WEIGHT[nbBinRFMaxDelta2__];

		fits_read_col(fitsptrSpec,TDOUBLE, 4,i+1,1,1,NULL,&ra,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 5,i+1,1,1,NULL,&de,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 6,i+1,1,1,NULL,&zz,   NULL,&sta);
		fits_read_col(fitsptrSpec,TINT,    7,i+1,1,1,NULL,&bit,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 8,i+1,1,1,NULL,&meanForestLambdaRF,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 9,i+1,1,1,NULL,&alpha,               NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 10,i+1,1,1,NULL,&beta,                NULL,&sta);

		fits_read_col(fitsptrSpec,TDOUBLE, 11,i+1,1,nbBinRFMaxDelta2__,NULL, &LAMBDA_OBS,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 12,i+1,1,nbBinRFMaxDelta2__,NULL, &NORM_FLUX,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMaxDelta2__,NULL, &NORM_FLUX_IVAR,  NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 14,i+1,1,nbBinRFMaxDelta2__,NULL, &FLUX_DLA,        NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMaxDelta2__,NULL, &DELTA,           NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 17,i+1,1,nbBinRFMaxDelta2__,NULL, &DELTA_WEIGHT,    NULL,&sta);

		if ( ra==0. && de==0. ) {
			//std::cout << "  Correlation::loadDataForest2::  ERROR:  ra==0. && de==0." << std::endl;
			nbCutted[5] ++;
			continue;
		}
		if (!mocksNoNoiseNoCont && ( (alpha == alphaStart__ && beta == betaStart__) || (fabs(alpha)>=maxAlpha__-0.5) || (fabs(beta)>=maxBeta__-0.05) ) ) {
			if (alpha == alphaStart__ && beta == betaStart__) nbCutted[0] ++;
			else if (fabs(alpha)>=maxAlpha__-0.5) nbCutted[1] ++;
			else if (fabs(beta)>=maxBeta__-0.05) nbCutted[2] ++;
			continue;
		}

		bool templateHasNegative = false;
		long double meanNeeded[5] = {0.};
		long double tmp_meanDelta[3] = {0.};
		const double oneOverOnePlusZ = 1./(1.+zz);
		std::vector< double > v_tmp_r;
		std::vector< double > v_tmp_d;
		std::vector< double > v_tmp_w;
		std::vector< double > v_tmp_z;
		std::vector< double > v_tmp_lRF;
		std::vector< double > v_tmp_lObs;

		for (unsigned int j=0; j<nbBinRFMaxDelta2__; j++) {

			if (DELTA_WEIGHT[j]<=0. || NORM_FLUX_IVAR[j]<=0. || FLUX_DLA[j]<C_DLACORR || LAMBDA_OBS[j]<lambdaObsMin__ || LAMBDA_OBS[j]>=lambdaObsMax__) continue;

			const double lambdaRF = LAMBDA_OBS[j]*oneOverOnePlusZ;
			if (lambdaRF<lambdaRFMinDelta2__ || lambdaRF>=lambdaRFMaxDelta2__) continue;

			//// Remove veto lines
			bool isLine = false;
			if (doVetoLines__) {
				for (unsigned int k=0; k<nbVetoLines__; k++) {
					 if (LAMBDA_OBS[j]>=vetoLine__[2*k] && LAMBDA_OBS[j]<vetoLine__[2*k+1]) {
						isLine = true;
						break;
					}
				}
			}
			if (isLine) continue;

			//// Get if the template is even one time negative
			if ( alpha+beta*(lambdaRF-meanForestLambdaRF) <= 0.) templateHasNegative = true;

			const double zi = LAMBDA_OBS[j]/lambdaRFLineDelta2__ -1.;
			v_tmp_r.push_back(hConvertRedshDist->Interpolate(zi));
			v_tmp_d.push_back(DELTA[j]);
			v_tmp_w.push_back(DELTA_WEIGHT[j]);
			v_tmp_z.push_back(zi);
			v_tmp_lRF.push_back(lambdaRF);
			v_tmp_lObs.push_back(LAMBDA_OBS[j]);

			tmp_meanDelta[0] += DELTA_WEIGHT[j]*DELTA[j];
			tmp_meanDelta[1] += DELTA_WEIGHT[j];
			tmp_meanDelta[2] ++;

			//// For Nicolas's estimator
			if (nicolasEstimator__) {
				meanNeeded[0] += DELTA_WEIGHT[j]*DELTA[j];
				meanNeeded[1] += DELTA_WEIGHT[j]*lambdaRF;
				meanNeeded[2] += DELTA_WEIGHT[j]*lambdaRF*lambdaRF;
				meanNeeded[3] += DELTA_WEIGHT[j]*lambdaRF*DELTA[j];
				meanNeeded[4] += DELTA_WEIGHT[j];
			}

			tp_flux_vs_lambdaRF->Fill(lambdaRF,NORM_FLUX[j],DELTA_WEIGHT[j]);
			tp_flux_vs_lambdaOBS->Fill(LAMBDA_OBS[j],NORM_FLUX[j],DELTA_WEIGHT[j]);
			h_fluxDLA->Fill(FLUX_DLA[j]);
		}

		const unsigned int tmp_nb = v_tmp_r.size();
		if (tmp_nb<=C_MIN_NB_PIXEL || templateHasNegative ) {
			if (tmp_nb<=C_MIN_NB_PIXEL) nbCutted[3] ++;
			else if (templateHasNegative) nbCutted[4] ++;
			continue;
		}

		//// Get Nb Pixel in forest
		meanDelta[0] += tmp_meanDelta[0];
		meanDelta[1] += tmp_meanDelta[1];
		meanDelta[2] += tmp_meanDelta[2];

		//// For Nicolas's estimator
		const long double meanDelta   = meanNeeded[0]/meanNeeded[4];
		const long double meanLambda  = meanNeeded[1]/meanNeeded[4];
		const long double numerator   = meanNeeded[3]-meanLambda*meanNeeded[0];
		const long double denominator = meanNeeded[2]-meanLambda*meanLambda*meanNeeded[4];
		const long double coef        = numerator/denominator;
		for (unsigned int j=0; j<tmp_nb; j++) {

			h_delta->Fill(v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_lambdaRF->Fill(v_tmp_lRF[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_lambdaOBS->Fill(v_tmp_lObs[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_z->Fill(v_tmp_z[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_vs_r->Fill(v_tmp_r[j], v_tmp_d[j],v_tmp_w[j]);
			
			if (nicolasEstimator__) v_tmp_d[j] -= meanDelta + (v_tmp_lRF[j]-meanLambda)*coef;

			h_delta_projected->Fill(v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_lambdaRF->Fill(v_tmp_lRF[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_lambdaOBS->Fill(v_tmp_lObs[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_z->Fill(v_tmp_z[j], v_tmp_d[j],v_tmp_w[j]);
			tp_delta_projected_vs_r->Fill(v_tmp_r[j], v_tmp_d[j],v_tmp_w[j]);

			//// Get Nb Pixel in forest
			meanDelta_Nicolas[0] += v_tmp_w[j]*v_tmp_d[j];
			meanDelta_Nicolas[1] += v_tmp_w[j]*meanDelta;
			meanDelta_Nicolas[2] += v_tmp_w[j]*coef;
			meanDelta_Nicolas[3] += v_tmp_w[j]*(v_tmp_lRF[j]-meanLambda)*coef;
			meanDelta_Nicolas[4] += v_tmp_w[j]*(meanLambda-meanForestLambdaRF);
			meanDelta_Nicolas[5] += v_tmp_w[j];
			meanDelta_Nicolas[6] ++;
		}

		h_meandelta->Fill(meanDelta);

		if (!mockJMC__) {
			ra = ra*C_DEGTORAD;
			de = de*C_DEGTORAD;
		}

		v_raDelta2__.push_back(ra);
		v_deDelta2__.push_back(de);
		v_CosDeDelta2__.push_back(cos(de));
		v_SinDeDelta2__.push_back(sin(de));
		v_zzDelta2__.push_back(zz);
		v_nbPixelDelta2__.push_back(tmp_nb);
			
		v_rDelta2__.push_back(v_tmp_r);
		v_dDelta2__.push_back(v_tmp_d);
		v_wDelta2__.push_back(v_tmp_w);
		v_zDelta2__.push_back(v_tmp_z);
		v_lRFDelta2__.push_back(v_tmp_lRF);
		v_lObsDelta2__.push_back(v_tmp_lObs);
		v_residual2_delta_vs_lRF__.push_back(v_tmp_lRF);
		v_residual2_delta_vs_lObs__.push_back(v_tmp_lObs);
	}

	fits_close_file(fitsptrSpec,&sta);
	
	nbForest2__ = v_raDelta2__.size();
	std::cout << "  number of good forest   = " << nbForest2__ << std::endl;
	for (unsigned int i=0; i<6; i++) std::cout << "  lost n°" << i << " = " << nbCutted[i] << std::endl;

	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;

	if (nicolasEstimator__) {
		std::cout << "  < delta_Nicolas >              = " << meanDelta_Nicolas[0]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < <delta> >                    = " << meanDelta_Nicolas[1]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < coef >                       = " << meanDelta_Nicolas[2]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < (l_i-l).coef >               = " << meanDelta_Nicolas[3]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  < lambda > - < lambda >_before = " << meanDelta_Nicolas[4]/meanDelta_Nicolas[5] << std::endl;
		std::cout << "  sum(w_i)_Nicolas               = " << meanDelta_Nicolas[5] << std::endl;
		std::cout << "  nb pixel_Nicolas               = " << (long long unsigned int)meanDelta_Nicolas[6] << std::endl;
	}

	delete hConvertRedshDist;

	//// Replace the values of delta by the one of the stack of delta vs. lambda_RF 
	for ( unsigned int i=0; i<nbForest2__; i++ ) {
		const unsigned int nbDelta = v_dDelta2__[i].size();
		for ( unsigned int j=0; j<nbDelta; j++ ) {
			if (nicolasEstimator__) {
				v_residual2_delta_vs_lRF__[i][j]  = tp_delta_projected_vs_lambdaRF->Interpolate( v_lRFDelta2__[i][j]  );
				v_residual2_delta_vs_lObs__[i][j] = tp_delta_projected_vs_lambdaOBS->Interpolate( v_lObsDelta2__[i][j]  );
			}
			else {
				v_residual2_delta_vs_lRF__[i][j]  = tp_delta_vs_lambdaRF->Interpolate( v_lRFDelta2__[i][j]  );
                                v_residual2_delta_vs_lObs__[i][j] = tp_delta_vs_lambdaOBS->Interpolate( v_lObsDelta2__[i][j]  );
			}
		}
	}

	// TH1D
	R_plot1D(h_meandelta,"<\\delta>");
	R_plot1D(h_fluxDLA,"< f_{DLA} >");
	R_plot1D(h_delta,"\\delta");
	R_plot1D(h_delta_projected,"\\delta_{proj.}");
	// TProfile
	R_plot1D(tp_flux_vs_lambdaRF,"\\lambda_{R.F.}","flux");
	R_plot1D(tp_flux_vs_lambdaOBS,"\\lambda_{Obs.}","flux");
	R_plot1D(tp_delta_vs_lambdaRF,"\\lambda_{R.F.}","\\delta");
	R_plot1D(tp_delta_vs_lambdaOBS,"\\lambda_{Obs.}","\\delta");
	R_plot1D(tp_delta_vs_z,"z_{pixel}","\\delta");
	R_plot1D(tp_delta_vs_r,"r_{pixel}","\\delta");
	R_plot1D(tp_delta_projected_vs_lambdaRF,"\\lambda_{R.F.}","\\delta_{proj.}");
	R_plot1D(tp_delta_projected_vs_lambdaOBS,"\\lambda_{Obs.}","\\delta_{proj.}");
	R_plot1D(tp_delta_projected_vs_z,"z_{pixel}","\\delta_{proj.}");
	R_plot1D(tp_delta_projected_vs_r,"r_{pixel}","\\delta_{proj.}");
	if (saveInRootFile__) {
		storeFile->Write();
		storeFile->Close();
	}
	delete storeFile;

	return;
}
void Correlation::loadDataForest_Raw(unsigned int bootIdx/*=0*/) {

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);
	// Get the minimal distance of pixels
	distMinPixel__ = cosmo->GetMinDistPixels(lambdaObsMin__, lambdaRFLine__);
	delete cosmo;

	// TProfile
	TProfile* tp_flux_vs_lambdaRF  = new TProfile("tp_flux_vs_lambdaRF","",(int)(20+lambdaRFMax__-lambdaRFMin__),lambdaRFMin__-10.,lambdaRFMax__+10.);
	TProfile* tp_flux_vs_lambdaOBS = new TProfile("tp_flux_vs_lambdaOBS","",(int)(20+lambdaObsMax__-lambdaObsMin__),lambdaObsMin__-10.,lambdaObsMax__+10.);

	/// index of forest
	unsigned int forestIdx = 0;
	unsigned int nbGoodForest = 0;
	bool breakk = false;

	/// Get the list 
	FILE *fp;
	char path[PATH_MAX];
	std::string command = "ls " + pathToRaw__;
	std::vector< std::string > listFiles;
	fp = popen(command.c_str(), "r");
	while (fgets(path, PATH_MAX, fp) != NULL) listFiles.push_back(path);

	if (mock_version==0 || mock_version==2) {
		listFiles.clear();
		listFiles.push_back(pathToRaw__);
	}

	/// Get nb of files
	const unsigned int nbFiles = listFiles.size();
	long double meanDelta[3] = {0.};

	for (unsigned int fileIdx=0; fileIdx<nbFiles; fileIdx++) {

		if (breakk) break;

		std::string path = listFiles[fileIdx];
		path.erase(std::remove(path.begin(), path.end(), '\n'), path.end());
		path.erase(std::remove(path.begin(), path.end(), ' '), path.end());

		/// Get the file
		std::cout << "  file = " << path << std::endl;
		const TString TSfitsnameSpec = path;
		int sta = 0;
		fitsfile* fitsptrSpec;
		fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);

		/// Get the number of spectra
		int tmp_nbSpectra = 0;
		unsigned int nbSpectra = 0;
		fits_get_num_hdus(fitsptrSpec, &tmp_nbSpectra, &sta);
		nbSpectra = tmp_nbSpectra-1;
		//std::cout << "  nbSpectra = " << nbSpectra << std::endl;
	
		/// Set to the first HDU
		fits_movabs_hdu(fitsptrSpec, 2,  NULL, &sta);
		
		/// Get the data
		for (unsigned int f=0; f<nbSpectra; f++) {

			forestIdx ++;
	
			/// Move to next HDU
			if (f!=0) fits_movrel_hdu(fitsptrSpec, 1,  NULL, &sta);
	
			/// Get the number of pixels in this HDU
			long tmp_nbPixels = 0;
			fits_get_num_rows(fitsptrSpec, &tmp_nbPixels, &sta);
			const unsigned int tmp_nbPixels2 = tmp_nbPixels;
			if (tmp_nbPixels2<=C_MIN_NB_PIXEL) continue;

			/// Variable from old FITS
			float X = 0.;
			float Y = 0.;
			float Z = 0.;
			float LAMBDA_OBS[tmp_nbPixels2];
			float FLUX[tmp_nbPixels2];
			float CONTINUUM[tmp_nbPixels2];
			fits_read_key(fitsptrSpec,TFLOAT,"X",   &X,NULL,&sta);
			fits_read_key(fitsptrSpec,TFLOAT,"Y",   &Y,NULL,&sta);
			fits_read_key(fitsptrSpec,TFLOAT,"ZQSO",&Z,NULL,&sta);
			fits_read_col(fitsptrSpec,TFLOAT, 1,1,1,tmp_nbPixels2,NULL, &LAMBDA_OBS,NULL,&sta);
			if (mock_version==0) {
				fits_read_col(fitsptrSpec,TFLOAT, 4,1,1,tmp_nbPixels2,NULL, &CONTINUUM,NULL,&sta);
				fits_read_col(fitsptrSpec,TFLOAT, 5,1,1,tmp_nbPixels2,NULL, &FLUX,NULL,&sta);
			}
			else if (mock_version==1) fits_read_col(fitsptrSpec,TFLOAT, 6,1,1,tmp_nbPixels2,NULL, &FLUX,NULL,&sta);
			else if (mock_version==2 || mock_version==4 || mock_version==5) fits_read_col(fitsptrSpec,TFLOAT, 2,1,1,tmp_nbPixels2,NULL, &FLUX,NULL,&sta);
			else if (mock_version==3) {
				fits_read_col(fitsptrSpec,TFLOAT, 2,1,1,tmp_nbPixels2,NULL, &FLUX,NULL,&sta);
				fits_read_col(fitsptrSpec,TFLOAT, 5,1,1,tmp_nbPixels2,NULL, &CONTINUUM,NULL,&sta);
			}

	
			/// Variables for new FITS
			const double oneOverOnePlusZ = 1./(1.+Z);
			std::vector< double > v_tmp_r;
			std::vector< double > v_tmp_d;
			std::vector< double > v_tmp_w;
			std::vector< double > v_tmp_z;
			std::vector< double > v_tmp_lRF;
			std::vector< double > v_tmp_lObs;

			/// Get the mean over 3 pixels
			double mean_over3Pixels[5] = {0.};
	
			for (unsigned int p=0; p<tmp_nbPixels2; p++) {

				/// Remove 1/3 of pixels because too heavy
				//if ( (mock_version==2 || mock_version==4) && !meanOver3pixels__ && p%3!=0) continue;

				/// bad pixels
				if (LAMBDA_OBS[p]<=0. || std::isnan(FLUX[p]) ) continue;
				

				if (mock_version==1 || mock_version==3 ) LAMBDA_OBS[p] = pow(10.,LAMBDA_OBS[p]);
				const double lambdaRFd = LAMBDA_OBS[p]*oneOverOnePlusZ;
	
				/// Pixel outside working region and remove pixels because of CCD and Sky lines 
				if (lambdaRFd<lambdaRFMin__ || lambdaRFd>lambdaRFMax__ || LAMBDA_OBS[p]<lambdaObsMin__ || LAMBDA_OBS[p]>=lambdaObsMax__) continue;

				if (mock_version==0 || mock_version==3 ) {
					if (CONTINUUM[p]==0.) continue;
					else FLUX[p] /= CONTINUUM[p];
				}

				/// Get the mean over 3 pixels
				double zi = LAMBDA_OBS[p]/lambdaRFLine__ -1.;
				double wi = pow( (zi+1.)/onePlusZ0__, halfGama__);
				mean_over3Pixels[0] += wi*lambdaRFd;
				mean_over3Pixels[1] += wi*LAMBDA_OBS[p];
				mean_over3Pixels[2] += wi*FLUX[p];
				mean_over3Pixels[3] += wi;
				mean_over3Pixels[4] ++;

				if (!meanOver3pixels__) {
					v_tmp_r.push_back(hConvertRedshDist->Interpolate(zi));
					v_tmp_d.push_back(FLUX[p]);
					v_tmp_w.push_back(wi);
					v_tmp_z.push_back(zi);
					v_tmp_lRF.push_back(lambdaRFd);
					v_tmp_lObs.push_back(LAMBDA_OBS[p]);

					meanDelta[0] += wi*FLUX[p];
					meanDelta[1] += wi;
					meanDelta[2] ++;

					tp_flux_vs_lambdaRF->Fill( lambdaRFd, FLUX[p], wi );
					tp_flux_vs_lambdaOBS->Fill( LAMBDA_OBS[p], FLUX[p], wi );
				}
				else if (mean_over3Pixels[4]==3.) {
					mean_over3Pixels[0] /= mean_over3Pixels[3];
					mean_over3Pixels[1] /= mean_over3Pixels[3];
					mean_over3Pixels[2] /= mean_over3Pixels[3];
					zi = mean_over3Pixels[1]/lambdaRFLine__ -1.;
					wi = pow( (zi+1.)/onePlusZ0__, halfGama__);
					v_tmp_r.push_back(hConvertRedshDist->Interpolate(zi));
					v_tmp_d.push_back(mean_over3Pixels[2]);
					v_tmp_w.push_back(wi);
					v_tmp_z.push_back(zi);
					v_tmp_lRF.push_back(mean_over3Pixels[0]);
					v_tmp_lObs.push_back(mean_over3Pixels[1]);
					mean_over3Pixels[0] = 0.;
					mean_over3Pixels[1] = 0.;
					mean_over3Pixels[2] = 0.;
					mean_over3Pixels[3] = 0.;
					mean_over3Pixels[4] = 0.;

					meanDelta[0] += wi*mean_over3Pixels[2];
					meanDelta[1] += wi;
					meanDelta[2] ++;

					tp_flux_vs_lambdaRF->Fill( mean_over3Pixels[0], mean_over3Pixels[2], wi );
					tp_flux_vs_lambdaOBS->Fill( mean_over3Pixels[1], mean_over3Pixels[2], wi );
				}
			}

			if ( v_tmp_r.size()==0. ) continue;

			nbGoodForest ++;
			v_ra__.push_back(X);
			v_de__.push_back(Y);
			v_zz__.push_back(Z);
			v_nbPixelDelta1__.push_back(v_tmp_r.size());
			v_idx__.push_back(forestIdx-1);
			v_r__.push_back(v_tmp_r);
			v_d__.push_back(v_tmp_d);
			v_w__.push_back(v_tmp_w);
			v_z__.push_back(v_tmp_z);
			v_lRF__.push_back(v_tmp_lRF);
			v_lObs__.push_back(v_tmp_lObs);

			std::vector< unsigned int > v_tmp_nb( v_tmp_r.size() ,0);
			v_nb__.push_back(v_tmp_nb);
			v_residual_delta_vs_lRF__.push_back(v_tmp_lRF);
			v_residual_delta_vs_lObs__.push_back(v_tmp_lObs);

			if (nbForest_!=0 && nbGoodForest==nbForest_) {
				breakk = true;
				fits_close_file(fitsptrSpec,&sta);
				break;
			}
		}
		fits_close_file(fitsptrSpec,&sta);
	}

        nbForest_ = v_ra__.size();
        std::cout << "  number of forest        = " << forestIdx << std::endl;
        std::cout << "  number of good forest   = " << nbForest_ << std::endl;

	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;

	//// Replace the values of delta by the one of the stack of delta vs. lambda_RF 
	for ( unsigned int i=0; i<nbForest_; i++ ) {
		const unsigned int nbDelta = v_d__[i].size();
		for ( unsigned int j=0; j<nbDelta; j++ ) {
			v_residual_delta_vs_lRF__[i][j]  = 0.;
			v_residual_delta_vs_lObs__[i][j] = 0.;
		}
	}

	delete tp_flux_vs_lambdaRF;
	delete tp_flux_vs_lambdaOBS;
	delete hConvertRedshDist;

	return;
}




// ---------------------------------------------------------------------
//
//		Change data
//
// ---------------------------------------------------------------------
void Correlation::removeFalseCorrelations(bool firstPass/*=true*/) {

	std::cout << "\n\n\n  ------ removeFalseCorrelations ------ " << std::endl;

	//// Constants
	const unsigned int nbLoop = 10;

	TH1D* hDeltaVsLambdaObs[nbLoop+1];
	TH1D* hDeltaVsLambdaObs_residual[nbLoop+1];
	for (unsigned int i=0; i<nbLoop+1; i++) {

		// Mean transmission flux
		TString name = "hDeltaVsLambdaObs_";
		name += i;
		hDeltaVsLambdaObs[i] = new TH1D(name,"",nbBinlambdaObs__+200,lambdaObsMin__-100.,lambdaObsMax__+100.);
		R_dealWithPlots_1D(hDeltaVsLambdaObs[i], "#lambda_{Obs.} (A)", "Mean transmission flux", "Method2: mean transmission flux");
		for (unsigned int j=0; j<nbBinlambdaObs__+200; j++) {
			hDeltaVsLambdaObs[i]->SetBinContent(j+1,1.);
			hDeltaVsLambdaObs[i]->SetBinError(j+1,0.);
		}
		// Mean transmission flux
		name = "hDeltaVsLambdaObs_residual_";
		name += i;
		hDeltaVsLambdaObs_residual[i] = new TH1D(name,"",nbBinlambdaObs__+200,lambdaObsMin__-100.,lambdaObsMax__+100.);
		R_dealWithPlots_1D(hDeltaVsLambdaObs_residual[i], "#lambda_{Obs.} (A)", "Mean transmission flux", "residual");
	}

	//// Start the loop
	for (unsigned int lp=0; lp<nbLoop; lp++) {

		double data[nbBinlambdaObs__][4];
		for (unsigned int i=0; i<nbBinlambdaObs__; i++) {
			for (unsigned int j=0; j<4; j++) {
				data[i][j] = 0.;
			}
		}

		//// Find < delta >
		for (unsigned int f=0; f<nbForest_; f++) {

			const unsigned int nb = v_nbPixelDelta1__[f];	
			for (unsigned int i=0; i<nb; i++) {

				const double l  = v_lObs__[f][i];
				double mean_flux = 1.;
				if (lp>0) mean_flux = hDeltaVsLambdaObs[lp-1]->Interpolate(l);
				const double d = v_d__[f][i]/mean_flux -1.;
				const double w = v_w__[f][i]*mean_flux*mean_flux;
				const double wd = w*d;

				const unsigned int idx = int( l-lambdaObsMin__);
				data[idx][0] += wd;
				data[idx][1] += w;
				data[idx][2] += wd*d;
				data[idx][3] ++;
			}
		}


		double meanDelta[4] = {0.};
		//// Get the number of pixels and the mean delta
		for (unsigned int i=0; i<nbBinlambdaObs__; i++) {
			meanDelta[0] += data[i][0];
			meanDelta[1] += data[i][1];
			meanDelta[2] += data[i][3];
		}
		std::cout << "  step  " << lp << std::endl;
		std::cout << "  < delta >        = " << meanDelta[0]/meanDelta[1] << std::endl;
		std::cout << "  sum(w_i)         = " << meanDelta[1]              << std::endl;
		std::cout << "  nb pixel         = " << (long long unsigned int)meanDelta[2]              << std::endl;

	
		unsigned int loopIdxForHist = lp-1;
		if (lp==0) loopIdxForHist = nbLoop;
	
		//// Delta vs. lambda_Obs
		for (unsigned int i=0; i<nbBinlambdaObs__; i++) {
			if (data[i][3]>1.) {
				double mean = data[i][0]/data[i][1];
				double err  = sqrt( (data[i][2]/data[i][1]-mean*mean)/data[i][3] );
				hDeltaVsLambdaObs[lp]->SetBinContent(i+100+1,(mean+1.)*hDeltaVsLambdaObs[loopIdxForHist]->GetBinContent(i+100+1) );
				hDeltaVsLambdaObs[lp]->SetBinError(i+100+1,err);
				hDeltaVsLambdaObs_residual[lp]->SetBinContent(i+100+1,mean);
				hDeltaVsLambdaObs_residual[lp]->SetBinError(i+100+1,err);
			}
		}
	
		// Method2: Set "hMeanFlux"
		double xxx1_value_m2 = 0.;
		double xxx2_value_m2 = 0.;
		double yyy1_value_m2 = 0.;
		double yyy2_value_m2 = 0.;
		double value_m2 = 0.;
		unsigned int nbEmptyPixels_m2 = 0;
		double valueFirstNotEmptyPixels = 0.;
		unsigned int idxFirstNotEmptyPixels = 0;
		for (unsigned int i=0; i<nbBinlambdaObs__+200; i++) {
	
			if (valueFirstNotEmptyPixels==0. && hDeltaVsLambdaObs[lp]->GetBinError(i+1)!=0.) {
				valueFirstNotEmptyPixels = hDeltaVsLambdaObs[lp]->GetBinContent(i+1);
				idxFirstNotEmptyPixels   = i;
			}
			
			// If empty and not the last pixel
			if (hDeltaVsLambdaObs[lp]->GetBinError(i+1) == 0. && i!=nbBinlambdaObs__+200-1) {
				if (nbEmptyPixels_m2 == 0 && i!=0) {
					xxx1_value_m2 = hDeltaVsLambdaObs[lp]->GetBinCenter(i);
					yyy1_value_m2 = hDeltaVsLambdaObs[lp]->GetBinContent(i);
					value_m2 = hDeltaVsLambdaObs[lp]->GetBinContent(i);
				}
				nbEmptyPixels_m2 ++;
			}
			else {
				
				// If first pixel not empty after some empty pixels
				if (nbEmptyPixels_m2!=0) {
	
					double a = 0.;
					double b = 0.;					
					
					// Find the mean value between the two not empty pixels
					if (value_m2!=0. && hDeltaVsLambdaObs[lp]->GetBinError(i+1)!=0.) {
						// Linear extrapolation
						xxx2_value_m2 = hDeltaVsLambdaObs[lp]->GetBinCenter(i+1);
						yyy2_value_m2 = hDeltaVsLambdaObs[lp]->GetBinContent(i+1);
						a = yyy1_value_m2 - yyy2_value_m2;
						b = yyy2_value_m2*xxx1_value_m2 - yyy1_value_m2*xxx2_value_m2;
						if (xxx1_value_m2-xxx2_value_m2 != 0.) {
							a /= (xxx1_value_m2-xxx2_value_m2);
							b /= (xxx1_value_m2-xxx2_value_m2);
						}
						value_m2 = (value_m2 + hDeltaVsLambdaObs[lp]->GetBinContent(i+1))/2.;
					}
					
					// Set all the empty pixels 
					for (unsigned int j=0; j<nbEmptyPixels_m2; j++) {
						double tmp_value = a*hDeltaVsLambdaObs[lp]->GetBinCenter(i-j) + b;
						if (hDeltaVsLambdaObs[lp]->GetBinError(i+1)==0.) tmp_value = value_m2;
						hDeltaVsLambdaObs[lp]->SetBinContent(i-j, tmp_value); //value_m2
						//hDeltaVsLambdaObs[lp]->SetBinError(i-j, maxError_m2);
					}
					
					// If the last pixel is also empty
					if (i==nbBinlambdaObs__+200-1) {
						hDeltaVsLambdaObs[lp]->SetBinContent(i+1, value_m2);
						//hDeltaVsLambdaObs[lp]->SetBinError(i+1, maxError_m2);
					}
	
					// Set all the values to zero
					nbEmptyPixels_m2 = 0;
					value_m2 = 0.;
					xxx1_value_m2 = 0.;
					xxx2_value_m2 = 0.;
					yyy1_value_m2 = 0.;
					yyy2_value_m2 = 0.;
				}
			}
		}
	
		//// Set values for the last pixels
		for (unsigned int i=0; i<idxFirstNotEmptyPixels; i++) {
			hDeltaVsLambdaObs[lp]->SetBinContent(i+1, valueFirstNotEmptyPixels);
		}
	}

	double meanDelta[4] = {0.};
	//// Set the new delta
	for (unsigned int f=0; f<nbForest_; f++) {
		const unsigned int nb = v_nbPixelDelta1__[f];	
		for (unsigned int i=0; i<nb; i++) {
			const double l = v_lObs__[f][i];
			const double mean_flux = hDeltaVsLambdaObs[nbLoop-1]->Interpolate(l);
			v_d__[f][i] = v_d__[f][i]/mean_flux -1.;
			v_w__[f][i] = v_w__[f][i]*mean_flux*mean_flux;
			v_residual_delta_vs_lObs__[f][i] = hDeltaVsLambdaObs_residual[nbLoop-1]->Interpolate(l);

			meanDelta[0] += v_w__[f][i]*v_d__[f][i];
			meanDelta[1] += v_w__[f][i];
			meanDelta[2] ++;
		}
	}
	std::cout << "\n  step  LAST" << std::endl;
	std::cout << "  < delta >        = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)         = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel         = " << (long long unsigned int)meanDelta[2]              << std::endl;

	for (unsigned int i=0; i<nbLoop+1; i++) {
		delete hDeltaVsLambdaObs[i];
		delete hDeltaVsLambdaObs_residual[i];
	}
}


















