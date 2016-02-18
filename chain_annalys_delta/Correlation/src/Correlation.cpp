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
#include "fitsio.h"
#include "TString.h"

#include "Correlation.h"
#include "LymanAlphaForestRegion.h"
#include "Cosmology.h"
#include "../../../Root/Library/RootHistoFunctions.h"
#include "../../../Cpp/Library/mathFunctions.h"
#include "../../../Constants/constants.h"
#include "../../../Constants/globalValues.h"

/*
std::string alQSOName[9] = {"QSO_ALL_LowZ", "QSO_ALL_highZ", "QSO_ALL_SGC", "QSO_ALL_NGC", "QSO_ALL_DR7", "QSO_ALL_DR12", "QSO_ALL_CORE", "QSO_ALL_LowMg", "QSO_ALL_HighMg"};
std::string alQSO[9]     = {"QSO_ALL_LowZ_2_397.fits", "QSO_ALL_highZ_2_397.fits", "QSO_ALL_SGC.fits", "QSO_ALL_NGC.fits", "QSO_ALL_DR7.fits", "QSO_ALL_DR12.fits", "QSO_ALL_CORE.fits", "QSO_ALL_LowMg.fits", "QSO_ALL_HighMg.fits"};
*/
std::string QSO__    = "";
std::string pathQ1__ = "";
std::string alQSOName[19] = {"QSO",
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
				"ALL_EVERY_OBJECTS"
};
std::string alQSO[19] = {"QSO_ALL_TESTS",
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
				"QSO_EBOSS",
				"QSO_DR7_DR12_EBOSS_2016_01_08",
				"ALL_OBJECTS",
				"ALL_EVERY_OBJECTS_2016_01_08"
};

///// Constants
const unsigned int nbBinlambdaObs__  = int(lambdaObsMax__-lambdaObsMin__);
double distMinPixel__ = 0.;
double distMinPixelDelta2__ = 0.;
unsigned int idxCommand_[6] = {0};
const std::string pathToMockJMC__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1563/";
std::string pathToSave__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_nicolasEstimator/";  //_nicolasEstimator
//std::string pathToSave__ = "";

///// Flags for Jean-Marc's simulations
const bool mocks              = false;
const bool mockJMC__          = true;
const bool mockBox__          = false;
const bool mocksNoNoiseNoCont = false;
const double randomPositionOfQSOInCellNotBeforeCorrelation__ = true;
//// Flags for covariance matrix estimation
const bool shuffleQSO     = false;
const bool shuffleForest  = false;
const bool randomQSO      = true;
const bool randomForest   = false;
const bool doBootstraps__ = false;


const bool doVetoLines__ = true;
const bool nicolasEstimator__ = true;
const bool doingCoAddSoRemoving__ = false;


Correlation::Correlation(int argc, char **argv) {

	std::cout << std::scientific;
	std::cout.precision(14);

	// command

	///// 0: correlation, 1: mockCatalogue, 2: bootstraps, 3: Wick, 4: mockChunck, 5: mockSimul
	for (unsigned int i=0; i<6; i++) {
		std::string string = argv[i+1];
		idxCommand_[i] = atoi(string.c_str());
	}

	
	QSO__    = alQSOName[idxCommand_[1]];
	pathQ1__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/";
	pathQ1__ += alQSO[idxCommand_[1]];
	pathQ1__ += ".fits";
	

	///// Set the number of forest to work on
	nbForest_   = 0;
	nbForest2__ = 0;
	nbQ1__      = 0;
	nbQ2__      = 0;
	if (!mocks && !mockJMC__) {
		pathForest__  = "/home/gpfs/manip/mnt/bao/hdumasde/Data/";
		pathForest__  += forest__;
		pathForest__  += "/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits";  //_method1
//		pathForest__  += "/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits";
//		pathForest__  += "/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits";
//		pathForest__   = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/FitsFile_DR12_reOBS_eBOSS_Guy/DR12_primery/DR12_primery_reOBS_eBOSS.fits";
//		pathForest__   = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/FitsFile_DR12_reOBS_Guy/DR12_primery/DR12_primery_reOBS.fits";
//		pathForest__   = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy/FitsFile_DR12_reOBS_eBOSS_noCoADD_Guy/DR12_primery/DR12_primery_reOBS_eBOSS_noCoADD.fits";
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
		pathQ1__      = pathForest__;
		pathToSave__  = pathForest__;

		if (!mocksNoNoiseNoCont) pathForest__ += "/Data/delta.fits";
		else pathForest__ += "/Data/delta.fits";
		
		pathQ1__ += "/Data/QSO_withRSD.fits";
		
		if (nicolasEstimator__) pathToSave__ += "/Results_NicolasDistortion/";
		else pathToSave__ += "/Results/";
		
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

	//
	if      (command == 0)  xi_1D_delta_delta();
	else if (command == 1)  xi_1DlRF_delta_delta();
	else if (command == 2)  xi_1DlRFDevide_delta_delta();
	else if (command == 3)  xi_1D_delta_delta2();
	else if (command == 4)  xi_1DlRFDevide_delta_delta2();
	//
	else if (command == 5)  xi_A_delta_delta();
	else if (command == 6)  xi_A_delta_delta_lambda();
	else if (command == 7)  xi_A_delta_delta2( idxCommand_[2] );
	else if (command == 8)  xi_A_delta_delta2_lambda();
	//
	else if (command == 9)  xi_delta_QSO(doBootstraps__,idxCommand_[2]);
	else if (command == 10) xi_delta_QSO_theta(doBootstraps__,idxCommand_[2]);
	else if (command == 11) xi_delta_QSO_lambda(doBootstraps__,idxCommand_[2]);
	else if (command == 12) xi_delta_QSO_distortionMatrix();
	else if (command == 13) xi_delta_QSO_distortionMatrix_1D();
	else if (command == 14) xi_delta_QSO_Metals_Models(absorber__[idxCommand_[2]],absorberName__[idxCommand_[2]]);
	else if (command == 15) xi_delta_QSO_Wick(idxCommand_[3]);
	//
	else if (command == 16) xi_QSO_QSO(doBootstraps__,idxCommand_[2]);
	else if (command == 17) xi_Q1_Q2();
	//
	else if (command == 18) xi_A_delta_delta_MockJMc();
	else if (command == 19) xi_delta_QSO_MockJMc(doBootstraps__,idxCommand_[2]);
	else if (command == 20) xi_delta_QSO_MockJMc_distortionMatrix();
	else if (command == 21) xi_delta_QSO_MockJMc_distortionMatrix_1D();
	else if (command == 22) xi_delta_QSO_Metals_Models_MockJMc(absorber__[idxCommand_[2]],absorberName__[idxCommand_[2]]);
	else if (command == 23) xi_QSO_QSO_MockJMc(doBootstraps__,idxCommand_[2]);

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

	std::cout << "\n\n\n\n  ------ xi_1D_delta_delta ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_1D_delta_delta.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	loadDataForest(pathForest__);

	if (mocksNoNoiseNoCont) {
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
	fFile.precision(14);

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

	std::cout << "\n\n\n\n  ------ xi_1DlRF_delta_delta ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_1DlRF_delta_delta.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	loadDataForest(pathForest__);


	if (mocksNoNoiseNoCont) {
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
	fFile.precision(14);

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

	std::cout << "\n\n\n\n  ------ xi_1DlRFDevide_delta_delta ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_1DlRFDevide_delta_delta.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	loadDataForest(pathForest__);


	if (mocksNoNoiseNoCont) {
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
	fFile.precision(14);

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
void Correlation::xi_1D_delta_delta2(void) {

	std::cout << "\n\n\n\n  ------ xi_1D_delta_delta2 ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_1D_delta_delta2.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	///// Forest 1
	loadDataForest(pathForest__);
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	///// Forest 2
	loadDataDelta2();
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

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {
		
		const double ra1       = v_ra__[f1];
		const double de1       = v_de__[f1];
		const double zz1       = v_zz__[f1];
		const unsigned int nb1 = v_nbPixelDelta1__[f1];
		
		for (unsigned int f2=0; f2<nbForest2__; f2++) {

			if (ra1!=v_raDelta2__[f2] || de1!=v_deDelta2__[f2] || zz1!=v_zzDelta2__[f2]) continue;
			
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
	fFile.precision(14);

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

	std::cout << "\n\n\n\n  ------ xi_1DlRF_delta_delta2 ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_1DlRF_delta_delta2.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	///// Forest 1
	loadDataForest(pathForest__);
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_idx__.clear();
	v_r__.clear();
	v_lObs__.clear();
	v_nb__.clear();
	///// Forest 2
	loadDataDelta2();
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

	std::cout << "\n  Starting\n" << std::endl;

	for (unsigned int f1=0; f1<nbForest_; f1++) {
		
		const double ra1       = v_ra__[f1];
		const double de1       = v_de__[f1];
		const double zz1       = v_zz__[f1];
		const unsigned int nb1 = v_nbPixelDelta1__[f1];

		for (unsigned int f2=0; f2<nbForest2__; f2++) {

			if (ra1!=v_raDelta2__[f2] || de1!=v_deDelta2__[f2] || zz1!=v_zzDelta2__[f2]) continue;
			
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
	fFile.precision(14);

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
void Correlation::xi_A_delta_delta(void) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_A_delta_delta.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	///// Load forest
	loadDataForest(pathForest__);
	v_zz__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Constants:
	const unsigned int max    = 200.;
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
	double dataMu[nbBin][nbBinM][7];
	double data2D[nbBin][nbBin][7];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<7; k++) {
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

				const double w1 = v_w__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double z1 = v_z__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-tmp_rParral);
					if (rParral>=max) continue;
					const double distTotPow2 = rPerp*rPerp + rParral*rParral;

					const double w1w2 = w1*v_w__[f2][i2];
					const double d1d2 = d1*v_d__[f2][i2];
					const double z1z2 = z1+v_z__[f2][i2];

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int(mu*50.);

						dataMu[idx][idxM][0] += w1w2*d1d2;
						dataMu[idx][idxM][1] += w1w2*d1d2*d1d2;
						dataMu[idx][idxM][2] += w1w2*distTot;
						dataMu[idx][idxM][3] += w1w2*mu;
						dataMu[idx][idxM][4] += w1w2*z1z2;
						dataMu[idx][idxM][5] += w1w2;
						dataMu[idx][idxM][6] ++;
					}
					
					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += w1w2*d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][1] += w1w2*d1d2*d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2*z1z2;
					data2D[rPerpBinIdx][rParralBinIdx][5] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][6] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::string prefix = forest__;
	prefix += ".txt";

	std::ofstream fFile;
	long double sumOne = 0.;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;

	///// Save the 2D cross-correlation
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_2D_";
	pathToSave += prefix;
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);
	
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
	fFile.precision(14);

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
void Correlation::xi_A_delta_delta_lambda(void) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_lambda ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_A_delta_delta_lambda.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	///// Load forest
	loadDataForest(pathForest__);
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
	fFile.precision(14);

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
void Correlation::xi_A_delta_delta2( unsigned int bootIdx/*=0*/ ) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta2 ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_A_delta_delta2.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	///// Forest 1
	loadDataForest(pathForest__);
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
	const unsigned int max = 200.;
	const unsigned int nbBin = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100.;

	const double maxPow2      = max*max;
	const double distMax      = sqrt(2.)*max+1.;
	const double distMaxPow2  = 2.*max*max;

	///// Find the maximal angle of seperation
	const double maxTheta = 2.*asin( distMax/(std::min(distMinPixel__,distMinPixelDelta2__)*2.) );
	const double cosmaxTheta = cos(maxTheta);
	std::cout << "  maxTheta = " << maxTheta*180./M_PI << " degree" << std::endl;

	///// Arrays for data
	double dataMu[nbBin][nbBinM][7];
	double data2D[nbBinX][nbBinY][7];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<7; k++) {
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

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = v_rDelta2__[f2][i2]-tmp_rParral;
					if (fabs(rParral)>=max) continue;
					const double distTotPow2 = rPerp*rPerp + rParral*rParral;

					const double w1w2     = w1*v_wDelta2__[f2][i2];
					const double d1d2     = d1*v_dDelta2__[f2][i2];
					const double w1w2z1z2 = w1w2*(z1+v_zDelta2__[f2][i2]);

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int((1.+mu)*50.);

						dataMu[idx][idxM][0] += w1w2*d1d2;
						dataMu[idx][idxM][1] += w1w2*d1d2*d1d2;
						dataMu[idx][idxM][2] += w1w2*distTot;
						dataMu[idx][idxM][3] += w1w2*mu;
						dataMu[idx][idxM][4] += w1w2z1z2;
						dataMu[idx][idxM][5] += w1w2;
						dataMu[idx][idxM][6] ++;
					}
					
					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(max+rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += w1w2*d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][1] += w1w2*d1d2*d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2z1z2;
					data2D[rPerpBinIdx][rParralBinIdx][5] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][6] ++;

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
	fFile.precision(14);

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
	pathToSave += prefix;
	if (doBootstraps__ || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

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
		

		/// find the index of each forest among the C_NBSUBSAMPLES sub-samples
		LymanForest* lymanForestObject = new LymanForest(forests, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		pathToSave = pathToSave__;
		pathToSave += "xi_A_delta_delta2_map_";
		pathToSave += prefix;
		lymanForestObject->SaveRegionMap(pathToSave);
		delete lymanForestObject;
	}


	return;
}
void Correlation::xi_A_delta_delta2_lambda(void) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta2_lambda ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_A_delta_delta2_lambda.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	///// Load forest 1
	loadDataForest(pathForest__);
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
	fFile.precision(14);

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
void Correlation::xi_delta_QSO(bool doBootstraps/*=False*/, unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	if (doBootstraps)  std::cout << "  subsampling N° "     << bootIdx << std::endl;
	if ( shuffleForest || shuffleQSO || randomQSO ) std::cout << "  shuffleForest seed " << bootIdx*10 << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__,doBootstraps,bootIdx);
	if (!nicolasEstimator__ && (doBootstraps || doingCoAddSoRemoving__) ) removeFalseCorrelations();

	///// Empty useless vectors
	v_zz__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

	///// Constants:
	///// The space between bins is of 10 Mpc.h^-1
	const double max          = 200.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100.;

	const double maxPow2      = max*max;

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
		std::srand (bootIdx*10);
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
		std::vector<double> tmp_r(v_rQ1__);
		///// Put the new data
		for (unsigned int i=0; i<nbQ1__; i++) {
			const unsigned int ii = randomIdx[i];
			v_zzQ1__[i] = tmp_zz[ii];
			v_rQ1__[i]  = tmp_r[ii];
		}
	}
	///// Needed to randomize the QSO
	if (randomQSO) {

		std::cout << "\n\n  Random QSO " << bootIdx << std::endl;

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

		std::srand (bootIdx*10);
		const double coefRA = 1./(raMax-raMin);
		const double coefDE = 1./(deMax-deMin);
		for (unsigned int i=0; i<nbQ1__; i++) {
			const double ra = (double)rand()/(double)(RAND_MAX*coefRA) +raMin;
			const double de = (double)rand()/(double)(RAND_MAX*coefDE) +deMin;
			v_raQ1__[i]    = ra;
			v_deQ1__[i]    = de;
			v_CosDeQ1__[i] = cos(de);
			v_SinDeQ1__[i] = sin(de);
		}
	}

	///// get an array for nb of pairs for the forest
	double a_nbPairs[nbForest_];
	for (unsigned int i=0; i<nbForest_; i++) {
		a_nbPairs[i] = 0.;
	}

	///// Arrays for data
	double data2D[nbBinX][nbBinY][7];
	double dataMu[nbBin][nbBinM][7];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<7; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}

	std::cout << "\n  Starting\n" << std::endl;

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
			const double distTransQsoLyaPow2 = v_rQ1__[q]*v_rQ1__[q]*(1.-cosTheta*cosTheta);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			///// Parrallel distance between the qso and the lya
			const double distParalQsoLya = v_rQ1__[q]*cosTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya );
			const double zQSO = v_zzQ1__[q];

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
				const double w   = v_w__[f][ii];
				const double d   = v_d__[f][ii];
				const double wd  = w*d;
				const double wdd = wd*d;
	
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
				}

				///// Fill the histogramm of xi(r_{perp}, r_{paral}
				const unsigned int rParralBinIdx = int(distP+max);
				data2D[rPerpBinIdx][rParralBinIdx][0] += wd;
				data2D[rPerpBinIdx][rParralBinIdx][1] += wdd;
				data2D[rPerpBinIdx][rParralBinIdx][2] += w*distTransQsoLya;
				data2D[rPerpBinIdx][rParralBinIdx][3] += w*distP;
				data2D[rPerpBinIdx][rParralBinIdx][4] += w*(zQSO+v_z__[f][ii]);
				data2D[rPerpBinIdx][rParralBinIdx][5] += w;
				data2D[rPerpBinIdx][rParralBinIdx][6] ++;

				///// Get the number of pairs
				nbPairs += w;
			}
		}
	
		///// Put the number of pairs comming with this forest
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
	if (doBootstraps)  prefix += "subsampling";
	if (shuffleForest) prefix += "shuffleForest";
	if (shuffleQSO)    prefix += "shuffleQSO";
	if (randomQSO)     prefix += "randomQSO";
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
	if (doBootstraps || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

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
	if (doBootstraps || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);
	
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

	
	if (!doBootstraps && !shuffleForest && !shuffleQSO && !randomQSO) {

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
void Correlation::xi_delta_QSO_theta(bool doBootstraps/*=False*/, unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_theta ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_theta.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	///// Constants
	const double radToDeg = 180.*M_PI;

	///// QSO
	loadDataQ1();
	v_zzQ1__.clear();
	v_rQ1__.clear();
	v_nb__.clear();
	///// Forest
	loadDataForest(pathForest__,doBootstraps,bootIdx);
	v_zz__.clear();
	v_r__.clear();
	v_z__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();

	///// Constants:
	///// The space between bins is of 10 Mpc.h^-1
	const double max         = 0.1;
	const unsigned int nbBin = 100;
	const double inverseBinSize = nbBin/max;

	///// Arrays for data
	double data[nbBin][7];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int k=0; k<5; k++) {
			data[i][k] = 0.;
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
		
		for (unsigned int q=0; q<nbQ1__; q++) {

			///// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT) continue;

			///// Angle between the two directions of the qso and the lya
			const double theta = radToDeg*acos( cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q] );
			if (theta>max) continue;
			const unsigned int idx  = int(theta*inverseBinSize);

			///// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {
				
				const double w   = v_w__[f][i];
				const double d   = v_d__[f][i];
				const double wd  = w*d;
				const double wdd = wd*d;

				data[idx][0] += wd;
				data[idx][1] += wdd;
				data[idx][2] += w*theta;
				data[idx][3] += w;
				data[idx][4] ++;
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_theta_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

	long double sumOne = 0.;
	long double sumWeg = 0.;

	///// Set the values of data
	///// [0] for value, [1] for error, [2] for bin center
	for (unsigned int i=0; i<nbBin; i++) {

		fFile << data[i][0];
		fFile << " " << data[i][1];
		fFile << " " << data[i][2];
		fFile << " " << data[i][3];
		fFile << " " << data[i][4];
		fFile << std::endl;

		sumWeg += data[i][3];
		sumOne += data[i][4];
	}
	fFile.close();

	std::cout << "  number of pairs      = " << sumOne          << std::endl;
	std::cout << "  number of w          = " << sumWeg          << std::endl;

	return;
}
void Correlation::xi_delta_QSO_lambda(bool doBootstraps/*=False*/, unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_lambda ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_lambda.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	if (doBootstraps)  std::cout << "  subsampling N° "     << bootIdx << std::endl;
	if (shuffleForest) std::cout << "  shuffleForest seed " << bootIdx*10 << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__,doBootstraps,bootIdx);
	if (doBootstraps) removeFalseCorrelations();
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
	fFile.precision(14);

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
void Correlation::xi_delta_QSO_distortionMatrix(void) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_distortionMatrix ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_distortionMatrix.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);

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
			const double distTransQsoLyaPow2 = v_rQ1__[q]*v_rQ1__[q]*(1.-cosTheta*cosTheta);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			///// Parrallel distance between the qso and the lya
			const double distParalQsoLya = v_rQ1__[q]*cosTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya*fromValToIdx );

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
	
				const double distP = v_r__[f][ii] - distParalQsoLya;				

				///// Look at the position of the Lya forest regarding the qso
				if (fabs(distP) >= max) {
					if (infPosBool || supPosBool) break;
					else continue;
				}

				const double w    = v_w__[f][ii];
				const double val0 = w*invSumWeight;
	
				///// Fill the histogramm of xi(r_{perp}, r_{paral}
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
	fFile.precision(14);
	
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

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_distortionMatrix_1D ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_distortionMatrix.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	///// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	///// Forest
	loadDataForest(pathForest__);

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
			const double distTransQsoLyaPow2 = v_rQ1__[q]*v_rQ1__[q]*(1.-cosTheta*cosTheta);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			///// Parrallel distance between the qso and the lya
			const double distParalQsoLya = v_rQ1__[q]*cosTheta;

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
	
				const double distP = v_r__[f][ii] - distParalQsoLya;				

				///// Look at the position of the Lya forest regarding the qso
				if (fabs(distP) >= max) {
					if (infPosBool || supPosBool) break;
					else continue;
				}

				const double distTotPow2 = distTransQsoLyaPow2 + distP*distP;
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
	fFile.precision(14);
	
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

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	//// Get the distance if it was this metal
	std::vector< std::vector< double > > v_r_metal(v_r__);
	for (unsigned int i=0; i<v_r_metal.size(); i++) {
		for (unsigned int j=0; j<v_r_metal[i].size(); j++) {
			v_r_metal[i][j] = hConvertRedshDist->Interpolate( v_lObs__[i][j]/lambdaFrMetal-1. );
		}
	}

	//// Empty useless vectors
	delete cosmo;
	v_zz__.clear();
	v_d__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();


	//// Get the monopole, quadrupol, exadecapol
	std::string pathToLoad = PATHTOWORK;
	pathToLoad += PATHTOCAMB;
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
	const unsigned int nbBinM = 100.;

	const double maxPow2      = max*max;


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
			//if (fabs(ra-v_raQ1__[q])<1.e-9 && fabs(dec-v_deQ1__[q])<1.e-9 ) continue;
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT ) continue;  // && fabs(z-v_zzQ1__[q])<0.001

			///// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];

			///// reject QSO with a distance too large
			const double distTransQsoLyaPow2 = v_rQ1__[q]*v_rQ1__[q]*(1.-cosTheta*cosTheta);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			///// Parrallel distance between the qso and the lya
			const double distParalQsoLya = v_rQ1__[q]*cosTheta;

			///// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel - distParalQsoLya;
			if ( distParalQsoFirstPixel >= max) continue;

			///// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			///// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya );
			const double zQSO = v_zzQ1__[q];

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
				const double distTot    = sqrt(distTotPow2);
				const double distTotMetal = sqrt(distTransQsoLyaPow2 + (v_r_metal[f][ii]-distParalQsoLya)*(v_r_metal[f][ii]-distParalQsoLya) );

				const double w   = v_w__[f][ii];
				unsigned int idxBinCAMB = nbBinCAMB;
				if (distTotMetal<maxDistCAMB) idxBinCAMB = int(distTotMetal);
				const double wxi0 = w*data_xi0[idxBinCAMB];
				const double wxi2 = w*data_xi2[idxBinCAMB];
				const double wxi4 = w*data_xi4[idxBinCAMB];
				const double wz   = w*(zQSO+v_z__[f][ii]);

				if (distTot < max) {
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
				const unsigned int rParralBinIdx = int(distP+max);
				data2D[rPerpBinIdx][rParralBinIdx][0] += wxi0;
				data2D[rPerpBinIdx][rParralBinIdx][1] += wxi2;
				data2D[rPerpBinIdx][rParralBinIdx][2] += wxi4;
				data2D[rPerpBinIdx][rParralBinIdx][3] += w*distTransQsoLya;
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
	fFile.precision(14);
	
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
	fFile.precision(14);
	
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
void Correlation::xi_delta_QSO_Wick(unsigned int diagramIdx) {

	//// Convert the diagram index to string
	std::stringstream convert;
	convert << diagramIdx;
	std::string strConvert = convert.str();

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_Wick ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO.py";
	command += commandEnd__;
	command += " ";
	command += strConvert;
	std::cout << command << "\n" << std::endl;
	
	//// QSO
	loadDataQ1();
	//// Forest
	loadDataForest(pathForest__);
	v_zz__.clear();
	v_idx__.clear();
	v_lRF__.clear();
	v_lObs__.clear();





	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max          = 200.;

	const unsigned int nbBin   = 50;
	const double inverseBinSize = nbBin/max;
	//const unsigned int nbBin2D = nbBin*nbBin*2;
	const double maxPow2       = max*max;

	//// Array with the covaraince
	double covVar1D[nbBin][3];
	double cov1D[nbBin][nbBin][3];
	for (unsigned int i=0; i<nbBin; i++) {
		covVar1D[i][0] = 0.;
		covVar1D[i][1] = 0.;
		covVar1D[i][2] = 0.;
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<3; k++) {
				cov1D[i][j][k] = 0.;
			}
		}
	}
	/*
	double data2D[nbBin2D][nbBin2D][2];
	for (unsigned int i=0; i<nbBin2D; i++) {
		for (unsigned int j=0; j<nbBin2D; j++) {
			for (unsigned int k=0; k<2; k++) {
				data1D[i][j][k] = 0.;
			}
		}
	}*/

	//// ------------------------------------------------------------------
	//// Load correlation delta-delta same forest
	std::string pathToFile = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/xi_1D_delta_delta_LYA.txt";
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


	//// ------------------------------------------------------------------
	//// Load correlation delta-delta different forest
	pathToFile = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/Tests/xi_A_delta_delta_2D_LYA.txt";
	std::ifstream fileData2(pathToFile.c_str());

	double xidd3D[nbBin][nbBin];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			xidd3D[i][j] = 0.;
		}
	}

	save0 = 0.;
	save1 = 0.;
	save2 = 0.;
	save3 = 0.;
	save4 = 0.;
	save5 = 0.;
	long double save6 = 0.;
	unsigned int i = 0;

	while (fileData2) {
		fileData2>>save0>>save1>>save2>>save3>>save4>>save5>>save6;
		if (fileData2==0) break;

		const unsigned int iX = i/int(max);
		const unsigned int iY = i%int(max);
		xidd3D[iX][iY] = save0/save5;
		i++;
	}
	fileData2.close();




	//// Get the covariance
	//// ------------------------------------------------------------------
	std::cout << "\n\n  Starting" << std::endl;

	for (unsigned int f=0; f<nbForest_; f++) {

		//// Get number of pixels in forest
		const unsigned int nbPixel = v_nbPixelDelta1__[f];
		const double firstPixel    = v_r__[f][0];
		const double lastPixel     = v_r__[f][nbPixel-1];
		const double cosDe         = v_CosDe__[f];
		const double sinDe         = v_SinDe__[f];
		const double ra            = v_ra__[f];
		const double dec           = v_de__[f];

		
		//// For T1 + T2 + T3 + T4
		std::vector<unsigned int> vect_idx;
		std::vector<double> vect_www;
		std::vector<double> vect_rrr;
		std::vector<unsigned int> vect_idxPix;
		std::vector<unsigned int> vect_idxQSO;
		

		for (unsigned int q=0; q<nbQ1__; q++) {

/*			
			//// For T1 + T2
			std::vector<unsigned int> vect_idx;
			std::vector<double> vect_www;
			std::vector<double> vect_rrr;
*/
			//// Remove the correlation qso-ownForest
			if (fabs(ra-v_raQ1__[q])<C_AUTOCORRCRIT && fabs(dec-v_deQ1__[q])<C_AUTOCORRCRIT) continue;

			//// Angle between the two directions of the qso and the lya
			const double cosTheta = cosDe*v_CosDeQ1__[q]*cos(ra-v_raQ1__[q]) + sinDe*v_SinDeQ1__[q];
			
			//// reject QSO with a distance too large
			const double distTransQsoLyaPow2 = v_rQ1__[q]*v_rQ1__[q]*(1.-cosTheta*cosTheta);
			if (distTransQsoLyaPow2 >= maxPow2) continue;
			
			//// Parrallel distance between the qso and the lya
			const double distParalQsoLya = v_rQ1__[q]*cosTheta;
			
			//// Distance between the qso and the first pixel
			const double distParalQsoFirstPixel = firstPixel - distParalQsoLya;
			if ( firstPixel - distParalQsoLya >= max) continue;
			
			//// Distance between the qso and the last pixel
			const double distParalQsoLastPixel = lastPixel - distParalQsoLya;
			if ( distParalQsoLastPixel <= -max) continue;

			//// Transvers distance between the qso and the lya
			//const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			//const unsigned int rPerpBinIdx = int( distTransQsoLya );
			//const double zQSO = v_zzQ1__[q];
			
			//// 'true' if Ly\alpha forest further than qso
			const bool infPosBool = (distParalQsoFirstPixel > 0.);
			///
			const bool supPosBool = (distParalQsoLastPixel < 0.);
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {
				
				unsigned int ii = i;
				if (supPosBool) ii = nbPixel-1-i;
				
				const double distP = v_r__[f][ii] - distParalQsoLya;				
				
				//// Look at the position of the Lya forest regarding the qso
				if (fabs(distP) >= max) {
					if (infPosBool) break;
					else if (supPosBool) break;
					else continue;
				}

				const double distTotPow2 = distTransQsoLyaPow2 + distP*distP;
				const double w           = v_w__[f][ii];

				if (distTotPow2 < maxPow2) {
					const unsigned int idx1 = sqrt(distTotPow2)*inverseBinSize;
					vect_idx.push_back(idx1);
					vect_www.push_back(w);
					vect_rrr.push_back(v_r__[f][ii]);

					//// For T1 + T2 + T3 + T4
					vect_idxPix.push_back(ii);
					vect_idxQSO.push_back(q);

					//// Variance
					covVar1D[idx1][0] += w*w*xidd1D[0];
					covVar1D[idx1][1] += w*w;
					covVar1D[idx1][2] ++;
				}
			}
	
		}  //// If T1 + T2 + T3 + T4

			//// Calculate the covariance
			const unsigned int nbPairs = vect_idx.size();
			
			for (unsigned int i=0; i<nbPairs; i++) {

				const unsigned int idx1    = vect_idx[i];
				//const unsigned int idxPix1 = vect_idxPix[i];
				//const unsigned int idxQSO1 = vect_idxQSO[i];
				const double       rrr1    = vect_rrr[i];
				const double       www1    = vect_www[i];

				//// Covariance
				for (unsigned int j=0; j<i; j++) {

					//// If T1 + T2
					//if (idxQSO1 != vect_idxQSO[j]) continue;
					//// If T1 + T2 + T3
					//if (idxQSO1 != vect_idxQSO[j] && idxPix1 != vect_idxPix[j]) continue;

					const unsigned int idx2 = vect_idx[j];
					const unsigned int idx  = 1 + int( fabs(rrr1-vect_rrr[j]) );
					const double ww         = www1*vect_www[j];

					cov1D[idx1][idx2][0] += ww*xidd1D[idx];
					cov1D[idx1][idx2][1] += ww;
					cov1D[idx1][idx2][2] ++;
				}
			}
		//}  ///

	}

	std::cout << "\n\n\n" << std::endl;
	std::ofstream fFile;
	std::string pathToSave;

	//// 1D
	pathToSave = pathToSave__;
	pathToSave += "xi_delta_QSO_1D_Wick";
	pathToSave += strConvert;
	pathToSave += "_";
	pathToSave += forest__;
	pathToSave += "_";
	pathToSave += QSO__;
	pathToSave += "_T1234_"; // T12  // T1234  // T124
	std::stringstream convert2;
	convert2 << nbForest_;
	pathToSave += convert2.str();
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

	for (unsigned int i=0; i<nbBin; i++) {
		fFile << i << " " << i;
		fFile << " " << covVar1D[i][0];
		fFile << " " << covVar1D[i][1];
		fFile << " " << covVar1D[i][2];
		fFile << " " << std::endl;
	}

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<i+1; j++) {
			fFile << i << " " << j;
			fFile << " " << cov1D[i][j][0]+cov1D[j][i][0];
			fFile << " " << cov1D[i][j][1]+cov1D[j][i][1];
			fFile << " " << cov1D[i][j][2]+cov1D[j][i][2];
			fFile << " " << std::endl;
		}
	}
	fFile.close();
	
	return;
}



// ---------------------------------------------------------------------
//
//		3D QSO 1 - QSO 2 correlation
//
// ---------------------------------------------------------------------
void Correlation::xi_QSO_QSO(bool doBootstraps/*=false*/, unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_QSO_QSO ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_Q_Q.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
	loadDataQ1();

	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();


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
		const double rr1 = v_rQ1__[q1];

		for (unsigned int q2=0; q2<q1; q2++) {

			//// Angle between the two directions of thes QSOs
			const double cosTheta = cosDe1*v_CosDeQ1__[q2]*cos(ra1-v_raQ1__[q2]) + sinDe1*v_SinDeQ1__[q2];
			
			//// Transvers distance between the two QSOs
			const double distTransQsoQsoPow2 = rr1*rr1*(1.-cosTheta*cosTheta);
			if (distTransQsoQsoPow2 >= maxPow2) continue;


			//// Parrallel distance between the two QSOs
			const double rPara = fabs(rr1*cosTheta-v_rQ1__[q2]);
			if (rPara >= max) continue;
			const unsigned int idxPara = int(rPara);

			//// Transvers distance between the two QSOs
			const double rPerp   = sqrt(distTransQsoQsoPow2);
			const unsigned int idxPerp = int(rPerp);

			const double distTotPow2 = distTransQsoQsoPow2 + rPara*rPara;

			//// 1D
			if (distTotPow2<maxPow2) {
				const double distTot    = sqrt(distTotPow2);
				const double mu         = rPara/distTot;
				const unsigned int idx  = int(distTot);
				const unsigned int idxM = int(50.*mu);

				dataMu[idx][idxM][0] ++;
				dataMu[idx][idxM][1] += distTot;
				dataMu[idx][idxM][2] += mu;
				dataMu[idx][idxM][3] += zz1+v_zzQ1__[q2];
			}

			//// 2D
			data2D[idxPerp][idxPara][0] ++;
			data2D[idxPerp][idxPara][1] += rPerp;
			data2D[idxPerp][idxPara][2] += rPara;
			data2D[idxPerp][idxPara][3] += zz1+v_zzQ1__[q2];
		}
	}
	
	std::cout << "\n  Saving\n" << std::endl;

	// 1D
	std::ofstream fFile;
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_QSO_QSO_1D_";
	pathToSave += QSO__;
	if (randomQSO) {
		pathToSave += "_randomQSO_";
		pathToSave += strBootIdx;;
	}
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

	long double sumOne = 0.;
	long double meanZ = 0.;

	//// Set the values of data
	for (unsigned int i=0; i<nbBin; i++) {

		long double save0 = 0.;
		long double save1 = 0.;
		long double save2 = 0.;
		long double save3 = 0.;

		for (unsigned int j=0; j<nbBinM; j++) {
			save0 += dataMu[i][j][0];
			save1 += dataMu[i][j][1];
			save2 += dataMu[i][j][2];
			save3 += dataMu[i][j][3];
		}

		fFile << save0;
		fFile << " " << save1;
		fFile << " " << save2;
		fFile << " " << save3/2.;
		fFile << std::endl;

		sumOne += save0;
		meanZ  += save3/2.;
	}
	fFile.close();

	//// Get mean redshift of pairs
	std::cout << "  < z >           = " << meanZ/sumOne << " +- " << 0. << std::endl;
	std::cout << "  number of pairs = " << sumOne          << std::endl;


	//// 2D
	pathToSave = pathToSave__;
	pathToSave += "xi_QSO_QSO_2D_";
	pathToSave += QSO__;
	if (randomQSO) {
		pathToSave += "_randomQSO_";
		pathToSave += strBootIdx;;
	}
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

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
		pathToSave += "_randomQSO_";
		pathToSave += strBootIdx;;
	}
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

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

	return;
}
void Correlation::xi_Q1_Q2(void) {

	std::cout << "\n\n\n\n  ------ xi_Q1_Q2 ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_Q1_Q2.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;
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
		const double r1    = v_rQ1__[q1];
		const double z1    = v_zzQ1__[q1];

		for (unsigned int q2=0; q2<nbQ2__; q2++) {

			//// Angle between the two directions of the qso and the lya
			double cosTheta = cosDe1*v_CosDeQ2__[q2]*cos(ra1-v_raQ2__[q2]) + sinDe1*v_SinDeQ2__[q2];
			if (fabs(1.-cosTheta)<1.e-10) cosTheta=1.;
			
			//// Transvers distance between the qso and the lya
			const double distTransQsoQsoPow2 = r1*r1*(1.-cosTheta*cosTheta);
			if (distTransQsoQsoPow2 >= maxPow2) continue;

			//// Parrallel distance between the qso and the lya
			const double distParalQsoQso = r1*cosTheta-v_rQ2__[q2];
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
	fFile.precision(14);

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
	fFile.precision(14);
	
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
void Correlation::xi_A_delta_delta_MockJMc(void) {

	std::cout << "\n\n\n\n  ------ xi_A_delta_delta_MockJMc ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_A_delta_delta.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	loadDataForest(pathForest__);
	if (mocksNoNoiseNoCont) {
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
	const unsigned int max = 200.;
	const unsigned int nbBin = int(max);
	const unsigned int nbBinM = 50.;

	const double maxPow2      = max*max;

	///// Arrays for data
	double dataMu[nbBin][nbBinM][7];
	double data2D[nbBin][nbBin][7];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<7; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}
	

	///// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::srand (42);

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

				const double r1 = v_r__[f1][i1];
				const double w1 = v_w__[f1][i1];
				const double d1 = v_d__[f1][i1];
				const double z1 = v_z__[f1][i1];

				for (unsigned int i2=0; i2<nb2; i2++) {

					const double rParral = fabs(v_r__[f2][i2]-r1);
					if (rParral>=max) continue;
					const double distTotPow2 = distTransPow2 + rParral*rParral;

					const double w1w2     = w1*v_w__[f2][i2];
					const double d1d2     = d1*v_d__[f2][i2];
					const double w1w2z1z2 = w1w2*(z1+v_z__[f2][i2]);

					if (distTotPow2 < maxPow2) {
						const double distTot    = sqrt(distTotPow2);
						const double mu         = rParral/distTot;
						const unsigned int idx  = int(distTot);
						const unsigned int idxM = int(mu*50.);

						dataMu[idx][idxM][0] += w1w2*d1d2;
						dataMu[idx][idxM][1] += w1w2*d1d2*d1d2;
						dataMu[idx][idxM][2] += w1w2*distTot;
						dataMu[idx][idxM][3] += w1w2*mu;
						dataMu[idx][idxM][4] += w1w2z1z2;
						dataMu[idx][idxM][5] += w1w2;
						dataMu[idx][idxM][6] ++;
					}
					
					///// Fill the histogramm of xi(r_{perp}, r_{paral}
					const unsigned int rParralBinIdx = int(rParral);
					data2D[rPerpBinIdx][rParralBinIdx][0] += w1w2*d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][1] += w1w2*d1d2*d1d2;
					data2D[rPerpBinIdx][rParralBinIdx][2] += w1w2*rPerp;
					data2D[rPerpBinIdx][rParralBinIdx][3] += w1w2*rParral;
					data2D[rPerpBinIdx][rParralBinIdx][4] += w1w2z1z2;
					data2D[rPerpBinIdx][rParralBinIdx][5] += w1w2;
					data2D[rPerpBinIdx][rParralBinIdx][6] ++;
				}
			}
		}
	}

	std::cout << "\n  Saving\n" << std::endl;

	std::ofstream fFile;
	long double sumZZZ = 0.;
	long double sumWeg = 0.;
	long double sumOne = 0.;

	///// Save the 2D cross-correlation
	std::string pathToSave = pathToSave__;
	pathToSave += "xi_A_delta_delta_2D_";
	pathToSave += forest__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

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
	pathToSave += forest__;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

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
void Correlation::xi_delta_QSO_MockJMc(bool doBootstraps/*=False*/, unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	//// QSO
	loadDataQ1();
	//// Forest
	loadDataForest(pathForest__,doBootstraps,bootIdx);
	if ( (doBootstraps || mocksNoNoiseNoCont) && !nicolasEstimator__) {
		removeFalseCorrelations();
	}
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();

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
	const double max          =  200.; // 200.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100.;

	const double maxPow2      = max*max;
	//// get an array for nb of pairs for the forest
	double a_nbPairs[nbForest_];
	for (unsigned int i=0; i<nbForest_; i++) {
		a_nbPairs[i] = 0.;
	}


	//// Arrays for data
	double dataMu[nbBin][nbBinM][7];
	double data2D[nbBinX][nbBinY][7];

	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBinM; j++) {
			for (unsigned int k=0; k<7; k++) {
				dataMu[i][j][k] = 0.;
			}
		}
	}
	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			for (unsigned int k=0; k<7; k++) {
				data2D[i][j][k] = 0.;
			}
		}
	}


	//// Vectors of randomized positions in cell
	std::vector<double> v_raRandForest;
	std::vector<double> v_deRandForest;
	std::vector<double> v_raRandQSO;
	std::vector<double> v_deRandQSO;
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::srand (42);

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
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}
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
			const double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
			}
			
			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// distance of QSO
			const double rQSO = v_rQ1__[q];
			
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

			//// Transvers distance between the qso and the lya
			const double distTransQsoLya = sqrt(distTransQsoLyaPow2);
			const unsigned int rPerpBinIdx = int( distTransQsoLya );
			const double zQSO = z2;
			
			//// Loops over all pixels of the forest
			for (unsigned int i=0; i<nbPixel; i++) {

				const double distP = v_r__[f][i] - rQSO;

				if (fabs(distP) >= max) continue;

				const double distTotPow2 = distTransQsoLyaPow2 + distP*distP;
				const double w   = v_w__[f][i];
				const double d   = v_d__[f][i];
				const double wd  = w*d;
				const double wdd = wd*d;
				const double wZpair = w*(zQSO+v_z__[f][i]);

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
	if (doBootstraps)  prefix += "subsampling";
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
	if (doBootstraps || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);

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
	if (doBootstraps || shuffleForest || shuffleQSO || randomQSO) pathToSave += prefix;
	pathToSave += ".txt";
	std::cout << "\n  " << pathToSave << std::endl;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(14);
	
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
			fFile << std::endl;
		}
	}
	fFile.close();

	
	if (!doBootstraps && !shuffleForest && !shuffleQSO && !randomQSO) {

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
	
	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_distortionMatrix ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_distortionMatrix.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	//// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	//// Forest
	loadDataForest(pathForest__);

	//// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_d__.clear();
	v_z__.clear();
	v_idx__.clear();
	v_lObs__.clear();
	v_nb__.clear();

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
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::srand (42);

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
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
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

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
			}
			
			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// distance of QSO
			const double rQSO = v_rQ1__[q];
			
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

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
	fFile.precision(14);
	
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
	
	std::cout << "\n\n\n\n  ------ xi_delta_QSO_MockJMc_distortionMatrix_1D ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_delta_QSO_distortionMatrix.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	//// QSO
	loadDataQ1();
	if (nbQ1__==0) return;
	//// Forest
	loadDataForest(pathForest__);

	//// Empty useless vectors
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_d__.clear();
	v_z__.clear();
	v_idx__.clear();
	v_lObs__.clear();
	v_nb__.clear();

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
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::srand (42);

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
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
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

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
			}
			
			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// distance of QSO
			const double rQSO = v_rQ1__[q];
			
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;
			
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
	fFile.precision(14);
	
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
	loadDataForest(pathForest__);

	//// Create the conversion table from redshift to distance
	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	//// Get the distance if it was this metal
	std::vector< std::vector< double > > v_r_metal(v_r__);
	for (unsigned int i=0; i<v_r_metal.size(); i++) {
		for (unsigned int j=0; j<v_r_metal[i].size(); j++) {
			v_r_metal[i][j] = hConvertRedshDist->Interpolate( v_lObs__[i][j]/lambdaFrMetal-1. );
		}
	}

	delete cosmo;
	v_CosDe__.clear();
	v_SinDe__.clear();
	v_zz__.clear();
	v_d__.clear();
	v_lRF__.clear();
	v_lObs__.clear();
	v_nb__.clear();


	//// Get the monopole, quadrupol, exadecapol
	std::string pathToLoad = PATHTOWORK;
	pathToLoad += PATHTOCAMB;
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
	
	const unsigned int nbBinCAMB = data_x.size();
	//// Set to zero the last pixel
	const double maxDistCAMB = data_x[nbBinCAMB-1];
	data_x.push_back(data_x[nbBinCAMB-1]+1.);
	data_xi0.push_back(0.);
	data_xi2.push_back(0.);
	data_xi4.push_back(0.);



	//// Constants:
	//// The space between bins is of 10 Mpc.h^-1
	const double max          =  200.; // 200.;

	const unsigned int nbBin  = int(max);
	const unsigned int nbBinX = nbBin;
	const unsigned int nbBinY = 2*nbBin;
	const unsigned int nbBinM = 100.;
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
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::srand (42);

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
		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raRandQSO[i] = v_raQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deRandQSO[i] = v_deQ1__[i] + sizeCell__*(1.*rand()/RAND_MAX-0.5);
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

		for (unsigned int q=0; q<nbQ1__; q++) {

			double x2 = v_raQ1__[q];
			double y2 = v_deQ1__[q];
			const double z2 = v_zzQ1__[q];

			//// Not in the same line of sight
			if (x==x2 && y==y2) continue;

			//// If random position in the cell
			if (randomPositionOfQSOInCellNotBeforeCorrelation__) {
				x2 = v_raRandQSO[q];
				y2 = v_deRandQSO[q];
			}
			
			//// Get the r_perp distance at the poxer of two
			const double distTransQsoLyaPow2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			if (distTransQsoLyaPow2 >= maxPow2) continue;

			//// distance of QSO
			const double rQSO = v_rQ1__[q];
			
			//// Distance between the qso and the first pixel
			if ( firstPixel-rQSO >= max) continue;
			
			//// Distance between the qso and the last pixel
			if ( lastPixel-rQSO <= -max) continue;

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
				if (distTotMetal<maxDistCAMB) idxBinCAMB = int(distTotMetal);
				const double wxi0 = w*data_xi0[idxBinCAMB];
				const double wxi2 = w*data_xi2[idxBinCAMB];
				const double wxi4 = w*data_xi4[idxBinCAMB];

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
	fFile.precision(14);
	
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
	fFile.precision(14);
	
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
void Correlation::xi_QSO_QSO_MockJMc(bool doBootstraps/*=false*/, unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n\n\n  ------ xi_QSO_QSO_MockJMc ------" << std::endl;
	std::string command = "  python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/xi_Q_Q.py";
	command += commandEnd__;
	std::cout << command << "\n" << std::endl;

	loadDataQ1();

	//// If doing with random
	std::stringstream convert;
	convert << bootIdx;
	const std::string strBootIdx = convert.str();

	std::vector<double> dataX(v_raQ1__);
	std::vector<double> dataY(v_deQ1__);
	std::vector<double> dataR(v_rQ1__);
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
			const double r = v_rQ1__[nbQSO];
			const double z = v_zzQ1__[nbQSO];

			bool cont = false;

			//// Pick a new one if already picked
			for (unsigned int i=0; i<nbQSO; i++) {
				if (v_raQ1__[i]==x && v_deQ1__[i]==y && v_rQ1__[i]==r) {
					cont = true;
					continue;
				}
			}
			if (cont) continue;

			v_raQ1__[nbQSO] = x;
			v_deQ1__[nbQSO] = y;
			v_rQ1__[nbQSO]  = r;
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


	///// Vectors of randomized positions in cell
	if (randomPositionOfQSOInCellNotBeforeCorrelation__) {

		std::srand (42);

		for (unsigned int i=0; i<nbQ1__; i++) {
			v_raQ1__[i] += sizeCell__*(1.*rand()/RAND_MAX-0.5);
			v_deQ1__[i] += sizeCell__*(1.*rand()/RAND_MAX-0.5);
		}
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
		
		const double x1 = v_raQ1__[q1];
		const double y1 = v_deQ1__[q1];
		const double r1 = v_rQ1__[q1];
		const double z1 = v_zzQ1__[q1];

		for (unsigned int q2=0; q2<q1; q2++) {

			//// s_perp
			const double rPerpPow2 = (x1-v_raQ1__[q2])*(x1-v_raQ1__[q2]) + (y1-v_deQ1__[q2])*(y1-v_deQ1__[q2]);
			if (rPerpPow2>=maxPow2) continue;

			//// s_paral
			const double rParaPow2 = (r1-v_rQ1__[q2])*(r1-v_rQ1__[q2]);
			if (rParaPow2>=maxPow2) continue;

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
	fFile.precision(14);

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
	fFile.precision(14);

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
			
				//// s_perp
				const double rPerpPow2 = (x1-v_raQ1__[q2])*(x1-v_raQ1__[q2]) + (y1-v_deQ1__[q2])*(y1-v_deQ1__[q2]);
				if (rPerpPow2>=maxPow2) continue;
				
				//// s_paral
				const double rParaPow2 = (r1-v_rQ1__[q2])*(r1-v_rQ1__[q2]);
				if (rParaPow2>=maxPow2) continue;
				
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
		fFile.precision(14);

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
		fFile.precision(14);

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
	if (idxCommand_[0]==15 || idxCommand_[0]==16 || idxCommand_[0]==21 || mockBox__) {
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

		if ( !mockBox__ && ra==0. && de==0. ) {
			std::cout << "  Correlation::loadDataQ1::  ERROR:  ra==0. && de==0. , RA = " << ra << " , Dec = " << de << " , idx = " << i << std::endl;
			continue;
		}
		if ( selectWindowRedshift && (zz<minRedshiftQSO || zz>maxRedshiftQSO) ) continue;

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
		if (mockBox__) v_rQ1__.push_back( zz );
		else v_rQ1__.push_back( hConvertRedshDist->Interpolate(zz) );
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
		v_rQ2__.push_back( hConvertRedshDist->Interpolate(zz) );
	}

	nbQ2__ = v_raQ2__.size();
	std::cout << "  number of good   Q2    = " << nbQ2__ << std::endl;

	delete hConvertRedshDist;

	return;
}
void Correlation::loadDataForest(std::string pathToFits,bool doBootstraps/*=false*/, unsigned int bootIdx/*=0*/) {

	std::cout << "\n\n  ------ load Data Forest ------" << std::endl;

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

	unsigned int regionMap[nrows];
	//// If doing sub-sampling, geting the array to know if the forest is in the given region
	if (doBootstraps) {

		//// Init the values
		for (unsigned int i=0; i<nrows; i++) {
			regionMap[i] = 100000;
		}

		//// For the cross-correlation
		std::string pathToLoad = pathToSave__;
		pathToLoad += "xi_delta_QSO_map_";
		pathToLoad += forest__;
		pathToLoad += "_";
		pathToLoad += QSO__;
		pathToLoad += ".txt";

		LymanForest* lymanForestObject = new LymanForest(pathToLoad, C_NBSUBSAMPLES, C_RA_SEPERATION_NGC_SGC, mockJMC__);
		lymanForestObject->GetRegionArray(regionMap);
		lymanForestObject->PrintRegionDetail(bootIdx);
		delete lymanForestObject;
		std::cout << "\n\n" << std::endl;
	}

	long unsigned int nbCutted[5] = {0}; 
	double meanDelta[3] = {0.};
	double meanDelta_Nicolas[3] = {0.};
	double meanDelta_Nicolas_overForests[4] = {0.};

	//// Load data
	for (unsigned int i=0; i<nbForest_; i++) {

		if (doBootstraps && regionMap[i]!=bootIdx) continue;

		//// Variables for data in FITS
		double ra = 0.;
		double de = 0.;
		double zz = 0.;
		double alpha1, alpha2, beta2;
		double meanForestLambdaRF = 0.;

		double DELTA[nbBinRFMax__];
		double DELTA_WEIGHT[nbBinRFMax__];
		double LAMBDA_OBS[nbBinRFMax__];
		double LAMBDA_RF[nbBinRFMax__];
		double NORM_FLUX[nbBinRFMax__];
		double NORM_FLUX_IVAR[nbBinRFMax__];
		double FLUX_DLA[nbBinRFMax__];

		fits_read_col(fitsptrSpec,TDOUBLE, 4,i+1,1,1,NULL,&ra,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 5,i+1,1,1,NULL,&de,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 6,i+1,1,1,NULL,&zz,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 8,i+1,1,1,NULL,&meanForestLambdaRF,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 9,i+1,1,1,NULL,&alpha1,                NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 11,i+1,1,1,NULL,&alpha2,               NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 12,i+1,1,1,NULL,&beta2,                NULL,&sta);

		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMax__,NULL, &LAMBDA_OBS,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 14,i+1,1,nbBinRFMax__,NULL, &LAMBDA_RF,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 15,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX_IVAR,  NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 17,i+1,1,nbBinRFMax__,NULL, &FLUX_DLA,        NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 18,i+1,1,nbBinRFMax__,NULL, &DELTA,           NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 20,i+1,1,nbBinRFMax__,NULL, &DELTA_WEIGHT,    NULL,&sta);

		if (!mocksNoNoiseNoCont && ( (alpha2 == alphaStart__ && beta2 == betaStart__) || (fabs(alpha2)>=maxAlpha__-0.5) || (fabs(beta2)>=maxBeta__-0.05) ) ) {
			if (alpha2 == alphaStart__ && beta2 == betaStart__) nbCutted[0] ++;
			else if (fabs(alpha2)>=maxAlpha__-0.5) nbCutted[1] ++;
			else if (fabs(beta2)>=maxBeta__-0.05) nbCutted[2] ++;
			continue;
		}

		//// If a reobs
		if (alpha1==isReobsFlag__) continue;

		bool templateHasNegative = false;
		double tmp_meanDelta[3] = {0.};
		long double meanNeeded[5] = {0.};
		std::vector< double > v_tmp_r;
		std::vector< double > v_tmp_d;
		std::vector< double > v_tmp_w;
		std::vector< double > v_tmp_z;
		std::vector< double > v_tmp_lRF;
		std::vector< double > v_tmp_lObs;

		for (unsigned int j=0; j<nbBinRFMax__; j++) {
			if (DELTA_WEIGHT[j]>0. && NORM_FLUX_IVAR[j]>0. && FLUX_DLA[j]>=C_DLACORR && LAMBDA_OBS[j]>=lambdaObsMin__ && LAMBDA_OBS[j]<lambdaObsMax__ && LAMBDA_RF[j]>=lambdaRFMin__ && LAMBDA_RF[j]<lambdaRFMax__) {

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
				v_tmp_lRF.push_back(LAMBDA_RF[j]);
				v_tmp_lObs.push_back(LAMBDA_OBS[j]);

				//// Get if the template is even one time negative
				if ( alpha2+beta2*(LAMBDA_RF[j]-meanForestLambdaRF) <= 0.) templateHasNegative = true;

				//// Get Nb Pixel in forest
				tmp_meanDelta[0] += DELTA_WEIGHT[j]*DELTA[j];
				tmp_meanDelta[1] += DELTA_WEIGHT[j];
				tmp_meanDelta[2] ++;

				//// For Nicolas's estimator
				if (nicolasEstimator__) {
					meanNeeded[0] += DELTA_WEIGHT[j]*DELTA[j];
					meanNeeded[1] += DELTA_WEIGHT[j]*LAMBDA_RF[j];
					meanNeeded[2] += DELTA_WEIGHT[j]*LAMBDA_RF[j]*LAMBDA_RF[j];
					meanNeeded[3] += DELTA_WEIGHT[j]*LAMBDA_RF[j]*DELTA[j];
					meanNeeded[4] += DELTA_WEIGHT[j];
				}
			}
		}

		const unsigned int tmp_nb = v_tmp_r.size();
		if (tmp_nb<C_MIN_NB_PIXEL || templateHasNegative ) {
			if (tmp_nb<C_MIN_NB_PIXEL) nbCutted[3] ++;
			else if (templateHasNegative) nbCutted[4] ++;
			continue;
		}

		//// Get Nb Pixel in forest
		meanDelta[0] += tmp_meanDelta[0];
		meanDelta[1] += tmp_meanDelta[1];
		meanDelta[2] += tmp_meanDelta[2];

		//// For Nicolas's estimator
		if (nicolasEstimator__) {
			const long double meanDelta   = meanNeeded[0]/meanNeeded[4];
			const long double meanLambda  = meanNeeded[1]/meanNeeded[4];
			const long double numerator   = meanNeeded[3]-meanLambda*meanNeeded[0];
			const long double denominator = meanNeeded[2]-meanLambda*meanLambda*meanNeeded[4];
			const long double coef        = numerator/denominator;
			for (unsigned int j=0; j<tmp_nb; j++) {
				v_tmp_d[j] -= meanDelta+(v_tmp_lRF[j]-meanLambda)*coef;


				//// Get Nb Pixel in forest
				meanDelta_Nicolas[0] += v_tmp_w[j]*v_tmp_d[j];
				meanDelta_Nicolas[1] += v_tmp_w[j];
				meanDelta_Nicolas[2] ++;
			}
			meanDelta_Nicolas_overForests[0] += meanDelta;
			meanDelta_Nicolas_overForests[1] += coef;
			meanDelta_Nicolas_overForests[2] += meanLambda-meanForestLambdaRF;
			meanDelta_Nicolas_overForests[3] ++;
		}

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

		std::vector< unsigned int > v_tmp_nb(tmp_nb,0);
		v_nb__.push_back(v_tmp_nb);
	}

	nbForest_ = v_ra__.size();
	std::cout << "  number of good forest   = " << nbForest_ << std::endl;
	for (unsigned int i=0; i<5; i++) std::cout << "  lost n°" << i << " = " << nbCutted[i] << std::endl;

	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;

	if (nicolasEstimator__) {
		std::cout << "  < delta_Nicolas >              = " << meanDelta_Nicolas[0]/meanDelta_Nicolas[1] << std::endl;
		std::cout << "  sum(w_i)_Nicolas               = " << meanDelta_Nicolas[1]              << std::endl;
		std::cout << "  nb pixel_Nicolas               = " << (long long unsigned int)meanDelta_Nicolas[2]              << std::endl;
		std::cout << "  < <delta> >                    = " << meanDelta_Nicolas_overForests[0]/meanDelta_Nicolas_overForests[3] << std::endl;
		std::cout << "  < coef >                       = " << meanDelta_Nicolas_overForests[1]/meanDelta_Nicolas_overForests[3] << std::endl;
		std::cout << "  < lambda > - < lambda >_before = " << meanDelta_Nicolas_overForests[2]/meanDelta_Nicolas_overForests[3] << std::endl;
	}

	fits_close_file(fitsptrSpec,&sta);

	delete hConvertRedshDist;

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

	long unsigned int nbCutted[5] = {0};
	double meanDelta[3] = {0.};
	double meanDelta_Nicolas[3] = {0.};
	double meanDelta_Nicolas_overForests[4] = {0.};

	//// Load data
	for (unsigned int i=0; i<nbForest2__; i++) {

		//// Variables for data in FITS
		double ra = 0.;
		double de = 0.;
		double zz = 0.;
		double alpha2, beta2;
		double meanForestLambdaRF = 0.;

		double LAMBDA_OBS[nbBinRFMaxDelta2__];
		double LAMBDA_RF[nbBinRFMaxDelta2__];
		double NORM_FLUX[nbBinRFMaxDelta2__];
		double NORM_FLUX_IVAR[nbBinRFMaxDelta2__];
		double FLUX_DLA[nbBinRFMaxDelta2__];
		double DELTA[nbBinRFMaxDelta2__];
		double DELTA_WEIGHT[nbBinRFMaxDelta2__];

		fits_read_col(fitsptrSpec,TDOUBLE, 4,i+1,1,1,NULL,&ra,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 5,i+1,1,1,NULL,&de,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 6,i+1,1,1,NULL,&zz,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 8, i+1,1,1,NULL,&meanForestLambdaRF,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 11,i+1,1,1,NULL,&alpha2,               NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 12,i+1,1,1,NULL,&beta2,                NULL,&sta);

		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMaxDelta2__,NULL, &LAMBDA_OBS,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 14,i+1,1,nbBinRFMaxDelta2__,NULL, &LAMBDA_RF,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 15,i+1,1,nbBinRFMaxDelta2__,NULL, &NORM_FLUX,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMaxDelta2__,NULL, &NORM_FLUX_IVAR,  NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 17,i+1,1,nbBinRFMaxDelta2__,NULL, &FLUX_DLA,        NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 18,i+1,1,nbBinRFMaxDelta2__,NULL, &DELTA,           NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 20,i+1,1,nbBinRFMaxDelta2__,NULL, &DELTA_WEIGHT,    NULL,&sta);

		if (!mocksNoNoiseNoCont && ( (alpha2 == alphaStart__ && beta2 == betaStart__) || (fabs(alpha2)>=maxAlpha__-0.5) || (fabs(beta2)>=maxBeta__-0.05) ) ) {
			if (alpha2 == alphaStart__ && beta2 == betaStart__) nbCutted[0] ++;
			else if (fabs(alpha2)>=maxAlpha__-0.5) nbCutted[1] ++;
			else if (fabs(beta2)>=maxBeta__-0.05) nbCutted[2] ++;
			continue;
		}

		bool templateHasNegative = false;
		long double meanNeeded[5] = {0.};
		double tmp_meanDelta[3] = {0.};
		std::vector< double > v_tmp_r;
		std::vector< double > v_tmp_d;
		std::vector< double > v_tmp_w;
		std::vector< double > v_tmp_z;
		std::vector< double > v_tmp_lRF;
		std::vector< double > v_tmp_lObs;

		for (unsigned int j=0; j<nbBinRFMaxDelta2__; j++) {

			if (DELTA_WEIGHT[j]>0. && NORM_FLUX_IVAR[j]>0. && FLUX_DLA[j]>=C_DLACORR && LAMBDA_OBS[j]>=lambdaObsMin__ && LAMBDA_OBS[j]<lambdaObsMax__ && LAMBDA_RF[j]>=lambdaRFMinDelta2__ && LAMBDA_RF[j]<lambdaRFMaxDelta2__) {

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
				if ( alpha2+beta2*(LAMBDA_RF[j]-meanForestLambdaRF) <= 0.) templateHasNegative = true;

				const double zi = LAMBDA_OBS[j]/lambdaRFLineDelta2__ -1.;
				v_tmp_r.push_back(hConvertRedshDist->Interpolate(zi));
				v_tmp_d.push_back(DELTA[j]);
				v_tmp_w.push_back(DELTA_WEIGHT[j]);
				v_tmp_z.push_back(zi);
				v_tmp_lRF.push_back(LAMBDA_RF[j]);
				v_tmp_lObs.push_back(LAMBDA_OBS[j]);

				tmp_meanDelta[0] += DELTA_WEIGHT[j]*DELTA[j];
				tmp_meanDelta[1] += DELTA_WEIGHT[j];
				tmp_meanDelta[2] ++;

				//// For Nicolas's estimator
				if (nicolasEstimator__) {
					meanNeeded[0] += DELTA_WEIGHT[j]*DELTA[j];
					meanNeeded[1] += DELTA_WEIGHT[j]*LAMBDA_RF[j];
					meanNeeded[2] += DELTA_WEIGHT[j]*LAMBDA_RF[j]*LAMBDA_RF[j];
					meanNeeded[3] += DELTA_WEIGHT[j]*LAMBDA_RF[j]*DELTA[j];
					meanNeeded[4] += DELTA_WEIGHT[j];
				}
			}
		}

		const unsigned int tmp_nb = v_tmp_r.size();
		if (tmp_nb<C_MIN_NB_PIXEL || templateHasNegative ) {
			if (tmp_nb<C_MIN_NB_PIXEL) nbCutted[3] ++;
			else if (templateHasNegative) nbCutted[4] ++;
			continue;
		}

		//// Get Nb Pixel in forest
		meanDelta[0] += tmp_meanDelta[0];
		meanDelta[1] += tmp_meanDelta[1];
		meanDelta[2] += tmp_meanDelta[2];

		//// For Nicolas's estimator
		if (nicolasEstimator__) {
			const long double meanDelta   = meanNeeded[0]/meanNeeded[4];
			const long double meanLambda  = meanNeeded[1]/meanNeeded[4];
			const long double numerator   = meanNeeded[3]-meanLambda*meanNeeded[0];
			const long double denominator = meanNeeded[2]-meanLambda*meanLambda*meanNeeded[4];
			const long double coef        = numerator/denominator;
			for (unsigned int j=0; j<tmp_nb; j++) {
				v_tmp_d[j] -= meanDelta+(v_tmp_lRF[j]-meanLambda)*coef;


				//// Get Nb Pixel in forest
				meanDelta_Nicolas[0] += v_tmp_w[j]*v_tmp_d[j];
				meanDelta_Nicolas[1] += v_tmp_w[j];
				meanDelta_Nicolas[2] ++;
			}
			meanDelta_Nicolas_overForests[0] += meanDelta;
			meanDelta_Nicolas_overForests[1] += coef;
			meanDelta_Nicolas_overForests[2] += meanLambda-meanForestLambdaRF;
			meanDelta_Nicolas_overForests[3] ++;
		}

		ra = ra*C_DEGTORAD;
		de = de*C_DEGTORAD;

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
	}

	fits_close_file(fitsptrSpec,&sta);
	
	nbForest2__ = v_raDelta2__.size();
	std::cout << "  number of good forest   = " << nbForest2__ << std::endl;
	for (unsigned int i=0; i<5; i++) std::cout << "  lost n°" << i << " = " << nbCutted[i] << std::endl;

	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;

	if (nicolasEstimator__) {
		std::cout << "  < delta_Nicolas >              = " << meanDelta_Nicolas[0]/meanDelta_Nicolas[1] << std::endl;
		std::cout << "  sum(w_i)_Nicolas               = " << meanDelta_Nicolas[1]              << std::endl;
		std::cout << "  nb pixel_Nicolas               = " << (long long unsigned int)meanDelta_Nicolas[2]              << std::endl;
		std::cout << "  < <delta> >                    = " << meanDelta_Nicolas_overForests[0]/meanDelta_Nicolas_overForests[3] << std::endl;
		std::cout << "  < coef >                       = " << meanDelta_Nicolas_overForests[1]/meanDelta_Nicolas_overForests[3] << std::endl;
		std::cout << "  < lambda > - < lambda >_before = " << meanDelta_Nicolas_overForests[2]/meanDelta_Nicolas_overForests[3] << std::endl;
	}

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
	const unsigned int nbLoop = 50;

	//// Remove <delta> vs. lambda_OBS
	const unsigned int nbBin = nbBinlambdaObs__;
	const double min = lambdaObsMin__;
	const double max = lambdaObsMax__;
	
	/*
	//// Remove <delta> vs. lambda_RF
	const unsigned int nbBin = int(lambdaRFMax__-lambdaRFMin__)+6;
	const double min = lambdaRFMin__-3.;
	const double max = lambdaRFMax__+3.;
	*/
	TH1D* histo = new TH1D("histo","",nbBin,min,max);

	//// Start the loop
	for (unsigned int lp=0; lp<nbLoop; lp++) {

		double data[nbBin][4];
		for (unsigned int i=0; i<nbBin; i++) {
			for (unsigned int j=0; j<4; j++) {
				data[i][j] = 0.;
			}
		}

		//// Find < delta >
		for (unsigned int f=0; f<nbForest_; f++) {

			const unsigned int nb = v_nbPixelDelta1__[f];	
			for (unsigned int i=0; i<nb; i++) {

				double coef = 1.;
				if (!firstPass) coef = v_nb__[f][i];
				if (coef==0.) continue; 

				const double l  = v_lObs__[f][i];
				//const double l    = v_lRF__[f][i];
				const double w    = v_w__[f][i];
				const double d    = v_d__[f][i];
				const double wd   = w*d;

				const unsigned int idx = int(l-min);
				data[idx][0] += coef*wd;
				data[idx][1] += coef*w;
				data[idx][2] += coef*wd*d;
				data[idx][3] ++;
			}
		}

		long double mean0[2] = {0.};

		//// Put it into an histo
		for (unsigned int i=0; i<nbBin; i++) {

			if (data[i][3]==0.) {
				histo->SetBinContent(i+1,0.);
				histo->SetBinError(i+1,0.);
				continue;
			}

			mean0[0] += data[i][0];
			mean0[1] += data[i][1];

			const double mean = data[i][0]/data[i][1];
			const double err  = sqrt( (data[i][2]/data[i][1]-mean*mean)/data[i][3] );
			
			histo->SetBinContent(i+1,mean);
			histo->SetBinError(i+1,err);
		}

		if (lp==0 || lp==nbLoop-1) std::cout << "  " << lp << " " << mean0[0]/mean0[1] << std::endl;


		// Fill empty bins
		double xxx1_value_m2 = 0.;
		double xxx2_value_m2 = 0.;
		double yyy1_value_m2 = 0.;
		double yyy2_value_m2 = 0.;
		double value_m2 = 0.;
		unsigned int nbEmptyPixels_m2 = 0;
		for (unsigned int i=0; i<nbBin; i++) {
			
			// If empty and not the last pixel
			if (histo->GetBinError(i+1) == 0. && i!=nbBin-1) {
				if (nbEmptyPixels_m2 == 0 && i!=0) {
					xxx1_value_m2 = histo->GetBinCenter(i);
					yyy1_value_m2 = histo->GetBinContent(i);
					value_m2 = histo->GetBinContent(i);
				}
				nbEmptyPixels_m2 ++;
			}
			else {
				
				// If first pixel not empty after some empty pixels
				if (nbEmptyPixels_m2!=0) {
	
					double a = 0.;
					double b = 0.;					
					
					// Find the mean value between the two not empty pixels
					if (value_m2!=0. && histo->GetBinError(i+1)!=0.) {
						// Linear extrapolation
						xxx2_value_m2 = histo->GetBinCenter(i+1);
						yyy2_value_m2 = histo->GetBinContent(i+1);
						a = yyy1_value_m2 - yyy2_value_m2;
						b = yyy2_value_m2*xxx1_value_m2 - yyy1_value_m2*xxx2_value_m2;
						if (xxx1_value_m2-xxx2_value_m2 != 0.) {
							a /= (xxx1_value_m2-xxx2_value_m2);
							b /= (xxx1_value_m2-xxx2_value_m2);
						}
						value_m2 = (value_m2 + histo->GetBinContent(i+1))/2.;
					}
					
					// Set all the empty pixels 
					for (unsigned int j=0; j<nbEmptyPixels_m2; j++) {
						double tmp_value = a*histo->GetBinCenter(i-j) + b;
						if (histo->GetBinError(i+1)==0.) tmp_value = value_m2;
						histo->SetBinContent(i-j, tmp_value); //value_m2
						//histo->SetBinError(i-j, maxError_m2);
					}
					
					// If the last pixel is also empty
					if (i==nbBin-1) {
						histo->SetBinContent(i+1, value_m2);
						//histo->SetBinError(i+1, maxError_m2);
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

		//// Remove the correlation
		for (unsigned int f=0; f<nbForest_; f++) {
			
			const unsigned int nb = v_nbPixelDelta1__[f];	
			for (unsigned int i=0; i<nb; i++) {
	
				const double l = v_lObs__[f][i];
				//const double l = v_lRF__[f][i];
	
				v_d__[f][i] -= histo->Interpolate(l);
			}
		}
	
		//std::stringstream convert;
		//convert << lp;
		//R_plot1D(histo, convert.str());

	}

	delete histo;
}





