//===================================================================================
//
//         FILE: GetDelta.cpp
//
//        USAGE: 
//
//  DESCRIPTION: ---
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

#include "GetDelta.h"

#include "../../../Root/Library/RootHistoFunctions.h"
#include "../../../Cpp/Library/mathFunctions.h"
#include "../../../Constants/constants.h"
#include "../../../Constants/globalValues.h"

#include <fstream>
#include <iostream>	// std::cout
#include <sstream>	//stringstream
#include <stdlib.h>	//atoi
#include <unistd.h>	//getopt
#include <cstdlib>	//std::rand, std::srand
#include <stdio.h>	//send a command
#include <algorithm>	//remove_if
#include <algorithm>
#include <string>
#include <iostream>
#include <cctype>



/// ROOT
#include "fitsio.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TProfile.h"

const unsigned int nbBinlambdaObs__  = int(lambdaObsMax__-lambdaObsMin__);
const double lambdaRFTemplateMin__ = lambdaRFMin__-3.;
const double lambdaRFTemplateMax__ = lambdaRFMax__+3.;
long  Bit16 = 1;
const double findPDF__   = false;
const bool doVetoLines__ = true;


GetDelta::GetDelta(int argc, char** argv) {

	isTest__   = false;
	withRSD__  = true;
	noMetals__ = true;

	std::cout << std::scientific;
	std::cout.precision(std::numeric_limits<double>::digits10);
	if (findPDF__) defineHistos();

	// Init the bit
	Bit16 = Bit16 <<16;

	/// 0: Box, 1: Simulation
	if (argc<3) return;
	const std::string box_idx  = argv[1];
	const std::string sim_idx  = argv[2];


	///
	pathToDataQSO__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1575/fits/spectra-785";
	pathToDataQSO__ += box_idx;
	pathToDataQSO__ += "-";
	pathToDataQSO__ += sim_idx;
	pathToDataQSO__ += ".fits";
	///
	pathToDataForest__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_00";
	pathToDataForest__ += box_idx;
	pathToDataForest__ += "/Simu_00";
	pathToDataForest__ += sim_idx;
	pathToDataForest__ += "/Raw/mocks-*";
	///
	pathToSave__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1575_with_good_metals/Box_00";
	pathToSave__ += box_idx;
	pathToSave__ += "/Simu_00";
	pathToSave__ += sim_idx;
	pathToSave__ += "/Data";
	if (noMetals__) pathToSave__ += "_no_metals";
	pathToSave__ += "/";

	std::cout << "\n"   << std::endl;
	std::cout << "  " << pathToDataQSO__ << std::endl;
	std::cout << "  " << pathToDataForest__ << std::endl;
	std::cout << "  " << pathToSave__ << std::endl;
	std::cout << "\n"   << std::endl;

	/// QSO
	//GetQSO();

	std::cout << "\n"   << std::endl;

	/// Delta
	GetData();

	if (findPDF__) {

		for (unsigned int i=0; i<1; i++) {
			for (unsigned int j=0; j<1; j++) {
				TString a = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1575/fits/spectra-785";
				//TString a = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1588/spectra-785";
				a += i;
				a += "-";
				a += j;
				a += ".fits";
				pathToDataQSO__ = a;
				//GetPDF(0);
				GetPDF(2);
			}
		}

		//pathToDataQSO__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/old/fev16/spectra-highz.fits";
		//GetPDF(1);
		//pathToDataQSO__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1573/fits/spectra-7800-0.fits";
		//GetPDF(2);

		saveHistos();
	}
	std::cout << "\n\n" << std::endl;

}
GetDelta::~GetDelta(void) {
	if (findPDF__) {
		delete hFluxVsLambdaObs__;
		delete hFluxPDF__;
		delete hRedshift__;
		delete hFlux__;
	}
}

void GetDelta::GetData(void) {

	/// Fits file where to save
	TString TSfitsnameSpec2 = pathToSave__;
	TSfitsnameSpec2 += "delta.fits";
	std::cout << "  " << TSfitsnameSpec2 << std::endl;
	int sta2    = 0;
	fitsfile* fitsptrSpec2;
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READWRITE, &sta2); //READWRITE

	/// index of forest
	unsigned int forestIdx = 0;
	unsigned int nb_QSO_same_XY = 0;
	float X_before = 0.;
	float Y_before = 0.;

	/// Get the list 
	FILE *fp;
	char path[PATH_MAX];
	std::string command = "ls -tr " + pathToDataForest__;
	std::vector< std::string > listFiles;
	fp = popen(command.c_str(), "r");
	while (fgets(path, PATH_MAX, fp) != NULL) listFiles.push_back(path);

	/// Get nb of files
	const unsigned int nbFiles = listFiles.size();

	for (unsigned int fileIdx=0; fileIdx<nbFiles; fileIdx++) {

		std::string path = listFiles[fileIdx];
		path.erase(std::remove(path.begin(), path.end(), '\n'), path.end());
		path.erase(std::remove(path.begin(), path.end(), ' '), path.end());
		std::cout << path << std::endl;

		/// Get the file
		const TString TSfitsnameSpec = path;
		int sta = 0;
		fitsfile* fitsptrSpec;
		fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);

		/// Get the number of spectra
		int tmp_nbSpectra = 0;
		unsigned int nbSpectra = 0;
		if (isTest__) nbSpectra = 30000;
		else {
			fits_get_num_hdus(fitsptrSpec, &tmp_nbSpectra, &sta);
			nbSpectra = tmp_nbSpectra-1;
		}
	
		/// Set to the first HDU
		fits_movabs_hdu(fitsptrSpec, 2,  NULL, &sta);
		
		/// Get the data
		for (unsigned int f=0; f<nbSpectra; f++) {
	
			/// Move to next HDU
			if (f!=0) fits_movrel_hdu(fitsptrSpec, 1,  NULL, &sta);
	
			/// Get the number of pixels in this HDU
			long tmp_nbPixels = 0;
			fits_get_num_rows(fitsptrSpec, &tmp_nbPixels, &sta);
			const unsigned int tmp_nbPixels2 = tmp_nbPixels;
	
			/// Variable from old FITS
			float X = 0.;
			float Y = 0.;
			fits_read_key(fitsptrSpec,TFLOAT,"X",   &X,NULL,&sta);
			fits_read_key(fitsptrSpec,TFLOAT,"Y",   &Y,NULL,&sta);

			if (X==X_before && Y_before==Y) {
				nb_QSO_same_XY ++;
			}
			X_before = X;
			Y_before = Y;

			if (tmp_nbPixels2<C_MIN_NB_PIXEL) continue;

			/// Variable from old FITS
			unsigned int plate;
			unsigned int mjd;
			unsigned int fiberid;
			float Z = 0.;
			float LAMBDA_OBS[tmp_nbPixels2];
			float FLUX[tmp_nbPixels2];
			float IVAR[tmp_nbPixels2];
			long  AND_MASK[tmp_nbPixels2];
			float CONTINUUM[tmp_nbPixels2];
			float METAL_FLUX[tmp_nbPixels2];
			fits_read_key(fitsptrSpec,TINT,"PLATE",  &plate,NULL,&sta);
			fits_read_key(fitsptrSpec,TINT,"MJD",    &mjd,NULL,&sta);
			fits_read_key(fitsptrSpec,TINT,"FIBERID",&fiberid,NULL,&sta);
			fits_read_key(fitsptrSpec,TFLOAT,"ZQSO",&Z,NULL,&sta);
			fits_read_col(fitsptrSpec,TFLOAT, 1,1,1,tmp_nbPixels2,NULL, &LAMBDA_OBS,NULL,&sta);
			fits_read_col(fitsptrSpec,TFLOAT, 2,1,1,tmp_nbPixels2,NULL, &FLUX,      NULL,&sta);
			fits_read_col(fitsptrSpec,TFLOAT, 3,1,1,tmp_nbPixels2,NULL, &IVAR,      NULL,&sta);
			fits_read_col(fitsptrSpec,TLONG,  4,1,1,tmp_nbPixels2,NULL, &AND_MASK,  NULL,&sta);
			fits_read_col(fitsptrSpec,TFLOAT, 5,1,1,tmp_nbPixels2,NULL, &CONTINUUM, NULL,&sta);
			if (noMetals__) fits_read_col(fitsptrSpec,TFLOAT, 7,1,1,tmp_nbPixels2,NULL, &METAL_FLUX, NULL,&sta);

			/// Variables for new FITS
			double XX = X;
                        double YY = Y;
                        double ZZ = Z;
			double lambdaOBS[nbPixelsTemplate__];
			double norm_flux[nbPixelsTemplate__];
			double norm_flux_ivar[nbPixelsTemplate__];
			double continuum[nbPixelsTemplate__];
			double delta[nbPixelsTemplate__];
			double delta_weight[nbPixelsTemplate__];
			unsigned int nbPixel = 0;
			double normFactor[2] = {0.};
			double meanForestLRF[2] = {0.};
			double meanFluxForest[2] = {0.};
			const double oneOverOnePlusZ = 1./(1.+Z);
	
			for (unsigned int p=0; p<tmp_nbPixels2; p++) {
				/// bad pixels
				if (LAMBDA_OBS[p]<=0. || IVAR[p]<=0. || AND_MASK[p]>=Bit16 || std::isnan(IVAR[p]) ) continue;
	
				LAMBDA_OBS[p] = pow(10.,LAMBDA_OBS[p]);
				const double lambdaRFd = LAMBDA_OBS[p]*oneOverOnePlusZ;
	
				/// Pixel outside working region and remove pixels because of CCD and Sky lines 
				if (lambdaRFd<lambdaRFTemplateMin__ || lambdaRFd>lambdaRFNormaMax__ || LAMBDA_OBS[p]<lambdaObsMin__ || LAMBDA_OBS[p]>=lambdaObsMax__) continue;

				//// Remove veto lines
				bool isLine = false;
				if (doVetoLines__) {
					for (unsigned int k=0; k<nbVetoLines__; k++) {
						 if (LAMBDA_OBS[p]>=vetoLine__[2*k] && LAMBDA_OBS[p]<vetoLine__[2*k+1]) {
							isLine = true;
							break;
						}
					}
				}
				if (isLine) continue;
				if (noMetals__) FLUX[p] /= METAL_FLUX[p];
	
				/// Normalisation zone
				if (lambdaRFd>=lambdaRFNormaMin__ && lambdaRFd<lambdaRFNormaMax__) {
					normFactor[0] += FLUX[p];
					normFactor[1] ++;
				}
				/// Keep only the template, i.e. the zone with the forest plus some pixels
				else if (lambdaRFd>=lambdaRFTemplateMin__ && lambdaRFd<lambdaRFTemplateMax__) {
	
					lambdaOBS[nbPixel]      = LAMBDA_OBS[p];
					norm_flux[nbPixel]      = FLUX[p];
					norm_flux_ivar[nbPixel] = IVAR[p];
					continuum[nbPixel]      = CONTINUUM[p];
					nbPixel ++;
	
					if (lambdaRFd>=lambdaRFMin__ && lambdaRFd<lambdaRFMax__) {
						meanForestLRF[0] += lambdaRFd;
						meanForestLRF[1] ++;
						meanFluxForest[0] += FLUX[p];
						meanFluxForest[1] ++;
					}
				}
			}

			if (nbPixel>nbPixelsTemplate__)  std::cout << "  GetDelta::GetData::  ERROR: nbPixel>nbPixelsTemplate__ " << nbPixel << " " << nbPixelsTemplate__ << std::endl;
	
			/// Normalisation zone
			if (meanForestLRF[1]<C_MIN_NB_PIXEL) continue;
	
			/// Normalize the forest
			double norm = 1.;
			if (normFactor[1]>0. && normFactor[0]>0.) norm = normFactor[0]/normFactor[1];
	
			/// Normalise the data
			for (unsigned int i=0; i<nbPixel; i++) {
				delta[i]           = norm_flux[i]/continuum[i] -1.;
				delta_weight[i]    = norm_flux_ivar[i]*continuum[i]*continuum[i];
				norm_flux[i]      /= norm;
				norm_flux_ivar[i] *= norm*norm;
				continuum[i]      /= norm;
			}

			/// Get the mean lambda_RF
			meanForestLRF[0] /= meanForestLRF[1];
			/// Get the mean_flux_in_forest
			double alpha = CONVERT_FROM_FLUX_TO_ALPHA*meanFluxForest[0]/(meanFluxForest[1]*norm);

			/// Save data in second file
			fits_write_col(fitsptrSpec2,TINT,    1,forestIdx+1,1,1, &plate, &sta2);
			fits_write_col(fitsptrSpec2,TINT,    2,forestIdx+1,1,1, &mjd, &sta2);
			fits_write_col(fitsptrSpec2,TINT,    3,forestIdx+1,1,1, &fiberid, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 4,forestIdx+1,1,1, &XX, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 5,forestIdx+1,1,1, &YY, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 6,forestIdx+1,1,1, &ZZ, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 8,forestIdx+1,1,1, &meanForestLRF[0], &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 9,forestIdx+1,1,1, &alpha, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 11,forestIdx+1,1,nbPixel, &lambdaOBS, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 12,forestIdx+1,1,nbPixel, &norm_flux,       &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 13,forestIdx+1,1,nbPixel, &norm_flux_ivar,   &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 15,forestIdx+1,1,nbPixel, &continuum,   &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 16,forestIdx+1,1,nbPixel, &delta,       &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 17,forestIdx+1,1,nbPixel, &delta_weight,   &sta2);

			forestIdx ++;
		}
		fits_close_file(fitsptrSpec,&sta);
		std::cout << "  nb forest      = " << forestIdx      << std::endl;
	}

	fits_close_file(fitsptrSpec2,&sta2);
	std::cout << "  nb forest      = " << forestIdx      << std::endl;
	std::cout << "  nb_QSO_same_XY = " << nb_QSO_same_XY << std::endl;

}
void GetDelta::GetPDF(unsigned int version) {

	/// index of forest
	unsigned int forestIdx = 0;

	/// Get the file
	std::cout << "  file = " << pathToDataQSO__ << std::endl;
	const TString TSfitsnameSpec = pathToDataQSO__;
	int sta = 0;
	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);

	/// Get the number of spectra
	int tmp_nbSpectra = 0;
	unsigned int nbSpectra = 0;
	if (isTest__) nbSpectra = 10000;
	else {
		fits_get_num_hdus(fitsptrSpec, &tmp_nbSpectra, &sta);
		nbSpectra = tmp_nbSpectra-1;
	}
	std::cout << "  nbSpectra = " << nbSpectra << std::endl;
	
	/// Set to the first HDU
	fits_movabs_hdu(fitsptrSpec, 2,  NULL, &sta);
		
	/// Get the data
	for (unsigned int f=0; f<nbSpectra; f++) {
	
		/// Move to next HDU
		if (f!=0) fits_movrel_hdu(fitsptrSpec, 1,  NULL, &sta);
	
		/// Get the number of pixels in this HDU
		long tmp_nbPixels = 0;
		fits_get_num_rows(fitsptrSpec, &tmp_nbPixels, &sta);
		const unsigned int tmp_nbPixels2 = tmp_nbPixels;
		if (!findPDF__ && tmp_nbPixels2<C_MIN_NB_PIXEL) continue;
	
		/// Variable from old FITS
		float Z = 0.;
		float LAMBDA_OBS[tmp_nbPixels2];
		float FLUX[tmp_nbPixels2];
		float CONT[tmp_nbPixels2];
		fits_read_key(fitsptrSpec,TFLOAT,"ZQSO",&Z,NULL,&sta);
		fits_read_col(fitsptrSpec,TFLOAT, 1,1,1,tmp_nbPixels2,NULL, &LAMBDA_OBS,NULL,&sta);
		if (version==0 || version==1) {
			fits_read_col(fitsptrSpec,TFLOAT, 5,1,1,tmp_nbPixels2,NULL, &FLUX,      NULL,&sta);
			fits_read_col(fitsptrSpec,TFLOAT, 4,1,1,tmp_nbPixels2,NULL, &CONT,      NULL,&sta);
		}
		else {
			fits_read_col(fitsptrSpec,TFLOAT, 2,1,1,tmp_nbPixels2,NULL, &FLUX,      NULL,&sta);
		}

		/// Variables for new FITS
		const double oneOverOnePlusZ = 1./(1.+Z);

		for (unsigned int p=0; p<tmp_nbPixels2; p++) {

			/// bad pixels
			if (LAMBDA_OBS[p]<=0.) continue;
			//if (LAMBDA_OBS[p]<lambdaObsMin__ || LAMBDA_OBS[p]>=lambdaObsMax__) continue;

			const double lambdaRFd = LAMBDA_OBS[p]*oneOverOnePlusZ;
			if (lambdaRFd<lambdaRFMin__ || lambdaRFd>=lambdaRFMax__) continue;

			
			if (version==0) {
				if (CONT[p]==0.) continue;
				FLUX[p] /= CONT[p];
			}

			/// Get the PDF
			hFluxPDF__->Fill(FLUX[p],LAMBDA_OBS[p]/lambdaRFLine__-1.);
			hRedshift__->Fill(LAMBDA_OBS[p]/lambdaRFLine__-1.);
			hFlux__->Fill(FLUX[p]);
			hFluxVsLambdaObs__->Fill(LAMBDA_OBS[p], FLUX[p]);
			hFluxVsRedshift__->Fill(LAMBDA_OBS[p]/lambdaRFLine__-1., FLUX[p]);
			hFluxPow2VsLambdaObs__->Fill(LAMBDA_OBS[p], FLUX[p]*FLUX[p]);
                        hFluxPow2VsRedshift__->Fill(LAMBDA_OBS[p]/lambdaRFLine__-1., FLUX[p]*FLUX[p]);

		}
		forestIdx ++;
	}

	fits_close_file(fitsptrSpec,&sta);
	std::cout << "  nb forest =  " << forestIdx << std::endl;
}
void GetDelta::GetQSO(void) {
	const TString TSfitsnameSpec = pathToDataQSO__;

	/// Constants
	const double oneOverc_speedOfLight = 1./(C_C_LIGHT/1000);

	/// Fits file where to save
	TString TSfitsnameSpec2 = pathToSave__;
	TSfitsnameSpec2 += "QSO";
	if (withRSD__) {
		TSfitsnameSpec2 += "_withRSD";
	}
	else {
		TSfitsnameSpec2 += "_noRSD";
	}
	TSfitsnameSpec2 += ".fits";
	std::cout << "  " << TSfitsnameSpec2 << std::endl;

	/// Variables for FITS
	int sta    = 0;

	/// Get the file where there is data
	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);

	/// Get the number of spectra
	int tmp_nbSpectra = 0;
	unsigned int nbSpectra = 0;
	if (isTest__) nbSpectra = 10000;
	else {
		fits_get_num_hdus(fitsptrSpec, &tmp_nbSpectra, &sta);
		nbSpectra = tmp_nbSpectra-1;
	}
	std::cout << "  nbSpectra = " << nbSpectra << std::endl;


	/// Fits file where to save
	int sta2    = 0;
	fitsfile* fitsptrSpec2;
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READWRITE, &sta2); //READWRITE  //READONLY

	/// Number of QSO
	unsigned int nbQSOs = 0;
	unsigned int nb_QSO_same_XY = 0;
	float X_before = 0.;
	float Y_before = 0.;
	/// mean of velocity
	double velocity_mean[3] = {0.};

	/// Set to the first HDU
	fits_movabs_hdu(fitsptrSpec, 2,  NULL, &sta);

	/// Get the parameters of the QSO
	for (unsigned int i=0; i<nbSpectra; i++) {

		/// Get the data
		if (i!=0) fits_movrel_hdu(fitsptrSpec, 1,  NULL, &sta);

		float X = 0.;
		float Y = 0.;
		float Z = 0.;
		float V = 0.;

		fits_read_key(fitsptrSpec,TFLOAT,"X",   &X,NULL,&sta);
		fits_read_key(fitsptrSpec,TFLOAT,"Y",   &Y,NULL,&sta);
		fits_read_key(fitsptrSpec,TFLOAT,"ZQSO",&Z,NULL,&sta);
		fits_read_key(fitsptrSpec,TFLOAT,"VELOCITY",&V,NULL,&sta);

		double XX = X;
                double YY = Y;
                double ZZ = Z;

		if (X==X_before && Y_before==Y) {
			nb_QSO_same_XY ++;
		}
		X_before = X;
		Y_before = Y;

		/// Get mean velocity
		velocity_mean[0] += V*oneOverc_speedOfLight;
		velocity_mean[1] += V*oneOverc_speedOfLight*(1.+ZZ);
		velocity_mean[2] ++;

		/// If with RSD
		if (withRSD__) ZZ += V*oneOverc_speedOfLight*(1.+ZZ);

		/// Save the data
		fits_write_col(fitsptrSpec2,TDOUBLE, 1,nbQSOs+1,1,1, &XX, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 2,nbQSOs+1,1,1, &YY, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 3,nbQSOs+1,1,1, &ZZ, &sta2);

		nbQSOs ++;
	}

	/// Close files
	fits_close_file(fitsptrSpec,&sta);
	fits_close_file(fitsptrSpec2,&sta2);

	std::cout << "  nbQSOs         = "  << nbQSOs << std::endl;
	std::cout << "  nb_QSO_same_XY = "  << nb_QSO_same_XY << std::endl;

	std::cout << "  < v/c >        = "  << velocity_mean[0] / velocity_mean[2] << std::endl;
	std::cout << "  < delta z >    = "  << velocity_mean[1] / velocity_mean[2] << std::endl;
}







void GetDelta::defineHistos() {

	hFluxVsLambdaObs__ = new TProfile("hFluxVsLambdaObs__","",nbBinlambdaObs__+200,lambdaObsMin__-100.,lambdaObsMax__+100.);
	hFluxVsRedshift__  = new TProfile("hFluxVsRedshift__","",nbBinsRedshift__,minRedshift__,maxRedshift__);

	hFluxPow2VsLambdaObs__ = new TProfile("hFluxPow2VsLambdaObs__","",nbBinlambdaObs__+200,lambdaObsMin__-100.,lambdaObsMax__+100.);
        hFluxPow2VsRedshift__  = new TProfile("hFluxPow2VsRedshift__","",nbBinsRedshift__,minRedshift__,maxRedshift__);

	hFluxPDF__         = new TH2D("hFluxPDF__",            "",nbBinsFlux__,minFlux__,maxFlux__,nbBinsRedshift__,minRedshift__,maxRedshift__);
	hRedshift__        = new TH1D("hRedshift__",           "",nbBinsRedshift__,minRedshift__,maxRedshift__);
	hFlux__            = new TH1D("hFlux__",               "",nbBinsFlux__,minFlux__,maxFlux__);
}
void GetDelta::saveHistos() {

	/// Get the df
	const double df = 1.*(maxFlux__-minFlux__)/nbBinsFlux__;

	/// Normalize the PDF
	for (unsigned int j=0; j<nbBinsRedshift__; j++) {

		const double integral = hRedshift__->GetBinContent(j+1);
		if (integral==0.) continue; 
		const double invIntegral = 1./(integral*df);

		for (unsigned int i=0; i<nbBinsFlux__; i++){
			hFluxPDF__->SetBinContent(i+1,j+1, hFluxPDF__->GetBinContent(i+1,j+1)*invIntegral );
		}
	}
	
        TH1D* h_sigma_LSS_VsRedshift  = new TH1D("h_sigma_LSS_VsRedshift","",nbBinsRedshift__,minRedshift__,maxRedshift__);
	for (unsigned int i=0; i<nbBinsRedshift__; i++) {
		const double val1  = hFluxVsRedshift__->GetBinContent(i+1);
		if (val1==0.) continue;
		const double value = hFluxPow2VsRedshift__->GetBinContent(i+1)/(val1*val1) -1;
		h_sigma_LSS_VsRedshift->SetBinContent( i+1, value);
		h_sigma_LSS_VsRedshift->SetBinError(i+1,0.000000001);
	}
	
	/// Plots histos
	R_plot1D(hRedshift__,"z", "#");
	R_plot1D(hFlux__,"flux", "#");
	R_plot2D(hFluxPDF__,"flux noNosie noCont", "z","pdf");
	R_plot1D(hFluxVsLambdaObs__,"lambda_{Obs.}", "flux");
	R_plot1D(hFluxVsRedshift__,"z", "flux");
	R_plot1D(h_sigma_LSS_VsRedshift,"z", "sigma_2_LSS");

/*
	/// Save PDF
	std::ofstream fFile;
	fFile.open("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/test_FluxPDF_mocksJMC.txt");
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	for (unsigned int i=0; i<nbBinsFlux__; i++){
		for (unsigned int j=0; j<nbBinsRedshift__; j++){
			fFile << hFluxPDF__->GetBinContent(i+1,j+1) << " ";
		}
	}
	fFile.close();


	/// Save mean transmission flux
	fFile.open("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/test_hDeltaVsLambdaObs_LYA_JMC.txt");
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	for (unsigned int i=0; i<nbBinlambdaObs__; i++){
		fFile << i << " " << hFluxVsLambdaObs__->GetBinCenter(i+1) << " " << hFluxVsLambdaObs__->GetBinContent(i+1) << std::endl;
	}
	fFile.close();
*/

	std::ofstream fFile;
	fFile.open("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/v1575_h_sigma_LSS_VsRedshift_LYA_JMC.txt");
        fFile << std::scientific;
        fFile.precision(std::numeric_limits<double>::digits10);
	for (unsigned int i=0; i<nbBinsRedshift__; i++) {
		const double center = h_sigma_LSS_VsRedshift->GetBinCenter( i+1 );
                const double value  = h_sigma_LSS_VsRedshift->GetBinContent( i+1 );
		if (value==0.) continue;
                fFile << i << " " << center << " " << value << " " << std::endl;
        }
	fFile.close();


	/// Get Nicolas' PDF
	TH2D* hFluxPDF2 = new TH2D("hFluxPDF2","",100,0.0,1.0,50,1.96,3.5);
	std::string pathToPDF = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Resources/PDF/FluxPDF_mocksAuto_from_Nicolas.txt";
	ifstream filePDF(pathToPDF.c_str());
	///
	TH1D* hFluxVsLambdaObs2 = new TProfile("hFluxVsLambdaObs2","",nbBinlambdaObs__+200,lambdaObsMin__-100.,lambdaObsMax__+100.);
	TH1D* hRedshift2        = new TH1D("hRedshift2",           "",50,1.96,3.5);

	for (unsigned int j=0; j<46; j++){
		for (unsigned int i=0; i<100; i++){
			double pdf = 0.;
			filePDF>>pdf;
			hFluxPDF2->SetBinContent(i+1,j+1,pdf);
		}
	}
	// last rows
	for (unsigned int j=46; j<50; j++){
		for (unsigned int i=0; i<100; i++){
			const double pdf = hFluxPDF2->GetBinContent(i+1,45+1);
			hFluxPDF2->SetBinContent(i+1,j+1,pdf);
		}
	}
	R_plot2D(hFluxPDF2,"flux noNosie noCont", "z","pdf");

	/// Fill 1D histo
	for (unsigned int i=0; i<50; i++) {
		double value = 0.;
		for (unsigned int j=0; j<100; j++) {
			const double T     = hFluxPDF2->GetXaxis()->GetBinCenter(j+1);
			const double proba = hFluxPDF2->GetBinContent(j+1,i+1);
			value += T*proba;
		}
		value /= 100.;
		hRedshift2->SetBinContent(i+1, value);
		hRedshift2->SetBinError(i+1, value/1000.);
	}

	R_plot1D(hFluxVsLambdaObs2,"lambda_{Obs.}", "flux Nicolas");
	R_plot1D(hRedshift2,"z", "flux Nicolas");

/*
	/// Save mean transmission flux
	std::ofstream fFile;
	fFile.open("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/mocksAutoNicolas_hDeltaVsLambdaObs.txt");
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);
	for (unsigned int i=0; i<50; i++){
		fFile << i << " " << (hRedshift2->GetBinCenter(i+1)+1.)*lambdaRFLine__ << " " << hRedshift2->GetBinContent(i+1) << std::endl;
	}
	fFile.close();
*/
}



























