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

#include <fstream>
#include <iostream>	// std::cout
#include <sstream>	//stringstream
#include <stdlib.h>	//atoi
#include <unistd.h>	//getopt
#include <cstdlib>	// std::rand, std::srand

/// ROOT
#include "fitsio.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TProfile.h"

const double findPDF__ = true; 
bool hasCont = true;

GetDelta::GetDelta(int argc, char** argv) {

	std::cout << std::scientific;
	std::cout.precision(15);
	if (findPDF__) defineHistos();

	/// 0: Box, 1: Simulation
	if (argc<3) return;
	const std::string box_idx  = argv[1];
	const std::string sim_idx  = argv[2];


	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-78";
	pathToData__ += box_idx;
	pathToData__ += "-";
	pathToData__ += sim_idx;
	pathToData__ += ".fits";

	pathToSave__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/noNoisenoCont/Box_00";
	pathToSave__ += box_idx;
	pathToSave__ += "/Simu_00";
	pathToSave__ += sim_idx;
	pathToSave__ += "/Data/";


	std::cout << "\n"   << std::endl;
	std::cout << "  " << pathToData__ << std::endl;
	std::cout << "  " << pathToSave__ << std::endl;
	std::cout << "\n"   << std::endl;

	isTest__  = false;

	/// QSO
	withRSD__ = true;
	//GetQSO();

	/// Delta
	noCont__  = true;
	noNoise__ = true;
	//GetData();

	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-0.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-1.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-2.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-3.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-4.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-5.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-6.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-7.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-8.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-780-9.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-0.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-1.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-2.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-3.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-4.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-5.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-6.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-7.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-8.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/v1547/fits/spectra-781-9.fits";
	GetData();
	pathToData__ = "/home/gpfs/manip/mnt0607/bao/jmlg/QSOlyaMocks/spectra-highz.fits";
	hasCont = false;
	GetData();

	if (findPDF__) saveHistos();
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
void GetDelta::defineHistos() {

	hFluxVsLambdaObs__ = new TProfile("hFluxVsLambdaObs__","",nbBinlambdaObs__,lambdaObsMin__,lambdaObsMax__);
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
	
/*
	/// Plots histos
	R_plot1D(hRedshift__,"z", "#");
	R_plot1D(hFlux__,"flux", "#");
	R_plot2D(hFluxPDF__,"flux noNosie noCont", "z","pdf");
	R_plot1D(hFluxVsLambdaObs__,"lambda_{Obs.}", "flux");
*/

	/// Save PDF
	std::ofstream fFile;
	fFile.open("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/FluxPDF_mocksJMC.txt");
	fFile << std::scientific;
	fFile.precision(17);
	for (unsigned int i=0; i<nbBinsFlux__; i++){
		for (unsigned int j=0; j<nbBinsRedshift__; j++){
			fFile << hFluxPDF__->GetBinContent(i+1,j+1) << " ";
		}
	}
	fFile.close();


	/// Save mean transmission flux
	fFile.open("/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/hDeltaVsLambdaObs_LYA_JMC.txt");
	fFile << std::scientific;
	fFile.precision(17);
	for (unsigned int i=0; i<nbBinlambdaObs__; i++){
		fFile << i << " " << hFluxVsLambdaObs__->GetBinCenter(i+1) << " " << hFluxVsLambdaObs__->GetBinContent(i+1) << std::endl;
	}
	fFile.close();


/*	
	/// Get Nicolas' PDF
	TH2D* hFluxPDF2 = new TH2D("hFluxPDF2","",100,0.0,1.0,50,1.96,3.5);
	std::string pathToPDF = "/home/gpfs/manip/mnt0607/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/RootFile/FluxPDF.txt";
	ifstream filePDF(pathToPDF.c_str());

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
*/
	
}
void GetDelta::GetData(void) {


	/// Where are the data
	const TString TSfitsnameSpec = pathToData__;
	std::cout << "  pathToData__ :  " << pathToData__ << std::endl;

	/// Fits file where to save
	TString TSfitsnameSpec2 = pathToSave__;
	TSfitsnameSpec2 += "delta";
	if (noNoise__) {
		TSfitsnameSpec2 += "_noNoise";
	}
	if (noCont__) {
		TSfitsnameSpec2 += "_noCont";
	}
	TSfitsnameSpec2 += ".fits";
	std::cout << "  " << TSfitsnameSpec2 << std::endl;

	/// Variables for FITS
	int sta    = 0;

	/// Get the file
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
	std::cout << " nbSpectra = " << nbSpectra << std::endl;

	
	/// Get the number of pixels
	unsigned int nbMinPixels = 100000;
	unsigned int nbPixels    = 0;

	/// Set to the first HDU
	fits_movabs_hdu(fitsptrSpec, 2,  NULL, &sta);

	for (unsigned int i=0; i<nbSpectra; i++) {

		/// Move to next HDU
		if (i!=0) fits_movrel_hdu(fitsptrSpec, 1,  NULL, &sta);

		long tmp_nbPixels = 0;
		fits_get_num_rows(fitsptrSpec, &tmp_nbPixels, &sta);
		const unsigned int tmp_nbPixels2 = tmp_nbPixels;

		nbMinPixels = std::min(nbMinPixels,tmp_nbPixels2);
		nbPixels    = std::max(nbPixels,tmp_nbPixels2);
	}
	std::cout << " nbMinPixels  = " << nbMinPixels << std::endl;
	std::cout << " nbPixels     = " << nbPixels << std::endl;

	/// Fits file where to save
	int sta2    = 0;
	fitsfile* fitsptrSpec2;
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READONLY, &sta2); //READWRITE

	std::cout << "\n" << std::endl;



	/// Variable from fit
	double X = 0.;
	double Y = 0.;
	double Z = 0.;
	float LAMBDA_OBS[nbPixels];
	float FLUX[nbPixels];
	float FLUX_ERR[nbPixels];
	float CONTINUUM[nbPixels];
	float FLUXNON[nbPixels];
	for (unsigned int j=0; j<nbPixels; j++) {
		LAMBDA_OBS[j] = 0.;
		FLUX[j]       = 0.;
		FLUX_ERR[j]   = 0.;
		CONTINUUM[j]  = 0.;
		FLUXNON[j]    = 0.;
	}


	/// Variables for new fit
	double LAMBDA_OBS2[nbPixelsTemplate__];
	double LAMBDA_RF2[nbPixelsTemplate__];
	double FLUX2[nbPixelsTemplate__];
	double FLUX_ERR2[nbPixelsTemplate__];
	double DELTA[nbPixelsTemplate__];
	double DELTA_ERR[nbPixelsTemplate__];
	double TEMPLATE2[nbPixelsTemplate__];
	double NORM_FACTOR    = 0.;
	unsigned int NORM_FACTOR_nb = 0;
	double MEAN_FOREST_LAMBDA_RF = 0.;
	unsigned int NB_PIXEL = 0;
	for (unsigned int j=0; j<nbPixelsTemplate__; j++) {
		LAMBDA_OBS2[j] = 0.;
		LAMBDA_RF2[j]  = 0.;
		FLUX2[j]       = 0.;
		FLUX_ERR2[j]   = 0.;
		DELTA[j]       = 0.;
		DELTA_ERR[j]   = 1.;
		TEMPLATE2[j]   = 0.;
	}


	/// index of forest
	unsigned int forestIdx = 0;


	/// Set to the first HDU
	fits_movabs_hdu(fitsptrSpec, 2,  NULL, &sta);
	
	/// Get the data
	for (unsigned int i=0; i<nbSpectra; i++) {

		/// Move to next HDU
		if (i!=0) fits_movrel_hdu(fitsptrSpec, 1,  NULL, &sta);

		/// Get the number of pixels in this HDU
		long tmp_nbPixels = 0;
		fits_get_num_rows(fitsptrSpec, &tmp_nbPixels, &sta);
		const unsigned int tmp_nbPixels2 = tmp_nbPixels;
		if (tmp_nbPixels2 < nbBinRFMin__ && !findPDF__) continue;
		

		/// Variables for new fit
		NORM_FACTOR    = 0.;
		NORM_FACTOR_nb = 0;
		MEAN_FOREST_LAMBDA_RF = 0.;
		NB_PIXEL = 0;


		fits_read_key_dbl(fitsptrSpec,"X",&X,NULL,&sta);
		fits_read_key_dbl(fitsptrSpec,"Y",&Y,NULL,&sta);
		fits_read_key_dbl(fitsptrSpec,"ZQSO",&Z,NULL,&sta);
		fits_read_col(fitsptrSpec,TFLOAT, 1,1,1,tmp_nbPixels2,NULL, &LAMBDA_OBS,NULL,&sta);
		fits_read_col(fitsptrSpec,TFLOAT, 2,1,1,tmp_nbPixels2,NULL, &FLUX,      NULL,&sta);
		fits_read_col(fitsptrSpec,TFLOAT, 3,1,1,tmp_nbPixels2,NULL, &FLUX_ERR,  NULL,&sta);
		fits_read_col(fitsptrSpec,TFLOAT, 4,1,1,tmp_nbPixels2,NULL, &CONTINUUM, NULL,&sta);
		fits_read_col(fitsptrSpec,TFLOAT, 5,1,1,tmp_nbPixels2,NULL, &FLUXNON,   NULL,&sta);

		const double oneOverOnePlusZ = 1./(1.+Z);
		unsigned int nbPixelTemplate = 0;
		for (unsigned int j=0; j<tmp_nbPixels2; j++) {

			/// If no noise
			if (noNoise__) {
				FLUX[j]     = FLUXNON[j];
				FLUX_ERR[j] = 1.;
			}
			///  If no continuum
			if (hasCont && noCont__) FLUX[j] /= CONTINUUM[j];
			const double lambdaRF = LAMBDA_OBS[j]*oneOverOnePlusZ;

			/// Get the PDF
			if (lambdaRF>=lambdaRFMin__ && lambdaRF<lambdaRFMax__ && findPDF__) {
				hFluxPDF__->Fill(FLUX[j],LAMBDA_OBS[j]/lambdaRFLine__-1.);
				hRedshift__->Fill(LAMBDA_OBS[j]/lambdaRFLine__-1.);
				hFlux__->Fill(FLUX[j]);
				hFluxVsLambdaObs__->Fill(LAMBDA_OBS[j], FLUX[j]);
			}

			/// CCD and too many sky lines and bad pixels
			if (LAMBDA_OBS[j]<lambdaObsMin__ || LAMBDA_OBS[j]>=lambdaObsMax__ || FLUX_ERR[j]<=0.) continue;

			/// Normalisation zone
			if (lambdaRF>=lambdaRFNormaMin__ && lambdaRF<lambdaRFNormaMax__) {
				NORM_FACTOR    += FLUX[j];
				NORM_FACTOR_nb ++;
			}
			/// Keep only the template, i.e. the zone with the forest plus some pixels
			else if (lambdaRF>=lambdaRFTemplateMin__ && lambdaRF<lambdaRFTemplateMax__) {

				LAMBDA_OBS2[nbPixelTemplate] = LAMBDA_OBS[j];
				LAMBDA_RF2[nbPixelTemplate]  = lambdaRF;
				FLUX2[nbPixelTemplate]       = FLUX[j];
				FLUX_ERR2[nbPixelTemplate]   = FLUX_ERR[j];
				DELTA[nbPixelTemplate]       = FLUX[j];
				TEMPLATE2[nbPixelTemplate]   = CONTINUUM[j];
				nbPixelTemplate ++;

				if (lambdaRF>=lambdaRFMin__ && lambdaRF<lambdaRFMax__) {
					MEAN_FOREST_LAMBDA_RF += lambdaRF;
					NB_PIXEL ++;
				}
			}
		}

		if (nbPixelTemplate>nbPixelsTemplate__)  std::cout << nbPixelTemplate << " " << nbPixelsTemplate__ << std::endl;

		/// Normalisation zone
		if (NORM_FACTOR_nb == 0 || NORM_FACTOR <= 0. || nbPixelTemplate == 0 || NB_PIXEL < nbBinRFMin__) {
			//std::cout << Z << " " << NORM_FACTOR_nb << " " << NORM_FACTOR << " " << nbPixelTemplate << " " << NB_PIXEL << std::endl;
			continue;
		}

		/// Get the mean lambda_RF
		MEAN_FOREST_LAMBDA_RF /= NB_PIXEL;

		/// Normalize or not the forest
		if (noCont__) {
			NORM_FACTOR = 1.;
		}
		else {
			NORM_FACTOR /= NORM_FACTOR_nb;

			/// Normalise the data
			for (unsigned int j=0; j<nbPixelTemplate; j++) {
				FLUX2[j]     /= NORM_FACTOR;
				FLUX_ERR2[j]  = 1./((FLUX_ERR2[j]/NORM_FACTOR)*(FLUX_ERR2[j]/NORM_FACTOR));
				TEMPLATE2[j] /= NORM_FACTOR;

				DELTA[j] /= TEMPLATE2[j];
			}
		}
	
		/// Save data in second file
		fits_write_col(fitsptrSpec2,TDOUBLE, 4,forestIdx+1,1,1, &X, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 5,forestIdx+1,1,1, &Y, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 6,forestIdx+1,1,1, &Z, &sta2);
		fits_write_col(fitsptrSpec2,TINT,    7,forestIdx+1,1,1, &NB_PIXEL, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 8,forestIdx+1,1,1, &MEAN_FOREST_LAMBDA_RF, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 13,forestIdx+1,1,nbPixelTemplate, &LAMBDA_OBS2, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 14,forestIdx+1,1,nbPixelTemplate, &LAMBDA_RF2, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 15,forestIdx+1,1,nbPixelTemplate, &FLUX2,       &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 16,forestIdx+1,1,nbPixelTemplate, &FLUX_ERR2,   &sta2);

		/// Delta
		fits_write_col(fitsptrSpec2,TDOUBLE, 18,forestIdx+1,1,nbPixelTemplate, &DELTA,       &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 19,forestIdx+1,1,nbPixelTemplate, &DELTA_ERR,   &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 20,forestIdx+1,1,nbPixelTemplate, &DELTA_ERR,   &sta2);

		fits_write_col(fitsptrSpec2,TDOUBLE, 21,forestIdx+1,1,nbPixelTemplate, &TEMPLATE2,   &sta2);
		forestIdx ++;
	}

	fits_close_file(fitsptrSpec,&sta);
	fits_close_file(fitsptrSpec2,&sta2);

}
void GetDelta::GetQSO(void) {
	const TString TSfitsnameSpec = pathToData__;

	/// Constants
	const double oneOverc_speedOfLight = 1./c_speedOfLight__;

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
	fits_get_num_hdus(fitsptrSpec, &tmp_nbSpectra, &sta);
	const unsigned int nbSpectra = tmp_nbSpectra-1;
	std::cout << "  nbSpectra = " << nbSpectra << std::endl;

	/// Fits file where to save
	int sta2    = 0;
	fitsfile* fitsptrSpec2;
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READWRITE, &sta2); //READWRITE

	/// Number of QSO
	unsigned int nbQSOs = 0;

	/// Set to the first HDU
	fits_movabs_hdu(fitsptrSpec, 2,  NULL, &sta);

	/// Get the parameters of the QSO
	for (unsigned int i=0; i<nbSpectra; i++) {

		/// Get the data
		if (i!=0) fits_movrel_hdu(fitsptrSpec, 1,  NULL, &sta);

		double X = 0.;
		double Y = 0.;
		double Z = 0.;
		double V = 0.;

		fits_read_key_dbl(fitsptrSpec,"X",   &X,NULL,&sta);
		fits_read_key_dbl(fitsptrSpec,"Y",   &Y,NULL,&sta);
		fits_read_key_dbl(fitsptrSpec,"ZQSO",&Z,NULL,&sta);
		fits_read_key_dbl(fitsptrSpec,"VELOCITY",&V,NULL,&sta);

		/// If with RSD
		if (withRSD__) Z += V*oneOverc_speedOfLight*(1.+Z);

		/// Save the data
		fits_write_col(fitsptrSpec2,TDOUBLE, 1,nbQSOs+1,1,1, &X, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 2,nbQSOs+1,1,1, &Y, &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 3,nbQSOs+1,1,1, &Z, &sta2);

		nbQSOs ++;
	}

	/// Close files
	fits_close_file(fitsptrSpec,&sta);
	fits_close_file(fitsptrSpec2,&sta2);

	std::cout << "  nbQSOs = "  << nbQSOs << std::endl;
}












