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

#include C_MYROOTLIB
#include C_MYCPPLIB

#include <fstream>
#include <iostream> // std::cout
#include <sstream>	//stringstream
#include <stdlib.h>     /* atoi */
#include <unistd.h> //getopt


/// ROOT
#include "fitsio.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TFile.h"


extern void Chi2(int &npar, double *gin, double &f, double *par, int iflag);
TMinuit *mygMinuit;
unsigned int nbPix;
double lambdaObs[4844];
double lambdaRF[4844];
double fluxModel[4844];
double flux[4844];
double fluxIvar[4844];


GetDelta::GetDelta(int argc, char** argv)
{
	get_flux_vs_lambdaObs();
}
GetDelta::~GetDelta(void)
{	
}

void GetDelta::get_flux_vs_lambdaObs(void)
{

	/// Constants 
	const unsigned int nbBinRFMax = 4844;
	const double lambdaRFNormaMax = 1285.;
	//const double lambdaRFMin      = 1040.;
	//const double lambdaRFMax      = 1200.;
	const unsigned int nbPixelTemplate = 20000;
	const unsigned int nbBinlambdaObs  = 20000;
	//const double lambdaObsMin = 3600.;
	//const double lambdaObsMax = 7235.;
	const unsigned int nbLoop = 20;
	const unsigned int nbFile = 54; //54

	const std::string pathToFolder = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_all/";
	std::vector< TString> path(nbFile,"");
	unsigned int nbSpectra = 0;
	for (unsigned int f=0; f<path.size(); f++) {
		path[f] += pathToFolder;
		if (f!=53) {
			path[f] += "DR12_Guy_";
			path[f] += f*5000;
			path[f] += "_";
			path[f] += (f+1)*5000;
			path[f] += ".fits";
		}
		else path[f] += "DR12_Guy_265000_267722.fits";

		std::cout << path[f] << std::endl;

		/// Variables for FITS
		int sta = 0;
		long nrows = 0;

		const TString TSfitsnameSpec = path[f];
		fitsfile* fitsptrSpec;
		fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
		fits_get_num_rows(fitsptrSpec, &nrows, &sta);
		fits_close_file(fitsptrSpec,&sta);
		nbSpectra += nrows;
	}

	/// vector to put the result of fit
	std::vector<double> fitCoef(nbSpectra,1.);
	std::cout << "  nbSpectra = " << nbSpectra << std::endl;

	/// Init the fit
	initFitCont();

	///
	const TString tPathToSave = "histos.root";
	TFile* storeFile = new TFile(tPathToSave,"RECREATE","");
	storeFile->cd();
	///
	TH1D* hDeltaVsLambdaRF = new TH1D("hDeltaVsLambdaRF","",nbPixelTemplate, 0., 20000.);
	R_dealWithPlots_1D(hDeltaVsLambdaRF, "#lambda_{R.F.} (A)", "<#delta+1>", "hDeltaVsLambdaRF");
	for (unsigned int i=0; i<nbPixelTemplate; i++) {
		hDeltaVsLambdaRF->SetBinContent(i+1,1.);
		hDeltaVsLambdaRF->SetBinError(i+1,0.);
	}
	///
	TH1D* hDeltaVsLambdaRFResiduals[nbLoop];
	for (unsigned int i=0; i<nbLoop; i++) {
		TString name = "hDeltaVsLambdaRFResiduals_";
		name += i;
		hDeltaVsLambdaRFResiduals[i] = new TH1D(name,"",nbPixelTemplate, 0., 20000.);
		R_dealWithPlots_1D(hDeltaVsLambdaRFResiduals[i], "#lambda_{R.F.} (A)", "<#delta+1>", "hDeltaVsLambdaRFResiduals");
	}
	///
	TProfile* hDeltaVsLambdaObsInRed[nbLoop];
	for (unsigned int i=0; i<nbLoop; i++) {
		TString name = "hDeltaVsLambdaObsInRed_";
		name += i;
		hDeltaVsLambdaObsInRed[i] = new TProfile(name,"",nbBinlambdaObs, 0., 20000.);
		R_dealWithPlots_1D(hDeltaVsLambdaObsInRed[i], "#lambda_{Obs.} (A)", "", "hDeltaVsLambdaObsInRed");
	}
	///
	TH1D* hParamFit = new TH1D("hParamFit","",1200, -6., 6.);
	R_dealWithPlots_1D(hParamFit, "", "", "hParamFit");




	for (unsigned int loop=0; loop<nbLoop; loop++) {

		std::cout << "  loop = " << loop << std::endl;
		unsigned int idxSpectrum = 0;

		/// Delta vs. lambda_Obs
		double deltaVSLambdaRF[nbPixelTemplate][2];
		for (unsigned int i=0; i<nbPixelTemplate; i++) {
			deltaVSLambdaRF[i][0] = 0.;
			deltaVSLambdaRF[i][1] = 0.;
		}


		for (unsigned int f=0; f<path.size(); f++) {

			/// Variables for FITS
			int sta = 0;
			long nrows = 0;

			const TString TSfitsnameSpec = path[f];
			fitsfile* fitsptrSpec;
			fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
			fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
			/// Load data
			for (unsigned int i=0; i<nrows; i++) {  //nrows
		
				double LAMBDA_OBS[nbBinRFMax];
				double LAMBDA_RF[nbBinRFMax];
				double FLUX[nbBinRFMax];
				double FLUX_IVAR[nbBinRFMax];
				nbPix = 0;

				fits_read_col(fitsptrSpec,TDOUBLE, 8,  i+1,1,nbBinRFMax,NULL,&LAMBDA_OBS,NULL,&sta);
				fits_read_col(fitsptrSpec,TDOUBLE, 9,  i+1,1,nbBinRFMax,NULL,&LAMBDA_RF,NULL,&sta);
				fits_read_col(fitsptrSpec,TDOUBLE, 10, i+1,1,nbBinRFMax,NULL,&FLUX,NULL,&sta);
				fits_read_col(fitsptrSpec,TDOUBLE, 11, i+1,1,nbBinRFMax,NULL,&FLUX_IVAR,NULL,&sta);

				/// Put values in the array
				for (unsigned int p=0; p<nbBinRFMax; p++) {
					if (FLUX_IVAR[p]<=0. || LAMBDA_RF[p]<lambdaRFNormaMax) continue;
					else {
						lambdaObs[nbPix]  = LAMBDA_OBS[p];
						lambdaRF[nbPix]   = LAMBDA_RF[p];
						fluxModel[nbPix]  = hDeltaVsLambdaRF->Interpolate(LAMBDA_RF[p]);
						flux[nbPix]       = FLUX[p];
						fluxIvar[nbPix]   = FLUX_IVAR[p];
						nbPix++;
					}
				}

				/// Fit
				int iflag;
				double arglist[10];
				mygMinuit->mnparm(0, "alpha", 1., fitCoef[idxSpectrum], -1000.0, 1000.0, iflag);
				// Minimization
				arglist[0] = 1000.;
				arglist[1] = 0.1;
				mygMinuit->mnexcm("MIGRAD", arglist,2, iflag);
				double err;
				double param;
				mygMinuit->GetParameter(0,param,err);

				/// Save parameter
				fitCoef[idxSpectrum] = param;
				if (param>0.) hParamFit->Fill(log10(param));

				/// Get residuals
				for (unsigned int i=0; i<nbPix; i++) {
					const double coef = param*fluxModel[i];
					const unsigned int idx = lambdaRF[i];
					deltaVSLambdaRF[idx][0] += fluxIvar[i]*coef*flux[i];
					deltaVSLambdaRF[idx][1] += fluxIvar[i]*coef*coef;
					hDeltaVsLambdaObsInRed[loop]->Fill( lambdaObs[i], flux[i]/coef, fluxIvar[i]*coef*coef);
				}

				idxSpectrum++;
			}
			fits_close_file(fitsptrSpec,&sta);
		}





		/// Get the new template
		for (unsigned int i=0; i<nbPixelTemplate; i++) {

			if (deltaVSLambdaRF[i][1]==0.) {
				hDeltaVsLambdaRF->SetBinContent(i+1,0.);
				hDeltaVsLambdaRF->SetBinError(i+1,0.);
			}
			else {
				deltaVSLambdaRF[i][0] /= deltaVSLambdaRF[i][1];
				const double value = hDeltaVsLambdaRF->GetBinContent(i+1)*deltaVSLambdaRF[i][0];
				hDeltaVsLambdaRF->SetBinContent(i+1,value);
				hDeltaVsLambdaRF->SetBinError(i+1,0.0001);
				hDeltaVsLambdaRFResiduals[loop]->SetBinContent(i+1,deltaVSLambdaRF[i][0]);
				hDeltaVsLambdaRFResiduals[loop]->SetBinError(i+1,0.0001);
			}
		}




		// Method2: Set "hMeanFlux"
		double xxx1_value_m2 = 0.;
		double xxx2_value_m2 = 0.;
		double yyy1_value_m2 = 0.;
		double yyy2_value_m2 = 0.;
		double value_m2 = 0.;
		unsigned int nbEmptyPixels_m2 = 0;
		for (unsigned int i=0; i<nbPixelTemplate; i++) {
		
			// If empty and not the last pixel
			if (hDeltaVsLambdaRF->GetBinError(i+1) == 0. && i!=nbPixelTemplate-1) {
				if (nbEmptyPixels_m2 == 0 && i!=0) {
					xxx1_value_m2 = hDeltaVsLambdaRF->GetBinCenter(i);
					yyy1_value_m2 = hDeltaVsLambdaRF->GetBinContent(i);
					value_m2 = hDeltaVsLambdaRF->GetBinContent(i);
				}
				nbEmptyPixels_m2 ++;
			}
			else {
			
				// If first pixel not empty after some empty pixels
				if (nbEmptyPixels_m2!=0) {

					double a = 0.;
					double b = 0.;					
				
					// Find the mean value between the two not empty pixels
					if (value_m2!=0. && hDeltaVsLambdaRF->GetBinError(i+1)!=0.) {
						// Linear extrapolation
						xxx2_value_m2 = hDeltaVsLambdaRF->GetBinCenter(i+1);
						yyy2_value_m2 = hDeltaVsLambdaRF->GetBinContent(i+1);
						a = yyy1_value_m2 - yyy2_value_m2;
						b = yyy2_value_m2*xxx1_value_m2 - yyy1_value_m2*xxx2_value_m2;
						if (xxx1_value_m2-xxx2_value_m2 != 0.) {
							a /= (xxx1_value_m2-xxx2_value_m2);
							b /= (xxx1_value_m2-xxx2_value_m2);
						}
						value_m2 = (value_m2 + hDeltaVsLambdaRF->GetBinContent(i+1))/2.;
					}
					
					// Set all the empty pixels 
					for (unsigned int j=0; j<nbEmptyPixels_m2; j++) {
						double tmp_value = a*hDeltaVsLambdaRF->GetBinCenter(i-j) + b;
						if (hDeltaVsLambdaRF->GetBinError(i+1)==0.) tmp_value = value_m2;
						hDeltaVsLambdaRF->SetBinContent(i-j, tmp_value); //value_m2
						//hDeltaVsLambdaRF->SetBinError(i-j, maxError_m2);
					}
				
					// If the last pixel is also empty
					if (i==nbPixelTemplate-1) {
						hDeltaVsLambdaRF->SetBinContent(i+1, value_m2);
						//hDeltaVsLambdaRF->SetBinError(i+1, maxError_m2);
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
	}

	storeFile->Write();
	storeFile->Close();
}
void write_histo_to_ascii(void) {

	

}
extern void Chi2(int &npar, double *gin, double &f, double *par, int iflag) {

	f = 0.;
	for (unsigned int i=0; i<nbPix; i++) {
		const double tmp = par[0]*fluxModel[i]-flux[i];
		f += tmp*tmp*fluxIvar[i];
	}
}
void GetDelta::initFitCont(void) {

	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "\n Initialisation of the continuum fit with minuit (method 2) " << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	mygMinuit = new TMinuit(1);
	mygMinuit->Clear();
	mygMinuit->SetFCN(Chi2);

	int iflag;
	double arglist[1];
	arglist[0] = 0.;
	mygMinuit->mnexcm("SET ERR", arglist, 1, iflag);

	// Set starting values and step sizes for parameters

	mygMinuit->mnparm(0, "a", 1.0, 0.1, -1000.0, 1000.0, iflag);

	mygMinuit->SetPrintLevel(-1);
}






















