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
#include "../../../Constants/globalValues.h"

#include <cmath>
#include <fstream>
#include <iostream>	// std::cout
#include <sstream>	//stringstream
#include <stdlib.h>	// atoi
#include <unistd.h>	//getopt


//// ROOT
#include "fitsio.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMath.h"


//// Constants
std::string pathToTxt__    = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/chain_annalys_delta/";
std::string pathToDLACat__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DLA_all.fits";
const unsigned int nbPixelTemplate__ = int(lambdaRFMax__-lambdaRFMin__)+6;
const unsigned int nbBinlambdaObs__  = int(lambdaObsMax__-lambdaObsMin__);
const double onePlusZ0__ = 1.+z0__;
const double halfGama__  = gama__/2.;
const unsigned int nbBinPowErrorDelta__ = 20*int(log10(maxPowErrorDelta__) - log10(minPowErrorDelta__));



//// Global variables
extern void Chi2_method1(int &npar, double *gin, double &f, double *par, int iflag);
extern void Chi2_method2(int &npar, double *gin, double &f, double *par, int iflag);
double ProbPixel(double cont,double flux, double sig, unsigned int idxPixel);
TH2D* hFluxPDF__;
TMinuit *mygMinuit;
double LambdaMeth2[nbBinRFMax__];
double FluxMeth2[nbBinRFMax__];
double FluxErrMeth2[nbBinRFMax__];
double FluxMeanMeth2[nbBinRFMax__];
unsigned int NbPixMeth2;
double LambdaMeanMeth2;
double a_PDF[nbBinsFlux__][nbBinRFMax__];
std::string pathForest__ = "";
std::string pathMoreForHist__ = "";


//// Flags
unsigned int stepDefinition = 2;
unsigned int stepAnnalyse   = 0;
unsigned int methodIndex__ = 2
const bool doVetoLines__          = false;
const bool setDLA__               = false;
const bool cutNotFittedSpectra__  = true;
const bool putReobsTogether__     = false;
const bool mocksColab__           = false;
const bool mockJMC__              = false;

GetDelta::GetDelta(int argc, char** argv) {

	std::cout << std::scientific;
	std::cout.precision(15);

	if (!mocksColab__ && !mockJMC__) {
		pathForest__   = "/home/gpfs/manip/mnt/bao/hdumasde/Data/";
		pathForest__  += forest__;
//		pathForest__  += "/FitsFile_DR12_Guy/DR12_reObs/DR12_reObs.fits";
//		pathForest__  += "/FitsFile_DR12_Guy/DR12_primery/DR12_primery_test_PDFMocksJMC_meanLambda_testNoCap.fits";
//		pathForest__  += "/FitsFile_eBOSS_Guy/all_eBOSS_primery/eBOSS_primery.fits";
		pathForest__  += "/FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits";
//pathForest__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/chain_annalys_delta/Get_delta/src/DR12_primery_test_1000.fits";
	}
	else {
	
		const std::string box = argv[1];
		const std::string sim = argv[2];

		pathMoreForHist__  = "_";
		pathMoreForHist__ += box;
		pathMoreForHist__ += "_";
		pathMoreForHist__ += sim;

		if (!mockJMC__) {
			std::cout << "  Mocks colab : box = " << box << " ,  simulation = " << sim << std::endl;

			pathForest__   = "/home/gpfs/manip/mnt0607/bao/hdumasde/MockV4/M3_0_";
			pathForest__  += box;
			pathForest__  += "/00";
			pathForest__  += sim;
			pathForest__  += "/mock.fits";

			pathMoreForHist__ += "_MocksColab";
		}
		else if (mockJMC__) {
			std::cout << "  Mocks Jean-Marc : box = " << box << " ,  simulation = " << sim << std::endl;

			pathForest__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_00";
			pathForest__ += box;
			pathForest__ += "/Simu_00";
			pathForest__ += sim;
			pathForest__ += "/Data/delta.fits";
		}
	}


	std::cout << "\n\n" << std::endl;

	//// Set Vectors and histos
	defineHistos();

	//// Set the DLA values
	if (setDLA__) updateDLA(pathForest__,0,0);

	if (stepAnnalyse==0) {

		loadDataForest(pathForest__,0,0,cutNotFittedSpectra__);

		//// Find the template
		for (unsigned int i=0; i<nbLoop__; i++) {
			std::cout << "\n  --- Step = " << i << std::endl;
			getHisto(i);
			if (i!=nbLoop__-1) updateDeltaVector(i);
		}

		//// Save changes
		updateDelta(pathForest__,nbLoop__-1,0,0);

		// put re-obs together
		if (putReobsTogether__) putReobsTogether(pathForest__,nbLoop__-1);
	}
	
	//// Fit the forest
	if (stepAnnalyse==1) {
		std::string str = argv[3];
		unsigned int start = atoi( str.c_str() );
		str = argv[4];
		unsigned int end   = atoi( str.c_str() );

		std::cout << "  " << start << " " << end << std::endl;
		
		loadDataForest(pathForest__,start,end);
		if (v_zz__.size()>=0) fitForests(start,end);
	}


	std::cout << "  Finished " << std::endl;

}
GetDelta::~GetDelta(void) {
	
	if (stepAnnalyse==0) {
		for (unsigned int i=0; i<nbLoop__+1; i++) {
			delete hDeltaVsLambdaObs__[i];
			delete hDeltaVsLambdaRF__[i];

			if (i<nbLoop__) {
				delete hTemplate__[i];
				delete grEta__[i];
				delete grSig__[i];
			}
		}
		delete hFluxPDF__;
	}
	else if (stepAnnalyse==1) delete mygMinuit;

}


void GetDelta::defineHistos() {

	for (unsigned int i=0; i<nbLoop__; i++) {
		// Template
		TString name = "hTemplate_";
		name += i;

		hTemplate__[i] = new TH1D(name,"",nbPixelTemplate__,lambdaRFMin__-3.,lambdaRFMax__+3.);
	}
	for (unsigned int i=0; i<nbLoop__+1; i++) {

		// Mean transmission flux
		TString name = "hDeltaVsLambdaObs_";
		name += i;

		hDeltaVsLambdaObs__[i] = new TH1D(name,"",nbBinlambdaObs__,lambdaObsMin__,lambdaObsMax__);
		R_dealWithPlots_1D(hDeltaVsLambdaObs__[i], "#lambda_{Obs.} (A)", "Mean transmission flux", "Method2: mean transmission flux");
		for (unsigned int j=0; j<nbBinlambdaObs__; j++) {
			hDeltaVsLambdaObs__[i]->SetBinContent(j+1,1.);
			hDeltaVsLambdaObs__[i]->SetBinError(j+1,0.);
		}
	}
	for (unsigned int i=0; i<nbLoop__+1; i++) {
		//// <delta+1> vs. lambda_RF
		TString name = "hDeltaVsLambdaRF_";
		name += i;

		hDeltaVsLambdaRF__[i] = new TH1D(name,"",nbPixelTemplate__,lambdaRFMin__-3.,lambdaRFMax__+3.);
		R_dealWithPlots_1D(hDeltaVsLambdaRF__[i], "#lambda_{R.F.} (A)", "<#delta+1>", "");
		for (unsigned int j=0; j<nbPixelTemplate__; j++) {
			hDeltaVsLambdaRF__[i]->SetBinContent(j+1,1.);
			hDeltaVsLambdaRF__[i]->SetBinError(j+1,0.);
		}
	}


	if ( ((mocksColab__ || !mockJMC__) && stepDefinition >= 1) || (mockJMC__ && stepDefinition >= 2)  ) {

		//// Put the value of the init of mean transmision flux
		std::string path = pathToTxt__;
		path += "hDeltaVsLambdaObs_";
		path += forest__;
		path += pathMoreForHist__;
		path += ".txt";
		std::ifstream file(path.c_str());

		unsigned int idx;
		double val;
		double err;
		while (file) {
			file>>idx>>val>>err;
			if(file==0) break;

			hDeltaVsLambdaObs__[nbLoop__]->SetBinContent(idx+1,val);
			hDeltaVsLambdaObs__[nbLoop__]->SetBinError(idx+1,err);
		}
		file.close();
	}
	if ( ((mocksColab__ || !mockJMC__) && stepDefinition >= 1) || (mockJMC__ && stepDefinition >= 2)  ) {

		//// Put the value of the init of mean transmision flux
		std::string path = pathToTxt__;
		path += "hDeltaVsLambdaRF_";
		path += forest__;
		path += pathMoreForHist__;
		path += ".txt";
		std::ifstream file(path.c_str());

		unsigned int idx;
		double val;
		double err;
		while (file) {
			file>>idx>>val>>err;
			if(file==0) break;

			hDeltaVsLambdaRF__[nbLoop__]->SetBinContent(idx+1,val);
			hDeltaVsLambdaRF__[nbLoop__]->SetBinError(idx+1,err);
		}
		file.close();
	}


	// Flux PDF
	//#define PATHTOCODE "/home/gpfs/manip/mnt0607/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation";
	//std::string pathToPDF = "/home/gpfs/manip/mnt0607/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/RootFile/FluxPDF.txt";
//	hFluxPDF__ = new TH2D("hFluxPDF__","",100,0.0,1.0,50,1.96,3.5);

	hFluxPDF__ = new TH2D("hFluxPDF__","",nbBinsFlux__,minFlux__,maxFlux__,nbBinsRedshift__,minRedshift__,maxRedshift__);
	hFluxPDF__->SetXTitle("flux");
	hFluxPDF__->SetYTitle("z");
	hFluxPDF__->SetZTitle("pdf");
	hFluxPDF__->SetTitle("Flux PDF: method 2");
	std::string pathToPDF = pathToTxt__;
	pathToPDF += "FluxPDF_mocksJMC.txt";
	ifstream filePDF(pathToPDF.c_str());

/*
	for (unsigned int j=0; j<46; j++){
		for (unsigned int i=0; i<100; i++){
			double pdf = 0.;
			filePDF>>pdf;
			hFluxPDF__->SetBinContent(i+1,j+1,pdf);
		}
	}
	// last rows
	for (unsigned int j=46; j<50; j++){
		for (unsigned int i=0; i<100; i++){
			const double pdf = hFluxPDF__->GetBinContent(i+1,45+1);
			hFluxPDF__->SetBinContent(i+1,j+1,pdf);
		}
	}
*/
	for (unsigned int i=0; i<nbBinsFlux__; i++){
		for (unsigned int j=0; j<nbBinsRedshift__; j++){
			double pdf = 0.;
			filePDF>>pdf;
			hFluxPDF__->SetBinContent(i+1,j+1,pdf);
		}
	}

	if (stepAnnalyse==1) {
		initFitCont();
		mygMinuit->SetPrintLevel(-1);
	}

}
void GetDelta::getHisto(unsigned int loopIdx) {

	//// Define the histos
	
	//// Template (m2)
	double m2_template[nbPixelTemplate__][4];

	//// Delta vs. lambda_RF
	double deltaVSLambdaRF[nbPixelTemplate__][4];

	//// Delta vs. lambda_Obs
	double deltaVSLambdaObs[nbBinlambdaObs__][4];


	for (unsigned int i=0; i<nbPixelTemplate__; i++) {
		for (unsigned int j=0; j<4; j++) {
			m2_template[i][j]     = 0.;
			deltaVSLambdaRF[i][j] = 0.;
		}
	}
	for (unsigned int i=0; i<nbBinlambdaObs__; i++) {
		for (unsigned int j=0; j<4; j++) {
			deltaVSLambdaObs[i][j] = 0.;
		}
	}

	//// For method 2
	TH1D* hVarDeltaVsError[C_NBBINZ+2];
	for (unsigned int j=0; j<C_NBBINZ+2; j++) {	
		TString name = "hVarDeltaVsErrorMethod2_";
		name += j;
			
		hVarDeltaVsError[j] = new TH1D(name,"",nbBinPowErrorDelta__, log10(minPowErrorDelta__), log10(maxPowErrorDelta__));
		R_dealWithPlots_1D(hVarDeltaVsError[j], "#sigma^{-2}_{pipeline}", "Var^{-1}_{#delta}", "Method2: ");
		R_binLogX(hVarDeltaVsError[j]);
	}

	//// Bin center for eta and sigma2LSS
	double binCenterZZ[C_NBBINZ+2][4] = {};
	//// Array for weight
	double arrayWeight[C_NBBINZ+2][nbBinPowErrorDelta__][4];
	for (unsigned int i=0; i<C_NBBINZ+2; i++) {
		for (unsigned int ii=0; ii<nbBinPowErrorDelta__; ii++) {
			for (unsigned int iii=0; iii<4; iii++) arrayWeight[i][ii][iii] = 0.;
		}
	}

	// Fit function for the variance of delta
	TF1* fitWeight = new TF1("fitWeight", "1./(pow(x, -1.)/[0] + [1])", minPowErrorDeltaForFit__, maxPowErrorDeltaForFit__);
	fitWeight->SetParName(0, "#eta");
	fitWeight->SetParameter(0, etaStart);
	fitWeight->SetParLimits(0, 0., 1000.);
	fitWeight->SetParName(1, "#sigma^{2}_{LSS}");
	fitWeight->SetParameter(1, sigma2LSSStart);
	fitWeight->SetParLimits(1, 0., 1000.);

	const unsigned int nbForest = v_zz__.size();
	double meanDelta_Zero[3] = {0.};
	
	for (unsigned int f=0; f<nbForest; f++) {

		const unsigned int nb = v_nbPixel__[f];

		for (unsigned int p=0; p<nb; p++) {

			const double flux      = v_NORM_FLUX__[f][p];
			const double z         = v_ZZZ__[f][p];
			const double d         = v_DELTA__[f][p];
			const double deltaIvar = v_DELTA_IVAR__[f][p];

			const double w         = v_DELTA_WEIGHT__[f][p];
			if (w<=0.) {
				if (w<0) std::cout << "  GetDelta::getHisto::  ERROR::  w<0 , w = " << w << "  , forest = " << f << "  , pixel = " << p << std::endl;
				if (v_LAMBDA_RF__[f][p] >= lambdaRFMin__ && v_LAMBDA_RF__[f][p] < lambdaRFMax__) meanDelta_Zero[2] ++;
				continue;
			}
			const double wf        = w*flux;
			const double wd        = w*d;
			const double wdd       = wd*d;
			const double wz        = w*z;
		
			//// Template (m2)
			unsigned int idx = int( v_LAMBDA_RF__[f][p]-lambdaRFMin__+3.);
			m2_template[idx][0] += wf;
			m2_template[idx][1] += w;
			m2_template[idx][2] += wf*flux;
			m2_template[idx][3] ++;

			//// Delta vs. lambda_RF
			deltaVSLambdaRF[idx][0] += wd;
			deltaVSLambdaRF[idx][1] += w;
			deltaVSLambdaRF[idx][2] += wdd;
			deltaVSLambdaRF[idx][3] ++;

			// Only in the forest
			if (v_LAMBDA_RF__[f][p] >= lambdaRFMin__ && v_LAMBDA_RF__[f][p] < lambdaRFMax__) {

				//// Delta vs. lambda_Obs
				idx = int( v_LAMBDA_OBS__[f][p]-lambdaObsMin__);
				deltaVSLambdaObs[idx][0] += wd;
				deltaVSLambdaObs[idx][1] += w;
				deltaVSLambdaObs[idx][2] += wdd;
				deltaVSLambdaObs[idx][3] ++;

				//// Weight
			
				idx = 0;
				if (z>=C_ZEXTREMA0 && z<C_ZEXTREMA1) idx = int( (z-C_ZEXTREMA0)/C_BINSIZEZ )+1;
				else if ( z>=C_ZEXTREMA1) idx = C_NBBINZ+1;

				// Histos for the weight
				binCenterZZ[idx][0] += wz;
				binCenterZZ[idx][1] += w;
				binCenterZZ[idx][2] += wz*z;
				binCenterZZ[idx][3] ++;
	
				if ( deltaIvar >= minPowErrorDelta__ && deltaIvar < maxPowErrorDelta__ ) {
					const unsigned int idxWeight = hVarDeltaVsError[0]->GetXaxis()->FindBin(deltaIvar)-1;
					arrayWeight[idx][idxWeight][0] += wd;
					arrayWeight[idx][idxWeight][1] += wdd;
					arrayWeight[idx][idxWeight][2] += w;
					arrayWeight[idx][idxWeight][3] ++;
				}
			}
		}
	}


	//// Get the number of pixels and the mean delta
	double meanDelta[3] = {0.};
	for (unsigned int i=0; i<nbBinlambdaObs__; i++) {
		meanDelta[0] += deltaVSLambdaObs[i][0];
		meanDelta[1] += deltaVSLambdaObs[i][1];
		meanDelta[2] += deltaVSLambdaObs[i][3];
	}
	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;
	std::cout << "  nb pixel (w==0) = " << (long long unsigned int)meanDelta_Zero[2]         << std::endl;
	std::cout << "  all pixel       = " << (long long unsigned int)(meanDelta[2]+meanDelta_Zero[2]) << "\n" << std::endl;




	//// Save the histos

	std::ofstream fFile;
	std::ofstream fFile1;
	std::ofstream fFile2;
	std::string tmp_pathToSave;





	//// Template (m2)
	tmp_pathToSave  = pathToTxt__;
	tmp_pathToSave += "template_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += pathMoreForHist__;
	tmp_pathToSave += ".txt";
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(17);

	//// Template (m2)
	for (unsigned int i=0; i<nbPixelTemplate__; i++) {
		if (m2_template[i][1]!=0.) {
			fFile << lambdaRFMin__-3.+0.5 + 1.*i;
			fFile << " " << m2_template[i][0]/m2_template[i][1];
			fFile << " " << m2_template[i][0];
			fFile << " " << m2_template[i][1] << std::endl;

			const double mean = m2_template[i][0]/m2_template[i][1];
			const double err  = sqrt( (m2_template[i][2]/m2_template[i][1]-mean*mean)/m2_template[i][3] );
			hTemplate__[loopIdx]->SetBinContent(i+1,mean);
			hTemplate__[loopIdx]->SetBinError(i+1,err);
		}
	}
	fFile.close();


	//// Delta vs. lambda_RF
	tmp_pathToSave  = pathToTxt__;
	tmp_pathToSave += "deltaVSLambdaRF_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += pathMoreForHist__;
	tmp_pathToSave += ".txt";
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(17);

	unsigned int loopIdxForHist1 = loopIdx-1;
	if (loopIdx==0) loopIdxForHist1 = nbLoop__;

	//// Delta vs. lambda_RF
	for (unsigned int i=0; i<nbPixelTemplate__; i++) {
		if (deltaVSLambdaRF[i][1]!=0.) {
			fFile << lambdaRFMin__-3.+0.5 + 1.*i;
			fFile << " " << deltaVSLambdaRF[i][0]/deltaVSLambdaRF[i][1];
			fFile << " " << deltaVSLambdaRF[i][0];
			fFile << " " << deltaVSLambdaRF[i][1] << std::endl;

			double mean = deltaVSLambdaRF[i][0]/deltaVSLambdaRF[i][1];
			double err  = sqrt( (deltaVSLambdaRF[i][2]/deltaVSLambdaRF[i][1]-mean*mean)/deltaVSLambdaRF[i][3] );
			if (stepDefinition == 0) {
				mean = 0.;
				err  = 0.001;
			}

			hDeltaVsLambdaRF__[loopIdx]->SetBinContent(i+1,(mean+1.)*hDeltaVsLambdaRF__[loopIdxForHist1]->GetBinContent(i+1) );
			hDeltaVsLambdaRF__[loopIdx]->SetBinError(i+1,err);
		}
	}
	fFile.close();

	
	//// Save hDeltaVsLambdaRF
	tmp_pathToSave  = pathToTxt__;
	tmp_pathToSave += "hDeltaVsLambdaRF_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += pathMoreForHist__;
	tmp_pathToSave += ".txt";
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(17);

	for (unsigned int i=0; i<nbPixelTemplate__; i++) {
		fFile << i;
		fFile << " " << hDeltaVsLambdaRF__[loopIdx]->GetBinContent(i+1);
		fFile << " " << hDeltaVsLambdaRF__[loopIdx]->GetBinError(i+1) << std::endl;
	}
	fFile.close();



	//// Delta vs. lambda_Obs
	tmp_pathToSave  = pathToTxt__;
	tmp_pathToSave += "deltaVSLambdaObs_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += pathMoreForHist__;
	tmp_pathToSave += ".txt";
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(17);

	unsigned int loopIdxForHist = loopIdx-1;
	if (loopIdx==0) loopIdxForHist = nbLoop__;

	//// Delta vs. lambda_Obs
	for (unsigned int i=0; i<nbBinlambdaObs__; i++) {
		if (deltaVSLambdaObs[i][1]!=0.) {
			fFile << lambdaObsMin__+0.5 + i;
			fFile << " " << deltaVSLambdaObs[i][0]/deltaVSLambdaObs[i][1];
			fFile << " " << deltaVSLambdaObs[i][0];
			fFile << " " << deltaVSLambdaObs[i][1] << std::endl;

			double mean = deltaVSLambdaObs[i][0]/deltaVSLambdaObs[i][1];
			double err  = sqrt( (deltaVSLambdaObs[i][2]/deltaVSLambdaObs[i][1]-mean*mean)/deltaVSLambdaObs[i][3] );
			if (stepDefinition == 0) {
				mean = 0.;
				err  = 0.001;
			}
			
			hDeltaVsLambdaObs__[loopIdx]->SetBinContent(i+1,(mean+1.)*hDeltaVsLambdaObs__[loopIdxForHist]->GetBinContent(i+1) );
			hDeltaVsLambdaObs__[loopIdx]->SetBinError(i+1,err);
		}
	}
	fFile.close();


	// Method2: Set "hMeanFlux"
	double xxx1_value_m2 = 0.;
	double xxx2_value_m2 = 0.;
	double yyy1_value_m2 = 0.;
	double yyy2_value_m2 = 0.;
	double value_m2 = 0.;
	unsigned int nbEmptyPixels_m2 = 0;
	double valueFirstNotEmptyPixels = 0.;
	unsigned int idxFirstNotEmptyPixels = 0;
	for (unsigned int i=0; i<nbBinlambdaObs__; i++) {

		if (valueFirstNotEmptyPixels==0. && hDeltaVsLambdaObs__[loopIdx]->GetBinError(i+1)!=0.) {
			valueFirstNotEmptyPixels = hDeltaVsLambdaObs__[loopIdx]->GetBinContent(i+1);
			idxFirstNotEmptyPixels   = i;
		}
		
		// If empty and not the last pixel
		if (hDeltaVsLambdaObs__[loopIdx]->GetBinError(i+1) == 0. && i!=nbBinlambdaObs__-1) {
			if (nbEmptyPixels_m2 == 0 && i!=0) {
				xxx1_value_m2 = hDeltaVsLambdaObs__[loopIdx]->GetBinCenter(i);
				yyy1_value_m2 = hDeltaVsLambdaObs__[loopIdx]->GetBinContent(i);
				value_m2 = hDeltaVsLambdaObs__[loopIdx]->GetBinContent(i);
			}
			nbEmptyPixels_m2 ++;
		}
		else {
			
			// If first pixel not empty after some empty pixels
			if (nbEmptyPixels_m2!=0) {

				double a = 0.;
				double b = 0.;					
				
				// Find the mean value between the two not empty pixels
				if (value_m2!=0. && hDeltaVsLambdaObs__[loopIdx]->GetBinError(i+1)!=0.) {
					// Linear extrapolation
					xxx2_value_m2 = hDeltaVsLambdaObs__[loopIdx]->GetBinCenter(i+1);
					yyy2_value_m2 = hDeltaVsLambdaObs__[loopIdx]->GetBinContent(i+1);
					a = yyy1_value_m2 - yyy2_value_m2;
					b = yyy2_value_m2*xxx1_value_m2 - yyy1_value_m2*xxx2_value_m2;
					if (xxx1_value_m2-xxx2_value_m2 != 0.) {
						a /= (xxx1_value_m2-xxx2_value_m2);
						b /= (xxx1_value_m2-xxx2_value_m2);
					}
					value_m2 = (value_m2 + hDeltaVsLambdaObs__[loopIdx]->GetBinContent(i+1))/2.;
				}
				
				// Set all the empty pixels 
				for (unsigned int j=0; j<nbEmptyPixels_m2; j++) {
					double tmp_value = a*hDeltaVsLambdaObs__[loopIdx]->GetBinCenter(i-j) + b;
					if (hDeltaVsLambdaObs__[loopIdx]->GetBinError(i+1)==0.) tmp_value = value_m2;
					hDeltaVsLambdaObs__[loopIdx]->SetBinContent(i-j, tmp_value); //value_m2
					//hDeltaVsLambdaObs__[loopIdx]->SetBinError(i-j, maxError_m2);
				}
				
				// If the last pixel is also empty
				if (i==nbBinlambdaObs__-1) {
					hDeltaVsLambdaObs__[loopIdx]->SetBinContent(i+1, value_m2);
					//hDeltaVsLambdaObs__[loopIdx]->SetBinError(i+1, maxError_m2);
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
		hDeltaVsLambdaObs__[loopIdx]->SetBinContent(i+1, valueFirstNotEmptyPixels);
	}


	tmp_pathToSave  = pathToTxt__;
	tmp_pathToSave += "hDeltaVsLambdaObs_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += pathMoreForHist__;
	tmp_pathToSave += ".txt";
	fFile.open(tmp_pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(17);

	for (unsigned int i=0; i<nbBinlambdaObs__; i++) {
		fFile << i;
		fFile << " " << hDeltaVsLambdaObs__[loopIdx]->GetBinContent(i+1);
		fFile << " " << hDeltaVsLambdaObs__[loopIdx]->GetBinError(i+1) << std::endl;
	}
	fFile.close();


	//// For the weight
	tmp_pathToSave  = pathToTxt__;
	tmp_pathToSave += "eta_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += pathMoreForHist__;
	tmp_pathToSave += ".txt";
	fFile1.open(tmp_pathToSave.c_str());
	fFile1 << std::scientific;
	fFile1.precision(17);
	tmp_pathToSave  = pathToTxt__;
	tmp_pathToSave += "sigma2LSS_";
	tmp_pathToSave += forest__;
	tmp_pathToSave += pathMoreForHist__;
	tmp_pathToSave += ".txt";
	fFile2.open(tmp_pathToSave.c_str());
	fFile2 << std::scientific;
	fFile2.precision(17);



	std::vector<double> xxxRed;
	std::vector<double> xxxRedError;
	std::vector<double> yyyEta;
	std::vector<double> yyyEtaError;
	std::vector<double> yyySig;
	std::vector<double> yyySigError;

	for (unsigned int z=0; z<C_NBBINZ+2; z++) {

		unsigned int toBeSaved = 0;
		
		if (binCenterZZ[z][1] == 0.) continue;
		binCenterZZ[z][0] /= binCenterZZ[z][1];
			
		for (unsigned int i=0; i<nbBinPowErrorDelta__; i++) {
			
			if (arrayWeight[z][i][3] <= 2.) continue;
			arrayWeight[z][i][0] /= arrayWeight[z][i][2];
			arrayWeight[z][i][1] /= arrayWeight[z][i][2];

			const double var2 = arrayWeight[z][i][1]-arrayWeight[z][i][0]*arrayWeight[z][i][0];
			const double nb   = arrayWeight[z][i][3];
				
			if (var2 < 0.) std::cout << "  GetDelta.cpp::GetDelta::  var2Method2 < 0. ::: " << z << " " << i << " " << nb << " " << var2 << " " << arrayWeight[z][i][2] << std::endl;
			if (var2 <= 0.) continue;
			
			toBeSaved ++;
				
			hVarDeltaVsError[z]->SetBinContent(i+1, 1./var2);
			hVarDeltaVsError[z]->SetBinError(i+1, 1./(var2*sqrt(nb)) );
		}

		if (toBeSaved < 10) continue;

		// Get the two parameters for this redshift
		hVarDeltaVsError[z]->Fit(fitWeight,"QR");
		Double_t* tmp_paramFit      = fitWeight->GetParameters();
		Double_t* tmp_paramErrorFit = fitWeight->GetParErrors();

		fFile1 << binCenterZZ[z][0];
		fFile1 << " " << tmp_paramFit[0];
		fFile1 << " " << tmp_paramErrorFit[0] << std::endl;

		fFile2 << binCenterZZ[z][0];
		fFile2 << " " << tmp_paramFit[1];
		fFile2 << " " << tmp_paramErrorFit[1] << std::endl;


		xxxRed.push_back(binCenterZZ[z][0]);
		xxxRedError.push_back( sqrt( (binCenterZZ[z][2]/binCenterZZ[z][1]- binCenterZZ[z][0]*binCenterZZ[z][0])/binCenterZZ[z][3]) );
		yyyEta.push_back(tmp_paramFit[0]);
		yyyEtaError.push_back(tmp_paramErrorFit[0]);
		yyySig.push_back(tmp_paramFit[1]);
		yyySigError.push_back(tmp_paramErrorFit[1]);

		//R_plot1D(hVarDeltaVsError[z]);
		//R_plot1D(hVarDeltaVsError[z]);

	}
	fFile1.close();
	fFile2.close();

	/*
	//// Print the fit parameters
	for (unsigned int i=0; i<xxxRed.size(); i++) {
		std::cout << "  " << xxxRed[i] << " " <<  yyyEta[i] << " " << yyySig[i] << std::endl;
	}*/

	grEta__[loopIdx] = new TGraphErrors(xxxRed.size(), &xxxRed[0], &yyyEta[0], &xxxRedError[0], &yyyEtaError[0]);
	R_dealWithPlots(grEta__[loopIdx], "z", "#eta", "Method1: #eta(z)");	
	grSig__[loopIdx] = new TGraphErrors(xxxRed.size(), &xxxRed[0], &yyySig[0], &xxxRedError[0], &yyySigError[0]);
	R_dealWithPlots(grSig__[loopIdx], "z", "#sigma^{2}_{LSS}", "Method1: #sigma^{2}_{LSS}(z)");
	
	// delete pointers
	for (unsigned int j=0; j<C_NBBINZ+2; j++) delete hVarDeltaVsError[j];
	delete fitWeight;

	return;
}
void GetDelta::updateDeltaVector(unsigned int loopIdx) {

	const unsigned int nbForest = v_zz__.size();

	if (loopIdx==0 && ( stepDefinition == 0  || (mockJMC__ && stepDefinition<2)) ) {
		for (unsigned int i=0; i<nbForest; i++) {

			const unsigned int nb = v_nbPixel__[i];
			double tmpMeanForestLambdaRF[2] = {0.};

			for (unsigned int j=0; j<nb; j++) {

				v_TEMPLATE__[i][j]          = hTemplate__[loopIdx]->Interpolate(v_LAMBDA_RF__[i][j]);
				const double tmp_template2  = (v_alpha2__[i] + v_beta2__[i]*(v_LAMBDA_RF__[i][j]-v_meanForestLambdaRF__[i]))*v_TEMPLATE__[i][j]*hDeltaVsLambdaObs__[loopIdx]->Interpolate(v_LAMBDA_OBS__[i][j]);
				v_DELTA__[i][j]             = ( v_NORM_FLUX__[i][j]/tmp_template2 -1. )/v_FLUX_DLA__[i][j];
				v_DELTA_IVAR__[i][j]        = v_NORM_FLUX_IVAR__[i][j]*tmp_template2*tmp_template2*v_FLUX_DLA__[i][j]*v_FLUX_DLA__[i][j];

				v_DELTA_WEIGHT__[i][j]  = std::max(0.,v_FACTORWEIGHT__[i][j]/(sigma2LSSStart+1./(etaStart*v_DELTA_IVAR__[i][j])));
				if (v_DELTA_WEIGHT__[i][j]==0.) continue;

				tmpMeanForestLambdaRF[0] += v_DELTA_WEIGHT__[i][j]*v_LAMBDA_RF__[i][j];
				tmpMeanForestLambdaRF[1] += v_DELTA_WEIGHT__[i][j];
			}

			if (tmpMeanForestLambdaRF[1]!=0.) v_meanForestLambdaRF__[i] = tmpMeanForestLambdaRF[0]/tmpMeanForestLambdaRF[1];
			else std::cout << "  GetDelta::updateDeltaVector::  ERROR::  tmpMeanForestLambdaRF[1]==0. , forest = " << i << std::endl;
		}
	}
	else {
		for (unsigned int i=0; i<nbForest; i++) {

			const unsigned int nb = v_nbPixel__[i];
			double tmpMeanForestLambdaRF[2] = {0.};

			for (unsigned int j=0; j<nb; j++) {

				v_TEMPLATE__[i][j]          = hTemplate__[loopIdx]->Interpolate(v_LAMBDA_RF__[i][j]);
				const double tmp_template2  = (v_alpha2__[i] + v_beta2__[i]*(v_LAMBDA_RF__[i][j]-v_meanForestLambdaRF__[i]))*v_TEMPLATE__[i][j]*hDeltaVsLambdaObs__[loopIdx]->Interpolate(v_LAMBDA_OBS__[i][j]);
				v_DELTA__[i][j]             = ( v_NORM_FLUX__[i][j]/tmp_template2 -1. )/v_FLUX_DLA__[i][j];
				v_DELTA_IVAR__[i][j]        = v_NORM_FLUX_IVAR__[i][j]*tmp_template2*tmp_template2*v_FLUX_DLA__[i][j]*v_FLUX_DLA__[i][j];

				const double eta        = grEta__[loopIdx]->Eval(v_ZZZ__[i][j]);
				const double sigma2LSS  = grSig__[loopIdx]->Eval(v_ZZZ__[i][j]);
				v_DELTA_WEIGHT__[i][j]  = std::max(0.,v_FACTORWEIGHT__[i][j]/(sigma2LSS+1./(eta*v_DELTA_IVAR__[i][j])));
				if (v_DELTA_WEIGHT__[i][j]==0.) continue;

				tmpMeanForestLambdaRF[0] += v_DELTA_WEIGHT__[i][j]*v_LAMBDA_RF__[i][j];
				tmpMeanForestLambdaRF[1] += v_DELTA_WEIGHT__[i][j];
			}

			if (tmpMeanForestLambdaRF[1]!=0.) v_meanForestLambdaRF__[i] = tmpMeanForestLambdaRF[0]/tmpMeanForestLambdaRF[1];
			else std::cout << "  GetDelta::updateDeltaVector::  ERROR::  tmpMeanForestLambdaRF[1]==0. , forest = " << i << std::endl;
		}
	}

	return;
}






// ---------------------------------------------------------------------
//
//		Fit forests
//
// ---------------------------------------------------------------------
void GetDelta::fitForests(unsigned int begin, unsigned int end) {

	//// Constants
	const double sizeBinsZ = (maxRedshift__-minRedshift__)/nbBinsRedshift__;

	//// Vectors with the alpha, beta, chi^{2}
	std::vector<double> v_alpha(end-begin);
	std::vector<double> v_beta(end-begin);
	std::vector<double> v_chi2(end-begin);
	std::vector<double> v_alphaErr(end-begin);
	std::vector<double> v_betaErr(end-begin);
	std::vector<int> v_iflag(end-begin);

	//// Arrays with bin center of X axis of 'hFluxPDF__'
	double binCenterX[nbBinsFlux__];
	for (unsigned int i=0; i<nbBinsFlux__; i++) {
		binCenterX[i] = hFluxPDF__->GetXaxis()->GetBinCenter(i+1);
	}

	const unsigned int nbForest = v_zz__.size();
	if (nbForest==0) return;

	for (unsigned int f=0; f<nbForest; f++) {

		//  start minimization
		int iflag = 0;
		double arglist[10] = {0.};

		// fill arrays
		LambdaMeanMeth2 = v_meanForestLambdaRF__[f];
		const unsigned int nb = v_nbPixel__[f];
		NbPixMeth2 = 0;

		for (unsigned int p=0; p<nb; p++) {
			if (v_LAMBDA_RF__[f][p] >= lambdaRFMin__ && v_LAMBDA_RF__[f][p] < lambdaRFMax__) {

				
				FluxMeth2[NbPixMeth2]     = v_NORM_FLUX__[f][p];
				FluxErrMeth2[NbPixMeth2]  = 1./sqrt(v_NORM_FLUX_IVAR__[f][p]);
				LambdaMeth2[NbPixMeth2]   = v_LAMBDA_RF__[f][p];
				FluxMeanMeth2[NbPixMeth2] = v_TEMPLATE__[f][p];

				double zpix = v_ZZZ__[f][p];
				if (zpix < minRedshift__)       zpix = minRedshift__+sizeBinsZ;
				else if (zpix >= maxRedshift__) zpix = maxRedshift__-sizeBinsZ;

				for (unsigned int idp=0; idp<nbBinsFlux__; idp++) {
					a_PDF[idp][NbPixMeth2] = hFluxPDF__->Interpolate(binCenterX[idp],zpix);
				}

				NbPixMeth2++;
			}
		}

		mygMinuit->mnparm(0, "alpha", v_alpha2__[f], 0.1, minAlpha__, maxAlpha__, iflag);
		mygMinuit->mnparm(1, "beta" , v_beta2__[f],  0.1, minBeta__,  maxBeta__,  iflag);

		// Minimization
		arglist[0] = 1000.;
		arglist[1] = 0.1;
		mygMinuit->mnexcm("MIGRAD", arglist,2, iflag);
		double err[2]   = {0.};
		double param[2] = {0.};
		for (unsigned int ivar=0; ivar<2; ivar++) mygMinuit->GetParameter(ivar,param[ivar],err[ivar]);
		v_alpha[f] = param[0];
		v_beta[f]  = param[1];
		v_chi2[f]  = mygMinuit->fAmin;
		v_alphaErr[f] = err[0];
		v_betaErr[f]  = err[1];
		v_iflag[f] = iflag
	}
	
	//// alpha_beta
	std::ofstream fFile;
	std::string pathToSave  = pathToTxt__;
	pathToSave += "alphaAndBeta_";
	pathToSave += forest__;
	pathToSave += "_";

	std::stringstream convert0;
	convert0 << begin;
	pathToSave += convert0.str();

	pathToSave += "_";

	std::stringstream convert1;
	convert1 << end;
	pathToSave += convert1.str();

	pathToSave += pathMoreForHist__;
	pathToSave += ".txt";
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(17);

	std::cout << "  pathToSave = " << pathToSave << std::endl;

	//// index  alpha  beta  chi^{2}  error_alpha   error_beta   iflag
	for (unsigned int i=0; i<nbForest; i++) {
		fFile << i+begin;
		fFile << " " << v_alpha[i];
		fFile << " " << v_beta[i];
		fFile << " " << v_chi2[i];
		fFile << " " << v_alphaErr[i];
		fFile << " " << v_betaErr[i];
		fFile << " " << v_iflag[i];
		fFile << std::endl;
	}
	fFile.close();

	return;
}
extern void Chi2_method1(int &npar, double *gin, double &f, double *par, int iflag) {

	f = 0.;
	for (unsigned int i=0; i<NbPixMeth2; i++) {
		const double cont = (par[0]+par[1]*(LambdaMeth2[i]-LambdaMeanMeth2))*FluxMeanMeth2[i];
		const double tmp = (FluxMeth2[i]-cont)/FluxErrMeth2[i];
		f += tmp*tmp;
	}
}
extern void Chi2_method2(int &npar, double *gin, double &f, double *par, int iflag) {

	f = 0.;
	for (unsigned int i=0; i<NbPixMeth2; i++) {
		const double cont = (par[0]+par[1]*(LambdaMeth2[i]-LambdaMeanMeth2))*FluxMeanMeth2[i];
		f += -2.*ProbPixel(cont,FluxMeth2[i],FluxErrMeth2[i],i);
	}
}
double ProbPixel(double cont,double flux, double sig, unsigned int idxPixel) {
	
	const double dF = 1./nbBinsFlux__;
	double prb = 0.;
	cont *= dF/sig;
	flux /= sig;
	//const double fact1 = -0.5/(sig*sqrt(2.*M_PI));
	const double fact1 = -0.5/sqrt(2.*M_PI);

	for (unsigned int i=0; i<nbBinsFlux__; i++) {
		const double F   = i+0.5;
		const double PDF = a_PDF[i][idxPixel];
		prb += exp( fact1*pow(cont*F-flux,2.) )*PDF;
	}
	return log( prb*dF );
}
void GetDelta::initFitCont(void) {


	std::cout << "\n\n\n  Initialisation of the continuum fit with minuit (method 2) \n\n" << std::endl;

	mygMinuit = new TMinuit(10);
	mygMinuit->Clear();
	if (methodIndex__==1) mygMinuit->SetFCN(Chi2_method1);
	else if (methodIndex__==2) mygMinuit->SetFCN(Chi2_method2);

	int iflag = 0;
	double arglist[10] = {0.};
	arglist[0] = 1.;
	mygMinuit->mnexcm("SET ERR", arglist, 1, iflag);

	// Set starting values and step sizes for parameters

	mygMinuit->mnparm(0, "alpha", alphaStart__, 0.1, minAlpha__, maxAlpha__, iflag);
	mygMinuit->mnparm(1, "beta",  betaStart__,  0.1, minBeta__,  maxBeta__, iflag);
}








// ---------------------------------------------------------------------
//
//		Load data
//
// ---------------------------------------------------------------------
void GetDelta::loadDataForest(std::string fitsnameSpec, unsigned int start, unsigned int end, bool cutNotFittedSpectra /*=false*/) {

	//// Variables for FITS
	const TString TSfitsnameSpec = fitsnameSpec;
	std::cout << "  " << fitsnameSpec << std::endl;
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
	unsigned int nLines = end;
	if (end == 0 || end>=nrows) nLines = nrows;
	if (start >= nrows) return;

	std::cout << "  number of loaded forest = " << nLines-start << std::endl;

	double meanDelta[3] = {0.};
	v_fromFitsIndexToVectorIndex__.resize(nrows,-1);

	//// Load data
	for (unsigned int i=start; i<nLines; i++) { //nLines

		double tmp_meanDelta[3] = {0.};

		//// Variables for data in FITS
		double zz = 0.;
		double meanForestLambdaRF = 0.;
		double alpha2 = 0.;
		double beta2  = 0.;

		double LAMBDA_OBS[nbBinRFMax__];  
		double LAMBDA_RF[nbBinRFMax__];
		double NORM_FLUX[nbBinRFMax__];
		double NORM_FLUX_IVAR[nbBinRFMax__];
		double FLUX_DLA[nbBinRFMax__];
		double DELTA[nbBinRFMax__];
		double DELTA_IVAR[nbBinRFMax__];
		double DELTA_WEIGHT[nbBinRFMax__];
		double TEMPLATE[nbBinRFMax__];

		fits_read_col(fitsptrSpec,TDOUBLE, 6, i+1,1,1,NULL,&zz,                   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 8, i+1,1,1,NULL,&meanForestLambdaRF,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 11,i+1,1,1,NULL,&alpha2,               NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 12,i+1,1,1,NULL,&beta2,                NULL,&sta);

		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMax__,NULL, &LAMBDA_OBS,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 14,i+1,1,nbBinRFMax__,NULL, &LAMBDA_RF,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 15,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX_IVAR,  NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 17,i+1,1,nbBinRFMax__,NULL, &FLUX_DLA,        NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 18,i+1,1,nbBinRFMax__,NULL, &DELTA,           NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 19,i+1,1,nbBinRFMax__,NULL, &DELTA_IVAR,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 20,i+1,1,nbBinRFMax__,NULL, &DELTA_WEIGHT,    NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 21,i+1,1,nbBinRFMax__,NULL, &TEMPLATE,        NULL,&sta);

		if (stepDefinition >= 2 && cutNotFittedSpectra==true && ( (alpha2 == alphaStart__ && beta2 == betaStart__) || (fabs(alpha2)>=maxAlpha__-0.5) || (fabs(beta2)>=maxBeta__-0.05) ) ) continue;

		//// Vector with data
		bool templateHasNegative = false;
		std::vector< double > v_tmp_13;
		std::vector< double > v_tmp_14;
		std::vector< double > v_tmp_15;
		std::vector< double > v_tmp_16;
		std::vector< double > v_tmp_17;
		std::vector< double > v_tmp_18;
		std::vector< double > v_tmp_19;
		std::vector< double > v_tmp_20;
		std::vector< double > v_tmp_21;
		std::vector< double > v_tmp_zz;
		std::vector< double > v_tmp_factorWeight;
		unsigned int tmp_nb = 0;

		for (unsigned int j=0; j<nbBinRFMax__; j++) {
			if (NORM_FLUX_IVAR[j]>0. && FLUX_DLA[j]>=C_DLACORR && LAMBDA_OBS[j]>=lambdaObsMin__ && LAMBDA_OBS[j]<lambdaObsMax__ && DELTA_WEIGHT[j]>0.) {

				//// Remove veto lines
				bool isLine = false;
				if (doVetoLines__) {
					for (unsigned int k=0; k<nbVetoLines__; k++) {
						 if (LAMBDA_OBS[j]>=vetoLine__[2*k] && LAMBDA_OBS[j]<vetoLine__[2*k+1]) {
							isLine = true;
							continue;
						}
					}
				}
				if (isLine) continue;

				v_tmp_13.push_back(LAMBDA_OBS[j]);
				v_tmp_14.push_back(LAMBDA_RF[j]);
				v_tmp_15.push_back(NORM_FLUX[j]);
				v_tmp_16.push_back(NORM_FLUX_IVAR[j]);
				v_tmp_17.push_back(FLUX_DLA[j]);
				
				if (!mockJMC__ && stepDefinition == 0) {
					v_tmp_18.push_back(0.);
					v_tmp_19.push_back(0.);
					v_tmp_20.push_back(1.);
					v_tmp_21.push_back(1.);
				}
				else {
					v_tmp_18.push_back(DELTA[j]);
					v_tmp_19.push_back(DELTA_IVAR[j]);
					v_tmp_20.push_back(DELTA_WEIGHT[j]);
					v_tmp_21.push_back(TEMPLATE[j]);
				}
				
				const double zi = LAMBDA_OBS[j]/lambdaRFLine__-1.;
				v_tmp_zz.push_back(zi);
				v_tmp_factorWeight.push_back(pow((zi+1.)/onePlusZ0__, halfGama__) );

				///// Get Nb Pixel in forest
				if (DELTA_WEIGHT[j]>0. && LAMBDA_RF[j]>=lambdaRFMin__ && LAMBDA_RF[j]<lambdaRFMax__) {
					tmp_meanDelta[0] += DELTA_WEIGHT[j]*DELTA[j];
					tmp_meanDelta[1] += DELTA_WEIGHT[j];
					tmp_meanDelta[2] ++;

					tmp_nb++;
					if ( alpha2+beta2*(LAMBDA_RF[j]-meanForestLambdaRF) <= 0.) templateHasNegative = true;
				}

			}
		}

		if (cutNotFittedSpectra==true && (tmp_nb<C_MIN_NB_PIXEL || templateHasNegative ) ) continue;

		///// Get Nb Pixel in forest
		meanDelta[0] += tmp_meanDelta[0];
		meanDelta[1] += tmp_meanDelta[1];
		meanDelta[2] += tmp_meanDelta[2];

		v_zz__.push_back(zz);
		v_meanForestLambdaRF__.push_back(meanForestLambdaRF);
		v_alpha2__.push_back(alpha2);
		v_beta2__.push_back(beta2);
		v_nbPixel__.push_back(v_tmp_18.size());
		v_fromFitsIndexToVectorIndex__[i] = v_zz__.size()-1;

		v_LAMBDA_OBS__.push_back(     v_tmp_13);
		v_LAMBDA_RF__.push_back(      v_tmp_14);
		v_NORM_FLUX__.push_back(      v_tmp_15);
		v_NORM_FLUX_IVAR__.push_back( v_tmp_16);
		v_FLUX_DLA__.push_back(       v_tmp_17);
		v_DELTA__.push_back(          v_tmp_18);
		v_DELTA_IVAR__.push_back(     v_tmp_19);
		v_DELTA_WEIGHT__.push_back(   v_tmp_20);
		v_TEMPLATE__.push_back(       v_tmp_21);

		v_ZZZ__.push_back(v_tmp_zz);
		v_FACTORWEIGHT__.push_back(v_tmp_factorWeight);
	}

	fits_close_file(fitsptrSpec,&sta);

	std::cout << "  number of good   forest = " << v_zz__.size() << std::endl;

	std::cout << "  \n\n--- Step = -1 " << std::endl;
	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;

	return;
}













// ---------------------------------------------------------------------
//
//		Change data
//
// ---------------------------------------------------------------------
void GetDelta::updateDLA(std::string fitsnameSpec, unsigned int start, unsigned int end)
{

	//// Get the list of DLA
	int sta2 = 0;
	long nrows2 = 0;

	const TString TSfitsnameSpec2 = pathToDLACat__;
	std::cout << "  " << TSfitsnameSpec2 << std::endl;

	fitsfile* fitsptrSpec2;
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READONLY, &sta2); //READONLY // READWRITE
	fits_get_num_rows(fitsptrSpec2, &nrows2, &sta2);

	std::cout << "  number of DLA = " << nrows2 << std::endl;

	// list of DLA
	const unsigned int NbDLACat = nrows2;
	//unsigned int DLAplate[NbDLACat], DLAmjd[NbDLACat], DLAfiber[NbDLACat];
	double raDLACat[NbDLACat], deDLACat[NbDLACat], zabsDLACat[NbDLACat], NHIDLACat[NbDLACat];

	for (unsigned int i=0; i<NbDLACat; i++) {
		fits_read_col(fitsptrSpec2,TDOUBLE, 1, i+1,1,1,NULL,&raDLACat[i],    NULL,&sta2);
		fits_read_col(fitsptrSpec2,TDOUBLE, 2, i+1,1,1,NULL,&deDLACat[i],    NULL,&sta2);
		fits_read_col(fitsptrSpec2,TDOUBLE, 3, i+1,1,1,NULL,&zabsDLACat[i],  NULL,&sta2);
		//fits_read_col(fitsptrSpec2,TINT,    4, i+1,1,1,NULL,&DLAplate[i],    NULL,&sta2);
		//fits_read_col(fitsptrSpec2,TINT,    5, i+1,1,1,NULL,&DLAmjd[i],      NULL,&sta2);
		//fits_read_col(fitsptrSpec2,TINT,    6, i+1,1,1,NULL,&DLAfiber[i],    NULL,&sta2);
		fits_read_col(fitsptrSpec2,TDOUBLE, 7, i+1,1,1,NULL,&NHIDLACat[i],   NULL,&sta2);
	}
	fits_close_file(fitsptrSpec2,&sta2);

	const TString TSfitsnameSpec = fitsnameSpec;
	std::cout << "  " << fitsnameSpec << std::endl;

	//// Variables for FITS
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READWRITE, &sta); //READONLY // READWRITE
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
	unsigned int nLines = end;
	if (end == 0) nLines = nrows;

	std::cout << "  number of forest = " << nLines << std::endl;
	
	//// Get number of detected DLA
	unsigned int nbDLAAllCat = 0;

	//// Load data
	for (unsigned int i=start; i<nLines; i++) {

		//// Variables for data in FITS
		//unsigned int plate = 0;
		//unsigned int mjd   = 0;
		//unsigned int fiber = 0;
		double ra, de;

		double LAMBDA_OBS[nbBinRFMax__];
		double FluxDLA[nbBinRFMax__];
		
		//fits_read_col(fitsptrSpec,TINT, 1, i+1,1,1,NULL,&plate, NULL,&sta);
		//fits_read_col(fitsptrSpec,TINT, 2, i+1,1,1,NULL,&mjd,   NULL,&sta);
		//fits_read_col(fitsptrSpec,TINT, 3, i+1,1,1,NULL,&fiber, NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 4,i+1,1,1,NULL,&ra,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 5,i+1,1,1,NULL,&de,   NULL,&sta);

		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMax__,NULL, &LAMBDA_OBS,      NULL,&sta);

		for (unsigned int p=0; p<nbBinRFMax__; p++) {
			FluxDLA[p] = 1.;
		}

		// match with DLA catalog
		unsigned int NbDLA=0;
		float zabsDLA[NbDLA0__], NHIDLA[NbDLA0__];
		for (unsigned int iii=0; iii<NbDLA0__; iii++){
			zabsDLA[iii]=0.0;
			NHIDLA[iii]=0.0;
		}
		for (unsigned int iii=0; iii<NbDLACat; iii++) {

			/*
			if (plate != DLAplate[iii]) continue;
			if (mjd   != DLAmjd[iii])   continue;
			if (fiber != DLAfiber[iii]) continue;
			*/
			if (ra!=raDLACat[iii] || de!=deDLACat[iii]) continue;

			zabsDLA[NbDLA]   = zabsDLACat[iii];
			NHIDLA[NbDLA]    = NHIDLACat[iii];

			// compute DLA flux 
			for (unsigned int p=0; p<nbBinRFMax__; p++) {
				FluxDLA[p] *= VoigtProfile(NHIDLA[NbDLA],LAMBDA_OBS[p],zabsDLA[NbDLA]);
			}
			NbDLA++;
			if (NbDLA>NbDLA0__) {
				std::cout << "  GetDelta::updateDLA::  ERROR:  NbDLA>NbDLA0__" << NbDLA << " " << NbDLA0__ << std::endl;
				break;
			}
		}
		//// Save
		if (NbDLA!=0) fits_write_col(fitsptrSpec,TDOUBLE, 17,i+1,1,nbBinRFMax__, &FluxDLA, &sta);
		nbDLAAllCat += NbDLA;
	}

	fits_close_file(fitsptrSpec,&sta);

	std::cout << "  nbDLAAllCat = " << nbDLAAllCat << std::endl;

	return;
}
void GetDelta::updateDelta(std::string fitsnameSpec, unsigned int loopIdx, unsigned int start, unsigned int end) {

	//// Variables for FITS
	const TString TSfitsnameSpec = fitsnameSpec;
	std::cout << "  " << fitsnameSpec << std::endl;
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READWRITE, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
	unsigned int nLines = end;
	if (end == 0) nLines = nrows;

	std::cout << "  number of loaded forest = " << nLines << std::endl;

	double meanDelta[3] = {0.};
	double meanDelta_Zero[3] = {0.};

	//// Load data
	for (unsigned int i=start; i<nLines; i++) { //nLines

		//// Variables for data in FITS
		double meanForestLambdaRF = 0.;
		double alpha2 = 0.;
		double beta2  = 0.;

		double LAMBDA_OBS[nbBinRFMax__];
		double LAMBDA_RF[nbBinRFMax__];
		double NORM_FLUX[nbBinRFMax__];
		double NORM_FLUX_IVAR[nbBinRFMax__];
		double FLUX_DLA[nbBinRFMax__];
		double DELTA[nbBinRFMax__];
		double DELTA_IVAR[nbBinRFMax__];
		double DELTA_WEIGHT[nbBinRFMax__];
		double TEMPLATE[nbBinRFMax__];

		fits_read_col(fitsptrSpec,TDOUBLE, 8, i+1,1,1,NULL,&meanForestLambdaRF,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 11,i+1,1,1,NULL,&alpha2,               NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 12,i+1,1,1,NULL,&beta2,                NULL,&sta);

		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMax__,NULL, &LAMBDA_OBS,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 14,i+1,1,nbBinRFMax__,NULL, &LAMBDA_RF,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 15,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX_IVAR,  NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 17,i+1,1,nbBinRFMax__,NULL, &FLUX_DLA,        NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 18,i+1,1,nbBinRFMax__,NULL, &DELTA,           NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 19,i+1,1,nbBinRFMax__,NULL, &DELTA_IVAR,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 20,i+1,1,nbBinRFMax__,NULL, &DELTA_WEIGHT,    NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 21,i+1,1,nbBinRFMax__,NULL, &TEMPLATE,        NULL,&sta);

		if (v_fromFitsIndexToVectorIndex__[i]!=-1) meanForestLambdaRF =  v_meanForestLambdaRF__[ v_fromFitsIndexToVectorIndex__[i] ];
		double tmp_template[nbBinRFMax__];
		double delta[nbBinRFMax__];
		double delta_ivar[nbBinRFMax__];
		double delta_weight[nbBinRFMax__];

		for (unsigned int j=0; j<nbBinRFMax__; j++) {

			if (NORM_FLUX_IVAR[j]<=0. || FLUX_DLA[j]<=C_FLUX_DLA_IS_ZERO) continue; 

			tmp_template[j]         = hTemplate__[loopIdx]->Interpolate(LAMBDA_RF[j]);
			const double tmp_template2 = (alpha2+beta2*(LAMBDA_RF[j]-meanForestLambdaRF))*tmp_template[j]*hDeltaVsLambdaObs__[loopIdx]->Interpolate(LAMBDA_OBS[j]);
			delta[j]                = ( NORM_FLUX[j]/tmp_template2 -1. )/FLUX_DLA[j];
			delta_ivar[j]           = NORM_FLUX_IVAR[j]*tmp_template2*tmp_template2*FLUX_DLA[j]*FLUX_DLA[j];

			const double zi         = LAMBDA_OBS[j]/lambdaRFLine__-1.;
			const double eta        = grEta__[loopIdx]->Eval(zi);
			const double sigma2LSS  = grSig__[loopIdx]->Eval(zi);
			delta_weight[j]         = std::max(0.,pow( (zi+1.)/onePlusZ0__, halfGama__)/(sigma2LSS+1./(eta*delta_ivar[j])));
			if (delta_weight[j]<=0.) {
				std::cout << "  GetDelta::updateDelta:: ERROR:: delta_weight[j]<=0.  , " << i << " " << j << " " << zi  << " " << LAMBDA_OBS[j] << " " << sigma2LSS << " " << eta << " " << delta_ivar[j] << " " << FLUX_DLA[j] << " " << NORM_FLUX_IVAR[j] << std::endl;
				if (LAMBDA_RF[j] >= lambdaRFMin__ && LAMBDA_RF[j] < lambdaRFMax__) meanDelta_Zero[2] ++;
				continue;
			}

			if (v_fromFitsIndexToVectorIndex__[i]!=-1 && delta_weight[j]>0. && NORM_FLUX_IVAR[j]>0. && FLUX_DLA[j]>=C_DLACORR && LAMBDA_RF[j]>=lambdaRFMin__ && LAMBDA_RF[j]<lambdaRFMax__ && LAMBDA_OBS[j]>=lambdaObsMin__ && LAMBDA_OBS[j]<lambdaObsMax__ && DELTA_WEIGHT[j]>0.) {
				meanDelta[0] += delta_weight[j]*delta[j];
				meanDelta[1] += delta_weight[j];
				meanDelta[2] ++;
			}

		}

		//// Set the new value of 'meanForestLambdaRF' if there are some pixels
		if (v_fromFitsIndexToVectorIndex__[i]!=-1) {
			fits_write_col(fitsptrSpec,TDOUBLE, 8,i+1,1,1, &meanForestLambdaRF, &sta);
		}

		fits_write_col(fitsptrSpec,TDOUBLE, 18,i+1,1,nbBinRFMax__, &delta,        &sta);
		fits_write_col(fitsptrSpec,TDOUBLE, 19,i+1,1,nbBinRFMax__, &delta_ivar,   &sta);
		fits_write_col(fitsptrSpec,TDOUBLE, 20,i+1,1,nbBinRFMax__, &delta_weight, &sta);
		fits_write_col(fitsptrSpec,TDOUBLE, 21,i+1,1,nbBinRFMax__, &tmp_template, &sta);


	}

	fits_close_file(fitsptrSpec,&sta);

	std::cout << "  \n\n--- Step = " << loopIdx+1 << std::endl;
	std::cout << "  < delta >       = " << meanDelta[0]/meanDelta[1] << std::endl;
	std::cout << "  sum(w_i)        = " << meanDelta[1]              << std::endl;
	std::cout << "  nb pixel        = " << (long long unsigned int)meanDelta[2]              << std::endl;
	std::cout << "  nb pixel (w==0) = " << (long long unsigned int)meanDelta_Zero[2]         << std::endl;
	std::cout << "  all pixel       = " << (long long unsigned int)(meanDelta[2]+meanDelta_Zero[2]) << "\n" << std::endl;


	return;
}
void GetDelta::updateFlux(std::string fitsnameSpec, unsigned int start, unsigned int end) {

	/*

	Load data into vectors

	*/

	//// Load CIV histo
	TH1D* hDeltaVsLambdaObsCIV = new TH1D("hDeltaVsLambdaObsCIV","",nbBinlambdaObs__,lambdaObsMin__,lambdaObsMax__);
	R_dealWithPlots_1D(hDeltaVsLambdaObsCIV, "#lambda_{Obs.} (A)", "Mean transmission flux", "Method2: mean transmission flux");
	for (unsigned int j=0; j<nbBinlambdaObs__; j++) {
		hDeltaVsLambdaObsCIV->SetBinContent(j+1,1.);
		hDeltaVsLambdaObsCIV->SetBinError(j+1,0.);
	}
	//// Put the value of the init of mean transmision flux
	std::string path = pathToTxt__;
	path += "hDeltaVsLambdaObs_CIV.txt";
	std::ifstream file(path.c_str());

	unsigned int idx;
	double val;
	double err;
	while (file) {
		file>>idx>>val>>err;
		if(file==0) break;

		hDeltaVsLambdaObsCIV->SetBinContent(idx+1,val);
		hDeltaVsLambdaObsCIV->SetBinError(idx+1,err);
	}
	file.close();
	
	R_plot1D(hDeltaVsLambdaObsCIV);
	R_plot1D(hDeltaVsLambdaObsCIV);




	const TString TSfitsnameSpec = fitsnameSpec;
	std::cout << "  " << fitsnameSpec << std::endl;

	//// Variables for FITS
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READWRITE, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
	unsigned int nLines = end;
	if (end == 0) nLines = nrows;

	std::cout << "  number of loaded forest = " << nLines << std::endl;

	double LAMBDA_OBS[nbBinRFMax__];
	double NORM_FLUX[nbBinRFMax__];
	double NORM_FLUX_IVAR[nbBinRFMax__];

	//// Load data
	for (unsigned int i=start; i<nLines; i++) { //nLines
		
		fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMax__,NULL, &LAMBDA_OBS,      NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 15,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX,       NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMax__,NULL, &NORM_FLUX_IVAR,  NULL,&sta);

		for (unsigned int j=0; j<nbBinRFMax__; j++) {

			if (NORM_FLUX_IVAR[j]<=0.) continue;
			const double devident = hDeltaVsLambdaObsCIV->Interpolate(LAMBDA_OBS[j]);
			NORM_FLUX[j] /= devident;
			NORM_FLUX_IVAR[j] *= devident*devident;
		}

		fits_write_col(fitsptrSpec,TDOUBLE, 15,i+1,1,nbBinRFMax__, &NORM_FLUX,        &sta);
		fits_write_col(fitsptrSpec,TDOUBLE, 16,i+1,1,nbBinRFMax__, &NORM_FLUX_IVAR,   &sta);


	}

	fits_close_file(fitsptrSpec,&sta);


	delete hDeltaVsLambdaObsCIV;

	return;
}
void GetDelta::putReobsTogether(std::string fitsnameSpec, unsigned int loopIdx)
{

	double isReobsFlag = -100.;

	const TString TSfitsnameSpec = fitsnameSpec;
	std::cout << "  " << fitsnameSpec << std::endl;

	//// Variables for FITS
	int sta    = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READWRITE, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	std::cout << "  number of loaded forest = " << nrows << std::endl;

	//// Array of Ra DE and bool
	std::vector<double> ra(nrows,0.);
	std::vector<double> de(nrows,0.);
	std::vector<bool> isReobs(nrows,false);

	//// Get RA and DE
	for (unsigned int i=0; i<nrows; i++) {
		fits_read_col(fitsptrSpec,TDOUBLE, 4, i+1,1,1,NULL,&ra[i], NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 5, i+1,1,1,NULL,&de[i], NULL,&sta);
	}

	std::cout << "  Start " << std::endl;

	//// Load data
	for (unsigned int i=0; i<nrows; i++) {

		if (isReobs[i]) continue;

		//// Get the primery
		double ra1 = ra[i];
		double de1 = de[i];
		std::vector<unsigned int> idxObservation;

		//// Get other observations
		for (unsigned int j=i+1; j<nrows; j++) {

			double tmpRa = ra[j];
			double tmpDe = de[j];

			if ( tmpRa==ra1 && tmpDe==de1) {
				idxObservation.push_back(j);
				isReobs[j] = true;
			}

		}

		//// If there is no other observations
		if (idxObservation.size()==0) continue;
		//// If there are other observations
		else {
			unsigned int NbLambda = 0;
			double LAMBDA_OBS[nbBinRFMax__];
			double DELTA[nbBinRFMax__];
			double DELTA_IVAR[nbBinRFMax__];
			double DELTA_WEIGHT[nbBinRFMax__];

			fits_read_col(fitsptrSpec,TDOUBLE, 13,i+1,1,nbBinRFMax__,NULL, &LAMBDA_OBS,      NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 18,i+1,1,nbBinRFMax__,NULL, &DELTA,           NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 19,i+1,1,nbBinRFMax__,NULL, &DELTA_IVAR,      NULL,&sta);

			//// Get number of pixels
			for (unsigned int p=0; p<nbBinRFMax__; p++) {
				if (LAMBDA_OBS[p]==0.) break;
				else NbLambda ++;
			}

			//// Get other observations
			for (unsigned int j=0; j<idxObservation.size(); j++) {

				unsigned int tmpNbLambda = 0;
				double tmpLAMBDA_OBS[nbBinRFMax__];
				double tmpDELTA[nbBinRFMax__];
				double tmpDELTA_IVAR[nbBinRFMax__];
				const unsigned int idx =  idxObservation[j];

				fits_read_col(fitsptrSpec,TDOUBLE, 13,idx+1,1,nbBinRFMax__,NULL, &tmpLAMBDA_OBS,      NULL,&sta);
				fits_read_col(fitsptrSpec,TDOUBLE, 18,idx+1,1,nbBinRFMax__,NULL, &tmpDELTA,           NULL,&sta);
				fits_read_col(fitsptrSpec,TDOUBLE, 19,idx+1,1,nbBinRFMax__,NULL, &tmpDELTA_IVAR,      NULL,&sta);

				//// Get number of pixels
				for (unsigned int p=0; p<nbBinRFMax__; p++) {
					if (tmpLAMBDA_OBS[p]==0.) break;
					else tmpNbLambda ++;
				}

				makeCoAdd(NbLambda, DELTA, LAMBDA_OBS, DELTA_IVAR, tmpNbLambda, tmpDELTA, tmpLAMBDA_OBS, tmpDELTA_IVAR);

				//// Set a useless value to tell that it is a re-obs
				fits_write_col(fitsptrSpec,TDOUBLE, 9,idx+1,1,1, &isReobsFlag, &sta);
			}

			//// Find the new DELTA_WEIGHT
			for (unsigned int p=0; p<nbBinRFMax__; p++) {

				const double zi         = LAMBDA_OBS[p]/lambdaRFLine__-1.;
				const double eta        = grEta__[loopIdx]->Eval(zi);
				const double sigma2LSS  = grSig__[loopIdx]->Eval(zi);
				DELTA_WEIGHT[p]         = pow( (zi+1.)/onePlusZ0__, halfGama__)/(sigma2LSS+1./(eta*DELTA_IVAR[p]));
			}

			fits_write_col(fitsptrSpec,TDOUBLE, 18,i+1,1,nbBinRFMax__, &DELTA,        &sta);
			fits_write_col(fitsptrSpec,TDOUBLE, 19,i+1,1,nbBinRFMax__, &DELTA_IVAR,   &sta);
			fits_write_col(fitsptrSpec,TDOUBLE, 20,i+1,1,nbBinRFMax__, &DELTA_WEIGHT, &sta);

		}		
	}

	fits_close_file(fitsptrSpec,&sta);

	return;
}
void GetDelta::makeCoAdd(unsigned int NbLambda, double* Flux, double* Lambda, double* ErrFlux, unsigned int NbLambdaObs, double* FluxObs, double* LambdaObs, double* ErrFluxObs) {

	// Compute shift
	const double LambdaRef = Lambda[0];
	const double LambdaNew = LambdaObs[0];

	const int tmp_NbLambda = std::min(NbLambda,NbLambdaObs);

	int ishift=0;
	if (LambdaNew>LambdaRef) {
		for (int ilam=0; ilam<tmp_NbLambda; ilam++) {
			if(Lambda[ilam]>LambdaNew){ 
				ishift=ilam-1;
				break;
			}
		}
		for (int ilam=ishift; ilam<tmp_NbLambda; ilam++) {
			sumFlux(Flux[ilam],ErrFlux[ilam],FluxObs[ilam-ishift],ErrFluxObs[ilam-ishift]);
		}
	}
	else {
		for (int ilam=0; ilam<tmp_NbLambda; ilam++) {
			if(LambdaObs[ilam]>LambdaRef){ 
				ishift=ilam-1;
				break;
			}
		}

		for (int ilam=0; ilam<tmp_NbLambda-ishift; ilam++) {
			sumFlux(Flux[ilam],ErrFlux[ilam],FluxObs[ilam+ishift],ErrFluxObs[ilam+ishift]);
		}
	}

	return;
}
void GetDelta::sumFlux (double& f1, double& s1, double f2, double s2) {

	if (s1<=0. && s2<=0.) return;

	if (s2<=0. || isnan(f2)) return;
	else {
		if ( s1>0. && !isnan(f1) ) {
			const double s = s1+s2;
			f1 = (f1*s1 + f2*s2 )/s;
			s1 = s;
		}
		else {
			f1 = f2;
			s1 = s2;
		}
	}
	return;
}






















double GetDelta::VoigtProfile(float nhi, float lamb, float z_abs) {

	const double gamma = 6.625e8;
	const double f = 0.4164;

	const double c1000 = 299792458.0; //m/s
	double b = 30.*1000.; // 30 km/s parametre Doppler
	nhi += 0.1;  // correction to tune distribution
	const double NN = pow(10,nhi);
	const double larf = lamb/(1+z_abs);

	const double u = (c1000/b)*(lambdaRFLine__/larf-1);

	const double a = lambdaRFLine__*1e-10*gamma/(4*M_PI*b);
	const double sig = sqrt(2.0);
	const double H = TMath::Voigt(u,sig,a*2.0);  // in root factor 2....
	b/=1000.;
	const double tau = 1.497e-15*NN*f*larf*H/b;

	double prof=exp(-tau);
	if(prof>0.999)prof=1.0;

	return prof;
}





