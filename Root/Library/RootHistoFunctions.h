//===================================================================================
//
//         FILE: RootHistoFunctions.h
//
//        USAGE: #include "/home/helion/Documents/Thèse/Code/Root/Library/RootHistoFunctions.h"
//               In a file where you need these functions
//
//  DESCRIPTION: Gathering of functions dealing with root objects
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

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#ifndef ROOTHISTOFUNCTIONS_H
#define ROOTHISTOFUNCTIONS_H
	
	void R_dealWithPlots(TCanvas* canvas);
	void R_dealWithPlots(TH1* histo);
	void R_dealWithPlots_1D(TH1* histo, std::string xTitle="", std::string yTitle="", std::string title="");
	void R_dealWithPlots_2D(TH2* histo, std::string xTitle="", std::string yTitle="", std::string zTitle="", std::string title="");
	void R_dealWithPlots_2D(TGraph2D* graph, std::string xTitle="", std::string yTitle="", std::string zTitle="", std::string title="");
	void R_dealWithPlots(TLegend* legend);
	void R_dealWithPlots(TF1* fonct);
	void R_dealWithPlots(TGraph* graph, std::string xTitle="", std::string yTitle="", std::string title="");
	
	void R_binLogX(TH1* h);
	
	void R_plot1D(TH1* histo,    std::string xTitle = "", std::string yTitle = "", std::string title = "");
	void R_plot1DArray(unsigned int nbplot, TH1D** histo, std::vector<std::string> legendPlot, std::string xTitle = "", std::string yTitle = "", std::string title = "");
	void R_plot1D(TGraph* graph, std::string xTitle = "", std::string yTitle = "", std::string title = "");
	void R_plot1DArray(unsigned int nbplot, TGraphErrors** graph, std::vector<std::string> legendPlot, std::string xTitle="", std::string yTitle="", std::string title="");
	void R_plot2D(TH2* histo,    std::string xTitle = "", std::string yTitle = "", std::string zTitle = "", std::string title = "");
	void R_plot2D(TProfile* histo,    std::string xTitle = "", std::string yTitle = "", std::string zTitle = "", std::string title = "");
	void R_plot2D(TMatrixD* matrix, std::string xTitle = "", std::string yTitle = "", std::string zTitle = "", std::string title = "");
	void R_plot2D(TMatrixDSym* matrix, std::string xTitle = "", std::string yTitle = "", std::string zTitle = "", std::string title = "");
	
	void R_rescaleHisto1D(TH1* histo);
	void R_rescaleHisto2D(TH2* histo);
	void R_rescaleHisto2DPolar(TH2* histo);
	void R_rescaleCovariance(TH2* histo);
	
	TMatrixDSym* R_getCorrelatMatrix(TMatrixDSym* mCovar);
	
	TH1D* R_getHisto1D(std::string pathToFile, std::string histoName);
	TH2D* R_getHisto2D(std::string pathToFile, std::string histoName);
	TGraphErrors* R_getTGraphErrors(std::string pathToFile, std::string histoName);
	TMatrixDSym* R_getMatrixFromFile(std::string fileName, std::string matrixName);
	
	TH2D* R_getTH2FromMatrix(TMatrixDSym* matrix);
	TH1D* R_getMeanHisto(std::vector<TH1D*> histos);
	void R_fillCorrelation(TH1D* histo, double mean[][4]);
	void R_fill1DHisto(TH1* histo, double** mean,  bool finMean=true, bool bError=false);
	void R_fill2DHisto(TH2* histo, double*** mean, bool finMean=true, bool bError=false);
	
	TH1D* R_getCrossCorrelatConstantR(TH2* histo, int idxBinR, double minVar);
	TH1D** R_getCorrelatMultipol(TH2* histo, std::string variable, double minVar);
	
	const int R_HCOLOR[240]={kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2, 
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2,
							   kBlue,kRed,kViolet,kBlack,kGreen,kOrange,kCyan,kMagenta-5,kYellow, kBlue-3, kGray, kGreen+2};

#endif


