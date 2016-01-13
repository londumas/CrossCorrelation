//===================================================================================
//
//         FILE: RootHistoFunctions.cpp
//
//        USAGE: ---
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

#include "RootHistoFunctions.h"

#include "math.h"
#include <vector>

#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h" // Needed because doing 'graph->GetXaxis()->CenterTitle()'
#include "TString.h"
#include "TFile.h"
#include "TROOT.h" //->SetDirectory(gROOT);
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#include "../../Cpp/Library/mathFunctions.h"

void R_dealWithPlots(TCanvas* canvas)
{
	// Sets the canvas parameters
	//
	// TCanvas* canvas : canvas to deal with
	//
	// return void     :

	gStyle->SetPaperSize(20,10);
	gStyle->SetCanvasColor(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	gStyle->SetPadColor(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBorderSize(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(0);
	gStyle->SetStatColor(0);
	gStyle->SetStatW(0.35);
	gStyle->SetStatH(0.15);
	gStyle->SetTitleColor(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0);

	canvas->SetTickx();
	canvas->SetTicky();
	canvas->SetGridx();
	canvas->SetGridy();
}
void R_dealWithPlots(TH1* histo)
{
	// Sets the histo parameters
	//
	// TH1* histo  : TH1 to deal with
	//
	// return void :

	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->CenterTitle();
	histo->GetXaxis()->SetTitleSize(0.045);
	histo->GetYaxis()->SetTitleSize(0.045);
	histo->GetXaxis()->SetTitleColor(1);
	histo->GetYaxis()->SetTitleColor(1);
}
void R_dealWithPlots_1D(TH1* histo, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string title /* = "" */)
{
	// Sets the 1D histo parameters
	//
	// TH1* histo                    : TH1 to deal with
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :
	
	R_dealWithPlots(histo);
	
	histo->SetMarkerStyle(20);
	histo->SetLineColor(kBlack);
	histo->SetMarkerColor(kBlue);
	
	const TString txTitle = xTitle;
	const TString tyTitle = yTitle;
	const TString tTitle = title;
	
	histo->SetXTitle(txTitle);
	histo->SetYTitle(tyTitle);
	histo->SetTitle(tTitle);	
}
void R_dealWithPlots_2D(TH2* histo, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string zTitle /* = "" */, std::string title /* = "" */)
{
	// Sets the 2D histo parameters
	//
	// TH2* histo                    : TH2 to deal with
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string zTitle /* = "" */ : z axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :
	
	R_dealWithPlots(histo);

	histo->GetZaxis()->CenterTitle();
	histo->GetZaxis()->SetTitleSize(0.045);
	histo->GetZaxis()->SetTitleColor(1);
	histo->GetZaxis()->SetTitleOffset(0.75);
	
	const TString txTitle = xTitle;
	const TString tyTitle = yTitle;
	const TString tzTitle = zTitle;
	const TString tTitle = title;
	histo->SetXTitle(txTitle);
	histo->SetYTitle(tyTitle);
	histo->SetZTitle(tzTitle);
	histo->SetTitle(tTitle);
}
void R_dealWithPlots_2D(TGraph2D* graph, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string zTitle /* = "" */, std::string title /* = "" */)
{
	// Sets the 2D histo parameters
	//
	// TH2* histo                    : TH2 to deal with
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string zTitle /* = "" */ : z axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :
	
	

	graph->GetXaxis()->CenterTitle();
	graph->GetYaxis()->CenterTitle();
	graph->GetZaxis()->CenterTitle();
	graph->GetXaxis()->SetTitleSize(0.045);
	graph->GetYaxis()->SetTitleSize(0.045);
	graph->GetZaxis()->SetTitleSize(0.045);
	graph->GetXaxis()->SetTitleColor(1);
	graph->GetYaxis()->SetTitleColor(1);
	graph->GetZaxis()->SetTitleColor(1);
	
	graph->SetMarkerColor(kBlue);
	graph->SetMarkerStyle(20);
	
	const TString txTitle = xTitle;
	const TString tyTitle = yTitle;
	const TString tzTitle = zTitle;
	const TString tTitle = title;
	graph->GetXaxis()->SetTitle(txTitle);
	graph->GetYaxis()->SetTitle(tyTitle);
	graph->GetZaxis()->SetTitle(tzTitle);
	graph->SetTitle(tTitle);
}
void R_dealWithPlots(TLegend* legend)
{
	// Sets the legend parameters
	//
	// TLegend * legend : legend to deal with
	//
	// return void      :

	legend->SetLineColor(0);
	legend->SetFillColor(0);
	//legend->SetMargin(0.);
}
void R_dealWithPlots(TF1* fonct)
{
	// Sets the TF1 parameters
	//
	// TF1* fonct  : pointer to the function
	//
	// return void :
	
	fonct->SetLineStyle(1);
	fonct->SetLineWidth(2);
	fonct->SetLineColor(kRed);
}
void R_dealWithPlots(TGraph* graph, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string title /* = "" */)
{
	// Sets the TGraph parameters
	//
	// TGraph* graph                 : pointer to the TGraph
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :
	
	graph->GetXaxis()->CenterTitle();
	graph->GetYaxis()->CenterTitle();
	graph->GetXaxis()->SetTitleSize(0.045);
	graph->GetYaxis()->SetTitleSize(0.045);
	graph->GetXaxis()->SetTitleColor(1);
	graph->GetYaxis()->SetTitleColor(1);
	
	graph->SetMarkerColor(kBlue);
	graph->SetMarkerStyle(20);
	
	const TString txTitle = xTitle;
	const TString tyTitle = yTitle;
	const TString tTitle = title;
	graph->GetXaxis()->SetTitle(txTitle);
	graph->GetYaxis()->SetTitle(tyTitle);
	graph->SetTitle(tTitle);
}

void R_binLogX(TH1* h)
{
	// Rebin an histo in log scale
	// 
	// TProfile* h : pointer to the TProfile 
	//
	// return      : void
	
	TAxis *axis = h->GetXaxis();
	const unsigned int bins = axis->GetNbins();
	
	Axis_t from = axis->GetXmin();
	Axis_t to = axis->GetXmax();
	Axis_t width = (to - from) / bins;
	Axis_t *new_bins = new Axis_t[bins + 1];
	
	for (unsigned int i = 0; i <= bins; i++) {
		new_bins[i] = pow(10, from + i * width);
	}
	axis->Set(bins, new_bins);
	delete[] new_bins;
}

void R_plot1D(TH1* histo, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 1d histogram
	//
	// TH1 *histo                    : histogram to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	TCanvas* canvas = new TCanvas("canvas"," ",2,2,600,450);
	R_dealWithPlots(canvas);
	canvas->cd();

	R_dealWithPlots_1D(histo, xTitle, yTitle, title);

	histo->Draw("ep");

	canvas->DrawClone();
	delete canvas;

}
void R_plot1DArray(unsigned int nbplot, TH1D** histo, std::vector<std::string> legendPlot, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 1d histogram
	//
	// TH1 *histo                    : histogram to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	if (nbplot < 1) return;
	if (legendPlot.size() < nbplot) legendPlot.resize(nbplot, "");


	TCanvas* canvas = new TCanvas("canvas"," ",2,2,600,450);
	R_dealWithPlots(canvas);
	canvas->cd();
	
	TLegend* legend = new TLegend(0.5, 0.5, 0.4, 0.4);
	R_dealWithPlots(legend);

	for (unsigned int i=0; i<nbplot; i++) {
		
		const TString plotName = legendPlot[i];
		
		R_dealWithPlots_1D(histo[i], xTitle, yTitle, title);
		histo[i]->SetMarkerColor(R_HCOLOR[i%10]);
		histo[i]->SetLineColor(R_HCOLOR[i%10]);
		histo[i]->SetLineWidth(3);
		if (plotName != "") legend->AddEntry(histo[i], plotName,"lep");
		
		if (i == 0) histo[i]->Draw("C");
		else histo[i]->Draw("C SAME");
	}
		
	legend->Draw();
	canvas->DrawClone();
	delete canvas;

}
void R_plot1D(TGraph* graph, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 1d TGraph
	//
	// TGraph* graph                 : graph to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	TCanvas* canvas = new TCanvas("canvas"," ",2,2,600,450);
	R_dealWithPlots(canvas);
	canvas->cd();

	R_dealWithPlots(graph, xTitle, yTitle, title);

	graph->DrawClone("AP");

	canvas->DrawClone();
	delete canvas;
}
void R_plot1DArray(unsigned int nbplot, TGraphErrors** graph, std::vector<std::string> legendPlot, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 1d histogram
	//
	// TH1 *histo                    : histogram to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	if (nbplot < 1) return;
	if (legendPlot.size() < nbplot) legendPlot.resize(nbplot, "");


	TCanvas* canvas = new TCanvas("canvas"," ",2,2,600,450);
	R_dealWithPlots(canvas);
	canvas->cd();
	
	TLegend* legend = new TLegend(0.5, 0.5, 0.4, 0.4);
	R_dealWithPlots(legend);

	for (unsigned int i=0; i<nbplot; i++) {
		
		const TString plotName = legendPlot[i];
		
		R_dealWithPlots(graph[i], xTitle, yTitle, title);
		graph[i]->SetMarkerColor(R_HCOLOR[i]);
		if (plotName != "") legend->AddEntry(graph[i], plotName,"lep");
		
		if (i == 0) graph[i]->Draw("ap");
		else graph[i]->Draw("ep SAME");
	}
		
	legend->Draw();
	canvas->DrawClone();
	delete canvas;

}
void R_plot2D(TH2* histo,    std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string zTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 2D histogram
	//
	// TH2 *histo                    : histogram to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string zTitle /* = "" */ : z axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	//TCanvas* canvas = new TCanvas("canvas"," ",2,2,1200,600);
	TCanvas* canvas = new TCanvas("canvas"," ",2,2,600,600);
	R_dealWithPlots(canvas);
	canvas->cd();

	R_dealWithPlots_2D(histo, xTitle, yTitle, zTitle, title);

	histo->DrawCopy("col z"); //pollego2z //cont4 pol z //col z cont4

	canvas->DrawClone();
	delete canvas;
}
void R_plot2D(TProfile* histo,    std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string zTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 2D histogram
	//
	// TH2 *histo                    : histogram to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string zTitle /* = "" */ : z axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	//TCanvas* canvas = new TCanvas("canvas"," ",2,2,1200,600);
	TCanvas* canvas = new TCanvas("canvas"," ",2,2,600,600);
	R_dealWithPlots(canvas);
	canvas->cd();

	//R_dealWithPlots_2D(histo, xTitle, yTitle, zTitle, title);

	histo->DrawCopy(); //pollego2z //cont4 pol z //col z cont4

	canvas->DrawClone();
	delete canvas;
}
void R_plot2D(TMatrixD* matrix, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string zTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 2D matrix
	//
	// TMatrixD *matrix              : matrix to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string zTitle /* = "" */ : z axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	TCanvas* canvas = new TCanvas("canvas"," ",2,2,600,600);
	R_dealWithPlots(canvas);
	canvas->cd();

	//R_dealWithPlots_2D(matrix, xTitle, yTitle, zTitle, title);

	matrix->Draw("col z"); //pollego2z //cont4 pol z //col z

	canvas->DrawClone();
	delete canvas;
}
void R_plot2D(TMatrixDSym* matrix, std::string xTitle /* = "" */, std::string yTitle /* = "" */, std::string zTitle /* = "" */, std::string title /* = "" */)
{
	// Plot a 2D matrix
	//
	// TMatrixD *matrix              : matrix to plot
	// std::string xTitle /* = "" */ : x axis title
	// std::string yTitle /* = "" */ : y axis title
	// std::string zTitle /* = "" */ : z axis title
	// std::string title  /* = "" */ : histo title
	//
	// return void                   :

	// Convert the matrix into a TH2D
	const unsigned int nbBin = matrix->GetNcols();
	TH2D* h2D = new TH2D("h2D", "", nbBin, 0., nbBin, nbBin, 0., nbBin);
	for (unsigned int i=0; i<nbBin; i++){
		for (unsigned int j=0; j<nbBin; j++){
			h2D->SetBinContent(i+1, j+1, (*matrix)(i,j));
			h2D->SetBinError(i+1, j+1, 0.);
		}	
	}
	
	R_plot2D(h2D, xTitle, yTitle, zTitle, title);
	delete h2D;
}

void R_rescaleHisto1D(TH1* histo)
{
	// Rescales a 1D histo by multiplying the bin content by the square
	// of the bin center value.
	//
	// TH1 *histo  : TH1 to rescale
	//
	// return void :

	const unsigned int nbBin = histo->GetNbinsX();

	for (unsigned int idxBin=0; idxBin<nbBin; idxBin++) {
		const double coeff = pow( histo->GetBinCenter(idxBin+1), 2. );
		histo->SetBinContent(idxBin+1, coeff*histo->GetBinContent(idxBin+1));
		histo->SetBinError(idxBin+1, coeff*histo->GetBinError(idxBin+1));
	}
}
void R_rescaleHisto2D(TH2* histo)
{
	// Rescales a 2D histo by multiplying the (x,y) bin content by the square
	// of the bin center value.
	//
	// TH2 *histo  : TH2 to rescale
	//
	// return void :

	const unsigned int nbOfBinX = histo->GetNbinsX();
	const unsigned int nbOfBinY = histo->GetNbinsY();

	for (unsigned int idxBinX=0; idxBinX<nbOfBinX; idxBinX++) {
		for (unsigned int idxBinY=0; idxBinY<nbOfBinY; idxBinY++) {
			const double coeff = pow(histo->GetXaxis()->GetBinCenter(idxBinX+1),2.) + pow(histo->GetYaxis()->GetBinCenter(idxBinY+1),2.);
			histo->SetBinContent(idxBinX+1, idxBinY+1, histo->GetBinContent(idxBinX+1, idxBinY+1)*coeff);
			histo->SetBinError(idxBinX+1, idxBinY+1, histo->GetBinError(idxBinX+1, idxBinY+1)*coeff);
		}
	}
}
void R_rescaleHisto2DPolar(TH2* histo)
{
	// Rescales a 2D histo in polar coord by multiplying the (theta, r)
	// bin content by the square of the bin 'r' value.
	//
	// TH2 *histo  : TH2 to rescale (theta, r)
	//
	// return void :

	const unsigned int nbOfBinX = histo->GetNbinsX();
	const unsigned int nbOfBinY = histo->GetNbinsY();

	for (unsigned int idxBinY=0; idxBinY<nbOfBinY; idxBinY++) {
		const double coeff = pow(histo->GetYaxis()->GetBinCenter(idxBinY+1),2.);
		for (unsigned int idxBinX=0; idxBinX<nbOfBinX; idxBinX++) {
			histo->SetBinContent(idxBinX+1, idxBinY+1, histo->GetBinContent(idxBinX+1, idxBinY+1)*coeff);
			histo->SetBinError(idxBinX+1, idxBinY+1, histo->GetBinError(idxBinX+1, idxBinY+1)*coeff);
		}
	}
}
void R_rescaleCovariance(TH2* histo)
{
	// Rescales a 2D histo containing the covariance matrix
	// between the bins of the cross correlation
	//
	// TH2* histo  : covariance matrix
	//
	// return void :

	const unsigned int nbBin = histo->GetNbinsX();

	// Divides the values inside a cell by the value inside the two diagonals
	for (unsigned int idxBinX=0; idxBinX<nbBin; idxBinX++) {
		for (unsigned int idxBinY=0; idxBinY<nbBin; idxBinY++) { //idxBinX or nbBin
			
			// Not on the diagonal
			if (idxBinX == idxBinY) continue;
			
			const double divide = sqrt(fabs(histo->GetBinContent(idxBinX+1, idxBinX+1)*histo->GetBinContent(idxBinY+1, idxBinY+1)));
			if (divide == 0) continue;
			const double value = histo->GetBinContent(idxBinX+1, idxBinY+1);
			if (value == 0) continue;
			histo->SetBinContent(idxBinX+1, idxBinY+1, value/divide);
		}
	}
	/*
	// Put the values to the symetrics
	for (int idxBinX=0; idxBinX<nbBin; idxBinX++) {
		for (int idxBinY=idxBinX+1; idxBinY<nbBin; idxBinY++) {
			histo->SetBinContent(idxBinX+1, idxBinY+1, histo->GetBinContent(idxBinY+1, idxBinX+1));
		}
	}*/
	
	// Put the diagonal to one
	for (unsigned int idxBin=0; idxBin<nbBin; idxBin++) {
		if (histo->GetBinContent(idxBin+1, idxBin+1) == 0) continue;
		histo->SetBinContent(idxBin+1, idxBin+1, 1.);
	}

}

TMatrixDSym* R_getCorrelatMatrix(TMatrixDSym* mCovar)
{
	// Get the correlation matrix
	//
	//
	
	const unsigned int nbBin = mCovar->GetNcols();
	
	// Weighting matrix to pass from covariance to correlation
	TMatrixDSym* mWeight = new TMatrixDSym(nbBin);
	for (unsigned int i=0; i<nbBin; i++) {
		if ((*mCovar)(i,i) != 0.) (*mWeight)(i,i) = 1./sqrt( (*mCovar)(i,i) );
	}
	
	// Pass from the covariance matrix 'mCovar' to
	//  the correlation matrix 'mCorr' using 
	//  mCorr = mWeight*mCovar*mWeight;
	TMatrixDSym* mCorr = new TMatrixDSym(nbBin);
	for (unsigned int i=0; i<nbBin; i++) {
		
		(*mCorr)(i,i) = (*mCovar)(i,i)*pow((*mWeight)(i,i), 2.);
		
		for (unsigned int j=0; j<i; j++) {
			(*mCorr)(i,j) = (*mCovar)(i,j)*(*mWeight)(i,i)*(*mWeight)(j,j);
			(*mCorr)(j,i) = (*mCorr)(i,j);
		}
	}
	
	
	return mCorr;
}

TH1D* R_getHisto1D(std::string pathToFile, std::string histoName)
{
	// Get a 1D histo from a root file
	//
	// std::string pathToFile : path to the file
	// std::string histoName  : name of the histogram
	//
	// return TH1D*           : histo
	
	// Load the file
	const TString tpathToFile = pathToFile;
	TFile* loadFile = new TFile(tpathToFile);
	
	// Get the histo
	const TString thistoName = histoName;
	//TH1D* histo = new TH1D();
	TH1D* histo = (TH1D*)loadFile->Get(thistoName);
	
	// Have the histo availeble after closing the TFile
	histo->SetDirectory(gROOT);

	delete loadFile;
	
	return histo;
}
TH2D* R_getHisto2D(std::string pathToFile, std::string histoName)
{
	// Get a 2D histo from a root file
	//
	// std::string pathToFile : path to the file
	// std::string histoName  : name of the histogram
	//
	// return TH2D*           : histo
	
	// Load the file
	const TString tpathToFile = pathToFile;
	TFile* loadFile = new TFile(tpathToFile);
	
	// Get the histo
	const TString thistoName = histoName;
	//TH2D* histo = new TH2D();
	TH2D* histo = (TH2D*)loadFile->Get(thistoName);
	
	// Have the histo availeble after closing the TFile
	histo->SetDirectory(gROOT);
	
	delete loadFile;
	
	return histo;
}
TGraphErrors* R_getTGraphErrors(std::string pathToFile, std::string histoName)
{
	// Get a 2D histo from a root file
	//
	// std::string pathToFile : path to the file
	// std::string histoName  : name of the histogram
	//
	// return TH2D*           : histo
	
	// Load the file
	const TString tpathToFile = pathToFile;
	TFile* loadFile = new TFile(tpathToFile);
	
	// Get the histo
	const TString thistoName = histoName;
	//TH2D* histo = new TH2D();
	TGraphErrors* histo = (TGraphErrors*)loadFile->Get(thistoName);
	delete loadFile;
	
	return histo;
}
TMatrixDSym* R_getMatrixFromFile(std::string fileName, std::string matrixName)
{
	// Get a matrix in a file
	//
	//
	
	const TString tFileName = fileName;
	TFile* loadFile = new TFile(tFileName, "READ");
	
	const TString tmatrixName = matrixName;
	TMatrixDSym* tmp_matrix = (TMatrixDSym*)loadFile->Get(tmatrixName);
	
	const unsigned int nbBin = tmp_matrix->GetNcols();
	
	TH2D* tmp_histo = new TH2D("tmp_histo", "", nbBin, 0., nbBin, nbBin, 0., nbBin);
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			tmp_histo->SetBinContent(i+1, j+1, (*tmp_matrix)(i,j));
		}
	}
	tmp_histo->SetDirectory(gROOT);
	
	delete tmp_matrix;
	delete loadFile;
	
	TMatrixDSym* matrix = new TMatrixDSym(nbBin);
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			(*matrix)(i,j) = tmp_histo->GetBinContent(i+1, j+1);
		}
	}
	delete tmp_histo;
	
	return matrix;
}

TH2D* R_getTH2FromMatrix(TMatrixDSym* matrix)
{
	
	const unsigned int nbBin = matrix->GetNcols();
	
	TH2D* histo = new TH2D("histo", "", nbBin, 0., nbBin, nbBin, 0., nbBin);
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			histo->SetBinContent(i+1, j+1, (*matrix)(i,j));
		}
	}
	
	return histo;
	
}
TH1D* R_getMeanHisto(std::vector<TH1D*> histos)
{
	// Returns the mean histo from a vector of histos
	//
	// std::vector<TH1D*> histos : vector of histos
	//
	// return TH1D*              : pointer to the mean histo
	
	// Get number of histos
	const unsigned int nbHisto = histos.size();
	
	// Get number of bins of histos
	const unsigned int nbBin = histos[0]->GetNbinsX();
	
	// Array to find the mean histo
	double mean[nbBin][4];
	//setArrayValuesToZero(nbBin, mean);
	
	// Fill the mean array
	for (unsigned int hIdx=0; hIdx<nbHisto; hIdx++) {
		
		const unsigned int nbBinHIdx = histos[hIdx]->GetNbinsX();
		if (nbBin != nbBinHIdx) std::cout << "RootHistoFunctions.cpp : getMeanHisto : nbBin != nbBinHIdx" << std::endl;
		
		for (unsigned int binIdx=0; binIdx<nbBinHIdx; binIdx++) {
			
			const double value = histos[hIdx]->GetBinContent(binIdx+1);
			const double error = histos[hIdx]->GetBinError(binIdx+1);
			if (error != 0.) {
				const double weight = 1./pow(error,2.);
				M_computeMean(mean[binIdx], value, weight);
			}
			else {
				std::cout << "RootHistoFunctions.cpp : getMeanHisto : error == 0., value = " << value << std::endl;
			}
		}
	}
	
	// Parameters of the histos
	const double binSize  = histos[0]->GetXaxis()->GetBinCenter(2) - histos[0]->GetXaxis()->GetBinCenter(1);
	const double lowEdge  = histos[0]->GetXaxis()->GetBinCenter(1) - binSize/2.;
	const double highEdge = histos[0]->GetXaxis()->GetBinCenter(nbBin) + binSize/2.;
	
	// Mean Histo
	TH1D* tmp_hMean = new TH1D("tmp_hMean", "", nbBin, lowEdge, highEdge);
	// Fill the histo
	R_fillCorrelation(tmp_hMean, mean);
	
	return tmp_hMean;
}
void R_fillCorrelation(TH1D* histo, double mean[][4])
{
	// Fill a 1D histogram "histo" with the values inside "mean"
	//
	// TH1D* histo       : histogram to fill
	// double array[][4] : array with value
	//                     (index X, index of mean array)
	//
	// return            : void

	const unsigned int nbBin = histo->GetNbinsX();

	for (unsigned int binIdx=0; binIdx<nbBin; binIdx++) {
		if (mean[binIdx][2] != 0.) {
			
			// Finish the array
			M_computeMean(mean[binIdx], 0., 0., true);
			
			// Fill the histo
			const double value = mean[binIdx][0];
			const double error = mean[binIdx][1];
			histo->SetBinContent(binIdx+1, value);
			histo->SetBinError(binIdx+1, error);
		}
	}
}
void R_fill1DHisto(TH1* histo, double** mean, bool finMean /*=true*/, bool bError /* = false */)
{
	// Fill a 1D histogram "histo" with the values inside "mean"
	//
	// TH1* histo     : histogram to fill
	// double** array : array with value
	//                   (index X, index of mean array)
	// bool finMean /*=true*/ : should we finish the mean calculation ?
	// bool bError     : tell if the error is the spread of the distrib 
	//                    or the invers sum of square of errors
	//
	// return         : void

	const unsigned int nbBin = histo->GetNbinsX();

	for (unsigned int binIdx=0; binIdx<nbBin; binIdx++) {
		if (mean[binIdx][2] != 0.) {
			
			// Finish the array
			if (finMean) M_computeMean(mean[binIdx], 0., 0., true);
			
			// Fill the histo
			const double value = mean[binIdx][0];
			double error = mean[binIdx][1];
			if (bError) error = mean[binIdx][2];
			histo->SetBinContent(binIdx+1, value);
			histo->SetBinError(binIdx+1, error);
		}
	}
}
void R_fill2DHisto(TH2* histo, double*** mean, bool finMean /*=true*/, bool bError /* = false */)
{
	// Fill a 2D histogram "histo" with the values inside "mean"
	//
	// TH2* histo      : histogram to fill
	// double*** array : array with value
	//                   (index X, index of mean array)
	// bool finMean /*=true*/ : should we finish the mean calculation ?
	// bool bError     : tell if the error is the spread of the distrib 
	//                    or the invers sum of square of errors
	//
	// return          : void

	const unsigned int nbBinX = histo->GetNbinsX();
	const unsigned int nbBinY = histo->GetNbinsY();

	for (unsigned int binIdx=0; binIdx<nbBinX; binIdx++) {
		for (unsigned int binIdx2=0; binIdx2<nbBinY; binIdx2++) {
			if (mean[binIdx][binIdx2][2] != 0.) {
				
				// Finish the array
				if (finMean) M_computeMean(mean[binIdx][binIdx2], 0., 0., true);
				
				// Fill the histo
				const double value = mean[binIdx][binIdx2][0];
				double error = mean[binIdx][binIdx2][1];
				if (bError) error = mean[binIdx][binIdx2][2];
				histo->SetBinContent(binIdx+1, binIdx2+1, value);
				histo->SetBinError(binIdx+1, binIdx2+1, error);
			}
		}
	}
}


TH1D* R_getCrossCorrelatConstantR(TH2* histo, int idxBinR, double minVar)
{
	// Get \xi(r) = f(\mu^{2}) in a TH1D from a polar TH2D
	//
	// TH2D* histo   : \xi= f(variable,r)
	// int idxBinR   : index of the bin in 'r'
	// double minVar : minimum of the variable as a range to fit
	//
	// return TH1D*  : pointer to the histo with the slice
	
	
	// range of the X axis of the histogramm
	const unsigned int nbPixel = histo->GetNbinsX();
	const double binSize = histo->GetXaxis()->GetBinCenter(2) - histo->GetXaxis()->GetBinCenter(1);
	const double lowEdge = histo->GetXaxis()->GetBinCenter(1) - binSize/2.;
	const double highEdge = histo->GetXaxis()->GetBinCenter(nbPixel) + binSize/2.;
	
	TH1D* histo1D = new TH1D("histo1D", "", nbPixel , lowEdge, highEdge);
	
	for (unsigned int idxBinX=0; idxBinX<nbPixel; idxBinX++) {
		//if (fabs(histo->GetXaxis()->GetBinCenter(idxBinX+1)) < minVar) continue;
		histo1D->SetBinContent(idxBinX+1, histo->GetBinContent(idxBinX+1, idxBinR+1));
		histo1D->SetBinError(idxBinX+1, histo->GetBinError(idxBinX+1, idxBinR+1));
	}
	
	return histo1D;
}
TH1D** R_getCorrelatMultipol(TH2* histo, std::string variable, double minVar)
{
	// Get [\xi_{i}(r)] in a TH1D from a polar TH2D
	//
	// TH2D* histo             : \xi= f(variable,r)
	// TH1D* correlatMulti[3]  : \xi_{i} = f(r)
	// std::string variable    : variable used : "theta", "mu", or "muSquare"
	// double minVar           : minimum of the variable as a range to fit
	//
	// return void             :
	
	const unsigned int nbOfBinX = histo->GetNbinsX();
	const double binSizeX = histo->GetXaxis()->GetBinCenter(2) - histo->GetXaxis()->GetBinCenter(1);
	const double lowEdgeX = histo->GetXaxis()->GetBinCenter(1) - binSizeX/2.;
	const double highEdgeX = histo->GetXaxis()->GetBinCenter(nbOfBinX) + binSizeX/2.;
	
	const unsigned int nbOfBinY = histo->GetNbinsY();
	const double binSizeY = histo->GetYaxis()->GetBinCenter(2) - histo->GetYaxis()->GetBinCenter(1);
	const double lowEdgeY = histo->GetYaxis()->GetBinCenter(1) - binSizeY/2.;
	const double highEdgeY = histo->GetYaxis()->GetBinCenter(nbOfBinY) + binSizeY/2.;

	// stores the resulted xi_{0}, xi_{2} and xi_{4} into an array
	//  of TH1D
	const unsigned int nbMulti = 3;
	TH1D** correlatMulti = new TH1D*[nbMulti];
	for (unsigned int idxMulti=0;  idxMulti<nbMulti; idxMulti++) {
		TString name = "correlatMulti";
		name += idxMulti;
		correlatMulti[idxMulti] = new TH1D(name, "", nbOfBinY, lowEdgeY, highEdgeY);
	}

	
	// Function fitting the multipols
	TF1* fitFunc;
	if (variable == "theta") {
		fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*pow(cos(x),2.)-0.5)", lowEdgeX, highEdgeX);
		//fitFunc = new TF1("fitFunc","[0] + [1]*cos(x) + [2]*(1.5*pow(cos(x),2.)-0.5)", lowEdgeX, highEdgeX);
		//fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*pow(cos(x),2.)-0.5) + [2]*((35./8.)*pow(cos(x),4.) - (15./4.)*pow(cos(x),2.) + 3./8.) ", lowEdgeX, highEdgeX);
		//fitFunc->FixParameter(2,0.);
	}
	else if (variable == "mu") {
		//fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*pow(x,2.)-0.5)", lowEdgeX, highEdgeX);
		//fitFunc = new TF1("fitFunc","[0] + [1]*x + [2]*(1.5*pow(x,2.)-0.5)", lowEdgeX, highEdgeX);
		//fitFunc = new TF1("fitFunc","[0] + [1]*x + [2]*(1.5*pow(x,2.)-0.5)", lowEdgeX, highEdgeX);
		fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*pow(x,2.)-0.5) + [2]*((35./8.)*pow(x,4.) - (15./4.)*pow(x,2.) + 3./8.) ", lowEdgeX, highEdgeX);
		//fitFunc->FixParameter(1,0.);
		//fitFunc->FixParameter(2,0.);
	}
	else if (variable == "muSquare") {
		fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*x-0.5)", lowEdgeX, highEdgeX);
		//fitFunc = new TF1("fitFunc","[0] + [1]*sqrt(x) + [2]*(1.5*x-0.5)", lowEdgeX, highEdgeX);
		//fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*x-0.5) + [2]*((35./8.)*pow(x,2.) - (15./4.)*x + 3./8.) ", lowEdgeX, highEdgeX);
		//fitFunc->FixParameter(2,0.);
	}
	else return correlatMulti;
	
	//TF1* fitFunc = new TF1("fitFunc","[0]-0.5*[1]+(3./8.)*[2] + ((3./2.)*[1]-(30./8.)*[2])*x + (35./8.)*[2]*pow(x,2.)",0.,1.);
	//TF1* fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*pow(x,2.)-0.5)",-1.,1.); //+ [1]*x
	/*TF1* fLegendre0 = new TF1("fLegendre0","[0]",0.,1.);
	TF1* fLegendre1 = new TF1("fLegendre1","[0]*x",0.,1.);
	TF1* fLegendre2 = new TF1("fLegendre2","[0]*(1.5*pow(x,2.)-0.5)",0.,1.);*/
	//fitFunc->FixParameter(1,0.);
	//fitFunc->FixParameter(2,0.);
	
	for (unsigned int idxBinY=0; idxBinY<nbOfBinY; idxBinY++) {
		TH1D* histo1D = R_getCrossCorrelatConstantR(histo, idxBinY, minVar);
		histo1D->Fit(fitFunc,"QRN");
		delete histo1D;
		
		Double_t* paramFit = fitFunc->GetParameters();
		Double_t* paramFitError = fitFunc->GetParErrors();
		for (unsigned int idxParam=0; idxParam<3; idxParam++) {
			correlatMulti[idxParam]->SetBinContent(idxBinY+1, paramFit[idxParam]);
			correlatMulti[idxParam]->SetBinError(idxBinY+1, paramFitError[idxParam]);
		}
	}
	delete fitFunc;
	
	return correlatMulti;
}


