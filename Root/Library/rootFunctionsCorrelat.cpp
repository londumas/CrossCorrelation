//===================================================================================
//
//         FILE: rootFunctionsCorrelat.cpp
//
//        USAGE: ---
//
//  DESCRIPTION: Gathering of functions dealing with root objects
//               for the cross correlations
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

#include "rootFunctionsCorrelat.h"

#include "math.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
//#include "TAxis.h" // Needed because doing 'graph->GetXaxis()->CenterTitle()'

/*
#include "/home/helion/Documents/Thèse/Code/Root/Library/RootHistoFunctions.h"
#include "/home/helion/Documents/Thèse/Code/Cpp/Library/mathFunctions.h"
*/

#include "RootHistoFunctions.h"
#include "/../../Cpp/Library/mathFunctions.h"


#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"

#include <Riostream.h>
#include "TLegend.h"
#include "TLegendEntry.h"

#include "Math/IFunction.h"
#include <cmath>
#include "TSystem.h"




void RC_rescale1D(TH1* histo, std::vector<double> meanDist)
{
	// Rescales a 1D histo by multiplying the bin content by the square
	// of the bin center value.
	//
	// TH1 *histo  : TH1 to rescale
	//
	// return void :

	const unsigned int nbBin = histo->GetNbinsX();

	bool coeffFromVector = true;
	if (meanDist.size() != nbBin) coeffFromVector = false;

	for (unsigned int i=0; i<nbBin; i++) {
		
		double coeff = 0.;
		if (coeffFromVector) coeff = pow( meanDist[i], 2.);
		else coeff = pow( histo->GetBinCenter(i+1), 2. );
		
		histo->SetBinContent(i+1, histo->GetBinContent(i+1)*coeff);
		histo->SetBinError(i+1, histo->GetBinError(i+1)*coeff);
	}
}
void RC_rescale2D(TH2* histo, std::vector<double> meanRparal, std::vector<double> meanRperp)
{
	// Rescales a 2D histo by multiplying the (x,y) bin content by the square
	// of the bin center value.
	//
	// TH2 *histo  : TH2 to rescale
	//
	// return void :

	const unsigned int nbOfBinX = histo->GetNbinsX();
	const unsigned int nbOfBinY = histo->GetNbinsY();
	
	bool coeffFromVector = true;
	if (meanRparal.size() != nbOfBinX) coeffFromVector = false;
	if (meanRperp.size() != nbOfBinY) coeffFromVector = false;

	for (unsigned int i=0; i<nbOfBinX; i++) {
		for (unsigned int j=0; j<nbOfBinY; j++) {
			double coeff = 0.;
			if (coeffFromVector) coeff = pow( meanRparal[i], 2.) + pow( meanRperp[j], 2.);
			else coeff = pow(histo->GetXaxis()->GetBinCenter(i+1),2.) + pow(histo->GetYaxis()->GetBinCenter(j+1),2.);
			histo->SetBinContent(i+1, j+1, histo->GetBinContent(i+1, j+1)*coeff);
			histo->SetBinError(i+1, j+1, histo->GetBinError(i+1, j+1)*coeff);
		}
	}
}
void RC_rescale2DPolar(TH2* histo)
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
void RC_rescaleCovariance(TH2* histo)
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

TGraphErrors* RC_getXCorrelatTGraphFromTH1(TH1D* histo, std::vector<double> meanDist)
{
	// Returns a pointer to a TGraphErrors with the cross correlation
	//  with the abscisse corresponding to the meanDist values.
	//
	// TH1D* histo                  : cross correlation
	// std::vector<double> meanDist : <|s|>
	//
	// return TGraphErrors*
	
	const unsigned int nbBin = histo->GetXaxis()->GetNbins();
	
	bool coeffFromVector = true;
	if (meanDist.size() != nbBin) coeffFromVector = false;
	
	// Arrays for a TGraphErrors
	double xxx[nbBin];
	double xxxError[nbBin];
	double yyy[nbBin];
	double yyyError[nbBin];
	
	for (unsigned int i=0; i<nbBin; i++) {
		
		double xValue = 0.;
		if (coeffFromVector) xValue = meanDist[i];
		else xValue = histo->GetBinCenter(i+1);
		
		xxx[i] = xValue;
		xxxError[i] = 0.;
		yyy[i] = histo->GetBinContent(i+1);
		yyyError[i] = histo->GetBinError(i+1);
	}
	
	TGraphErrors* graph = new TGraphErrors(nbBin, xxx, yyy, xxxError, yyyError);

	return graph;
}
TH1D* RC_getXCorrelat1DFrom2D(TH2* histo)
{
	// Returns a pointer to a TH1D with #xi(|s|)
	//
	// TH2* histo   : #xi(s_{parallel}, s_{perp})
	//
	// return TH1D* :
	
	const unsigned int nbBinX = histo->GetXaxis()->GetNbins();
	const unsigned int nbBinY = histo->GetYaxis()->GetNbins();
	
	const double minY = histo->GetYaxis()->GetXmin();
	const double maxY = histo->GetYaxis()->GetXmax();
	
	double** mean = M_init2DArray(nbBinY, 4);
	
	TH1D* h1XCorrelat = new TH1D("h1XCorrelat", "", nbBinY, minY, maxY);
	//R_dealWithPlots_1D(h1XCorrelat, "|s| (Mpc.h^{-1})", "#xi(|s|)", "Cross Correlation from 2D : #xi(|s|)");

	for (unsigned int i=0; i<nbBinX; i++) {
		for (unsigned int j=0; j<nbBinY; j++) {
			
			const double xCorrelatError = histo->GetBinError(i+1, j+1);
			if (xCorrelatError == 0.) continue;
			const double weight = 1./pow(xCorrelatError, 2.);
			
			const double value = histo->GetBinContent(i+1, j+1);
			
			const double dist = sqrt( pow(histo->GetXaxis()->GetBinCenter(i+1),2.) + pow(histo->GetYaxis()->GetBinCenter(j+1),2.) );
			unsigned int binIdx = h1XCorrelat->GetXaxis()->FindBin(dist)-1;
			if (binIdx >= nbBinY) binIdx = nbBinY-1;
			
			M_computeMean(mean[binIdx], value, weight);
		}
	}
	
	// Fill h1XCorrelat
	for (unsigned int i=0; i<nbBinY; i++) {
		
		M_computeMean(mean[i], 0., 0., true);
		
		h1XCorrelat->SetBinContent(i+1, mean[i][0]);
		h1XCorrelat->SetBinError(i+1, mean[i][2]);
	}
	M_delete2DArray(mean, nbBinY, 4);
	
	return h1XCorrelat;
}

TH1D* RC_getCrossCorrelatConstantR(TH2* histo, unsigned int idxBinR, double minVar)
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
TH1D** RC_getCorrelatMultipol(TH2* histo, std::string variable, double minVar)
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

void RC_getLegendrePoly(TH2* histo)
{
	const unsigned int nbBin = 100;
	const double dist = 10.;
	TH1D* h1ConstR = new TH1D("h1ConstR", "", nbBin, -M_PI/2., M_PI/2.);
	
	for (unsigned int i=0; i<nbBin; i++) {
		const double theta = h1ConstR->GetBinCenter(i+1);
		const double rParal = dist*sin(theta);
		const double rPerp  = dist*cos(theta);
		const double value = histo->Interpolate(rParal, rPerp);
		h1ConstR->SetBinContent(i+1, value);
		h1ConstR->SetBinError(i+1, 0.);
	}
	
	TF1* fitFunc = new TF1("fitFunc","[0] + [1]*(1.5*pow(cos(x),2.)-0.5)", -M_PI/2., M_PI/2.);
	h1ConstR->Fit(fitFunc, "R");
	
	R_plot1D(h1ConstR);
	delete h1ConstR;
	
}




