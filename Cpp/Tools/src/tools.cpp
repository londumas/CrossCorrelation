//===================================================================================
//
//         FILE: tools.cpp
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

#include "tools.h"
#include "../../../Root/Library/RootHistoFunctions.h"
#include "../../../Cpp/Library/mathFunctions.h"
#include "../../../Constants/constants.h"
#include "../../../Constants/globalValues.h"
#include "../../../chain_annalys_delta/Correlation/src/Cosmology.h"


#include <fstream>
#include <iostream>	// std::cout
#include <sstream>	//stringstream
#include <stdlib.h>	//atoi
#include <unistd.h>	//getopt


/// ROOT
#include "fitsio.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMath.h"


Tools::Tools(int argc, char** argv) {

	//get_Ra_Dec_DLA();
	//look_DLA();
	//get_delta_nicolas();
	get_Catalogue();
	//get_flux_vs_lambdaObs();
	//get_weigted_covar_matrix();

}
Tools::~Tools(void) {	
}

void Tools::get_Ra_Dec_DLA(void) {

	std::string forest_file = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits";
	std::string DLA_file    = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DLA_all.fits";

	/// Get the list of DLA
	const TString TSfitsnameSpec2 = DLA_file;
	std::cout << TSfitsnameSpec2 << std::endl;
	int sta2 = 0;
	long nrows2 = 0;
	fitsfile* fitsptrSpec2;
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READONLY, &sta2); //READONLY // READWRITE
	fits_get_num_rows(fitsptrSpec2, &nrows2, &sta2);

	std::cout << "  number of DLA = " << nrows2 << std::endl;

	// list of DLA
	const unsigned int NbDLACat = nrows2;
	unsigned int DLAplate[NbDLACat], DLAmjd[NbDLACat], DLAfiber[NbDLACat], nbFound[NbDLACat];  
	double zabsDLACat[NbDLACat], NHIDLACat[NbDLACat], raDLACat[NbDLACat], deDLACat[NbDLACat];

	for (unsigned int i=0; i<NbDLACat; i++) {
		fits_read_col(fitsptrSpec2,TDOUBLE, 1, i+1,1,1,NULL,&raDLACat[i],    NULL,&sta2);
		fits_read_col(fitsptrSpec2,TDOUBLE, 2, i+1,1,1,NULL,&deDLACat[i],    NULL,&sta2);
		fits_read_col(fitsptrSpec2,TDOUBLE, 3, i+1,1,1,NULL,&zabsDLACat[i],  NULL,&sta2);
		fits_read_col(fitsptrSpec2,TINT,    4, i+1,1,1,NULL,&DLAplate[i],    NULL,&sta2);
		fits_read_col(fitsptrSpec2,TINT,    5, i+1,1,1,NULL,&DLAmjd[i],      NULL,&sta2);
		fits_read_col(fitsptrSpec2,TINT,    6, i+1,1,1,NULL,&DLAfiber[i],    NULL,&sta2);
		fits_read_col(fitsptrSpec2,TDOUBLE, 7, i+1,1,1,NULL,&NHIDLACat[i],   NULL,&sta2);
		nbFound[i] = 0.;
	}
	fits_close_file(fitsptrSpec2,&sta2);



	/// Variables for FITS
	const TString TSfitsnameSpec = forest_file;
	std::cout << TSfitsnameSpec << std::endl;
	int sta = 0;
	long nrows = 0;
	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta); //READONLY // READWRITE
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
	unsigned int nLines = nrows;
	std::cout << "  number of forest = " << nLines << std::endl;

	
	/// Load data
	for (unsigned int i=0; i<nrows; i++) {

		/// Variables for data in FITS
		unsigned int NSPEC_BOSS = 0;
		double ra = 0.;
		double de = 0.;
		
		fits_read_col(fitsptrSpec,TDOUBLE, 2, i+1,1,1,NULL,&ra,           NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 3, i+1,1,1,NULL,&de,           NULL,&sta);
		fits_read_col(fitsptrSpec,TINT, 23,i+1,1,1,NULL,&NSPEC_BOSS,   NULL,&sta);

		const unsigned int nbObs = NSPEC_BOSS+1;
		unsigned int plate[nbObs];
		unsigned int mjd[nbObs];
		unsigned int fiber[nbObs];


		if (NSPEC_BOSS>0) {
			/// Get other observations
			fits_read_col(fitsptrSpec,TINT, 24, i+1,1,NSPEC_BOSS,NULL,&plate, NULL,&sta);
			fits_read_col(fitsptrSpec,TINT, 25, i+1,1,NSPEC_BOSS,NULL,&mjd,   NULL,&sta);
			fits_read_col(fitsptrSpec,TINT, 26, i+1,1,NSPEC_BOSS,NULL,&fiber, NULL,&sta);
		}

		/// Get best observations
		fits_read_col(fitsptrSpec,TINT, 5, i+1,1,1,NULL,&plate[NSPEC_BOSS], NULL,&sta);
		fits_read_col(fitsptrSpec,TINT, 6, i+1,1,1,NULL,&mjd[NSPEC_BOSS],   NULL,&sta);
		fits_read_col(fitsptrSpec,TINT, 7, i+1,1,1,NULL,&fiber[NSPEC_BOSS], NULL,&sta);		

		// match with DLA catalog
		for (unsigned int ii=0; ii<nbObs; ii++) {

//std::cout << i << " " << ii << " " << plate[ii] << " " << mjd[ii] << " " << fiber[ii] << std::endl;

			const unsigned int tmpplate = plate[ii];
			const unsigned int tmpmjd   = mjd[ii];
			const unsigned int tmpfiber = fiber[ii];

			for (unsigned int iii=0; iii<NbDLACat; iii++) {

				if (tmpplate != DLAplate[iii]) continue;
				if (tmpmjd   != DLAmjd[iii])   continue;
				if (tmpfiber != DLAfiber[iii]) continue;

				raDLACat[iii] = ra;
				deDLACat[iii] = de;
				nbFound[iii] ++;
			}
		}
	}

	fits_close_file(fitsptrSpec,&sta);

	/// Look how many DLA have not been found
	for (unsigned int i=0; i<NbDLACat; i++) {
		if (nbFound[i]==0) {
			std::cout << zabsDLACat[i] << " " << DLAplate[i] << " " << DLAmjd[i] << " " << DLAfiber[i] << std::endl;
		}
	}

	/// Store the RA and DEC
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READWRITE, &sta2); //READONLY // READWRITE
	fits_get_num_rows(fitsptrSpec2, &nrows2, &sta2);

	for (unsigned int i=0; i<NbDLACat; i++) {
		fits_write_col(fitsptrSpec2,TDOUBLE, 1,i+1,1,1, &raDLACat[i], &sta2);
		fits_write_col(fitsptrSpec2,TDOUBLE, 2,i+1,1,1, &deDLACat[i], &sta2);
	}
	fits_close_file(fitsptrSpec2,&sta2);



	return;
}
void Tools::look_DLA(void) {

	const std::string pathToDLACat__ = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DLA_all.fits";
	const std::string pathForest__   = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA//FitsFile_DR12_Guy/DR12_primery/DR12_primery.fits";
	unsigned int start = 0;
	unsigned int end   = 0;

	Cosmology* cosmo = new Cosmology(C_H, C_OMEGAM, C_OMEGAB);
	TH1D* hConvertRedshDist = cosmo->createHistoConvertRedshDist(C_NBBINREDSH, C_ZEXTREMABINCONVERT0, C_ZEXTREMABINCONVERT1);

	TProfile* h_norm_flux            = new TProfile("h_norm_flux","" ,8000,-8000.,8000.);
	TProfile* h_delta                = new TProfile("h_delta",""     ,8000,-8000.,8000.);
	TProfile* h_delta_vs_dist        = new TProfile("h_delta_vs_dist",""     ,8000,-8000.,8000.);
	TProfile* h_correction           = new TProfile("h_correction","",8000,-8000.,8000.);

	TProfile* h_delta_vs_ratio       = new TProfile("h_delta_vs_ratio","",8000,0.,5.);

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
		fits_read_col(fitsptrSpec2,TDOUBLE, 7, i+1,1,1,NULL,&NHIDLACat[i],   NULL,&sta2);
	}
	fits_close_file(fitsptrSpec2,&sta2);

	const TString TSfitsnameSpec = pathForest__;
	std::cout << "  " << pathForest__ << std::endl;

	//// Variables for FITS
	int sta = 0;
	long nrows = 0;

	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
	unsigned int nLines = end;
	if (end == 0) nLines = nrows;

	std::cout << "  number of forest = " << nLines << std::endl;
	
	//// Get number of detected DLA
	unsigned int nbDLAAllCat = 0;

	//// Load data
	for (unsigned int i=start; i<nLines; i++) {

		//// Variables for data in FITS
		double ra, de, zz;
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
		
		fits_read_col(fitsptrSpec,TDOUBLE, 4,i+1,1,1,NULL,&ra,   NULL,&sta);
		fits_read_col(fitsptrSpec,TDOUBLE, 5,i+1,1,1,NULL,&de,   NULL,&sta);
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

		if ((alpha2 == alphaStart__ && beta2 == betaStart__) || (fabs(alpha2)>=maxAlpha__-0.5) || (fabs(beta2)>=maxBeta__-0.05) ) {
			continue;
		}

		// match with DLA catalog
		unsigned int NbDLA=0;
		double zabsDLA[NbDLA0__], NHIDLA[NbDLA0__];
		for (unsigned int iii=0; iii<NbDLA0__; iii++){
			zabsDLA[iii]=0.0;
			NHIDLA[iii]=0.0;
		}
		for (unsigned int iii=0; iii<NbDLACat; iii++) {

			if (ra!=raDLACat[iii] || de!=deDLACat[iii]) continue;

			zabsDLA[NbDLA]   = zabsDLACat[iii];
			NHIDLA[NbDLA]    = NHIDLACat[iii];

			// compute DLA flux 
			for (unsigned int p=0; p<nbBinRFMax__; p++) {  //C_DLACORR

				if (DELTA_WEIGHT[p]>0. && NORM_FLUX_IVAR[p]>0. && FLUX_DLA[p]>C_DLACORR && LAMBDA_OBS[p]>=lambdaObsMin__ && LAMBDA_OBS[p]<lambdaObsMax__ && LAMBDA_RF[p]>=lambdaRFMin__ && LAMBDA_RF[p]<lambdaRFMax__) {
					// Fill histos
					const double l = (zabsDLA[NbDLA]+1.)*lambdaRFLine__;
					h_norm_flux->Fill(LAMBDA_OBS[p]-l, NORM_FLUX[p],DELTA_WEIGHT[p]);
					h_delta->Fill(LAMBDA_OBS[p]-l, DELTA[p],DELTA_WEIGHT[p]);
					h_delta_vs_dist->Fill(hConvertRedshDist->Interpolate(LAMBDA_OBS[p]/lambdaRFLine__-1.)-hConvertRedshDist->Interpolate(zabsDLA[NbDLA]), DELTA[p],DELTA_WEIGHT[p]);
					h_correction->Fill(LAMBDA_OBS[p]-l,  VoigtProfile(NHIDLA[NbDLA],LAMBDA_OBS[p],zabsDLA[NbDLA]),DELTA_WEIGHT[p]);
					h_delta_vs_ratio->Fill(LAMBDA_OBS[p]/l , DELTA[p],DELTA_WEIGHT[p]);
				}
			}
			NbDLA++;
			if (NbDLA>NbDLA0__) {
				std::cout << "  GetDelta::updateDLA::  ERROR:  NbDLA>NbDLA0__" << NbDLA << " " << NbDLA0__ << std::endl;
				break;
			}
		}
		nbDLAAllCat += NbDLA;
	}

	fits_close_file(fitsptrSpec,&sta);

	std::cout << "  nbDLAAllCat = " << nbDLAAllCat << std::endl;


	R_plot1D(h_norm_flux,"h_norm_flux");
	R_plot1D(h_delta,"h_delta");
	R_plot1D(h_delta_vs_dist,"h_delta_vs_dist");
	R_plot1D(h_correction,"h_correction");
	R_plot1D(h_delta_vs_ratio,"h_delta_vs_ratio");

	return;

	return;
}
double Tools::VoigtProfile(double nhi, double lamb, double z_abs) {

	const double gamma = 6.625e8;
	const double f = 0.4164;

	const double c1000 = 299792458.0; //m/s
	double b = 30.*1000.; // 30 km/s parametre Doppler
	nhi += 0.1;  // correction to tune distribution
	const double NN = pow(10.,nhi);
	const double larf = lamb/(1.+z_abs);

	const double u = (c1000/b)*(C_C_LYA_LINE/larf-1.);

	const double a = C_C_LYA_LINE*1.e-10*gamma/(4.*M_PI*b);
	const double sig = sqrt(2.);
	const double H = TMath::Voigt(u,sig,a*2.);  // in root factor 2....
	b/=1000.;
	const double tau = 1.497e-15*NN*f*larf*H/b;

	double prof=exp(-tau);

	return prof;
}
void Tools::get_delta_nicolas(void) {

	std::string pathToDataForest = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/DR12_Nicolas/deltas/plate-*";
	std::string pathToSave       = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/LYA/DR12_Nicolas/delta.fits";

	/// Fits file where to save
	TString TSfitsnameSpec2 = pathToSave;
	std::cout << "  " << TSfitsnameSpec2 << std::endl;
	int sta2    = 0;
	fitsfile* fitsptrSpec2;
	fits_open_table(&fitsptrSpec2,TSfitsnameSpec2, READWRITE, &sta2);

	/// index of forest
	unsigned int forestIdx = 0;
	/// Get the mean_flux_in_forest
	double meanForestLRF = 1120.;

	/// Get the list 
	FILE *fp;
	char path[PATH_MAX];
	std::string command = "ls " + pathToDataForest;
	std::vector< std::string > listFiles;
	fp = popen(command.c_str(), "r");
	while (fgets(path, PATH_MAX, fp) != NULL) listFiles.push_back(path);

	/// Get nb of files
	const unsigned int nbFiles = listFiles.size();

	for (unsigned int fileIdx=0; fileIdx<nbFiles; fileIdx++) {

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
			if (tmp_nbPixels2<C_MIN_NB_PIXEL) continue;

			//PLATE   =                 4063
			//MJD     =                55364
			//FIBERID =                  360
			//RA      =    4.343726224078581
			//DEC     =   0.2856228179296067
			//Z       =   2.4519999027252197
			//## R_COMOV  LOGLAM  DELTA  WEIGHT  CONT  MSHA
			/// Variable from old FITS
			unsigned int plate;
			unsigned int mjd;
			unsigned int fiberid;
			double ra;
			double de;
			double zz;
			double alpha;
			double beta;
			double LAMBDA_OBS[tmp_nbPixels2];
			double DELTA[tmp_nbPixels2];
			double DELTA_WEIGHT[tmp_nbPixels2];
			double CONTINUUM[tmp_nbPixels2];
			fits_read_key(fitsptrSpec,TINT,"PLATE",  &plate,NULL,&sta);
			fits_read_key(fitsptrSpec,TINT,"MJD",    &mjd,NULL,&sta);
			fits_read_key(fitsptrSpec,TINT,"FIBERID",&fiberid,NULL,&sta);
			fits_read_key(fitsptrSpec,TDOUBLE,"RA",  &ra,NULL,&sta);
			fits_read_key(fitsptrSpec,TDOUBLE,"DEC", &de,NULL,&sta);
			fits_read_key(fitsptrSpec,TDOUBLE,"Z",   &zz,NULL,&sta);
			fits_read_key(fitsptrSpec,TDOUBLE,"P0", &alpha,NULL,&sta);
			fits_read_key(fitsptrSpec,TDOUBLE,"P1", &beta,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 2,1,1,tmp_nbPixels2,NULL, &LAMBDA_OBS,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 3,1,1,tmp_nbPixels2,NULL, &DELTA,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 4,1,1,tmp_nbPixels2,NULL, &DELTA_WEIGHT,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 5,1,1,tmp_nbPixels2,NULL, &CONTINUUM,NULL,&sta);
	
			/// Variables for new FITS
			double lambdaOBS[nbBinRFMax__ ];
			double lambdaRF[nbBinRFMax__ ];
			double norm_flux[nbBinRFMax__];
			double norm_flux_ivar[nbBinRFMax__];
			double delta[nbBinRFMax__];
			double delta_ivar[nbBinRFMax__];
			double delta_weight[nbBinRFMax__];
			double continuum[nbBinRFMax__];
			unsigned int nbPixel = 0;
			const double oneOverOnePlusZ = 1./(1.+zz);

			ra *= C_RADTODEG;
			de *= C_RADTODEG;
	
			for (unsigned int p=0; p<tmp_nbPixels2; p++) {

				LAMBDA_OBS[p] = pow(10.,LAMBDA_OBS[p]);
				const double lambdaRFd = LAMBDA_OBS[p]*oneOverOnePlusZ;
	
				lambdaOBS[nbPixel]      = LAMBDA_OBS[p];
				lambdaRF[nbPixel]       = lambdaRFd;
				norm_flux[nbPixel]      = (DELTA[p]+1.)*CONTINUUM[p];
				norm_flux_ivar[nbPixel] = 1.;
				delta[nbPixel]          = DELTA[p];
				delta_ivar[nbPixel]     = 1.;
				delta_weight[nbPixel]   = DELTA_WEIGHT[p];
				continuum[nbPixel]      = CONTINUUM[p];
				nbPixel ++;
			}

			/// Save data in second file
			fits_write_col(fitsptrSpec2,TINT,    1,forestIdx+1,1,1, &plate, &sta2);
			fits_write_col(fitsptrSpec2,TINT,    2,forestIdx+1,1,1, &mjd, &sta2);
			fits_write_col(fitsptrSpec2,TINT,    3,forestIdx+1,1,1, &fiberid, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 4,forestIdx+1,1,1, &ra, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 5,forestIdx+1,1,1, &de, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 6,forestIdx+1,1,1, &zz, &sta2);
			fits_write_col(fitsptrSpec2,TINT,    7,forestIdx+1,1,1, &nbPixel, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 8,forestIdx+1,1,1, &meanForestLRF, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 11,forestIdx+1,1,1, &alpha, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 12,forestIdx+1,1,1, &beta, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 13,forestIdx+1,1,nbPixel, &lambdaOBS, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 14,forestIdx+1,1,nbPixel, &lambdaRF, &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 15,forestIdx+1,1,nbPixel, &norm_flux,       &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 16,forestIdx+1,1,nbPixel, &norm_flux_ivar,   &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 18,forestIdx+1,1,nbPixel, &delta,       &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 19,forestIdx+1,1,nbPixel, &delta_ivar,   &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 20,forestIdx+1,1,nbPixel, &delta_weight,   &sta2);
			fits_write_col(fitsptrSpec2,TDOUBLE, 21,forestIdx+1,1,nbPixel, &continuum,   &sta2);
			forestIdx ++;

		}

		fits_close_file(fitsptrSpec,&sta);
		std::cout << "  nb forest =  " << forestIdx << std::endl;
	}

	fits_close_file(fitsptrSpec2,&sta2);
	std::cout << "  nb forest =  " << forestIdx << std::endl;

}
void Tools::get_Catalogue(void) {

	std::cout << "\n\n" << std::endl;

	/// Path to save
	std::string pathToSave = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR7_DR12_EBOSS_2016_05_24.fits";
	/// Number of catalogues
	const unsigned int nbCat = 3;

	/// PathToCat
	std::vector<std::string> nameCat(nbCat);
	std::vector<std::string> pathToCat(nbCat);
	nameCat[0]   = "DR12";
	pathToCat[0] = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits";
	nameCat[1]   = "EBOSS";
	pathToCat[1] = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_EBOSS_updated_2016_05_24.fits";
	nameCat[2]   = "DR7";
	pathToCat[2] = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/knownquasarstar.060910.fits";

	/// Index of data in FITS
	unsigned int idxData[nbCat][3];
	/// DR12
	idxData[0][0] = 2;
	idxData[0][1] = 3;
	idxData[0][2] = 8;
	/// EBOSS
	idxData[1][0] = 1;
	idxData[1][1] = 2;
	idxData[1][2] = 3;
	/// DR7
	idxData[2][0] = 2;
	idxData[2][1] = 3;
	idxData[2][2] = 4;

	/// Load data
	std::vector<std::vector<double > > ra(nbCat);
	std::vector<std::vector<double > > de(nbCat);
	std::vector<std::vector<double > > zz(nbCat);
	std::vector<std::vector<bool   > > bo(nbCat);

	for (unsigned int c=0; c<nbCat; c++) {

		std::vector<double> tmp_ra;
		std::vector<double> tmp_de;
		std::vector<double> tmp_zz;

		int sta    = 0;
		long nrows = 0;
		const TString TSfitsnameSpec = pathToCat[c];
		fitsfile* fitsptrSpec;
		fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
		fits_get_num_rows(fitsptrSpec, &nrows, &sta);
		std::cout << nameCat[c]  << " " << nrows << std::endl;

		for (unsigned int q=0; q<nrows; q++) { /// nrows

			double ra = 0.;
			double de = 0.;
			double zz = 0.;

			fits_read_col(fitsptrSpec,TDOUBLE, idxData[c][0], q+1,1,1,NULL,&ra, NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, idxData[c][1], q+1,1,1,NULL,&de, NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, idxData[c][2], q+1,1,1,NULL,&zz, NULL,&sta);

			//if (zz<=1.5) continue; /// 0.001
			tmp_ra.push_back(ra);
			tmp_de.push_back(de);
			tmp_zz.push_back(zz);
		}

		fits_close_file(fitsptrSpec,&sta);

		std::vector<bool> tmp_bo(tmp_ra.size(),false);

		/// Save data
		ra[c] = tmp_ra;
		de[c] = tmp_de;
		zz[c] = tmp_zz;
		bo[c] = tmp_bo;
	}

	/// Print numbers
	for (unsigned int c=0; c<nbCat; c++) {
		std::cout << "  " << c << " : " << ra[c].size() << std::endl;
	}
	std::cout << "\n" << std::endl;

	/// Find re-observations inside own catalogues (Don't do for DR12Q)
	for (unsigned int c=1; c<nbCat; c++) {
		const unsigned int nb = ra[c].size();

		for (unsigned int q1=0; q1<nb; q1++) {

			if (bo[c][q1]) continue;
			if (q1%5000 == 0) std::cout << c << " " << q1 << "  /  " << nb << std::endl;

			const double ra1 = ra[c][q1];
			const double de1 = de[c][q1];

			for (unsigned int q2=q1+1; q2<nb; q2++) {

				const double ra2 = ra[c][q2];
				const double de2 = de[c][q2];

				if ( (fabs(ra1-ra2) < C_AUTOCORRCRITDEG) && (fabs(de1-de2) < C_AUTOCORRCRITDEG) ) {
					bo[c][q2] = true;
				}
			}
		}
	}
	std::cout << "\n" << std::endl;

	/// Print number of re-observed
	for (unsigned int c=0; c<nbCat; c++) {
		const unsigned int nb = ra[c].size();

		unsigned int n = 0;
		for (unsigned int q=0; q<nb; q++) {
			if (bo[c][q]) n++;
		}
		std::cout << "  Auto reobs " << c << " : " << n << std::endl;
	}
	std::cout << "\n" << std::endl;

	/// Find re-observations
	for (unsigned int c1=0; c1<nbCat; c1++) {
		const unsigned int nb1 = ra[c1].size();

		for (unsigned int c2=c1+1; c2<nbCat; c2++) {
			const unsigned int nb2 = ra[c2].size();

			for (unsigned int q1=0; q1<nb1; q1++) {

				if (bo[c1][q1]) continue;
				if (q1%5000 == 0) std::cout << c1 << " " << c2 << " " << q1 << "  /  " << nb1 << " " << nb2 << std::endl;

				const double ra1 = ra[c1][q1];
				const double de1 = de[c1][q1];

				for (unsigned int q2=0; q2<nb2; q2++) {

					if (bo[c2][q2]) continue;
					const double ra2 = ra[c2][q2];
					const double de2 = de[c2][q2];

					if ( (fabs(ra1-ra2) < C_AUTOCORRCRITDEG) && (fabs(de1-de2) < C_AUTOCORRCRITDEG) ) {
						bo[c2][q2] = true;
					}
				}
			}
		}
	}
	std::cout << "\n" << std::endl;

	/// Print number of re-observed
	for (unsigned int c=0; c<nbCat; c++) {
		const unsigned int nb = ra[c].size();

		unsigned int n = 0;
		for (unsigned int q=0; q<nb; q++) {
			if (bo[c][q]) n++;
		}
		std::cout << "  Cross reobs " << c << " : " << n << std::endl;
	}
	std::cout << "\n" << std::endl;

	/// Store the data
	unsigned int nbFinal = 0;
	int sta    = 0;
	long nrows = 0;
	const TString TSfitsnameSpec = pathToSave;
	fitsfile* fitsptrSpec;
	fits_open_table(&fitsptrSpec,TSfitsnameSpec, READWRITE, &sta);
	fits_get_num_rows(fitsptrSpec, &nrows, &sta);

	for (unsigned int c=0; c<nbCat; c++) {
		const unsigned int nb = ra[c].size();
		for (unsigned int q=0; q<nb; q++) {

			if (bo[c][q]) continue;

			fits_write_col(fitsptrSpec,TDOUBLE, 1,nbFinal+1,1,1, &ra[c][q], &sta);
			fits_write_col(fitsptrSpec,TDOUBLE, 2,nbFinal+1,1,1, &de[c][q], &sta);
			fits_write_col(fitsptrSpec,TDOUBLE, 3,nbFinal+1,1,1, &zz[c][q], &sta);

			nbFinal ++;
		}
	}
	fits_close_file(fitsptrSpec,&sta);

	std::cout << "  nbFinal = " << nbFinal << std::endl;
}
void Tools::get_flux_vs_lambdaObs(void) {

	/// Constants 
	const unsigned int nbBinRFMax = 4844;
	const double lambdaRFNormaMax = 1285.;
	const double lambdaRFMin      = 1040.;
	const double lambdaRFMax      = 1200.;

	const std::string pathToFolder = "/home/gpfs/manip/mnt/bao/hdumasde/Data/LYA/FitsFile_DR12_Guy_all/";
	std::vector< TString> path(20,"");
	for (unsigned int f=0; f<path.size(); f++) {
		path[f] += pathToFolder;
		path[f] += "DR12_Guy_";
		path[f] += f*5000;
		path[f] += "_";
		path[f] += (f+1)*5000;
		path[f] += ".fits";
	}

	/// vector to put the result of fit
	std::vector<double> fitCoef;


	///
	TProfile* hDeltaVsLambdaObs = new TProfile("hDeltaVsLambdaObs","",10000, 0., 20000.);
	R_dealWithPlots_1D(hDeltaVsLambdaObs, "#lambda_{Obs.} (A)", "Mean transmission flux", "hDeltaVsLambdaObs");
	///
	TProfile* hDeltaVsLambdaObsInForest = new TProfile("hDeltaVsLambdaObsInForest","",10000, 0., 20000.);
	R_dealWithPlots_1D(hDeltaVsLambdaObsInForest, "#lambda_{Obs.} (A)", "Mean transmission flux in forest", "hDeltaVsLambdaObsInForest");
	///
	TProfile* hDeltaVsLambdaObsInForestDivided = new TProfile("hDeltaVsLambdaObsInForestDivided","",10000, 0., 20000.);
	R_dealWithPlots_1D(hDeltaVsLambdaObsInForestDivided, "#lambda_{Obs.} (A)", "Mean transmission flux in forest", "hDeltaVsLambdaObsInForestDivided");
	///
	TProfile* hDeltaVsLambdaRF = new TProfile("hDeltaVsLambdaRF","",10000, 0., 10000.);
	R_dealWithPlots_1D(hDeltaVsLambdaRF, "#lambda_{R.F.} (A)", "<#delta+1>", "hDeltaVsLambdaRF");
	///
	TProfile* hDeltaVsLambdaRFDivided = new TProfile("hDeltaVsLambdaRFDivided","",10000, 0., 10000.);
	R_dealWithPlots_1D(hDeltaVsLambdaRFDivided, "#lambda_{R.F.} (A)", "<#delta+1>", "hDeltaVsLambdaRFDivided");
	///
	TH1D* hCoef = new TH1D("hCoef","",2000, -100., 100.);
	R_dealWithPlots_1D(hCoef, "coef", "", "hCoef");



	for (unsigned int f=0; f<path.size(); f++) {

		/// Variables for FITS
		int sta = 0;
		long nrows = 0;

		const TString TSfitsnameSpec = path[f];
		std::cout << TSfitsnameSpec << std::endl;
		fitsfile* fitsptrSpec;
		fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
		fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
		/// Load data
		for (unsigned int i=0; i<nrows; i++) {
		
			double LAMBDA_RF[nbBinRFMax];
			double FLUX[nbBinRFMax];
			double FLUX_IVAR[nbBinRFMax];
			double mean[2] = {};

			fits_read_col(fitsptrSpec,TDOUBLE, 9,  i+1,1,nbBinRFMax,NULL,&LAMBDA_RF,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 10, i+1,1,nbBinRFMax,NULL,&FLUX,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 11, i+1,1,nbBinRFMax,NULL,&FLUX_IVAR,NULL,&sta);

			for (unsigned int p=0; p<nbBinRFMax; p++) {
				if (FLUX_IVAR[p]<=0.) continue;
				else {
					hDeltaVsLambdaRF->Fill( LAMBDA_RF[p], FLUX[p],FLUX_IVAR[p]);
					mean[0] += FLUX_IVAR[p]*FLUX[p];
					mean[1] += FLUX_IVAR[p];
				}
			}

			mean[0] /= mean[1];
			fitCoef.push_back(mean[0]);
			hCoef->Fill(mean[0]);
		}
		fits_close_file(fitsptrSpec,&sta);
	}

	unsigned int idxSpec = 0;
	/// Divide by the callibration
	for (unsigned int f=0; f<path.size(); f++) {

		/// Variables for FITS
		int sta = 0;
		long nrows = 0;

		TString TSfitsnameSpec = path[f];
		std::cout << TSfitsnameSpec << std::endl;
		fitsfile* fitsptrSpec;
		fits_open_table(&fitsptrSpec,TSfitsnameSpec, READONLY, &sta);
		fits_get_num_rows(fitsptrSpec, &nrows, &sta);
	
		/// Load data
		for (unsigned int i=0; i<nrows; i++) {
		
			double LAMBDA_OBS[nbBinRFMax];
			double LAMBDA_RF[nbBinRFMax];
			double FLUX[nbBinRFMax];
			double FLUX_IVAR[nbBinRFMax];
			const double mean = fitCoef[idxSpec];

			fits_read_col(fitsptrSpec,TDOUBLE, 8,  i+1,1,nbBinRFMax,NULL,&LAMBDA_OBS,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 9,  i+1,1,nbBinRFMax,NULL,&LAMBDA_RF,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 10, i+1,1,nbBinRFMax,NULL,&FLUX,NULL,&sta);
			fits_read_col(fitsptrSpec,TDOUBLE, 11, i+1,1,nbBinRFMax,NULL,&FLUX_IVAR,NULL,&sta);

			for (unsigned int p=0; p<nbBinRFMax; p++) {
				if (FLUX_IVAR[p]<=0. || LAMBDA_RF[p]<lambdaRFNormaMax) continue;
				else {
					const double coef = mean*hDeltaVsLambdaRF->Interpolate(LAMBDA_RF[p]);
					hDeltaVsLambdaObs->Fill(LAMBDA_OBS[p], FLUX[p]/coef, FLUX_IVAR[p]*coef*coef);
				}				
			}

			idxSpec ++;
		}
		fits_close_file(fitsptrSpec,&sta);
	}


	R_plot1D(hCoef, "coef", "", "hCoef");
	R_plot1D(hDeltaVsLambdaRF, "#lambda_{R.F.} (A)", "", "hDeltaVsLambdaRF");
	R_plot1D(hDeltaVsLambdaObs, "#lambda_{Obs.} (A)", "", "hDeltaVsLambdaObs");
}
void Tools::get_weigted_covar_matrix(void) {

	/// Constants
	const unsigned int nbBoot = 2000;
	const unsigned int nbBin1D  = 200;
	const unsigned int nbBinX = nbBin1D;
	const unsigned int nbBinY = 2*nbBin1D;
	//const unsigned int nbBin  = 50*100;
	const unsigned int nbBin  = 50;
	const double binSize = nbBin1D/nbBin;

	/// All
	double all[nbBoot][nbBin][4];
	for (unsigned int i=0; i<nbBoot; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<4; k++) {
				all[i][j][k] = 0.;
			}
		}
	}
	/// Mean
	double mean[nbBin][4];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<4; j++) {
			mean[i][j] = 0.;
		}
	}
	/// Covariance and Correlation and diagCov
	double cov[nbBin][nbBin][4];
	double cor[nbBin][nbBin];
	double diagCov[nbBin];
	for (unsigned int i=0; i<nbBin; i++) {
		diagCov[i] = 0.;
		for (unsigned int j=0; j<nbBin; j++) {
			for (unsigned int k=0; k<4; k++) {
				cov[i][j][k] = 0.;
			}
			cor[i][j] = 0.;
		}
	}

	/// Get path to data
	std::string pathToFile[nbBoot];
	for (unsigned int i=0; i<nbBoot; i++) {

		std::stringstream convert;
		convert << i;

		//pathToFile[i]  = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_correctedBadFit3/xi_delta_QSO_1D_LYA_QSO_bootstrap_";
		pathToFile[i]  = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/Txt/FitsFile_DR12_Guy_primery_nicolasEstimator/xi_delta_QSO_1D_LYA_QSO_subsampling_";
		//pathToFile[i]  = "/home/gpfs/manip/mnt0607/bao/hdumasde/Mock_JMLG/v1547/Box_000/Simu_000/Results/xi_delta_QSO_1D_LYA_QSO_subsampling_";
		pathToFile[i] += convert.str();
		pathToFile[i] += ".txt";
	}

	/// All: Get the data
	for (unsigned int i=0; i<nbBoot; i++) {

		unsigned int idx = 0;
		double save0 = 0.;
		double save1 = 0.;
		double save2 = 0.;
		double save3 = 0.;
		double save4 = 0.;
		double save5 = 0.;
		double save6 = 0.;

		ifstream fileData(pathToFile[i].c_str());
		while (fileData) {
			fileData >> save0 >> save1 >> save2 >> save3 >> save4 >> save5 >> save6;
			if (fileData==0) break;

			const unsigned int idx2 = idx/binSize;
			/*
			const unsigned int iX = i/int(nbBinY);
			const unsigned int iY = i%int(nbBinY);
			const unsigned int idX = iX/int(binSize);
			const unsigned int idY = iY/int(binSize);
			const unsigned int idx2 = idX*100 + idY;
			*/
			all[i][idx2][0] += save0;
			all[i][idx2][1] += save1;
			all[i][idx2][2] += save5;
			all[i][idx2][3] += save6;
			idx ++;
		}

		for (unsigned int j=0; j<nbBin; j++) {

			if (all[i][j][3]<=1.) continue;
			save0 = all[i][j][0];
			save1 = all[i][j][1];
			save2 = all[i][j][2];
			save3 = all[i][j][3];

			all[i][j][0] = save0/save2;
			all[i][j][1] = sqrt( (save1/save2-all[i][j][0]*all[i][j][0])/save3);
			all[i][j][2] = save3/(save1/save2-all[i][j][0]*all[i][j][0]);
		}
	}

	/// Mean
	for (unsigned int i=0; i<nbBoot; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			if (all[i][j][3]<=1.) continue;
			mean[j][0] += all[i][j][2]*all[i][j][0];
			mean[j][1] += all[i][j][2]*all[i][j][0]*all[i][j][0];
			mean[j][2] += all[i][j][2];
			mean[j][3] ++;
		}
	}
	/// Mean
	for (unsigned int i=0; i<nbBin; i++) {
		mean[i][0] /= mean[i][2];
		mean[i][1] = sqrt( (mean[i][1]/mean[i][2]-mean[i][0]*mean[i][0])/mean[i][3]);
	}
	

	/// Covariance
	for (unsigned int i=0; i<nbBoot; i++) {
		for (unsigned int j=0; j<nbBin; j++) {

			if (all[i][j][3]<=1.) continue;
			const double m1 = all[i][j][0]-mean[j][0];
			const double w1 = all[i][j][2];

			for (unsigned int k=0; k<nbBin; k++) {

				if (all[i][k][3]<=1.) continue;
				const double m2 = all[i][k][0]-mean[k][0];
				const double w2 = all[i][k][2];

				cov[j][k][0] += w1*w2*m1*m2;
				cov[j][k][1] += w1*w2;
				cov[j][k][2] += w1*w2*w1*w2;
				cov[j][k][3] ++;
			}
		}
	}

	/// Covariance: Find the values
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {

			const double save0 = cov[i][j][0];
			const double save1 = cov[i][j][1];
			const double save2 = cov[i][j][2];
			cov[i][j][0] = save0*save1/( save1*save1-save2 );
		}
		diagCov[i] = cov[i][i][0];
	}



	/// Plot mean xi
	TH1D* histo = new TH1D("histo","",nbBin,0.,200.);
	for (unsigned int k=0; k<3; k++) {
		
		for (unsigned int i=0; i<nbBin; i++) {
			const double coef = pow(histo->GetBinCenter(i+1),k);
			histo->SetBinContent(i+1, coef*mean[i][0]);
			histo->SetBinError(i+1, coef*mean[i][1]);
		}
		R_dealWithPlots_1D(histo);
		R_plot1D(histo);
	}



	/// Plot covariance and correlation matrix
	TH2D* histo2 = new TH2D("histo2","",nbBin,0.,200.,nbBin,0.,200.);
	TH2D* histo3 = new TH2D("histo3","",nbBin,0.,200.,nbBin,0.,200.);
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<nbBin; j++) {
			histo2->SetBinContent(i+1,j+1, cov[i][j][0]);
			histo2->SetBinError(i+1,j+1, 0.00001);
			histo3->SetBinContent(i+1,j+1, cov[i][j][0]/sqrt(diagCov[i]*diagCov[j])  );
			histo3->SetBinError(i+1,j+1, 0.00001);
		}
	}
	R_plot2D(histo2);
	R_plot2D(histo3);

	/// Plot the variogram
	TH1D* histo4 = new TH1D("histo4","",nbBin,0.,200.);
	double var[nbBin][2];
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<2; j++) {
			var[i][j] = 0.;
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		for (unsigned int j=0; j<i+1; j++) {
			var[i-j][0] += histo3->GetBinContent(i+1,j+1);
			var[i-j][1] ++; 
		}
	}
	for (unsigned int i=0; i<nbBin; i++) {
		var[i][0] /= var[i][1];
		histo4->SetBinContent(i+1, var[i][0]);
		histo4->SetBinError(i+1, 0.00001);
	}
	R_dealWithPlots_1D(histo4);
	R_plot1D(histo4);






}


















