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


GetDelta::GetDelta(int argc, char** argv) {

	//get_Ra_Dec_DLA();
	//get_Catalogue();
	get_flux_vs_lambdaObs();
	//get_weigted_covar_matrix();

}
GetDelta::~GetDelta(void) {	
}

void GetDelta::get_Ra_Dec_DLA(void)
{

	std::string fitsnameSpec = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits";

	/// Get the list of DLA
	int sta2 = 0;
	long nrows2 = 0;

	const TString TSfitsnameSpec2 = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/DLA_all.fits";
	std::cout << TSfitsnameSpec2 << std::endl;

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

	const TString TSfitsnameSpec = fitsnameSpec;
	std::cout << fitsnameSpec << std::endl;

	/// Variables for FITS
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
void GetDelta::get_Catalogue(void)
{
	/// Constants
	const double degToArcsec = 3600.;
	const double crit  = 2./degToArcsec;
	//const double critZ = 0.001;

	/// Path to save
	std::string pathToSave = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_DR7_DR12_EBOSS.fits";

	/// Number of catalogues
	const unsigned int nbCat = 3;

	/// PathToCat
	std::vector<std::string> nameCat(nbCat);
	std::vector<std::string> pathToCat(nbCat);
	nameCat[0]   = "DR12";
	pathToCat[0] = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/DR12Q_v2_10.fits";
	nameCat[1]   = "DR7";
	pathToCat[1] = "/home/gpfs/manip/mnt/bao/hdumasde/Data/Catalogue/knownquasarstar.060910.fits";
	nameCat[2]   = "EBOSS";
	pathToCat[2] = "/home/gpfs/manip/mnt0607/bao/hdumasde/Data/Catalogue/QSO_EBOSS.fits";

	/// Index of data in FITS
	unsigned int idxData[nbCat][3];
	/// DR12
	idxData[0][0] = 2;
	idxData[0][1] = 3;
	idxData[0][2] = 8;
	/// DR7
	idxData[1][0] = 2;
	idxData[1][1] = 3;
	idxData[1][2] = 4;
	/// EBOSS
	idxData[2][0] = 1;
	idxData[2][1] = 2;
	idxData[2][2] = 3;

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

	/// Find re-observations inside own catalogues
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

				if ( (fabs(ra1-ra2) < crit) && (fabs(de1-de2) < crit) ) {
					bo[c][q2] = true;
				}
			}
		}
	}

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

					if ( (fabs(ra1-ra2) < crit) && (fabs(de1-de2) < crit) ) {
						bo[c2][q2] = true;
					}
				}
			}
		}
	}

	/// Print number of re-observed
	for (unsigned int c=0; c<nbCat; c++) {
		const unsigned int nb = ra[c].size();

		unsigned int n = 0;
		for (unsigned int q=0; q<nb; q++) {
			if (bo[c][q]) n++;
		}
		std::cout << "  Reobs " << c << " : " << n << std::endl;
	}


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
void GetDelta::get_flux_vs_lambdaObs(void)
{

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
void GetDelta::get_weigted_covar_matrix(void)
{

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


















