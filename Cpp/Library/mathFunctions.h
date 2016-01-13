//===================================================================================
//
//         FILE: mathFunctions.h
//
//        USAGE: #include "/home/helion/Documents/Thèse/Code/Cpp/Library/mathFunctions.h"
//               In a file where you need these functions
//
//  DESCRIPTION: Gathering of functions doing mathematical stuff
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
#include <string>
#include <math.h>

#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

	int M_findBinIndex(const double& value, const double& min, const double& max, const int& nbBin, const double& binSize);
	void M_computeMean(double* mean, const double& value, const double& weight = 1., bool end=false);
	void M_computeCovariance(double* mean, const double& value, const double& weight1=1., const double& weight2=1., bool end=false);
	void M_setArrayValuesToZero(const int& nbBin, double mean[][4]);
	
	double* M_init1DArray(int nbBinX);
	double** M_init2DArray(int nbBinX, int nbBinY);
	double*** M_init3DArray(int nbBinX, int nbBinY, int nbBinZ);
	double**** M_init4DArray(int nbBinX, int nbBinY, int nbBinZ, int nbBinT);
	double***** M_init5DArray(int nbBinX, int nbBinY, int nbBinZ, int nbBinT, int nbBinU);
	
	void M_set2DArrayValuesToZero(double** array, int nbBinX, int nbBinY);
	void M_set3DArrayValuesToZero(double*** array, int nbBinX, int nbBinY, int nbBinZ);
	
	void M_delete2DArray(double** array, int nbBinX, int nbBinY);
	void M_delete3DArray(double*** array, int nbBinX, int nbBinY, int nbBinZ);
	void M_delete5DArray(double***** array, int nbBinX, int nbBinY, int nbBinZ, int nbBinT, int nbBinU);
	
	std::vector<std::vector<double> > M_getVectorFromFile(std::string pathToTxt, const int nbColumn);
	std::vector<std::string> M_getPathFromFile(std::string pathToTxt, std::string pathStart);
	
	inline void M_computeMeanFastest(double* mean, double value, double weight= 1.) {
		mean[2] += weight;
		mean[3] += pow(weight,2.);
		double meanBefore = mean[0];
		double weightNorma = weight/mean[2];
		mean[0] += weightNorma*(value-mean[0]);
		mean[1] = (1.-weightNorma)*mean[1] + weightNorma*(value-meanBefore)*(value-mean[0]);
	}

#endif


