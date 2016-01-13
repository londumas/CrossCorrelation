//===================================================================================
//
//         FILE: mathFunctions.cpp
//
//        USAGE: ---
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

#include "mathFunctions.h"

#include "math.h"
#include <vector>
#include <string>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream

int M_findBinIndex(const double& value, const double& min, const double& max, const int& nbBin, const double& binSize)
{
	// Compute the index of a value between "min" and "max" of an array
	// composed of "nbBin" with a bin size of "binSize"
	//
	// const double& value   : value for which the index is computed
	// const double& min     : min value of the array
	// const double& max     : max value of the array
	// const int& nbBin      : nb of bin the array is devided into
	// const double& binSize : size of one bin
	//
	// return int            : index of the bin (-1 if out of range )

	int idxBin = -1;
	if ( (value >= min) && (value <= max) ) {
		idxBin = int( (nbBin-1)*(value-min-binSize*0.5)/(max-min-binSize) + 0.5 );
		if (idxBin==nbBin) idxBin -= 1;
	}
	else {
		std::cout << "f::findIndex:: value out of range : value = " << value << " when min = " << min << " and max = " << max << std::endl;
	}

	return idxBin;
}
void M_computeMean(double* mean, const double& value, const double& weight /* = 1. */, bool end /* = false */)
{
	// Compute the mean and the error on the mean when the sum of weight
	//  is one.
	// "nfs-uxsup.csx.cam.ac.uk/~fanf2/hermes/doc/antiforgery/stats.pdf"
	//
	// double * mean[4] : mean[0] is the mean
	//                    mean[1] the error on the mean
	//                    mean[2] sum of weight
	//                    mean[3] sum of square of weight
	// double value     : the value to add to the mean
	// double weight    : the weight to give to the value (1. by default)
	// bool end         : calcul the error on the mean from the std (false by default)
	//
	// return void
	
	if (!end) {
		if (weight != 0.) {
			mean[2] += weight;
			mean[3] += pow(weight,2.);
			const double meanBefore = mean[0];
			const double weightNorma = weight/mean[2];
			mean[0] += weightNorma*(value-mean[0]);
			// Variance = sigma^2
			mean[1] = (1.-weightNorma)*mean[1] + weightNorma*(value-meanBefore)*(value-mean[0]);
		}
	}
	else {
		if (mean[2] == 0.) return;
		else {
			const double value = mean[1];
			// Error on the mean e = sigma/sqrt(N) = sqrt(Var)*sqrt(sum(w^2))/sum(w)
			mean[1] = sqrt(mean[1]*mean[3])/mean[2];
			mean[2] = 1./sqrt(mean[2]);
			// variance = sigma^2
			mean[3] = value;
		}
	}
}
void M_computeCovariance(double* mean, const double& value, const double& weight1 /* = 1. */, const double& weight2 /* = 1. */, bool end /*=false*/)
{
	// Compute the covariance and the error on the covariance
	//
	// double * mean[3] : mean[0] is the mean
	//                    mean[1]
	//                    mean[2] is the sum of weight1
	//                    mean[3] is the sum of weight2
	// double value     : the value to add to the mean
	// double weight1   : the first weight to give to the value (1. by default)
	// double weight2   : the second weight to give to the value (1. by default)
	// bool end         :
	//
	// return void

	if (weight1*weight2 != 0.) {
		const double weightSumTotBefore = mean[2]*mean[3];
		mean[2] += weight1;
		mean[3] += weight2;
		const double weightSumTotNow = mean[2]*mean[3];
		mean[0] = (weightSumTotBefore*mean[0] + weight1*weight2*value)/weightSumTotNow;
		mean[1] += 1;
		
	}
}
void M_setArrayValuesToZero(const int& nbBin, double mean[][4])
{
	// Set all the values of an array to zero
	//
	// const int& nbBin : number of bin of the array
	// double mean[][4] : array itself
	//
	// return void      : 

	for (int idxBin=0; idxBin<nbBin; idxBin++) {
		for (int idxInfo=0; idxInfo<4; idxInfo++) {
			mean[idxBin][idxInfo] = 0.;
		}
	}
}

double* M_init1DArray(int nbBinX)
{
	// Set all the values of an array to zero
	//
	// const int& nbBin : number of bin of the array
	// double mean[][4] : array itself
	//
	// return void      : 

	double* array = new double[nbBinX];
	for (int i=0; i<nbBinX; i++) {
		array[i] = 0.;
	}
	
	return array;
}
double** M_init2DArray(int nbBinX, int nbBinY)
{
	// Set all the values of an array to zero
	//
	// const int& nbBin : number of bin of the array
	// double mean[][4] : array itself
	//
	// return void      : 

	double** array = new double*[nbBinX];
	for (int i=0; i<nbBinX; i++) {
	    array[i] = new double[nbBinY];
	    for (int j=0; j<nbBinY; j++) {
			array[i][j] = 0.;
		}
	}
	
	return array;
}
double*** M_init3DArray(int nbBinX, int nbBinY, int nbBinZ)
{
	// Set all the values of an array to zero
	//
	// const int& nbBin : number of bin of the array
	// double mean[][4] : array itself
	//
	// return void      : 

	double*** array = new double**[nbBinX];
	for (int i=0; i<nbBinX; i++) {
	    array[i] = new double*[nbBinY];
	    for (int j=0; j<nbBinY; j++) {
			array[i][j] = new double[nbBinZ];
			for (int k=0; k<nbBinZ; k++) {
				array[i][j][k] = 0.;
			}
		}
	}
	
	return array;
}
double**** M_init4DArray(int nbBinX, int nbBinY, int nbBinZ, int nbBinT)
{
	// Init a 4D array of double with all values to zero
	//
	// int nbBinX
	// int nbBinY
	// int nbBinZ
	// int nbBinT
	//
	// return double*****

	double**** array = new double***[nbBinX];
	for (int i=0; i<nbBinX; i++) {
	    array[i] = new double**[nbBinY];
	    for (int j=0; j<nbBinY; j++) {
			array[i][j] = new double*[nbBinZ];
			for (int k=0; k<nbBinZ; k++) {
				array[i][j][k] = new double[nbBinT];
				for (int l=0; l<nbBinT; l++) {
					array[i][j][k][l] = 0.;
				}
			}
		}
	}
	
	return array;
}
double***** M_init5DArray(int nbBinX, int nbBinY, int nbBinZ, int nbBinT, int nbBinU)
{
	// Init a 5D array of double with all values to zero
	//
	// int nbBinX
	// int nbBinY
	// int nbBinZ
	// int nbBinT
	// int nbBinU
	//
	// return double*****

	double***** array = new double****[nbBinX];
	for (int i=0; i<nbBinX; i++) {
	    array[i] = new double***[nbBinY];
	    for (int j=0; j<nbBinY; j++) {
			array[i][j] = new double**[nbBinZ];
			for (int k=0; k<nbBinZ; k++) {
				array[i][j][k] = new double*[nbBinT];
				for (int l=0; l<nbBinT; l++) {
					array[i][j][k][l] = new double[nbBinU];
					for (int m=0; m<nbBinU; m++) {
						array[i][j][k][l][m] = 0.;
					}
				}
			}
		}
	}
	
	return array;
}
void M_set2DArrayValuesToZero(double** array, int nbBinX, int nbBinY)
{
	// Set all the values of an array to zero
	//
	// const int& nbBin : number of bin of the array
	// double mean[][4] : array itself
	//
	// return void      : 

	for (int ix=0; ix<nbBinX; ix++) {
		for (int iy=0; iy<nbBinY; iy++) {
			array[ix][iy] = 0.;
		}
	}
}
void M_set3DArrayValuesToZero(double*** array, int nbBinX, int nbBinY, int nbBinZ)
{
	// Set all the values of an array to zero
	//
	// const int& nbBin : number of bin of the array
	// double mean[][4] : array itself
	//
	// return void      : 

	for (int ix=0; ix<nbBinX; ix++) {
		for (int iy=0; iy<nbBinY; iy++) {
			for (int iz=0; iz<nbBinZ; iz++) {
				array[ix][iy][iz] = 0.;
			}
		}
	}
}
void M_delete2DArray(double** array, int nbBinX, int nbBinY)
{
	for (int i=0; i<nbBinX; i++) {
	    delete [] array[i];
	}
	delete [] array;
}
void M_delete3DArray(double*** array, int nbBinX, int nbBinY, int nbBinZ)
{
	for (int i=0; i<nbBinX; i++) {
		for (int j=0; j<nbBinY; j++) {
		    delete [] array[i][j];
		}
		delete [] array[i];
	}
	delete [] array;
}
void M_delete5DArray(double***** array, int nbBinX, int nbBinY, int nbBinZ, int nbBinT, int nbBinU)
{
	// Delete a 5D array of double.
	//
	// double***** array
	// int nbBinX
	// int nbBinY
	// int nbBinZ
	// int nbBinT
	// int nbBinU
	//
	// return void

	for (int i=0; i<nbBinX; i++) {
	    for (int j=0; j<nbBinY; j++) {
			for (int k=0; k<nbBinZ; k++) {
				for (int l=0; l<nbBinT; l++) {
					delete [] array[i][j][k][l];
				}
				delete [] array[i][j][k];
			}
			delete [] array[i][j];
		}
		delete [] array[i];
	}
	delete [] array;
}

std::vector<std::vector<double> > M_getVectorFromFile(std::string pathToTxt, const int nbColumn)
{
	// Returns a vector with the content of a file
	//
	// std::string pathToTxt                    : path to the file
	// const int nbCulumn                       : number of columns
	//
	// return std::vector<std::vector<double> > : vector with the info
	//											   vec[i] = row i
	
	// Open a file where there is a list of information
	std::ifstream file(pathToTxt.c_str());

	std::vector<std::vector<double> > values(nbColumn);
	std::vector<double> tmp_values(nbColumn);

	while (file) {
		
		if (nbColumn == 1)      file>>tmp_values[0];
		else if (nbColumn == 2) file>>tmp_values[0]>>tmp_values[1];
		else if (nbColumn == 3) file>>tmp_values[0]>>tmp_values[1]>>tmp_values[2];
		
		if(file==0) break;
		else {
			for (int i=0; i<nbColumn; i++) {
				values[i].push_back(tmp_values[i]);
			}
		}
	}
	
	return values;
}
std::vector<std::string> M_getPathFromFile(std::string pathToTxt, std::string pathStart)
{
	// Returns a vector with the content of a file
	//
	// std::string pathToTxt                    : path to the file
	//
	// return std::vector<std::string> : vector with the info
	
	// Open a file where there is a list of information
	std::ifstream file(pathToTxt.c_str());

	std::vector<std::string> vector;

	while (file) {
		std::string tmp_string;
		file>>tmp_string;
		std::string string = pathStart;
		string.append(tmp_string);
		
		if(file==0) break;
		else vector.push_back(string);
	}
	
	return vector;
}


