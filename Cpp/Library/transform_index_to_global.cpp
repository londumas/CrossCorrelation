//===================================================================================
//
//         FILE: cpp.cpp
//
//        USAGE: compile with " g++ cpp.cpp -o  cpp.o "
//
//  DESCRIPTION: Transform a file of the 2D cross correlation :
//                 < i1 i2 xi(i1, i2) >
//                 into a global index based list
//                 < j xi(j) >
//                 where 0 <= i1 < N1 , 0 <= i2 < N2
//                       j = (i1*N2+i2)*N3+i3
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

#include <iostream>
#include "vector.h"
#include "string.h"
#include <fstream>      // std::ifstream

std::vector<std::vector<double> > getVectorFromFile(std::string pathToTxt, const int nbColumn);

int main(int argc, char **argv)
{
	// Get the data inside a file
	std::string pathToData = "";
	std::vector<std::vector<double> > data = M_getVectorFromFile(pathToData, 3);
	
	const int nbPixel = data.size();
	
	for ()

	return 0;
}

std::vector<std::vector<double> > getVectorFromFile(std::string pathToTxt, const int nbColumn)
{
	// Returns a vector with the content of a file
	//
	// std::string pathToTxt                    : path to the file
	// const int nbCulumn                       : number of columns
	//
	// return std::vector<std::vector<double> > : vector with the info
	//											   vec[i] = column i
	
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
