//===================================================================================
//
//         FILE: LymanAlphaForestRegion.h
//
//        USAGE: ---
//
//  DESCRIPTION: Devide a sky of forests into 'nbRegion_' regions according to 
//		the statistical weight of this forets (i.e. the sum of weight over 
//		pairs forest-QSO)
//		The sky can be a RA-Dec (data) or X-Y (mocks)
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

#include <string>
#include <vector>

#ifndef LYMANALPHAFORESTREGION_H
#define LYMANALPHAFORESTREGION_H

class LymanForest
{

	private:

		// Flag: 
		//	* true if using a X-Y sky (euclidean coordinate)
		//	* false if using a RA-Dec sky (spherical coordinate)
		bool euclidean_;
		// The sky of DR12 is separeted into "North Galactic Cap" (NGC) and 
		//	"South Galactic Cap" (SGC). This double gives the value of RA
		//	that separated the two regions
		double raSeperationTwoRegions_;
		// Nb of regions to devide the sky into
		unsigned int nbRegion_;
		
		// In the following vector, the first index gives a region,
		//	the second index a forest
		// RA of the forest
		std::vector <std::vector <double> > ra_;
		// Dec of the forest
		std::vector <std::vector <double> > de_;
		// Sum of weight of the forest
		std::vector <std::vector <double> > pa_;
		// ID of the forest: (its index in the data file)
		std::vector <std::vector <unsigned int> > id_;
		
		// Coordonate: first index gives a region, second index gives:
		//	minRA, maxRA, minDec, maxDec
		std::vector <std::vector <double> > co_;
		
		// Attributes of the sky
		double minRa_;
		double maxRa_;
		double minDe_;
		double maxDe_;
		unsigned int nbForests_;
		double nbPairs_;
		
		// Load forest from a given Ascii file:
		//		- std::string pathToFile: path to file with the following
		// structure: idx  ra  de  nbPairs			
		void LoadForest(std::vector< std::vector<double> > forests);

		void LoadForestFromMap(std::string pathToLoad);
		
		// Performs the devision of the sky into regions
		bool DevideInRegions(void);
		
		// Print data
		void PrintData(void);

	public:
	
		// Class constructor:
		//		- std::string pathToFile: path to forest attributes
		//		- unsigned int nbRegions: number of regions to devide the sky into
		//		- This double gives the value of RA that separates NGC and SGC
		//		 (see private variable definition)
		//		- bool euclidean=false: true if eucldeen coordinate, false if spherical
		LymanForest(std::vector< std::vector<double> > forests, unsigned int nbRegions, double raSeperationTwoRegions, bool euclidean=false);


		LymanForest(std::string pathToLoad, unsigned int nbRegions, double raSeperationTwoRegions, bool euclidean=false);

		// Class destructor:
		~LymanForest();

		// Fill an array of size 4 with the coordinate of the region:
		//	minRA, maxRA, minDec, maxDec
		//		- int regionIdx: region to get the coordinate of
		//		- double* array: array of size 4 with: minRA, maxRA, minDec, maxDec
		void GetCoordRegion(int regionIdx, double* array);

		//
		void PrintRegionDetail(int regionIdx);
		
		// Fill an array of size superior to the number of forest in the
		//	all survey with the list of region's index of forests
		//		- unsigned int* array: array of size min the number
		//		of forest in the all survey, with the index of the 
		//		region correcponding to each forest
		void GetRegionArray(std::vector <unsigned int >& array);
		
		// Save into an Ascii file the list of forest with its
		//	corresponding region, file structure is:
		//	fFile << indexRegion  RA  Dec  nbPairs
		//		- std::string pathToSave: path where to save
		void SaveRegionMap(std::string pathToSave);
};

#endif











