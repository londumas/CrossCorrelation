//===================================================================================
//
//         FILE: LymanAlphaForestRegion.h
//
//        USAGE: ---
//
//  DESCRIPTION: Class of list of forest to do regions
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

//#include "/home/helion/Documents/Thèse/Code/Cpp/X_Correlat_firstStep/src/Constants.h"
#include "/home/usr201/mnt/hdumasde/Code/Cpp/X_Correlat_firstStep/src/Constants.h"

#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"

#include <string>
#include <vector>

#ifndef LYMANALPHAFORESTREGION_H
#define LYMANALPHAFORESTREGION_H

class LymanForest
{

	private:
	
		bool mockjMc_;
		double raSeperationTwoRegions_;
	
		unsigned int nbRegion_;
		
		std::vector <std::vector <double> > ra_;
		std::vector <std::vector <double> > de_;
		std::vector <std::vector <double> > pa_;
		std::vector <std::vector <unsigned int> > id_;
		std::vector <std::vector <double> > co_;
		
		double minRa_;
		double maxRa_;
		double minDe_;
		double maxDe_;
		
		unsigned int nbForests_;
		double nbPairs_;
		
		void LoadForest(std::string pathToFile);
		void DevideInRegions(void);
		void PrintData(void);

	public:
		LymanForest(std::string pathToFile, unsigned int nbRegions);
		void GetCoordRegion(int regionIdx, double* array);
		void GetRegionArray(unsigned int* array);
		void SaveRegionMap(void);
};

#endif


