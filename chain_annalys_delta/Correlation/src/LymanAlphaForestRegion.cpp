//===================================================================================
//
//         FILE: LymanAlphaForestRegion.cpp
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

#include "LymanAlphaForestRegion.h"
#include "../../../Constants/constants.h"

#include <fstream>
#include <iostream>	// std::cout
#include <sstream>	// stringstream
#include <cmath>
#include <limits>	// numeric_limits

// Nb of allowed tries in a 'while' loop before ending
long unsigned int nbAllowedTriesWhileLoop_ = 1000000;

LymanForest::LymanForest(std::vector< std::vector<double> > forests, unsigned int nbRegions, double raSeperationTwoRegions, bool euclidean/*=false*/) {
	
	euclidean_ = euclidean;
	if (euclidean_) raSeperationTwoRegions_ = 0.;
	else raSeperationTwoRegions_ = raSeperationTwoRegions;
	
	if (nbRegions < 1) {
		std::cout << "  LymanForest::LymanForest:: ERROR: nbRegions < 1" <<std::endl;
		return;
	}

	if (nbRegions>forests[0].size()) {
		nbRegion_ = forests[0].size();
		std::cout << "  LymanForest::LymanForest:: WARNING:  nbRegions>forests[0].size()" <<std::endl;
	}
	else nbRegion_ = nbRegions;

	std::cout << "\n  nbRegion_ = " << nbRegion_ << std::endl;
	std::cout << "  raSeperationTwoRegions_ = " << raSeperationTwoRegions_ << std::endl;
	std::cout << "  euclidean_ = " << euclidean_ << std::endl;

	ra_ = std::vector <std::vector <double> >(1);
	de_ = std::vector <std::vector <double> >(1);
	pa_ = std::vector <std::vector <double> >(1);
	id_ = std::vector <std::vector <unsigned int> >(1);
	
	co_ = std::vector <std::vector <double> >(1);
	co_[0].resize(4, 0.);
	nbPairs_   = 0.;
	nbForests_ = 0;

	// Load the list from 'forests'
	LoadForest(forests);
	
	// Devide the sky into 'nbRegion_' regions
	DevideInRegions();
	nbRegion_ = ra_.size();
	
	// Print the information on the screen
	PrintData();
}
LymanForest::LymanForest(std::string pathToLoad, unsigned int nbRegions, double raSeperationTwoRegions, bool euclidean/*=false*/) {

	euclidean_ = euclidean;
	if (euclidean_) raSeperationTwoRegions_ = 0.;
	else raSeperationTwoRegions_ = raSeperationTwoRegions;

	if (nbRegions < 1) {
		std::cout << "  ERROR: LymanForest::LymanForest: nbRegions < 1" <<std::endl;
		return;
	}
	nbRegion_ = nbRegions;

	std::cout << "\n  nbRegion_ = " << nbRegion_ << std::endl;
	std::cout << "  raSeperationTwoRegions_ = " << raSeperationTwoRegions_ << std::endl;
	std::cout << "  euclidean_ = " << euclidean_ << std::endl;

	ra_ = std::vector <std::vector <double> >(nbRegion_);
	de_ = std::vector <std::vector <double> >(nbRegion_);
	pa_ = std::vector <std::vector <double> >(nbRegion_);
	id_ = std::vector <std::vector <unsigned int> >(nbRegion_);

	co_ = std::vector <std::vector <double> >(nbRegion_, std::vector <double>(4,0.) );

	nbPairs_   = 0.;
	nbForests_ = 0;

	// Load the list from 'pathToLoad'
	LoadForestFromMap(pathToLoad);
}
LymanForest::~LymanForest() {
	ra_.clear();
	de_.clear();
	pa_.clear();
	id_.clear();
	co_.clear();
}

void LymanForest::LoadForest(std::vector< std::vector<double> > forests) {
	
	std::cout << "\n\n\n" <<  std::endl;
	std::cout << "***************************************************" <<  std::endl;
	
	const unsigned int nbForest = forests[0].size();

	for (unsigned int i=0; i<nbForest; i++) {

		double ra = forests[1][i];

		if (ra < M_PI/2. && !euclidean_) ra += 2.*M_PI;
		id_[0].push_back( forests[0][i] );
		ra_[0].push_back( ra );
		de_[0].push_back( forests[2][i] );
		pa_[0].push_back( forests[3][i] );
		
		if (nbForests_ == 0) {
			co_[0][0] = ra_[0][0];
			co_[0][1] = ra_[0][0];
			co_[0][2] = de_[0][0];
			co_[0][3] = de_[0][0];
		}
		else {
			co_[0][0] = std::min(co_[0][0], ra );
			co_[0][1] = std::max(co_[0][1], ra );
			co_[0][2] = std::min(co_[0][2], forests[2][i] );
			co_[0][3] = std::max(co_[0][3], forests[2][i] );
		}
		nbPairs_ += forests[3][i];
		nbForests_ ++;
	}

	if (nbRegion_ > nbForests_) {
		std::cout << "  ERROR: LymanForest::LoadForest: nbRegions > nbForest" <<std::endl;
		return;
	}
	
	minRa_ = co_[0][0];
	maxRa_ = co_[0][1];
	minDe_ = co_[0][2];
	maxDe_ = co_[0][3];
}
void LymanForest::LoadForestFromMap(std::string pathToLoad) {
	
	std::cout << "\n\n\n" <<  std::endl;
	std::cout << "***************************************************" <<  std::endl;
	std::cout << "  Read Lya list from : " << std::endl;
	std::cout << "  " << pathToLoad << std::endl;

	std::ifstream fileData(pathToLoad.c_str());
	while (fileData) {

		unsigned int id = 0;
		unsigned int re = 0;
		double ra = 0.;
		double de = 0.;
		double nbPairs = 0.;

		fileData>>id>>re>>ra>>de>>nbPairs;
		if (fileData==0) break;

		if (re>nbRegion_) std::cout << "  LymanForest::LoadForestFromMap::  ERROR:  re>nbRegion_, re = " << re << "  , nbRegion_ = " << nbRegion_ << std::endl;

		if (ra < M_PI/2. && !euclidean_) ra += 2.*M_PI;
		id_[re].push_back( id );
		ra_[re].push_back( ra );
		de_[re].push_back( de );
		pa_[re].push_back( nbPairs );

		if (nbForests_==0) {
			minRa_ = ra;
			maxRa_ = ra;
			minDe_ = de;
			maxDe_ = de;
		}
		else {
			minRa_ = std::min(minRa_, ra);
			maxRa_ = std::max(maxRa_, ra);
			minDe_ = std::min(minDe_, de);
			maxDe_ = std::max(maxDe_, de);
		}

		if (id_[re].size() == 0) {
			co_[re][0] = ra;
			co_[re][1] = ra;
			co_[re][2] = de;
			co_[re][3] = de;
		}
		else {
			co_[re][0] = std::min(co_[re][0], ra);
			co_[re][1] = std::max(co_[re][1], ra);
			co_[re][2] = std::min(co_[re][2], de);
			co_[re][3] = std::max(co_[re][3], de);
		}

		nbPairs_ += nbPairs;
		nbForests_++;
	}
	fileData.close();

	if (nbRegion_ > nbForests_) {
		std::cout << "  ERROR: LymanForest::LoadForest: nbRegions > nbForest" <<std::endl;
		return;
	}
}
bool LymanForest::DevideInRegions(void) {
	
	// Maximum number of steps we have to use to reach
	//  equally spaced strips
	const unsigned nbStepsDe = 10000;
	// Maximum number of steps we have to use to reach
	//  equally spaced squares
	const unsigned nbStepsRa = 10000;
	
	
	if (nbRegion_ == 1) return true;
	else {

		std::vector <std::vector <double> > tmp_ra(nbRegion_);
		std::vector <std::vector <double> > tmp_de(nbRegion_);
		std::vector <std::vector <double> > tmp_pa(nbRegion_);
		std::vector <std::vector <unsigned int> > tmp_id(nbRegion_);
		
		std::vector <std::vector <double> > tmp_co(nbRegion_);
		for (unsigned int i=0; i<nbRegion_; i++) {
			tmp_co[i].resize(4, 0.);
		}
		
		std::vector <std::vector <double> > tmp2_ra(2);
		std::vector <std::vector <double> > tmp2_de(2);
		std::vector <std::vector <double> > tmp2_pa(2);
		std::vector <std::vector <unsigned int> > tmp2_id(2);
		
		std::vector <std::vector <double> > tmp2_co(2);
		for (unsigned int i=0; i<2; i++) {
			tmp2_co[i].resize(4, 0.);
		}
		
		// Put the forest in their respective regions
		for (unsigned int i=0; i<ra_[0].size(); i++) {
		
			unsigned int regIdx = 0;
			if (ra_[0][i] < raSeperationTwoRegions_ && !euclidean_) regIdx = 1;
			
			tmp2_ra[regIdx].push_back(ra_[0][i]);
			tmp2_de[regIdx].push_back(de_[0][i]);
			tmp2_pa[regIdx].push_back(pa_[0][i]);
			tmp2_id[regIdx].push_back(id_[0][i]);
			
			if (i == 0) {
				tmp2_co[regIdx][0] = ra_[0][i];
				tmp2_co[regIdx][1] = ra_[0][i];
				tmp2_co[regIdx][2] = de_[0][i];
				tmp2_co[regIdx][3] = de_[0][i];
			}
			else {
				tmp2_co[regIdx][0] = std::min(tmp2_co[regIdx][0], ra_[0][i]);
				tmp2_co[regIdx][1] = std::max(tmp2_co[regIdx][1], ra_[0][i]);
				tmp2_co[regIdx][2] = std::min(tmp2_co[regIdx][2], de_[0][i]);
				tmp2_co[regIdx][3] = std::max(tmp2_co[regIdx][3], de_[0][i]);
			}
		}
		
		if (nbRegion_ == 2) {
			tmp_ra = tmp2_ra;
			tmp_de = tmp2_de;
			tmp_pa = tmp2_pa;
			tmp_id = tmp2_id;
			tmp_co = tmp2_co;
		}		
		else {			
			
			// Get the total number of pairs in the two regions
			double nbPairsInRegions[2] = {};
			for (unsigned int i=0; i<2; i++) {
				for (unsigned int j=0; j<tmp2_ra[i].size(); j++) {
					nbPairsInRegions[i] += tmp2_pa[i][j];
				}
			}
			std::cout << "  nb of pairs in region 0: " << nbPairsInRegions[0] << std::endl;
			std::cout << "  nb of pairs in region 1: " << nbPairsInRegions[1] << std::endl;

			// Find the best way to devide the two regions
			unsigned int nbRegionInRegions[2] = {};
			nbRegionInRegions[0] = nbRegion_*nbPairsInRegions[0]/(nbPairsInRegions[0]+nbPairsInRegions[1]);
			nbRegionInRegions[1] = nbRegion_-nbRegionInRegions[0];
			//std::cout << "  " << 1.*nbRegion_*nbPairsInRegions[0]/(nbPairsInRegions[0]+nbPairsInRegions[1]) << " " << 1.*nbRegion_*nbPairsInRegions[1]/(nbPairsInRegions[0]+nbPairsInRegions[1]) << std::endl;
			//std::cout << "  " << 1.*nbPairsInRegions[0]/(nbPairsInRegions[0]+nbPairsInRegions[1]) << " " << 1.*nbPairsInRegions[1]/(nbPairsInRegions[0]+nbPairsInRegions[1]) << std::endl;
			//std::cout << "  " << 1.*nbPairsInRegions[1]/nbPairsInRegions[0] << std::endl;
			std::cout << "  nb of sub-regions in region 0: " << nbRegionInRegions[0] << std::endl;
			std::cout << "  nb of sub-regions in region 1: " << nbRegionInRegions[1] << std::endl;
			
			// Get the size of the two regions
			double deltaRa[2] = {};
			double deltaDe[2] = {};
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				deltaRa[i] = tmp2_co[i][1]-tmp2_co[i][0];
				deltaDe[i] = tmp2_co[i][3]-tmp2_co[i][2];
			}
			
			// Get the number of strips in each region in order
			//  to have regions in approximetly squared shape
			unsigned int nbStrips[2] = {};
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				if (deltaDe[i]>0. && deltaRa[i]>0.) nbStrips[i] = sqrt(nbRegionInRegions[i]*deltaDe[i]/deltaRa[i]);
				else nbStrips[i] = 1;
			}
			std::cout << "  nb of strips per regions = " << nbStrips[0] << " " << nbStrips[1] << std::endl;
			
			// Get the number of regions per band
			unsigned int nbColumn[2] = {};
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				nbColumn[i] = nbRegionInRegions[i]/nbStrips[i];
			}
			std::cout << "  nb of columns per regions = " << nbColumn[0] << " " << nbColumn[1] << std::endl;
			
			// Get the spear regions to distribute over the sky
			unsigned int nbSpears[2] = {};
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				nbSpears[i] = nbRegionInRegions[i] - nbStrips[i]*nbColumn[i];
			}
			std::cout << "  nb of spears regions = " << nbSpears[0] << " " << nbSpears[1] << std::endl;
			
			// Get the number of regions in each strips
			std::vector<std::vector<unsigned int> > nbOfRegionsPerStrips(2);
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				nbOfRegionsPerStrips[i].resize(nbStrips[i]);
				for (unsigned int j=0; j<nbStrips[i]; j++) {
					const unsigned int idx = j; //nbStrips[i]-j-1;
					
					unsigned int nbSquares = nbColumn[i];
					if (nbSpears[i] != 0) {
						nbSquares ++;
						nbSpears[i] --;
					}
					
					nbOfRegionsPerStrips[i][idx] = nbSquares;
					std::cout << "  region =  " << i << " , strip = " << idx << " , nb of regions = " << nbOfRegionsPerStrips[i][idx] << std::endl;	
				}
			}

			// Get the number of pairs in each strips
			std::vector<std::vector<double> > nbOfPairsPerStrips(2);
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				nbOfPairsPerStrips[i].resize(nbStrips[i]);
				double nbTotPairs = 0;
				for (unsigned int j=0; j<nbStrips[i]; j++) {
					if (j == nbStrips[i]-1) nbOfPairsPerStrips[i][j] = nbPairsInRegions[i] - nbTotPairs;
					else nbOfPairsPerStrips[i][j] = nbOfRegionsPerStrips[i][j]*nbPairsInRegions[i]/nbRegionInRegions[i];
					nbTotPairs += nbOfRegionsPerStrips[i][j]*nbPairsInRegions[i]/nbRegionInRegions[i];
					std::cout << "  region =  " << i << " , strip = " << j << " , nb of pairs = " << nbOfPairsPerStrips[i][j] << std::endl;	
				}
			}


			// Find the intervals in declination of each strips
			std::vector<std::vector<std::vector<double> > > deIntervalsEachStrips(2);
			
			// Find also the true number of pairs in each strips
			std::vector<std::vector<double> > nbOfPairsPerStripsTrue(2);
			
			unsigned int stripIdx = 0;
			
			for (unsigned int i=0; i<2; i++) {

				if (nbRegionInRegions[i]==0) continue;
				deIntervalsEachStrips[i].resize(nbStrips[i]);
				nbOfPairsPerStripsTrue[i].resize(nbStrips[i]);

				double tmp_nbPairs = 0;
				const double deMin = tmp2_co[i][2];
				const double deMax = tmp2_co[i][3];
				const double sizeStrips = deltaDe[i]/nbStepsDe;

				for (unsigned int j=0; j<nbStrips[i]-1; j++) {

					deIntervalsEachStrips[i][j].resize(2);

					double deMinStrip = 0.;
					if (j == 0) {
						if (deMin < 0) deMinStrip = 1.1*deMin;
						else deMinStrip = 0.9*deMin;
					}
					else {
						deMinStrip = deIntervalsEachStrips[i][j-1][1];
					}

					// Find the best declination max for the strip
					double deMaxStrip = deMinStrip;
					double nbPairs = 0;
					bool neg = true;

					unsigned int nbLoop = 0;					
					while (neg && nbLoop<nbAllowedTriesWhileLoop_) {
						nbLoop++;

						deMaxStrip += sizeStrips;
						
						nbPairs = 0;
						for (unsigned int l=0; l<tmp2_ra[i].size(); l++) {
							if ( (tmp2_de[i][l] >= deMinStrip) && (tmp2_de[i][l] < deMaxStrip) ) {
								nbPairs += tmp2_pa[i][l];
							}
						}
						
						neg = (1.*nbPairs-1.*nbOfPairsPerStrips[i][j] < 0. );
						if (neg == false && stripIdx%2 == 0) deMaxStrip -= sizeStrips;
					}

					if (nbLoop>=nbAllowedTriesWhileLoop_) {
						std::cout << "  LymanForest::DevideInRegions::  ERROR:  nbLoop>=nbAllowedTriesWhileLoop_" << std::endl;
						return false;
					}
					
					// Find the final number of pairs
					nbPairs = 0;
					for (unsigned int l=0; l<tmp2_ra[i].size(); l++) {
						if ( (tmp2_de[i][l] >= deMinStrip) && (tmp2_de[i][l] < deMaxStrip) ) {
							nbPairs += tmp2_pa[i][l];
						}
					}
					
					//std::cout << i << " " << j << " " << nbPairs << " " << nbOfPairsPerStrips[i][j] << " " << (1.*nbPairs-1.*nbOfPairsPerStrips[i][j])/nbOfPairsPerStrips[i][j] << std::endl;
					nbOfPairsPerStripsTrue[i][j] = nbPairs;
					tmp_nbPairs += nbPairs;
					
					deIntervalsEachStrips[i][j][0] = deMinStrip;
					deIntervalsEachStrips[i][j][1] = deMaxStrip;

					stripIdx ++;
				}

				double deMinStrip = 0.;
				double deMaxStrip = 0.;

				if (nbStrips[i] > 1) deMinStrip = deIntervalsEachStrips[i][nbStrips[i]-2][1];
				if (deMax < 0) deMaxStrip = 0.9*deMax;
				else deMaxStrip = 1.1*deMax;

				deIntervalsEachStrips[i][nbStrips[i]-1].resize(2);
				deIntervalsEachStrips[i][nbStrips[i]-1][0] = deMinStrip;
				deIntervalsEachStrips[i][nbStrips[i]-1][1] = deMaxStrip;

				double nbPairs = nbPairsInRegions[i] - tmp_nbPairs;
				nbOfPairsPerStripsTrue[i][nbStrips[i]-1] = nbPairs;
				//std::cout << i << " " << nbStrips[i]-1 << " " << nbPairs << " " << nbOfPairsPerStrips[i][nbStrips[i]-1] << " " << (1.*nbPairs-1.*nbOfPairsPerStrips[i][nbStrips[i]-1])/nbOfPairsPerStrips[i][nbStrips[i]-1] << std::endl;
			}

			// Find the intervals in declination and right assencion
			//  of all regions
			unsigned int regIdx1 = 0;
			std::vector<std::vector<std::vector<std::vector<double> > > > intervalsEachRegions(2);
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				intervalsEachRegions[i].resize(nbStrips[i]);
				
				const double raMin = tmp2_co[i][0];
				const double raMax = tmp2_co[i][1];
				
				for (unsigned int j=0; j<nbStrips[i]; j++) {

					intervalsEachRegions[i][j].resize(nbOfRegionsPerStrips[i][j]);
					
					const double sizeColumn = deltaRa[i]/nbStepsRa;
					const double deMin = deIntervalsEachStrips[i][j][0];
					const double deMax = deIntervalsEachStrips[i][j][1];
					
					const double nbPairsPerRegions = 1.*nbOfPairsPerStripsTrue[i][j]/nbOfRegionsPerStrips[i][j];
					double tmp_nbPairs = 0;
					
					for (unsigned int k=0; k<nbOfRegionsPerStrips[i][j]-1; k++) {
						intervalsEachRegions[i][j][k].resize(4);
						
						double raMinSquare = 0.;
						if (k == 0) {
							if (raMin < 0) raMinSquare = 1.1*raMin;
							else raMinSquare = 0.9*raMin;
						}
						else {
							raMinSquare = intervalsEachRegions[i][j][k-1][1];
						}
						
						// Find the best right assencion max for the squares
						double raMaxSquare = raMinSquare;
						double nbPairs = 0;
						bool neg = true;

						unsigned int nbLoop = 0;
						while (neg && nbLoop<nbAllowedTriesWhileLoop_) {
							nbLoop++;
							raMaxSquare += sizeColumn;

							nbPairs = 0;
							for (unsigned int l=0; l<tmp2_ra[i].size(); l++) {
								if ( (tmp2_ra[i][l] >= raMinSquare) && (tmp2_ra[i][l] < raMaxSquare) && (tmp2_de[i][l] >= deMin) && (tmp2_de[i][l] < deMax) ) {
									nbPairs += tmp2_pa[i][l];
								}
							}
							
							neg = (1.*nbPairs-1.*nbPairsPerRegions < 0. );
							if (neg == false && regIdx1%2 == 0) raMaxSquare -= sizeColumn;
						}

						if (nbLoop>=nbAllowedTriesWhileLoop_) {
							std::cout << "  LymanForest::DevideInRegions::  ERROR:  nbLoop>=nbAllowedTriesWhileLoop_" << std::endl;
							return false;
						}
						
						nbPairs = 0;
						for (unsigned int l=0; l<tmp2_ra[i].size(); l++) {
							if ( (tmp2_ra[i][l] >= raMinSquare) && (tmp2_ra[i][l] < raMaxSquare) && (tmp2_de[i][l] >= deMin) && (tmp2_de[i][l] < deMax) ) {
								nbPairs += tmp2_pa[i][l];
							}
						}
						
						//std::cout << i << " " << j << " " << nbPairs << " " << nbPairsPerRegions << " " << (1.*nbPairs-1.*nbPairsPerRegions)/nbPairsPerRegions << std::endl;
						tmp_nbPairs += nbPairs;
						
						intervalsEachRegions[i][j][k][0] = raMinSquare;
						intervalsEachRegions[i][j][k][1] = raMaxSquare;
						intervalsEachRegions[i][j][k][2] = deMin;
						intervalsEachRegions[i][j][k][3] = deMax;
						
						regIdx1++;
					}
					
					//long long unsigned int nbPairs = nbOfPairsPerStripsTrue[i][j]-tmp_nbPairs;
					//std::cout << i << " " << j << " " << nbPairs << " " << nbPairsPerRegions << " " << (1.*nbPairs-1.*nbPairsPerRegions)/nbPairsPerRegions << std::endl;
					

					const double raMinSquare = intervalsEachRegions[i][j][nbOfRegionsPerStrips[i][j]-2][1];
					double raMaxSquare = 0.;
					if (raMax < 0) raMaxSquare = 0.9*raMax;
					else raMaxSquare = 1.1*raMax;
					
					intervalsEachRegions[i][j][nbOfRegionsPerStrips[i][j]-1].resize(4);
					intervalsEachRegions[i][j][nbOfRegionsPerStrips[i][j]-1][0] = raMinSquare;
					intervalsEachRegions[i][j][nbOfRegionsPerStrips[i][j]-1][1] = raMaxSquare;
					intervalsEachRegions[i][j][nbOfRegionsPerStrips[i][j]-1][2] = deMin;
					intervalsEachRegions[i][j][nbOfRegionsPerStrips[i][j]-1][3] = deMax;
				}
			}
			
			// Fill the vectors of each of the 'nbRegion_' regions
			unsigned int regIdx = 0;
			for (unsigned int i=0; i<2; i++) {
				if (nbRegionInRegions[i]==0) continue;
				for (unsigned int j=0; j<nbStrips[i]; j++) {
					for (unsigned int k=0; k<nbOfRegionsPerStrips[i][j]; k++) {
						for (unsigned int l=0; l<tmp2_ra[i].size(); l++) {
							if ( (tmp2_ra[i][l] >= intervalsEachRegions[i][j][k][0]) && (tmp2_ra[i][l] < intervalsEachRegions[i][j][k][1]) && (tmp2_de[i][l] >= intervalsEachRegions[i][j][k][2]) && (tmp2_de[i][l] < intervalsEachRegions[i][j][k][3]) ) {
								tmp_ra[regIdx].push_back(tmp2_ra[i][l]);
								tmp_de[regIdx].push_back(tmp2_de[i][l]);
								tmp_pa[regIdx].push_back(tmp2_pa[i][l]);
								tmp_id[regIdx].push_back(tmp2_id[i][l]);
							}
						}
						regIdx++;
					}
				}		
			}
			
			// Fill the 'tmp_co' vector.
			for (unsigned int i=0; i<nbRegion_; i++) {
				for (unsigned int j=0; j<tmp_ra[i].size(); j++) {
					if (j == 0) {
						tmp_co[i][0] = tmp_ra[i][0];
						tmp_co[i][1] = tmp_ra[i][0];
						tmp_co[i][2] = tmp_de[i][0];
						tmp_co[i][3] = tmp_de[i][0];
					}
					else {
						tmp_co[i][0] = std::min(tmp_co[i][0], tmp_ra[i][j] );
						tmp_co[i][1] = std::max(tmp_co[i][1], tmp_ra[i][j] );
						tmp_co[i][2] = std::min(tmp_co[i][2], tmp_de[i][j] );
						tmp_co[i][3] = std::max(tmp_co[i][3], tmp_de[i][j] );
					}
				}
			}
		}
		
		ra_ = tmp_ra;
		de_ = tmp_de;
		pa_ = tmp_pa;
		id_ = tmp_id;
		co_ = tmp_co;
		
	}
	return true;
}
void LymanForest::PrintData(void) {

	unsigned int nbForesTot = 0;
	double nbPairsTot = 0;

	for (unsigned int i=0; i<nbRegion_; i++) {
		
		nbForesTot += ra_[i].size();
		for (unsigned int j=0; j<ra_[i].size(); j++) {
			nbPairsTot += pa_[i][j];
		}
	}
	
	if (nbForests_ != nbForesTot) std::cout << " ERROR: nbForests_ != nbForesTot " << nbForests_ << " " << nbForesTot << " " << nbForesTot-nbForests_ << std::endl;
	if (nbPairs_ != nbPairsTot)  std::cout  << " ERROR: nbPairs_   != nbPairsTot " << nbPairs_   << " " << nbPairsTot << " " << nbPairsTot-nbPairs_   << std::endl;
	
	// Get the mean number of pairs per regions
	const double nbPairsPerRegions = 1.*nbPairsTot/nbRegion_;

	std::cout << "***************************************************" <<  std::endl;
	std::cout << "  Total nb of forest = " << nbForesTot << std::endl;
	std::cout << "  Total nb of pairs pixel-qso  = " << nbPairsTot << std::endl;
	std::cout << std::endl;
	
	for (unsigned int i=0; i<nbRegion_; i++) {
		
		double nbPairs = 0;
		for (unsigned int j=0; j<ra_[i].size(); j++) {
			nbPairs += pa_[i][j];
		}
		std::cout << "  Region N°" << i << " , nb forest = " << ra_[i].size() << " , nb pairs pixel-qso = " << nbPairs << " , percent with respect to mean = " << 100.*(nbPairs-nbPairsPerRegions)/nbPairsPerRegions << std::endl;
	}
}


void LymanForest::GetCoordRegion(int regionIdx, double* array) {
	
	int idx = regionIdx;
	if (idx == -1) idx = 0;
	
	for (unsigned int i=0; i<ra_[idx].size(); i++) {
		
		double ra = ra_[idx][i];
		if (ra >= 2.*M_PI) ra -= 2.*M_PI;
		const double de = de_[idx][i];
	
		if (i == 0) {
				array[0] = ra;
				array[1] = ra;
				array[2] = de;
				array[3] = de;
			}
		else {
			array[0] = std::min(array[0], ra);
			array[1] = std::max(array[1], ra);
			array[2] = std::min(array[2], de);
			array[3] = std::max(array[3], de);
		}
	}
	
}
void LymanForest::PrintRegionDetail(int regionIdx) {

	if (regionIdx<0 || (unsigned int)(regionIdx)>=nbRegion_) std::cout << "  LymanForest::PrintRegionDetail:: ERROR: regionIdx<0 || regionIdx>=nbRegion_" << std::endl;

	double nbPairs = 0;	
	for (unsigned int i=0; i<ra_[regionIdx].size(); i++) {
		nbPairs += pa_[regionIdx][i];
	}
	std::cout << "  Region N°" << regionIdx << " , nb forest = " << ra_[regionIdx].size() << " , nb pairs pixel-qso = " << nbPairs << std::endl;

	return;
}
void LymanForest::GetRegionArray(unsigned int* array) {
	
	for (unsigned int i=0; i<nbRegion_; i++) {
		for (unsigned int j=0; j<ra_[i].size(); j++) {
			array[ id_[i][j] ] = i;
		}
	}
	
}
void LymanForest::SaveRegionMap(std::string pathToSave) {
	
	std::cout << "\n  " << pathToSave << std::endl;
	
	std::ofstream fFile;
	fFile.open(pathToSave.c_str());
	fFile << std::scientific;
	fFile.precision(std::numeric_limits<double>::digits10);

	//const double piTimes2   = M_PI*2.;
	double radToDeg = 1.;
	if (!euclidean_) radToDeg = C_RADTODEG;

	for (unsigned int i=0; i<nbRegion_; i++) {
		const unsigned int nbForest = ra_[i].size();
		for (unsigned int j=0; j<nbForest; j++) {

			//double ra = ra_[i][j];
			//if (ra > piTimes2) ra -= piTimes2;
			//ra *= radToDeg;

			fFile << id_[i][j];
			fFile << " " << i;
			fFile << " " << ra_[i][j]*radToDeg;
			fFile << " " << de_[i][j]*radToDeg;
			fFile << " " << pa_[i][j];
			fFile << std::endl;
		}
	}
	fFile.close();

	std::cout << "\n  LymanForest::SaveRegionMap: Finished" << std::endl;
}





















