//===================================================================================
//
//         FILE: CrossCorrelation.h
//
//        USAGE: ---
//
//  DESCRIPTION: Do the cross correlation between Lya forest and Qso
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

#include "Constants.h"

#include <fstream>
#include <vector>
#include "TH1D.h"
#include "TGraphErrors.h"


#ifndef GETDELTA_H
#define GETDELTA_H



class GetDelta
{
	private:

		void initFitCont(void);
		void get_flux_vs_lambdaObs(void);

	public:
		GetDelta(int argc, char** argv);
		~GetDelta();
};


#endif


