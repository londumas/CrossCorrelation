//===================================================================================
//
//         FILE: main.cpp
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
#include "TRint.h"

int main(int argc, char** argv) {	
	// Needed to see the plots when running inside a program
	//TRint theApp("App", &argc, argv);

	Tools* xCorr = new Tools(argc, argv);
	delete xCorr;
	
	// Needed to see the plots when running inside a program
	//theApp.Run();

	return 0;
}







