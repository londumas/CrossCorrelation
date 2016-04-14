//===================================================================================
//
//         FILE: main.cpp
//
//        USAGE: inside '/src' do 'make'
//               then inside '/src' do
//                   < time ../bin/main >> ../../../../Results/TxtFiles/output.txt& >
//               On iclust do before :
//                'tcsh'
//                'setenv CFITSIO /home/nfs/manip/mnt/bao/EXTLibs/cfitsio'
//                              
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
#include "TRint.h"

int main(int argc, char** argv)
{	
	// Needed to see the plots when running inside a program
	TRint theApp("App", &argc, argv);

	GetDelta* xCorr = new GetDelta(argc, argv);
	delete xCorr;
	
	// Needed to see the plots when running inside a program
	theApp.Run();

	return 0;
}







