//===================================================================================
//
//         FILE: tools.h
//
//        USAGE: ---
//
//  DESCRIPTION: 
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

#ifndef TOOLS_H
#define TOOLS_H

class Tools
{
	private:

		void get_Ra_Dec_DLA(void);
		void get_Catalogue(void);
		void get_flux_vs_lambdaObs(void);
		void get_weigted_covar_matrix(void);

	public:
		Tools(int argc, char** argv);
		~Tools();
};


#endif


