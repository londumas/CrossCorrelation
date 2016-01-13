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

#ifndef GETDELTA_H
#define GETDELTA_H

class GetDelta
{
	private:

		void get_Ra_Dec_DLA(void);
		void get_Catalogue(void);
		void get_flux_vs_lambdaObs(void);
		void get_weigted_covar_matrix(void);

	public:
		GetDelta(int argc, char** argv);
		~GetDelta();
};


#endif


