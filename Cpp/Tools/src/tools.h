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
		void look_DLA(void);
		void get_delta_nicolas(void);
		void get_Catalogue(void);
		void get_flux_vs_lambdaObs(void);
		void get_weigted_covar_matrix(void);
		double VoigtProfile(double nhi, double lamb, double z_abs);

	public:
		Tools(int argc, char** argv);
		~Tools();
};


const double lambdaRFLine__  = 1215.67;
const double lambdaRFMin__   = 1040.;
const double lambdaRFMax__   = 1200.;
const unsigned int nbBinRFMax__   = 647;
const double alphaStart__ = 1.3;
#define C_NBBINZ 12
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5
const unsigned int NbDLA0__ = 15;
const double lambdaObsMin__   = 3600.; //3547.; //
const double lambdaObsMax__   = 7235.; //10326.; //
const double minAlpha__   = -100.;
const double maxAlpha__   = 100.;
const double betaStart__ = 0.;
const double minBeta__  = -0.6;
const double maxBeta__  =  0.6;

#endif


