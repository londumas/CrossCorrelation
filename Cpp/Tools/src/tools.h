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
		void get_template_function_redshift(void);

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


/*
/// If LYB
const std::string forest__  = "LYB";
const double lambdaRFLine__  = 1025.72;
const double lambdaRFMin__   = 800.;
const double lambdaRFMax__   = 1020.;
const unsigned int nbBinRFMax__ = 1085;
const double alphaStart__ = 1.;
#define C_NBBINZ 12
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5
*/

/// If LYA
const std::string forest__  = "LYA";
const double lambdaRFLine__  = 1215.67;
const double lambdaRFMin__   = 1040.;
const double lambdaRFMax__   = 1200.;
const unsigned int nbBinRFMax__   = 647;
const double alphaStart__ = 1.3;
#define C_NBBINZ 12
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5

/*
/// If SIIV
const std::string forest__  = "SIIV";
const double lambdaRFLine__  = 1393.76018;
const double lambdaRFMin__   = 1286.;
const double lambdaRFMax__   = 1380.;
const unsigned int nbBinRFMax__ = 326;
const double alphaStart__ = 1.;
#define C_NBBINZ 9
#define C_ZEXTREMA0 1.4
#define C_ZEXTREMA1 3.4
*/
/*
/// If CIV
const std::string forest__   = "CIV";
const double lambdaRFLine__  = 1548.2049;
const double lambdaRFMin__   = 1410.;
const double lambdaRFMax__   = 1530.;
const unsigned int nbBinRFMax__   = 373;
const double alphaStart__ = 1.;
#define C_NBBINZ 9
#define C_ZEXTREMA0 1.4
#define C_ZEXTREMA1 3.2
*/

/*
/// If MGII
const std::string forest__  = "MGII";
const double lambdaRFLine__  = 2804.;
const double lambdaRFMin__   = 1970.;
const double lambdaRFMax__   = 2760.;
const unsigned int nbBinRFMax__ = 1450;
const double alphaStart__ = 1.;
#define C_NBBINZ 9
#define C_ZEXTREMA0 0.4
#define C_ZEXTREMA1 1.7
*/




















#endif


