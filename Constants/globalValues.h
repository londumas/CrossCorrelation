//===================================================================================
//
//         FILE: globalValues.h
//
//        USAGE:
//
//  DESCRIPTION: global values that can change
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

#include <vector>

#ifndef GLOBALVALUES_H
#define GLOBALVALUES_H

// ---------------------------------------------------------------------
//
//		Forests
//
// ---------------------------------------------------------------------

//// const double C_AUTOCORRCRIT = 0.0005*M_PI/180.;
#define C_AUTOCORRCRIT 8.72664625997164823e-06
#define C_AUTOCORRCRITDEG 0.0005
//// The sky of DR12 is separeted into "North Galactic Cap" (NGC) and 
////	"South Galactic Cap" (SGC). This "global value" gives the value
////	of RA that separated the two regions
#define C_RA_SEPERATION_NGC_SGC 5.

//// Define the limit for which a FLUX_DLA is considered zero
#define C_FLUX_DLA_IS_ZERO 0.000000000000001
//// Define the max of DLA correction
#define C_DLACORR  0.8
//// Number of sub-samples the sky is split into:
#define C_NBSUBSAMPLES 80
//// Minimum number of pixels in a forest
#define C_MIN_NB_PIXEL 50
//// Conversion from mean_flux_in_forest to alpha
#define C_CONVERT_FROM_FLUX_TO_ALPHA 1.154

//// Lines in the sky to veto or not
const unsigned int nbVetoLines__ = 21;
const double vetoLine__[42] = {3615.,3619.,3932.,3937.,3966.,3972.,4042.,4050.,
        4357.,4362.,5458.,5467.,5573.,5585.,5682.,5695.,
        5885.,5902.,6235.,6241.,6256.,6263.,6296.,6311.,
        6320.,6334.,6362.,6369.,6498.,6502.,6554.,6557.,
        6825.,6840.,6862.,6870.,6922.,6928.,6948.,6954.,
        6977.,6982.};

//// List of absorbers in IGM:
const double absorber__[29] = {
2852.96, 2803.5324, 2796.3511, 2600.1724835,  2586.6495659,  2576.877, 2382.7641781,  2374.4603294,  2344.2129601,
1862.79113,  1854.71829, 1670.7886, 1608.4511, 1550.77845, 1548.2049, 1526.70698,
1402.77291, 1393.76018, 1335., 1304.3702,   1302.1685,  1260.4221,  1242.804, 1238.821,
1215.67 ,1206.500, 1193.2897, 1190.4158,
1025.7223
};
const std::string absorberName__[29] = {
"MgI",   "MgII_a",  "MgII_b",  "FeII_a",      "FeII_b",      "MnII",   "FeII_c",      "FeII_d",      "FeII_e",
"AlIII(1863)" ,  "AlIII(1855)" , "AlII(1671)" ,   "FeII(1609)",  "CIV(1551)" ,   "CIV(1548)",   "SiII(1527)",
"SiIV(1403)",   "SiIV(1394)",   "CII(1335)" ,"SiII(1304)",    "OI(1302)",     "SiII(1260)",   "NV(1243)",   "NV(1239)",
"LYA"   ,"SiIII(1207)",  "SiII(1193)",  "SiII(1190)",
"LYB"
};

// ---------------------------------------------------------------------
//
//		Cosmology
//
// ---------------------------------------------------------------------

//// Defines the number of bins for the interpolation dist=f(z)
#define C_NBBINREDSH 12000
//// Defines the extrema for the interpolation dist=f(z)
#define C_ZEXTREMABINCONVERT0 0.
#define C_ZEXTREMABINCONVERT1 12.

#define C_H 0.70 // (100 km/s/Mpc)^-1
#define C_OMEGAM 0.27
#define C_OMEGAB 0.0463



#endif


















/*

// const double C_DEGTORAD = M_PI/180;
#define C_DEGTORAD 1.74532925199432955e-02

// Number of regions the statistics is split into:
//  has to seen inside qso and Lya files. = 8
#define C_NBREGION 8

// Number of bootstraps the sky is split into :
//  can be for now '1' or '2' or '8' or '81' or '139'
//  can be for simulations : 78
#define C_NBBOOTSTRAP 80

// Defines if you are debuging or not
#define C_DEBUG false
#define C_DEBUGNUMBER 10000 //

#define C_REOBSERV false

// Which method to use : 1 or 2
#define C_METHOD 2

// Either to look at one plat or all
#define C_PLATEBOOL false

// What type of data
#define C_DATASET "DR12"
//#define C_DATASET "eBOSS"
//#define C_DATASET "Mock__DR11"
#define C_DATASET2 ""

#define C_NAMEOFSAVING "eBOSS"


#define C_LISTPATH       "" // /home/usr201/mnt/hdumasde/Data/List/QSO_DR12_DR7.txt"
//#define C_LISTPATH       "/home/usr201/mnt/hdumasde/Data/List/QSO_DR12_DR7_eBOSS.txt"
#define C_LYALIST        "" ///home/usr201/mnt/hdumasde/Data/myLists/lyaMap.txt"
#define C_PLATESLISTPATH "" ///home/usr201/mnt/hdumasde/Data/myLists/platesProdV7.txt"
#define C_SAVEFOLDERPATH "/home/gpfs/manip/mnt/bao/hdumasde/Results/RootFile"
#define C_PRODNUMBER     "/home/gpfs/manip/mnt/bao/hdumasde/CrossCorrelation/Prod_DR12_1"
#define C_NBOFDATAFILES  60

#define C_PATHTOWEIGHT   "/home/gpfs/manip/mnt0607/bao/hdumasde/CrossCorrelation_StartingAgainFrom1347/CrossCorrelation/RootFile/source_histos_DR12.root"
#define C_NAMEOFSOURCEHISTOS ""

#define C_MYCPPLIB       "/home/usr201/mnt/hdumasde/Code/Cpp/Library/mathFunctions.h"
#define C_MYROOTLIB      "/home/usr201/mnt/hdumasde/Code/Root/Library/RootHistoFunctions.h"
#define C_MYCOSMOLIB     "/home/usr201/mnt/hdumasde/Code/Cpp/X_Correlat_firstStep/src/Cosmology.h"
#define C_MYQUASARLIB    "/home/usr201/mnt/hdumasde/Code/Cpp/X_Correlat_firstStep/src/QuasarList.h"
#define C_MYSPECTRALIB   "/home/usr201/mnt/hdumasde/Code/Cpp/Fast_computation/X_Correlat_firstStep/src/TTreeSpectra_Method2.h"
#define C_MYLYALIB       "/home/usr201/mnt/hdumasde/Code/Cpp/X_Correlat_secondStep_notCommonMeanDelta/src/LymanAlphaForestRegion.h"
#define C_NBOFSOURCES    1
#define C_PRODNUMBER2     ""
#define C_NBOFDATAFILES2  100



// Defines the criterion [ra, dec] if two objects are in the same
// direction
// const double C_AUTOCORRCRIT = C_DEGTORAD*0.0005;
#define C_AUTOCORRCRIT 8.72664625997164823e-06

// Defines the size of an array with the number of lya that have been cutted
// cut + 1, because SelCut[0]=total
#define C_NBCUT 21
// cutParam={Meanz, SNRmin, SNRmax, MeanFluxLambdaMin, MeanFluxLambdaMax,
//			MaxErr, MeanErr, MeanRmin, MeanRmax, iBAL, iDLA, nbMinPixel, chi2Method2_Min, chi2Method2_Max, aMethod2_Min, aMethod2_Max, fabs(bMethod2)_Max, ra=dec=0., fluxRapport1225Min=1.}
const double C_CUTPARAM[C_NBCUT] = {1.96, .1, 20., 2., 0., 10., 0., 150., 20., 1., -1., 50., 0., 5000., 0., 3., .1, 0., -100.};

// Extrema of the cross-correlation (in Mpc.h^{-1})
#define C_DISTEXTREMA0 0.
#define C_DISTEXTREMA1 160. //160.
// C_DISTMAX = sqrt(2)*C_DISTEXTREMA1
#define C_DISTMAX 227. //(for 160)  //72 //290 (for 200)
// Number of bin of the cross-correlation
#define C_NBBINCORRELAT 160 //50
#define C_BINSIZEDIST1D 1.
// Number of bin for the 2D correlation in the X axis
//#define C_NBBINCORRELAT2DX 20
// Number of bin for the 2D correlation in the Y axis
//#define C_NBBINCORRELAT2DY 40
// Number of bin for theta
//#define C_NBBINTHETA 20
// Number of bin for \mu^{2}
//#define C_NBBINMU 20
// Size of a bin of the correlation
//const double C_BINSIZEDIST1D = (C_DISTEXTREMA1-C_DISTEXTREMA0)/C_NBBINCORRELAT;
//const double C_BINSIZEDIST2DX = (C_DISTEXTREMA1 - 0.)/C_NBBINCORRELAT2DX;
//const double C_BINSIZEDIST2DY = 2.*C_DISTEXTREMA1/C_NBBINCORRELAT2DY;
//const double C_BINSIZETHETA2D = M_PI/C_NBBINTHETA;
//const double C_BINSIZEMU = (1.-(-1))/C_NBBINMU;
//#define C_BINSIZEDIST2DX 10.
//#define C_BINSIZEDIST2DY 10.
//#define C_BINSIZETHETA2D 1.57079632679489656e-01
//#define C_BINSIZEMU 0.1

#define C_NBPASS 4
// Min and max of deltas
#define C_DELTAEXTREMA0 -300.
#define C_DELTAEXTREMA1 300.
// Min and max of deltaErrSquare
#define C_DSQUAREEXTREMA0 1.e-5 //0.
#define C_DSQUAREEXTREMA1 10000. //0.8
// Starting values for the fit
const double C_PARSETVALUE[6] = {6.84738e-02, -1.05170, 1.32220e+03, 9e-3, 2.62123, 5.51456e-01};
// Number of bins for the error plots
#define C_NBBINPERROR 100
// Coeff of the kminus correction (3.8/2)
#define C_HALFGAMAWEIGHT 1.9
// Create the array for drowing delta=f(lambda)
#define C_NBBINDELTA 700
// Extrema for the wavelength
#define C_LAMBDAEXTREMA0 3600.
#define C_LAMBDAEXTREMA1 7100.
//const double C_BINSIZELAMBDA = (C_LAMBDAEXTREMA1 - C_LAMBDAEXTREMA0)/C_NBBINDELTA;
#define C_BINSIZELAMBDA 5.
#define C_LAMBDARFEXTREMA0 1040.
#define C_LAMBDARFEXTREMA1 1200.

// Define the range of redshift
#define C_Z0 2.25
#define C_ZEXTREMA0 2.
#define C_ZEXTREMA1 4.5
// Number of bins of slice of redshift
#define C_NBBINZ 12
// const double C_BINSIZEZ = (C_ZEXTREMA1 - C_ZEXTREMA0)/C_NBBINZ;
#define C_BINSIZEZ 2.08333333333333343e-01
*/







