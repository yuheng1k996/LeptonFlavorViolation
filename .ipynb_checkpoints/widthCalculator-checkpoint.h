#ifndef _WIDTHCALCULATOR
#define _WIDTHCALCULATOR


#include <math.h>
#include <iostream>

namespace widthCalculator {

//physical constants
	constexpr double pi   = 3.1415926;
	constexpr double stw2 = 0.231224;
	constexpr double stw  = sqrt(stw2);
	constexpr double ctw  = sqrt(1-stw2);
	constexpr double ttw  = stw/ctw;
	constexpr double GF   = 1.16637876e-5;
	constexpr double MW   = 80.37912;
	constexpr double g2   = sqrt(GF*8*pow(MW,2)/sqrt(2.));

	constexpr double glL  = g2/sqrt(2)*ttw;
	constexpr double gnuL = glL;
	constexpr double guL  = -g2/(3*sqrt(2.))*ttw;
	constexpr double gdL  = guL;
	constexpr double gdR  = -2*g2/(3*sqrt(2.))*ttw;

	constexpr double GSnu = (gnuL - 0.5*gdL - 0.5*gdR);
	constexpr double GTnu = (0.25*gdL + 0.25*gdR);
	constexpr double GSL  = (0.5*guL + 0.5*gdR - glL);
	constexpr double GTL  = (0.25*guL+0.25*gdR);




//particle constants
	//mtau
	constexpr double mtau  = 1.77682;
	constexpr double tauSMGamma = 2.26578e-12;//GeV


	//pi
	constexpr double mpip  = 0.13957061;
	constexpr double mpi0  = 0.1349770;
	constexpr double I     = 1;
	constexpr double fpip  = 0.130;
	constexpr double mud   = 3.39E-3;
	constexpr double mu    = 2.3210e-3;
	constexpr double md    = 4.719e-3;
	constexpr double fSpip = I*pow(mpip,2)/(2*mud)*fpip;
	constexpr double fSpi0 = I*pow(mpi0,2)/(2*md)*(fpip/sqrt(2.));
	
	//kaon
	constexpr double ms = 0.0929;
	constexpr double fk = 0.1557;
	constexpr double fkV = 0.230;


	constexpr double mk0 = 0.49761;
	constexpr double fSk0 = I*pow(mk0,2)/(ms+md)*fk;

	constexpr double mkstar0 = 0.89594;
	constexpr double fTkstar0 = fkV;


	constexpr double mkp = 0.493677;
	constexpr double fSkp = I*pow(mkp,2)/(ms+mu)*fk;

	constexpr double mkstarp = 0.89555;
	constexpr double fTkstarp = fkV;





//functions

    // kinematic function
    const double lambdaHalf(double x, double y, double z);

    //tau decay widths
    double chiralTauWidth(double mNeutralino, double mtau, double mMeson, double fSMeson,double lambOverMSq);    
    double tensorTauWidth(double mNeutralino, double mtau, double mMeson, double fTMeson, double lambOverMSq);

    //neutralino decay widths
    double chiralNeutralinoWidth(double mNeutralino, double mMeson, double fSMeson,double lambOverMSq);
    double tensorNeutralinoWidth(double mNeutralino, double mMeson, double fTMeson, double lambOverMSq);



}
#endif
