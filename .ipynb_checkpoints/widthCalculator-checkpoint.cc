#include "widthCalculator.h"



namespace widthCalculator {

    // kinematic function
    const double lambdaHalf(double x, double y, double z) {
        double res2 = x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*y*z;
        if (res2 < 0) {
            std::cerr << "Warning: Your parameters demand for a decay which is kinematically forbidden!" << std::endl;
            return 0;
        }
        return sqrt(res2);
    }


//Tau decay widths:
    double chiralTauWidth(double mNeutralino, double mtau, double mMeson, double fSMeson,double lambOverMSq)
    {        
        if(mtau - mMeson - mNeutralino <= 0)  return 0;    
	double gammaForUnitRPVCoupling = 1./(128.*pi)*pow(GSL,2)*pow(abs(fSMeson),2)/pow(mtau,3)*(pow(mtau,2) - pow(mMeson,2) + pow(mNeutralino,2) ) * lambdaHalf(pow(mMeson,2), pow(mtau,2), pow(mNeutralino,2));
	return gammaForUnitRPVCoupling*pow(lambOverMSq,2);
    }


    double tensorTauWidth(double mNeutralino, double mtau, double mMeson, double fTMeson, double lambOverMSq)
    {        
        if(mtau - mMeson - mNeutralino <= 0)  return 0;    
	double gammaForUnitRPVCoupling = lambdaHalf(pow(mMeson,2), pow(mtau,2), pow(mNeutralino,2))/(2.*pi*pow(mtau,3))*pow(GTL,2)*pow(fTMeson,2)*(2*pow(pow(mtau,2) - pow(mNeutralino,2),2) - pow(mMeson,2)*(pow(mMeson,2)+pow(mtau,2)+pow(mNeutralino,2) ));
	return gammaForUnitRPVCoupling*pow(lambOverMSq,2);
    }


//Neutralino decay widths: //assuming Dirac neutrinos
    double chiralNeutralinoWidth(double mNeutralino, double mMeson, double fSMeson,double lambOverMSq)
    {
        if(mNeutralino - mMeson  <= 0)  return 0;    
	double gammaForUnitRPVCoupling = 2 *  lambdaHalf( pow(mNeutralino,2), pow(mMeson,2), 0  )/ (128.*pi* pow(mNeutralino,3)) * pow(GSnu , 2) *abs(pow(fSMeson,2)) * (pow(mNeutralino,2) - pow(mMeson,2) ); 
	return gammaForUnitRPVCoupling*pow(lambOverMSq,2);
    }


    double tensorNeutralinoWidth(double mNeutralino, double mMeson, double fTMeson, double lambOverMSq)
    {
        if(mNeutralino - mMeson  <= 0)  return 0;   
	double gammaForUnitRPVCoupling = 2 * lambdaHalf( pow(mNeutralino,2), pow(mMeson,2), 0  )/ (2*pi*pow(mNeutralino,3)) * pow(GTnu,2) * abs(pow(fTMeson,2)) * (2*pow(mNeutralino,4) - pow(mMeson,2)*( pow(mMeson,2) + pow(mNeutralino,2) )  ); 
	return gammaForUnitRPVCoupling*pow(lambOverMSq,2);
    }

}
