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
    double violationLepWidth(double mPseudoscalar, double mlAlpha, double mlBeta, double g_AlphaBetaOverLamb)
    {        
        if(mlAlpha - mlBeta - mPseudoscalar <= 0)  return 0;    
	double gammaForUnitLFVCoupling = lambdaHalf(pow(mlAlpha,2), pow(mlBeta,2), pow(mPseudoscalar,2))/(16.*pi*pow(mlAlpha,3))*(2*pow(pow(mlAlpha,2) - pow(mlBeta,2),2) - pow(mPseudoscalar,2)*(pow((mlAlpha-mlBeta),2)+pow(mlAlpha+mlBeta,2)));
	return gammaForUnitLFVCoupling*pow(g_AlphaBetaOverLamb,2);
    }
 
//Pseudoscalar decay widths:
    double pseudoscalarLepWidth(double mPseudoscalar, double mlAlpha, double g_AlphaAlphaOverLamb)
    {
        if(mPseudoscalar - mlAlpha - mlAlpha <= 0)  return 0;    
	double gammaForUnitLFVCoupling = 2 * lambdaHalf(pow(mPseudoscalar,2), pow(mlAlpha,2), pow(mlAlpha,2))/ (2.*pi)* pow(mlAlpha,1); 
	return gammaForUnitLFVCoupling*pow(g_AlphaAlphaOverLamb,2);
    }

    double pseudoscalarGmWidth(double mPseudoscalar, double g_GmGmOverLamb)
    {
    double gammaForUnitGmCoupling = pow(mPseudoscalar,3)/ (64.*pi);
    return gammaForUnitGmCoupling*pow(g_GmGmOverLamb,3);
    }
}