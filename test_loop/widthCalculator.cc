#include "widthCalculator.h"
#include <complex>

using namespace std;
typedef complex<double> dcomp;

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
    
    // loop function B1
    dcomp loopB1(double mlA, double mPs) {
        double t;
        dcomp i;
        dcomp f;
        
        t = pow(2.*mlA/mPs,2);
        i = -1;
        i = sqrt(i);
        
        if (t>=1.){
            f = asin(1./sqrt(t));
        }else if (1.>t>=0.){
            f = pi/2. + (i/2.)*(log((1.+sqrt(1.-t))/(1.-sqrt(1.-t))));
        }else{
            return 0.;
        }
        return 1.-t*f*f;
    }


//Tau decay widths:  
    double violationLepWidth(double mPs, double mlA, double mlB, double g_AB_L)
    {        
        if(mlA - mlB - mPs <= 0)  return 0;    
	double gammaForUnitLFVCoupling = lambdaHalf(pow(mlA,2), pow(mlB,2), pow(mPs,2))/(16.*pi*pow(mlA,3))*(2*pow(pow(mlA,2) - pow(mlB,2),2) - pow(mPs,2)*(pow((mlA-mlB),2)+pow(mlA+mlB,2)));
	return gammaForUnitLFVCoupling*pow(g_AB_L,2);
    }
 
//Pseudoscalar decay widths:
    double pseudoscalarLepWidth(double mPs, double mlA, double g_AA_L)
    {
        if(mPs- mlA - mlA <= 0)  return 0;    
	double gammaForUnitLFVCoupling = 2 * lambdaHalf(pow(mPs,2), pow(mlA,2), pow(mlA,2))/ (2.*pi)* pow(mlA,1); 
	return gammaForUnitLFVCoupling*pow(g_AA_L,2);
    }

    double pseudoscalarGmWidth(double mPs, double mlA, double mlB, double mlC, double g_AA_L, double g_BB_L, double g_CC_L)
    {
    double g_2;
    dcomp b;
        b = 2*g_AA_L*loopB1(mlA, mPs);
        b = b + 2*g_BB_L*loopB1(mlB, mPs);
        b = b + 2*g_CC_L*loopB1(mlC, mPs);
        g_2 = pow(afa/pi,2)*real(b*conj(b));
    double gammaForUnitGmCoupling = pow(mPs,3)/ (64*pi);
    return gammaForUnitGmCoupling*g_2;
    }
}