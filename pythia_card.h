#ifndef _PYTHIACARD
#define _PYTHIACARD


#include <exception>
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include "Pythia.h"
#include "widthCalculator.h"


class pythia_card {
public:
    pythia_card();
    ~pythia_card() {};    

    void setVerbose() {verbose = true;};
    
    //lepton flavor violation    
    void setMX(double mX_in) {mX = mX_in;};
    void setg_TM_Lamb(double g_TM_Lamb_in) {g_TM_Lamb = g_TM_Lamb_in;};
    void setg_TE_Lamb(double g_TE_Lamb_in) {g_TE_Lamb = g_TE_Lamb_in;};
    void setg_ME_Lamb(double g_ME_Lamb_in) {g_ME_Lamb = g_ME_Lamb_in;};
    void setg_TT_Lamb(double g_TT_Lamb_in) {g_TT_Lamb = g_TT_Lamb_in;};
    void setg_MM_Lamb(double g_MM_Lamb_in) {g_MM_Lamb = g_MM_Lamb_in;};
    void setg_EE_Lamb(double g_EE_Lamb_in) {g_EE_Lamb = g_EE_Lamb_in;};
    void setg_GG_Lamb(double g_GG_Lamb_in) {g_GG_Lamb = g_GG_Lamb_in;};
    
    bool doCalculations(); //< evaluates widths 
    bool initPythia(); //< Initialises Pythia, if needed. 
    bool runPythia(int nEventsMC);
    
        
private: 
     

    Pythia8::Pythia* pythia; //< Pythia8 object for simulation
    bool verbose; //< declares amount of Information Pythia writes
    
    double number_total_tau_tau;//total number of ee->tautau events 
     
    bool decayingInsideFidVol(Pythia8::Particle XXX); 
    double detectorEffi(Pythia8::Particle XXX);
    double decayProbabilityBelle2Part1(Pythia8::Particle XXX);
    double decayProbabilityBelle2Part2(Pythia8::Particle XXX);
   
    //lepton flavor violation
    double mX; // given in GeV
    double g_TM_Lamb; // g_taumu/Lambda in GeV^-1
    double g_TE_Lamb; // g_taue/Lambda in GeV^-1
    double g_ME_Lamb; // g_mue/Lambda in GeV^-1
    double g_TT_Lamb; // g_tautau/Lambda in GeV^-1
    double g_MM_Lamb; // g_mumu/Lambda in GeV^-1
    double g_EE_Lamb; // g_ee/Lambda in GeV^-1
    double g_GG_Lamb; // g_gmgm/Lambda in GeV^-1
    
    double Gammax2tautau(double mX, double g_TT_Lamb);
    double Gammax2mumu(double mX, double g_MM_Lamb);
    double Gammax2ee(double mX, double g_EE_Lamb);
    double Gammax2gmgm(double mX, double g_GG_Lamb);
    
    double totalGammaX(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb);
    double ctau(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb);
    
    double BRx2tautau(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb);
    double BRx2mumu(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb);
    double BRx2ee(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb);
    double BRx2gmgm(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb);
    
    double BRx2Visibles(double mX, double g_TT_lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb);
    
    double Gammatau2xmu(double mX, double g_TM_Lamb);
    double Gammatau2xe(double mX, double g_TE_Lamb);
    double NewTotalGammatau(double mX, double g_TM_Lamb, double g_TE_Lamb);
    double Gammamu2xe(double mX, double g_ME_Lamb);
    double NewTotalGammamu(double mX, double g_ME_Lamb);
    
    double BRtau2xmu(double mX, double g_TM_Lamb, double g_TE_Lamb);
    double BRtau2xe(double mX, double g_TM_Lamb, double g_TE_Lamb);
    double BRmu2xe(double mX, double g_ME_Lamb);
    
    double reallyProducedX;
    double reallyObservedX;
    
};
#endif
