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
    void setg_TM_L(double g_TM_L_in) {g_TM_L = g_TM_L_in;};
    void setg_TE_L(double g_TE_L_in) {g_TE_L = g_TE_L_in;};
    void setg_ME_L(double g_ME_L_in) {g_ME_L = g_ME_L_in;};
    void setg_TT_L(double g_TT_L_in) {g_TT_L = g_TT_L_in;};
    void setg_MM_L(double g_MM_L_in) {g_MM_L = g_MM_L_in;};
    void setg_EE_L(double g_EE_L_in) {g_EE_L = g_EE_L_in;};
    
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
    double g_TM_L; // g_taumu/Lambda in GeV^-1
    double g_TE_L; // g_taue/Lambda in GeV^-1
    double g_ME_L; // g_mue/Lambda in GeV^-1
    double g_TT_L; // g_tautau/Lambda in GeV^-1
    double g_MM_L; // g_mumu/Lambda in GeV^-1
    double g_EE_L; // g_ee/Lambda in GeV^-1
    
    double Gammax2tautau(double mX, double g_TT_L);
    double Gammax2mumu(double mX, double g_MM_L);
    double Gammax2ee(double mX, double g_EE_L);
    double Gammax2gmgm(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    
    double totalGammaX(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    double ctau(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    
    double BRx2tautau(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    double BRx2mumu(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    double BRx2ee(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    double BRx2gmgm(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    
    double BRx2Visibles(double mX, double g_TT_L, double g_MM_L, double g_EE_L);
    
    double Gammatau2xmu(double mX, double g_TM_L);
    double Gammatau2xe(double mX, double g_TE_L);
    double NewTotalGammatau(double mX, double g_TM_L, double g_TE_L);
    double Gammamu2xe(double mX, double g_ME_L);
    double NewTotalGammamu(double mX, double g_ME_L);
    
    double BRtau2xmu(double mX, double g_TM_L, double g_TE_L);
    double BRtau2xe(double mX, double g_TM_L, double g_TE_L);
    double BRmu2xe(double mX, double g_ME_L);
    
    double reallyProducedX;
    double reallyObservedX;
    
};
#endif
