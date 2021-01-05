
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
    void setg_TauMuOverLamb(double g_TauMuOverLamb_in) {g_TauMuOverLamb = g_TauMuOverLamb_in;};
    void setg_TauEleOverLamb(double g_TauEleOverLamb_in) {g_TauEleOverLamb = g_TauEleOverLamb_in;};
    void setg_MuEleOverLamb(double g_MuEleOverLamb_in) {g_MuEleOverLamb = g_MuEleOverLamb_in;};
    void setg_TauTauOverLamb(double g_TauTauOverLamb_in) {g_TauTauOverLamb = g_TauTauOverLamb_in;};
    void setg_MuMuOverLamb(double g_MuMuOverLamb_in) {g_MuMuOverLamb = g_MuMuOverLamb_in;};
    void setg_EleEleOverLamb(double g_EleEleOverLamb_in) {g_EleEleOverLamb = g_EleEleOverLamb_in;};
    void setg_GmGmOverLamb(double g_GmGmOverLamb_in) {g_GmGmOverLamb = g_GmGmOverLamb_in;};
    
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
    double g_TauMuOverLamb; // g_taumu/Lambda in GeV^-1
    double g_TauEleOverLamb; // g_taue/Lambda in GeV^-1
    double g_MuEleOverLamb; // g_mue/Lambda in GeV^-1
    double g_TauTauOverLamb; // g_tautau/Lambda in GeV^-1
    double g_MuMuOverLamb; // g_mumu/Lambda in GeV^-1
    double g_EleEleOverLamb; // g_ee/Lambda in GeV^-1
    double g_GmGmOverLamb; // g_gmgm/Lambda in GeV^-1
    
    double Gammax2tautau(double mX, double g_TauTauOverLamb);
    double Gammax2mumu(double mX, double g_MuMuOverLamb);
    double Gammax2ee(double mX, double g_EleEleOverLamb);
    double Gammax2gmgm(double mX, double g_GmGmOverLamb);
    
    double totalGammaX(double mX, double g_TauTauOverLamb, double g_MuMuOverLamb, double g_EleEleOverLamb, double g_GmGmOverlamb);
    double ctau(double mX, double g_TauTauOverLamb, double g_MuMuOverLamb, double g_EleEleOverLamb, double g_GmGmOverLamb);
    
    double BRx2tautau(double mX, double g_TauTauOverLamb, double g_MuMuOverLamb, double g_EleEleOverLamb, double g_GmGmOverlamb);
    double BRx2mumu(double mX, double g_TauTauOverLamb, double g_MuMuOverLamb, double g_EleEleOverLamb, double g_GmGmOverLamb);
    double BRx2ee(double mX, double g_TauTauOverLamb, double g_MuMuOverLamb, double g_EleEleOverLamb, double g_GmGmOverLamb);
    double BRx2gmgm(double mX, double g_TauTauOverLamb, double g_MuMuOverLamb, double g_EleEleOverLamb, double g_GmGmOverLamb);
    
    double BRx2Visibles(double mX, double g_TauTauOverlamb, double g_MuMuOverLamb, double g_EleEleOverLamb, double g_GmGmOverlamb);
    
    double Gammatau2xmu(double mX, double g_TauMuOverLamb);
    double Gammatau2xe(double mX, double g_TauEleOverLamb);
    double NewTotalGammatau(double mX, double g_TauMuOverLamb, double g_TauEleOverLamb);
    double Gammamu2xe(double mX, double g_MuEleOverLamb);
    double NewTotalGammamu(double mX, double g_MuEleOverLamb);
    
    double BRtau2xmu(double mX, double g_TauMuOverLamb, double g_TauEleOverLamb);
    double BRtau2xe(double mX, double g_TauMuOverLamb, double g_TauEleOverLamb);
    double BRmu2xe(double mX, double g_MuEleOverLamb);
    
    double reallyProducedX;
    double reallyObservedX;
    
};
#endif
