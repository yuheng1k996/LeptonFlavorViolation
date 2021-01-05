#include <iostream>
//#include "widthCalculator.h"
#include "pythia_card.h"

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cout << "./main scen mX g_AlphaBetaOverLamb g_AlphaAlphaOverLamb NMC" << std::endl;
        std::cout << "   - scen: scenario number" << std::endl;
        std::cout << "    --1: tau2xe, x2ee" << std::endl;
        std::cout << "    --2: tau2xe, x2mumu" << std::endl;
        std::cout << "    --3: tau2xmu, x2ee" << std::endl;
        std::cout << "    --4: tau2xmu, x2mumu" << std::endl;
        std::cout << "   - mX: X mass in GeV" << std::endl;
        std::cout << "   - g_AlphaBetaOverLamb: g_alphabeta/Lambda in GeV^-1" << std::endl;
        std::cout << "   - g_AlphaAlphaOverLamb: g_alphaalpha/Lambda in GeV^-1" << std::endl;
        std::cout << "   - NMC: number of MC simulation events" << std::endl;
        exit(1);
    }
    
    int scen = atof(argv[1]);
    double mX = atof(argv[2]);
    double g_TauMuOverLamb = 0;
    double g_TauEleOverLamb = 0;
    double g_MuEleOverLamb = 0;
    double g_TauTauOverLamb = 0;
    double g_MuMuOverLamb = 0;
    double g_EleEleOverLamb = 0;
    double g_GmGmOverLamb = 0;
    int nMC = atof(argv[5]);
    
    if(scen == 1){
        double g_TauEleOverLamb = atof(argv[3]);
        double g_EleEleOverLamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TauEleOverLamb:  " << g_TauEleOverLamb << std::endl;
    std::cout << "g_EleEleOverLamb:  " << g_EleEleOverLamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TauMuOverLamb(g_TauMuOverLamb);
    mychecker.setg_TauEleOverLamb(g_TauEleOverLamb);
    mychecker.setg_MuEleOverLamb(g_MuEleOverLamb);
    mychecker.setg_TauTauOverLamb(g_TauTauOverLamb);
    mychecker.setg_MuMuOverLamb(g_MuMuOverLamb);
    mychecker.setg_EleEleOverLamb(g_EleEleOverLamb);
    mychecker.setg_GmGmOverLamb(g_GmGmOverLamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
        
    }else if (scen == 2){
        double g_TauEleOverLamb = atof(argv[3]);
        double g_MuMuOverLamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TauEleOverLamb:  " << g_TauEleOverLamb << std::endl;
    std::cout << "g_MuMuOverLamb:  " << g_MuMuOverLamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TauMuOverLamb(g_TauMuOverLamb);
    mychecker.setg_TauEleOverLamb(g_TauEleOverLamb);
    mychecker.setg_MuEleOverLamb(g_MuEleOverLamb);
    mychecker.setg_TauTauOverLamb(g_TauTauOverLamb);
    mychecker.setg_MuMuOverLamb(g_MuMuOverLamb);
    mychecker.setg_EleEleOverLamb(g_EleEleOverLamb);
    mychecker.setg_GmGmOverLamb(g_GmGmOverLamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else if (scen == 3){
        double g_TauMuOverLamb = atof(argv[3]);
        double g_EleEleOverLamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TauMuOverLamb:  " << g_TauMuOverLamb << std::endl;
    std::cout << "g_EleEleOverLamb:  " << g_EleEleOverLamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TauMuOverLamb(g_TauMuOverLamb);
    mychecker.setg_TauEleOverLamb(g_TauEleOverLamb);
    mychecker.setg_MuEleOverLamb(g_MuEleOverLamb);
    mychecker.setg_TauTauOverLamb(g_TauTauOverLamb);
    mychecker.setg_MuMuOverLamb(g_MuMuOverLamb);
    mychecker.setg_EleEleOverLamb(g_EleEleOverLamb);
    mychecker.setg_GmGmOverLamb(g_GmGmOverLamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else if (scen == 4){
        double g_TauMuOverLamb = atof(argv[3]);
        double g_MuMuOverLamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TauMuOverLamb:  " << g_TauMuOverLamb << std::endl;
    std::cout << "g_MuMuOverLamb:  " << g_MuMuOverLamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
       
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TauMuOverLamb(g_TauMuOverLamb);
    mychecker.setg_TauEleOverLamb(g_TauEleOverLamb);
    mychecker.setg_MuEleOverLamb(g_MuEleOverLamb);
    mychecker.setg_TauTauOverLamb(g_TauTauOverLamb);
    mychecker.setg_MuMuOverLamb(g_MuMuOverLamb);
    mychecker.setg_EleEleOverLamb(g_EleEleOverLamb);
    mychecker.setg_GmGmOverLamb(g_GmGmOverLamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else{
        std::cout << "./main scen mX g_AlphaBetaOverLamb g_AlphaAlphaOverLamb NMC" << std::endl;
        std::cout << "   - scen: scenario number" << std::endl;
        std::cout << "    --1: g_TauEleOverLamb, g_EleEleOverLamb" << std::endl;
        std::cout << "    --2: g_TauEleOverLamb, g_MuMuOverLamb" << std::endl;
        std::cout << "    --3: g_TauMuOverLamb, g_EleEleOverLamb" << std::endl;
        std::cout << "    --4: g_TauMuOverLamb, g_MuMuOverLamb" << std::endl;
        std::cout << "   - mX: X mass in GeV" << std::endl;
        std::cout << "   - g_AlphaBetaOverLamb: g_alphabeta/Lambda in GeV^-1" << std::endl;
        std::cout << "   - g_AlphaAlphaOverLamb: g_alphaalpha/Lambda in GeV^-1" << std::endl;
        std::cout << "   - NMC: number of MC simulation events" << std::endl;
        exit(1);
    }

    return 0;
}
