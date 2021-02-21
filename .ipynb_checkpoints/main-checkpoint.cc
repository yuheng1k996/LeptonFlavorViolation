#include <iostream>
//#include "widthCalculator.h"
#include "pythia_card.h"

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cout << "./main scen mX g_AB_Lamb g_CC_Lamb NMC" << std::endl;
        std::cout << "   - scen: scenario number" << std::endl;
        std::cout << "    --1: A = tau, B = e , C = e" << std::endl;
        std::cout << "    --2: A = tau, B = e , C = mu" << std::endl;
        std::cout << "    --3: A = tau, B = mu, C = e" << std::endl;
        std::cout << "    --4: A = tau, B = mu, C = mu" << std::endl;
        std::cout << "   - mX: X mass in GeV" << std::endl;
        std::cout << "   - g_AB_Lamb: g_ab/Lambda in GeV^-1" << std::endl;
        std::cout << "   - g_CC_Lamb: g_cc/Lambda in GeV^-1" << std::endl;
        std::cout << "   - NMC: number of MC simulation events" << std::endl;
        exit(1);
    }
    
    int scen = atof(argv[1]);
    double mX = atof(argv[2]);
    double g_TM_Lamb = 0;
    double g_TE_Lamb = 0;
    double g_ME_Lamb = 0;
    double g_TT_Lamb = 0;
    double g_MM_Lamb = 0;
    double g_EE_Lamb = 0;
    double g_GG_Lamb = 0;
    int nMC = atof(argv[5]);
    
    if(scen == 1){
        double g_TE_Lamb = atof(argv[3]);
        double g_EE_Lamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TE_Lamb:  " << g_TE_Lamb << std::endl;
    std::cout << "g_EE_Lamb:  " << g_EE_Lamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_Lamb(g_TM_Lamb);
    mychecker.setg_TE_Lamb(g_TE_Lamb);
    mychecker.setg_ME_Lamb(g_ME_Lamb);
    mychecker.setg_TT_Lamb(g_TT_Lamb);
    mychecker.setg_MM_Lamb(g_MM_Lamb);
    mychecker.setg_EE_Lamb(g_EE_Lamb);
    mychecker.setg_GG_Lamb(g_GG_Lamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
        
    }else if (scen == 2){
        double g_TE_Lamb = atof(argv[3]);
        double g_MM_Lamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TE_Lamb:  " << g_TE_Lamb << std::endl;
    std::cout << "g_MM_Lamb:  " << g_MM_Lamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_Lamb(g_TM_Lamb);
    mychecker.setg_TE_Lamb(g_TE_Lamb);
    mychecker.setg_ME_Lamb(g_ME_Lamb);
    mychecker.setg_TT_Lamb(g_TT_Lamb);
    mychecker.setg_MM_Lamb(g_MM_Lamb);
    mychecker.setg_EE_Lamb(g_EE_Lamb);
    mychecker.setg_GG_Lamb(g_GG_Lamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else if (scen == 3){
        double g_TM_Lamb = atof(argv[3]);
        double g_EE_Lamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TM_Lamb:  " << g_TM_Lamb << std::endl;
    std::cout << "g_EE_Lamb:  " << g_EE_Lamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_Lamb(g_TM_Lamb);
    mychecker.setg_TE_Lamb(g_TE_Lamb);
    mychecker.setg_ME_Lamb(g_ME_Lamb);
    mychecker.setg_TT_Lamb(g_TT_Lamb);
    mychecker.setg_MM_Lamb(g_MM_Lamb);
    mychecker.setg_EE_Lamb(g_EE_Lamb);
    mychecker.setg_GG_Lamb(g_GG_Lamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else if (scen == 4){
        double g_TM_Lamb = atof(argv[3]);
        double g_MM_Lamb = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TM_Lamb:  " << g_TM_Lamb << std::endl;
    std::cout << "g_MM_Lamb:  " << g_MM_Lamb << std::endl;
    std::cout << "NMC " << nMC << std::endl;
       
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_Lamb(g_TM_Lamb);
    mychecker.setg_TE_Lamb(g_TE_Lamb);
    mychecker.setg_ME_Lamb(g_ME_Lamb);
    mychecker.setg_TT_Lamb(g_TT_Lamb);
    mychecker.setg_MM_Lamb(g_MM_Lamb);
    mychecker.setg_EE_Lamb(g_EE_Lamb);
    mychecker.setg_GG_Lamb(g_GG_Lamb);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else{
        std::cout << "./main scen mX g_AB_Lamb g_CC_Lamb NMC" << std::endl;
        std::cout << "   - scen: scenario number" << std::endl;
        std::cout << "    --1: A = tau, B = e , C = e" << std::endl;
        std::cout << "    --2: A = tau, B = e , C = mu" << std::endl;
        std::cout << "    --3: A = tau, B = mu, C = e" << std::endl;
        std::cout << "    --4: A = tau, B = mu, C = mu" << std::endl;
        std::cout << "   - mX: X mass in GeV" << std::endl;
        std::cout << "   - g_AB_Lamb: g_ab/Lambda in GeV^-1" << std::endl;
        std::cout << "   - g_CC_Lamb: g_cc/Lambda in GeV^-1" << std::endl;
        std::cout << "   - NMC: number of MC simulation events" << std::endl;
        exit(1);
    }

    return 0;
}
