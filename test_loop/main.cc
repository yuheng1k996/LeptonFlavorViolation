#include <iostream>
//#include "widthCalculator.h"
#include "pythia_card.h"

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cout << "./main scen mX g_AB_L g_CC_L NMC" << std::endl;
        std::cout << "   - scen: scenario number" << std::endl;
        std::cout << "    --0: A = a  , B = b , C = c" << std::endl;
        std::cout << "    --1: A = tau, B = e , C = e" << std::endl;
        std::cout << "    --2: A = tau, B = e , C = mu" << std::endl;
        std::cout << "    --3: A = tau, B = mu, C = e" << std::endl;
        std::cout << "    --4: A = tau, B = mu, C = mu" << std::endl;
        std::cout << "   - mX: X mass in GeV" << std::endl;
        std::cout << "   - g_AB_L: g_ab/Lambda in GeV^-1" << std::endl;
        std::cout << "   - g_CC_L: g_cc/Lambda in GeV^-1" << std::endl;
        std::cout << "   - NMC: number of MC simulation events" << std::endl;
        exit(1);
    }
    
    int scen = atof(argv[1]);
    double mX = atof(argv[2]);
    double g_TM_L = 0;
    double g_TE_L = 0;
    double g_ME_L = 0;
    double g_TT_L = 0;
    double g_MM_L = 0;
    double g_EE_L = 0;
    int nMC = atof(argv[5]);

   
    if(scen == 0){
        //off-diagonal coupling
        double g_TM_L = atof(argv[3]);
        double g_TE_L = g_TM_L;
        double g_ME_L = g_TM_L;
        //diagonal coulping
        double g_TT_L = atof(argv[4]);
        double g_MM_L = g_TT_L;
        double g_EE_L = g_TT_L;    
        
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_AB_L:  " << g_TE_L << std::endl;
    std::cout << "g_CC_L:  " << g_EE_L << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_L(g_TM_L);
    mychecker.setg_TE_L(g_TE_L);
    mychecker.setg_ME_L(g_ME_L);
    mychecker.setg_TT_L(g_TT_L);
    mychecker.setg_MM_L(g_MM_L);
    mychecker.setg_EE_L(g_EE_L);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
        
    }else if(scen == 1){
        double g_TE_L = atof(argv[3]);
        double g_EE_L = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TE_L:  " << g_TE_L << std::endl;
    std::cout << "g_EE_L:  " << g_EE_L << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_L(g_TM_L);
    mychecker.setg_TE_L(g_TE_L);
    mychecker.setg_ME_L(g_ME_L);
    mychecker.setg_TT_L(g_TT_L);
    mychecker.setg_MM_L(g_MM_L);
    mychecker.setg_EE_L(g_EE_L);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
        
    }else if (scen == 2){
        double g_TE_L = atof(argv[3]);
        double g_MM_L = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TE_L:  " << g_TE_L << std::endl;
    std::cout << "g_MM_L:  " << g_MM_L << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_L(g_TM_L);
    mychecker.setg_TE_L(g_TE_L);
    mychecker.setg_ME_L(g_ME_L);
    mychecker.setg_TT_L(g_TT_L);
    mychecker.setg_MM_L(g_MM_L);
    mychecker.setg_EE_L(g_EE_L);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else if (scen == 3){
        double g_TM_L = atof(argv[3]);
        double g_EE_L = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TM_L:  " << g_TM_L << std::endl;
    std::cout << "g_EE_L:  " << g_EE_L << std::endl;
    std::cout << "NMC " << nMC << std::endl;
        
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_L(g_TM_L);
    mychecker.setg_TE_L(g_TE_L);
    mychecker.setg_ME_L(g_ME_L);
    mychecker.setg_TT_L(g_TT_L);
    mychecker.setg_MM_L(g_MM_L);
    mychecker.setg_EE_L(g_EE_L);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else if (scen == 4){
        double g_TM_L = atof(argv[3]);
        double g_MM_L = atof(argv[4]);
    std::cout << "mX:  " << mX << std::endl;
    std::cout << "g_TM_L:  " << g_TM_L << std::endl;
    std::cout << "g_MM_L:  " << g_MM_L << std::endl;
    std::cout << "NMC " << nMC << std::endl;
       
    pythia_card mychecker;
    mychecker.setVerbose();
    mychecker.setMX(mX);
    mychecker.setg_TM_L(g_TM_L);
    mychecker.setg_TE_L(g_TE_L);
    mychecker.setg_ME_L(g_ME_L);
    mychecker.setg_TT_L(g_TT_L);
    mychecker.setg_MM_L(g_MM_L);
    mychecker.setg_EE_L(g_EE_L);
    
    if (!mychecker.doCalculations())
        return 1;
    if (!mychecker.initPythia())
        return 1;
    if (!mychecker.runPythia(nMC))
        return 1;
     
    }else{
        std::cout << "./main scen mX g_AB_L g_CC_L NMC" << std::endl;
        std::cout << "   - scen: scenario number" << std::endl;
        std::cout << "    --0: A = a  , B = b , C = c" << std::endl;
        std::cout << "    --1: A = tau, B = e , C = e" << std::endl;
        std::cout << "    --2: A = tau, B = e , C = mu" << std::endl;
        std::cout << "    --3: A = tau, B = mu, C = e" << std::endl;
        std::cout << "    --4: A = tau, B = mu, C = mu" << std::endl;
        std::cout << "   - mX: X mass in GeV" << std::endl;
        std::cout << "   - g_AB_L: g_ab/Lambda in GeV^-1" << std::endl;
        std::cout << "   - g_CC_L: g_cc/Lambda in GeV^-1" << std::endl;
        std::cout << "   - NMC: number of MC simulation events" << std::endl;
        exit(1);
    }

    return 0;
}
