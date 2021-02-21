#include "pythia_card.h"
#include <iostream>
#include <fstream>
using namespace std;

// consider two DVs:  1) x -> l lbar
//		              2) x -> gm gm

// Translates a number into a string
std::string floatToString(double x) {
    std::ostringstream ss;
    ss << x;
    std::string s(ss.str());
    return s;
}

// Translates a number into a string
std::string intToString(int x) {
    std::ostringstream ss;
    ss << x;
    std::string s(ss.str());
    return s;
}



pythia_card::pythia_card() {
    pythia = new Pythia8::Pythia("../xmldoc/", false);
    verbose = false;
    mX = 0;  //mass of X in GeV
    g_TM_Lamb = 0;
    g_TE_Lamb = 0;
    g_ME_Lamb = 0;
    g_TT_Lamb = 0;
    g_MM_Lamb = 0;
    g_EE_Lamb = 0;
    g_GG_Lamb = 0;
    reallyProducedX = 0;//total number of X produced at the real experiment
};    


bool pythia_card::doCalculations() {
    number_total_tau_tau = 4.6e10;//total number of ee->tautau events at Belle II with 50/ab
    double BrTau2Leptons = 0.3521;//PDG
  
    reallyProducedX = number_total_tau_tau * (2 * BrTau2Leptons * (BRtau2xmu(mX, g_TM_Lamb, g_TE_Lamb) + BRtau2xe(mX, g_TM_Lamb, g_TE_Lamb)));

    return true;
}



bool pythia_card::initPythia() {
    try {
        //set beam parameters
        pythia->readString("Beams:idA = 11");
        pythia->readString("Beams:idB = -11");
        pythia->readString("Beams:frameType = 2");
        pythia->readString("Beams:eA = 7.");
        pythia->readString("Beams:eB = 4.");

        pythia->readString("WeakSingleBoson:ffbar2ffbar(s:gm) = on"); //off-shell photon -> 15% efficiency (N_{tautau}/NMC~15%) but fast.

	//The following gives 100% efficiency to tautau, but very slow
        //pythia->readString("WeakSingleBoson:ffbar2ffbar(s:gmZ) = on"); //off-shell photon or Z0
        //pythia->readString("23:oneChannel =  1 1 103 15 -15 ");//Z->tau+ tau-
        

        //set random seed on 
        pythia->readString("Random:setSeed = on");
        pythia->readString("Random:seed = 0");//pick new random number seed for each run, based on clock     
        
        //to speed up the simulation.
        pythia->readString("PartonLevel:MPI = 0");//should switch off. small influence on pT of 1000022. speeds up alot.
        //pythia->readString("PartonLevel:ISR = 0");//should not switch off ISR. large influence
        //pythia->readString("PartonLevel:FSR = 0");//should not switch off ISR. large influence
        pythia->readString("HadronLevel:Hadronize = 0"); //should switch off, as it has no influence but speeds up   		  

        //add a new pseudoscalar particle
        pythia->readString("36:new = x");
        pythia->readString("36:spinType = 1");//scalar boson
        pythia->readString("36:chargeType = 0");
        pythia->readString("36:colType = 0");
        pythia->readString("36:m0 = "+floatToString(mX)); // pythia masses are in GeV

	    //set x decay
        pythia->readString("36:mayDecay = 1 ");
        //converting ctau in meter to mm by the factor 1.E3
        pythia->readString("36:tau0 = "+floatToString(ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)*1.E3));
        pythia->readString("36:tau0 = "+floatToString(ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)*1.E3));
        pythia->readString("36:tau0 = "+floatToString(ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)*1.E3));
        //x -> l lbar
        pythia->readString("36:oneChannel =  1 "+floatToString(BRx2tautau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))+ " 103 15 -15 ");
        pythia->readString("36:addChannel =  1 "+floatToString(BRx2mumu(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))+ " 103 13 -13 ");
        pythia->readString("36:addChannel =  1 "+floatToString(BRx2ee(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))+ " 103 11 -11 ");
        pythia->readString("36:addChannel =  1 "+floatToString(BRx2gmgm(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))+ " 103 22 22 ");//fermion loop diagram

	    //set tau decay channels
        pythia->readString("15:oneChannel =  1 "+floatToString(BRtau2xmu(mX,g_TM_Lamb,g_TE_Lamb))+" 103 36 13 ");//tau ->x mu
        pythia->readString("15:addChannel =  1 "+floatToString(BRtau2xe(mX,g_TM_Lamb,g_TE_Lamb))+" 103 36 11 ");//tau ->x e
	
        
        // Some general things regarding output
        if(verbose) {
            pythia->readString("Init:showProcesses = on");
            pythia->readString("Init:showChangedSettings = on");
            //pythia->readString("Main:showChangedSettings = on");
            pythia->readString("Stat:showProcessLevel = on");
            pythia->readString("Stat:showErrors = on");
            pythia->readString("Init:showProcesses = on");
            pythia->readString("Print:quiet = off");
            pythia->readString("Init:showChangedParticleData = on");
            pythia->readString("Next:numberCount = 10000");
        }
        else {
            pythia->readString("Init:showProcesses = off");
            pythia->readString("Init:showChangedSettings = off");
            //pythia->readString("Main:showChangedSettings = off");
            pythia->readString("Stat:showProcessLevel = off");
            pythia->readString("Stat:showErrors = off");
            pythia->readString("Init:showProcesses = off");
            pythia->readString("Print:quiet = on");
            pythia->readString("Init:showChangedSettings = off");
            pythia->readString("Init:showChangedParticleData = off");
        }
        pythia->init(); 
    }
    catch(std::exception& e) {
        std::cerr << "!!! Error occured while trying to initialise Pythia: " << e.what() << std::endl;
        return false;
    }
    return true;
}

bool pythia_card::runPythia(int nEventsMC) {

    double producedtau = 0;
    double producedX = 0;
    double observedX = 0;
    double visibleX = 0;
    
    double XdecayinsideFidVol = 0;
    
    int xForwardSphere = 0;
    int xBackwardSphere = 0;

//Creating files to scan
    try{
        ofstream myfile;
        myfile.open ("./scanning/data_g_TM_" +floatToString(g_TM_Lamb)+ "_g_TE_" +floatToString(g_TE_Lamb)+ "_g_MM_" +floatToString(g_MM_Lamb)+ "_g_EE_" +floatToString(g_EE_Lamb)+ ".txt");
        
        for (int iEvent = 0; iEvent < nEventsMC; ++iEvent) {
                if (!pythia->next()) continue;
                // Check the list of final state particles
                for (int i = 0; i < pythia->event.size(); ++i)
		{
			//count tau leptons, (whose daughters are not tau leptons)
                	if (abs(pythia->event[i].id()) == 15 && abs(pythia->event[pythia->event[i].daughter1()].id()) != 15)
			{
                        	producedtau += 1;
			}
                	if (abs(pythia->event[i].id()) == 36)
			{	//count x
				producedX += 1;
				double p1 = decayProbabilityBelle2Part1(pythia->event[i]);
				double p2 = decayProbabilityBelle2Part2(pythia->event[i]);
                        	observedX += p1+p2;
				visibleX += (p1+p2)*BRx2Visibles(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)*detectorEffi(pythia->event[i]);
				if (decayingInsideFidVol(pythia->event[i])) XdecayinsideFidVol+=1;
				if (pythia->event[i].p().eta()>0)xForwardSphere+=1;
				if (pythia->event[i].p().eta()<0)xBackwardSphere+=1;
			}
        }
         }
        if(verbose)
            pythia->stat();
        //Close the scanning files
        myfile << producedX << endl;
        myfile.close();
    }
    catch(std::exception& e) {
        std::cerr << "!!! Error occured while trying to run Pythia: " << e.what() << std::endl;
        return false;
    }


    double reallyobservedX = observedX/double(producedX) * reallyProducedX;
    double observedFiducialEfficiencyX = observedX/double(producedX);

    double reallyvisibleX = visibleX/double(producedX) * reallyProducedX;
    double visibleFiducialEfficiencyX = visibleX/double(producedX);

    // Results
    std::cout << "Gammatau2xmu: " << Gammatau2xmu(mX,g_TM_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "Gammatau2xe: " << Gammatau2xe(mX,g_TE_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "BRtau2xmu: " << BRtau2xmu(mX,g_TM_Lamb,g_TE_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "BRtau2xe: " << BRtau2xe(mX,g_TM_Lamb,g_TE_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "produced tau excluding those daughters are also tau lepton: " << producedtau << '\n';
    std::cout << "produced X: " << producedX << '\n';   
    std::cout << "produced X/NMC: " << producedX/double(nEventsMC) << '\n';
    std::cout << std::endl;
    std::cout << "Proportion of X flying forward (eta>0): " << xForwardSphere / producedX << '\n';
    std::cout << "Proportion of X flying backward (eta<0): " << xBackwardSphere / producedX << '\n';
    std::cout << '\n';
    std::cout << "Total Gamma [GeV]: " << totalGammaX(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb) << '\n';
    std::cout << "ctau [m]: " << ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "Gammax2tautau: " << Gammax2tautau(mX,g_TT_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "Gammax2mumu: " << Gammax2mumu(mX,g_MM_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "Gammax2ee: " << Gammax2ee(mX,g_EE_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "Gammax2gmgm: " << Gammax2gmgm(mX,g_GG_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "BRx2tautau: " << BRx2tautau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "BRx2mumu: " << BRx2mumu(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "BRx2ee: " << BRx2ee(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "BRx2gmgm: " << BRx2gmgm(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb) << '\n';
    std::cout << '\n';
    std::cout << "X: insideFidVol: " << XdecayinsideFidVol << '\n';
    std::cout << "observedX: " << observedX << '\n';
    std::cout << "visibleX: " << visibleX << '\n';
    std::cout << '\n';
    std::cout << "X Observed Fiducial efficiency: " << observedFiducialEfficiencyX << '\n';
    std::cout << "X Visible Fiducial efficiency: " << visibleFiducialEfficiencyX << '\n';
    std::cout << '\n';
    std::cout << "reallyProduced tau: " << 2*number_total_tau_tau << '\n';
    std::cout << "total number of tau-tau events: " << number_total_tau_tau << '\n';
    //std::cout << "Br(tau -> leptons): " << "35.21%" << '\n';
    std::cout << '\n';
    std::cout << "reallyProducedX: " << reallyProducedX << '\n';
    std::cout << "reallyobservedX: " << reallyobservedX << '\n';
    std::cout << "reallyvisibleX: " << reallyvisibleX << '\n';
    std::cout << '\n';
    
    return true;
}



double pythia_card::detectorEffi(Pythia8::Particle XXX)//detector efficiency
{
	return 0.1;
}


//lepton flavor violation
double pythia_card::Gammatau2xmu(double mX, double g_TM_Lamb)
{
	return widthCalculator::violationLepWidth(mX, widthCalculator::mtau, widthCalculator::mmu, g_TM_Lamb);
}


double pythia_card::Gammatau2xe(double mX, double g_TE_Lamb)
{
	return widthCalculator::violationLepWidth(mX, widthCalculator::mtau, widthCalculator::me, g_TE_Lamb);
}


double pythia_card::NewTotalGammatau(double mX, double g_TM_Lamb, double g_TE_Lamb)
{
	return widthCalculator::tauSMGamma + Gammatau2xmu(mX, g_TM_Lamb) + Gammatau2xe(mX, g_TE_Lamb);
}


double pythia_card::BRtau2xmu(double mX, double g_TM_Lamb, double g_TE_Lamb)
{
	return Gammatau2xmu(mX, g_TM_Lamb) / NewTotalGammatau(mX, g_TM_Lamb, g_TE_Lamb);
}

double pythia_card::BRtau2xe(double mX, double g_TM_Lamb, double g_TE_Lamb)
{
	return Gammatau2xe(mX, g_TE_Lamb / NewTotalGammatau(mX, g_TM_Lamb, g_TE_Lamb));
}


double pythia_card::Gammamu2xe(double mX, double g_ME_Lamb)
{
	return widthCalculator::violationLepWidth(mX, widthCalculator::mmu, widthCalculator::me, g_ME_Lamb);
}


double pythia_card::NewTotalGammamu(double mX, double g_ME_Lamb)
{
    return widthCalculator::muSMGamma + Gammamu2xe(mX, g_ME_Lamb);
}


double pythia_card::BRmu2xe(double mX, double g_ME_Lamb)
{
    return Gammamu2xe(mX, g_ME_Lamb) / NewTotalGammamu(mX, g_ME_Lamb);
}


double pythia_card::Gammax2tautau(double mX, double g_TT_Lamb)
{
	return widthCalculator::pseudoscalarLepWidth(mX,widthCalculator::mtau,g_TT_Lamb);
}

double pythia_card::Gammax2mumu(double mX, double g_MM_Lamb)
{
	return widthCalculator::pseudoscalarLepWidth(mX,widthCalculator::mmu,g_MM_Lamb);
}

double pythia_card::Gammax2ee(double mX, double g_EE_Lamb)
{
	return widthCalculator::pseudoscalarLepWidth(mX,widthCalculator::me,g_EE_Lamb);
}

double pythia_card::Gammax2gmgm(double mX, double g_GG_Lamb)
{
	return widthCalculator::pseudoscalarGmWidth(mX,g_GG_Lamb);
}


double pythia_card::totalGammaX(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb)
{
	return Gammax2tautau(mX,g_TT_Lamb) + Gammax2mumu(mX,g_MM_Lamb) + Gammax2ee(mX,g_EE_Lamb) + Gammax2gmgm(mX,g_GG_Lamb);
}


double pythia_card::ctau(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_lamb)//ctau in meter
{
	double chbar = 1.973269631e-16;//GeV m
	return chbar/ totalGammaX(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb);
}


double pythia_card::BRx2tautau(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb)
{
	return Gammax2tautau(mX,g_TT_Lamb)/totalGammaX(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb);
}

double pythia_card::BRx2mumu(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb)
{
	return Gammax2mumu(mX,g_MM_Lamb)/totalGammaX(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb);
}

double pythia_card::BRx2ee(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb)
{
	return Gammax2ee(mX,g_EE_Lamb)/totalGammaX(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb);
}

double pythia_card::BRx2gmgm(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb)
{
	return Gammax2gmgm(mX,g_GG_Lamb)/totalGammaX(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb);
}


double pythia_card::BRx2Visibles(double mX, double g_TT_Lamb, double g_MM_Lamb, double g_EE_Lamb, double g_GG_Lamb)
{
	constexpr double BRx2Visibles=1.;//need to change if charged (k*0 -> K^+ pi^-)
	return   (BRx2tautau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)+BRx2mumu(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)+BRx2ee(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)+BRx2gmgm(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))*BRx2Visibles;
}


bool pythia_card::decayingInsideFidVol(Pythia8::Particle XXX)
{
    constexpr double r_min = 100.;//in mm
    constexpr double r_max = 800.; 
    constexpr double z_min = -400.;
    constexpr double z_max = 1200.;

    double rDec = sqrt(pow(XXX.xDec(),2)+pow(XXX.yDec(),2));
    double zDec = XXX.zDec();
    
    if (rDec > r_min and rDec < r_max and zDec > z_min and zDec < z_max)
    return true;
    else return false;
}

double pythia_card::decayProbabilityBelle2Part1(Pythia8::Particle XXX) {//Part1 is a symmetric ATLAS-like detector
    constexpr double R_I = 0.1;//in meter
    constexpr double R_O = 0.8; 
    constexpr double L_D = 0.4;

    // Identify the kinematic properties of the pseudoscalar
    //Pythia always calculates in GeV
    double gamma = XXX.e()/(mX);
    double beta_z = fabs(XXX.pz()/XXX.e());
    double beta = sqrt(1. - pow(mX/XXX.e(), 2));
    double theta = XXX.p().theta();            
    double phi = XXX.p().phi();     
 

    if (tan(theta) == 0.0) {
        std::cout << "The impossible happened!" << std::endl;
        return 0;
    }      
    double L1 = fabs(R_I/tan(theta));
    if (L1 >= L_D)
        return 0;
    double L2 = std::min(L_D, fabs(R_O/tan(theta))) - L1;
    // The probability that the pseudoscalar would decay in the detector is then as follows (Decay law:            
    return exp(-L1/(beta_z*gamma*ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))) * (1. - exp(-L2/(beta_z*gamma*ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))));
}


double pythia_card::decayProbabilityBelle2Part2(Pythia8::Particle XXX) {//Part2 is AL3X-like
    constexpr double L_D = 0.4;//horizontal distance from the IP to the detector
    constexpr double L_L = 0.1;//vertical distance from the IP to the detector  or  equivalently the inner radius
    constexpr double L_d = 0.8;// length of the detector: 120 - 40 = 80 cm
    constexpr double L_H = 0.7;// transverse length = 80 - 10 = 70 cm

    // Identify the kinematic properties of the pseudoscalar
    //Pythia always calculates in GeV
    double gamma = XXX.e()/(mX);
    double beta_z = fabs(XXX.pz()/XXX.e());
    double beta = sqrt(1. - pow(mX/XXX.e(), 2));
    double theta = XXX.p().theta();            
    double phi = XXX.p().phi();     
 

    if (tan(theta) == 0.0) {
        std::cout << "The impossible happened!" << std::endl;
        return 0;
    }   
    double L1 = std::min(std::max(L_D,L_L/tan(theta)), L_D + L_d);
    double L2 = std::min(std::max(L_D,(L_L + L_H)/tan(theta)), L_D + L_d) - L1;
    if (L_L/tan(theta) >= L_D+L_d) // theta too small, pseudoscalar flying too lowly
        return 0;
    if ((L_L + L_H)/tan(theta) <= L_D) //theta too large, pseudoscalar flying too highly
        return 0; 
    // The probability that the pseudoscalar would decay in the detector is then as follows (Decay law:            
    return  exp(-L1/(beta_z*gamma*ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb))) * (1. - exp(-L2/(beta_z*gamma*ctau(mX,g_TT_Lamb,g_MM_Lamb,g_EE_Lamb,g_GG_Lamb)))); 
}


