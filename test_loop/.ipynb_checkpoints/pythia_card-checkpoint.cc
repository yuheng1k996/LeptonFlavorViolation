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
    g_TM_L = 0;
    g_TE_L = 0;
    g_ME_L = 0;
    g_TT_L = 0;
    g_MM_L = 0;
    g_EE_L = 0;
    reallyProducedX = 0;//total number of X produced at the real experiment
};    


bool pythia_card::doCalculations() {
    number_total_tau_tau = 4.6e10;//total number of ee->tautau events at Belle II with 50/ab
    double BrTau2OneProng = 0.85;//PDG
  
    reallyProducedX = number_total_tau_tau * (2 * BrTau2OneProng * (BRtau2xmu(mX, g_TM_L, g_TE_L) + BRtau2xe(mX, g_TM_L, g_TE_L)));

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
        //pythia->readString("NewGaugeBoson:ffbar2gmZZprime = on");
        //pythia->readString("Zprime:gmZmode = 1");
        //pythia->readString("22:oneChannel =  1 1 103 15 -15 ");//
        //pythia->readString("32:oneChannel =  1 1 103 15 -15 ");//

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
        pythia->readString("36:tau0 = "+floatToString(ctau(mX,g_TT_L,g_MM_L,g_EE_L)*1.E3));
        //x -> l lbar
        pythia->readString("36:oneChannel =  1 "+floatToString(BRx2tautau(mX,g_TT_L,g_MM_L,g_EE_L))+ " 103 15 -15 ");
        pythia->readString("36:addChannel =  1 "+floatToString(BRx2mumu(mX,g_TT_L,g_MM_L,g_EE_L))+ " 103 13 -13 ");
        pythia->readString("36:addChannel =  1 "+floatToString(BRx2ee(mX,g_TT_L,g_MM_L,g_EE_L))+ " 103 11 -11 ");
        pythia->readString("36:addChannel =  1 "+floatToString(BRx2gmgm(mX,g_TT_L,g_MM_L,g_EE_L))+ " 103 22 22 ");//fermion loop diagram

	    //set tau decay channels
        pythia->readString("15:oneChannel =  1 "+floatToString(BRtau2xmu(mX,g_TM_L,g_TE_L))+" 103 36 13 ");//tau ->x mu
        pythia->readString("15:addChannel =  1 "+floatToString(BRtau2xe(mX,g_TM_L,g_TE_L))+" 103 36 11 ");//tau ->x e
	
        
        //allow for very long-lived particles
        pythia->readString("ResonanceWidths:minWidth = 1.97e-30");


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
            pythia->readString("Next:numberCount = 100000");
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

    double producedtau{};
    double producedX{};
    double observedX{};
    
    double xForwardSphere{};
    double xBackwardSphere{};
    
    double n_prompt_smallr{};
    double n_prompt_d0z0{};
    double n_prompt_baselineeff{};

    double n_disp_fidvol{};
    double n_disp_baselineeff{};
    double n_disp_DispTrackEff{};
    
    double n_mod_DispTrackEffee{};
    double n_mod_DispTrackEffmumu{};
    double n_mod_DispTrackEff{};

    try{
        //creating files to scan
        //ofstream myfile;
        //myfile.open ("./scanning/data_mX_" +floatToString(mX)+ "_g_TM_" +floatToString(g_TM_L)+ "_g_TE_" +floatToString(g_TE_L)+ "_g_MM_" +floatToString(g_MM_L)+ "_g_EE_" +floatToString(g_EE_L)+ ".txt");
        
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
			{	//count X
				producedX += 1;
                        
                    		//1) prompt cuts
				//cut 1.1) enough hits in the vertex detector, realized by requiring small r
				double rDec = sqrt(pow(pythia->event[i].xDec(),2)+pow(pythia->event[i].yDec(),2));
                        	if (  rDec < 100. )//in mm
                {
					n_prompt_smallr += 1.;
					
					//cut 1.2) small d0 and z0 for the llbar produced from x decay
                    			bool l1 = d0z0constraints(pythia->event[pythia->event[i].daughter1()]);
                    			bool l2 = d0z0constraints(pythia->event[pythia->event[i].daughter2()]);
					if (l1*l2){
						n_prompt_d0z0 += 1.;
						n_prompt_baselineeff += 1.*baselineefficiency(i);
                    }
                }
                    
                        
                    
                    //myfile << sqrt(pow(x1Prod,2)+pow(y1Prod,2))  << " " << z1Prod << " " <<   d1 <<  " " << z1 << " " << d2 << " " << z2 << '\n';
                 

                        
                //2) displaced cuts
				//cut 2.1) fiducial volume
                     		if(decayingInsideFidVol(pythia->event[i])){
					n_disp_fidvol += 1.;
					//cut 2.2) detector efficiency
					n_disp_baselineeff += baselineefficiency(i);
					n_disp_DispTrackEff += baselineefficiency(i) * DispTrackEff(pythia->event[i]);
                }


				//check number of x traveling forward or backward
                 		if (pythia->event[i].p().eta()>0)xForwardSphere+=1;
				if (pythia->event[i].p().eta()<0)xBackwardSphere+=1;
			}
                    
                    
                //3.1) modified exponential functions tau -> xe
                    if(abs(pythia->event[pythia->event[i].daughter1()].id()) == 36 && abs(pythia->event[pythia->event[i].daughter2()].id()) == 11){
				double p1 = decayProbabilityBelle2Part1(pythia->event[pythia->event[i].daughter1()]);
				double p2 = decayProbabilityBelle2Part2(pythia->event[pythia->event[i].daughter1()]);
                        	observedX += p1+p2;
                            if (p2 < 0){
                                std::cout <<  "pythia->event[pythia->event[i].daughter2()].id(): " << pythia->event[pythia->event[i].daughter2()].id() << '\n';
                            }
                        //baselineeff (tau -> eee) and (tau -> emumu)
					    n_mod_DispTrackEffee += (p1+p2)*(0.060*BRx2ee(mX,g_TT_L,g_MM_L,g_EE_L)+0.061*BRx2mumu(mX,g_TT_L,g_MM_L,g_EE_L));

                    }
                //3.2) modified exponential functions tau -> xmu
                    if(abs(pythia->event[pythia->event[i].daughter1()].id()) == 36 && abs(pythia->event[pythia->event[i].daughter2()].id()) == 13){
				double p1 = decayProbabilityBelle2Part1(pythia->event[pythia->event[i].daughter1()]);
				double p2 = decayProbabilityBelle2Part2(pythia->event[pythia->event[i].daughter1()]);
                        	observedX += p1+p2;
                        //baselineeff (tau -> muee) and (tau -> mumumu)
					    n_mod_DispTrackEffmumu += (p1+p2)*(0.093*BRx2ee(mX,g_TT_L,g_MM_L,g_EE_L)+0.076*BRx2mumu(mX,g_TT_L,g_MM_L,g_EE_L));

                    }
            }
        }
        if(verbose)
            pythia->stat();
        //close the scanning files
        //myfile << producedX << endl;
        //myfile.close();
    }
    catch(std::exception& e) {
        std::cerr << "!!! Error occured while trying to run Pythia: " << e.what() << std::endl;
        return false;
    }


    double Prompt_reallyobservedX = n_prompt_smallr / producedX * reallyProducedX;
    double Prompt_reallyvisibleX = n_prompt_baselineeff * BRx2Visibles(mX,g_TT_L,g_MM_L,g_EE_L) / producedX * reallyProducedX;
    double Displaced_reallyobservedX = n_disp_fidvol / producedX * reallyProducedX;
    double Displaced_reallyvisibleX = n_disp_DispTrackEff * BRx2Visibles(mX,g_TT_L,g_MM_L,g_EE_L) / producedX * reallyProducedX;
    
    n_mod_DispTrackEff = n_mod_DispTrackEffee + n_mod_DispTrackEffmumu;
    double Modified_reallyobservedX = n_mod_DispTrackEff / producedX * reallyProducedX;
    double Modified_reallyvisibleX = n_mod_DispTrackEff * BRx2Visibles(mX,g_TT_L,g_MM_L,g_EE_L) / producedX * reallyProducedX;


    // Results
    std::cout << "Gammatau2xmu: " << Gammatau2xmu(mX,g_TM_L) << '\n';
    std::cout << "Gammatau2xe: " << Gammatau2xe(mX,g_TE_L) << '\n';
    std::cout << '\n';
    std::cout << "BRtau2xmu: " << BRtau2xmu(mX,g_TM_L,g_TE_L) << '\n';
    std::cout << "BRtau2xe: " << BRtau2xe(mX,g_TM_L,g_TE_L) << '\n';
    std::cout << "BRtau2xmu + BRtau2xe: " << BRtau2xmu(mX,g_TM_L,g_TE_L) + BRtau2xe(mX,g_TM_L,g_TE_L) << '\n';
    std::cout << '\n';
    std::cout << "Total Gamma of X [GeV]: " << totalGammaX(mX,g_TT_L,g_MM_L,g_EE_L) << '\n';
    std::cout << "ctau [m]: " << ctau(mX,g_TT_L,g_MM_L,g_EE_L) << '\n';
    std::cout << '\n';
    std::cout << "loopB1(tau,mX): " << widthCalculator::loopB1(widthCalculator::mtau,mX) << '\n';
    std::cout << "loopB1(mu ,mX): " << widthCalculator::loopB1(widthCalculator::mmu, mX) << '\n';
    std::cout << "loopB1(e  ,mX): " << widthCalculator::loopB1(widthCalculator::me,  mX) << '\n';
    std::cout << '\n';
    std::cout << "Gammax2tautau: " << Gammax2tautau(mX,g_TT_L) << '\n';
    std::cout << "Gammax2mumu: " << Gammax2mumu(mX,g_MM_L) << '\n';
    std::cout << "Gammax2ee: " << Gammax2ee(mX,g_EE_L) << '\n';
    std::cout << "Gammax2gmgm: " << Gammax2gmgm(mX,g_TT_L,g_MM_L,g_EE_L) << '\n';
    std::cout << '\n';
    std::cout << "BRx2tautau: " << BRx2tautau(mX,g_TT_L,g_MM_L,g_EE_L) << '\n';
    std::cout << "BRx2mumu: " << BRx2mumu(mX,g_TT_L,g_MM_L,g_EE_L) << '\n';
    std::cout << "BRx2ee: " << BRx2ee(mX,g_TT_L,g_MM_L,g_EE_L) << '\n';
    std::cout << "BRx2gmgm: " << BRx2gmgm(mX,g_TT_L,g_MM_L,g_EE_L) << '\n';
    std::cout << '\n';
    std::cout << "Proportion of X flying forward (eta>0): " << xForwardSphere / producedX << '\n';
    std::cout << "Proportion of X flying backward (eta<0): " << xBackwardSphere / producedX << '\n';
    std::cout << '\n';
    std::cout << "produced tau excluding those daughters are also tau lepton: " << producedtau << '\n';
    std::cout << "produced X: " << producedX << '\n';   
    std::cout << "produced X/NMC: " << producedX/double(nEventsMC) << '\n';
    std::cout << '\n';
    std::cout << "Prompt analysis results: \n";
    std::cout << "n_prompt_smallr: " << n_prompt_smallr << '\n';
    std::cout << "n_prompt_d0z0: " << n_prompt_d0z0 << '\n';
    std::cout << "n_prompt_baselineeff: " << n_prompt_baselineeff << '\n';
    std::cout << '\n';
    std::cout << "prompt small_r efficiencies (ind. and cum.): " << n_prompt_smallr/producedX << " " << n_prompt_smallr/producedX  << '\n';
    std::cout << "prompt d0z0 efficiencies (ind. and cum.): " << n_prompt_d0z0/n_prompt_smallr << " " << n_prompt_d0z0/producedX << '\n';
    std::cout << "prompt baseline efficiencies (ind. and cum.): " << n_prompt_baselineeff/n_prompt_d0z0 << " " << n_prompt_baselineeff/producedX << '\n';
    std::cout << '\n';
    std::cout << "Displaced analysis results: \n";
    std::cout << "n_disp_fidvol: " << n_disp_fidvol << '\n';
    std::cout << "n_disp_baselineeff: " << n_disp_baselineeff << '\n';
    std::cout << "n_disp_DispTrackEff: " << n_disp_DispTrackEff << '\n';
    std::cout << '\n';
    std::cout << "observedX: " << observedX <<'\n';
    std::cout << "n_mod_DispTrackEff: " << n_mod_DispTrackEff <<'\n';
    std::cout << '\n';
    std::cout << "displaced fiducial volume efficiencies (ind. and cum.): " << n_disp_fidvol/producedX << " " << n_disp_fidvol/producedX  << '\n';
    std::cout << "displaced baseline efficiencies (ind. and cum.): " << n_disp_baselineeff/n_disp_fidvol << " " << n_disp_baselineeff/producedX << '\n';
    std::cout << "displaced displaced-tracking efficiencies (ind. and cum.): " << n_disp_DispTrackEff/n_disp_baselineeff << " " << n_disp_DispTrackEff/producedX << '\n';
    std::cout << '\n';
    std::cout << "modified efficiencies: " << n_disp_DispTrackEff/producedX << '\n';
    std::cout << '\n';
    std::cout << "Expected total signal numbers at Belle II with 50/ab:	 \n";
    std::cout << "total number of tau-tau events: " << number_total_tau_tau << '\n';
    std::cout << "reallyProduced tau: " << 2*number_total_tau_tau << '\n';
    std::cout << "reallyProducedX: " << reallyProducedX << '\n';
    std::cout << "Prompt_reallyobservedX: " << Prompt_reallyobservedX << '\n';
    std::cout << "Prompt_reallyvisibleX: " << Prompt_reallyvisibleX << '\n';
    std::cout << "Displaced_reallyobservedX: " << Displaced_reallyobservedX << '\n';
    std::cout << "Displaced_reallyvisibleX: " << Displaced_reallyvisibleX << '\n';
    std::cout << '\n';
    std::cout << "Modified_reallyobservedX: " << Modified_reallyobservedX << '\n';
    std::cout << "Modified_reallyvisibleX: " << Modified_reallyvisibleX << '\n';
    
    return true;
}

																																		

//lepton flavor violation

//physical quantities

double pythia_card::Gammatau2xmu(double mX, double g_TM_L)
{
	return widthCalculator::violationLepWidth(mX, widthCalculator::mtau, widthCalculator::mmu, g_TM_L);
}

double pythia_card::Gammatau2xe(double mX, double g_TE_L)
{
	return widthCalculator::violationLepWidth(mX, widthCalculator::mtau, widthCalculator::me, g_TE_L);
}

double pythia_card::NewTotalGammatau(double mX, double g_TM_L, double g_TE_L)
{
	return widthCalculator::tauSMGamma + Gammatau2xmu(mX, g_TM_L) + Gammatau2xe(mX, g_TE_L);
}

double pythia_card::BRtau2xmu(double mX, double g_TM_L, double g_TE_L)
{
	return Gammatau2xmu(mX, g_TM_L) / NewTotalGammatau(mX, g_TM_L, g_TE_L);
}

double pythia_card::BRtau2xe(double mX, double g_TM_L, double g_TE_L)
{
	return Gammatau2xe(mX, g_TE_L) / NewTotalGammatau(mX, g_TM_L, g_TE_L);
}



double pythia_card::Gammax2tautau(double mX, double g_TT_L)
{
	return widthCalculator::pseudoscalarLepWidth(mX,widthCalculator::mtau,g_TT_L);
}

double pythia_card::Gammax2mumu(double mX, double g_MM_L)
{
	return widthCalculator::pseudoscalarLepWidth(mX,widthCalculator::mmu,g_MM_L);
}

double pythia_card::Gammax2ee(double mX, double g_EE_L)
{
	return widthCalculator::pseudoscalarLepWidth(mX,widthCalculator::me,g_EE_L);
}

double pythia_card::Gammax2gmgm(double mX, double g_TT_L, double g_MM_L, double g_EE_L)
{
	return widthCalculator::pseudoscalarGmWidth(mX, widthCalculator::mtau, widthCalculator::mmu, widthCalculator::me, g_TT_L, g_MM_L, g_EE_L);
}

double pythia_card::totalGammaX(double mX, double g_TT_L, double g_MM_L, double g_EE_L)
{
	return Gammax2tautau(mX,g_TT_L) + Gammax2mumu(mX,g_MM_L) + Gammax2ee(mX,g_EE_L) + Gammax2gmgm(mX,g_TT_L,g_MM_L,g_EE_L);
}


double pythia_card::ctau(double mX, double g_TT_L, double g_MM_L, double g_EE_L)//ctau in meter
{
	double chbar = 1.973269631e-16;//GeV m
	return chbar/ totalGammaX(mX,g_TT_L,g_MM_L,g_EE_L);
}


double pythia_card::BRx2tautau(double mX, double g_TT_L, double g_MM_L, double g_EE_L)
{
	return Gammax2tautau(mX,g_TT_L)/totalGammaX(mX,g_TT_L,g_MM_L,g_EE_L);
}

double pythia_card::BRx2mumu(double mX, double g_TT_L, double g_MM_L, double g_EE_L)
{
	return Gammax2mumu(mX,g_MM_L)/totalGammaX(mX,g_TT_L,g_MM_L,g_EE_L);
}

double pythia_card::BRx2ee(double mX, double g_TT_L, double g_MM_L, double g_EE_L)
{
	return Gammax2ee(mX,g_EE_L)/totalGammaX(mX,g_TT_L,g_MM_L,g_EE_L);
}

double pythia_card::BRx2gmgm(double mX, double g_TT_L, double g_MM_L, double g_EE_L)
{
	return Gammax2gmgm(mX,g_TT_L,g_MM_L,g_EE_L)/totalGammaX(mX,g_TT_L,g_MM_L,g_EE_L);
}

double pythia_card::BRx2Visibles(double mX, double g_TT_L, double g_MM_L, double g_EE_L)
{
	double BRx2Visibles=1.-BRx2gmgm(mX,g_TT_L,g_MM_L,g_EE_L);
	return   BRx2Visibles;
}


																																		
//functions of cuts


bool pythia_card::d0z0constraints(Pythia8::Particle XXX)
{
    double z0 = abs(XXX.zProd() - ((XXX.xProd()*XXX.px()+XXX.yProd()*XXX.py())*XXX.pz())/(XXX.pT2()));
    double rT = sqrt(pow(XXX.xProd(),2)+pow(XXX.yProd(),2));
    double d0 = abs(XXX.px()*XXX.yProd()-XXX.py()*XXX.xProd()) / XXX.pT();
    double rC = (XXX.pT()/0.45)*1.E3;//in mm
    d0 = sqrt(pow(rC+d0,2)+pow(rT,2)-pow(d0,2))-rC;
    
    //std::cout << "d0_" << d0 << "_z0_" << z0 << '\n';
    //require d0 < 0.5cm and z0 < 3cm
    if (d0 < 5 && z0 < 30){
        return true;
    }
    else {return false;}
}

bool pythia_card::decayingInsideFidVol(Pythia8::Particle XXX)
{
    constexpr double r_min = 10.;//in mm
    constexpr double r_max = 800.; 
    constexpr double z_min = -400.;
    constexpr double z_max = 1200.;

    double rDec = sqrt(pow(XXX.xDec(),2)+pow(XXX.yDec(),2));
    double zDec = XXX.zDec();
    
    if (rDec > r_min and rDec < r_max and zDec > z_min and zDec < z_max){
    return true;
    }
    else {return false;}
}


double pythia_card::DispTrackEff(Pythia8::Particle XXX)//detector efficiency
{
    //fit the efficiency linearly (in mm)
    double rDec = sqrt(pow(XXX.xDec(),2)+pow(XXX.yDec(),2));
    double rmin= 10.;
    double rmax= 800.;
    if(rDec<rmin || rDec>rmax){//for r smaller than rmin or larger than rmax, no detection
	    return 0.;
    }
    else{// for r between rmin and rmax, use the linear interpolation (eff=1 if r=10, eff=0 if r=800)
        return 1*(rDec-rmax)/(rmin-rmax);
    }
}



double pythia_card::baselineefficiency(int i)//i is the particle index of X in an event
{
	double epsilon{};
	if ( abs(pythia->event[i+1].id()) == 11  )//tau->X+e
	{
		if (abs(pythia->event[pythia->event[i].daughter1()].id())==11) //tau->X+e, X->ee
		{
			epsilon=0.060;
		}
		else if (abs(pythia->event[pythia->event[i].daughter1()].id())==13) //tau->X+e, X->mu mu
		{
			epsilon=0.061;
		}
	}
	else if ( abs(pythia->event[i+1].id()) == 13   )//tau->X+mu
	{
		if (abs(pythia->event[pythia->event[i].daughter1()].id())==11) //tau->X+mu, X->ee
		{
			epsilon=0.093;
		}
		else if (abs(pythia->event[pythia->event[i].daughter1()].id())==13) //tau->X+mu, X->mu mu
		{
			epsilon=0.076;
		}
	}
	return epsilon;
}


																																		
//calculate the efficiency in an exponential integral

double pythia_card::decayProbabilityBelle2Part1(Pythia8::Particle XXX) {//Part1 is a symmetric ATLAS-like detector
    constexpr double R_I = 0.01;//in meter
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
     //return exp(-L1/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (1. - exp(-L2/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))));
     
    // the efficiency is a function (1*(l-l_max)/(l_min-l_max)), l_min=0.01, eff=1, l_max=0.8, eff=0 ; therefore, for every dL, the modified dacay probability is
    //exp(-L/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))*(1/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)))*(1*(L*tan(theta)-0.8)/(0.01-0.8))*dL,from L1 to L1+L2
    
    return exp(-L1/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (1.01266 - 1.26582*fabs(tan(theta))*(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)+ L1)) + exp(-(L1+L2)/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (-1.01266 + 1.26582*fabs(tan(theta))*(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)+ L2+L1));
    //note that r should be non-negative, so it is necessary to take the absolute value of tan(theta)
}


double pythia_card::decayProbabilityBelle2Part2(Pythia8::Particle XXX) {//Part2 is AL3X-like
    constexpr double L_D = 0.4;//horizontal distance from the IP to the detector
    constexpr double L_L = 0.01;//vertical distance from the IP to the detector  or  equivalently the inner radius
    constexpr double L_d = 0.8;// length of the detector: 120 - 40 = 80 cm
    constexpr double L_H = 0.79;// transverse length = 80 - 1 = 79 cm

    // Identify the kinematic properties of the neutralino
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
    double efficiency = (20*(-4 + 5*(L1 + L2 + (beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)))*tan(theta) + 
       exp(L2/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)))*(4 - 5*(L1 + (beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)))*tan(theta))))/
   (79.*exp((L1 + L2)/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))));
        
        
        //exp(-L1/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (1.01266 - 1.26582*fabs(tan(theta))*(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)+ L1)) + exp(-(L1+L2)/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (-1.01266 + 1.26582*fabs(tan(theta))*(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)+ L2+L1));
        //exp(-L1/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (1. - exp(-L2/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))));
    if (efficiency < 0 ){
        std::cout << "L1: "<< L1 << ", L2: " << L2  <<", ctau: "<<  ctau(mX, g_TT_L, g_MM_L, g_EE_L) << ", beta_z: " << beta_z << ", gamma: "<< gamma << ", tantheta:" << tan(theta) <<", efficiency: " << std::setprecision(20) << -efficiency << '\n';
    }
     return efficiency; //  not work anymore!!!
    
    // The probability that the pseudoscalar would decay in the detector is then as follows (Decay law:            
    //return  exp(-L1/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (1. - exp(-L2/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)))); 
   
      //return exp(-L1/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (1.01266 - 1.26582*fabs(tan(theta))*(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)+ L1)) + exp(-(L1+L2)/(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L))) * (-1.01266 + 1.26582*fabs(tan(theta))*(beta_z*gamma*ctau(mX, g_TT_L, g_MM_L, g_EE_L)+ L2+L1));
   
}