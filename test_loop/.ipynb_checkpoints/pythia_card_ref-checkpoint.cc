#include "pythia_card.h"

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
    mN1 = 0;  //mass of N1 in GeV
    lambOverMSq = 0; //lambOverMSq in GeV^-2
    number_total_tau_tau = 0;
    reallyProducedN1 = 0;//total number of N1 produced at the real experiment
};    


bool pythia_card::doCalculations() {
    number_total_tau_tau=4.6e10;//total number of ee->tautau events at Belle II with 50/ab
    double BrTau2OneProng = 0.85;//PDG
    reallyProducedN1 = number_total_tau_tau * 2 * BrTau2OneProng *  (   BRTautoNpi(mN1,lambOverMSq) +   BRTautoNrho(mN1,lambOverMSq) +  BRTautoNK(mN1,lambOverMSq) + BRTautoNvK(mN1,lambOverMSq) +BRTautoe(mN1,lambOverMSq)  +BRTautomu(mN1,lambOverMSq) )    ;
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

        pythia->readString("WeakSingleBoson:ffbar2ffbar(s:gm) = on"); //off-shell photon    -> 15% efficiency (N_{tautau}/NMC~15%) but fast.

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

        //add a new fermion particle
        pythia->readString("1000022:new = n1");//no anti-name so that n1-> Majorana
        pythia->readString("1000022:spinType = 2");//fermion
        pythia->readString("1000022:chargeType = 0");
        pythia->readString("1000022:colType = 0");
        pythia->readString("1000022:m0 = "+floatToString(mN1)); // pythia masses are in GeV


	//set n1 decay 
	 pythia->readString("1000022:mayDecay = 1 ");
        pythia->readString("1000022:tau0 = "+floatToString(ctau(mN1,lambOverMSq)*1.E3)); //converting ctau in meter to mm by the factor 1.E3
        pythia->readString("1000022:oneChannel =  1 "+floatToString(BRNtopi0(mN1,lambOverMSq))+ " 103 16 111 ");//n -> pi0 nu
        pythia->readString("1000022:addChannel =  1 "+floatToString(BRNtorho0(mN1,lambOverMSq))+ " 103 16 113 ");//n -> rho0 nu
        pythia->readString("1000022:addChannel =  1 "+floatToString(BRNtoeta(mN1,lambOverMSq))+ " 103 16 221 ");//n -> eta nu
        pythia->readString("1000022:addChannel =  1 "+floatToString(BRNtoetap(mN1,lambOverMSq))+ " 103 16 331 ");//n1 -> eta' nu
        pythia->readString("1000022:addChannel =  1 "+floatToString(BRNtoomega(mN1,lambOverMSq))+ " 103 16 223 ");//n1 -> omega nu
        pythia->readString("1000022:addChannel =  1 "+floatToString(BRNtophi(mN1,lambOverMSq))+ " 103 16 333 ");//n1 -> phi nu



	//set tau decay channels
	pythia->readString("15:oneChannel =  1 "+floatToString(BRTautoNpi(mN1,lambOverMSq))+" 103 1000022 -211 ");//tau- ->n pi-
	pythia->readString("15:addChannel =  1 "+floatToString(BRTautoNrho(mN1,lambOverMSq))+" 103 1000022 -213 ");//tau- ->n rho-
	pythia->readString("15:oneChannel =  1 "+floatToString(BRTautoNK(mN1,lambOverMSq))+" 103 1000022 -321 ");//tau- ->n kaon-
	pythia->readString("15:addChannel =  1 "+floatToString(BRTautoNvK(mN1,lambOverMSq))+" 103 1000022 -323 ");//tau- ->n v kaon-
	pythia->readString("15:oneChannel =  1 "+floatToString(BRTautoe(mN1,lambOverMSq))+" 103 1000022 11 -12 ");//tau- ->n e nue
	pythia->readString("15:addChannel =  1 "+floatToString(BRTautomu(mN1,lambOverMSq))+" 103 1000022 13 -14 ");//tau- ->n muon nu_mu
	
	
	
	//pythia->readString("15:addChannel =  1 "+floatToString(BRtauTon1a1m(mN1,lambOverMSq))+" 103 1000022 -20213 ");//tau- ->n1 a1-    //no need to include a1 as it is kinematicaly irrelevant
  

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

    double producedtau{};
    double producedN1{};
    double observedN1{};
    double visibleN1{};
    double decayN{};
 

    int insideFidVol{};

    int n1ForwardSphere{};
    int n1BackwardSphere{};
    try{
        for (int iEvent = 0; iEvent < nEventsMC; ++iEvent) {
                if (!pythia->next()) continue;
                // Check the list of final state particles
                for (int i = 0; i < pythia->event.size(); ++i)
		{
			//count tau leptons, (whose daughters are not tau leptons)
                	if (abs(pythia->event[i].id()) == 15 && abs(pythia->event[pythia->event[i].daughter1()].id())!=15 )
			{
                        	producedtau += 1;
			}
                	if (abs(pythia->event[i].id()) == 1000022)
			{	//count n1
				producedN1 += 1;
				double p1 = decayProbabilityBelle2Part1(pythia->event[i]);
				double p2 = decayProbabilityBelle2Part2(pythia->event[i]);
				double p3 = decayProbabilityBelle2Part3(pythia->event[i]);
				double p4 = decayProbabilityBelle2Part4(pythia->event[i]);
				decayN +=p3+p4;
                        	observedN1 += p1+p2;
				visibleN1 += (p1+p2)*BRn1ToVisibles(mN1,lambOverMSq);
				if (decayingInsideFidVol(pythia->event[i])) insideFidVol+=1;
				if (pythia->event[i].p().eta()>0)n1ForwardSphere+=1;
				if (pythia->event[i].p().eta()<0)n1BackwardSphere+=1;
			}
                }
        }
        if(verbose)
            pythia->stat();   
    }
    catch(std::exception& e) {
        std::cerr << "!!! Error occured while trying to run Pythia: " << e.what() << std::endl;
        return false;
    }

    double reallyobservedN1 = observedN1/double(producedN1) * reallyProducedN1;			 // total N decay insider the chamber with linear efficiency 
    double reallyDecayN = decayN/double(producedN1) * reallyProducedN1;                              //total N decay insidide the chamber
    double observedFiducialEfficiencyN1 = observedN1/double(producedN1);

    double reallyvisibleN1 = visibleN1/double(producedN1) * reallyProducedN1;
    double visibleFiducialEfficiencyN1 = visibleN1/double(producedN1);

    // Results
    std::cout << "TautoNpi: " << TautoNpi(mN1,lambOverMSq) << '\n';
    std::cout << "TautoNrho: " << TautoNrho(mN1,lambOverMSq) << '\n';
    std::cout << "TautoNK: " << TautoNK(mN1,lambOverMSq) << '\n';
    std::cout << "TautoNvK: " << TautoNvK(mN1,lambOverMSq) << '\n';
    std::cout << "Tautoe: " << Tautoe(mN1,lambOverMSq) << '\n';
    std::cout << "Tautomu: " << Tautomu(mN1,lambOverMSq) << '\n';
    std::cout << "totalTau: " << totalTau(mN1,lambOverMSq) << '\n';
    
    std::cout << '\n';
    std::cout << "BRTautoNpi: " << BRTautoNpi(mN1,lambOverMSq) << '\n';
    std::cout << "BRTautoNrho: " << BRTautoNrho(mN1,lambOverMSq) << '\n';
    std::cout << "BRTautoNK: " << BRTautoNK(mN1,lambOverMSq) << '\n';
    std::cout << "BRTautoNvK: " << BRTautoNvK(mN1,lambOverMSq) << '\n';
    std::cout << "BRTautoe: " << BRTautoe(mN1,lambOverMSq) << '\n';
    std::cout << "BRTautomu: " << BRTautomu(mN1,lambOverMSq) << '\n';
    std::cout << '\n';
    std::cout << "produced tau excluding those whose daughters are also tau lepton: " << producedtau << '\n';
    std::cout << "produced N1: " << producedN1 << '\n';   
    std::cout << "produced N1/NMC: " << producedN1/double(nEventsMC) << '\n';
    std::cout << std::endl;
    std::cout << "Proportion of N1 flying forward (eta>0): " << n1ForwardSphere / producedN1 << '\n';
    std::cout << "Proportion of N1 flying backward (eta<0): " << n1BackwardSphere / producedN1 << '\n';
    std::cout << '\n';
    std::cout << "Total Gamma [GeV]: " << totalN(mN1,lambOverMSq) << '\n';
    std::cout << "ctau [m]: " << ctau(mN1,lambOverMSq) << '\n';
    std::cout << '\n';
    std::cout << "Ntopi0: " << Ntopi0(mN1, lambOverMSq) << '\n';
    std::cout << "Ntoeta: " << Ntoeta(mN1, lambOverMSq) << '\n';
    std::cout << "Ntoetap: " << Ntoetap(mN1, lambOverMSq) << '\n';
    std::cout << "Ntorho0: " << Ntorho0(mN1, lambOverMSq) << '\n';
    std::cout << "Ntoomega: " << Ntoomega(mN1, lambOverMSq) << '\n';
    std::cout << "Ntophi: " << Ntophi(mN1, lambOverMSq) << '\n';
    std::cout << "Ntoll: " << Ntoll(mN1, lambOverMSq) << '\n';
    std::cout << "Ntonunu: " << Ntonunu(mN1, lambOverMSq) << '\n';
    std::cout << '\n';
    std::cout << "BRNtopi0: " << BRNtopi0(mN1, lambOverMSq) << '\n';
    std::cout << "BRNtoeta: " << BRNtoeta(mN1, lambOverMSq) << '\n';
    std::cout << "BRNtoetap: " << BRNtoetap(mN1, lambOverMSq) << '\n';
    std::cout << "BRNtorho0: " << BRNtorho0(mN1, lambOverMSq) << '\n';
    std::cout << "BRNtoomega: " << BRNtoomega(mN1, lambOverMSq) << '\n';
    std::cout << "BRNtophi: " << BRNtophi(mN1, lambOverMSq) << '\n';
        std::cout << "BRNtoll: " << BRNtoll(mN1, lambOverMSq) << '\n';
    std::cout << "BRNtonunu: " << BRNtonunu(mN1, lambOverMSq) << '\n';
    std::cout << '\n';
    //std::cout << "insideFidVol: " << insideFidVol << '\n';    //needed only when we decay n1 in py8
    std::cout << "observedN: " << observedN1 << '\n';
    std::cout << "detector effi: " << "10%" << '\n';
    std::cout << "visibleN: " << visibleN1 << '\n';
    std::cout << '\n';
    std::cout << "Observed Fiducial efficiency: " << observedFiducialEfficiencyN1 << '\n';
    std::cout << "Visible Fiducial efficiency: " << visibleFiducialEfficiencyN1 << '\n';
    std::cout << '\n';
    std::cout << "reallyProduced tau: " << 2*number_total_tau_tau << '\n';
    std::cout << "total number of tau-tau events: " << number_total_tau_tau << '\n';
    std::cout << "Br(tau -> leptons): " << "35.21%" << '\n';
    std::cout << '\n';
    std::cout << "reallyProducedN1: " << reallyProducedN1 << '\n';
    std::cout << "reallyobservedN1: " << reallyobservedN1 << '\n';
   std::cout << "reallyDecayN: " << reallyDecayN << '\n';
    std::cout << "reallyvisibleN1: " << reallyvisibleN1 << '\n';
    std::cout << '\n';

    return true;
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

												//    double TautoP(double mN, double mtau, double mMeson, double fMeson,double mixingsq, double VUD);            

double pythia_card::TautoNpi(double mN1, double lambOverMSq)
{
	return widthCalculator::TautoP(mN1, widthCalculator::mtau, widthCalculator::mpip,  widthCalculator::fpip, lambOverMSq, widthCalculator::Vud);
}

double pythia_card::TautoNK(double mN1, double lambOverMSq)
{
	return widthCalculator::TautoP(mN1, widthCalculator::mtau, widthCalculator::mK,  widthCalculator::fK, lambOverMSq, widthCalculator::Vus) ;
}



												//  double TautoV(double mN, double mtau, double mMeson, double fMeson, double mixingsq, double VUD);


double pythia_card::TautoNrho(double mN1, double lambOverMSq)
{
	return widthCalculator::TautoV(mN1, widthCalculator::mtau, widthCalculator::mrho,  widthCalculator::frho, lambOverMSq, widthCalculator::Vud);
	
}

double pythia_card::TautoNvK(double mN1, double lambOverMSq)
{
	return widthCalculator::TautoV(mN1, widthCalculator::mtau, widthCalculator::mvK,  widthCalculator::fvK, lambOverMSq, widthCalculator::Vus)  ;
	
}

//tau to e 
double pythia_card::Tautoe(double x, double SinThetaSq)
{
	if (x <= 1.7 and x >= 0 ) return 0.001*(4.04858494948796e-13 +3.734086406832738e-16* x -1.0373545234574142e-12*pow(x,2)     + 2.0655301820395976e-13*pow(x,3) + 1.8058287474510143e-12*pow(x,4) - 
  2.732699307841911e-12*pow(x,5) + 2.516449794869419e-12*pow(x,6) - 
  1.8362565836029536e-12*pow(x,7) + 1.009010637026267e-12*pow(x,8) - 
   3.709073485459008e-13*pow(x,9) + 8.001383201273446e-14*pow(x,10)- 
   7.6692720872473e-15*pow(x,11)) * SinThetaSq ;
	else return 0;
}

//tau to muon 
double pythia_card::Tautomu(double x, double lambOverMSq)
{
	if (x <= 1.65 and x >= 0 ) return 0.001*(3.9368114930019283e-13 +3.4236155788150233e-16* x -1.0257872375875665e-12*pow(x,2)     + 2.00636149753519e-13*pow(x,3) + 1.8306065358548687e-12*pow(x,4) - 
   2.813169465733193e-12*pow(x,5) + 2.660953184199404e-12*pow(x,6) - 
   2.002982147402113e-12*pow(x,7) + 1.133297105220641e-12*pow(x,8) - 
   4.2853497866026385e-13*pow(x,9) + 9.512455852291831e-14*pow(x,10)- 
   9.380273795033076e-15*pow(x,11)) * lambOverMSq ;
	else return 0;
}









double pythia_card::totalTau(double mN1,double lambOverMSq)
{
	return widthCalculator::tauSMGamma+ TautoNpi(mN1,lambOverMSq)+ TautoNrho(mN1,lambOverMSq) +Tautoe( mN1,  lambOverMSq)+ Tautomu(mN1,lambOverMSq)+TautoNvK(mN1,lambOverMSq) +  TautoNK(mN1,lambOverMSq) ;
}




double pythia_card::BRTautoNpi(double mN1,double lambOverMSq)
{
	return TautoNpi(mN1,lambOverMSq) / totalTau(mN1,lambOverMSq);
}

double pythia_card::BRTautoNrho(double mN1,double lambOverMSq)
{
	return TautoNrho(mN1,lambOverMSq) / totalTau(mN1,lambOverMSq);
}

double pythia_card::BRTautoNK(double mN1,double lambOverMSq)
{
	return TautoNK(mN1,lambOverMSq) / totalTau(mN1,lambOverMSq);
}

double pythia_card::BRTautoNvK(double mN1,double lambOverMSq)
{
	return TautoNvK(mN1,lambOverMSq) / totalTau(mN1,lambOverMSq);
}
double pythia_card::BRTautoe(double mN1,double lambOverMSq)
{
	return Tautoe(mN1,lambOverMSq) / totalTau(mN1,lambOverMSq);
}
double pythia_card::BRTautomu(double mN1,double lambOverMSq)
{
	return Tautomu(mN1,lambOverMSq) / totalTau(mN1,lambOverMSq);
}




																	//(double mN, double mMeson, double fMeson,double mixingsq)

double pythia_card::Ntopi0(double mN1,double lambOverMSq)
{
	return widthCalculator::NtoP0(mN1,  widthCalculator::mpi0, widthCalculator::fpip, lambOverMSq) ;
}

double pythia_card::Ntoeta(double mN1,double lambOverMSq)
{
	return widthCalculator::NtoP0(mN1,  widthCalculator::meta, widthCalculator::feta, lambOverMSq);
}

double pythia_card::Ntoetap(double mN1,double lambOverMSq)
{
	return widthCalculator::NtoP0(mN1,  widthCalculator::metap, widthCalculator::fetap, lambOverMSq);
}



															// NtoV0(double mN, double mMeson, double fMeson, double gMeson, double mixingsq)


double pythia_card::Ntorho0(double mN1,double lambOverMSq)
{
	return widthCalculator::NtoV0(mN1,  widthCalculator::mrho0, widthCalculator::frho0, widthCalculator::grho0 ,lambOverMSq);
}
double pythia_card::Ntoomega(double mN1,double lambOverMSq)
{
	return widthCalculator::NtoV0(mN1,  widthCalculator::momega, widthCalculator::fomega, widthCalculator::gomega ,lambOverMSq);
}
double pythia_card::Ntophi(double mN1,double lambOverMSq)
{
	return widthCalculator::NtoV0(mN1,  widthCalculator::mphi, widthCalculator::fphi, widthCalculator::gphi ,lambOverMSq);
}


// N to e- e+ or mu- mu+
double pythia_card::Ntoll(double x, double SinThetaSq)
{
	if (x <= 1.7 and x >= 0.2 ) return (4.48116e-18 -1.04679e-16* x +8.88488e-16*pow(x,2)     - 3.42138e-15*pow(x,3) + 6.14543e-15*pow(x,4) - 
   3.74499e-15*pow(x,5) + 1.04271e-14*pow(x,6) - 
   7.98159e-15*pow(x,7) + 4.18499e-15*pow(x,8) - 
   1.43098e-15*pow(x,9) + 2.87294e-16*pow(x,10)- 
   2.56774e-17*pow(x,11)) * SinThetaSq;
	else return 0;
}
//N to 3 nu
double pythia_card::Ntonunu(double x, double SinThetaSq)
{
	if (x <= 1.7 and x >= 0 ) return   (pow(widthCalculator::GF,2)*4.* pow(x,5) * (1./(768. * pow(widthCalculator::pi,3))) ) * SinThetaSq;
	else return 0;
}





double pythia_card::totalN(double mN1,double lambOverMSq)
{
     return  Ntopi0(mN1,lambOverMSq) + Ntoeta(mN1,lambOverMSq) + Ntoetap(mN1,lambOverMSq) +Ntorho0(mN1,lambOverMSq) + Ntoomega(mN1,lambOverMSq)+ Ntophi(mN1,lambOverMSq) + Ntoll(mN1, lambOverMSq)+ Ntonunu(mN1, lambOverMSq);
}

double pythia_card::ctau(double mN1, double lambOverMSq)//ctau in meter
{
	double chbar=1.9732699e-16;//GeV m
	return chbar/ totalN(mN1,lambOverMSq);
}













double pythia_card::BRNtopi0(double mN1,double lambOverMSq)
{
	return Ntopi0(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}

double pythia_card::BRNtorho0(double mN1,double lambOverMSq)
{
	return Ntorho0(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}

double pythia_card::BRNtoeta(double mN1,double lambOverMSq)
{
	return Ntoeta(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}

double pythia_card::BRNtoetap(double mN1,double lambOverMSq)
{
	return Ntoetap(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}

double pythia_card::BRNtoomega(double mN1,double lambOverMSq)
{
	return Ntoomega(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}

double pythia_card::BRNtophi(double mN1,double lambOverMSq)
{
	return Ntophi(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}
double pythia_card::BRNtoll(double mN1,double lambOverMSq)
{
	return Ntoll(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}
double pythia_card::BRNtonunu(double mN1,double lambOverMSq)
{
	return Ntonunu(mN1,lambOverMSq)/totalN(mN1,lambOverMSq);
}





double pythia_card::BRn1ToVisibles (double mN1,double lambOverMSq)
{


	constexpr double BRpi0ToVisibles=0.;

	constexpr double BRrho0Topicpic=1.;

	constexpr double BRomegaTopicpicpi0 = 0.893;
	constexpr double BRomegaTopicpic    = 0.0153;


	constexpr double BRetaTopicpicpi0=0.2292;
	constexpr double BRetaTopicpicgamma=0.0422;
	constexpr double BRetaTogammagamma=0.3941;
	constexpr double BRetaTopi0pi0pi0=0.3268;

	constexpr double BRetapTopicpiceta=0.426;
	constexpr double BRetapTorhogamma = 0.289;
	constexpr double BRetapTopi0pi0eta = 0.228;
	constexpr double BRetapToomegagamma = 0.0262;
	constexpr double BRPhitokckc = 0.489;
	constexpr double BRPhitopicpicpi0 = 0.1532;
	

	constexpr double effi_picpic = 0.12;
	constexpr double effi_picpicpi0 = 0.084;
	constexpr double effi_picpicgamma = 0.084;
	constexpr double effi_picpicgammagamma = 0.084;
	constexpr double effi_picpicpicpicpi0  = 0.071;
	constexpr double effi_picpicpi0pi0pi0 = 0.041;
	constexpr double effi_picpicpicpicgamma = 0.071;
	constexpr double effi_picpicpi0pi0gamma = 0.041;
	constexpr double effi_picpicpi0gamma = 0.059;


	return (BRNtopi0(mN1,lambOverMSq)*BRpi0ToVisibles 
	     +  BRNtorho0(mN1,lambOverMSq)*BRrho0Topicpic*effi_picpic
  	     + BRNtoeta(mN1,lambOverMSq)*(BRetaTopicpicpi0*effi_picpicpi0  +  BRetaTopicpicgamma*effi_picpicgamma	    )
	     + BRNtoetap(mN1,lambOverMSq)*
	       (		   	       BRetapTopicpiceta*(BRetaTopicpicpi0*effi_picpicpicpicpi0  +  BRetaTopicpicgamma*effi_picpicpicpicgamma + BRetaTogammagamma*effi_picpicgammagamma  + BRetaTopi0pi0pi0*effi_picpicpi0pi0pi0  )  
					     + BRetapTorhogamma*(BRrho0Topicpic*effi_picpicgamma)
				             + BRetapTopi0pi0eta*(BRetaTopicpicpi0*effi_picpicpi0pi0pi0 + BRetaTopicpicgamma*effi_picpicpi0pi0gamma)
					     + BRetapToomegagamma*(BRomegaTopicpicpi0*effi_picpicpi0gamma  +  BRomegaTopicpic*effi_picpicgamma )
	       )
	     +  BRNtoomega(mN1,lambOverMSq)*(BRomegaTopicpicpi0*effi_picpicpi0  + BRomegaTopicpic*effi_picpic  ) +BRNtophi(mN1,lambOverMSq)*(BRPhitokckc*effi_picpic +BRPhitopicpicpi0 * effi_picpicpi0 )  +BRNtoll(mN1,lambOverMSq)*0.12 )*0.75;
}    // 0.75 is the background free additional efficiency by 1908.09719.  0.12 for BRNtoll from Abi.



double pythia_card::decayProbabilityBelle2Part1(Pythia8::Particle XXX) {//Part1 is a symmetric ATLAS-like detector
    constexpr double R_I = 0.1;//in meter
    constexpr double R_O = 0.8; 
    constexpr double L_D = 0.4;

    // Identify the kinematic properties of the neutralino
    //Pythia always calculates in GeV
    double gamma = XXX.e()/(mN1);
    double beta_z = fabs(XXX.pz()/XXX.e());
    double beta = sqrt(1. - pow(mN1/XXX.e(), 2));
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
    // The probability that the neutralino would decay in the detector is then as follows (Decay law:            
     // return exp(-L1/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (1. - exp(-L2/(beta_z*gamma*ctau(mN1,lambOverMSq)))); //  not work anymore!!!
    // the efficiency is a function (-10/7 *x + 8/7), x=0.1, eff=1, x=0.8, eff=0
    //so for every dL, the modified decay probability is e^(-L/(beta_z*gamma*ctau(mN1, lambOverMSq))) * 1/((beta_z*gamma*ctau(mN1, lambOverMSq))) * (-10/7 *L*tan(theta) + 8/7) *dL,from L1 to L1+L2
    return 2./7. * (exp(-L1/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (4.-5.* fabs(tan(theta)) * (beta_z*gamma*ctau(mN1, lambOverMSq)+ L1 ))    + exp(-(L1+L2)/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (-4.+5.* fabs(tan(theta)) * (beta_z*gamma*ctau(mN1, lambOverMSq)+ L2+L1 )) )  ;
}




// to get the total decay events inside the chamber, we can not use the modified decay probability
double pythia_card::decayProbabilityBelle2Part3(Pythia8::Particle XXX) {//Part1 is a symmetric ATLAS-like detector
    constexpr double R_I = 0.1;//in meter
    constexpr double R_O = 0.8; 
    constexpr double L_D = 0.4;

    // Identify the kinematic properties of the neutralino
    //Pythia always calculates in GeV
    double gamma = XXX.e()/(mN1);
    double beta_z = fabs(XXX.pz()/XXX.e());
    double beta = sqrt(1. - pow(mN1/XXX.e(), 2));
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
    // The probability that the neutralino would decay in the detector is then as follows (Decay law:            
      return exp(-L1/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (1. - exp(-L2/(beta_z*gamma*ctau(mN1,lambOverMSq)))); 
  
}












double pythia_card::decayProbabilityBelle2Part2(Pythia8::Particle XXX) {//Part2 is AL3X-like
    constexpr double L_D = 0.4;//horizontal distance from the IP to the detector
    constexpr double L_L = 0.1;//vertical distance from the IP to the detector  or  equivalently the inner radius
    constexpr double L_d = 0.8;// length of the detector: 120 - 40 = 80 cm
    constexpr double L_H = 0.7;// transverse length = 80 - 10 = 70 cm

    // Identify the kinematic properties of the neutralino
    //Pythia always calculates in GeV
    double gamma = XXX.e()/(mN1);
    double beta_z = fabs(XXX.pz()/XXX.e());
    double beta = sqrt(1. - pow(mN1/XXX.e(), 2));
    double theta = XXX.p().theta();            
    double phi = XXX.p().phi();     
 

    if (tan(theta) == 0.0) {
        std::cout << "The impossible happened!" << std::endl;
        return 0;
    }   
    double L1 = std::min(std::max(L_D,L_L/tan(theta)), L_D + L_d);
    double L2 = std::min(std::max(L_D,(L_L + L_H)/tan(theta)), L_D + L_d) - L1;
    if (L_L/tan(theta) >= L_D+L_d) // theta too small, neutralino flying too lowly
        return 0;
    if ((L_L + L_H)/tan(theta) <= L_D) //theta too large, neutralino flying too highly
        return 0; 
    // The probability that the neutralino would decay in the detector is then as follows (Decay law:            
    //return  exp(-L1/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (1. - exp(-L2/(beta_z*gamma*ctau(mN1, lambOverMSq)))); 
   
      return 2./7. * (exp(-L1/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (4.-5.* fabs(tan(theta)) * (beta_z*gamma*ctau(mN1, lambOverMSq)+ L1 ))    + exp(-(L1+L2)/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (-4.+5.* fabs(tan(theta)) * (beta_z*gamma*ctau(mN1, lambOverMSq)+ L2+L1 )) )  ;
   
}




// to get the decay events inside the chamber 
double pythia_card::decayProbabilityBelle2Part4(Pythia8::Particle XXX) {//Part2 is AL3X-like
    constexpr double L_D = 0.4;//horizontal distance from the IP to the detector
    constexpr double L_L = 0.1;//vertical distance from the IP to the detector  or  equivalently the inner radius
    constexpr double L_d = 0.8;// length of the detector: 120 - 40 = 80 cm
    constexpr double L_H = 0.7;// transverse length = 80 - 10 = 70 cm

    // Identify the kinematic properties of the neutralino
    //Pythia always calculates in GeV
    double gamma = XXX.e()/(mN1);
    double beta_z = fabs(XXX.pz()/XXX.e());
    double beta = sqrt(1. - pow(mN1/XXX.e(), 2));
    double theta = XXX.p().theta();            
    double phi = XXX.p().phi();     
 

    if (tan(theta) == 0.0) {
        std::cout << "The impossible happened!" << std::endl;
        return 0;
    }   
    double L1 = std::min(std::max(L_D,L_L/tan(theta)), L_D + L_d);
    double L2 = std::min(std::max(L_D,(L_L + L_H)/tan(theta)), L_D + L_d) - L1;
    if (L_L/tan(theta) >= L_D+L_d) // theta too small, neutralino flying too lowly
        return 0;
    if ((L_L + L_H)/tan(theta) <= L_D) //theta too large, neutralino flying too highly
        return 0; 
    // The probability that the neutralino would decay in the detector is then as follows (Decay law:            
    return  exp(-L1/(beta_z*gamma*ctau(mN1, lambOverMSq))) * (1. - exp(-L2/(beta_z*gamma*ctau(mN1, lambOverMSq)))); 
   

   
}



















