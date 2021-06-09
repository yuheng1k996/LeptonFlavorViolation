//////////////////////////////////////////////////
// Written by Giovanna Cottin (gfcottin@gmail.com)
//////////////////////////////////////////////////
// Inspired by ATLAS displaced vertex analysis in SUSY-2016-08
// Optimized for DV+elec
// Cuts: 
// Event Selection :
// One electron with pT>120 GeV within eta< 2.47
// Different triggers depending on analysis strategy
// DV Selection - Event pass if at least one displaced vertex with :
//   *   pT of tracks coming from DV  > 1 GeV
//   *   |d0| of tracks coming from DV  > 2 mm
//   This last two cuts define the signal region :
//   *   Number of tracks coming out of vertex > 5 
//   *   Invariant mass of vertex > 10 GeV 

#include "Pythia8/Pythia.h"
// #include "ToyDetector_leptons-ATLAS.h"

// ATLAS parametrized_truthEff at vertex level
#include "parametrized_truthEff.h"

#include <vector>
#include <string>
#include <iterator>
#include <algorithm>


// Translates a char* to a string
std::string charToString(char* x){
    std::ostringstream ss;
    ss << x;
    std::string s(ss.str());
    return s;
     
}


double neutrino2ejjgamma(double x){//partial decay width of N to ejj for unit mixing sqred
    double decaywidth{};
	if (1<x && x <= 6.275 )
	{
	decaywidth = 3.635039046420314e-12 - 1.0359367199120318e-11*x + 1.0567562718989518e-11*pow(x,2) - 
      4.226675885204053e-12*pow(x,3) + 2.361726847941139e-13*pow(x,4) + 3.4079490878956233e-13*pow(x,5) - 
      1.2704204641880472e-14*pow(x,6) + 1.2858274445647907e-15*pow(x,7) - 4.908486746146544e-17*pow(x,8);
	}
	else if (6.275 < x && x <= 82.06)
	{
	decaywidth = 1.7543515798612051e-9 - 5.449403377474963e-10*x + 
      6.885125520443618e-11*pow(x,2) - 6.551474712994119e-12*pow(x,3) + 2.777691533841118e-13*pow(x,4) + 
      2.9065709124055413e-13*pow(x,5) - 5.6079722222346335e-18*pow(x,6);
	}
	else if (82.06 < x )
	{
	decaywidth = (4.195908486e-7*(5.40215e11 - 1.25356e8*pow(x,2) + pow(x,6)))/pow(x,3);
	}
	else { 
        decaywidth = 0;
        }
	return decaywidth/2.;      //  /2. for Dirac
}

double neutrinolifetimefunction(double x) { //total width of nu4 in GeV for unit mixing sqred
    double decaywidth{};
	if (300<x && x <= 7000 )
	{
	decaywidth = -0.168099316322736 - 0.019176244649974233*x - 5.76312310809768e-8*pow(x,2) + 
      1.8673498869383768e-6*pow(x,3) - 5.367634497405086e-16*pow(x,4);
	}
	else if (120 < x && x <= 300)
	{
	decaywidth = 4.155725651877626 - 0.09043739071504049*x + 0.00044000314769598647*pow(x,2) + 
      6.560126640853136e-7*pow(x,3) + 1.2512719181228776e-9*pow(x,4);
	}
	else if (100 < x && x <= 120)
	{
	decaywidth = 18.77630370817525 - 0.6006856983548512*x + 0.006886162038434954*pow(x,2) - 
      0.000034378503914951415*pow(x,3) + 7.034971985716572e-8*pow(x,4);
	}
	else if (90 < x && x<=100)
	{
	decaywidth = -108.91502687400279 + 4.764103006740356*x - 0.07756810027258784*pow(x,2) + 
      0.0005560589056519785*pow(x,3) - 1.4764738630394293e-6*pow(x,4);
	}
	else if (82.31 < x && x <= 90)
	{
	decaywidth = 14.142625517703555 - 0.5479005663683709*x + 0.007913900147133817*pow(x,2) - 
      0.00005143345047679638*pow(x,3) + 1.3121272321849977e-7*pow(x,4);
	}
	else if (6.275 < x && x <= 82.31)
	{
	decaywidth = 2.2962852191318524e-9 - 7.069258842347951e-10*x + 9.348240881033966e-11*pow(x,2) - 
      9.217602129578459e-12*pow(x,3) + 3.6920457406342167e-13*pow(x,4) + 
      5.603338531128519e-13*pow(x,5) - 7.433911047747564e-18*pow(x,6);
	}
    	else if (0.8 < x && x <= 6.275)
	{
        decaywidth = -9.448531563556749e-11 + 5.022495628129828e-10*x - 1.1337774485856246e-9*pow(x,2) +       1.4307366794941334e-9*pow(x,3) - 1.1221845408247908e-9*pow(x,4) + 
      5.776653156928192e-10*pow(x,5) - 2.0004481252306443e-10*pow(x,6) + 
      4.704330325891135e-11*pow(x,7) - 7.40035741821473e-12*pow(x,8) + 
      7.456707016474243e-13*pow(x,9) - 4.3512612733257766e-14*pow(x,10) + 
      1.1183042344083286e-15*pow(x,11);
	}   
	else if (0.28 < x && x <= 0.8)
	{
        decaywidth =  -2.948649033457808e-11 + 7.004391570928488e-10*x - 7.46447899912376e-9*pow(x,2) + 
      4.7108479465027736e-8*pow(x,3) - 1.9564480245446248e-7*pow(x,4) + 
      5.615192973781638e-7*pow(x,5) - 1.1367146910727858e-6*pow(x,6) + 
      1.6234614593312558e-6*pow(x,7) - 1.6035895904413e-6*pow(x,8) + 
      1.0436679858269422e-6*pow(x,9) - 4.029484263538355e-7*pow(x,10) + 
      6.994315681784629e-8*pow(x,11);
	}
	else if ( 0.14 < x && x  <= 0.28)
	{
        decaywidth = 4.469971272478775e-16 - 6.7204618339090756e-15*x + 2.534573198113779e-14*pow(x,2) + 
      1.5711850175220388e-14*pow(x,3) - 2.2463913028976076e-13*pow(x,4) + 
      1.1935223602729758e-12*pow(x,5) - 8.085327618133332e-13*pow(x,6);
	}
	else if (0.05  < x && x  <= 0.14)
	{
        decaywidth =   4.6220833653384e-17 - 3.4865630734252773e-15*x + 
      1.0703197866730849e-13*pow(x,2) - 1.7126556021920763e-12*pow(x,3) + 
      1.507657556025811e-11*pow(x,4) - 6.897918816618373e-11*pow(x,5) + 
      1.3008401850849766e-10*pow(x,6);
	} else { 
        decaywidth = 0;
    }
	return decaywidth/2.;      //  /2. for Dirac
}  


double BrN2ejj(double x)
{
	return neutrino2ejjgamma(x) / neutrinolifetimefunction(x);
}

using namespace Pythia8;

//Calculate Delta R distance
// double CalcDR(double eta1, double eta2, double phi1, double phi2){
//     double dphi = abs(phi1-phi2);
//     if( dphi > M_PI) dphi = 2*M_PI- dphi;
//     double deta = eta1 - eta2;
//     return sqrt(dphi*dphi + deta*deta);
// }

// //Calculate Delta R cosmic (this is neglegible on signal efficiency)
// double CalcDRcosmic(double eta1, double eta2, double phi1, double phi2){
//     double dphi = abs(phi1-phi2);
//     if( dphi > M_PI) dphi = 2*M_PI- dphi;
//     double deta = eta1 - eta2;
//     return sqrt((M_PI-dphi)*(M_PI-dphi)-(eta1+eta2)*(eta1+eta2))
// }

int main(int argc, char *argv[]) {

if (argc != 8) {
        std::cout << "./main mLLP MixSq NMC LUMI Scenario XS LHE-path" << std::endl;
        std::cout << "   - mLLP: mass of LLP \n";
        std::cout << "   - MixSq:  mixing square \n";
        std::cout << "   - NMC: number of MC events for the Pythia simulation \n";       
        std::cout << "   - LUMI: luminosity in unit of fb^-1 \n";     
        std::cout << "   - Scenario: d11, q11, u11, diag_dqu, all_dqu, dqu11 \n";
        std::cout << "   - XS: XS of (pp > n1 n1~) in fb for this scenario and mass \n"; 
        std::cout << "   - LHE-path: path to the lhe file \n"; 
        exit(1);
    }

   double mLLP = atof(argv[1]);
   double MixSq = atof(argv[2]);
   int    NMC = atof(argv[3]);
   double lumi = atof(argv[4]);   
   string scenario= charToString(argv[5]);
   double sigma = atof(argv[6]);
   string lhepath= charToString(argv[7]); //  "/../../unweighted_events.lhe" as a string


   cout << "mLLP [GeV]: " << mLLP << '\n';
   cout << "MixSq : "<< MixSq <<'\n';
   cout << "NMC: " << NMC << '\n';
   cout << "Int Lumi [fb^{-1}]: " << lumi << '\n';

   Pythia pythia;
   Event& event = pythia.event;
   // Initialize Les Houches Event File run. List initialization information.
   pythia.readString("Beams:frameType = 4");
   /////////////////////////////////////////////////////////////////////////////////////

   pythia.readString("Beams:LHEF = "+lhepath);

   pythia.readString("Random:setSeed = on");
   pythia.readString("Random:seed = 0");//pick new random number seed for each run, based on clock     

   pythia.readString("PartonLevel:MPI = 0");//can switch off. speeds up a lot

   pythia.readString("Next:numberCount = 50000");

   cout<<lhepath<<endl;


   double ctauLLP{};
   int LLPPID{9990012};

   double chbar{1.9732699e-16};//GeV m

   int nEvent   = NMC;
  //int nEvent = 10000;//FROM LHE File
   int nAbort = 50;   

  ////// Initialize.
  pythia.init();

  //convert MixSq to string with scientific notation
  // Create an output string stream
  std::ostringstream streamObjMixSq;
  //Add double to stream
  streamObjMixSq << MixSq;
  // Get string from output string stream
  std::string MixSqString = streamObjMixSq.str();
  
  //save files
  //ofstream decays;
  //decays.open("decays_info/" + scenario + "/decays_GRID-"+std::to_string(mLLP)+"_"+MixSqString+".dat"); 
  

  // Enableing the following gives the same result, just checking
  // Disable SLHA for Neutrino. Neutrino is already decayed inside MadSpin
  // pythia.readString("SLHA:allowUserOverride=on");
  // pythia.readString("9910012:onMode = off");//  switch off all decay channels in pythia as N already decayed 

  // Consistency lifetime check
  cout<<"###################"<<endl;
  cout<<"pythia.particleData.tau0(" << LLPPID << ") = "<<pythia.particleData.tau0(LLPPID)<<endl;
  cout<<"pythia.particleData.m0(" << LLPPID << ")   = "<<pythia.particleData.m0(LLPPID)<<endl;


  int iAbort = 0;
  int nPromptElec=0;
  int nDVReco=0;
  int nMaterial=0;
  int nFidutial=0;
  int nTrk=0;
  int nMass=0;
  int nEvtEff=0;
  int nDVEff=0;

  // ToyDetector detector;

  // // Begin event loop; generate until none left in input file.
  // for (int iEvent = 0; ; ++iEvent) {
  // // for (int iEvent = 0; iEvent < 10; ++iEvent) {
  //   cout << "######## Event "<<iEvent<<" #########"<<endl;

  //   // Generate events, and check whether generation failed.
  //   if (!pythia.next()) {
  //     // If failure because reached end of file then exit event loop.
  //     if (pythia.info.atEndOfFile()) break;
  //     // First few failures write off as "acceptable" errors, then quit.
  //     if (++iAbort < nAbort) continue;
  //     break;
  //   }

    //Start event loop, need consistent nEvent 
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        // Generate events. Quit if many failure.
        if (!pythia.next()) {
          // event.list();
          cout<<" inside pythia!next"<<endl;
          if (++iAbort < nAbort) continue;
            cout << " Event generation aborted prematurely, owing to error!\n";
          break;
        }
        
    bool passEvtPromptElec  = false; // electron trigger
    bool passElecEff        = false; // electron efficiency
    /////////////////
    //Event Selection
    /////////////////
    std::vector<double> elecIndex;
    ///////////////////////////////////////
    // Electron trigger
    ///////////////////////////////////////
    for (int i= 0; i < event.size(); i++){
      if (abs(event[i].id()) == 11){
        // Apply electron efficiency
        float iRndNumber = (std::rand()/(float)RAND_MAX); 
        // cout<<"iRndNumber = "<<iRndNumber<<endl;
        passElecEff = 0.7 >iRndNumber;
        if((event[i].pT()>120.) && abs(event[i].eta())<2.47) {
          passEvtPromptElec=true;
          elecIndex.push_back(i);
        }
      }
    }
    if(passEvtPromptElec && passElecEff) nPromptElec++;
      else continue; 
      /////////////////////////////////////////////
      //Displaced Vertex truth identification
      //Find all particles coming from a Neutrino
      /////////////////////////////////////////////  
      std::vector<int> motherIndices;   
      //Get neutrinos from event
      for (int i= 0; i < event.size(); i++){
        if (abs(event[i].id()) == LLPPID and event[event[i].daughter1()].id() != LLPPID   ){
	  //mLLP = event[i].m();
	  ctauLLP = event[i].tau0();
          double trux = event[i].xDec() - event[i].xProd();
          double truy = event[i].yDec() - event[i].yProd();
          double truz = event[i].zDec() - event[i].zProd();    
          double vProd = sqrt(pow2(event[i].xProd())+pow2(event[i].yProd())+pow2(event[i].zProd()));     
           //If neutrino is displaced 
          if ((abs(trux) > 0.0001) || (abs(truy) > 0.0001)) {
            double decayLength = sqrt(pow2(trux) + pow2(truy) + pow2(truz));
            double p = sqrt(pow2(event[i].px())+pow2(event[i].py())+pow2(event[i].pz()));
            double mass = event[i].m();
            double betagamma =  p/mass;

            double dist = event[i].vDec().pAbs();//should be equal to decayLenght
            //In this case dcyLen and dist are the same. vDec() gives the decay vertex coordinates,  
            //calculated from the production vertex (at origin) and the proper lifetime 
             //cout<<"###################"<<endl;
             //cout<<"dist (from vDec) [mm] =  "<< event[i].vDec().pAbs()<<endl;
             //cout<<"decayLength [mm]         =  "<< decayLength<<endl;
             //cout<<"gamma                    =  "<< gamma<<endl;
            motherIndices.push_back(i);
            //decays<<pythia.particleData.tau0(LLPPID)<<" "<<betagamma<<" "<<decayLength<<endl;
          }
        }
      }
     // At least one DV in event 
     // Cuts
      bool vtxPassFidutial = false;
      bool vtxPassesNtrk   = false;
      bool vtxPassesMass   = false;
      bool passDVEff       = false;
      bool passEvtEff      = false; 
      bool passDVAll       = false;
      double Rdecay_largest = 0.0;
     // // Loop over DVs
       // cout<<"motherIndices.size() = "<<motherIndices.size()<<endl;
      for (int i = 0; i < motherIndices.size(); ++i) {
        int DVindex = motherIndices[i];
        double xDV=event[DVindex].xDec() - event[DVindex].xProd();  
        double yDV=event[DVindex].yDec() - event[DVindex].yProd();  
        double zDV=event[DVindex].zDec() - event[DVindex].zProd();  
        double rDV=sqrt(xDV*xDV+yDV*yDV);
        //cout<< "rDV = "<<rDV<<endl;
        Vec4 total4p;
        Vec4 total4pPion;  
        //Get all daughters
        vector<int> daughters= event[DVindex].daughterListRecursive();// Need this method to get all daughters
        std::vector<int> daughterIndices;
        for (int j = 0; j < daughters.size(); ++j){
          if (!event[daughters[j]].isFinal()) continue;
          if (std::find(daughterIndices.begin(),daughterIndices.end(),daughters[j]) != daughterIndices.end()) continue;
          daughterIndices.push_back(daughters[j]);
        }
        std::vector<int> chargedDaughters;
        // cout<<"daughterIndices.size()="<< daughterIndices.size()<<endl;         
        
        // find triggered electron on list of daughters from DV
        //exitst will return if the electron index is matched to track index 
        for(unsigned int trk =0;trk<daughterIndices.size();trk++){
          int trackIndex=daughterIndices[trk];
          std::vector<int> trackidxs;
          trackidxs.push_back(trackIndex);
            for(int idx=0;idx<elecIndex.size();idx++){
              int dvElecIndex=elecIndex[idx];
              bool exists = std::find(std::begin(trackidxs), std::end(trackidxs), dvElecIndex) != std::end(trackidxs);
              // if(exists){
              // cout<<"dvElecIndex =  "<<dvElecIndex<<endl;
              // cout<<"trackIndex = "<<trackIndex<<endl;
              // }
              if(!exists) break; // leave if no electron matched to dv
              // cout<<" exists = "<<exists<<endl;
          }// close electron index
          //track cuts inside DV  
          bool passTrackQuality  = true; // true unless probe otherwise
          double rTrk=sqrt(pow2((event[trackIndex].vProd()).px())+ pow2((event[trackIndex].vProd()).py()));
           // cout<<"rTrk = "<<rTrk<<endl;
          double zTrk=abs((event[trackIndex].vProd()).pz());   
          double phixy = event[trackIndex].vProd().phi();
          double deltaphi = phixy - event[trackIndex].phi();
          if (abs(deltaphi) > 3.1415) deltaphi = 2 * 3.1415 - abs(deltaphi);
          //d0 in the abscence of magnetic field
          //Make sure tracks are displaced
          double d0 = rTrk*sin(deltaphi);
          // cout<<"event[trackIndex].id()="<<event[trackIndex].id()<<endl;
          //select good daughter tracks
          // cout<<"event[trackIndex].isCharged()="<<event[trackIndex].isCharged()<<endl;
          // cout<<"d0 = "<<d0<<endl;
          if(!event[trackIndex].isCharged()) continue;
          if(abs(d0)<2.0 || event[trackIndex].pT()<1.0) continue;
          chargedDaughters.push_back(trackIndex);
          // DV invariant mass calculation with pionMass hypothesis  
          double particlepx = event[trackIndex].px();
          double particlepy = event[trackIndex].py();
          double particlepz = event[trackIndex].pz();
          double particleE = event[trackIndex].e();
          double particleEPion = sqrt(pow2(0.1395700)+pow2(particlepx)+pow2(particlepy)+pow2(particlepz));
          total4p+=Vec4(particlepx,particlepy,particlepz,particleE);
          total4pPion+=Vec4(particlepx,particlepy,particlepz,particleEPion);
        }// close track loop
        double totalmPion=total4pPion.mCalc();
        // Apply DV cuts. At leats ONE DV satisfying all criteria
        //cout<<"chargedDaughters.size() ="<<chargedDaughters.size()<<endl;
        //cout<<"totalmPion              = "<<totalmPion<<endl;
        vtxPassFidutial = (4.<rDV) && (rDV<300) && (abs(zDV)<300);
        vtxPassesNtrk   = chargedDaughters.size()>=4;
        vtxPassesMass   = totalmPion>=5.;
        // Apply DV efficiency
        float iRndNumber1 = (std::rand()/(float)RAND_MAX); 
        passDVEff         = vertexEff_Regions(rDV, totalmPion, chargedDaughters.size())>iRndNumber1;
        if(vtxPassFidutial && vtxPassesNtrk && vtxPassesMass && passDVEff) {
          //Find max of the one/two DVs passing all selection
          if(rDV>Rdecay_largest) {
            Rdecay_largest=rDV;
          }
          passDVAll=true;
          break;// leave vertex for for, as I have found at least one satisfying all criteria        
        }
      }//close vertex loop
      if(passDVAll){
        // cout<<"AT LEAST ONE VERTEX PASSES ALL CONDITIONS"<<endl;
      }
      // Apply event level efficiency (with MET trigger) using the maximum rDV present in the event.
      // float iRndNumber2 = (std::rand()/(float)RAND_MAX); 
      // if(eventEff_MET(Rdecay_largest,truthMET.pT())>iRndNumber2) passEvtEff=true;
      //CutFlow count
      if(vtxPassFidutial) nFidutial++;
        else continue;
      if(vtxPassesNtrk) nTrk++;
        else continue;
      if(vtxPassesMass) nMass++;
        else continue;        
      if(passDVEff) nDVEff++;
        else continue;    
     
    }//End of event loop.
   

    //sigma = pythia.info.sigmaGen()*1e12;//convert mb to fb    not using this, because it includes br(n1->ejj) as we use madspin, but the MG5 computation is not good for mN >  the EW scale,  use pow(BrN2ejj(mLLP),2) below for the computation of signal number

    pythia.stat();

    cout << "XS [fb]: " << sigma <<'\n'; //in fb
    cout << "Gamma [GeV]: " << neutrinolifetimefunction(mLLP)*MixSq << '\n';
    cout << "ctau [m]: "<< chbar/(neutrinolifetimefunction(mLLP)*MixSq) << '\n';   //when ctau is too large, pythia8 tau0() function returns 0
    cout << "Br(N->ejj): " << BrN2ejj(mLLP) << '\n';

    cout << "No. of Events,  step-wise efficiency in %, cumulative efficiency in % \n";
    cout << "All events "        << nEvent                << " " << nEvent*100./double(nEvent)            << " " << nEvent*100./double(nEvent)    << '\n';
    cout << "Electron trigger  " << nPromptElec           << " " << nPromptElec*100./double(nEvent)           << " " << nPromptElec*100./double(nEvent) << '\n';
    cout << "DV Fidutial "       << nFidutial             << " " << nFidutial*100./double(nPromptElec)    << " " << nFidutial*100./double(nEvent) << '\n';
    cout << "DV Ntrk "           << nTrk                  << " " << nTrk*100./double(nFidutial)     << " " << nTrk*100./double(nEvent)      << '\n';
    cout << "DV Mass "           << nMass                 << " " << nMass*100./double(nTrk)               << " " << nMass*100./double(nEvent)     << '\n';
    cout << "DV Eff "            << nDVEff                << " " << nDVEff*100./double(nMass)             << " " << nDVEff*100./double(nEvent)    << '\n';

    cout << "Final cut efficiency: " << nDVEff/double(nEvent) << '\n';
    cout << "Final signal number: " << sigma * lumi * ( 1. - (1. - BrN2ejj(mLLP)) * (1. - BrN2ejj(mLLP))  ) * (nDVEff/double(nEvent))  << '\n';



  //decays.close();
  return 0;
}

