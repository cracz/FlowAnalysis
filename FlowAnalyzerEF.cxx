//////
// Do not use vector::erase() to change any vectors. 
// This will invalidate any for() loops iterating over the vectors
// and make things much more complicated. For bad events after 
// creation of the "Event" vector, use the badEvent flag.
//////

// C++ headers
#include <iostream>
#include <vector>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TSystem.h"
#include "TKey.h"
#include "TMath.h"

// PicoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"


// EPD Util headers
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StBbcGeom.h"



//=========================================================
//          SOME CONTROLS
//=========================================================
const Int_t ORDER_N = 1;                // Order of anisotropic flow
TString ORDER_N_STR;                   // ORDER_N but as a TString for titles/labels

const Double_t PSI_BOUNDS = TMath::Pi()/ORDER_N + 1;
const Double_t Q_BOUNDS = 100;

const Int_t EPD_FORMAT       = 2;       // format=0/1/2 for StEpdHit/StMuEpdHit/StPicoEpdHit
const Int_t EPD_MAX_WEIGHT   = 2;      // max nMIP weight; recommended value, but variable
const Double_t EPD_THRESHOLD = 0.3;    // recommended value, but variable

const Double_t MIN_TPC_ETA_CUT = -2.0;
const Double_t AGAP_TPC_ETA_CUT = -1.1;//-1.6;
const Double_t GAPB_TPC_ETA_CUT = -1.0;//-0.5;
const Double_t MAX_TPC_ETA_CUT = 0.0;
/*
const Double_t MIN_ETA_CUT = -5.6;//-5.16;    // Cuts for the EPD subevents (A, B, C, D)
const Double_t AB_ETA_CUT  = -4.05;//-3.82;
const Double_t BC_ETA_CUT  = -3.3;//-3.28;
const Double_t CD_ETA_CUT  = -2.9;//-2.87;
const Double_t MAX_ETA_CUT = -2.30;
*/

const Double_t MIN_ETA_CUT = -5.6;
const Double_t EF_ETA_CUT  = -3.5;
const Double_t MAX_ETA_CUT = -2.4;

const Double_t R_VTX_CUT = 2.0;         // 2D r value, good vertices are within this value
const Double_t Z_VTX_CUT_LOW  = 199.5;
const Double_t Z_VTX_CUT_HIGH = 201.5;

const Int_t MIN_TRACKS = 5;             // Min number of tracks/hits in each sub-event

const Double_t PI_KA_MOM_CUT = 1.6;     // Reject those above these total momenta (GeV/c)
const Double_t PR_MOM_CUT    = 3.0;

const Int_t SHIFT_TERMS = 10;           // Number of terms to use when shifting event plane angles

const Int_t I_BAD_VALUE    = -99;
const Double_t D_BAD_VALUE = -99.0;

const Int_t CENT_BINS  = 16;             // Number of centrality bins (max 16)
const Int_t FIRST_CENT = 16 - CENT_BINS;            // Starting point for centrality dependent plots

const Double_t Y_MID = -1.05;       // Mid rapidity
/*
const Double_t Y_MID_PI_P = -0.82;       // Mid rapidity for pions (plus and minus), kaons, and protons
const Double_t Y_MID_PI_M = -0.82;
const Double_t Y_MID_KA_P = -0.72;
const Double_t Y_MID_KA_M = -0.77;
const Double_t Y_MID_PR   = -0.52;
*/
Int_t RUN_ITERATION = 0;
// 0 = No correction info yet; save raw (Xn,Yn) distributions
// 1 = Correction file found, but only <Xn> and <Yn> for re-centering.
//     Also save <sin> <cos> at this step for shifting in the next step.
// 2 = Correction file found, and <sin> <cos> values found so that shifting can be performed.
//=========================================================
//          
//=========================================================




// Custom type to hold important info for every good event
struct Event
{
  bool badEvent;      // Flag for marking events to ignore
  Int_t centID;
  Int_t primTracks;   // Number of primary tracks before track cuts (used for centrality)

  Int_t nTracksTpc;
  Double_t XnTpc;
  Double_t YnTpc;
  Double_t psiTpc; 

  Int_t nTracksTpcA;      // Number of GOOD tracks in the sub-event
  Double_t XnTpcA;
  Double_t YnTpcA;
  Double_t psiTpcA;       // Overall EP angle without removing autocorrelations
  std::vector<Double_t> phiValuesTpcA;   // Azimuthal angles for all TPC A particles in the event
  std::vector<Double_t> etaValuesTpcA;

  Int_t nTracksTpcB;
  Double_t XnTpcB;
  Double_t YnTpcB;
  Double_t psiTpcB;
  std::vector<Double_t> phiValuesTpcB;
  std::vector<Double_t> etaValuesTpcB;

  std::vector<Double_t> phiValuesPionP;
  std::vector<Double_t> phiValuesPionM;
  std::vector<Double_t> phiValuesKaonP;
  std::vector<Double_t> phiValuesKaonM;
  std::vector<Double_t> phiValuesProton;

  Int_t nhitsEpd;
  Double_t XnEpd;
  Double_t YnEpd;
  Double_t psiEpd;

  Int_t nhitsEpdE;
  Double_t XnEpdE;
  Double_t YnEpdE;
  Double_t psiEpdE;
  std::vector<Double_t> phiValuesEpdE;
  std::vector<Double_t> etaValuesEpdE;

  Int_t nhitsEpdF;
  Double_t XnEpdF;
  Double_t YnEpdF;
  Double_t psiEpdF;
  std::vector<Double_t> phiValuesEpdF;
  std::vector<Double_t> etaValuesEpdF;


  void reset()
  {
    badEvent  = false;  //Reset all values in the struct to reuse
    primTracks = 0;
    centID = I_BAD_VALUE;

    nTracksTpc = 0;
    XnTpc = 0;
    YnTpc = 0;
    psiTpc = D_BAD_VALUE;        //Just some number to use that is out of bounds

    nTracksTpcA = 0;
    XnTpcA = 0;
    YnTpcA = 0;
    psiTpcA = D_BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcA);
    std::vector<Double_t>().swap(etaValuesTpcA);

    nTracksTpcB = 0;
    XnTpcB = 0;
    YnTpcB = 0;
    psiTpcB = D_BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcB);
    std::vector<Double_t>().swap(etaValuesTpcB);

    std::vector<Double_t>().swap(phiValuesPionP);
    std::vector<Double_t>().swap(phiValuesPionM);
    std::vector<Double_t>().swap(phiValuesKaonP);
    std::vector<Double_t>().swap(phiValuesKaonM);
    std::vector<Double_t>().swap(phiValuesProton);

    nhitsEpd = 0;
    XnEpd = 0;
    YnEpd = 0;
    psiEpd = D_BAD_VALUE;

    nhitsEpdE = 0;
    XnEpdE = 0;
    YnEpdE = 0;
    psiEpdE = D_BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdE);
    std::vector<Double_t>().swap(etaValuesEpdE);

    nhitsEpdF = 0;
    XnEpdF = 0;
    YnEpdF = 0;
    psiEpdF = D_BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdF);
    std::vector<Double_t>().swap(etaValuesEpdF);
  }
};




////////
//   Moves event plane angles back into the correct period.
////////
Double_t angleShift(Double_t angle, Int_t order)
{
  if (angle < -TMath::Pi()/(Double_t)order) { angle += TMath::TwoPi()/(Double_t)order; }
  else if (angle >  TMath::Pi()/(Double_t)order) { angle -= TMath::TwoPi()/(Double_t)order; }
  return angle;
};

////////
//   Checks if an event has any flow vectors equal to zero. Updates the event's member variable "badEvent".
////////
void checkZeroQ(Event event)
{
  if (event.XnTpc == 0 && event.YnTpc == 0) { event.badEvent = true; }
  //else if (event.XnTpcA == 0 && event.YnTpcA == 0) { event.badEvent = true; }
  else if (event.XnTpcB == 0 && event.YnTpcB == 0) { event.badEvent = true; }
  else if (event.XnEpd == 0 && event.YnEpd == 0) { event.badEvent = true; }
  else if (event.XnEpdE == 0 && event.YnEpdE == 0) { event.badEvent = true; }
  else if (event.XnEpdF == 0 && event.YnEpdF == 0) { event.badEvent = true; }
};


////////
//   Using px, py, pz, and rest mass, return rapidity
////////
Double_t rapidity(Double_t px, Double_t py, Double_t pz, Double_t mass)
{
  Double_t rapidity, energy, momentum = 0;
  momentum = TMath::Sqrt(px*px + py*py + pz*pz);
  energy   = TMath::Sqrt(momentum*momentum + mass*mass);
  rapidity = TMath::ATanH(pz/energy);
  return rapidity;
};

////////
//   Using px, py, pz, and rest mass, return transverse mass
////////
Double_t transMass(Double_t px, Double_t py, Double_t mass) {return TMath::Sqrt(mass*mass + px*px + py*py);};


////////
//   Using px, py, pz, and rest mass, fill histograms raw of dN/dy,
// and dN/dmT (shifted left by m0), and a 2D histogram of mT-m0 vs y.
////////
void fillRawSpect(Double_t px, Double_t py, Double_t pz, Double_t mass, TH1D *dndy, TH1D *dndm, TH2D *MvsY)
{
  Double_t y  = rapidity(px, py, pz, mass);
  Double_t mT = transMass(px, py, mass);
  Double_t M  = mT - mass;
  dndy->Fill(y);
  dndm->Fill(M);
  MvsY->Fill(y,M, 1/(TMath::TwoPi() * mT));
};





void FlowAnalyzerEF(TString inFile, TString jobID)
{
  std::cout << "Initializing..." << std::endl;

  if (gSystem->AccessPathName(inFile)) { std::cout << "Error reading input file!" << std::endl; return;}


  ORDER_N_STR.Form("%d", ORDER_N);


  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("EpdHit",1);
  picoReader->SetStatus("BTofHit",1);
  picoReader->SetStatus("BTofPidTraits",1);
  if (!picoReader->chain()) { std::cout << "No chain found." << std::endl; return; }
  
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  Long64_t events2read  = picoReader->chain()->GetEntries();
    
  std::cout << "Number of events to read: " << events2read << std::endl;


  TClonesArray *epdHits = new TClonesArray("StPicoEpdHit");
  unsigned int found;

  TChain *picoChain = picoReader->chain();
  picoChain->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk                                                                                
  std::cout << "EpdHit Branch returned found = " << found << std::endl;
  picoChain->SetBranchAddress("EpdHit",&epdHits);


  /*
  StEpdEpFinder *epdEpFinder = new StEpdEpFinder(CENT_BINS, "StEpdEpFinderCorrectionHistograms_OUTPUT_"+jobID+".root", "StEpdEpFinderCorrectionHistograms_INPUT.root");
  epdEpFinder->SetEpdHitFormat(EPD_FORMAT);
  epdEpFinder->SetnMipThreshold(EPD_THRESHOLD);
  epdEpFinder->SetMaxTileWeight(EPD_MAX_WEIGHT);
  */

  // INPUT FILE FOR CORRECTION INFORMATION
  TString correctionInputName = "correctionInfo_INPUT.root";
  TFile *correctionInputFile = TFile::Open(correctionInputName, "READ");
  if (!correctionInputFile)
    {
      RUN_ITERATION = 0;
      std::cout << "No correction file found." << std::endl
		<< "Re-centering and shifting will not be performed." << std::endl;
    }
  else
    {
      TKey *key;
      TIter next(correctionInputFile->GetListOfKeys());
      TProfile *profile;

      while( (key = (TKey*)next()) )
	{
	  TClass *cl = gROOT->GetClass(key->GetClassName());
	  
	  if (cl->InheritsFrom("TProfile"))
	    {
	      profile = (TProfile*)key->ReadObj();
	      if (profile->GetEntries() == 0)
		{
		  std::cout << "TProfiles are empty!" << std::endl
			    << "Re-centering will be performed and TProfiles will be filled." << std::endl;
		  RUN_ITERATION = 1;
		  break;
		}
	      else if (profile->GetEntries() != 0)
		{
		  std::cout << "Non-empty TProfiles found!" << std::endl
			    << "Re-centering and event plane shifting will be performed." << std::endl;
		  RUN_ITERATION = 2;
		  break;
		}
	    }
	}
    }


  // OUTPUT FILE FOR CORRECTION INFORMATION
  TString correctionOutputName = "correctionInfo_OUTPUT_"+jobID+".root";
  TFile *correctionOutputFile;
  if (RUN_ITERATION == 0 || RUN_ITERATION == 1) { correctionOutputFile = new TFile(correctionOutputName, "RECREATE"); }

  // MAIN OUTPUT FILE
  TString outFile = jobID+".picoDst.result.root";
  TFile *outputFile = new TFile(outFile,"RECREATE");
  outputFile->cd();


  // HISTOGRAMS
  TH1D *h_event_check = new TH1D("h_event_check","Event number after each cut;;Events", 3, 0, 3);
  const char *eventSections[3] = {"No cuts", "Trigger cut", "Vertex cut"};  
  h_event_check->SetStats(0);

  TH1D *h_track_check = new TH1D("h_track_check","Track number after each cut;;Tracks", 4, 0, 4);
  const char *trackSections[4] = {"Event cuts only", "TOF #beta cut", "QA cuts", "PID cuts"};  
  h_track_check->SetStats(0);

  TH1D *h_nhits      = new TH1D("h_nhits", "nHits;Number of hits;Tracks", 50, 0, 50);
  TH1D *h_nhits_fit  = new TH1D("h_nhits_fit","nHitsFit;Number of hits;Tracks", 50, 0, 50);
  TH1D *h_nhits_dEdx = new TH1D("h_nhits_dEdx","nHitsdEdx;Number of hits;Tracks", 50, 0, 50);

  TH1D *h_primTracks = new TH1D("h_primTracks","Raw Number of Primary Tracks;Tracks;Events", 200, 0, 200);

  TH1D *h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", 100, 190, 210);

  TH1D *h_eta_s   = new TH1D("h_eta_s", "Particle #eta_{CM};#eta-#eta_{mid};Particles", 600, -6, 2);
  TH1D *h_eta_TPC_s = new TH1D("h_eta_TPC_s", "TPC tracks' #eta_{CM};#eta-#eta_{mid};Particles", 600, -2, 2);


  TH1D *h_tileWeights = new TH1D("h_tileWeights", "EPD Tile Weights;Hits;nMIP Weights", 5, -1, 4);
  TH1D *h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TH1D *h_pp_dndm = new TH1D("h_pp_dndm", "#pi^{+} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}", 60, 0, 3);//30, 0, 3);
  TH1D *h_pm_dndm = new TH1D("h_pm_dndm", "#pi^{-} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}", 60, 0, 3);//30, 0, 3);
  TH1D *h_kp_dndm = new TH1D("h_kp_dndm", "K^{+} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",   60, 0, 3);//30, 0, 3);
  TH1D *h_km_dndm = new TH1D("h_km_dndm", "K^{-} Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",   60, 0, 3);//30, 0, 3);
  TH1D *h_pr_dndm = new TH1D("h_pr_dndm", "Proton Raw m_{T} Spectrum;m_{T}-m_{0} (GeV);dN/dm_{T}",  60, 0, 3);//30, 0, 3);

  TH1D *h_pp_dndy = new TH1D("h_pp_dndy", "#pi^{+} Raw Rapidity Spectrum;y;dN/dy", 40, -2, 0);//50, -2, 0.5);
  TH1D *h_pm_dndy = new TH1D("h_pm_dndy", "#pi^{-} Raw Rapidity Spectrum;y;dN/dy", 40, -2, 0);//50, -2, 0.5);
  TH1D *h_kp_dndy = new TH1D("h_kp_dndy", "K^{+} Raw Rapidity Spectrum;y;dN/dy",   40, -2, 0);//50, -2, 0.5);
  TH1D *h_km_dndy = new TH1D("h_km_dndy", "K^{-} Raw Rapidity Spectrum;y;dN/dy",   40, -2, 0);//50, -2, 0.5);
  TH1D *h_pr_dndy = new TH1D("h_pr_dndy", "Proton Raw Rapidity Spectrum;y;dN/dy",  40, -2, 0);//50, -2, 0.5);

  TH1D *h_pp_pT = new TH1D("h_pp_pT", "#pi^{+} p_{T};p_{T} (GeV);", 100, 0, 5);
  TH1D *h_pm_pT = new TH1D("h_pm_pT", "#pi^{-} p_{T};p_{T} (GeV);", 100, 0, 5);
  TH1D *h_kp_pT = new TH1D("h_kp_pT", "K^{+} p_{T};p_{T} (GeV);",   100, 0, 5);
  TH1D *h_km_pT = new TH1D("h_km_pT", "K^{-} p_{T};p_{T} (GeV);",   100, 0, 5);
  TH1D *h_pr_pT = new TH1D("h_pr_pT", "Proton p_{T};p_{T} (GeV);",  100, 0, 5);

  TH1D *h_pp_mom = new TH1D("h_pp_mom", "#pi^{+} Total Momentum;|p| (GeV);", 100, 0, 5);
  TH1D *h_pm_mom = new TH1D("h_pm_mom", "#pi^{-} Total Momentum;|p| (GeV);", 100, 0, 5);
  TH1D *h_kp_mom = new TH1D("h_kp_mom", "K^{+} Total Momentum;|p| (GeV);",   100, 0, 5);
  TH1D *h_km_mom = new TH1D("h_km_mom", "K^{-} Total Momentum;|p| (GeV);",   100, 0, 5);
  TH1D *h_pr_mom = new TH1D("h_pr_mom", "Proton Total Momentum;|p| (GeV);",  100, 0, 5);

  TH2D *h2_pp_MvsY  = new TH2D("h2_pp_MvsY", "#pi^{+} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}", 16, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_pm_MvsY  = new TH2D("h2_pm_MvsY", "#pi^{-} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}", 16, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_kp_MvsY  = new TH2D("h2_kp_MvsY", "K^{+} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}",   16, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_km_MvsY  = new TH2D("h2_km_MvsY", "K^{-} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}",   16, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_pr_MvsY  = new TH2D("h2_pr_MvsY", "Proton m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}",  16, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);


  TH1D *h_psiTpc_RAW  = new TH1D("h_psiTpc_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcA_RAW = new TH1D("h_psiTpcA_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_RAW = new TH1D("h_psiTpcB_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpd_RAW  = new TH1D("h_psiEpd_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdE_RAW = new TH1D("h_psiEpdE_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD E);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdF_RAW = new TH1D("h_psiEpdF_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD F);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  TProfile *p_v2Plot = new TProfile("p_v2Plot", "v_{2} by Centrality;Centrality;<cos(2(#phi - #psi_{"+ORDER_N_STR+"}))>", CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);

  // Profiles for resolution terms
  TProfile *p_TpcAB = new TProfile("p_TpcAB","TPC A-B Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{TPC,B}_{"+ORDER_N_STR+"}))>",
				   CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcAEpdE = new TProfile("p_TpcAEpdE","TPC A EPD E Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,E}_{"+ORDER_N_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdF = new TProfile("p_TpcAEpdF","TPC A EPD F Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,F}_{"+ORDER_N_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcBEpdE = new TProfile("p_TpcBEpdE","TPC B EPD E Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,E}_{"+ORDER_N_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdF = new TProfile("p_TpcBEpdF","TPC B EPD F Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,F}_{"+ORDER_N_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_EpdEEpdF = new TProfile("p_EpdEEpdF","EPD E EPD F Correlations;Centrality;<cos(2(#psi^{EPD,E}_{"+ORDER_N_STR+"}-#psi^{EPD,F}_{"+ORDER_N_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);


  TH2D *h2_pp_vs_eta = new TH2D("h2_pp_vs_eta","Tile Weight for Supersectors vs #eta;#eta;Supersector", 400, -6, -2, 12, 0.5, 12.5);
  TH2D *h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  TH2D *h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);

  TH2D *h2_betap  = new TH2D("h2_betap","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_m2_p   = new TH2D("h2_m2_p", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  TH2D *h2_m2_vs_qpT  = new TH2D("h2_m2_vs_qpT", "m^{2} vs q*p_{T};q*p_{T} (GeV);m^{2} (GeV^{2})", 500, -3, 3, 500, -0.1, 1.2);
  TH2D *h2_dEdx_vs_qp = new TH2D("h2_dEdx_vs_qp", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_nSig_vs_qp_pi = new TH2D("h2_nSig_vs_qp_pi", "Pion n#sigma vs q|p|;q|p| (GeV); n#sigma_{#pi}", 500, -5, 5, 400, -8, 8);
  TH2D *h2_nSig_vs_qp_ka = new TH2D("h2_nSig_vs_qp_ka", "Kaon n#sigma vs q|p|;q|p| (GeV); n#sigma_{K}", 500, -5, 5, 400, -8, 8);
  TH2D *h2_nSig_vs_qp_pr = new TH2D("h2_nSig_vs_qp_pr", "Proton n#sigma vs q|p|;q|p| (GeV); n#sigma_{p}", 500, -5, 5, 400, -8, 8);
  TH2D *h2_pi_m2_vs_TPC_nsig = new TH2D("h2_pi_m2_vs_TPC_nsig", "m^{2} vs #pi TPC n#sigma;n#sigma_{#pi};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);
  TH2D *h2_ka_m2_vs_TPC_nsig = new TH2D("h2_ka_m2_vs_TPC_nsig", "m^{2} vs K TPC n#sigma;n#sigma_{K};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);
  TH2D *h2_pr_m2_vs_TPC_nsig = new TH2D("h2_pr_m2_vs_TPC_nsig", "m^{2} vs Proton TPC n#sigma;n#sigma_{pro};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);

  TH2D *h2_y_vs_eta = new TH2D("h2_y_vs_eta", "TPC All Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_pp = new TH2D("h2_y_vs_eta_pp", "TPC #pi^{+} y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_pp = new TH2D("h2_y_vs_eta_pt0p5to1_pp", "#pi^{+} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pp = new TH2D("h2_y_vs_eta_pt1to1p5_pp", "#pi^{+} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pp = new TH2D("h2_y_vs_eta_pt1p5to2_pp", "#pi^{+} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_pm = new TH2D("h2_y_vs_eta_pm", "TPC #pi^{-} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_pm = new TH2D("h2_y_vs_eta_pt0p5to1_pm", "#pi^{-} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pm = new TH2D("h2_y_vs_eta_pt1to1p5_pm", "#pi^{-} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pm = new TH2D("h2_y_vs_eta_pt1p5to2_pm", "#pi^{-} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_kp = new TH2D("h2_y_vs_eta_kp", "TPC K^{+} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_kp = new TH2D("h2_y_vs_eta_pt0p5to1_kp", "K^{+} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_kp = new TH2D("h2_y_vs_eta_pt1to1p5_kp", "K^{+} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_kp = new TH2D("h2_y_vs_eta_pt1p5to2_kp", "K^{+} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_km = new TH2D("h2_y_vs_eta_km", "TPC K^{-} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_km = new TH2D("h2_y_vs_eta_pt0p5to1_km", "K^{-} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_km = new TH2D("h2_y_vs_eta_pt1to1p5_km", "K^{-} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_km = new TH2D("h2_y_vs_eta_pt1p5to2_km", "K^{-} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_pr = new TH2D("h2_y_vs_eta_pr", "TPC Proton Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_pr = new TH2D("h2_y_vs_eta_pt0p5to1_pr", "Proton y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pr = new TH2D("h2_y_vs_eta_pt1to1p5_pr", "Proton y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pr = new TH2D("h2_y_vs_eta_pt1p5to2_pr", "Proton y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  /*
  TH2D *h2_phiSearchTpc = new TH2D("h2_phiSearchTpc", "Azimuthal Distribution by Centrality;#phi;Centrality (%)", 
				   200, -4, 4, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TH2D *h2_phiSearchEpd = new TH2D("h2_phiSearchEpd", "Azimuthal Distribution by Centrality;#phi;Centrality (%)", 
				   200, -4, 4, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  */
  // Here the name refers to the eta region that will be displayed/searched using the event plane angle from the opposite region

  TProfile2D *h2_v2SearchTpc = new TProfile2D("h2_v2SearchTpc", "<cos(2(#phi^{TPC} - #psi^{EPD}_{"+ORDER_N_STR+"}))>;#eta;Centrality (%)", 
					      12, -2, 0, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2SearchEpd = new TProfile2D("h2_v2SearchEpd", "<cos(2(#phi^{EPD} - #psi^{TPC}_{"+ORDER_N_STR+"}))>;#eta;Centrality (%)", 
					      12, -6.0, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2SearchEpdTpcA = new TProfile2D("h2_v2SearchEpdTpcA", "<cos(2(#phi^{EPD} - #psi^{TPC,A}_{"+ORDER_N_STR+"}))>;#eta;Centrality (%)", 
						  12, -6.0, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2SearchEpdTpcB = new TProfile2D("h2_v2SearchEpdTpcB", "<cos(2(#phi^{EPD} - #psi^{TPC,B}_{"+ORDER_N_STR+"}))>;#eta;Centrality (%)", 
						  12, -6.0, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  h2_v2SearchTpc->SetStats(0);
  h2_v2SearchEpd->SetStats(0);
  h2_v2SearchEpdTpcA->SetStats(0);
  h2_v2SearchEpdTpcB->SetStats(0);

  // The indices here are equivalent to the corresponding centrality ID
  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};

  // The indices here are opposite to the corresponding centrality ID (array is backward)
  //const char *centralityBins[16] = {"0-5", "5-10", "10-15" "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80"};


  TH2D *h2_psiEpdETpcA = new TH2D("h2_psiEpdETpcA", "#psi^{EPD,E} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{E}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdFTpcA = new TH2D("h2_psiEpdFTpcA", "#psi^{EPD,F} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{F}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiEpdETpcB = new TH2D("h2_psiEpdETpcB", "#psi^{EPD,E} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{E}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdFTpcB = new TH2D("h2_psiEpdFTpcB", "#psi^{EPD,F} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{F}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiTpcATpcB = new TH2D("h2_psiTpcATpcB", "#psi^{TPC,A} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{TPC}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);



  // CORRECTION HISTOGRAMS
  TProfile *p_sinAvgsTpc  = new TProfile("p_sinAvgsTpc", "Sin Averages (TPC);j (Correction term);<sin(jn#psi^{TPC}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpc  = new TProfile("p_cosAvgsTpc", "Cos Averages (TPC);j (Correction term);<sin(jn#psi^{TPC}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsTpcA = new TProfile("p_sinAvgsTpcA", "Sin Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC,A}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpcA = new TProfile("p_cosAvgsTpcA", "Cos Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC,A}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsTpcB = new TProfile("p_sinAvgsTpcB", "Sin Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC,B}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpcB = new TProfile("p_cosAvgsTpcB", "Cos Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC,B}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpd  = new TProfile("p_sinAvgsEpd", "Sin Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpd  = new TProfile("p_cosAvgsEpd", "Cos Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdE = new TProfile("p_sinAvgsEpdE", "Sin Averages (EPD E);j (Correction term);<sin(jn#psi^{EPD,E}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdE = new TProfile("p_cosAvgsEpdE", "Cos Averages (EPD E);j (Correction term);<sin(jn#psi^{EPD,E}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdF = new TProfile("p_sinAvgsEpdF", "Sin Averages (EPD F);j (Correction term);<sin(jn#psi^{EPD,F}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdF = new TProfile("p_cosAvgsEpdF", "Cos Averages (EPD F);j (Correction term);<sin(jn#psi^{EPD,F}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);

  TH1D *h_XnTpc  = new TH1D("h_XnTpc", "X_n Distribution (TPC);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpc  = new TH1D("h_YnTpc", "Y_n Distribution (TPC);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcA = new TH1D("h_XnTpcA", "X_n Distribution (TPC A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcA = new TH1D("h_YnTpcA", "Y_n Distribution (TPC A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcB = new TH1D("h_XnTpcB", "X_n Distribution (TPC B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcB = new TH1D("h_YnTpcB", "Y_n Distribution (TPC B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpd  = new TH1D("h_XnEpd", "X_n Distribution (EPD);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpd  = new TH1D("h_YnEpd", "Y_n Distribution (EPD);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdE = new TH1D("h_XnEpdE", "X_n Distribution (EPD E);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdE = new TH1D("h_YnEpdE", "Y_n Distribution (EPD E);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdF = new TH1D("h_XnEpdF", "X_n Distribution (EPD F);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdF = new TH1D("h_YnEpdF", "Y_n Distribution (EPD F);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);

  
  Int_t goodRunList[170] = { 19151031, 19151034, 19151036, 19151039, 19151041, 19151043, 19151044, 19151045,
			     19151046, 19151047, 19151048, 19151049, 19151050, 19151052, 19151053, 19151054,
			     19151055, 19151056, 19151066, 19151067, 19151068, 19151069, 19151070, 19151071,
			     19151072, 19151082, 19151083, 19151084, 19152002, 19152003, 19152008, 19152009,
			     19152010, 19152014, 19152016, 19152021, 19152023, 19152024, 19152025, 19152027,
			     19152028, 19152029, 19152030, 19152031, 19152032, 19152033, 19152034, 19152035,
			     19152036, 19152037, 19152038, 19152039, 19152040, 19152041, 19152042, 19152043,
			     19152044, 19152045, 19152046, 19152048, 19152051, 19152052, 19152053, 19152054,
			     19152055, 19152071, 19152073, 19152074, 19152075, 19152076, 19152081, 19153001,
			     19153002, 19153003, 19153004, 19153007, 19153009, 19153010, 19153011, 19153012,
			     19153013, 19153014, 19153015, 19153016, 19153017, 19153018, 19153019, 19153020,
			     19153021, 19153022, 19153024, 19153025, 19153027, 19153028, 19153029, 19153031,
			     19153033, 19153034, 19153035, 19153036, 19153037, 19153042, 19153043, 19153044,
			     19153050, 19153051, 19153052, 19153053, 19153054, 19153055, 19153056, 19153057,
			     19153058, 19153059, 19153061, 19153062, 19153063, 19153064, 19153066, 19154001,
			     19154002, 19154005, 19154007, 19154027, 19154028, 19154029, 19154030, 19154031,
			     19154032, 19154036, 19154037, 19154038, 19154039, 19154040, 19154041, 19154044,
			     19154045, 19154046, 19154047, 19154048, 19154049, 19154052, 19154053, 19154054,
			     19154055, 19154056, 19154057, 19154058, 19154061, 19154063, 19154064, 19154065,
			     19154066, 19154067, 19155001, 19155003, 19155004, 19155005, 19155006, 19155008,
			     19155009, 19155010, 19155011, 19155016, 19155017, 19155018, 19155019, 19155020,
			     19155021, 19155022 };

  Event eventInfo;
  std::vector<Event> v_events;    // Vector of all events and their info using the custom struct
  std::vector<UInt_t> triggerIDs;

  // EVENT LOOP
  for (Long64_t ievent = 0; ievent < events2read; ievent++)
    {
      eventInfo.reset();

      Bool_t readEvent = picoReader->readPicoEvent(ievent);
      if( !readEvent ) { std::cout << "Event could not be read; aborting analysis." << std::endl; break; }
        
      StPicoDst *dst = picoReader->picoDst();
      StPicoEvent *event = dst->event();
      if( !event ) { std::cout << "No event found; aborting analysis." << std::endl; break; }

      Bool_t b_bad_run = true;
      for (Int_t i = 0; i < 170; i++) { if (event->runId() == goodRunList[i]) { b_bad_run = false; break; } }
      if (b_bad_run) continue;

      h_event_check->Fill(eventSections[0], 1);

      //=========================================================
      //          Trigger Selection
      //=========================================================
      // loop for the trigger ids and see if any match minbias ID 620052

      triggerIDs.clear();
      triggerIDs = event->triggerIds();
      Bool_t b_bad_trig = true;

      for (UInt_t i = 0; i < triggerIDs.size(); i++) { if (triggerIDs[i] == 620052) {b_bad_trig = false;} } // minBias ID for 3 GeV: 620052

      if (b_bad_trig) continue;
      //=========================================================
      //      END Trigger Selection
      //=========================================================

      h_event_check->Fill(eventSections[1], 1);
            
      //=========================================================
      //          Z-VTX Selection
      //=========================================================
      // Fill vertex coordinates and check the z-vertex position

      TVector3 pVtx = event->primaryVertex();

      Double_t d_xvtx = pVtx.x();
      Double_t d_yvtx = pVtx.y();
      Double_t d_zvtx = pVtx.z();
      
      Double_t d_rvtx = TMath::Sqrt(d_xvtx * d_xvtx + (d_yvtx + 2) * (d_yvtx + 2));

      h_zvtx->Fill(d_zvtx);

      Bool_t b_bad_rvtx = ( d_rvtx >= R_VTX_CUT );
      Bool_t b_bad_zvtx = ( (d_zvtx <= Z_VTX_CUT_LOW)  || (d_zvtx >= Z_VTX_CUT_HIGH));

      if (b_bad_zvtx) continue;

      h2_trans_vtx->Fill(d_xvtx, d_yvtx);

      if (b_bad_rvtx) continue;

      h2_trans_vtx_cut->Fill(d_xvtx, d_yvtx);
      //=========================================================
      //      END Z-VTX Selection
      //=========================================================

      h_event_check->Fill(eventSections[2], 1);


      Int_t nTracks = dst->numberOfTracks();
      if (nTracks < MIN_TRACKS) continue;                // Preliminary cut to hopefully speed things up a bit. This cut repeated below also.

      // TRACK LOOP OVER PRIMARY TRACKS
      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
	{            
	  StPicoTrack *picoTrack = dst->track(iTrk);            
	  if(picoTrack == NULL) continue;
	  if(!picoTrack->isPrimary()) continue;  // Require primary tracks

	  h_track_check->Fill(trackSections[0], 1);

	  eventInfo.primTracks++;

	  //=========================================================
	  //          TOF Beta Cuts
	  //=========================================================
	  if (!picoTrack->isTofTrack()) continue;  // Only TOF tracks

	  Double_t d_tofBeta;
	      
	  StPicoBTofPidTraits *trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());
	      
	  d_tofBeta = trait->btofBeta();
	      
	  if (d_tofBeta < 0.01) continue;
	  //=========================================================
	  //          End TOF Beta Cuts
	  //=========================================================

	  h_track_check->Fill(trackSections[1], 1);

	  //=========================================================
	  //          Track QA Cuts
	  //=========================================================
	  h_nhits->Fill(picoTrack->nHits());
	  h_nhits_fit->Fill(picoTrack->nHitsFit());
	  h_nhits_dEdx->Fill(picoTrack->nHitsDedx());

	  unsigned short nHits = picoTrack->nHits();
	  
	  bool b_bad_hits     = ( nHits < 10 );
	  bool b_bad_dEdx     = ( picoTrack->nHitsDedx() <= 5 );
	  bool b_bad_tracking = ( ((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) <= 0.52 );
	  bool b_bad_DCA      = ( picoTrack->gDCA(pVtx.X(),pVtx.Y(),pVtx.Z()) >= 3.0 );

	  if (b_bad_hits || b_bad_dEdx || b_bad_tracking || b_bad_DCA) continue;
	  //=========================================================
	  //          End Track QA Cuts
	  //=========================================================

	  h_track_check->Fill(trackSections[2], 1);

	  Double_t d_TPCnSigmaPion   = picoTrack->nSigmaPion();
	  Double_t d_TPCnSigmaProton = picoTrack->nSigmaProton();
	  Double_t d_TPCnSigmaKaon   = picoTrack->nSigmaKaon();


	  TVector3 mom_vec  = picoTrack->pMom();
	  Double_t d_charge = picoTrack->charge();
	  Double_t d_dEdx   = picoTrack->dEdx();
	  Double_t d_mom    = picoTrack->pPtot();
	  Double_t d_m2     = d_mom*d_mom*( (1 / (d_tofBeta*d_tofBeta)) - 1);
	  Double_t d_px     = picoTrack->pMom().x();
	  Double_t d_py     = picoTrack->pMom().y();
	  Double_t d_pz     = picoTrack->pMom().z();
	  Double_t d_pT     = picoTrack->pPt();
	  Double_t d_phi    = mom_vec.Phi();
	  Double_t d_eta    = mom_vec.Eta();


	  // Fill histos and save important event info in the custom struct type
	  if (d_charge != 0)
	    {
	      h_eta_TPC_s->Fill(d_eta - Y_MID);

	      eventInfo.nTracksTpc++;
	      if (d_eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
		{
		  eventInfo.XnTpc += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
		  eventInfo.YnTpc += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
		}
	      else if (d_eta < Y_MID)
		{
		  eventInfo.XnTpc -= d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
		  eventInfo.YnTpc -= d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
		}


	      if (d_eta > MIN_TPC_ETA_CUT && d_eta < AGAP_TPC_ETA_CUT)          // TPC A  (sign change happens later
		{
		  eventInfo.nTracksTpcA++;
		  eventInfo.XnTpcA += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
		  eventInfo.YnTpcA += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
		  eventInfo.phiValuesTpcA.push_back(d_phi);
		  eventInfo.etaValuesTpcA.push_back(d_eta);
		}
	      else if (d_eta > GAPB_TPC_ETA_CUT && d_eta < MAX_TPC_ETA_CUT)     // TPC B
		{
		  eventInfo.nTracksTpcB++;
		  eventInfo.XnTpcB += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
		  eventInfo.YnTpcB += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
		  eventInfo.phiValuesTpcB.push_back(d_phi);
		  eventInfo.etaValuesTpcB.push_back(d_eta);
		}


	      h2_betap->Fill(d_charge * d_mom, 1/d_tofBeta);
	      h2_m2_p->Fill(d_charge * d_mom, d_mom * d_mom * (1/(d_tofBeta*d_tofBeta) - 1));
	      h2_m2_vs_qpT->Fill(d_charge * d_pT, d_m2);
	      h2_dEdx_vs_qp->Fill(d_charge * d_mom, d_dEdx);

	      h2_pi_m2_vs_TPC_nsig->Fill(d_TPCnSigmaPion, d_m2);
	      h2_ka_m2_vs_TPC_nsig->Fill(d_TPCnSigmaKaon, d_m2);
	      h2_pr_m2_vs_TPC_nsig->Fill(d_TPCnSigmaProton, d_m2);

	      //=========================================================
	      //          PID Cuts
	      //=========================================================
	      Bool_t pion   = ( (d_m2 > -0.1) && (d_m2 < 0.1)   && (d_TPCnSigmaPion > -2)   && (d_TPCnSigmaPion < 2) );    // PARTICLE TAGGING by TPC nSigma values
	      Bool_t kaon   = ( (d_m2 > 0.15) && (d_m2 < 0.34)  && (d_TPCnSigmaKaon > -2)   && (d_TPCnSigmaKaon < 2) );
	      Bool_t proton = ( (d_m2 > 0.65) && (d_m2 < 1.11)  && (d_TPCnSigmaProton > -2) && (d_TPCnSigmaProton < 2) );

	      if (!pion && !kaon && !proton) continue;

	      if (pion && kaon)   continue;   // Particle must be exclusively one particular type.
	      if (pion && proton) continue;
	      if (kaon && proton) continue;

	      if (pion)   h2_nSig_vs_qp_pi->Fill(d_charge * d_mom, d_TPCnSigmaPion);
	      if (kaon)   h2_nSig_vs_qp_ka->Fill(d_charge * d_mom, d_TPCnSigmaKaon);
	      if (proton) h2_nSig_vs_qp_pr->Fill(d_charge * d_mom, d_TPCnSigmaProton);
	      //=========================================================
	      //          END PID Cuts
	      //=========================================================

	      h_track_check->Fill(trackSections[3], 1);

	      Double_t d_m0_pi = 0.1396;   //Rest masses
	      Double_t d_m0_ka = 0.4937;
	      Double_t d_m0_pr = 0.9383;


	      if (pion)
		{
		  if (d_charge > 0) 
		    {
		      if (d_mom > PI_KA_MOM_CUT) continue;
		      
		      eventInfo.phiValuesPionP.push_back(d_phi);

		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h2_y_vs_eta_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h_pp_pT->Fill(d_pT);
		      h_pp_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_pi, h_pp_dndy, h_pp_dndm, h2_pp_MvsY);

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		    }
		  else if (d_charge < 0) 
		    {
		      if (d_mom > PI_KA_MOM_CUT) continue;

		      eventInfo.phiValuesPionM.push_back(d_phi);

		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h2_y_vs_eta_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h_pm_pT->Fill(d_pT);
		      h_pm_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_pi, h_pm_dndy, h_pm_dndm, h2_pm_MvsY);

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		    }
		}
	      else if (kaon)
		{
		  if (d_charge > 0) 
		    {
		      if (d_mom > PI_KA_MOM_CUT) continue;

		      eventInfo.phiValuesKaonP.push_back(d_phi);

		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h2_y_vs_eta_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h_kp_pT->Fill(d_pT);
		      h_kp_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_ka, h_kp_dndy, h_kp_dndm, h2_kp_MvsY);

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		    }
		  else if (d_charge < 0)		 
		    {
		      if (d_mom > PI_KA_MOM_CUT) continue;

		      eventInfo.phiValuesKaonM.push_back(d_phi);

		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h2_y_vs_eta_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h_km_pT->Fill(d_pT);
		      h_km_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_ka, h_km_dndy, h_km_dndm, h2_km_MvsY);

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		    }
		}
	      else if (proton)
		{
		  if (d_charge > 0) 
		    {
		      if (d_mom > PR_MOM_CUT) continue;

		      eventInfo.phiValuesProton.push_back(d_phi);

		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr));
		      h2_y_vs_eta_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr));
		      h_pr_pT->Fill(d_pT);
		      h_pr_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_pr, h_pr_dndy, h_pr_dndm, h2_pr_MvsY);

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr)); }
		    }
		}

	    }// End if(d_charge != 0)
	}//End track loop


      // ASSIGN CENTRALITY ID
      if     ( eventInfo.primTracks >=   3 && eventInfo.primTracks <=   4 ) eventInfo.centID =  0;  // 75% - 80% (Peripheral)
      else if( eventInfo.primTracks >=   5 && eventInfo.primTracks <=   6 ) eventInfo.centID =  1;
      else if( eventInfo.primTracks >=   7 && eventInfo.primTracks <=   9 ) eventInfo.centID =  2;
      else if( eventInfo.primTracks >=  10 && eventInfo.primTracks <=  13 ) eventInfo.centID =  3;
      else if( eventInfo.primTracks >=  14 && eventInfo.primTracks <=  17 ) eventInfo.centID =  4;
      else if( eventInfo.primTracks >=  18 && eventInfo.primTracks <=  23 ) eventInfo.centID =  5;
      else if( eventInfo.primTracks >=  24 && eventInfo.primTracks <=  29 ) eventInfo.centID =  6;
      else if( eventInfo.primTracks >=  30 && eventInfo.primTracks <=  37 ) eventInfo.centID =  7;
      else if( eventInfo.primTracks >=  38 && eventInfo.primTracks <=  46 ) eventInfo.centID =  8;
      else if( eventInfo.primTracks >=  47 && eventInfo.primTracks <=  56 ) eventInfo.centID =  9;
      else if( eventInfo.primTracks >=  57 && eventInfo.primTracks <=  68 ) eventInfo.centID = 10;
      else if( eventInfo.primTracks >=  69 && eventInfo.primTracks <=  82 ) eventInfo.centID = 11;
      else if( eventInfo.primTracks >=  83 && eventInfo.primTracks <=  98 ) eventInfo.centID = 12;
      else if( eventInfo.primTracks >=  99 && eventInfo.primTracks <= 117 ) eventInfo.centID = 13;
      else if( eventInfo.primTracks >= 118 && eventInfo.primTracks <= 140 ) eventInfo.centID = 14;
      else if( eventInfo.primTracks >= 141 && eventInfo.primTracks <= 195 ) eventInfo.centID = 15;  // 0% - 5% (Central)

      if (eventInfo.centID == I_BAD_VALUE) continue;
      
      if (eventInfo.centID < FIRST_CENT) continue;

      //=========================================================
      //                EPD STUFF
      //=========================================================
      //StEpdEpInfo result = epdEpFinder->Results(epdHits,pVtx,eventInfo.centID);
      StEpdGeom *epdGeom = new StEpdGeom();
      StPicoEpdHit *epdHit;
      int tileID;
      TVector3 tileVector;     // Vector from vertex to center of tile that was hit
      int tileSector;
      //int tileRow;
      double tileWeight;
      double tileEta;
      double tilePhi;

      for (int iEpdHit = 0; iEpdHit < epdHits->GetEntries(); iEpdHit++)
	{
	  epdHit = (StPicoEpdHit*)(epdHits->At(iEpdHit));

	  tileID = epdHit->id();
	  if (tileID > 0) continue;      // Exclude the West side
	  
	  tileVector = epdGeom->TileCenter(tileID) - pVtx;
	  tileSector = epdHit->position();
	  //tileRow = epdHit->row();
	  tileEta = tileVector.Eta();
	  tilePhi = tileVector.Phi();
	  tileWeight = (epdHit->nMIP() > EPD_THRESHOLD) ? ( (epdHit->nMIP() > EPD_MAX_WEIGHT)?EPD_MAX_WEIGHT:epdHit->nMIP() ) : 0;

	  h2_pp_vs_eta->Fill(tileEta, tileSector, tileWeight);
	  h_tileWeights->Fill(tileWeight);

	  h_eta_s->Fill(tileEta - Y_MID);

	  eventInfo.nhitsEpd++;
	  eventInfo.XnEpd += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	  eventInfo.YnEpd += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);

	  if (tileEta > MIN_ETA_CUT && tileEta < EF_ETA_CUT)  // Sub E
	    {
	      eventInfo.nhitsEpdE++;
	      eventInfo.XnEpdE += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdE += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdE.push_back(tilePhi);
	      eventInfo.etaValuesEpdE.push_back(tileEta);
	    }
	  else if (tileEta >= EF_ETA_CUT && tileEta < MAX_ETA_CUT)  // Sub F
	    {
	      eventInfo.nhitsEpdF++;
	      eventInfo.XnEpdF += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdF += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdF.push_back(tilePhi);
	      eventInfo.etaValuesEpdF.push_back(tileEta);
	    }

	}
      delete epdGeom;
      //=========================================================
      //            END EPD STUFF
      //=========================================================


      if (eventInfo.nTracksTpc  < MIN_TRACKS) continue;
      //if (eventInfo.nTracksTpcA < MIN_TRACKS) continue;
      if (eventInfo.nTracksTpcB < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpd    < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdE   < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdF   < MIN_TRACKS) continue;
      

      checkZeroQ(eventInfo);
      if (eventInfo.badEvent) continue;


      // RAW SUB-EVENT PLANE ANGLES //
      if (ORDER_N == 1)           // Q vectors must change sign past mid-rapidity; I think this is for 1st order event-planes only. Full TPC already takes this into account.
	{
	  eventInfo.XnTpcA *= -1.0;
	  eventInfo.YnTpcA *= -1.0;
	  eventInfo.XnEpd  *= -1.0;
	  eventInfo.YnEpd  *= -1.0;
	  eventInfo.XnEpdE *= -1.0;
	  eventInfo.YnEpdE *= -1.0;
	  eventInfo.XnEpdF *= -1.0;
	  eventInfo.YnEpdF *= -1.0;
	}
      eventInfo.psiTpc  = TMath::ATan2(eventInfo.YnTpc,  eventInfo.XnTpc)  / (Double_t)ORDER_N;
      eventInfo.psiTpcA = TMath::ATan2(eventInfo.YnTpcA, eventInfo.XnTpcA) / (Double_t)ORDER_N;
      eventInfo.psiTpcB = TMath::ATan2(eventInfo.YnTpcB, eventInfo.XnTpcB) / (Double_t)ORDER_N;
      eventInfo.psiEpd  = TMath::ATan2(eventInfo.YnEpd,  eventInfo.XnEpd)  / (Double_t)ORDER_N;
      eventInfo.psiEpdE = TMath::ATan2(eventInfo.YnEpdE, eventInfo.XnEpdE) / (Double_t)ORDER_N;
      eventInfo.psiEpdF = TMath::ATan2(eventInfo.YnEpdF, eventInfo.XnEpdF) / (Double_t)ORDER_N;


      // Fill eta distribution here since this is past all possible cuts.

      for (unsigned int i = 0; i < eventInfo.etaValuesTpcA.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcA.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcB.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcB.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdE.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdE.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdF.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdF.at(i) - Y_MID ); } 


      h_primTracks->Fill(eventInfo.primTracks);
      h_centralities->Fill(eventInfo.centID);

      h_XnTpc->Fill(eventInfo.XnTpc);
      h_YnTpc->Fill(eventInfo.YnTpc);
      h_XnTpcA->Fill(eventInfo.XnTpcA);
      h_YnTpcA->Fill(eventInfo.YnTpcA);
      h_XnTpcB->Fill(eventInfo.XnTpcB);
      h_YnTpcB->Fill(eventInfo.YnTpcB);
      h_XnEpd->Fill(eventInfo.XnEpd);
      h_YnEpd->Fill(eventInfo.YnEpd);
      h_XnEpdE->Fill(eventInfo.XnEpdE);
      h_YnEpdE->Fill(eventInfo.YnEpdE);
      h_XnEpdF->Fill(eventInfo.XnEpdF);
      h_YnEpdF->Fill(eventInfo.YnEpdF);

      h_psiTpc_RAW->Fill(eventInfo.psiTpc);
      h_psiTpcA_RAW->Fill(eventInfo.psiTpcA);
      h_psiTpcB_RAW->Fill(eventInfo.psiTpcB);
      h_psiEpd_RAW->Fill(eventInfo.psiEpd);
      h_psiEpdE_RAW->Fill(eventInfo.psiEpdE);
      h_psiEpdF_RAW->Fill(eventInfo.psiEpdF);

      v_events.push_back(eventInfo);   // Store this event with all of its attributes
    }//End event loop

  eventInfo.reset();

  TH1D *h_XnTpc_RC;   // Re-centered histograms
  TH1D *h_YnTpc_RC;
  TH1D *h_XnTpcA_RC;
  TH1D *h_YnTpcA_RC;
  TH1D *h_XnTpcB_RC;
  TH1D *h_YnTpcB_RC;
  TH1D *h_XnEpd_RC;
  TH1D *h_YnEpd_RC;
  TH1D *h_XnEpdE_RC;
  TH1D *h_YnEpdE_RC;
  TH1D *h_XnEpdF_RC;
  TH1D *h_YnEpdF_RC;

  TH1D *h_psiTpc_RC;
  TH1D *h_psiTpcA_RC;
  TH1D *h_psiTpcB_RC;
  TH1D *h_psiEpd_RC;
  TH1D *h_psiEpdE_RC;
  TH1D *h_psiEpdF_RC;

  TH1D *h_psiTpc_FLAT;
  TH1D *h_psiTpcA_FLAT;
  TH1D *h_psiTpcB_FLAT;
  TH1D *h_psiEpd_FLAT;
  TH1D *h_psiEpdE_FLAT;
  TH1D *h_psiEpdF_FLAT;


  //=========================================================
  //          Re-centering (Xn, Yn) Distributions
  //=========================================================

  if (RUN_ITERATION == 1 || RUN_ITERATION == 2)
    {
      std::cout << "Re-centering flow vectors and accumulating sin/cos averages..." << std::endl;

      h_XnTpc_RC  = new TH1D("h_XnTpc_RC", "Re-centered X_n Distribution (TPC);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnTpc_RC  = new TH1D("h_YnTpc_RC", "Re-centered Y_n Distribution (TPC);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnTpcA_RC = new TH1D("h_XnTpcA_RC", "Re-centered X_n Distribution (TPC A);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnTpcA_RC = new TH1D("h_YnTpcA_RC", "Re-centered Y_n Distribution (TPC A);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnTpcB_RC = new TH1D("h_XnTpcB_RC", "Re-centered X_n Distribution (TPC B);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnTpcB_RC = new TH1D("h_YnTpcB_RC", "Re-centered Y_n Distribution (TPC B);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpd_RC  = new TH1D("h_XnEpd_RC", "Re-centered X_n Distribution (EPD);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpd_RC  = new TH1D("h_YnEpd_RC", "Re-centered Y_n Distribution (EPD);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdE_RC = new TH1D("h_XnEpdE_RC", "Re-centered X_n Distribution (EPD E);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdE_RC = new TH1D("h_YnEpdE_RC", "Re-centered Y_n Distribution (EPD E);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdF_RC = new TH1D("h_XnEpdF_RC", "Re-centered X_n Distribution (EPD F);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdF_RC = new TH1D("h_YnEpdF_RC", "Re-centered Y_n Distribution (EPD F);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);

      h_psiTpc_RC  = new TH1D("h_psiTpc_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiTpcA_RC = new TH1D("h_psiTpcA_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiTpcB_RC = new TH1D("h_psiTpcB_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpd_RC  = new TH1D("h_psiEpd_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdE_RC = new TH1D("h_psiEpdE_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD E);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdF_RC = new TH1D("h_psiEpdF_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD F);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

      TH1D *h_XnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_XnTpc");
      TH1D *h_XnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcA");
      TH1D *h_XnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcB");
      TH1D *h_XnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_XnEpd");
      TH1D *h_XnEpdE_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdE");
      TH1D *h_XnEpdF_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdF");

      TH1D *h_YnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_YnTpc");
      TH1D *h_YnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcA");
      TH1D *h_YnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcB");
      TH1D *h_YnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_YnEpd");
      TH1D *h_YnEpdE_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdE");
      TH1D *h_YnEpdF_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdF");

      Double_t d_XnTpc_Avg  = h_XnTpc_INPUT->GetMean();
      Double_t d_XnTpcA_Avg = h_XnTpcA_INPUT->GetMean();
      Double_t d_XnTpcB_Avg = h_XnTpcB_INPUT->GetMean();
      Double_t d_XnEpd_Avg  = h_XnEpd_INPUT->GetMean();
      Double_t d_XnEpdE_Avg = h_XnEpdE_INPUT->GetMean();
      Double_t d_XnEpdF_Avg = h_XnEpdF_INPUT->GetMean();

      Double_t d_YnTpc_Avg  = h_YnTpc_INPUT->GetMean();
      Double_t d_YnTpcA_Avg = h_YnTpcA_INPUT->GetMean();
      Double_t d_YnTpcB_Avg = h_YnTpcB_INPUT->GetMean();
      Double_t d_YnEpd_Avg  = h_YnEpd_INPUT->GetMean();
      Double_t d_YnEpdE_Avg = h_YnEpdE_INPUT->GetMean();
      Double_t d_YnEpdF_Avg = h_YnEpdF_INPUT->GetMean();


      Int_t numOfEvents = v_events.size();
      Int_t badEvents = 0;

      for (int i = 0; i < numOfEvents; i++)
	{
	  v_events.at(i).XnTpc  -= d_XnTpc_Avg;
	  v_events.at(i).XnTpcA -= d_XnTpcA_Avg;
	  v_events.at(i).XnTpcB -= d_XnTpcB_Avg;
	  v_events.at(i).XnEpd  -= d_XnEpd_Avg;
	  v_events.at(i).XnEpdE -= d_XnEpdE_Avg;
	  v_events.at(i).XnEpdF -= d_XnEpdF_Avg;

	  v_events.at(i).YnTpc  -= d_YnTpc_Avg;
	  v_events.at(i).YnTpcA -= d_YnTpcA_Avg;
	  v_events.at(i).YnTpcB -= d_YnTpcB_Avg;
	  v_events.at(i).YnEpd  -= d_YnEpd_Avg;
	  v_events.at(i).YnEpdE -= d_YnEpdE_Avg;
	  v_events.at(i).YnEpdF -= d_YnEpdF_Avg;

	  checkZeroQ(v_events.at(i));

	  if ( v_events.at(i).badEvent ) { badEvents++; continue; }

	  h_XnTpc_RC->Fill(v_events.at(i).XnTpc);
	  h_XnTpcA_RC->Fill(v_events.at(i).XnTpcA);
	  h_XnTpcB_RC->Fill(v_events.at(i).XnTpcB);
	  h_XnEpd_RC->Fill(v_events.at(i).XnEpd);
	  h_XnEpdE_RC->Fill(v_events.at(i).XnEpdE);
	  h_XnEpdF_RC->Fill(v_events.at(i).XnEpdF);

	  h_YnTpc_RC->Fill(v_events.at(i).YnTpc);
	  h_YnTpcA_RC->Fill(v_events.at(i).YnTpcA);
	  h_YnTpcB_RC->Fill(v_events.at(i).YnTpcB);
	  h_YnEpd_RC->Fill(v_events.at(i).YnEpd);
	  h_YnEpdE_RC->Fill(v_events.at(i).YnEpdE);
	  h_YnEpdF_RC->Fill(v_events.at(i).YnEpdF);

	  // Recalculate the event plane angles after re-centering	  
	  v_events.at(i).psiTpc  = TMath::ATan2(v_events.at(i).YnTpc,  v_events.at(i).XnTpc)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiTpcA = TMath::ATan2(v_events.at(i).YnTpcA, v_events.at(i).XnTpcA) / (Double_t)ORDER_N; 
	  v_events.at(i).psiTpcB = TMath::ATan2(v_events.at(i).YnTpcB, v_events.at(i).XnTpcB) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpd  = TMath::ATan2(v_events.at(i).YnEpd,  v_events.at(i).XnEpd)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdE = TMath::ATan2(v_events.at(i).YnEpdE, v_events.at(i).XnEpdE) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdF = TMath::ATan2(v_events.at(i).YnEpdF, v_events.at(i).XnEpdF) / (Double_t)ORDER_N; 

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_N);
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_N);
	  v_events.at(i).psiEpdE = angleShift(v_events.at(i).psiEpdE, ORDER_N);
	  v_events.at(i).psiEpdF = angleShift(v_events.at(i).psiEpdF, ORDER_N);

	  h_psiTpc_RC->Fill(v_events.at(i).psiTpc);
	  h_psiTpcA_RC->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_RC->Fill(v_events.at(i).psiTpcB);
	  h_psiEpd_RC->Fill(v_events.at(i).psiEpd);
	  h_psiEpdE_RC->Fill(v_events.at(i).psiEpdE);
	  h_psiEpdF_RC->Fill(v_events.at(i).psiEpdF);


	  // Accumulate terms for averages over the re-centered angles for event plane angle shifting
	  for (int j = 1; j <= SHIFT_TERMS; j++)
	    {
	      p_sinAvgsTpc->Fill(j,  TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc));
	      p_cosAvgsTpc->Fill(j,  TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc));
	      p_sinAvgsTpcA->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      p_cosAvgsTpcA->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      p_sinAvgsTpcB->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      p_cosAvgsTpcB->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      p_sinAvgsEpd->Fill(j,  TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd));
	      p_cosAvgsEpd->Fill(j,  TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd));
	      p_sinAvgsEpdE->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdE));
	      p_cosAvgsEpdE->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdE));
	      p_sinAvgsEpdF->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdF));
	      p_cosAvgsEpdF->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdF));
	    }

	}// End loop over v_events
      std::cout << "Bad Events after re-centering: " << badEvents << std::endl;
    }
  //=========================================================
  //          End Re-centering
  //=========================================================


  //=========================================================
  //          Event Plane Angle Shifting
  //=========================================================

  if (RUN_ITERATION == 2)
    {
      std::cout << "Performing event plane angle shifting..." << std::endl;

      Int_t numOfEvents = v_events.size();

      h_psiTpc_FLAT  = new TH1D("h_psiTpc_FLAT", "Flattened Event Plane Angle (TPC, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiTpcA_FLAT = new TH1D("h_psiTpcA_FLAT", "Flattened Event Plane Angle (TPC A, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiTpcB_FLAT = new TH1D("h_psiTpcB_FLAT", "Flattened Event Plane Angle (TPC B, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiEpd_FLAT  = new TH1D("h_psiEpd_FLAT", "Flattened Event Plane Angle (EPD, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdE_FLAT = new TH1D("h_psiEpdE_FLAT", "Flattened Event Plane Angle (EPD E, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdF_FLAT = new TH1D("h_psiEpdF_FLAT", "Flattened Event Plane Angle (EPD F, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

      TProfile *p_sinAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsTpc");
      TProfile *p_cosAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsTpc");
      TProfile *p_sinAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcA");
      TProfile *p_cosAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcA");
      TProfile *p_sinAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcB");
      TProfile *p_cosAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcB");
      TProfile *p_sinAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsEpd");
      TProfile *p_cosAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsEpd");
      TProfile *p_sinAvgsEpdE_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdE");
      TProfile *p_cosAvgsEpdE_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdE");
      TProfile *p_sinAvgsEpdF_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdF");
      TProfile *p_cosAvgsEpdF_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdF");


      // Get corrected event plane angles //
      for (Int_t i = 0; i < numOfEvents; i++)  // Loop over v_events to correct their angles
	{
	  if ( v_events.at(i).badEvent == true) { continue; }

	  Double_t psiTpc_delta  = 0;
	  Double_t psiTpcA_delta = 0;
	  Double_t psiTpcB_delta = 0;
	  Double_t psiEpd_delta  = 0;
	  Double_t psiEpdE_delta = 0;
	  Double_t psiEpdF_delta = 0;

	  Double_t jthSinAvg_Tpc  = 0;
	  Double_t jthCosAvg_Tpc  = 0;
	  Double_t jthSinAvg_TpcA = 0;
	  Double_t jthCosAvg_TpcA = 0;
	  Double_t jthSinAvg_TpcB = 0;
	  Double_t jthCosAvg_TpcB = 0;
	  Double_t jthSinAvg_Epd  = 0;
	  Double_t jthCosAvg_Epd  = 0;
	  Double_t jthSinAvg_EpdE = 0;
	  Double_t jthCosAvg_EpdE = 0;
	  Double_t jthSinAvg_EpdF = 0;
	  Double_t jthCosAvg_EpdF = 0;


	  for (Int_t j = 1; j <= SHIFT_TERMS; j++)    // Build the correction sums
	    {
	      jthSinAvg_Tpc  = p_sinAvgsTpc_INPUT->GetBinContent(j);
	      jthCosAvg_Tpc  = p_cosAvgsTpc_INPUT->GetBinContent(j);
	      jthSinAvg_TpcA = p_sinAvgsTpcA_INPUT->GetBinContent(j);
	      jthCosAvg_TpcA = p_cosAvgsTpcA_INPUT->GetBinContent(j);
	      jthSinAvg_TpcB = p_sinAvgsTpcB_INPUT->GetBinContent(j);
	      jthCosAvg_TpcB = p_cosAvgsTpcB_INPUT->GetBinContent(j);
	      jthSinAvg_Epd  = p_sinAvgsEpd_INPUT->GetBinContent(j);
	      jthCosAvg_Epd  = p_cosAvgsEpd_INPUT->GetBinContent(j);
	      jthSinAvg_EpdE = p_sinAvgsEpdE_INPUT->GetBinContent(j);
	      jthCosAvg_EpdE = p_cosAvgsEpdE_INPUT->GetBinContent(j);
	      jthSinAvg_EpdF = p_sinAvgsEpdF_INPUT->GetBinContent(j);
	      jthCosAvg_EpdF = p_cosAvgsEpdF_INPUT->GetBinContent(j);

	      psiTpc_delta  += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_Tpc*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc) 
									+jthCosAvg_Tpc*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc));
	      psiTpcA_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcA*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA) 
									+jthCosAvg_TpcA*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      psiTpcB_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcB*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB) 
									+jthCosAvg_TpcB*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      psiEpd_delta  += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_Epd*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd)
									+jthCosAvg_Epd*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd));
	      psiEpdE_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdE*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdE)
									+jthCosAvg_EpdE*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdE));
	      psiEpdF_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdF*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdF)
									+jthCosAvg_EpdF*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdF));
	    }


	  v_events.at(i).psiTpc  += psiTpc_delta;
	  v_events.at(i).psiTpcA += psiTpcA_delta;
	  v_events.at(i).psiTpcB += psiTpcB_delta;
	  v_events.at(i).psiEpd  += psiEpd_delta;
	  v_events.at(i).psiEpdE += psiEpdE_delta;
	  v_events.at(i).psiEpdF += psiEpdF_delta;

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_N);
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_N);
	  v_events.at(i).psiEpdE = angleShift(v_events.at(i).psiEpdE, ORDER_N);
	  v_events.at(i).psiEpdF = angleShift(v_events.at(i).psiEpdF, ORDER_N);

	  h_psiTpc_FLAT->Fill(v_events.at(i).psiTpc);
	  h_psiTpcA_FLAT->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_FLAT->Fill(v_events.at(i).psiTpcB);
	  h_psiEpd_FLAT->Fill(v_events.at(i).psiEpd);
	  h_psiEpdE_FLAT->Fill(v_events.at(i).psiEpdE);
	  h_psiEpdF_FLAT->Fill(v_events.at(i).psiEpdF);


	  // 2D Correlations between event planes
	  h2_psiEpdETpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdE);
	  h2_psiEpdFTpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdF);

	  h2_psiEpdETpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdE);
	  h2_psiEpdFTpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdF);

	  h2_psiTpcATpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiTpcA);
	  //


	  // 1D correlation averages used in calculating resolution using the 3 sub-event method
	  p_TpcAB->Fill(v_events.at(i).centID,    TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiTpcB)));

	  p_TpcAEpdE->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdE)));
	  p_TpcAEpdF->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdF)));
	  p_TpcBEpdE->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdE)));
	  p_TpcBEpdF->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdF)));

	  p_EpdEEpdF->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiEpdE - v_events.at(i).psiEpdF)));
	  //

	  // 2D searches through eta and centrality for correlations between detectors
	  // ONLY USE THIS SECTION IF THE EPD REGIONS COVER THE WHOLE EPD!! MIGHT NOT MAKE SENSE OTHERWISE

	  Int_t tpcTracksA = v_events.at(i).phiValuesTpcA.size();
	  Int_t tpcTracksB = v_events.at(i).phiValuesTpcB.size();
	  Int_t epdHitsE   = v_events.at(i).phiValuesEpdE.size();
	  Int_t epdHitsF   = v_events.at(i).phiValuesEpdF.size();
	  Double_t phiTpc;
	  Double_t etaTpc;
	  Double_t phiEpd;
	  Double_t etaEpd;
	  Double_t psiTpc  = v_events.at(i).psiTpc;
	  Double_t psiEpd  = v_events.at(i).psiEpd;
	  Double_t psiTpcA = v_events.at(i).psiTpcA;
	  Double_t psiTpcB = v_events.at(i).psiTpcB;
	  Int_t centralityID = v_events.at(i).centID;

	  for (int k = 0; k < epdHitsE; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdE.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdE.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsF; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdF.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdF.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }

	  for (int k = 0; k < tpcTracksA; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcA.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcA.at(k);

	      h2_v2SearchTpc->Fill(etaTpc, centralityID, TMath::Cos((Double_t)ORDER_N * (phiTpc - psiEpd)));
	      //h2_phiSearchTpc->Fill(phiTpc, centralityID);
	    }
	  for (int k = 0; k < tpcTracksB; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcB.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcB.at(k);

	      h2_v2SearchTpc->Fill(etaTpc, centralityID, TMath::Cos((Double_t)ORDER_N * (phiTpc - psiEpd)));
	      //h2_phiSearchTpc->Fill(phiTpc, centralityID);
	    }

	  v_events.at(i).reset(); // Try to free up space?


	  //=========================================================
	  //          Get v_2 values
	  //=========================================================

	  Double_t cosTerm = 0;
	  Double_t phi = 0;
	  Double_t psi = 0;
	  /*
	  // CORRELATIONS AND FLOW COEFFICIENTS
	  for (Int_t i = 0; i < numOfEvents; i++) //Don't need this loop
	    {
	      psi = v_events.at(i).psi;


	    }
	  */
	}// End shift loop over events
    }
  //=========================================================
  //          End Event Plane Angle Shifting
  //=========================================================

  // Switch y-axis labels to centrality percentages

  Int_t labelIndex;
  for (int i = 1; i <= CENT_BINS; i++) 
    {
      labelIndex = FIRST_CENT + i - 1;
      h2_v2SearchTpc->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      h2_v2SearchEpd->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_v2SearchEpdTpcA->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_v2SearchEpdTpcB->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      //h2_phiSearchTpc->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      //h2_phiSearchEpd->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
    }



  outputFile->cd();
  outputFile->Write();

  if (RUN_ITERATION == 0 || RUN_ITERATION == 1)
    {
      correctionOutputFile->cd();

      p_sinAvgsTpc   ->Write();
      p_cosAvgsTpc   ->Write();
      p_sinAvgsTpcA  ->Write();
      p_cosAvgsTpcA  ->Write();
      p_sinAvgsTpcB  ->Write();
      p_cosAvgsTpcB  ->Write();
      p_sinAvgsEpd   ->Write();
      p_cosAvgsEpd   ->Write();
      p_sinAvgsEpdE  ->Write();
      p_cosAvgsEpdE  ->Write();
      p_sinAvgsEpdF  ->Write();
      p_cosAvgsEpdF  ->Write();
      h_XnTpc        ->Write();
      h_YnTpc        ->Write();
      h_XnTpcA       ->Write();
      h_YnTpcA       ->Write();
      h_XnTpcB       ->Write();
      h_YnTpcB       ->Write();
      h_XnEpd        ->Write();
      h_YnEpd        ->Write();
      h_XnEpdE       ->Write();
      h_YnEpdE       ->Write();
      h_XnEpdF       ->Write();
      h_YnEpdF       ->Write();

      gROOT->GetListOfFiles()->Remove(correctionOutputFile);
      correctionOutputFile->Close();
    }

  gROOT->GetListOfFiles()->Remove(outputFile);
  outputFile->Close();

  if (RUN_ITERATION == 1 || RUN_ITERATION == 2)
    {
      gROOT->GetListOfFiles()->Remove(correctionInputFile);
      correctionInputFile->Close();
    }

  //epdEpFinder->Finish();
  picoReader->Finish();

  std::cout << "Done!" << std::endl;
}//End FlowAnalyzer()





      /* REMOVING AUTO-CORRELATIONS */
      /*
	for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	{
	Double_t newXn = v_events.at(i).Xn - v_events.at(i).pTValues.at(j) * TMath::Cos((Double_t)ORDER_N * v_events.at(i).phiValues.at(j));
	Double_t newYn = v_events.at(i).Yn - v_events.at(i).pTValues.at(j) * TMath::Sin((Double_t)ORDER_N * v_events.at(i).phiValues.at(j));
	Double_t newPsi = TMath::ATan2(newYn, newXn) / (Double_t)ORDER_N;
	v_events.at(i).psiValues.push_back(newPsi);
	}	  

	      // GET V2 WITHOUT AUTOCORRELATIONS HERE
	      for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	      {
	      phi = v_events.at(i).phiValues.at(j);
	      psi = v_events.at(i).psiValues.at(j);  // This psi was calculated without particle 'j'
	  
	      cos = TMath::Cos(2 * (phi - psi));
	      h_v2Plot->Fill(cos);
	      } 

      */
