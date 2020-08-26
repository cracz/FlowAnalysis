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

const Double_t MIN_ETA_CUT = -5.6;//-5.16;    // Cuts for the EPD subevents (A, B, C, D)
const Double_t AB_ETA_CUT  = -4.05;//-3.82;
const Double_t BC_ETA_CUT  = -3.3;//-3.28;
const Double_t CD_ETA_CUT  = -2.9;//-2.87;
const Double_t MAX_ETA_CUT = -2.30;
/*
const Double_t C1_MIN_ETA_CUT = -4.0;       // Possble choices for sub-event C in the EPD
const Double_t C1_MAX_ETA_CUT = -3.0;
const Double_t C2_MIN_ETA_CUT = -4.0;
const Double_t C2_MAX_ETA_CUT = -2.3;
const Double_t C3_MIN_ETA_CUT = -4.4;
const Double_t C3_MAX_ETA_CUT = -2.3;
*/
const Double_t R_VTX_CUT = 2.0;         // 2D r value, good vertices are within this value
const Double_t Z_VTX_CUT_LOW  = 199.5;
const Double_t Z_VTX_CUT_HIGH = 201.5;

const Int_t MIN_TRACKS = 5;             // Min number of tracks/hits in each sub-event

const Int_t SHIFT_TERMS = 10;           // Number of terms to use when shifting event plane angles
const Int_t CENT_BINS = 8;             // Number of centrality bins (max 16)
const Int_t FIRST_CENT = 16 - CENT_BINS;            // Starting point for centrality dependent plots
const Int_t BAD_VALUE = -99;

const Double_t Y_MID    = -1.05;       // Mid rapidity
const Double_t Y_MID_PI_P = -0.82;       // Mid rapidity for pions (plus and minus), kaons, and protons
const Double_t Y_MID_PI_M = -0.82;
const Double_t Y_MID_KA_P = -0.76;
const Double_t Y_MID_KA_M = -0.80;
const Double_t Y_MID_PR = -0.60;

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

  Int_t nPionP;
  Double_t XnPionP;
  Double_t YnPionP;
  Double_t psiPionP;

  Int_t nPionM;
  Double_t XnPionM;
  Double_t YnPionM;
  Double_t psiPionM;

  Int_t nKaonP;
  Double_t XnKaonP;
  Double_t YnKaonP;
  Double_t psiKaonP;

  Int_t nKaonM;
  Double_t XnKaonM;
  Double_t YnKaonM;
  Double_t psiKaonM;

  Int_t nProton;
  Double_t XnProton;
  Double_t YnProton;
  Double_t psiProton;

  Int_t nhitsEpd;
  Double_t XnEpd;
  Double_t YnEpd;
  Double_t psiEpd;

  Int_t nhitsEpdA;
  Double_t XnEpdA;
  Double_t YnEpdA;
  Double_t psiEpdA;
  std::vector<Double_t> phiValuesEpdA;
  std::vector<Double_t> etaValuesEpdA;

  Int_t nhitsEpdB;
  Double_t XnEpdB;
  Double_t YnEpdB;
  Double_t psiEpdB;
  std::vector<Double_t> phiValuesEpdB;
  std::vector<Double_t> etaValuesEpdB;

  Int_t nhitsEpdC;
  Double_t XnEpdC;
  Double_t YnEpdC;
  Double_t psiEpdC;
  std::vector<Double_t> phiValuesEpdC;
  std::vector<Double_t> etaValuesEpdC;

  Int_t nhitsEpdD;
  Double_t XnEpdD;
  Double_t YnEpdD;
  Double_t psiEpdD;
  std::vector<Double_t> phiValuesEpdD;
  std::vector<Double_t> etaValuesEpdD;


  void reset()
  {
    badEvent  = false;  //Reset all values in the struct to reuse
    primTracks = 0;
    centID = BAD_VALUE;

    nTracksTpc = 0;
    XnTpc = 0;
    YnTpc = 0;
    psiTpc = BAD_VALUE;

    nTracksTpcA = 0;
    XnTpcA = 0;
    YnTpcA = 0;
    psiTpcA = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcA);
    std::vector<Double_t>().swap(etaValuesTpcA);

    nTracksTpcB = 0;
    XnTpcB = 0;
    YnTpcB = 0;
    psiTpcB = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcB);
    std::vector<Double_t>().swap(etaValuesTpcB);

    nPionP = 0;
    XnPionP = 0;
    YnPionP = 0;
    psiPionP = BAD_VALUE;

    nPionM = 0;
    XnPionM = 0;
    YnPionM = 0;
    psiPionM = BAD_VALUE;

    nKaonP = 0;
    XnKaonP = 0;
    YnKaonP = 0;
    psiKaonP = BAD_VALUE;

    nKaonM = 0;
    XnKaonM = 0;
    YnKaonM = 0;
    psiKaonM = BAD_VALUE;

    nProton = 0;
    XnProton = 0;
    YnProton = 0;
    psiProton = BAD_VALUE;

    nhitsEpd = 0;
    XnEpd = 0;
    YnEpd = 0;
    psiEpd = BAD_VALUE;

    nhitsEpdA = 0;
    XnEpdA = 0;
    YnEpdA = 0;
    psiEpdA = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdA);
    std::vector<Double_t>().swap(etaValuesEpdA);

    nhitsEpdB = 0;
    XnEpdB = 0;
    YnEpdB = 0;
    psiEpdB = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdB);
    std::vector<Double_t>().swap(etaValuesEpdB);

    nhitsEpdC = 0;
    XnEpdC = 0;
    YnEpdC = 0;
    psiEpdC = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdC);
    std::vector<Double_t>().swap(etaValuesEpdC);

    nhitsEpdD = 0;
    XnEpdD = 0;
    YnEpdD = 0;
    psiEpdD = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdD);
    std::vector<Double_t>().swap(etaValuesEpdD);
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
  else if (event.XnTpcA == 0 && event.YnTpcA == 0) { event.badEvent = true; }
  else if (event.XnTpcB == 0 && event.YnTpcB == 0) { event.badEvent = true; }
  else if (event.XnPionP == 0 && event.YnPionP == 0) { event.badEvent = true; }
  else if (event.XnPionM == 0 && event.YnPionM == 0) { event.badEvent = true; }
  //else if (event.XnKaonP == 0 && event.YnKaonP == 0) { event.badEvent = true; }
  //else if (event.XnKaonM == 0 && event.YnKaonM == 0) { event.badEvent = true; }
  else if (event.XnProton == 0 && event.YnProton == 0) { event.badEvent = true; }
  else if (event.XnEpd == 0 && event.YnEpd == 0) { event.badEvent = true; }
  else if (event.XnEpdA == 0 && event.YnEpdA == 0) { event.badEvent = true; }
  else if (event.XnEpdB == 0 && event.YnEpdB == 0) { event.badEvent = true; }
  else if (event.XnEpdC == 0 && event.YnEpdC == 0) { event.badEvent = true; }
  else if (event.XnEpdD == 0 && event.YnEpdD == 0) { event.badEvent = true; }
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
void fillRawSpect(Double_t px, Double_t py, Double_t pz, Double_t mass, TH1D *dndy, TH1D *dndm)/*, TH2D *MvsY)*/
{
  Double_t y  = rapidity(px, py, pz, mass);
  Double_t mT = transMass(px, py, mass);
  Double_t M  = mT - mass;
  dndy->Fill(y);
  dndm->Fill(M);
  //MvsY->Fill(y,M, 1/(TMath::TwoPi() * mT));
};




void FlowAnalyzer(TString inFile, TString jobID)
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

  //TH1D *h_pTA   = new TH1D("h_pTA", "TPC A Region p_{T};p_{T} (GeV);Particles", 100, 0, 5);
  //TH1D *h_pTB   = new TH1D("h_pTB", "TPC B Region p_{T};p_{T} (GeV);Particles", 100, 0, 5);

  TH1D *h_eta_s   = new TH1D("h_eta_s", "Particle #eta_{CM};#eta-#eta_{mid};Particles", 600, -6, 2);

  TH1D *h_tileWeights = new TH1D("h_tileWeights", "EPD Tile Weights;Hits;nMIP Weights", 5, -1, 4);
  TH1D *h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  //TH1D *h_pi_TPC_nsig = new TH1D("h_pi_TPC_nsig", "#pi TPC n#sigma;n#sigma;Tracks",    200, -10, 10);
  //TH1D *h_ka_TPC_nsig = new TH1D("h_ka_TPC_nsig", "K TPC n#sigma;n#sigma;Tracks",      200, -10, 10);
  //TH1D *h_pr_TPC_nsig = new TH1D("h_pr_TPC_nsig", "Proton TPC n#sigma;n#sigma;Tracks", 200, -10, 10);

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

  TH1D *h_psiTpc_RAW  = new TH1D("h_psiTpc_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcA_RAW = new TH1D("h_psiTpcA_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_RAW = new TH1D("h_psiTpcB_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiPionP_RAW  = new TH1D("h_psiPionP_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", #pi^{+});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiPionM_RAW  = new TH1D("h_psiPionM_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", #pi^{-});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiKaonP_RAW  = new TH1D("h_psiKaonP_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", K^{+});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiKaonM_RAW  = new TH1D("h_psiKaonM_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", K^{-});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiProton_RAW = new TH1D("h_psiProton_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC Protons);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpd_RAW  = new TH1D("h_psiEpd_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdA_RAW = new TH1D("h_psiEpdA_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdB_RAW = new TH1D("h_psiEpdB_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdC_RAW = new TH1D("h_psiEpdC_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD C);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdD_RAW = new TH1D("h_psiEpdD_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD D);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  //TH1D *h_v2Plot = new TH1D("h_v2Plot", "Plot to Retrieve v_{2};cos(2(#phi - #psi_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);

  // Profiles for resolution terms
  TProfile *p_TpcAB = new TProfile("p_TpcAB","TPC A-B Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{TPC,B}_{"+ORDER_N_STR+"}))>",
				   CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  
  TProfile *p_TpcAEpdA = new TProfile("p_TpcAEpdA","TPC A EPD A Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,A}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdB = new TProfile("p_TpcAEpdB","TPC A EPD B Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,B}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdC = new TProfile("p_TpcAEpdC","TPC A EPD C Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,C}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdD = new TProfile("p_TpcAEpdD","TPC A EPD D Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,D}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcBEpdA = new TProfile("p_TpcBEpdA","TPC B EPD A Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,A}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdB = new TProfile("p_TpcBEpdB","TPC B EPD B Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,B}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdC = new TProfile("p_TpcBEpdC","TPC B EPD C Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,C}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdD = new TProfile("p_TpcBEpdD","TPC B EPD D Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,D}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_EpdAEpdB = new TProfile("p_EpdAEpdB","EPD A EPD B Correlations;Centrality;<cos(2(#psi^{EPD,A}_{"+ORDER_N_STR+"}-#psi^{EPD,B}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_EpdAEpdC = new TProfile("p_EpdAEpdC","EPD A EPD C Correlations;Centrality;<cos(2(#psi^{EPD,A}_{"+ORDER_N_STR+"}-#psi^{EPD,C}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_EpdAEpdD = new TProfile("p_EpdAEpdD","EPD A EPD D Correlations;Centrality;<cos(2(#psi^{EPD,A}_{"+ORDER_N_STR+"}-#psi^{EPD,D}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_EpdBEpdC = new TProfile("p_EpdBEpdC","EPD B EPD C Correlations;Centrality;<cos(2(#psi^{EPD,B}_{"+ORDER_N_STR+"}-#psi^{EPD,C}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_EpdBEpdD = new TProfile("p_EpdBEpdD","EPD B EPD D Correlations;Centrality;<cos(2(#psi^{EPD,B}_{"+ORDER_N_STR+"}-#psi^{EPD,D}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_EpdCEpdD = new TProfile("p_EpdCEpdD","EPD C EPD D Correlations;Centrality;<cos(2(#psi^{EPD,C}_{"+ORDER_N_STR+"}-#psi^{EPD,D}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  //TH1D *h_abCorr_n2 = new TH1D("h_abCorr_n2", "2nd Order Subevent a-b Correlations;cos(2(#psi^{a}_{"+ORDER_N_STR+"} - #psi^{b}_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);
  //TH1D *h_acCorr_n2 = new TH1D("h_acCorr_n2", "2nd Order Subevent a-c Correlations;cos(2(#psi^{a}_{"+ORDER_N_STR+"} - #psi^{c}_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);
  //TH1D *h_bcCorr_n2 = new TH1D("h_bcCorr_n2", "2nd Order Subevent b-c Correlations;cos(2(#psi^{b}_{"+ORDER_N_STR+"} - #psi^{c}_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);


  TH2D *h2_pp_vs_eta = new TH2D("h2_pp_vs_eta","Tile Weight for Supersectors vs #eta;#eta;Supersector", 400, -6, -2, 12, 0.5, 12.5);
  TH2D *h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  TH2D *h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);

  TH2D *h2_betap = new TH2D("h2_betap","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_m2_p   = new TH2D("h2_m2_p", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 500, -3, 3, 500, -0.1, 15);
  TH2D *h2_m2_vs_qpT  = new TH2D("h2_m2_vs_qpT", "m^{2} vs q*p_{T};q*p_{T} (GeV);m^{2} (GeV^{2})", 500, -3, 3, 500, -0.1, 1.2);
  TH2D *h2_dEdx_vs_qp = new TH2D("h2_dEdx_vs_qp", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  TH2D *h2_nSig_vs_qp_pi = new TH2D("h2_nSig_vs_qp_pi", "Pion n#sigma vs q|p|;q|p| (GeV); n#sigma_{#pi}", 400, -2, 2, 400, -8, 8);
  TH2D *h2_nSig_vs_qp_ka = new TH2D("h2_nSig_vs_qp_ka", "Kaon n#sigma vs q|p|;q|p| (GeV); n#sigma_{K}", 400, -2, 2, 400, -8, 8);
  TH2D *h2_nSig_vs_qp_pr = new TH2D("h2_nSig_vs_qp_pr", "Proton n#sigma vs q|p|;q|p| (GeV); n#sigma_{p}", 400, -2, 2, 400, -8, 8);
  TH2D *h2_pi_m2_vs_TPC_nsig = new TH2D("h2_pi_m2_vs_TPC_nsig", "m^{2} vs #pi TPC n#sigma;n#sigma_{#pi};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);
  TH2D *h2_ka_m2_vs_TPC_nsig = new TH2D("h2_ka_m2_vs_TPC_nsig", "m^{2} vs K TPC n#sigma;n#sigma_{K};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);
  TH2D *h2_pr_m2_vs_TPC_nsig = new TH2D("h2_pr_m2_vs_TPC_nsig", "m^{2} vs Proton TPC n#sigma;n#sigma_{pro};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);

  TH2D *h2_y_vs_eta = new TH2D("h2_y_vs_eta", "TPC A All Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_pp = new TH2D("h2_y_vs_eta_pp", "TPC A #pi^{+} y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_pp = new TH2D("h2_y_vs_eta_pt0p5to1_pp", "#pi^{+} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pp = new TH2D("h2_y_vs_eta_pt1to1p5_pp", "#pi^{+} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pp = new TH2D("h2_y_vs_eta_pt1p5to2_pp", "#pi^{+} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_pm = new TH2D("h2_y_vs_eta_pm", "TPC A #pi^{-} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_pm = new TH2D("h2_y_vs_eta_pt0p5to1_pm", "#pi^{-} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pm = new TH2D("h2_y_vs_eta_pt1to1p5_pm", "#pi^{-} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pm = new TH2D("h2_y_vs_eta_pt1p5to2_pm", "#pi^{-} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_kp = new TH2D("h2_y_vs_eta_kp", "TPC A K^{+} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_kp = new TH2D("h2_y_vs_eta_pt0p5to1_kp", "K^{+} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_kp = new TH2D("h2_y_vs_eta_pt1to1p5_kp", "K^{+} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_kp = new TH2D("h2_y_vs_eta_pt1p5to2_kp", "K^{+} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_km = new TH2D("h2_y_vs_eta_km", "TPC A K^{-} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt0p5to1_km = new TH2D("h2_y_vs_eta_pt0p5to1_km", "K^{-} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_km = new TH2D("h2_y_vs_eta_pt1to1p5_km", "K^{-} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_km = new TH2D("h2_y_vs_eta_pt1p5to2_km", "K^{-} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_pr = new TH2D("h2_y_vs_eta_pr", "TPC A Proton Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
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


  TH2D *h2_psiEpdATpcA = new TH2D("h2_psiEpdATpcA", "#psi^{EPD,A} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcA = new TH2D("h2_psiEpdBTpcA", "#psi^{EPD,B} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdCTpcA = new TH2D("h2_psiEpdCTpcA", "#psi^{EPD,C} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{C}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdDTpcA = new TH2D("h2_psiEpdDTpcA", "#psi^{EPD,D} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{D}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiEpdATpcB = new TH2D("h2_psiEpdATpcB", "#psi^{EPD,A} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcB = new TH2D("h2_psiEpdBTpcB", "#psi^{EPD,B} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdCTpcB = new TH2D("h2_psiEpdCTpcB", "#psi^{EPD,C} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{C}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdDTpcB = new TH2D("h2_psiEpdDTpcB", "#psi^{EPD,D} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{D}", 
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
  TProfile *p_sinAvgsPionP = new TProfile("p_sinAvgsPionP", "Sin Averages (#pi^{+});j (Correction term);<sin(jn#psi^{#pi^{+}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsPionP = new TProfile("p_cosAvgsPionP", "Cos Averages (#pi^{+});j (Correction term);<sin(jn#psi^{#pi^{+}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsPionM = new TProfile("p_sinAvgsPionM", "Sin Averages (#pi^{-});j (Correction term);<sin(jn#psi^{#pi^{-}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsPionM = new TProfile("p_cosAvgsPionM", "Cos Averages (#pi^{-});j (Correction term);<sin(jn#psi^{#pi^{-}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  //TProfile *p_sinAvgsKaonP = new TProfile("p_sinAvgsKaonP", "Sin Averages (K^{+});j (Correction term);<sin(jn#psi^{K^{+}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  //TProfile *p_cosAvgsKaonP = new TProfile("p_cosAvgsKaonP", "Cos Averages (K^{+});j (Correction term);<sin(jn#psi^{K^{+}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  //TProfile *p_sinAvgsKaonM = new TProfile("p_sinAvgsKaonM", "Sin Averages (K^{-});j (Correction term);<sin(jn#psi^{K^{-}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  //TProfile *p_cosAvgsKaonM = new TProfile("p_cosAvgsKaonM", "Cos Averages (K^{-});j (Correction term);<sin(jn#psi^{K^{-}}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsProton = new TProfile("p_sinAvgsProton", "Sin Averages (Protons);j (Correction term);<sin(jn#psi^{pro}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsProton = new TProfile("p_cosAvgsProton", "Cos Averages (Protons);j (Correction term);<sin(jn#psi^{pro}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpd  = new TProfile("p_sinAvgsEpd", "Sin Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpd  = new TProfile("p_cosAvgsEpd", "Cos Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdA = new TProfile("p_sinAvgsEpdA", "Sin Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD,A}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdA = new TProfile("p_cosAvgsEpdA", "Cos Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD,A}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdB = new TProfile("p_sinAvgsEpdB", "Sin Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD,B}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdB = new TProfile("p_cosAvgsEpdB", "Cos Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD,B}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdC = new TProfile("p_sinAvgsEpdC", "Sin Averages (EPD C);j (Correction term);<sin(jn#psi^{EPD,C}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdC = new TProfile("p_cosAvgsEpdC", "Cos Averages (EPD C);j (Correction term);<sin(jn#psi^{EPD,C}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdD = new TProfile("p_sinAvgsEpdD", "Sin Averages (EPD D);j (Correction term);<sin(jn#psi^{EPD,D}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdD = new TProfile("p_cosAvgsEpdD", "Dos Averages (EPD D);j (Correction term);<sin(jn#psi^{EPD,D}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);

  TH1D *h_XnTpc  = new TH1D("h_XnTpc", "X_n Distribution (TPC);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpc  = new TH1D("h_YnTpc", "Y_n Distribution (TPC);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcA = new TH1D("h_XnTpcA", "X_n Distribution (TPC A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcA = new TH1D("h_YnTpcA", "Y_n Distribution (TPC A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcB = new TH1D("h_XnTpcB", "X_n Distribution (TPC B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcB = new TH1D("h_YnTpcB", "Y_n Distribution (TPC B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnPionP = new TH1D("h_XnPionP", "X_n Distribution (#pi^{+});X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnPionP = new TH1D("h_YnPionP", "Y_n Distribution (#pi^{+});Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnPionM = new TH1D("h_XnPionM", "X_n Distribution (#pi^{-});X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnPionM = new TH1D("h_YnPionM", "Y_n Distribution (#pi^{-});Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnKaonP = new TH1D("h_XnKaonP", "X_n Distribution (K^{+});X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnKaonP = new TH1D("h_YnKaonP", "Y_n Distribution (K^{+});Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnKaonM = new TH1D("h_XnKaonM", "X_n Distribution (K^{-});X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnKaonM = new TH1D("h_YnKaonM", "Y_n Distribution (K^{-});Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnProton = new TH1D("h_XnProton", "X_n Distribution (Protons);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnProton = new TH1D("h_YnProton", "Y_n Distribution (Protons);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpd  = new TH1D("h_XnEpd", "X_n Distribution (EPD);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpd  = new TH1D("h_YnEpd", "Y_n Distribution (EPD);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdA = new TH1D("h_XnEpdA", "X_n Distribution (EPD A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdA = new TH1D("h_YnEpdA", "Y_n Distribution (EPD A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdB = new TH1D("h_XnEpdB", "X_n Distribution (EPD B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdB = new TH1D("h_YnEpdB", "Y_n Distribution (EPD B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdC = new TH1D("h_XnEpdC", "X_n Distribution (EPD C);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdC = new TH1D("h_YnEpdC", "Y_n Distribution (EPD C);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdD = new TH1D("h_XnEpdD", "X_n Distribution (EPD D);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdD = new TH1D("h_YnEpdD", "Y_n Distribution (EPD D);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);

  
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

	  StPicoBTofPidTraits *trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());

	  Double_t d_tofBeta = trait->btofBeta();

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
	      //h_eta_s->Fill(d_eta - Y_MID);

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


	      if (d_eta > MIN_TPC_ETA_CUT && d_eta < AGAP_TPC_ETA_CUT)          // TPC A
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

	      //h_pi_TPC_nsig->Fill(d_TPCnSigmaPion);
	      //h_ka_TPC_nsig->Fill(d_TPCnSigmaProton);
	      //h_pr_TPC_nsig->Fill(d_TPCnSigmaKaon);

	      h2_nSig_vs_qp_pi->Fill(d_charge * d_mom, d_TPCnSigmaPion);
	      h2_nSig_vs_qp_ka->Fill(d_charge * d_mom, d_TPCnSigmaKaon);
	      h2_nSig_vs_qp_pr->Fill(d_charge * d_mom, d_TPCnSigmaProton);

	      h2_pi_m2_vs_TPC_nsig->Fill(d_TPCnSigmaPion, d_m2);
	      h2_ka_m2_vs_TPC_nsig->Fill(d_TPCnSigmaKaon, d_m2);
	      h2_pr_m2_vs_TPC_nsig->Fill(d_TPCnSigmaProton, d_m2);

	      h2_m2_vs_qpT->Fill(d_charge * d_pT, d_m2);
	      h2_dEdx_vs_qp->Fill(d_charge * d_mom, d_dEdx);

	      //=========================================================
	      //          PID Cuts
	      //=========================================================

	      Bool_t pion   = ( (d_TPCnSigmaPion > -2)   && (d_TPCnSigmaPion < 2)   && (d_m2 > 0.0) && (d_m2 < 0.05) );    // PARTICLE TAGGING
	      Bool_t kaon   = ( (d_TPCnSigmaKaon > -2)   && (d_TPCnSigmaKaon < 2)   && (d_m2 > 0.2) && (d_m2 < 0.3)  );
	      Bool_t proton = ( (d_TPCnSigmaProton > -2) && (d_TPCnSigmaProton < 2) && (d_m2 > 0.8) && (d_m2 < 1.0)  );

	      if (!pion && !kaon && !proton) continue;

	      if (pion && kaon)   continue;   // Particle must be exclusively one particular type.
	      if (pion && proton) continue;
	      if (kaon && proton) continue;

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
		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h2_y_vs_eta_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h_pp_pT->Fill(d_pT);
		      h_pp_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_pi, h_pp_dndy, h_pp_dndm);

		      eventInfo.nPionP++;
		      if (d_eta > Y_MID_PI_P)        // Account for Q vector sign change past mid-rapidity.
			{
			  eventInfo.XnPionP += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
			  eventInfo.YnPionP += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
			}
		      else if (d_eta < Y_MID_PI_P)
			{
			  eventInfo.XnPionP -= d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
			  eventInfo.YnPionP -= d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
			}

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		    }
		  else if (d_charge < 0) 
		    {
		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h2_y_vs_eta_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi));
		      h_pm_pT->Fill(d_pT);
		      h_pm_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_pi, h_pm_dndy, h_pm_dndm);

		      eventInfo.nPionM++;
		      if (d_eta > Y_MID_PI_M)        // Account for Q vector sign change past mid-rapidity.
			{
			  eventInfo.XnPionM += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
			  eventInfo.YnPionM += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
			}
		      else if (d_eta < Y_MID_PI_M)
			{
			  eventInfo.XnPionM -= d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
			  eventInfo.YnPionM -= d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
			}


		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		    }
		}
	      else if (kaon)
		{
		  if (d_charge > 0) 
		    {
		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h2_y_vs_eta_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h_kp_pT->Fill(d_pT);
		      h_kp_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_ka, h_kp_dndy, h_kp_dndm);

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		    }
		  else if (d_charge < 0)		 
		    {
		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h2_y_vs_eta_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka));
		      h_km_pT->Fill(d_pT);
		      h_km_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_ka, h_km_dndy, h_km_dndm);

		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		    }
		}
	      else if (proton)
		{
		  if (d_charge > 0) 
		    {
		      h2_y_vs_eta->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr));
		      h2_y_vs_eta_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr));
		      h_pr_pT->Fill(d_pT);
		      h_pr_mom->Fill(d_mom);
		      fillRawSpect(d_px, d_py, d_pz, d_m0_pr, h_pr_dndy, h_pr_dndm);

		      eventInfo.nProton++;
		      if (d_eta > Y_MID_PR)        // Account for Q vector sign change past mid-rapidity.
			{
			  eventInfo.XnProton += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
			  eventInfo.YnProton += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
			}
		      else if (d_eta < Y_MID_PR)
			{
			  eventInfo.XnProton -= d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
			  eventInfo.YnProton -= d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
			}

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

      if (eventInfo.centID == BAD_VALUE) continue;
      
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

	  /*
	  if (tileEta > A_MIN_ETA_CUT && tileEta < A_MAX_ETA_CUT)
	    {
	    }
	  else if (tileEta > B_MIN_ETA_CUT && tileEta < B_MAX_ETA_CUT)
	    {
	    }
	  else if (tileEta > C_MIN_ETA_CUT && tileEta < C_MAX_ETA_CUT)
	    {
	    }
	  */

	  eventInfo.nhitsEpd++;
	  eventInfo.XnEpd += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	  eventInfo.YnEpd += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);

	  if (tileEta > MIN_ETA_CUT && tileEta < AB_ETA_CUT)  // Sub A
	    {
	      eventInfo.nhitsEpdA++;
	      eventInfo.XnEpdA += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdA += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdA.push_back(tilePhi);
	      eventInfo.etaValuesEpdA.push_back(tileEta);
	    }
	  else if (tileEta >= AB_ETA_CUT && tileEta < BC_ETA_CUT)  // Sub B
	    {
	      eventInfo.nhitsEpdB++;
	      eventInfo.XnEpdB += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdB += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdB.push_back(tilePhi);
	      eventInfo.etaValuesEpdB.push_back(tileEta);
	    }
	  else if (tileEta >= BC_ETA_CUT && tileEta < CD_ETA_CUT)  // Sub C
	    {
	      eventInfo.nhitsEpdC++;
	      eventInfo.XnEpdC += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdC += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdC.push_back(tilePhi);
	      eventInfo.etaValuesEpdC.push_back(tileEta);
	    }
	  else if (tileEta >= CD_ETA_CUT && tileEta < MAX_ETA_CUT)  // Sub D
	    {
	      eventInfo.nhitsEpdD++;
	      eventInfo.XnEpdD += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdD += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdD.push_back(tilePhi);
	      eventInfo.etaValuesEpdD.push_back(tileEta);
	    }

	}
      delete epdGeom;

      //=========================================================
      //            END EPD STUFF
      //=========================================================


      if (eventInfo.nTracksTpc  < MIN_TRACKS) continue;
      if (eventInfo.nTracksTpcA < MIN_TRACKS) continue;
      if (eventInfo.nTracksTpcB < MIN_TRACKS) continue;
      if (eventInfo.nPionP      < MIN_TRACKS) continue;
      if (eventInfo.nPionM      < MIN_TRACKS) continue;
      if (eventInfo.nProton     < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpd    < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdA   < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdB   < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdC   < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdD   < MIN_TRACKS) continue;
      
      checkZeroQ(eventInfo);

      if (eventInfo.badEvent) continue;

      // RAW SUB-EVENT PLANE ANGLES //
      if (ORDER_N % 2 == 1)           // Q vectors must change sign past mid-rapidity; I think this is for 1st order event-planes only. Full TPC already takes this into account.
	{
	  eventInfo.XnTpcA *= -1.0;
	  eventInfo.YnTpcA *= -1.0;
	  eventInfo.XnEpd  *= -1.0;
	  eventInfo.YnEpd  *= -1.0;
	  eventInfo.XnEpdA *= -1.0;
	  eventInfo.YnEpdA *= -1.0;
	  eventInfo.XnEpdB *= -1.0;
	  eventInfo.YnEpdB *= -1.0;
	  eventInfo.XnEpdC *= -1.0;
	  eventInfo.YnEpdC *= -1.0;
	}
      eventInfo.psiTpc  = TMath::ATan2(eventInfo.YnTpc,  eventInfo.XnTpc)  / (Double_t)ORDER_N;
      eventInfo.psiTpcA = TMath::ATan2(eventInfo.YnTpcA, eventInfo.XnTpcA) / (Double_t)ORDER_N;
      eventInfo.psiTpcB = TMath::ATan2(eventInfo.YnTpcB, eventInfo.XnTpcB) / (Double_t)ORDER_N;
      eventInfo.psiPionP  = TMath::ATan2(eventInfo.YnPionP,  eventInfo.XnPionP)  / (Double_t)ORDER_N;
      eventInfo.psiPionM  = TMath::ATan2(eventInfo.YnPionM,  eventInfo.XnPionM)  / (Double_t)ORDER_N;
      eventInfo.psiProton = TMath::ATan2(eventInfo.YnProton, eventInfo.XnProton) / (Double_t)ORDER_N;
      eventInfo.psiEpd  = TMath::ATan2(eventInfo.YnEpd,  eventInfo.XnEpd)  / (Double_t)ORDER_N;
      eventInfo.psiEpdA = TMath::ATan2(eventInfo.YnEpdA, eventInfo.XnEpdA) / (Double_t)ORDER_N;
      eventInfo.psiEpdB = TMath::ATan2(eventInfo.YnEpdB, eventInfo.XnEpdB) / (Double_t)ORDER_N;
      eventInfo.psiEpdC = TMath::ATan2(eventInfo.YnEpdC, eventInfo.XnEpdC) / (Double_t)ORDER_N;
      eventInfo.psiEpdD = TMath::ATan2(eventInfo.YnEpdD, eventInfo.XnEpdD) / (Double_t)ORDER_N;


      // Fill eta distribution here since this is past all possible cuts.

      for (unsigned int i = 0; i < eventInfo.etaValuesTpcA.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcA.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcB.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcB.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdA.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdA.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdB.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdB.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdC.at(i) - Y_MID ); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdC.at(i) - Y_MID ); } 


      h_primTracks->Fill(eventInfo.primTracks);
      h_centralities->Fill(eventInfo.centID);

      h_XnTpc->Fill(eventInfo.XnTpc);
      h_YnTpc->Fill(eventInfo.YnTpc);
      h_XnTpcA->Fill(eventInfo.XnTpcA);
      h_YnTpcA->Fill(eventInfo.YnTpcA);
      h_XnTpcB->Fill(eventInfo.XnTpcB);
      h_YnTpcB->Fill(eventInfo.YnTpcB);
      h_XnPionP->Fill(eventInfo.XnPionP);
      h_YnPionP->Fill(eventInfo.YnPionP);
      h_XnPionM->Fill(eventInfo.XnPionM);
      h_YnPionM->Fill(eventInfo.YnPionM);
      h_XnProton->Fill(eventInfo.XnProton);
      h_YnProton->Fill(eventInfo.YnProton);
      h_XnEpd->Fill(eventInfo.XnEpd);
      h_YnEpd->Fill(eventInfo.YnEpd);
      h_XnEpdA->Fill(eventInfo.XnEpdA);
      h_YnEpdA->Fill(eventInfo.YnEpdA);
      h_XnEpdB->Fill(eventInfo.XnEpdB);
      h_YnEpdB->Fill(eventInfo.YnEpdB);
      h_XnEpdC->Fill(eventInfo.XnEpdC);
      h_YnEpdC->Fill(eventInfo.YnEpdC);
      h_XnEpdD->Fill(eventInfo.XnEpdD);
      h_YnEpdD->Fill(eventInfo.YnEpdD);

      h_psiTpc_RAW->Fill(eventInfo.psiTpc);
      h_psiTpcA_RAW->Fill(eventInfo.psiTpcA);
      h_psiTpcB_RAW->Fill(eventInfo.psiTpcB);
      h_psiEpd_RAW->Fill(eventInfo.psiEpd);
      h_psiPionP_RAW->Fill(eventInfo.psiPionP);
      h_psiPionM_RAW->Fill(eventInfo.psiPionM);
      h_psiProton_RAW->Fill(eventInfo.psiProton);
      h_psiEpdA_RAW->Fill(eventInfo.psiEpdA);
      h_psiEpdB_RAW->Fill(eventInfo.psiEpdB);
      h_psiEpdC_RAW->Fill(eventInfo.psiEpdC);
      h_psiEpdD_RAW->Fill(eventInfo.psiEpdD);

      v_events.push_back(eventInfo);   // Store this event with all of its attributes
    }//End event loop

  eventInfo.reset();

  TH1D *h_XnTpc_RC;   // Re-centered histograms
  TH1D *h_YnTpc_RC;
  TH1D *h_XnTpcA_RC;
  TH1D *h_YnTpcA_RC;
  TH1D *h_XnTpcB_RC;
  TH1D *h_YnTpcB_RC;
  TH1D *h_XnPionP_RC;
  TH1D *h_YnPionP_RC;
  TH1D *h_XnPionM_RC;
  TH1D *h_YnPionM_RC;
  TH1D *h_XnProton_RC;
  TH1D *h_YnProton_RC;
  TH1D *h_XnEpd_RC;
  TH1D *h_YnEpd_RC;
  TH1D *h_XnEpdA_RC;
  TH1D *h_YnEpdA_RC;
  TH1D *h_XnEpdB_RC;
  TH1D *h_YnEpdB_RC;
  TH1D *h_XnEpdC_RC;
  TH1D *h_YnEpdC_RC;
  TH1D *h_XnEpdD_RC;
  TH1D *h_YnEpdD_RC;

  TH1D *h_psiTpc_RC;
  TH1D *h_psiTpcA_RC;
  TH1D *h_psiTpcB_RC;
  TH1D *h_psiPionP_RC;
  TH1D *h_psiPionM_RC;
  TH1D *h_psiProton_RC;
  TH1D *h_psiEpd_RC;
  TH1D *h_psiEpdA_RC;
  TH1D *h_psiEpdB_RC;
  TH1D *h_psiEpdC_RC;
  TH1D *h_psiEpdD_RC;

  TH1D *h_psiTpc_FLAT;
  TH1D *h_psiTpcA_FLAT;
  TH1D *h_psiTpcB_FLAT;
  TH1D *h_psiPionP_FLAT;
  TH1D *h_psiPionM_FLAT;
  TH1D *h_psiProton_FLAT;
  TH1D *h_psiEpd_FLAT;
  TH1D *h_psiEpdA_FLAT;
  TH1D *h_psiEpdB_FLAT;
  TH1D *h_psiEpdC_FLAT;
  TH1D *h_psiEpdD_FLAT;


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
      h_XnPionP_RC  = new TH1D("h_XnPionP_RC", "Re-centered X_n Distribution (#pi^{+});X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnPionP_RC  = new TH1D("h_YnPionP_RC", "Re-centered Y_n Distribution (#pi^{+});Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnPionM_RC  = new TH1D("h_XnPionM_RC", "Re-centered X_n Distribution (#pi^{-});X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnPionM_RC  = new TH1D("h_YnPionM_RC", "Re-centered Y_n Distribution (#pi^{-});Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnProton_RC = new TH1D("h_XnProton_RC", "Re-centered X_n Distribution (Protons);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnProton_RC = new TH1D("h_YnProton_RC", "Re-centered Y_n Distribution (Protons);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpd_RC  = new TH1D("h_XnEpd_RC", "Re-centered X_n Distribution (EPD);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpd_RC  = new TH1D("h_YnEpd_RC", "Re-centered Y_n Distribution (EPD);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdA_RC = new TH1D("h_XnEpdA_RC", "Re-centered X_n Distribution (EPD A);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdA_RC = new TH1D("h_YnEpdA_RC", "Re-centered Y_n Distribution (EPD A);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdB_RC = new TH1D("h_XnEpdB_RC", "Re-centered X_n Distribution (EPD B);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdB_RC = new TH1D("h_YnEpdB_RC", "Re-centered Y_n Distribution (EPD B);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdC_RC = new TH1D("h_XnEpdC_RC", "Re-centered X_n Distribution (EPD C);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdC_RC = new TH1D("h_YnEpdC_RC", "Re-centered Y_n Distribution (EPD C);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdD_RC = new TH1D("h_XnEpdD_RC", "Re-centered X_n Distribution (EPD D);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdD_RC = new TH1D("h_YnEpdD_RC", "Re-centered Y_n Distribution (EPD D);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);

      h_psiTpc_RC  = new TH1D("h_psiTpc_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiTpcA_RC = new TH1D("h_psiTpcA_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiTpcB_RC = new TH1D("h_psiTpcB_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiPionP_RC  = new TH1D("h_psiPionP_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", #pi^{+});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiPionM_RC  = new TH1D("h_psiPionM_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", #pi^{-});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiProton_RC = new TH1D("h_psiProton_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", Protons);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpd_RC  = new TH1D("h_psiEpd_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdA_RC = new TH1D("h_psiEpdA_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdB_RC = new TH1D("h_psiEpdB_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC_RC = new TH1D("h_psiEpdC_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD C);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdD_RC = new TH1D("h_psiEpdD_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD D);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

      TH1D *h_XnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_XnTpc");
      TH1D *h_XnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcA");
      TH1D *h_XnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcB");
      TH1D *h_XnPionP_INPUT = (TH1D*)correctionInputFile->Get("h_XnPionP");
      TH1D *h_XnPionM_INPUT = (TH1D*)correctionInputFile->Get("h_XnPionM");
      TH1D *h_XnProton_INPUT = (TH1D*)correctionInputFile->Get("h_XnProton");
      TH1D *h_XnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_XnEpd");
      TH1D *h_XnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdA");
      TH1D *h_XnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdB");
      TH1D *h_XnEpdC_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdC");
      TH1D *h_XnEpdD_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdD");

      TH1D *h_YnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_YnTpc");
      TH1D *h_YnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcA");
      TH1D *h_YnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcB");
      TH1D *h_YnPionP_INPUT = (TH1D*)correctionInputFile->Get("h_YnPionP");
      TH1D *h_YnPionM_INPUT = (TH1D*)correctionInputFile->Get("h_YnPionM");
      TH1D *h_YnProton_INPUT = (TH1D*)correctionInputFile->Get("h_YnProton");
      TH1D *h_YnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_YnEpd");
      TH1D *h_YnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdA");
      TH1D *h_YnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdB");
      TH1D *h_YnEpdC_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdC");
      TH1D *h_YnEpdD_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdD");

      Double_t d_XnTpc_Avg  = h_XnTpc_INPUT->GetMean();
      Double_t d_XnTpcA_Avg = h_XnTpcA_INPUT->GetMean();
      Double_t d_XnTpcB_Avg = h_XnTpcB_INPUT->GetMean();
      Double_t d_XnPionP_Avg  = h_XnPionP_INPUT->GetMean();
      Double_t d_XnPionM_Avg  = h_XnPionM_INPUT->GetMean();
      Double_t d_XnProton_Avg = h_XnProton_INPUT->GetMean();
      Double_t d_XnEpd_Avg  = h_XnEpd_INPUT->GetMean();
      Double_t d_XnEpdA_Avg = h_XnEpdA_INPUT->GetMean();
      Double_t d_XnEpdB_Avg = h_XnEpdB_INPUT->GetMean();
      Double_t d_XnEpdC_Avg = h_XnEpdC_INPUT->GetMean();
      Double_t d_XnEpdD_Avg = h_XnEpdD_INPUT->GetMean();

      Double_t d_YnTpc_Avg  = h_YnTpc_INPUT->GetMean();
      Double_t d_YnTpcA_Avg = h_YnTpcA_INPUT->GetMean();
      Double_t d_YnTpcB_Avg = h_YnTpcB_INPUT->GetMean();
      Double_t d_YnPionP_Avg  = h_YnPionP_INPUT->GetMean();
      Double_t d_YnPionM_Avg  = h_YnPionM_INPUT->GetMean();
      Double_t d_YnProton_Avg = h_YnProton_INPUT->GetMean();
      Double_t d_YnEpd_Avg  = h_YnEpd_INPUT->GetMean();
      Double_t d_YnEpdA_Avg = h_YnEpdA_INPUT->GetMean();
      Double_t d_YnEpdB_Avg = h_YnEpdB_INPUT->GetMean();
      Double_t d_YnEpdC_Avg = h_YnEpdC_INPUT->GetMean();
      Double_t d_YnEpdD_Avg = h_YnEpdD_INPUT->GetMean();


      Int_t numOfEvents = v_events.size();
      Int_t badEvents = 0;

      for (int i = 0; i < numOfEvents; i++)
	{
	  v_events.at(i).XnTpc  -= d_XnTpc_Avg;
	  v_events.at(i).XnTpcA -= d_XnTpcA_Avg;
	  v_events.at(i).XnTpcB -= d_XnTpcB_Avg;
	  v_events.at(i).XnPionP  -= d_XnPionP_Avg;
	  v_events.at(i).XnPionM  -= d_XnPionM_Avg;
	  v_events.at(i).XnProton -= d_XnProton_Avg;
	  v_events.at(i).XnEpd  -= d_XnEpd_Avg;
	  v_events.at(i).XnEpdA -= d_XnEpdA_Avg;
	  v_events.at(i).XnEpdB -= d_XnEpdB_Avg;
	  v_events.at(i).XnEpdC -= d_XnEpdC_Avg;
	  v_events.at(i).XnEpdD -= d_XnEpdD_Avg;

	  v_events.at(i).YnTpc  -= d_YnTpc_Avg;
	  v_events.at(i).YnTpcA -= d_YnTpcA_Avg;
	  v_events.at(i).YnTpcB -= d_YnTpcB_Avg;
	  v_events.at(i).YnPionP  -= d_YnPionP_Avg;
	  v_events.at(i).YnPionM  -= d_YnPionM_Avg;
	  v_events.at(i).YnProton -= d_YnProton_Avg;
	  v_events.at(i).YnEpd  -= d_YnEpd_Avg;
	  v_events.at(i).YnEpdA -= d_YnEpdA_Avg;
	  v_events.at(i).YnEpdB -= d_YnEpdB_Avg;
	  v_events.at(i).YnEpdC -= d_YnEpdC_Avg;
	  v_events.at(i).YnEpdD -= d_YnEpdD_Avg;

	  checkZeroQ(v_events.at(i));

	  if ( v_events.at(i).badEvent ) { badEvents++; continue; }

	  h_XnTpc_RC->Fill(v_events.at(i).XnTpc);
	  h_XnTpcA_RC->Fill(v_events.at(i).XnTpcA);
	  h_XnTpcB_RC->Fill(v_events.at(i).XnTpcB);
	  h_XnPionP_RC->Fill(v_events.at(i).XnPionP);
	  h_XnPionM_RC->Fill(v_events.at(i).XnPionM);
	  h_XnProton_RC->Fill(v_events.at(i).XnProton);
	  h_XnEpd_RC->Fill(v_events.at(i).XnEpd);
	  h_XnEpdA_RC->Fill(v_events.at(i).XnEpdA);
	  h_XnEpdB_RC->Fill(v_events.at(i).XnEpdB);
	  h_XnEpdC_RC->Fill(v_events.at(i).XnEpdC);
	  h_XnEpdD_RC->Fill(v_events.at(i).XnEpdD);

	  h_YnTpc_RC->Fill(v_events.at(i).YnTpc);
	  h_YnTpcA_RC->Fill(v_events.at(i).YnTpcA);
	  h_YnTpcB_RC->Fill(v_events.at(i).YnTpcB);
	  h_YnPionP_RC->Fill(v_events.at(i).YnPionP);
	  h_YnPionM_RC->Fill(v_events.at(i).YnPionM);
	  h_YnProton_RC->Fill(v_events.at(i).YnProton);
	  h_YnEpd_RC->Fill(v_events.at(i).YnEpd);
	  h_YnEpdA_RC->Fill(v_events.at(i).YnEpdA);
	  h_YnEpdB_RC->Fill(v_events.at(i).YnEpdB);
	  h_YnEpdC_RC->Fill(v_events.at(i).YnEpdC);
	  h_YnEpdD_RC->Fill(v_events.at(i).YnEpdD);

	  // Recalculate the event plane angles after re-centering	  
	  v_events.at(i).psiTpc  = TMath::ATan2(v_events.at(i).YnTpc,  v_events.at(i).XnTpc)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiTpcA = TMath::ATan2(v_events.at(i).YnTpcA, v_events.at(i).XnTpcA) / (Double_t)ORDER_N; 
	  v_events.at(i).psiTpcB = TMath::ATan2(v_events.at(i).YnTpcB, v_events.at(i).XnTpcB) / (Double_t)ORDER_N; 
	  v_events.at(i).psiPionP  = TMath::ATan2(v_events.at(i).YnPionP, v_events.at(i).XnPionP)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiPionM  = TMath::ATan2(v_events.at(i).YnPionM, v_events.at(i).XnPionM)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiProton = TMath::ATan2(v_events.at(i).YnProton, v_events.at(i).XnProton)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpd  = TMath::ATan2(v_events.at(i).YnEpd,  v_events.at(i).XnEpd)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdA = TMath::ATan2(v_events.at(i).YnEpdA, v_events.at(i).XnEpdA) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdB = TMath::ATan2(v_events.at(i).YnEpdB, v_events.at(i).XnEpdB) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdC = TMath::ATan2(v_events.at(i).YnEpdC, v_events.at(i).XnEpdC) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdD = TMath::ATan2(v_events.at(i).YnEpdD, v_events.at(i).XnEpdD) / (Double_t)ORDER_N; 

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_N);
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiPionP  = angleShift(v_events.at(i).psiPionP, ORDER_N);
	  v_events.at(i).psiPionM  = angleShift(v_events.at(i).psiPionM, ORDER_N);
	  v_events.at(i).psiProton = angleShift(v_events.at(i).psiProton, ORDER_N);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_N);
	  v_events.at(i).psiEpdA = angleShift(v_events.at(i).psiEpdA, ORDER_N);
	  v_events.at(i).psiEpdB = angleShift(v_events.at(i).psiEpdB, ORDER_N);
	  v_events.at(i).psiEpdC = angleShift(v_events.at(i).psiEpdC, ORDER_N);
	  v_events.at(i).psiEpdD = angleShift(v_events.at(i).psiEpdD, ORDER_N);

	  h_psiTpc_RC->Fill(v_events.at(i).psiTpc);
	  h_psiTpcA_RC->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_RC->Fill(v_events.at(i).psiTpcB);
	  h_psiPionP_RC->Fill(v_events.at(i).psiPionP);
	  h_psiPionM_RC->Fill(v_events.at(i).psiPionM);
	  h_psiProton_RC->Fill(v_events.at(i).psiProton);
	  h_psiEpd_RC->Fill(v_events.at(i).psiEpd);
	  h_psiEpdA_RC->Fill(v_events.at(i).psiEpdA);
	  h_psiEpdB_RC->Fill(v_events.at(i).psiEpdB);
	  h_psiEpdC_RC->Fill(v_events.at(i).psiEpdC);
	  h_psiEpdD_RC->Fill(v_events.at(i).psiEpdD);


	  // Accumulate terms for averages over the re-centered angles for event plane angle shifting
	  for (int j = 1; j <= SHIFT_TERMS; j++)
	    {
	      p_sinAvgsTpc->Fill(j,  TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc));
	      p_cosAvgsTpc->Fill(j,  TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc));
	      p_sinAvgsTpcA->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      p_cosAvgsTpcA->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      p_sinAvgsTpcB->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      p_cosAvgsTpcB->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      p_sinAvgsPionP->Fill(j,  TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionP));
	      p_cosAvgsPionP->Fill(j,  TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionP));
	      p_sinAvgsPionM->Fill(j,  TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionM));
	      p_cosAvgsPionM->Fill(j,  TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionM));
	      p_sinAvgsProton->Fill(j,  TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiProton));
	      p_cosAvgsProton->Fill(j,  TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiProton));
	      p_sinAvgsEpd->Fill(j,  TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd));
	      p_cosAvgsEpd->Fill(j,  TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd));
	      p_sinAvgsEpdA->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdA));
	      p_cosAvgsEpdA->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdA));
	      p_sinAvgsEpdB->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdB));
	      p_cosAvgsEpdB->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdB));
	      p_sinAvgsEpdC->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC));
	      p_cosAvgsEpdC->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC));
	      p_sinAvgsEpdD->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdD));
	      p_cosAvgsEpdD->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdD));
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
      h_psiPionP_FLAT  = new TH1D("h_psiPionP_FLAT", "Flattened Event Plane Angles (n = "+ORDER_N_STR+", #pi^{+});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiPionM_FLAT  = new TH1D("h_psiPionM_FLAT", "Flattened Event Plane Angles (n = "+ORDER_N_STR+", #pi^{-});#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiProton_FLAT = new TH1D("h_psiProton_FLAT", "Flattened Event Plane Angles (n = "+ORDER_N_STR+", Protons);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpd_FLAT  = new TH1D("h_psiEpd_FLAT", "Flattened Event Plane Angle (EPD, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdA_FLAT = new TH1D("h_psiEpdA_FLAT", "Flattened Event Plane Angle (EPD A, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdB_FLAT = new TH1D("h_psiEpdB_FLAT", "Flattened Event Plane Angle (EPD B, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC_FLAT = new TH1D("h_psiEpdC_FLAT", "Flattened Event Plane Angle (EPD C, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdD_FLAT = new TH1D("h_psiEpdD_FLAT", "Flattened Event Plane Angle (EPD D, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

      TProfile *p_sinAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsTpc");
      TProfile *p_cosAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsTpc");
      TProfile *p_sinAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcA");
      TProfile *p_cosAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcA");
      TProfile *p_sinAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcB");
      TProfile *p_cosAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcB");
      TProfile *p_sinAvgsPionP_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsPionP");
      TProfile *p_cosAvgsPionP_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsPionP");
      TProfile *p_sinAvgsPionM_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsPionM");
      TProfile *p_cosAvgsPionM_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsPionM");
      TProfile *p_sinAvgsProton_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsProton");
      TProfile *p_cosAvgsProton_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsProton");
      TProfile *p_sinAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsEpd");
      TProfile *p_cosAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsEpd");
      TProfile *p_sinAvgsEpdA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdA");
      TProfile *p_cosAvgsEpdA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdA");
      TProfile *p_sinAvgsEpdB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdB");
      TProfile *p_cosAvgsEpdB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdB");
      TProfile *p_sinAvgsEpdC_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdC");
      TProfile *p_cosAvgsEpdC_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdC");
      TProfile *p_sinAvgsEpdD_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdD");
      TProfile *p_cosAvgsEpdD_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdD");


      // Get corrected event plane angles //
      for (Int_t i = 0; i < numOfEvents; i++)  // Loop over v_events to correct their angles
	{
	  if ( v_events.at(i).badEvent == true) { continue; }

	  Double_t psiTpc_delta  = 0;
	  Double_t psiTpcA_delta = 0;
	  Double_t psiTpcB_delta = 0;
	  Double_t psiPionP_delta  = 0;
	  Double_t psiPionM_delta  = 0;
	  Double_t psiProton_delta = 0;
	  Double_t psiEpd_delta  = 0;
	  Double_t psiEpdA_delta = 0;
	  Double_t psiEpdB_delta = 0;
	  Double_t psiEpdC_delta = 0;
	  Double_t psiEpdD_delta = 0;

	  Double_t jthSinAvg_Tpc  = 0;
	  Double_t jthCosAvg_Tpc  = 0;
	  Double_t jthSinAvg_TpcA = 0;
	  Double_t jthCosAvg_TpcA = 0;
	  Double_t jthSinAvg_TpcB = 0;
	  Double_t jthCosAvg_TpcB = 0;
	  Double_t jthSinAvg_PionP  = 0;
	  Double_t jthCosAvg_PionP  = 0;
	  Double_t jthSinAvg_PionM  = 0;
	  Double_t jthCosAvg_PionM  = 0;
	  Double_t jthSinAvg_Proton = 0;
	  Double_t jthCosAvg_Proton = 0;
	  Double_t jthSinAvg_Epd  = 0;
	  Double_t jthCosAvg_Epd  = 0;
	  Double_t jthSinAvg_EpdA = 0;
	  Double_t jthCosAvg_EpdA = 0;
	  Double_t jthSinAvg_EpdB = 0;
	  Double_t jthCosAvg_EpdB = 0;
	  Double_t jthSinAvg_EpdC = 0;
	  Double_t jthCosAvg_EpdC = 0;
	  Double_t jthSinAvg_EpdD = 0;
	  Double_t jthCosAvg_EpdD = 0;


	  for (Int_t j = 1; j <= SHIFT_TERMS; j++)    // Build the correction sums
	    {
	      jthSinAvg_Tpc  = p_sinAvgsTpc_INPUT->GetBinContent(j);
	      jthCosAvg_Tpc  = p_cosAvgsTpc_INPUT->GetBinContent(j);
	      jthSinAvg_TpcA = p_sinAvgsTpcA_INPUT->GetBinContent(j);
	      jthCosAvg_TpcA = p_cosAvgsTpcA_INPUT->GetBinContent(j);
	      jthSinAvg_TpcB = p_sinAvgsTpcB_INPUT->GetBinContent(j);
	      jthCosAvg_TpcB = p_cosAvgsTpcB_INPUT->GetBinContent(j);
	      jthSinAvg_PionP  = p_sinAvgsPionP_INPUT->GetBinContent(j);
	      jthCosAvg_PionP  = p_cosAvgsPionP_INPUT->GetBinContent(j);
	      jthSinAvg_PionM  = p_sinAvgsPionM_INPUT->GetBinContent(j);
	      jthCosAvg_PionM  = p_cosAvgsPionM_INPUT->GetBinContent(j);
	      jthSinAvg_Proton = p_sinAvgsProton_INPUT->GetBinContent(j);
	      jthCosAvg_Proton = p_cosAvgsProton_INPUT->GetBinContent(j);
	      jthSinAvg_Epd  = p_sinAvgsEpd_INPUT->GetBinContent(j);
	      jthCosAvg_Epd  = p_cosAvgsEpd_INPUT->GetBinContent(j);
	      jthSinAvg_EpdA = p_sinAvgsEpdA_INPUT->GetBinContent(j);
	      jthCosAvg_EpdA = p_cosAvgsEpdA_INPUT->GetBinContent(j);
	      jthSinAvg_EpdB = p_sinAvgsEpdB_INPUT->GetBinContent(j);
	      jthCosAvg_EpdB = p_cosAvgsEpdB_INPUT->GetBinContent(j);
	      jthSinAvg_EpdC = p_sinAvgsEpdC_INPUT->GetBinContent(j);
	      jthCosAvg_EpdC = p_cosAvgsEpdC_INPUT->GetBinContent(j);
	      jthSinAvg_EpdD = p_sinAvgsEpdD_INPUT->GetBinContent(j);
	      jthCosAvg_EpdD = p_cosAvgsEpdD_INPUT->GetBinContent(j);

	      psiTpc_delta  += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_Tpc*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc) 
									+jthCosAvg_Tpc*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc));
	      psiTpcA_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcA*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA) 
									+jthCosAvg_TpcA*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      psiTpcB_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcB*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB) 
									+jthCosAvg_TpcB*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      psiPionP_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_PionP*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionP) 
									 +jthCosAvg_PionP*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionP));
	      psiPionM_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_PionM*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionM) 
									 +jthCosAvg_PionM*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiPionM));
	      psiProton_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_Proton*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiProton) 
									  +jthCosAvg_Proton*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiProton));
	      psiEpd_delta  += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_Epd*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd)
									+jthCosAvg_Epd*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd));
	      psiEpdA_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdA*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdA)
									+jthCosAvg_EpdA*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdA));
	      psiEpdB_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdB*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdB)
									+jthCosAvg_EpdB*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdB));
	      psiEpdC_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdC*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC)
									+jthCosAvg_EpdC*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC));
	      psiEpdD_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdD*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdD)
									+jthCosAvg_EpdD*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdD));
	    }


	  v_events.at(i).psiTpc  += psiTpc_delta;
	  v_events.at(i).psiTpcA += psiTpcA_delta;
	  v_events.at(i).psiTpcB += psiTpcB_delta;
	  v_events.at(i).psiPionP  += psiPionP_delta;
	  v_events.at(i).psiPionM  += psiPionM_delta;
	  v_events.at(i).psiProton += psiProton_delta;
	  v_events.at(i).psiEpd  += psiEpd_delta;
	  v_events.at(i).psiEpdA += psiEpdA_delta;
	  v_events.at(i).psiEpdB += psiEpdB_delta;
	  v_events.at(i).psiEpdC += psiEpdC_delta;
	  v_events.at(i).psiEpdD += psiEpdD_delta;

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_N);
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiPionP  = angleShift(v_events.at(i).psiPionP, ORDER_N);
	  v_events.at(i).psiPionM  = angleShift(v_events.at(i).psiPionM, ORDER_N);
	  v_events.at(i).psiProton = angleShift(v_events.at(i).psiProton, ORDER_N);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_N);
	  v_events.at(i).psiEpdA = angleShift(v_events.at(i).psiEpdA, ORDER_N);
	  v_events.at(i).psiEpdB = angleShift(v_events.at(i).psiEpdB, ORDER_N);
	  v_events.at(i).psiEpdC = angleShift(v_events.at(i).psiEpdC, ORDER_N);
	  v_events.at(i).psiEpdD = angleShift(v_events.at(i).psiEpdD, ORDER_N);

	  h_psiTpc_FLAT->Fill(v_events.at(i).psiTpc);
	  h_psiTpcA_FLAT->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_FLAT->Fill(v_events.at(i).psiTpcB);
	  h_psiPionP_FLAT->Fill(v_events.at(i).psiPionP);
	  h_psiPionM_FLAT->Fill(v_events.at(i).psiPionM);
	  h_psiProton_FLAT->Fill(v_events.at(i).psiProton);
	  h_psiEpd_FLAT->Fill(v_events.at(i).psiEpd);
	  h_psiEpdA_FLAT->Fill(v_events.at(i).psiEpdA);
	  h_psiEpdB_FLAT->Fill(v_events.at(i).psiEpdB);
	  h_psiEpdC_FLAT->Fill(v_events.at(i).psiEpdC);
	  h_psiEpdD_FLAT->Fill(v_events.at(i).psiEpdD);


	  // 2D Correlations between event planes
	  h2_psiEpdATpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdA);
	  h2_psiEpdBTpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdB);
	  h2_psiEpdCTpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdC);
	  h2_psiEpdDTpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdD);

	  h2_psiEpdATpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdA);
	  h2_psiEpdBTpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdB);
	  h2_psiEpdCTpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdC);
	  h2_psiEpdDTpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdD);

	  h2_psiTpcATpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiTpcA);
	  //


	  // 1D correlation averages used in calculating resolution using the 3 sub-event method
	  p_TpcAB->Fill(v_events.at(i).centID,    TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiTpcB)));

	  p_TpcAEpdA->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdA)));
	  p_TpcAEpdB->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdB)));
	  p_TpcAEpdC->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdC)));
	  p_TpcAEpdD->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdD)));
	  p_TpcBEpdA->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdA)));
	  p_TpcBEpdB->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdB)));
	  p_TpcBEpdC->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdC)));
	  p_TpcBEpdD->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdD)));

	  p_EpdAEpdB->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiEpdA - v_events.at(i).psiEpdB)));
	  p_EpdAEpdC->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiEpdA - v_events.at(i).psiEpdC)));
	  p_EpdAEpdD->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiEpdA - v_events.at(i).psiEpdD)));
	  p_EpdBEpdC->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiEpdB - v_events.at(i).psiEpdC)));
	  p_EpdBEpdD->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiEpdB - v_events.at(i).psiEpdD)));
	  p_EpdCEpdD->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiEpdC - v_events.at(i).psiEpdD)));
	  //


	  // 2D searches through eta and centrality for correlations between detectors
	  // ONLY USE THIS SECTION IF THE EPD REGIONS COVER THE WHOLE EPD!! MIGHT NOT MAKE SENSE OTHERWISE

	  Int_t tpcTracksA = v_events.at(i).phiValuesTpcA.size();
	  Int_t tpcTracksB = v_events.at(i).phiValuesTpcB.size();
	  Int_t epdHitsA   = v_events.at(i).phiValuesEpdA.size();
	  Int_t epdHitsB   = v_events.at(i).phiValuesEpdB.size();
	  Int_t epdHitsC   = v_events.at(i).phiValuesEpdC.size();
	  Int_t epdHitsD   = v_events.at(i).phiValuesEpdD.size();
	  Double_t phiTpc;
	  Double_t etaTpc;
	  Double_t phiEpd;
	  Double_t etaEpd;
	  Double_t psiTpc  = v_events.at(i).psiTpc;
	  Double_t psiEpd  = v_events.at(i).psiEpd;
	  Double_t psiTpcA = v_events.at(i).psiTpcA;
	  Double_t psiTpcB = v_events.at(i).psiTpcB;
	  Int_t centralityID = v_events.at(i).centID;

	  for (int k = 0; k < epdHitsA; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdA.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdA.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsB; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdB.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdB.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsC; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdC.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdC.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsD; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdD.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdD.at(k);

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


  //=========================================================
  //          Use A Histogram To Get v_2
  //=========================================================
  /*
  Double_t cosTerm = 0;
  Double_t phi = 0;
  Double_t psi = 0;
  Double_t AB_n2 = 0;
  Double_t AC_n2 = 0;
  Double_t AD_n2 = 0;
  Double_t BC_n2 = 0;
  Double_t BD_n2 = 0;
  Double_t CD_n2 = 0;
  Double_t ATpc_n2 = 0;
  Double_t BTpc_n2 = 0;
  Double_t CTpc_n2 = 0;
  Double_t DTpc_n2 = 0;

  Double_t psiEpdA = 0;
  Double_t psiEpdB = 0;
  Double_t psiEpdC = 0;
  Double_t psiEpdC = 0;
  Double_t psiTpc  = 0;

  // CORRELATIONS AND FLOW COEFFICIENTS
  for (Int_t i = 0; i < numOfEvents; i++)
    {
      // SUBEVENT CORRELATIONS
      psiTpc  = v_events.at(i).psiTpc;
      psiEpdA = v_events.at(i).psiEpdA;
      psiEpdB = v_events.at(i).psiEpdB;
      psiEpdC = v_events.at(i).psiEpdC;
      psiEpdC = v_events.at(i).psiEpdC;

      AB_n2   = TMath::Cos(2.0 * (psiEpdA - psiEpdB));
      AC_n2   = TMath::Cos(2.0 * (psiEpdA - psiEpdC));
      AD_n2   = TMath::Cos(2.0 * (psiEpdA - psiEpdC));
      BC_n2   = TMath::Cos(2.0 * (psiEpdB - psiEpdC));
      BD_n2   = TMath::Cos(2.0 * (psiEpdB - psiEpdC));
      CD_n2   = TMath::Cos(2.0 * (psiEpdC - psiEpdC));
      ATpc_n2 = TMath::Cos(2.0 * (psiEpdA - psiTpc));
      BTpc_n2 = TMath::Cos(2.0 * (psiEpdB - psiTpc));
      CTpc_n2 = TMath::Cos(2.0 * (psiEpdC - psiTpc));
      DTpc_n2 = TMath::Cos(2.0 * (psiEpdC - psiTpc));
      
      //h_AB_n2->Fill(AB_n2);
      //h_AC_n2->Fill(AC_n2);
      //h_AD_n2->Fill(AD_n2);
      //h_BC_n2->Fill(BC_n2);
      //h_BD_n2->Fill(BD_n2);
      //h_CD_n2->Fill(CD_n2);
      

      switch (v_events.at(i).centID)
	{
	case 0: 
	  h_CD_n2_c00->Fill(CD_n2);
	  break; 
	case 1: 
	  h_CD_n2_c01->Fill(CD_n2);
	  break;
	case 2: 
	  h_CD_n2_c02->Fill(CD_n2);
	  break;
	case 3: 
	  h_CD_n2_c03->Fill(CD_n2);
	  break;
	case 4: 
	  h_CD_n2_c04->Fill(CD_n2);
	  break;
	case 5: 
	  h_CD_n2_c05->Fill(CD_n2);
	  break;
	case 6: 
	  h_CD_n2_c06->Fill(CD_n2);
	  break;
	case 7: 
	  h_CD_n2_c07->Fill(CD_n2);
	  break;
	case 8: 
	  h_CD_n2_c08->Fill(CD_n2);
	  break;
	case 9: 
	  h_CD_n2_c09->Fill(CD_n2);
	  break;
	case 10: 
	  h_CD_n2_c10->Fill(CD_n2);
	  break;
	case 11: 
	  h_CD_n2_c11->Fill(CD_n2);
	  break;
	case 12: 
	  h_CD_n2_c12->Fill(CD_n2);
	  break;
	case 13: 
	  h_CD_n2_c13->Fill(CD_n2);
	  break;
	case 14: 
	  h_CD_n2_c14->Fill(CD_n2);
	  break;
	case 15: 
	  h_CD_n2_c15->Fill(CD_n2);
	  break;
	}


      // GET V2 WITHOUT AUTOCORRELATIONS HERE
      for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	{
	  phi = v_events.at(i).phiValues.at(j);
	  psi = v_events.at(i).psiValues.at(j);  // This psi was calculated without particle 'j'
	  
	  cos = TMath::Cos(2 * (phi - psi));
	  h_v2Plot->Fill(cos);
	}
 
	}
  */
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
      p_sinAvgsPionP   ->Write();
      p_cosAvgsPionP   ->Write();
      p_sinAvgsPionM   ->Write();
      p_cosAvgsPionM   ->Write();
      p_sinAvgsProton   ->Write();
      p_cosAvgsProton   ->Write();
      p_sinAvgsEpd   ->Write();
      p_cosAvgsEpd   ->Write();
      p_sinAvgsEpdA  ->Write();
      p_cosAvgsEpdA  ->Write();
      p_sinAvgsEpdB  ->Write();
      p_cosAvgsEpdB  ->Write();
      p_sinAvgsEpdC  ->Write();
      p_cosAvgsEpdC  ->Write();
      p_sinAvgsEpdD  ->Write();
      p_cosAvgsEpdD  ->Write();
      h_XnTpc        ->Write();
      h_YnTpc        ->Write();
      h_XnTpcA       ->Write();
      h_YnTpcA       ->Write();
      h_XnTpcB       ->Write();
      h_YnTpcB       ->Write();
      h_XnPionP      ->Write();
      h_YnPionP      ->Write();
      h_XnPionM      ->Write();
      h_YnPionM      ->Write();
      h_XnProton     ->Write();
      h_YnProton     ->Write();
      h_XnEpd        ->Write();
      h_YnEpd        ->Write();
      h_XnEpdA       ->Write();
      h_YnEpdA       ->Write();
      h_XnEpdB       ->Write();
      h_YnEpdB       ->Write();
      h_XnEpdC       ->Write();
      h_YnEpdC       ->Write();
      h_XnEpdD       ->Write();
      h_YnEpdD       ->Write();

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
}//End EventPlane





      /* REMOVING AUTO-CORRELATIONS */
      /*
	for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	{
	Double_t newXn = v_events.at(i).Xn - v_events.at(i).pTValues.at(j) * TMath::Cos((Double_t)ORDER_N * v_events.at(i).phiValues.at(j));
	Double_t newYn = v_events.at(i).Yn - v_events.at(i).pTValues.at(j) * TMath::Sin((Double_t)ORDER_N * v_events.at(i).phiValues.at(j));
	Double_t newPsi = TMath::ATan2(newYn, newXn) / (Double_t)ORDER_N;
	v_events.at(i).psiValues.push_back(newPsi);
	}	  

      */

  /*
  TH1D *h_AB_n2_c00 = new TH1D("h_AB_n2_c00","2nd Order A-B Correlations (75%-80%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c01 = new TH1D("h_AB_n2_c01","2nd Order A-B Correlations (70%-75%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c02 = new TH1D("h_AB_n2_c02","2nd Order A-B Correlations (65%-70%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c03 = new TH1D("h_AB_n2_c03","2nd Order A-B Correlations (60%-65%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c04 = new TH1D("h_AB_n2_c04","2nd Order A-B Correlations (55%-60%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c05 = new TH1D("h_AB_n2_c05","2nd Order A-B Correlations (50%-55%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c06 = new TH1D("h_AB_n2_c06","2nd Order A-B Correlations (45%-50%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c07 = new TH1D("h_AB_n2_c07","2nd Order A-B Correlations (40%-45%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c08 = new TH1D("h_AB_n2_c08","2nd Order A-B Correlations (35%-40%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c09 = new TH1D("h_AB_n2_c09","2nd Order A-B Correlations (30%-35%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c10 = new TH1D("h_AB_n2_c10","2nd Order A-B Correlations (25%-30%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c11 = new TH1D("h_AB_n2_c11","2nd Order A-B Correlations (20%-25%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c12 = new TH1D("h_AB_n2_c12","2nd Order A-B Correlations (15%-20%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c13 = new TH1D("h_AB_n2_c13","2nd Order A-B Correlations (10%-15%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c14 = new TH1D("h_AB_n2_c14","2nd Order A-B Correlations (5%-10%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_AB_n2_c15 = new TH1D("h_AB_n2_c15","2nd Order A-B Correlations (0%-5%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);

  TH1D *h_AC_n2_c00 = new TH1D("h_AC_n2_c00","2nd Order A-C Correlations (75%-80%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c01 = new TH1D("h_AC_n2_c01","2nd Order A-C Correlations (70%-75%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c02 = new TH1D("h_AC_n2_c02","2nd Order A-C Correlations (65%-70%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c03 = new TH1D("h_AC_n2_c03","2nd Order A-C Correlations (60%-65%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c04 = new TH1D("h_AC_n2_c04","2nd Order A-C Correlations (55%-60%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c05 = new TH1D("h_AC_n2_c05","2nd Order A-C Correlations (50%-55%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c06 = new TH1D("h_AC_n2_c06","2nd Order A-C Correlations (45%-50%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c07 = new TH1D("h_AC_n2_c07","2nd Order A-C Correlations (40%-45%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c08 = new TH1D("h_AC_n2_c08","2nd Order A-C Correlations (35%-40%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c09 = new TH1D("h_AC_n2_c09","2nd Order A-C Correlations (30%-35%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c10 = new TH1D("h_AC_n2_c10","2nd Order A-C Correlations (25%-30%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c11 = new TH1D("h_AC_n2_c11","2nd Order A-C Correlations (20%-25%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c12 = new TH1D("h_AC_n2_c12","2nd Order A-C Correlations (15%-20%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c13 = new TH1D("h_AC_n2_c13","2nd Order A-C Correlations (10%-15%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c14 = new TH1D("h_AC_n2_c14","2nd Order A-C Correlations (5%-10%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AC_n2_c15 = new TH1D("h_AC_n2_c15","2nd Order A-C Correlations (0%-5%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);

  TH1D *h_AD_n2_c00 = new TH1D("h_AD_n2_c00","2nd Order A-D Correlations (75%-80%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c01 = new TH1D("h_AD_n2_c01","2nd Order A-D Correlations (70%-75%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c02 = new TH1D("h_AD_n2_c02","2nd Order A-D Correlations (65%-70%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c03 = new TH1D("h_AD_n2_c03","2nd Order A-D Correlations (60%-65%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c04 = new TH1D("h_AD_n2_c04","2nd Order A-D Correlations (55%-60%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c05 = new TH1D("h_AD_n2_c05","2nd Order A-D Correlations (50%-55%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c06 = new TH1D("h_AD_n2_c06","2nd Order A-D Correlations (45%-50%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c07 = new TH1D("h_AD_n2_c07","2nd Order A-D Correlations (40%-45%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c08 = new TH1D("h_AD_n2_c08","2nd Order A-D Correlations (35%-40%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c09 = new TH1D("h_AD_n2_c09","2nd Order A-D Correlations (30%-35%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c10 = new TH1D("h_AD_n2_c10","2nd Order A-D Correlations (25%-30%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c11 = new TH1D("h_AD_n2_c11","2nd Order A-D Correlations (20%-25%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c12 = new TH1D("h_AD_n2_c12","2nd Order A-D Correlations (15%-20%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c13 = new TH1D("h_AD_n2_c13","2nd Order A-D Correlations (10%-15%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c14 = new TH1D("h_AD_n2_c14","2nd Order A-D Correlations (5%-10%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_AD_n2_c15 = new TH1D("h_AD_n2_c15","2nd Order A-D Correlations (0%-5%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);

  TH1D *h_BC_n2_c00 = new TH1D("h_BC_n2_c00","2nd Order B-C Correlations (75%-80%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c01 = new TH1D("h_BC_n2_c01","2nd Order B-C Correlations (70%-75%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c02 = new TH1D("h_BC_n2_c02","2nd Order B-C Correlations (65%-70%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c03 = new TH1D("h_BC_n2_c03","2nd Order B-C Correlations (60%-65%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c04 = new TH1D("h_BC_n2_c04","2nd Order B-C Correlations (55%-60%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c05 = new TH1D("h_BC_n2_c05","2nd Order B-C Correlations (50%-55%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c06 = new TH1D("h_BC_n2_c06","2nd Order B-C Correlations (45%-50%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c07 = new TH1D("h_BC_n2_c07","2nd Order B-C Correlations (40%-45%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c08 = new TH1D("h_BC_n2_c08","2nd Order B-C Correlations (35%-40%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c09 = new TH1D("h_BC_n2_c09","2nd Order B-C Correlations (30%-35%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c10 = new TH1D("h_BC_n2_c10","2nd Order B-C Correlations (25%-30%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c11 = new TH1D("h_BC_n2_c11","2nd Order B-C Correlations (20%-25%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c12 = new TH1D("h_BC_n2_c12","2nd Order B-C Correlations (15%-20%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c13 = new TH1D("h_BC_n2_c13","2nd Order B-C Correlations (10%-15%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c14 = new TH1D("h_BC_n2_c14","2nd Order B-C Correlations (5%-10%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BC_n2_c15 = new TH1D("h_BC_n2_c15","2nd Order B-C Correlations (0%-5%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{b}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);

  TH1D *h_BD_n2_c00 = new TH1D("h_BD_n2_c00","2nd Order B-D Correlations (75%-80%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c01 = new TH1D("h_BD_n2_c01","2nd Order B-D Correlations (70%-75%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c02 = new TH1D("h_BD_n2_c02","2nd Order B-D Correlations (65%-70%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c03 = new TH1D("h_BD_n2_c03","2nd Order B-D Correlations (60%-65%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c04 = new TH1D("h_BD_n2_c04","2nd Order B-D Correlations (55%-60%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c05 = new TH1D("h_BD_n2_c05","2nd Order B-D Correlations (50%-55%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c06 = new TH1D("h_BD_n2_c06","2nd Order B-D Correlations (45%-50%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c07 = new TH1D("h_BD_n2_c07","2nd Order B-D Correlations (40%-45%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c08 = new TH1D("h_BD_n2_c08","2nd Order B-D Correlations (35%-40%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c09 = new TH1D("h_BD_n2_c09","2nd Order B-D Correlations (30%-35%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c10 = new TH1D("h_BD_n2_c10","2nd Order B-D Correlations (25%-30%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c11 = new TH1D("h_BD_n2_c11","2nd Order B-D Correlations (20%-25%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c12 = new TH1D("h_BD_n2_c12","2nd Order B-D Correlations (15%-20%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c13 = new TH1D("h_BD_n2_c13","2nd Order B-D Correlations (10%-15%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c14 = new TH1D("h_BD_n2_c14","2nd Order B-D Correlations (5%-10%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_BD_n2_c15 = new TH1D("h_BD_n2_c15","2nd Order B-D Correlations (0%-5%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);

  TH1D *h_CD_n2_c00 = new TH1D("h_CD_n2_c00","2nd Order C-D Correlations (75%-80%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c01 = new TH1D("h_CD_n2_c01","2nd Order C-D Correlations (70%-75%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c02 = new TH1D("h_CD_n2_c02","2nd Order C-D Correlations (65%-70%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c03 = new TH1D("h_CD_n2_c03","2nd Order C-D Correlations (60%-65%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c04 = new TH1D("h_CD_n2_c04","2nd Order C-D Correlations (55%-60%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c05 = new TH1D("h_CD_n2_c05","2nd Order C-D Correlations (50%-55%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c06 = new TH1D("h_CD_n2_c06","2nd Order C-D Correlations (45%-50%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c07 = new TH1D("h_CD_n2_c07","2nd Order C-D Correlations (40%-45%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c08 = new TH1D("h_CD_n2_c08","2nd Order C-D Correlations (35%-40%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c09 = new TH1D("h_CD_n2_c09","2nd Order C-D Correlations (30%-35%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c10 = new TH1D("h_CD_n2_c10","2nd Order C-D Correlations (25%-30%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c11 = new TH1D("h_CD_n2_c11","2nd Order C-D Correlations (20%-25%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c12 = new TH1D("h_CD_n2_c12","2nd Order C-D Correlations (15%-20%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c13 = new TH1D("h_CD_n2_c13","2nd Order C-D Correlations (10%-15%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c14 = new TH1D("h_CD_n2_c14","2nd Order C-D Correlations (5%-10%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_n2_c15 = new TH1D("h_CD_n2_c15","2nd Order C-D Correlations (0%-5%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);

  TH1D *h_ATpc_n2_c00 = new TH1D("h_ATpc_n2_c00","2nd Order A-TPC Correlations (75%-80%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c01 = new TH1D("h_ATpc_n2_c01","2nd Order A-TPC Correlations (70%-75%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c02 = new TH1D("h_ATpc_n2_c02","2nd Order A-TPC Correlations (65%-70%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c03 = new TH1D("h_ATpc_n2_c03","2nd Order A-TPC Correlations (60%-65%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c04 = new TH1D("h_ATpc_n2_c04","2nd Order A-TPC Correlations (55%-60%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c05 = new TH1D("h_ATpc_n2_c05","2nd Order A-TPC Correlations (50%-55%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c06 = new TH1D("h_ATpc_n2_c06","2nd Order A-TPC Correlations (45%-50%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c07 = new TH1D("h_ATpc_n2_c07","2nd Order A-TPC Correlations (40%-45%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c08 = new TH1D("h_ATpc_n2_c08","2nd Order A-TPC Correlations (35%-40%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c09 = new TH1D("h_ATpc_n2_c09","2nd Order A-TPC Correlations (30%-35%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c10 = new TH1D("h_ATpc_n2_c10","2nd Order A-TPC Correlations (25%-30%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c11 = new TH1D("h_ATpc_n2_c11","2nd Order A-TPC Correlations (20%-25%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c12 = new TH1D("h_ATpc_n2_c12","2nd Order A-TPC Correlations (15%-20%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c13 = new TH1D("h_ATpc_n2_c13","2nd Order A-TPC Correlations (10%-15%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c14 = new TH1D("h_ATpc_n2_c14","2nd Order A-TPC Correlations (5%-10%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_ATpc_n2_c15 = new TH1D("h_ATpc_n2_c15","2nd Order A-TPC Correlations (0%-5%);cos(2(#psi^{a}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);

  TH1D *h_BTpc_n2_c00 = new TH1D("h_BTpc_n2_c00","2nd Order B-TPC Correlations (75%-80%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c01 = new TH1D("h_BTpc_n2_c01","2nd Order B-TPC Correlations (70%-75%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c02 = new TH1D("h_BTpc_n2_c02","2nd Order B-TPC Correlations (65%-70%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c03 = new TH1D("h_BTpc_n2_c03","2nd Order B-TPC Correlations (60%-65%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c04 = new TH1D("h_BTpc_n2_c04","2nd Order B-TPC Correlations (55%-60%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c05 = new TH1D("h_BTpc_n2_c05","2nd Order B-TPC Correlations (50%-55%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c06 = new TH1D("h_BTpc_n2_c06","2nd Order B-TPC Correlations (45%-50%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c07 = new TH1D("h_BTpc_n2_c07","2nd Order B-TPC Correlations (40%-45%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c08 = new TH1D("h_BTpc_n2_c08","2nd Order B-TPC Correlations (35%-40%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c09 = new TH1D("h_BTpc_n2_c09","2nd Order B-TPC Correlations (30%-35%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c10 = new TH1D("h_BTpc_n2_c10","2nd Order B-TPC Correlations (25%-30%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c11 = new TH1D("h_BTpc_n2_c11","2nd Order B-TPC Correlations (20%-25%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c12 = new TH1D("h_BTpc_n2_c12","2nd Order B-TPC Correlations (15%-20%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c13 = new TH1D("h_BTpc_n2_c13","2nd Order B-TPC Correlations (10%-15%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c14 = new TH1D("h_BTpc_n2_c14","2nd Order B-TPC Correlations (5%-10%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_BTpc_n2_c15 = new TH1D("h_BTpc_n2_c15","2nd Order B-TPC Correlations (0%-5%);cos(2(#psi^{b}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);

  TH1D *h_CTpc_n2_c00 = new TH1D("h_CTpc_n2_c00","2nd Order C-TPC Correlations (75%-80%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c01 = new TH1D("h_CTpc_n2_c01","2nd Order C-TPC Correlations (70%-75%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c02 = new TH1D("h_CTpc_n2_c02","2nd Order C-TPC Correlations (65%-70%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c03 = new TH1D("h_CTpc_n2_c03","2nd Order C-TPC Correlations (60%-65%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c04 = new TH1D("h_CTpc_n2_c04","2nd Order C-TPC Correlations (55%-60%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c05 = new TH1D("h_CTpc_n2_c05","2nd Order C-TPC Correlations (50%-55%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c06 = new TH1D("h_CTpc_n2_c06","2nd Order C-TPC Correlations (45%-50%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c07 = new TH1D("h_CTpc_n2_c07","2nd Order C-TPC Correlations (40%-45%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c08 = new TH1D("h_CTpc_n2_c08","2nd Order C-TPC Correlations (35%-40%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c09 = new TH1D("h_CTpc_n2_c09","2nd Order C-TPC Correlations (30%-35%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c10 = new TH1D("h_CTpc_n2_c10","2nd Order C-TPC Correlations (25%-30%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c11 = new TH1D("h_CTpc_n2_c11","2nd Order C-TPC Correlations (20%-25%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c12 = new TH1D("h_CTpc_n2_c12","2nd Order C-TPC Correlations (15%-20%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c13 = new TH1D("h_CTpc_n2_c13","2nd Order C-TPC Correlations (10%-15%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c14 = new TH1D("h_CTpc_n2_c14","2nd Order C-TPC Correlations (5%-10%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_CTpc_n2_c15 = new TH1D("h_CTpc_n2_c15","2nd Order C-TPC Correlations (0%-5%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);

  TH1D *h_DTpc_n2_c00 = new TH1D("h_DTpc_n2_c00","2nd Order D-TPC Correlations (75%-80%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c01 = new TH1D("h_DTpc_n2_c01","2nd Order D-TPC Correlations (70%-75%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c02 = new TH1D("h_DTpc_n2_c02","2nd Order D-TPC Correlations (65%-70%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c03 = new TH1D("h_DTpc_n2_c03","2nd Order D-TPC Correlations (60%-65%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c04 = new TH1D("h_DTpc_n2_c04","2nd Order D-TPC Correlations (55%-60%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c05 = new TH1D("h_DTpc_n2_c05","2nd Order D-TPC Correlations (50%-55%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c06 = new TH1D("h_DTpc_n2_c06","2nd Order D-TPC Correlations (45%-50%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c07 = new TH1D("h_DTpc_n2_c07","2nd Order D-TPC Correlations (40%-45%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c08 = new TH1D("h_DTpc_n2_c08","2nd Order D-TPC Correlations (35%-40%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c09 = new TH1D("h_DTpc_n2_c09","2nd Order D-TPC Correlations (30%-35%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c10 = new TH1D("h_DTpc_n2_c10","2nd Order D-TPC Correlations (25%-30%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c11 = new TH1D("h_DTpc_n2_c11","2nd Order D-TPC Correlations (20%-25%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c12 = new TH1D("h_DTpc_n2_c12","2nd Order D-TPC Correlations (15%-20%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c13 = new TH1D("h_DTpc_n2_c13","2nd Order D-TPC Correlations (10%-15%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c14 = new TH1D("h_DTpc_n2_c14","2nd Order D-TPC Correlations (5%-10%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  TH1D *h_DTpc_n2_c15 = new TH1D("h_DTpc_n2_c15","2nd Order D-TPC Correlations (0%-5%);cos(2(#psi^{d}_{"+ORDER_N_STR+"}-#psi^{TPC}_{"+ORDER_N_STR+"}));Events",150,-1.5,1.5);
  */



      //if (eventInfo.psiTpc < 0) eventInfo.psiTpc += TMath::TwoPi();

      // Dividing TPC into subevents
      /*
      for (UInt_t i = 0; i < eventInfo.etaValues.size(); i++)
	{
	  if (eventInfo.etaValues.at(i) < 0 && eventInfo.etaValues.at(i) >= etaCut1)
	    {
	      eventInfo.phiValuesA.push_back(eventInfo.phiValues.at(i));
	      eventInfo.pTValuesA.push_back(eventInfo.pTValues.at(i));

	      eventInfo.XnA += eventInfo.pTValues.at(i) * TMath::Cos((Double_t)ORDER_N * eventInfo.phiValues.at(i));
	      eventInfo.YnA += eventInfo.pTValues.at(i) * TMath::Sin((Double_t)ORDER_N * eventInfo.phiValues.at(i));

	      h_etaA->Fill(eventInfo.etaValues.at(i));
	      eventInfo.nTracksA++;
	    }
	  else if (eventInfo.etaValues.at(i) < etaCut1 && eventInfo.etaValues.at(i) >= etaCut2)
	    {
	      eventInfo.phiValuesB.push_back(eventInfo.phiValues.at(i));
	      eventInfo.pTValuesB.push_back(eventInfo.pTValues.at(i));

	      eventInfo.XnB += eventInfo.pTValues.at(i) * TMath::Cos((Double_t)ORDER_N * eventInfo.phiValues.at(i));
	      eventInfo.YnB += eventInfo.pTValues.at(i) * TMath::Sin((Double_t)ORDER_N * eventInfo.phiValues.at(i));

	      h_etaB->Fill(eventInfo.etaValues.at(i));
	      eventInfo.nTracksB++;
	    }
	  else if (eventInfo.etaValues.at(i) < etaCut2 && eventInfo.etaValues.at(i) >= -1.5)
	    {
	      eventInfo.phiValuesC.push_back(eventInfo.phiValues.at(i));
	      eventInfo.pTValuesC.push_back(eventInfo.pTValues.at(i));

	      eventInfo.XnC += eventInfo.pTValues.at(i) * TMath::Cos((Double_t)ORDER_N * eventInfo.phiValues.at(i));
	      eventInfo.YnC += eventInfo.pTValues.at(i) * TMath::Sin((Double_t)ORDER_N * eventInfo.phiValues.at(i));

	      h_etaC->Fill(eventInfo.etaValues.at(i));
	      eventInfo.nTracksC++;
	    }
	}

      if (eventInfo.nTracks < 18) continue;                  // Make sure there are at least 6 GOOD tracks in each subevent
      if (eventInfo.nTracksA < 6) continue;
      if (eventInfo.nTracksB < 6) continue;
      if (eventInfo.nTracksC < 6) continue;

      eventInfo.psi  = TMath::ATan2(eventInfo.Yn,  eventInfo.Xn)  / (Double_t)ORDER_N;
      eventInfo.psiA = TMath::ATan2(eventInfo.YnA, eventInfo.XnA) / (Double_t)ORDER_N;
      eventInfo.psiB = TMath::ATan2(eventInfo.YnB, eventInfo.XnB) / (Double_t)ORDER_N;
      eventInfo.psiC = TMath::ATan2(eventInfo.YnC, eventInfo.XnC) / (Double_t)ORDER_N;

      h_Xn->Fill(eventInfo.Xn);
      h_Yn->Fill(eventInfo.Yn);
      h_psi->Fill(eventInfo.psi);

      h_XnA->Fill(eventInfo.XnA);
      h_YnA->Fill(eventInfo.YnA);
      h_psiA->Fill(eventInfo.psiA);

      h_XnB->Fill(eventInfo.XnB);
      h_YnB->Fill(eventInfo.YnB);
      h_psiB->Fill(eventInfo.psiB);

      h_XnC->Fill(eventInfo.XnC);
      h_YnC->Fill(eventInfo.YnC);
      h_psi_C->Fill(eventInfo.psiC);
      */
