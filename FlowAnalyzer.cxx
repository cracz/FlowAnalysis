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
const Int_t ORDER_N = 2;                // Order of anisotropic flow
TString ORDER_N_STR;                   // ORDER_N but as a TString for titles/labels

//Double_t PSI_BOUNDS;
const Double_t PSI_BOUNDS = TMath::Pi()/ORDER_N + 1;
const Double_t Q_BOUNDS = 100;

const Int_t EPD_FORMAT       = 2;       // format=0/1/2 for StEpdHit/StMuEpdHit/StPicoEpdHit
const Int_t EPD_MAX_WEIGHT   = 2;      // max nMIP weight; recommended value, but variable
const Double_t EPD_THRESHOLD = 0.3;    // recommended value, but variable

const Double_t MIN_TPC_ETA_CUT = -2.0;
const Double_t AGAP_TPC_ETA_CUT = -1.0;//-1.1;
const Double_t GAPB_TPC_ETA_CUT = -1.0;
const Double_t MAX_TPC_ETA_CUT = 0.0;

const Double_t MIN_ETA_CUT = -5.3;//-5.16;    // Cuts for the EPD subevents (A, B, C, D)
const Double_t AB_ETA_CUT  = -4.05;//-3.82;
const Double_t BC_ETA_CUT  = -3.3;//-3.28;
const Double_t CD_ETA_CUT  = -2.9;//-2.87;
const Double_t MAX_ETA_CUT = -2.60;

const Double_t R_VTX_CUT = 2.0;         // 2D r value, good vertices are within this value
const Double_t Z_VTX_CUT_LOW  = 200.0;
const Double_t Z_VTX_CUT_HIGH = 201.5;

const Int_t MIN_TRACKS = 5;             // Min number of tracks/hits in each sub-event

const Int_t SHIFT_TERMS = 10;           // Number of terms to use when shifting event plane angles
const Int_t CENT_BINS = 8;             // Number of centrality bins (max 16)
const Int_t FIRST_CENT = 16 - CENT_BINS;            // Starting point for centrality dependent plots
const Int_t BAD_VALUE = -99;
const Double_t Y_MID = -1.05;          // Mid rapidity

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
  
  Int_t nTracksTpc;      // Number of GOOD tracks in the sub-event
  Double_t XnTpc;
  Double_t YnTpc;
  Double_t psiTpc;       // Overall EP angle without removing autocorrelations
  //std::vector<Double_t> phiValuesTpc;   // Azimuthal angles for all TPC particles in the event
  //std::vector<Double_t> etaValuesTpc;  // Eta values for TPC particles in the event
  //std::vector<Double_t> pTValuesTpc;   // pT values for all TPC particles in the event
  //std::vector<Double_t> psiValuesTpc;  // EP angles after removing autocorrelations

  Int_t nTracksTpcA;      // Number of GOOD tracks in the sub-event
  Double_t XnTpcA;
  Double_t YnTpcA;
  Double_t psiTpcA;       // Overall EP angle without removing autocorrelations
  std::vector<Double_t> phiValuesTpcA;   // Azimuthal angles for all TPC A particles in the event
  std::vector<Double_t> etaValuesTpcA;  // Eta values for TPC A particles in the event
  //std::vector<Double_t> pTValuesTpcA;   // pT values for all TPC A particles in the event

  Int_t nTracksTpcB;
  Double_t XnTpcB;
  Double_t YnTpcB;
  Double_t psiTpcB;
  std::vector<Double_t> phiValuesTpcB;
  std::vector<Double_t> etaValuesTpcB;
  //std::vector<Double_t> pTValuesTpcB;

  Int_t nhitsEpd;      // Number of GOOD tracks in the sub-event
  Double_t XnEpd;
  Double_t YnEpd;
  Double_t psiEpd;       // Overall EP angle without removing autocorrelations
  //std::vector<Double_t> phiValuesEpd;   // Azimuthal angles for all EPD particles in the event
  //std::vector<Double_t> etaValuesEpd;  // Eta values for EPD particles in the event
  //std::vector<Double_t> psiValuesEpd;  // EP angles after removing autocorrelations


  Int_t nhitsEpdA;
  Double_t XnEpdA;
  Double_t YnEpdA;
  Double_t psiEpdA;       // Subevent A EP angle without removing autocorrelations
  std::vector<Double_t> phiValuesEpdA;  // Subevent A azimuthal angles
  std::vector<Double_t> etaValuesEpdA;  // Subevent A eta values

  Int_t nhitsEpdB;
  Double_t XnEpdB;
  Double_t YnEpdB;
  Double_t psiEpdB;       // Subevent B EP angle without removing autocorrelations  
  std::vector<Double_t> phiValuesEpdB;  // Subevent B azimuthal angles
  std::vector<Double_t> etaValuesEpdB;  // Subevent B eta values

  Int_t nhitsEpdC;
  Double_t XnEpdC;
  Double_t YnEpdC;
  Double_t psiEpdC;       // Subevent C EP angle without removing autocorrelations  
  std::vector<Double_t> phiValuesEpdC;  // Subevent C azimuthal angles
  std::vector<Double_t> etaValuesEpdC;  // Subevent C eta values

  Int_t nhitsEpdD;
  Double_t XnEpdD;
  Double_t YnEpdD;
  Double_t psiEpdD;       // Subevent D EP angle without removing autocorrelations  
  std::vector<Double_t> phiValuesEpdD;  // Subevent D azimuthal angles
  std::vector<Double_t> etaValuesEpdD;  // Subevent D eta values

  void reset()
  {
    badEvent  = false;  //Reset all values in the struct to reuse
    primTracks = 0;
    centID = BAD_VALUE;

    nTracksTpc = 0;
    XnTpc = 0;
    YnTpc = 0;
    psiTpc = BAD_VALUE;
    //phiValuesTpc.clear();
    //etaValuesTpc.clear();
    //pTValuesTpc.clear();
    //psiValuesTpc.clear();

    nTracksTpcA = 0;
    XnTpcA = 0;
    YnTpcA = 0;
    psiTpcA = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcA);
    std::vector<Double_t>().swap(etaValuesTpcA);
    //pTValuesTpcA.clear();

    nTracksTpcB = 0;
    XnTpcB = 0;
    YnTpcB = 0;
    psiTpcB = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcB);
    std::vector<Double_t>().swap(etaValuesTpcB);
    //pTValuesTpcB.clear();

    nhitsEpd = 0;
    XnEpd = 0;
    YnEpd = 0;
    psiEpd = BAD_VALUE;
    //phiValuesEpd.clear();
    //etaValuesEpd.clear();
    //psiValuesEpd.clear();


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


// Moves event plane angles back into the correct period.
Double_t angleShift(Double_t angle, Int_t order)
{
  if (angle < -TMath::Pi()/(Double_t)order) { angle += TMath::TwoPi()/(Double_t)order; }
  else if (angle >  TMath::Pi()/(Double_t)order) { angle -= TMath::TwoPi()/(Double_t)order; }
  return angle;
};


// Checks if an event has any flow vectors equal to zero. Updates the event's member variable "badEvent".
void checkZeroQ(Event event)
{
  if (event.XnTpc == 0 && event.YnTpc == 0) { event.badEvent = true; }
  else if (event.XnEpd == 0 && event.YnEpd == 0) { event.badEvent = true; }
  else if (event.XnTpcA == 0 && event.YnTpcA == 0) { event.badEvent = true; }
  else if (event.XnTpcB == 0 && event.YnTpcB == 0) { event.badEvent = true; }
  else if (event.XnEpdA == 0 && event.YnEpdA == 0) { event.badEvent = true; }
  else if (event.XnEpdB == 0 && event.YnEpdB == 0) { event.badEvent = true; }
  else if (event.XnEpdC == 0 && event.YnEpdC == 0) { event.badEvent = true; }
  else if (event.XnEpdD == 0 && event.YnEpdD == 0) { event.badEvent = true; }
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
  /*
  TH1D *h_eta   = new TH1D("h_eta", "Particle #eta;#eta;Particles", 600, -6, 0);
  TH1D *h_etaTpc = new TH1D("h_etaTpc", "Particle #eta (TPC);#eta;Particles", 600, -6, 0);
  TH1D *h_etaEpdA = new TH1D("h_etaEpdA", "Particle #eta (EPD A);#eta;Particles", 600, -6, 0);
  TH1D *h_etaEpdB = new TH1D("h_etaEpdB", "Particle #eta (EPD B);#eta;Particles", 600, -6, 0);
  TH1D *h_etaEpdC = new TH1D("h_etaEpdC", "Particle #eta (EPD C);#eta;Particles", 600, -6, 0);
  TH1D *h_etaEpdD = new TH1D("h_etaEpdD", "Particle #eta (EPD D);#eta;Particles", 600, -6, 0);
  */
  TH1D *h_tileWeights = new TH1D("h_tileWeights", "EPD Tile Weights;Hits;nMIP Weights", 5, -1, 4);
  TH1D *h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TH1D *h_psiTpc_RAW  = new TH1D("h_psiTpc_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcA_RAW = new TH1D("h_psiTpcA_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_RAW = new TH1D("h_psiTpcB_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpd_RAW  = new TH1D("h_psiEpd_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdA_RAW = new TH1D("h_psiEpdA_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdB_RAW = new TH1D("h_psiEpdB_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdC_RAW = new TH1D("h_psiEpdC_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD C);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdD_RAW = new TH1D("h_psiEpdD_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD D);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  //TH1D *h_v2Plot = new TH1D("h_v2Plot", "Plot to Retrieve v_{2};cos(2(#phi - #psi_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);

  TH1D *h_CD_c00 = new TH1D("h_CD_c00","C-D Correlations (75%-80%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c01 = new TH1D("h_CD_c01","C-D Correlations (70%-75%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c02 = new TH1D("h_CD_c02","C-D Correlations (65%-70%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c03 = new TH1D("h_CD_c03","C-D Correlations (60%-65%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c04 = new TH1D("h_CD_c04","C-D Correlations (55%-60%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c05 = new TH1D("h_CD_c05","C-D Correlations (50%-55%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c06 = new TH1D("h_CD_c06","C-D Correlations (45%-50%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c07 = new TH1D("h_CD_c07","C-D Correlations (40%-45%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c08 = new TH1D("h_CD_c08","C-D Correlations (35%-40%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c09 = new TH1D("h_CD_c09","C-D Correlations (30%-35%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c10 = new TH1D("h_CD_c10","C-D Correlations (25%-30%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c11 = new TH1D("h_CD_c11","C-D Correlations (20%-25%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c12 = new TH1D("h_CD_c12","C-D Correlations (15%-20%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c13 = new TH1D("h_CD_c13","C-D Correlations (10%-15%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c14 = new TH1D("h_CD_c14","C-D Correlations (5%-10%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);
  TH1D *h_CD_c15 = new TH1D("h_CD_c15","C-D Correlations (0%-5%);cos(2(#psi^{c}_{"+ORDER_N_STR+"}-#psi^{d}_{"+ORDER_N_STR+"}));Particles",150,-1.5,1.5);


  //TH1D *h_abCorr_n2 = new TH1D("h_abCorr_n2", "2nd Order Subevent a-b Correlations;cos(2(#psi^{a}_{"+ORDER_N_STR+"} - #psi^{b}_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);
  //TH1D *h_acCorr_n2 = new TH1D("h_acCorr_n2", "2nd Order Subevent a-c Correlations;cos(2(#psi^{a}_{"+ORDER_N_STR+"} - #psi^{c}_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);
  //TH1D *h_bcCorr_n2 = new TH1D("h_bcCorr_n2", "2nd Order Subevent b-c Correlations;cos(2(#psi^{b}_{"+ORDER_N_STR+"} - #psi^{c}_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);

  /*
  TH2D *h2_betap = new TH2D("h2_betap","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_m2_p   = new TH2D("h2_m2_p", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 500, -3, 3, 500, -0.1, 15);
  */
  TH2D *h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  TH2D *h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  TH2D *h2_pp_vs_eta = new TH2D("h2_pp_vs_eta","Tile Weight for Supersectors vs #eta;#eta;Supersector", 400, -6, -2, 12, 0.5, 12.5);


  // Here the name refers to the eta region that will be displayed/searched using the event plane angle from the opposite region
  TProfile2D *h2_v2SearchTpc = new TProfile2D("h2_v2SearchTpc", "<cos(2(#phi^{TPC} - #psi^{EPD}_{"+ORDER_N_STR+"}))>;#eta;Centrality (%)", 
					      12, -2, 0, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2SearchEpd = new TProfile2D("h2_v2SearchEpd", "<cos(2(#phi^{EPD} - #psi^{TPC}_{"+ORDER_N_STR+"}))>;#eta;Centrality (%)", 
					      12, -5.3, -2.6, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2SearchEpdTpcB = new TProfile2D("h2_v2SearchEpdTpcB", "<cos(2(#phi^{EPD} - #psi^{TPC,B}_{"+ORDER_N_STR+"}))>;#eta;Centrality (%)", 
					      12, -5.3, -2.6, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  h2_v2SearchTpc->SetStats(0);
  h2_v2SearchEpd->SetStats(0);
  h2_v2SearchEpdTpcB->SetStats(0);
  // The indices here are equivalent to the corresponding centrality ID
  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};

  // The indices here are opposite to the corresponding centrality ID (array is backward)
  //const char *centralityBins[16] = {"0-5", "5-10", "10-15" "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80"};


  TH2D *h2_psiEpdAB = new TH2D("h2_psiEpdAB", "#psi^{EPD}_{A} vs #psi^{EPD}_{B} (Order "+ORDER_N_STR+");#psi^{EPD}_{B};#psi^{EPD}_{A}", 
			       200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdAC = new TH2D("h2_psiEpdAC", "#psi^{EPD}_{A} vs #psi^{EPD}_{C} (Order "+ORDER_N_STR+");#psi^{EPD}_{C};#psi^{EPD}_{A}", 
			       200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdAD = new TH2D("h2_psiEpdAD", "#psi^{EPD}_{A} vs #psi^{EPD}_{D} (Order "+ORDER_N_STR+");#psi^{EPD}_{D};#psi^{EPD}_{A}", 
			       200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBC = new TH2D("h2_psiEpdBC", "#psi^{EPD}_{B} vs #psi^{EPD}_{C} (Order "+ORDER_N_STR+");#psi^{EPD}_{C};#psi^{EPD}_{B}", 
			       200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBD = new TH2D("h2_psiEpdBD", "#psi^{EPD}_{B} vs #psi^{EPD}_{D} (Order "+ORDER_N_STR+");#psi^{EPD}_{D};#psi^{EPD}_{B}", 
			       200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdCD = new TH2D("h2_psiEpdCD", "#psi^{EPD}_{C} vs #psi^{EPD}_{D} (Order "+ORDER_N_STR+");#psi^{EPD}_{D};#psi^{EPD}_{C}", 
			       200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  
  TH2D *h2_psiEpdATpcA = new TH2D("h2_psiEpdATpcA", "#psi^{EPD}_{A} vs #psi^{TPC}_{A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcA = new TH2D("h2_psiEpdBTpcA", "#psi^{EPD}_{B} vs #psi^{TPC}_{A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdCTpcA = new TH2D("h2_psiEpdCTpcA", "#psi^{EPD}_{C} vs #psi^{TPC}_{A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{C}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdDTpcA = new TH2D("h2_psiEpdDTpcA", "#psi^{EPD}_{D} vs #psi^{TPC}_{A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{D}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiEpdATpcB = new TH2D("h2_psiEpdATpcB", "#psi^{EPD}_{A} vs #psi^{TPC}_{B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcB = new TH2D("h2_psiEpdBTpcB", "#psi^{EPD}_{B} vs #psi^{TPC}_{B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdCTpcB = new TH2D("h2_psiEpdCTpcB", "#psi^{EPD}_{C} vs #psi^{TPC}_{B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{C}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdDTpcB = new TH2D("h2_psiEpdDTpcB", "#psi^{EPD}_{D} vs #psi^{TPC}_{B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{D}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  
  TH2D *h2_psiTpcATpcB = new TH2D("h2_psiTpcATpcB", "#psi^{TPC}_{A} vs #psi^{TPC}_{B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{TPC}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);



  // CORRECTION HISTOGRAMS
  TProfile *p_sinAvgsTpc = new TProfile("p_sinAvgsTpc", "Sin Averages (TPC);j (Correction term);<sin(jn#psi^{TPC}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpc = new TProfile("p_cosAvgsTpc", "Cos Averages (TPC);j (Correction term);<sin(jn#psi^{TPC}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpd = new TProfile("p_sinAvgsEpd", "Sin Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpd = new TProfile("p_cosAvgsEpd", "Cos Averages (EPD);j (Correction term);<sin(jn#psi^{EPD}_{n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);

  TProfile *p_sinAvgsTpcA = new TProfile("p_sinAvgsTpcA", "Sin Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpcA = new TProfile("p_cosAvgsTpcA", "Cos Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsTpcB = new TProfile("p_sinAvgsTpcB", "Sin Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpcB = new TProfile("p_cosAvgsTpcB", "Cos Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdA = new TProfile("p_sinAvgsEpdA", "Sin Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdA = new TProfile("p_cosAvgsEpdA", "Cos Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdB = new TProfile("p_sinAvgsEpdB", "Sin Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdB = new TProfile("p_cosAvgsEpdB", "Cos Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdC = new TProfile("p_sinAvgsEpdC", "Sin Averages (EPD C);j (Correction term);<sin(jn#psi^{EPD}_{C,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdC = new TProfile("p_cosAvgsEpdC", "Cos Averages (EPD C);j (Correction term);<sin(jn#psi^{EPD}_{C,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdD = new TProfile("p_sinAvgsEpdD", "Sin Averages (EPD D);j (Correction term);<sin(jn#psi^{EPD}_{D,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdD = new TProfile("p_cosAvgsEpdD", "Cos Averages (EPD D);j (Correction term);<sin(jn#psi^{EPD}_{D,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);

  TH1D *h_XnTpc = new TH1D("h_XnTpc", "X_n Distribution (TPC);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpc = new TH1D("h_YnTpc", "Y_n Distribution (TPC);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpd = new TH1D("h_XnEpd", "X_n Distribution (EPD);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpd = new TH1D("h_YnEpd", "Y_n Distribution (EPD);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);

  TH1D *h_XnTpcA = new TH1D("h_XnTpcA", "X_n Distribution (TPC A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcA = new TH1D("h_YnTpcA", "Y_n Distribution (TPC A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcB = new TH1D("h_XnTpcB", "X_n Distribution (TPC B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcB = new TH1D("h_YnTpcB", "Y_n Distribution (TPC B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdA = new TH1D("h_XnEpdA", "X_n Distribution (EPD A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdA = new TH1D("h_YnEpdA", "Y_n Distribution (EPD A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdB = new TH1D("h_XnEpdB", "X_n Distribution (EPD B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdB = new TH1D("h_YnEpdB", "Y_n Distribution (EPD B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdC = new TH1D("h_XnEpdC", "X_n Distribution (EPD C);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdC = new TH1D("h_YnEpdC", "Y_n Distribution (EPD C);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdD = new TH1D("h_XnEpdD", "X_n Distribution (EPD D);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdD = new TH1D("h_YnEpdD", "Y_n Distribution (EPD D);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);


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

      //Bool_t b_bad_xvtx = ( (d_xvtx < -1.0) || (d_xvtx > 1.0)  );
      //Bool_t b_bad_yvtx = ( (d_yvtx < -2.6) || (d_yvtx > -1.4) );
      Bool_t b_bad_rvtx = ( d_rvtx >= R_VTX_CUT );
      Bool_t b_bad_zvtx = ( (d_zvtx <= Z_VTX_CUT_LOW)  || (d_zvtx >= Z_VTX_CUT_HIGH));

      if (b_bad_zvtx) continue;

      h2_trans_vtx->Fill(d_xvtx, d_yvtx);

      //if (b_bad_xvtx || b_bad_yvtx) continue;
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
	  /*
	  if (!picoTrack->isTofTrack()) continue;  // Only TOF tracks

	  StPicoBTofPidTraits *trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());

	  Double_t d_tofBeta = trait->btofBeta();

	  if (d_tofBeta < 0.01) continue;
	  */
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

	  TVector3 mom_vec  = picoTrack->pMom();
	  Double_t d_charge = picoTrack->charge();
	  Double_t d_mom    = picoTrack->pPtot();
	  Double_t d_pT     = picoTrack->pPt();
	  Double_t d_phi    = mom_vec.Phi();
	  Double_t d_eta    = mom_vec.Eta();

	  // Fill histos and save important event info in the custom struct type
	  if (d_charge != 0)
	    {
	      //h2_betap->Fill(d_charge*d_mom, 1/d_tofBeta);
	      //h2_m2_p->Fill(d_charge*d_mom, d_mom*d_mom*(1/(d_tofBeta*d_tofBeta) - 1));

	      eventInfo.nTracksTpc++;
	      //eventInfo.phiValuesTpc.push_back(d_phi);
	      //eventInfo.etaValuesTpc.push_back(d_eta);
	      eventInfo.XnTpc += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
	      eventInfo.YnTpc += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);


	      if (d_eta > MIN_TPC_ETA_CUT && d_eta < AGAP_TPC_ETA_CUT)
		{
		  eventInfo.phiValuesTpcA.push_back(d_phi);
		  eventInfo.etaValuesTpcA.push_back(d_eta);

		  eventInfo.nTracksTpcA++;
		  //eventInfo.pTValuesTpcA.push_back(d_pT);
		  eventInfo.XnTpcA += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
		  eventInfo.YnTpcA += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
		}
	      else if (d_eta /*>*/>= GAPB_TPC_ETA_CUT && d_eta < MAX_TPC_ETA_CUT)
		{
		  eventInfo.phiValuesTpcB.push_back(d_phi);
		  eventInfo.etaValuesTpcB.push_back(d_eta);

		  eventInfo.nTracksTpcB++;
		  //eventInfo.pTValuesTpcB.push_back(d_pT);
		  eventInfo.XnTpcB += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
		  eventInfo.YnTpcB += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
		}

	    }
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

	  eventInfo.nhitsEpd++;
	  eventInfo.XnEpd += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	  eventInfo.YnEpd += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	  //eventInfo.etaValuesEpd.push_back(tileEta);
	  //eventInfo.phiValuesEpd.push_back(tilePhi);

	  if (tileEta > MIN_ETA_CUT && tileEta < AB_ETA_CUT)  // Sub A
	    {
	      eventInfo.nhitsEpdA++;
	      eventInfo.XnEpdA += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdA += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.etaValuesEpdA.push_back(tileEta);
	      eventInfo.phiValuesEpdA.push_back(tileEta);
	    }
	  else if (tileEta /*>*/>= AB_ETA_CUT && tileEta < BC_ETA_CUT)  // Sub B
	    {
	      eventInfo.nhitsEpdB++;
	      eventInfo.XnEpdB += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdB += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.etaValuesEpdB.push_back(tileEta);
	      eventInfo.phiValuesEpdB.push_back(tileEta);
	    }
	  else if (tileEta /*>*/>= BC_ETA_CUT && tileEta < CD_ETA_CUT)  // Sub C
	    {
	      eventInfo.nhitsEpdC++;
	      eventInfo.XnEpdC += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdC += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.etaValuesEpdC.push_back(tileEta);
	      eventInfo.phiValuesEpdC.push_back(tileEta);
	    }
	  else if (tileEta /*>*/>= CD_ETA_CUT && tileEta < MAX_ETA_CUT)  // Sub D
	    {
	      eventInfo.nhitsEpdD++;
	      eventInfo.XnEpdD += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdD += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.etaValuesEpdD.push_back(tileEta);
	      eventInfo.phiValuesEpdD.push_back(tileEta);
	    }

	}
      delete epdGeom;

      //=========================================================
      //            END EPD STUFF
      //=========================================================



      if (eventInfo.nTracksTpc < MIN_TRACKS) continue;
      if (eventInfo.nTracksTpcA < MIN_TRACKS) continue;
      if (eventInfo.nTracksTpcB < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpd < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdA < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdB < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdC < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdD < MIN_TRACKS) continue;
      if (eventInfo.XnTpc == 0 && eventInfo.YnTpc == 0) continue;
      if (eventInfo.XnEpd == 0 && eventInfo.YnEpd == 0) continue;
      if (eventInfo.XnTpcA == 0 && eventInfo.YnTpcA == 0) continue;
      if (eventInfo.XnTpcB == 0 && eventInfo.YnTpcB == 0) continue;
      if (eventInfo.XnEpdA == 0 && eventInfo.YnEpdA == 0) continue;
      if (eventInfo.XnEpdB == 0 && eventInfo.YnEpdB == 0) continue;
      if (eventInfo.XnEpdC == 0 && eventInfo.YnEpdC == 0) continue;
      if (eventInfo.XnEpdD == 0 && eventInfo.YnEpdD == 0) continue;


      // RAW SUB-EVENT PLANE ANGLES //
      if (ORDER_N % 2 == 1)                 // Q vectors of EPD East are opposite sign of West for odd harmonics. Switch back to compare to TPC.
	{
	  eventInfo.XnEpd  *= -1.0;
	  eventInfo.YnEpd  *= -1.0;
	  eventInfo.XnEpdA *= -1.0;
	  eventInfo.YnEpdA *= -1.0;
	  eventInfo.XnEpdB *= -1.0;
	  eventInfo.YnEpdB *= -1.0;
	  eventInfo.XnEpdC *= -1.0;
	  eventInfo.YnEpdC *= -1.0;
	  eventInfo.XnEpdD *= -1.0;
	  eventInfo.YnEpdD *= -1.0;
	}
      eventInfo.psiTpc  = TMath::ATan2(eventInfo.YnTpc,  eventInfo.XnTpc)  / (Double_t)ORDER_N;
      eventInfo.psiEpd  = TMath::ATan2(eventInfo.YnEpd,  eventInfo.XnEpd)  / (Double_t)ORDER_N;
      eventInfo.psiTpcA = TMath::ATan2(eventInfo.YnTpcA, eventInfo.XnTpcA) / (Double_t)ORDER_N;
      eventInfo.psiTpcB = TMath::ATan2(eventInfo.YnTpcB, eventInfo.XnTpcB) / (Double_t)ORDER_N;
      eventInfo.psiEpdA = TMath::ATan2(eventInfo.YnEpdA, eventInfo.XnEpdA) / (Double_t)ORDER_N;
      eventInfo.psiEpdB = TMath::ATan2(eventInfo.YnEpdB, eventInfo.XnEpdB) / (Double_t)ORDER_N;
      eventInfo.psiEpdC = TMath::ATan2(eventInfo.YnEpdC, eventInfo.XnEpdC) / (Double_t)ORDER_N;
      eventInfo.psiEpdD = TMath::ATan2(eventInfo.YnEpdD, eventInfo.XnEpdD) / (Double_t)ORDER_N;


      // Filling histos here since this is past all possible cuts
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcA.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcA.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcB.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcB.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdA.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdA.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdB.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdB.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdC.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdD.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdD.at(i) - Y_MID ); }
      //for (unsigned int i = 0; i < eventInfo.etaValues.size();   i++)  { h_eta->Fill(eventInfo.etaValues.at(i)); }
      //for (unsigned int i = 0; i < eventInfo.pTValuesTpcA.size(); i++)  { h_pTA->Fill(eventInfo.pTValuesTpcA.at(i)); }
      //for (unsigned int i = 0; i < eventInfo.pTValuesTpcB.size(); i++)  { h_pTB->Fill(eventInfo.pTValuesTpcB.at(i)); }

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
      h_XnEpdA->Fill(eventInfo.XnEpdA);
      h_YnEpdA->Fill(eventInfo.YnEpdA);
      h_XnEpdB->Fill(eventInfo.XnEpdB);
      h_YnEpdB->Fill(eventInfo.YnEpdB);
      h_XnEpdC->Fill(eventInfo.XnEpdC);
      h_YnEpdC->Fill(eventInfo.YnEpdC);
      h_XnEpdD->Fill(eventInfo.XnEpdD);
      h_YnEpdD->Fill(eventInfo.YnEpdD);

      h_psiTpc_RAW->Fill(eventInfo.psiTpc);
      h_psiEpd_RAW->Fill(eventInfo.psiEpd);
      h_psiTpcA_RAW->Fill(eventInfo.psiTpcA);
      h_psiTpcB_RAW->Fill(eventInfo.psiTpcB);
      h_psiEpdA_RAW->Fill(eventInfo.psiEpdA);
      h_psiEpdB_RAW->Fill(eventInfo.psiEpdB);
      h_psiEpdC_RAW->Fill(eventInfo.psiEpdC);
      h_psiEpdD_RAW->Fill(eventInfo.psiEpdD);

      v_events.push_back(eventInfo);   // Store this event with all of its attributes
    }//End event loop

  eventInfo.reset();

  TH1D *h_xvtx = h2_trans_vtx->ProjectionX();
  TH1D *h_yvtx = h2_trans_vtx->ProjectionY();


  TH1D *h_XnTpc_RC;   // Re-centered histograms
  TH1D *h_YnTpc_RC;
  TH1D *h_XnTpcA_RC;
  TH1D *h_YnTpcA_RC;
  TH1D *h_XnTpcB_RC;
  TH1D *h_YnTpcB_RC;
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
  TH1D *h_psiEpd_RC;
  TH1D *h_psiEpdA_RC;
  TH1D *h_psiEpdB_RC;
  TH1D *h_psiEpdC_RC;
  TH1D *h_psiEpdD_RC;

  TH1D *h_psiTpc_FLAT;
  TH1D *h_psiTpcA_FLAT;
  TH1D *h_psiTpcB_FLAT;
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
      h_psiEpd_RC  = new TH1D("h_psiEpd_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdA_RC = new TH1D("h_psiEpdA_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdB_RC = new TH1D("h_psiEpdB_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC_RC = new TH1D("h_psiEpdC_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD C);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdD_RC = new TH1D("h_psiEpdD_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD D);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

      TH1D *h_XnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_XnTpc");
      TH1D *h_XnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcA");
      TH1D *h_XnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcB");
      TH1D *h_XnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_XnEpd");
      TH1D *h_XnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdA");
      TH1D *h_XnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdB");
      TH1D *h_XnEpdC_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdC");
      TH1D *h_XnEpdD_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdD");

      TH1D *h_YnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_YnTpc");
      TH1D *h_YnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcA");
      TH1D *h_YnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcB");
      TH1D *h_YnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_YnEpd");
      TH1D *h_YnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdA");
      TH1D *h_YnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdB");
      TH1D *h_YnEpdC_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdC");
      TH1D *h_YnEpdD_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdD");

      Double_t d_XnTpc_Avg  = h_XnTpc_INPUT->GetMean();
      Double_t d_XnTpcA_Avg = h_XnTpcA_INPUT->GetMean();
      Double_t d_XnTpcB_Avg = h_XnTpcB_INPUT->GetMean();
      Double_t d_XnEpd_Avg  = h_XnEpd_INPUT->GetMean();
      Double_t d_XnEpdA_Avg = h_XnEpdA_INPUT->GetMean();
      Double_t d_XnEpdB_Avg = h_XnEpdB_INPUT->GetMean();
      Double_t d_XnEpdC_Avg = h_XnEpdC_INPUT->GetMean();
      Double_t d_XnEpdD_Avg = h_XnEpdD_INPUT->GetMean();

      Double_t d_YnTpc_Avg  = h_YnTpc_INPUT->GetMean();
      Double_t d_YnTpcA_Avg = h_YnTpcA_INPUT->GetMean();
      Double_t d_YnTpcB_Avg = h_YnTpcB_INPUT->GetMean();
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
	  v_events.at(i).XnEpd  -= d_XnEpd_Avg;
	  v_events.at(i).XnEpdA -= d_XnEpdA_Avg;
	  v_events.at(i).XnEpdB -= d_XnEpdB_Avg;
	  v_events.at(i).XnEpdC -= d_XnEpdC_Avg;
	  v_events.at(i).XnEpdD -= d_XnEpdD_Avg;

	  v_events.at(i).YnTpc  -= d_YnTpc_Avg;
	  v_events.at(i).YnTpcA -= d_YnTpcA_Avg;
	  v_events.at(i).YnTpcB -= d_YnTpcB_Avg;
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
	  h_XnEpd_RC->Fill(v_events.at(i).XnEpd);
	  h_XnEpdA_RC->Fill(v_events.at(i).XnEpdA);
	  h_XnEpdB_RC->Fill(v_events.at(i).XnEpdB);
	  h_XnEpdC_RC->Fill(v_events.at(i).XnEpdC);
	  h_XnEpdD_RC->Fill(v_events.at(i).XnEpdD);

	  h_YnTpc_RC->Fill(v_events.at(i).YnTpc);
	  h_YnTpcA_RC->Fill(v_events.at(i).YnTpcA);
	  h_YnTpcB_RC->Fill(v_events.at(i).YnTpcB);
	  h_YnEpd_RC->Fill(v_events.at(i).YnEpd);
	  h_YnEpdA_RC->Fill(v_events.at(i).YnEpdA);
	  h_YnEpdB_RC->Fill(v_events.at(i).YnEpdB);
	  h_YnEpdC_RC->Fill(v_events.at(i).YnEpdC);
	  h_YnEpdD_RC->Fill(v_events.at(i).YnEpdD);

	  // Recalculate the event plane angles after re-centering	  
	  v_events.at(i).psiTpc  = TMath::ATan2(v_events.at(i).YnTpc,  v_events.at(i).XnTpc)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiTpcA = TMath::ATan2(v_events.at(i).YnTpcA, v_events.at(i).XnTpcA) / (Double_t)ORDER_N; 
	  v_events.at(i).psiTpcB = TMath::ATan2(v_events.at(i).YnTpcB, v_events.at(i).XnTpcB) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpd  = TMath::ATan2(v_events.at(i).YnEpd,  v_events.at(i).XnEpd)  / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdA = TMath::ATan2(v_events.at(i).YnEpdA, v_events.at(i).XnEpdA) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdB = TMath::ATan2(v_events.at(i).YnEpdB, v_events.at(i).XnEpdB) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdC = TMath::ATan2(v_events.at(i).YnEpdC, v_events.at(i).XnEpdC) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdD = TMath::ATan2(v_events.at(i).YnEpdD, v_events.at(i).XnEpdD) / (Double_t)ORDER_N; 

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_N);  // Maintain the correct periodicity of angles
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_N);
	  v_events.at(i).psiEpdA = angleShift(v_events.at(i).psiEpdA, ORDER_N);
	  v_events.at(i).psiEpdB = angleShift(v_events.at(i).psiEpdB, ORDER_N);
	  v_events.at(i).psiEpdC = angleShift(v_events.at(i).psiEpdC, ORDER_N);
	  v_events.at(i).psiEpdD = angleShift(v_events.at(i).psiEpdD, ORDER_N);

	  h_psiTpc_RC->Fill(v_events.at(i).psiTpc);
	  h_psiTpcA_RC->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_RC->Fill(v_events.at(i).psiTpcB);
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
	      jthSinAvg_Tpc = p_sinAvgsTpc_INPUT->GetBinContent(j);
	      jthCosAvg_Tpc = p_cosAvgsTpc_INPUT->GetBinContent(j);
	      jthSinAvg_TpcA = p_sinAvgsTpcA_INPUT->GetBinContent(j);
	      jthCosAvg_TpcA = p_cosAvgsTpcA_INPUT->GetBinContent(j);
	      jthSinAvg_TpcB = p_sinAvgsTpcB_INPUT->GetBinContent(j);
	      jthCosAvg_TpcB = p_cosAvgsTpcB_INPUT->GetBinContent(j);
	      jthSinAvg_Epd = p_sinAvgsEpd_INPUT->GetBinContent(j);
	      jthCosAvg_Epd = p_cosAvgsEpd_INPUT->GetBinContent(j);
	      jthSinAvg_EpdA = p_sinAvgsEpdA_INPUT->GetBinContent(j);
	      jthCosAvg_EpdA = p_cosAvgsEpdA_INPUT->GetBinContent(j);
	      jthSinAvg_EpdB = p_sinAvgsEpdB_INPUT->GetBinContent(j);
	      jthCosAvg_EpdB = p_cosAvgsEpdB_INPUT->GetBinContent(j);
	      jthSinAvg_EpdC = p_sinAvgsEpdC_INPUT->GetBinContent(j);
	      jthCosAvg_EpdC = p_cosAvgsEpdC_INPUT->GetBinContent(j);
	      jthSinAvg_EpdD = p_sinAvgsEpdD_INPUT->GetBinContent(j);
	      jthCosAvg_EpdD = p_cosAvgsEpdD_INPUT->GetBinContent(j);

	      psiTpc_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_Tpc*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc) 
								       +jthCosAvg_Tpc*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpc));
	      psiTpcA_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcA*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA) 
									+jthCosAvg_TpcA*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      psiTpcB_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcB*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB) 
									+jthCosAvg_TpcB*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      psiEpd_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_Epd*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpd)
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
	  v_events.at(i).psiEpd  += psiEpd_delta;
	  v_events.at(i).psiEpdA += psiEpdA_delta;
	  v_events.at(i).psiEpdB += psiEpdB_delta;
	  v_events.at(i).psiEpdC += psiEpdC_delta;
	  v_events.at(i).psiEpdD += psiEpdD_delta;

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_N);
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_N);
	  v_events.at(i).psiEpdA = angleShift(v_events.at(i).psiEpdA, ORDER_N);
	  v_events.at(i).psiEpdB = angleShift(v_events.at(i).psiEpdB, ORDER_N);
	  v_events.at(i).psiEpdC = angleShift(v_events.at(i).psiEpdC, ORDER_N);
	  v_events.at(i).psiEpdD = angleShift(v_events.at(i).psiEpdD, ORDER_N);

	  h_psiTpc_FLAT->Fill(v_events.at(i).psiTpc);
	  h_psiTpcA_FLAT->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_FLAT->Fill(v_events.at(i).psiTpcB);
	  h_psiEpd_FLAT->Fill(v_events.at(i).psiEpd);
	  h_psiEpdA_FLAT->Fill(v_events.at(i).psiEpdA);
	  h_psiEpdB_FLAT->Fill(v_events.at(i).psiEpdB);
	  h_psiEpdC_FLAT->Fill(v_events.at(i).psiEpdC);
	  h_psiEpdD_FLAT->Fill(v_events.at(i).psiEpdD);


	  // 2D Correlations between angles
	  h2_psiEpdAB->Fill(v_events.at(i).psiEpdB,v_events.at(i).psiEpdA);
	  h2_psiEpdAC->Fill(v_events.at(i).psiEpdC,v_events.at(i).psiEpdA);
	  h2_psiEpdAD->Fill(v_events.at(i).psiEpdD,v_events.at(i).psiEpdA);
	  h2_psiEpdBC->Fill(v_events.at(i).psiEpdC,v_events.at(i).psiEpdB);
	  h2_psiEpdBD->Fill(v_events.at(i).psiEpdD,v_events.at(i).psiEpdB);
	  h2_psiEpdCD->Fill(v_events.at(i).psiEpdD,v_events.at(i).psiEpdC);

	  h2_psiEpdATpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdA);
	  h2_psiEpdBTpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdB);
	  h2_psiEpdCTpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdC);
	  h2_psiEpdDTpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdD);

	  h2_psiEpdATpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdA);
	  h2_psiEpdBTpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdB);
	  h2_psiEpdCTpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdC);
	  h2_psiEpdDTpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdD);

	  h2_psiTpcATpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiTpcA);


	  Int_t tpcHitsA = v_events.at(i).phiValuesTpcA.size();
	  Int_t tpcHitsB = v_events.at(i).phiValuesTpcB.size();
	  Int_t epdHitsA = v_events.at(i).phiValuesEpdA.size();
	  Int_t epdHitsB = v_events.at(i).phiValuesEpdB.size();
	  Int_t epdHitsC = v_events.at(i).phiValuesEpdC.size();
	  Int_t epdHitsD = v_events.at(i).phiValuesEpdD.size();
	  Double_t phiTpc;
	  Double_t etaTpc;
	  Double_t phiEpd;
	  Double_t etaEpd;
	  Double_t psiTpc  = v_events.at(i).psiTpc;
	  Double_t psiTpcB = v_events.at(i).psiTpcB;
	  Double_t psiEpd  = v_events.at(i).psiEpd;
	  Int_t centralityID = v_events.at(i).centID;

	  for (int k = 0; k < epdHitsA; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdA.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdA.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	    }
	  for (int k = 0; k < epdHitsB; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdB.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdB.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	    }
	  for (int k = 0; k < epdHitsC; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdC.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdC.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	    }
	  for (int k = 0; k < epdHitsD; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdD.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdD.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	    }

	  for (int k = 0; k < tpcHitsA; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcA.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcA.at(k);

	      h2_v2SearchTpc->Fill(etaTpc, centralityID, TMath::Cos((Double_t)ORDER_N * (phiTpc - psiEpd)));
	    }
	  for (int k = 0; k < tpcHitsB; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcB.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcB.at(k);

	      h2_v2SearchTpc->Fill(etaTpc, centralityID, TMath::Cos((Double_t)ORDER_N * (phiTpc - psiEpd)));
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
  Double_t psiEpdD = 0;
  Double_t psiTpc  = 0;

  // CORRELATIONS AND FLOW COEFFICIENTS
  for (Int_t i = 0; i < numOfEvents; i++)
    {
      // SUBEVENT CORRELATIONS
      psiTpc  = v_events.at(i).psiTpc;
      psiEpdA = v_events.at(i).psiEpdA;
      psiEpdB = v_events.at(i).psiEpdB;
      psiEpdC = v_events.at(i).psiEpdC;
      psiEpdD = v_events.at(i).psiEpdD;

      AB_n2   = TMath::Cos(2.0 * (psiEpdA - psiEpdB));
      AC_n2   = TMath::Cos(2.0 * (psiEpdA - psiEpdC));
      AD_n2   = TMath::Cos(2.0 * (psiEpdA - psiEpdD));
      BC_n2   = TMath::Cos(2.0 * (psiEpdB - psiEpdC));
      BD_n2   = TMath::Cos(2.0 * (psiEpdB - psiEpdD));
      CD_n2   = TMath::Cos(2.0 * (psiEpdC - psiEpdD));
      ATpc_n2 = TMath::Cos(2.0 * (psiEpdA - psiTpc));
      BTpc_n2 = TMath::Cos(2.0 * (psiEpdB - psiTpc));
      CTpc_n2 = TMath::Cos(2.0 * (psiEpdC - psiTpc));
      DTpc_n2 = TMath::Cos(2.0 * (psiEpdD - psiTpc));
      
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
