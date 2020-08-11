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
const Double_t MIN_ETA_CUT = -6;//-5.16;    // Cuts for the EPD subevents (A, B, C, D)
const Double_t AB_ETA_CUT  = -4.05;//-3.82;
const Double_t BC_ETA_CUT  = -3.3;//-3.28;
const Double_t CD_ETA_CUT  = -2.9;//-2.87;
const Double_t MAX_ETA_CUT = -2.30;
*/
const Double_t C1_MIN_ETA_CUT = -4.0;       // Possble choices for sub-event C in the EPD
const Double_t C1_MAX_ETA_CUT = -3.0;
const Double_t C2_MIN_ETA_CUT = -4.0;
const Double_t C2_MAX_ETA_CUT = -2.3;
const Double_t C3_MIN_ETA_CUT = -4.4;
const Double_t C3_MAX_ETA_CUT = -2.3;

const Double_t R_VTX_CUT = 2.0;         // 2D r value, good vertices are within this value
const Double_t Z_VTX_CUT_LOW  = 199.5;
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

  Int_t nTracksTpcA;      // Number of GOOD tracks in the sub-event
  Double_t XnTpcA;
  Double_t YnTpcA;
  Double_t psiTpcA;       // Overall EP angle without removing autocorrelations
  std::vector<Double_t> phiValuesTpcA;   // Azimuthal angles for all TPC A particles in the event

  Int_t nTracksTpcB;
  Double_t XnTpcB;
  Double_t YnTpcB;
  Double_t psiTpcB;
  std::vector<Double_t> phiValuesTpcB;

  Int_t nhitsEpdC1;
  Double_t XnEpdC1;
  Double_t YnEpdC1;
  Double_t psiEpdC1;       // Subevent A EP angle without removing autocorrelations
  std::vector<Double_t> phiValuesEpdC1;  // Subevent A azimuthal angles

  Int_t nhitsEpdC2;
  Double_t XnEpdC2;
  Double_t YnEpdC2;
  Double_t psiEpdC2;       // Subevent B EP angle without removing autocorrelations  
  std::vector<Double_t> phiValuesEpdC2;  // Subevent B azimuthal angles

  Int_t nhitsEpdC3;
  Double_t XnEpdC3;
  Double_t YnEpdC3;
  Double_t psiEpdC3;       // Subevent D EP angle without removing autocorrelations  
  std::vector<Double_t> phiValuesEpdC3;  // Subevent D azimuthal angles

  void reset()
  {
    badEvent  = false;  //Reset all values in the struct to reuse
    primTracks = 0;
    centID = BAD_VALUE;

    nTracksTpcA = 0;
    XnTpcA = 0;
    YnTpcA = 0;
    psiTpcA = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcA);

    nTracksTpcB = 0;
    XnTpcB = 0;
    YnTpcB = 0;
    psiTpcB = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesTpcB);

    nhitsEpdC1 = 0;
    XnEpdC1 = 0;
    YnEpdC1 = 0;
    psiEpdC1 = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdC1);

    nhitsEpdC2 = 0;
    XnEpdC2 = 0;
    YnEpdC2 = 0;
    psiEpdC2 = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdC2);

    nhitsEpdC3 = 0;
    XnEpdC3 = 0;
    YnEpdC3 = 0;
    psiEpdC3 = BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdC3);
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
  if (event.XnTpcA == 0 && event.YnTpcA == 0) { event.badEvent = true; }
  else if (event.XnTpcB == 0 && event.YnTpcB == 0) { event.badEvent = true; }
  else if (event.XnEpdC1 == 0 && event.YnEpdC1 == 0) { event.badEvent = true; }
  else if (event.XnEpdC2 == 0 && event.YnEpdC2 == 0) { event.badEvent = true; }
  else if (event.XnEpdC3 == 0 && event.YnEpdC3 == 0) { event.badEvent = true; }
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
  TH1D *h_etaEpdC1 = new TH1D("h_etaEpdC1", "Particle #eta (EPD A);#eta;Particles", 600, -6, 0);
  TH1D *h_etaEpdC2 = new TH1D("h_etaEpdC2", "Particle #eta (EPD B);#eta;Particles", 600, -6, 0);
  TH1D *h_etaEpdC = new TH1D("h_etaEpdC", "Particle #eta (EPD C);#eta;Particles", 600, -6, 0);
  TH1D *h_etaEpdC3 = new TH1D("h_etaEpdC3", "Particle #eta (EPD D);#eta;Particles", 600, -6, 0);
  */
  TH1D *h_tileWeights = new TH1D("h_tileWeights", "EPD Tile Weights;Hits;nMIP Weights", 5, -1, 4);
  TH1D *h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TH1D *h_psiTpcA_RAW = new TH1D("h_psiTpcA_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_RAW = new TH1D("h_psiTpcB_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", TPC B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdC1_RAW = new TH1D("h_psiEpdC1_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdC2_RAW = new TH1D("h_psiEpdC2_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdC3_RAW = new TH1D("h_psiEpdC3_RAW", "Raw Event Plane Angles (n = "+ORDER_N_STR+", EPD D);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  //TH1D *h_v2Plot = new TH1D("h_v2Plot", "Plot to Retrieve v_{2};cos(2(#phi - #psi_{"+ORDER_N_STR+"}));Particles", 150, -1.5, 1.5);

  TProfile *p_TpcAB = new TProfile("p_TpcAB","TPC A-B Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{TPC,B}_{"+ORDER_N_STR+"}))>",
				   CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  
  TProfile *p_TpcAEpdC1 = new TProfile("p_TpcAEpdC1","TPC A EPD C1 Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,C1}_{"+ORDER_N_STR+"}))>", 
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdC2 = new TProfile("p_TpcAEpdC2","TPC A EPD C2 Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,C2}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdC3 = new TProfile("p_TpcAEpdC3","TPC A EPD C3 Correlations;Centrality;<cos(2(#psi^{TPC,A}_{"+ORDER_N_STR+"}-#psi^{EPD,C3}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcBEpdC1 = new TProfile("p_TpcBEpdC1","TPC B EPD C1 Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,C1}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdC2 = new TProfile("p_TpcBEpdC2","TPC B EPD C2 Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,C2}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdC3 = new TProfile("p_TpcBEpdC3","TPC B EPD C3 Correlations;Centrality;<cos(2(#psi^{TPC,B}_{"+ORDER_N_STR+"}-#psi^{EPD,C3}_{"+ORDER_N_STR+"}))>",
				       CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

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
  /*
  TH2D *h2_phiSearchTpc = new TH2D("h2_phiSearchTpc", "Azimuthal Distribution by Centrality;#phi;Centrality (%)", 
				   200, -4, 4, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TH2D *h2_phiSearchEpd = new TH2D("h2_phiSearchEpd", "Azimuthal Distribution by Centrality;#phi;Centrality (%)", 
				   200, -4, 4, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

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
  */
  // The indices here are equivalent to the corresponding centrality ID
  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};

  // The indices here are opposite to the corresponding centrality ID (array is backward)
  //const char *centralityBins[16] = {"0-5", "5-10", "10-15" "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80"};


  TH2D *h2_psiEpdC1TpcA = new TH2D("h2_psiEpdC1TpcA", "#psi^{EPD,C1} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdC2TpcA = new TH2D("h2_psiEpdC2TpcA", "#psi^{EPD,C2} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdC3TpcA = new TH2D("h2_psiEpdC3TpcA", "#psi^{EPD,C3} vs #psi^{TPC,A} (Order "+ORDER_N_STR+");#psi^{TPC}_{A};#psi^{EPD}_{D}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiEpdC1TpcB = new TH2D("h2_psiEpdC1TpcB", "#psi^{EPD,C1} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdC2TpcB = new TH2D("h2_psiEpdC2TpcB", "#psi^{EPD,C2} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdC3TpcB = new TH2D("h2_psiEpdC3TpcB", "#psi^{EPD,C3} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{EPD}_{D}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  
  TH2D *h2_psiTpcATpcB = new TH2D("h2_psiTpcATpcB", "#psi^{TPC,A} vs #psi^{TPC,B} (Order "+ORDER_N_STR+");#psi^{TPC}_{B};#psi^{TPC}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);



  // CORRECTION HISTOGRAMS
  TProfile *p_sinAvgsTpcA = new TProfile("p_sinAvgsTpcA", "Sin Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpcA = new TProfile("p_cosAvgsTpcA", "Cos Averages (TPC A);j (Correction term);<sin(jn#psi^{TPC}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsTpcB = new TProfile("p_sinAvgsTpcB", "Sin Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsTpcB = new TProfile("p_cosAvgsTpcB", "Cos Averages (TPC B);j (Correction term);<sin(jn#psi^{TPC}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdC1 = new TProfile("p_sinAvgsEpdC1", "Sin Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdC1 = new TProfile("p_cosAvgsEpdC1", "Cos Averages (EPD A);j (Correction term);<sin(jn#psi^{EPD}_{A,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdC2 = new TProfile("p_sinAvgsEpdC2", "Sin Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdC2 = new TProfile("p_cosAvgsEpdC2", "Cos Averages (EPD B);j (Correction term);<sin(jn#psi^{EPD}_{B,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_sinAvgsEpdC3 = new TProfile("p_sinAvgsEpdC3", "Sin Averages (EPD D);j (Correction term);<sin(jn#psi^{EPD}_{D,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);
  TProfile *p_cosAvgsEpdC3 = new TProfile("p_cosAvgsEpdC3", "Cos Averages (EPD D);j (Correction term);<sin(jn#psi^{EPD}_{D,n})>", SHIFT_TERMS, 1, SHIFT_TERMS + 1);

  TH1D *h_XnTpcA = new TH1D("h_XnTpcA", "X_n Distribution (TPC A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcA = new TH1D("h_YnTpcA", "Y_n Distribution (TPC A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnTpcB = new TH1D("h_XnTpcB", "X_n Distribution (TPC B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnTpcB = new TH1D("h_YnTpcB", "Y_n Distribution (TPC B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdC1 = new TH1D("h_XnEpdC1", "X_n Distribution (EPD A);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdC1 = new TH1D("h_YnEpdC1", "Y_n Distribution (EPD A);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdC2 = new TH1D("h_XnEpdC2", "X_n Distribution (EPD B);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdC2 = new TH1D("h_YnEpdC2", "Y_n Distribution (EPD B);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_XnEpdC3 = new TH1D("h_XnEpdC3", "X_n Distribution (EPD D);X_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);
  TH1D *h_YnEpdC3 = new TH1D("h_YnEpdC3", "Y_n Distribution (EPD D);Y_n;Events", 250, -Q_BOUNDS, Q_BOUNDS);

  
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

	      h_eta_s->Fill(d_eta - Y_MID);

	      if (d_eta > MIN_TPC_ETA_CUT && d_eta < AGAP_TPC_ETA_CUT)
		{
		  eventInfo.phiValuesTpcA.push_back(d_phi);

		  eventInfo.nTracksTpcA++;
		  eventInfo.XnTpcA += d_pT * TMath::Cos((Double_t)ORDER_N * d_phi);
		  eventInfo.YnTpcA += d_pT * TMath::Sin((Double_t)ORDER_N * d_phi);
		}
	      else if (d_eta > GAPB_TPC_ETA_CUT && d_eta < MAX_TPC_ETA_CUT)
		{
		  eventInfo.phiValuesTpcB.push_back(d_phi);

		  eventInfo.nTracksTpcB++;
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

	  h_eta_s->Fill(tileEta - Y_MID);

	  if (tileEta > C1_MIN_ETA_CUT && tileEta < C1_MAX_ETA_CUT)
	    {
	      eventInfo.nhitsEpdC1++;
	      eventInfo.XnEpdC1 += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdC1 += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdC1.push_back(tilePhi);	      
	    }
	  else if (tileEta > C2_MIN_ETA_CUT && tileEta < C2_MAX_ETA_CUT)
	    {
	      eventInfo.nhitsEpdC2++;
	      eventInfo.XnEpdC2 += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdC2 += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdC2.push_back(tilePhi);
	    }
	  else if (tileEta > C3_MIN_ETA_CUT && tileEta < C3_MAX_ETA_CUT)
	    {
	      eventInfo.nhitsEpdC3++;
	      eventInfo.XnEpdC3 += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdC3 += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.phiValuesEpdC3.push_back(tilePhi);
	    }
	  
	  /*
	  if (tileEta > MIN_ETA_CUT && tileEta < AB_ETA_CUT)  // Sub A
	    {
	    }
	  else if (tileEta >= AB_ETA_CUT && tileEta < BC_ETA_CUT)  // Sub B
	    {
	    }
	  else if (tileEta >= BC_ETA_CUT && tileEta < CD_ETA_CUT)  // Sub C
	    {
	      eventInfo.nhitsEpdC++;
	      eventInfo.XnEpdC += tileWeight * TMath::Cos((Double_t)ORDER_N * tilePhi);
	      eventInfo.YnEpdC += tileWeight * TMath::Sin((Double_t)ORDER_N * tilePhi);
	      eventInfo.etaValuesEpdC.push_back(tileEta);
	      eventInfo.phiValuesEpdC.push_back(tilePhi);
	    }
	  else if (tileEta >= CD_ETA_CUT && tileEta < MAX_ETA_CUT)  // Sub D
	    {
	    }
	  */

	}
      delete epdGeom;

      //=========================================================
      //            END EPD STUFF
      //=========================================================



      if (eventInfo.nTracksTpcA < MIN_TRACKS) continue;
      if (eventInfo.nTracksTpcB < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdC1 < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdC2 < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdC3 < MIN_TRACKS) continue;
      if (eventInfo.XnTpcA == 0 && eventInfo.YnTpcA == 0) continue;
      if (eventInfo.XnTpcB == 0 && eventInfo.YnTpcB == 0) continue;
      if (eventInfo.XnEpdC1 == 0 && eventInfo.YnEpdC1 == 0) continue;
      if (eventInfo.XnEpdC2 == 0 && eventInfo.YnEpdC2 == 0) continue;
      if (eventInfo.XnEpdC3 == 0 && eventInfo.YnEpdC3 == 0) continue;


      // RAW SUB-EVENT PLANE ANGLES //
      if (ORDER_N % 2 == 1)                 // Q vectors of EPD East are opposite sign of West for odd harmonics. Switch back to compare to TPC.
	{
	  eventInfo.XnEpdC1 *= -1.0;
	  eventInfo.YnEpdC1 *= -1.0;
	  eventInfo.XnEpdC2 *= -1.0;
	  eventInfo.YnEpdC2 *= -1.0;
	  eventInfo.XnEpdC3 *= -1.0;
	  eventInfo.YnEpdC3 *= -1.0;
	}
      eventInfo.psiTpcA = TMath::ATan2(eventInfo.YnTpcA, eventInfo.XnTpcA) / (Double_t)ORDER_N;
      eventInfo.psiTpcB = TMath::ATan2(eventInfo.YnTpcB, eventInfo.XnTpcB) / (Double_t)ORDER_N;
      eventInfo.psiEpdC1 = TMath::ATan2(eventInfo.YnEpdC1, eventInfo.XnEpdC1) / (Double_t)ORDER_N;
      eventInfo.psiEpdC2 = TMath::ATan2(eventInfo.YnEpdC2, eventInfo.XnEpdC2) / (Double_t)ORDER_N;
      eventInfo.psiEpdC3 = TMath::ATan2(eventInfo.YnEpdC3, eventInfo.XnEpdC3) / (Double_t)ORDER_N;


      // Filling histos here since this is past all possible cuts
      /*
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcA.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcA.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcB.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesTpcB.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC1.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdC1.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC2.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdC2.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdC.at(i) - Y_MID ); }
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC3.size(); i++) { h_eta_s->Fill( eventInfo.etaValuesEpdC3.at(i) - Y_MID ); }
      */
      //for (unsigned int i = 0; i < eventInfo.etaValues.size();   i++)  { h_eta->Fill(eventInfo.etaValues.at(i)); }
      //for (unsigned int i = 0; i < eventInfo.pTValuesTpcA.size(); i++)  { h_pTA->Fill(eventInfo.pTValuesTpcA.at(i)); }
      //for (unsigned int i = 0; i < eventInfo.pTValuesTpcB.size(); i++)  { h_pTB->Fill(eventInfo.pTValuesTpcB.at(i)); }

      h_primTracks->Fill(eventInfo.primTracks);
      h_centralities->Fill(eventInfo.centID);

      h_XnTpcA->Fill(eventInfo.XnTpcA);
      h_YnTpcA->Fill(eventInfo.YnTpcA);
      h_XnTpcB->Fill(eventInfo.XnTpcB);
      h_YnTpcB->Fill(eventInfo.YnTpcB);
      h_XnEpdC1->Fill(eventInfo.XnEpdC1);
      h_YnEpdC1->Fill(eventInfo.YnEpdC1);
      h_XnEpdC2->Fill(eventInfo.XnEpdC2);
      h_YnEpdC2->Fill(eventInfo.YnEpdC2);
      h_XnEpdC3->Fill(eventInfo.XnEpdC3);
      h_YnEpdC3->Fill(eventInfo.YnEpdC3);

      h_psiTpcA_RAW->Fill(eventInfo.psiTpcA);
      h_psiTpcB_RAW->Fill(eventInfo.psiTpcB);
      h_psiEpdC1_RAW->Fill(eventInfo.psiEpdC1);
      h_psiEpdC2_RAW->Fill(eventInfo.psiEpdC2);
      h_psiEpdC3_RAW->Fill(eventInfo.psiEpdC3);

      v_events.push_back(eventInfo);   // Store this event with all of its attributes
    }//End event loop

  eventInfo.reset();

  //TH1D *h_xvtx = h2_trans_vtx->ProjectionX();
  //TH1D *h_yvtx = h2_trans_vtx->ProjectionY();


  TH1D *h_XnTpcA_RC;   // Re-centered histograms
  TH1D *h_YnTpcA_RC;
  TH1D *h_XnTpcB_RC;
  TH1D *h_YnTpcB_RC;
  TH1D *h_XnEpdC1_RC;
  TH1D *h_YnEpdC1_RC;
  TH1D *h_XnEpdC2_RC;
  TH1D *h_YnEpdC2_RC;
  TH1D *h_XnEpdC3_RC;
  TH1D *h_YnEpdC3_RC;

  TH1D *h_psiTpcA_RC;
  TH1D *h_psiTpcB_RC;
  TH1D *h_psiEpdC1_RC;
  TH1D *h_psiEpdC2_RC;
  TH1D *h_psiEpdC3_RC;

  TH1D *h_psiTpcA_FLAT;
  TH1D *h_psiTpcB_FLAT;
  TH1D *h_psiEpdC1_FLAT;
  TH1D *h_psiEpdC2_FLAT;
  TH1D *h_psiEpdC3_FLAT;


  //=========================================================
  //          Re-centering (Xn, Yn) Distributions
  //=========================================================

  if (RUN_ITERATION == 1 || RUN_ITERATION == 2)
    {
      std::cout << "Re-centering flow vectors and accumulating sin/cos averages..." << std::endl;

      h_XnTpcA_RC = new TH1D("h_XnTpcA_RC", "Re-centered X_n Distribution (TPC A);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnTpcA_RC = new TH1D("h_YnTpcA_RC", "Re-centered Y_n Distribution (TPC A);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnTpcB_RC = new TH1D("h_XnTpcB_RC", "Re-centered X_n Distribution (TPC B);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnTpcB_RC = new TH1D("h_YnTpcB_RC", "Re-centered Y_n Distribution (TPC B);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdC1_RC = new TH1D("h_XnEpdC1_RC", "Re-centered X_n Distribution (EPD A);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdC1_RC = new TH1D("h_YnEpdC1_RC", "Re-centered Y_n Distribution (EPD A);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdC2_RC = new TH1D("h_XnEpdC2_RC", "Re-centered X_n Distribution (EPD B);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdC2_RC = new TH1D("h_YnEpdC2_RC", "Re-centered Y_n Distribution (EPD B);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_XnEpdC3_RC = new TH1D("h_XnEpdC3_RC", "Re-centered X_n Distribution (EPD D);X_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);
      h_YnEpdC3_RC = new TH1D("h_YnEpdC3_RC", "Re-centered Y_n Distribution (EPD D);Y_n;Events", 200, -Q_BOUNDS, Q_BOUNDS);

      h_psiTpcA_RC = new TH1D("h_psiTpcA_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiTpcB_RC = new TH1D("h_psiTpcB_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", TPC B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC1_RC = new TH1D("h_psiEpdC1_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD A);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC2_RC = new TH1D("h_psiEpdC2_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD B);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC3_RC = new TH1D("h_psiEpdC3_RC", "Re-centered Event Plane Angles (n = "+ORDER_N_STR+", EPD D);#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

      TH1D *h_XnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcA");
      TH1D *h_XnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcB");
      TH1D *h_XnEpdC1_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdC1");
      TH1D *h_XnEpdC2_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdC2");
      TH1D *h_XnEpdC3_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdC3");

      TH1D *h_YnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcA");
      TH1D *h_YnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcB");
      TH1D *h_YnEpdC1_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdC1");
      TH1D *h_YnEpdC2_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdC2");
      TH1D *h_YnEpdC3_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdC3");

      Double_t d_XnTpcA_Avg = h_XnTpcA_INPUT->GetMean();
      Double_t d_XnTpcB_Avg = h_XnTpcB_INPUT->GetMean();
      Double_t d_XnEpdC1_Avg = h_XnEpdC1_INPUT->GetMean();
      Double_t d_XnEpdC2_Avg = h_XnEpdC2_INPUT->GetMean();
      Double_t d_XnEpdC3_Avg = h_XnEpdC3_INPUT->GetMean();

      Double_t d_YnTpcA_Avg = h_YnTpcA_INPUT->GetMean();
      Double_t d_YnTpcB_Avg = h_YnTpcB_INPUT->GetMean();
      Double_t d_YnEpdC1_Avg = h_YnEpdC1_INPUT->GetMean();
      Double_t d_YnEpdC2_Avg = h_YnEpdC2_INPUT->GetMean();
      Double_t d_YnEpdC3_Avg = h_YnEpdC3_INPUT->GetMean();


      Int_t numOfEvents = v_events.size();
      Int_t badEvents = 0;

      for (int i = 0; i < numOfEvents; i++)
	{
	  v_events.at(i).XnTpcA -= d_XnTpcA_Avg;
	  v_events.at(i).XnTpcB -= d_XnTpcB_Avg;
	  v_events.at(i).XnEpdC1 -= d_XnEpdC1_Avg;
	  v_events.at(i).XnEpdC2 -= d_XnEpdC2_Avg;
	  v_events.at(i).XnEpdC3 -= d_XnEpdC3_Avg;

	  v_events.at(i).YnTpcA -= d_YnTpcA_Avg;
	  v_events.at(i).YnTpcB -= d_YnTpcB_Avg;
	  v_events.at(i).YnEpdC1 -= d_YnEpdC1_Avg;
	  v_events.at(i).YnEpdC2 -= d_YnEpdC2_Avg;
	  v_events.at(i).YnEpdC3 -= d_YnEpdC3_Avg;

	  checkZeroQ(v_events.at(i));

	  if ( v_events.at(i).badEvent ) { badEvents++; continue; }

	  h_XnTpcA_RC->Fill(v_events.at(i).XnTpcA);
	  h_XnTpcB_RC->Fill(v_events.at(i).XnTpcB);
	  h_XnEpdC1_RC->Fill(v_events.at(i).XnEpdC1);
	  h_XnEpdC2_RC->Fill(v_events.at(i).XnEpdC2);
	  h_XnEpdC3_RC->Fill(v_events.at(i).XnEpdC3);

	  h_YnTpcA_RC->Fill(v_events.at(i).YnTpcA);
	  h_YnTpcB_RC->Fill(v_events.at(i).YnTpcB);
	  h_YnEpdC1_RC->Fill(v_events.at(i).YnEpdC1);
	  h_YnEpdC2_RC->Fill(v_events.at(i).YnEpdC2);
	  h_YnEpdC3_RC->Fill(v_events.at(i).YnEpdC3);

	  // Recalculate the event plane angles after re-centering	  
	  v_events.at(i).psiTpcA = TMath::ATan2(v_events.at(i).YnTpcA, v_events.at(i).XnTpcA) / (Double_t)ORDER_N; 
	  v_events.at(i).psiTpcB = TMath::ATan2(v_events.at(i).YnTpcB, v_events.at(i).XnTpcB) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdC1 = TMath::ATan2(v_events.at(i).YnEpdC1, v_events.at(i).XnEpdC1) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdC2 = TMath::ATan2(v_events.at(i).YnEpdC2, v_events.at(i).XnEpdC2) / (Double_t)ORDER_N; 
	  v_events.at(i).psiEpdC3 = TMath::ATan2(v_events.at(i).YnEpdC3, v_events.at(i).XnEpdC3) / (Double_t)ORDER_N; 

	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiEpdC1 = angleShift(v_events.at(i).psiEpdC1, ORDER_N);
	  v_events.at(i).psiEpdC2 = angleShift(v_events.at(i).psiEpdC2, ORDER_N);
	  v_events.at(i).psiEpdC3 = angleShift(v_events.at(i).psiEpdC3, ORDER_N);

	  h_psiTpcA_RC->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_RC->Fill(v_events.at(i).psiTpcB);
	  h_psiEpdC1_RC->Fill(v_events.at(i).psiEpdC1);
	  h_psiEpdC2_RC->Fill(v_events.at(i).psiEpdC2);
	  h_psiEpdC3_RC->Fill(v_events.at(i).psiEpdC3);


	  // Accumulate terms for averages over the re-centered angles for event plane angle shifting
	  for (int j = 1; j <= SHIFT_TERMS; j++)
	    {
	      p_sinAvgsTpcA->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      p_cosAvgsTpcA->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      p_sinAvgsTpcB->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      p_cosAvgsTpcB->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      p_sinAvgsEpdC1->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC1));
	      p_cosAvgsEpdC1->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC1));
	      p_sinAvgsEpdC2->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC2));
	      p_cosAvgsEpdC2->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC2));
	      p_sinAvgsEpdC3->Fill(j, TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC3));
	      p_cosAvgsEpdC3->Fill(j, TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC3));
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

      h_psiTpcA_FLAT = new TH1D("h_psiTpcA_FLAT", "Flattened Event Plane Angle (TPC A, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiTpcB_FLAT = new TH1D("h_psiTpcB_FLAT", "Flattened Event Plane Angle (TPC B, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiEpdC1_FLAT = new TH1D("h_psiEpdC1_FLAT", "Flattened Event Plane Angle (EPD A, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC2_FLAT = new TH1D("h_psiEpdC2_FLAT", "Flattened Event Plane Angle (EPD B, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC3_FLAT = new TH1D("h_psiEpdC3_FLAT", "Flattened Event Plane Angle (EPD D, order "+ORDER_N_STR+");#psi_{"+ORDER_N_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

      TProfile *p_sinAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcA");
      TProfile *p_cosAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcA");
      TProfile *p_sinAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcB");
      TProfile *p_cosAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcB");
      TProfile *p_sinAvgsEpdC1_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdC1");
      TProfile *p_cosAvgsEpdC1_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdC1");
      TProfile *p_sinAvgsEpdC2_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdC2");
      TProfile *p_cosAvgsEpdC2_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdC2");
      TProfile *p_sinAvgsEpdC3_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdC3");
      TProfile *p_cosAvgsEpdC3_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdC3");


      // Get corrected event plane angles //
      for (Int_t i = 0; i < numOfEvents; i++)  // Loop over v_events to correct their angles
	{
	  if ( v_events.at(i).badEvent == true) { continue; }

	  Double_t psiTpcA_delta = 0;
	  Double_t psiTpcB_delta = 0;
	  Double_t psiEpdC1_delta = 0;
	  Double_t psiEpdC2_delta = 0;
	  Double_t psiEpdC3_delta = 0;

	  Double_t jthSinAvg_TpcA = 0;
	  Double_t jthCosAvg_TpcA = 0;
	  Double_t jthSinAvg_TpcB = 0;
	  Double_t jthCosAvg_TpcB = 0;
	  Double_t jthSinAvg_EpdC1 = 0;
	  Double_t jthCosAvg_EpdC1 = 0;
	  Double_t jthSinAvg_EpdC2 = 0;
	  Double_t jthCosAvg_EpdC2 = 0;
	  Double_t jthSinAvg_EpdC3 = 0;
	  Double_t jthCosAvg_EpdC3 = 0;


	  for (Int_t j = 1; j <= SHIFT_TERMS; j++)    // Build the correction sums
	    {
	      jthSinAvg_TpcA = p_sinAvgsTpcA_INPUT->GetBinContent(j);
	      jthCosAvg_TpcA = p_cosAvgsTpcA_INPUT->GetBinContent(j);
	      jthSinAvg_TpcB = p_sinAvgsTpcB_INPUT->GetBinContent(j);
	      jthCosAvg_TpcB = p_cosAvgsTpcB_INPUT->GetBinContent(j);
	      jthSinAvg_EpdC1 = p_sinAvgsEpdC1_INPUT->GetBinContent(j);
	      jthCosAvg_EpdC1 = p_cosAvgsEpdC1_INPUT->GetBinContent(j);
	      jthSinAvg_EpdC2 = p_sinAvgsEpdC2_INPUT->GetBinContent(j);
	      jthCosAvg_EpdC2 = p_cosAvgsEpdC2_INPUT->GetBinContent(j);
	      jthSinAvg_EpdC3 = p_sinAvgsEpdC3_INPUT->GetBinContent(j);
	      jthCosAvg_EpdC3 = p_cosAvgsEpdC3_INPUT->GetBinContent(j);

	      psiTpcA_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcA*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA) 
									+jthCosAvg_TpcA*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcA));
	      psiTpcB_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_TpcB*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB) 
									+jthCosAvg_TpcB*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiTpcB));
	      psiEpdC1_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdC1*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC1)
									+jthCosAvg_EpdC1*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC1));
	      psiEpdC2_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdC2*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC2)
									+jthCosAvg_EpdC2*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC2));
	      psiEpdC3_delta += (2.0/((Double_t)j*(Double_t)ORDER_N)) * (-jthSinAvg_EpdC3*TMath::Cos((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC3)
									+jthCosAvg_EpdC3*TMath::Sin((Double_t)j * (Double_t)ORDER_N * v_events.at(i).psiEpdC3));
	    }


	  v_events.at(i).psiTpcA += psiTpcA_delta;
	  v_events.at(i).psiTpcB += psiTpcB_delta;
	  v_events.at(i).psiEpdC1 += psiEpdC1_delta;
	  v_events.at(i).psiEpdC2 += psiEpdC2_delta;
	  v_events.at(i).psiEpdC3 += psiEpdC3_delta;

	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_N);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_N);
	  v_events.at(i).psiEpdC1 = angleShift(v_events.at(i).psiEpdC1, ORDER_N);
	  v_events.at(i).psiEpdC2 = angleShift(v_events.at(i).psiEpdC2, ORDER_N);
	  v_events.at(i).psiEpdC3 = angleShift(v_events.at(i).psiEpdC3, ORDER_N);

	  h_psiTpcA_FLAT->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_FLAT->Fill(v_events.at(i).psiTpcB);
	  h_psiEpdC1_FLAT->Fill(v_events.at(i).psiEpdC1);
	  h_psiEpdC2_FLAT->Fill(v_events.at(i).psiEpdC2);
	  h_psiEpdC3_FLAT->Fill(v_events.at(i).psiEpdC3);


	  // 2D Correlations between angles
	  h2_psiEpdC1TpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdC1);
	  h2_psiEpdC2TpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdC2);
	  h2_psiEpdC3TpcA->Fill(v_events.at(i).psiTpcA,v_events.at(i).psiEpdC3);

	  h2_psiEpdC1TpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdC1);
	  h2_psiEpdC2TpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdC2);
	  h2_psiEpdC3TpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiEpdC3);

	  h2_psiTpcATpcB->Fill(v_events.at(i).psiTpcB,v_events.at(i).psiTpcA);
	  //


	  // 1D correlation averages used in calculating resolution using the 3 sub-event method
	  p_TpcAB->Fill(v_events.at(i).centID,     TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiTpcB)));
	  p_TpcAEpdC1->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdC1)));
	  p_TpcAEpdC2->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdC2)));
	  p_TpcAEpdC3->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdC3)));
	  p_TpcBEpdC1->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdC1)));
	  p_TpcBEpdC2->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdC2)));
	  p_TpcBEpdC3->Fill(v_events.at(i).centID, TMath::Cos(2.0 * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdC3)));
	  //


	  // 2D searches through eta and centrality for correlations between detectors
	  /*
	  Int_t tpcHitsA  = v_events.at(i).phiValuesTpcA.size();
	  Int_t tpcHitsB  = v_events.at(i).phiValuesTpcB.size();
	  Int_t epdHitsC1 = v_events.at(i).phiValuesEpdC1.size();
	  Int_t epdHitsC2 = v_events.at(i).phiValuesEpdC2.size();
	  Int_t epdHitsC3 = v_events.at(i).phiValuesEpdC3.size();
	  Double_t phiTpc;
	  Double_t etaTpc;
	  Double_t phiEpd;
	  Double_t etaEpd;
	  Double_t psiTpc  = v_events.at(i).psiTpc;
	  Double_t psiTpcA = v_events.at(i).psiTpcA;
	  Double_t psiTpcB = v_events.at(i).psiTpcB;
	  //Double_t psiEpd  = v_events.at(i).psiEpd;
	  Int_t centralityID = v_events.at(i).centID;

	  for (int k = 0; k < epdHitsA; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdC1.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdC1.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsB; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdC2.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdC2.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsC; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdC.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdC.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsD; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdC3.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdC3.at(k);

	      h2_v2SearchEpd->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpc)));
	      h2_v2SearchEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcA)));
	      h2_v2SearchEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos((Double_t)ORDER_N * (phiEpd - psiTpcB)));
	      h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }

	  for (int k = 0; k < tpcHitsA; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcA.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcA.at(k);

	      h2_v2SearchTpc->Fill(etaTpc, centralityID, TMath::Cos((Double_t)ORDER_N * (phiTpc - psiEpd)));
	      h2_phiSearchTpc->Fill(phiTpc, centralityID);
	    }
	  for (int k = 0; k < tpcHitsB; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcB.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcB.at(k);

	      h2_v2SearchTpc->Fill(etaTpc, centralityID, TMath::Cos((Double_t)ORDER_N * (phiTpc - psiEpd)));
	      h2_phiSearchTpc->Fill(phiTpc, centralityID);
	    }
	  */
	  v_events.at(i).reset(); // Try to free up space?

	}// End shift loop over events

    }
  //=========================================================
  //          End Event Plane Angle Shifting
  //=========================================================

  // Switch y-axis labels to centrality percentages
  /*
  Int_t labelIndex;
  for (int i = 1; i <= CENT_BINS; i++) 
    {
      labelIndex = FIRST_CENT + i - 1;
      h2_v2SearchTpc->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      h2_v2SearchEpd->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_v2SearchEpdTpcA->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_v2SearchEpdTpcB->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_phiSearchTpc->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      h2_phiSearchEpd->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
    }
  */

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

  Double_t psiEpdC1 = 0;
  Double_t psiEpdC2 = 0;
  Double_t psiEpdC = 0;
  Double_t psiEpdC3 = 0;
  Double_t psiTpc  = 0;

  // CORRELATIONS AND FLOW COEFFICIENTS
  for (Int_t i = 0; i < numOfEvents; i++)
    {
      // SUBEVENT CORRELATIONS
      psiTpc  = v_events.at(i).psiTpc;
      psiEpdC1 = v_events.at(i).psiEpdC1;
      psiEpdC2 = v_events.at(i).psiEpdC2;
      psiEpdC = v_events.at(i).psiEpdC;
      psiEpdC3 = v_events.at(i).psiEpdC3;

      AB_n2   = TMath::Cos(2.0 * (psiEpdC1 - psiEpdC2));
      AC_n2   = TMath::Cos(2.0 * (psiEpdC1 - psiEpdC));
      AD_n2   = TMath::Cos(2.0 * (psiEpdC1 - psiEpdC3));
      BC_n2   = TMath::Cos(2.0 * (psiEpdC2 - psiEpdC));
      BD_n2   = TMath::Cos(2.0 * (psiEpdC2 - psiEpdC3));
      CD_n2   = TMath::Cos(2.0 * (psiEpdC - psiEpdC3));
      ATpc_n2 = TMath::Cos(2.0 * (psiEpdC1 - psiTpc));
      BTpc_n2 = TMath::Cos(2.0 * (psiEpdC2 - psiTpc));
      CTpc_n2 = TMath::Cos(2.0 * (psiEpdC - psiTpc));
      DTpc_n2 = TMath::Cos(2.0 * (psiEpdC3 - psiTpc));
      
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

      p_sinAvgsTpcA  ->Write();
      p_cosAvgsTpcA  ->Write();
      p_sinAvgsTpcB  ->Write();
      p_cosAvgsTpcB  ->Write();
      p_sinAvgsEpdC1  ->Write();
      p_cosAvgsEpdC1  ->Write();
      p_sinAvgsEpdC2  ->Write();
      p_cosAvgsEpdC2  ->Write();
      p_sinAvgsEpdC3  ->Write();
      p_cosAvgsEpdC3  ->Write();
      h_XnTpcA       ->Write();
      h_YnTpcA       ->Write();
      h_XnTpcB       ->Write();
      h_YnTpcB       ->Write();
      h_XnEpdC1       ->Write();
      h_YnEpdC1       ->Write();
      h_XnEpdC2       ->Write();
      h_YnEpdC2       ->Write();
      h_XnEpdC3       ->Write();
      h_YnEpdC3       ->Write();

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
