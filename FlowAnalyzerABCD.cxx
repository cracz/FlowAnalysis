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
const Double_t ORDER_N = 2.0;          // Order of anisotropic flow (v_n)
TString ORDER_N_STR;                   // ORDER_N but as a TString for titles/labels

const Double_t ORDER_M = 1.0;          // Order of event plane angle (psi_m)
TString ORDER_M_STR;                   // ORDER_M but as a TString for titles/labels

const Double_t PSI_BOUNDS = TMath::Pi()/ORDER_M + 1;
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
/*
const Double_t MIN_ETA_CUT = -5.1;
const Double_t EF_ETA_CUT  = -3.3;
const Double_t MAX_ETA_CUT = -2.4;
*/
const Double_t R_VTX_CUT = 2.0;         // 2D radius, good vertices are within this value
const Double_t Z_VTX_CUT_LOW  = 199.5;
const Double_t Z_VTX_CUT_HIGH = 201.5;

const Int_t MIN_TRACKS = 8;             // Min number of tracks/hits in each sub-event

const Double_t PI_KA_MOM_CUT = 1.6;     // Reject those above these total momenta (GeV/c)
const Double_t PR_MOM_CUT    = 3.0;

const Int_t SHIFT_TERMS = 10;           // Number of terms to use when shifting event plane angles

const Int_t I_BAD_VALUE    = -999;
const Double_t D_BAD_VALUE = -999.0;

const Int_t CENT_BINS  = 16;             // Number of centrality bins to show (max 16)  LEAVE AT 16 FOR NOW, BEST FOR RESOLUTION STUFF
const Int_t FIRST_CENT = 16 - CENT_BINS;            // Starting point for centrality dependent plots

const Double_t Y_MID = -1.05;       // Mid rapidity
/*
const Double_t Y_MID_PP = -0.83;       // Mid rapidity for pions (plus and minus), kaons, and protons
const Double_t Y_MID_PM = -0.83;
const Double_t Y_MID_KP = -0.73;
const Double_t Y_MID_KM = -0.80;
const Double_t Y_MID_PR = -0.61;
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

  Int_t nHitsEpd;
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

    nHitsEpd = 0;
    XnEpd = 0;
    YnEpd = 0;
    psiEpd = D_BAD_VALUE;

    nhitsEpdA = 0;
    XnEpdA = 0;
    YnEpdA = 0;
    psiEpdA = D_BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdA);
    std::vector<Double_t>().swap(etaValuesEpdA);

    nhitsEpdB = 0;
    XnEpdB = 0;
    YnEpdB = 0;
    psiEpdB = D_BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdB);
    std::vector<Double_t>().swap(etaValuesEpdB);

    nhitsEpdC = 0;
    XnEpdC = 0;
    YnEpdC = 0;
    psiEpdC = D_BAD_VALUE;
    std::vector<Double_t>().swap(phiValuesEpdC);
    std::vector<Double_t>().swap(etaValuesEpdC);

    nhitsEpdD = 0;
    XnEpdD = 0;
    YnEpdD = 0;
    psiEpdD = D_BAD_VALUE;
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
  //if (event.XnTpc == 0 && event.YnTpc == 0) { event.badEvent = true; }
  //else if (event.XnTpcA == 0 && event.YnTpcA == 0) { event.badEvent = true; }
  if (event.XnTpcB == 0 && event.YnTpcB == 0) { event.badEvent = true; }
  else if (event.XnEpd  == 0 && event.YnEpd  == 0) { event.badEvent = true; }
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
void fillRawSpect(Double_t px, Double_t py, Double_t pz, Double_t mass, TH1D *dndy, TH1D *dndm, TH2D *MvsY)
{
  Double_t y  = rapidity(px, py, pz, mass);
  Double_t mT = transMass(px, py, mass);
  Double_t M  = mT - mass;
  dndy->Fill(y);
  dndm->Fill(M);
  MvsY->Fill(y,M, 1/(TMath::TwoPi() * mT));
};





void FlowAnalyzerABCD(TString inFile, TString jobID)
{
  std::cout << "Initializing..." << std::endl;

  if (gSystem->AccessPathName(inFile)) { std::cout << "Error reading input file!" << std::endl; return;}


  ORDER_N_STR.Form("%d", (Int_t)ORDER_N);
  ORDER_M_STR.Form("%d", (Int_t)ORDER_M);


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


  // INPUT FILE FOR EVENT PLANE RESOLUTION INFORMATION
  /*
  TH1D *h_resolutions;
  Bool_t resolutionsFound = false;
  TString resolutionInputName = "resolutionInfo_INPUT.root";
  TFile *resolutionInputFile;
  if (RUN_ITERATION == 2) 
    { 
      resolutionInputFile = TFile::Open(resolutionInputName, "READ"); 
      if (!resolutionInputFile) { std::cout << "No resolution file was found!" << std::endl; }
      else 
	{ 
	  resolutionsFound = true; 
	  std::cout << "Resolution file found!" << std::endl; 
	}
    }
  */

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
  const char *trackSections[4] = {"Event cuts only", "QA Cuts", "TOF #beta cut", "PID cuts"};  
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

  TH2D *h2_pp_MvsY  = new TH2D("h2_pp_MvsY", "#pi^{+} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}", 32, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_pm_MvsY  = new TH2D("h2_pm_MvsY", "#pi^{-} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}", 32, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_kp_MvsY  = new TH2D("h2_kp_MvsY", "K^{+} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}",   32, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_km_MvsY  = new TH2D("h2_km_MvsY", "K^{-} m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}",   32, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);
  TH2D *h2_pr_MvsY  = new TH2D("h2_pr_MvsY", "Proton m_{T} vs. Rapidity (Weighted);y;m_{T} - m_{0}",  32, -1.6, 0, 60, 0, 3);//50, -2, 0.5, 60, 0, 3);


  TH1D *h_psiTpc_RAW  = new TH1D("h_psiTpc_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", TPC);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcA_RAW = new TH1D("h_psiTpcA_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", TPC A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiTpcB_RAW = new TH1D("h_psiTpcB_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", TPC B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpd_RAW  = new TH1D("h_psiEpd_RAW", "Raw Event Plane Angles (m = "+ORDER_M_STR+", EPD);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdA_RAW = new TH1D("h_psiEpdA_RAW", "Raw Event Plane Angles (n = "+ORDER_M_STR+", EPD A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdB_RAW = new TH1D("h_psiEpdB_RAW", "Raw Event Plane Angles (n = "+ORDER_M_STR+", EPD B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdC_RAW = new TH1D("h_psiEpdC_RAW", "Raw Event Plane Angles (n = "+ORDER_M_STR+", EPD C);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
  TH1D *h_psiEpdD_RAW = new TH1D("h_psiEpdD_RAW", "Raw Event Plane Angles (n = "+ORDER_M_STR+", EPD D);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  TProfile *p_vn_EpdE = new TProfile("p_vn_EpdE", "v_{"+ORDER_N_STR+"} by Centrality (EPD E);Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_EpdF = new TProfile("p_vn_EpdF", "v_{"+ORDER_N_STR+"} by Centrality (EPD F);Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_TpcA = new TProfile("p_vn_TpcA", "v_{"+ORDER_N_STR+"} by Centrality (TPC A);Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_TpcB = new TProfile("p_vn_TpcB", "v_{"+ORDER_N_STR+"} by Centrality (TPC B);Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);

  TProfile *p_vn_pp = new TProfile("p_vn_pp", "#pi^{+} v_{"+ORDER_N_STR+"} by Centrality;Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pm = new TProfile("p_vn_pm", "#pi^{-} v_{"+ORDER_N_STR+"} by Centrality;Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_kp = new TProfile("p_vn_kp", "K^{+} v_{"+ORDER_N_STR+"} by Centrality;Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_km = new TProfile("p_vn_km", "K^{-} v_{"+ORDER_N_STR+"} by Centrality;Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);
  TProfile *p_vn_pr = new TProfile("p_vn_pr", "Proton v_{"+ORDER_N_STR+"} by Centrality;Centrality;<cos("+ORDER_N_STR+"(#phi - #psi_{"+ORDER_M_STR+"}))>", 
				    CENT_BINS, FIRST_CENT, FIRST_CENT+CENT_BINS);


  TH1D *h_psiEpdE_NoAuto = new TH1D("h_psiEpdE_NoAuto", "EP Angles, No Auto-Correlations (m = "+ORDER_M_STR+", EPD E);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

  // Profiles for resolution terms
  TProfile *p_TpcAB = new TProfile("p_TpcAB","TPC A-B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{TPC,B}_{"+ORDER_M_STR+"}))>",
				   CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcAEpdA = new TProfile("p_TpcAEpdA","TPC A EPD A Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{EPD,A}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdB = new TProfile("p_TpcAEpdB","TPC A EPD B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{EPD,B}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdC = new TProfile("p_TpcAEpdC","TPC A EPD C Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{EPD,C}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcAEpdD = new TProfile("p_TpcAEpdD","TPC A EPD D Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,A}_{"+ORDER_M_STR+"}-#psi^{EPD,D}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_TpcBEpdA = new TProfile("p_TpcBEpdA","TPC B EPD A Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,B}_{"+ORDER_M_STR+"}-#psi^{EPD,A}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdB = new TProfile("p_TpcBEpdB","TPC B EPD B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,B}_{"+ORDER_M_STR+"}-#psi^{EPD,B}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdC = new TProfile("p_TpcBEpdC","TPC B EPD C Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,B}_{"+ORDER_M_STR+"}-#psi^{EPD,C}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_TpcBEpdD = new TProfile("p_TpcBEpdD","TPC B EPD D Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{TPC,B}_{"+ORDER_M_STR+"}-#psi^{EPD,D}_{"+ORDER_M_STR+"}))>",
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_EpdAEpdB = new TProfile("p_EpdAEpdB","EPD A EPD B Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{EPD,A}_{"+ORDER_M_STR+"}-#psi^{EPD,B}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_EpdAEpdC = new TProfile("p_EpdAEpdC","EPD A EPD C Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{EPD,A}_{"+ORDER_M_STR+"}-#psi^{EPD,C}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_EpdAEpdD = new TProfile("p_EpdAEpdD","EPD A EPD D Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{EPD,A}_{"+ORDER_M_STR+"}-#psi^{EPD,D}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_EpdBEpdC = new TProfile("p_EpdBEpdC","EPD B EPD C Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{EPD,B}_{"+ORDER_M_STR+"}-#psi^{EPD,C}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile *p_EpdBEpdD = new TProfile("p_EpdBEpdD","EPD B EPD D Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{EPD,B}_{"+ORDER_M_STR+"}-#psi^{EPD,D}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile *p_EpdCEpdD = new TProfile("p_EpdCEpdD","EPD C EPD D Correlations;Centrality;<cos("+ORDER_N_STR+"(#psi^{EPD,C}_{"+ORDER_M_STR+"}-#psi^{EPD,D}_{"+ORDER_M_STR+"}))>", 
				      CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);

  TProfile2D *p2_pp_vs_eta = new TProfile2D("p2_pp_vs_eta","<TnMIP> for Supersectors vs #eta;#eta;Supersector", 400, -6, -2, 12, 0.5, 12.5);
  TH2D *h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  TH2D *h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);

  TH2D *h2_betap  = new TH2D("h2_betap","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_m2_qp   = new TH2D("h2_m2_qp", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  TH2D *h2_m2_vs_qpT  = new TH2D("h2_m2_vs_qpT", "m^{2} vs q*p_{T};q*p_{T} (GeV);m^{2} (GeV^{2})", 300, -3, 3, 300, -0.1, 1.2);
  TH2D *h2_dEdx_vs_qp = new TH2D("h2_dEdx_vs_qp", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);

  TH2D *h2_nSig_vs_qp_pi = new TH2D("h2_nSig_vs_qp_pi", "Pion n#sigma vs q|p|;q|p| (GeV); n#sigma_{#pi}", 500, -5, 5, 400, -8, 8);
  TH2D *h2_nSig_vs_qp_ka = new TH2D("h2_nSig_vs_qp_ka", "Kaon n#sigma vs q|p|;q|p| (GeV); n#sigma_{K}", 500, -5, 5, 400, -8, 8);
  TH2D *h2_nSig_vs_qp_pr = new TH2D("h2_nSig_vs_qp_pr", "Proton n#sigma vs q|p|;q|p| (GeV); n#sigma_{p}", 500, -5, 5, 400, -8, 8);

  TH2D *h2_nSig_vs_qpT_pi = new TH2D("h2_nSig_vs_qpT_pi", "Pion n#sigma vs q*p_{T};q*p_{T} (GeV); n#sigma_{#pi}", 500, -5, 5, 400, -8, 8);
  TH2D *h2_nSig_vs_qpT_ka = new TH2D("h2_nSig_vs_qpT_ka", "Kaon n#sigma vs q*p_{T};q*p_{T} (GeV); n#sigma_{K}", 500, -5, 5, 400, -8, 8);
  TH2D *h2_nSig_vs_qpT_pr = new TH2D("h2_nSig_vs_qpT_pr", "Proton n#sigma vs q*p_{T};q*p_{T} (GeV); n#sigma_{p}", 500, -5, 5, 400, -8, 8);

  TH2D *h2_pi_m2_vs_TPC_nsig = new TH2D("h2_pi_m2_vs_TPC_nsig", "m^{2} vs #pi TPC n#sigma;n#sigma_{#pi};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);
  TH2D *h2_ka_m2_vs_TPC_nsig = new TH2D("h2_ka_m2_vs_TPC_nsig", "m^{2} vs K TPC n#sigma;n#sigma_{K};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);
  TH2D *h2_pr_m2_vs_TPC_nsig = new TH2D("h2_pr_m2_vs_TPC_nsig", "m^{2} vs Proton TPC n#sigma;n#sigma_{pro};m^{2} (GeV^{2})", 500, -5, 5, 500, -0.1, 1.2);

  TH2D *h2_phi_vs_eta_TPC = new TH2D("h2_phi_vs_eta_TPC", "TPC;#eta;#phi", 300, -2.2, 0.2, 300, -4, 4);
  TH2D *h2_phi_vs_eta_EPD = new TH2D("h2_phi_vs_eta_EPD", "EPD;#eta;#phi", 300, -6, -2.5, 300, -4, 4);

  TH2D *h2_y_vs_eta = new TH2D("h2_y_vs_eta", "TPC All Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);

  TH2D *h2_y_vs_eta_pp = new TH2D("h2_y_vs_eta_pp", "TPC #pi^{+} y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  /*
  TH2D *h2_y_vs_eta_pt0p5to1_pp = new TH2D("h2_y_vs_eta_pt0p5to1_pp", "#pi^{+} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pp = new TH2D("h2_y_vs_eta_pt1to1p5_pp", "#pi^{+} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pp = new TH2D("h2_y_vs_eta_pt1p5to2_pp", "#pi^{+} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  */
  TH2D *h2_y_vs_eta_pm = new TH2D("h2_y_vs_eta_pm", "TPC #pi^{-} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  /*
  TH2D *h2_y_vs_eta_pt0p5to1_pm = new TH2D("h2_y_vs_eta_pt0p5to1_pm", "#pi^{-} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pm = new TH2D("h2_y_vs_eta_pt1to1p5_pm", "#pi^{-} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pm = new TH2D("h2_y_vs_eta_pt1p5to2_pm", "#pi^{-} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  */
  TH2D *h2_y_vs_eta_kp = new TH2D("h2_y_vs_eta_kp", "TPC K^{+} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  /*
  TH2D *h2_y_vs_eta_pt0p5to1_kp = new TH2D("h2_y_vs_eta_pt0p5to1_kp", "K^{+} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_kp = new TH2D("h2_y_vs_eta_pt1to1p5_kp", "K^{+} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_kp = new TH2D("h2_y_vs_eta_pt1p5to2_kp", "K^{+} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  */
  TH2D *h2_y_vs_eta_km = new TH2D("h2_y_vs_eta_km", "TPC K^{-} Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  /*
  TH2D *h2_y_vs_eta_pt0p5to1_km = new TH2D("h2_y_vs_eta_pt0p5to1_km", "K^{-} y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_km = new TH2D("h2_y_vs_eta_pt1to1p5_km", "K^{-} y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_km = new TH2D("h2_y_vs_eta_pt1p5to2_km", "K^{-} y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  */
  TH2D *h2_y_vs_eta_pr = new TH2D("h2_y_vs_eta_pr", "TPC Proton Charged y vs #eta;#eta;y", 40, -2, 0, 40, -2, 0);
  /*
  TH2D *h2_y_vs_eta_pt0p5to1_pr = new TH2D("h2_y_vs_eta_pt0p5to1_pr", "Proton y vs #eta (0.5 < p_{T} < 1 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1to1p5_pr = new TH2D("h2_y_vs_eta_pt1to1p5_pr", "Proton y vs #eta (1 < p_{T} < 1.5 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  TH2D *h2_y_vs_eta_pt1p5to2_pr = new TH2D("h2_y_vs_eta_pt1p5to2_pr", "Proton y vs #eta (1.5 < p_{T} < 2 GeV);#eta;y", 40, -2, 0, 40, -2, 0);
  */

  TH2D *h2_pT_vs_yCM_pp = new TH2D("h2_pT_vs_yCM_pp", "#pi^{+};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 2.5);
  TH2D *h2_pT_vs_yCM_pm = new TH2D("h2_pT_vs_yCM_pm", "#pi^{-};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 2.5);
  TH2D *h2_pT_vs_yCM_kp = new TH2D("h2_pT_vs_yCM_kp", "K^{+};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 2.5);
  TH2D *h2_pT_vs_yCM_km = new TH2D("h2_pT_vs_yCM_km", "K^{-};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 2.5);
  TH2D *h2_pT_vs_yCM_pr = new TH2D("h2_pT_vs_yCM_pr", "Proton;y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 3.0);


  

  /*
  TH2D *h2_phiSearchTpc = new TH2D("h2_phiSearchTpc", "Azimuthal Distribution by Centrality;#phi;Centrality (%)", 
				   200, -4, 4, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TH2D *h2_phiSearchEpd = new TH2D("h2_phiSearchEpd", "Azimuthal Distribution by Centrality;#phi;Centrality (%)", 
				   200, -4, 4, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  */
  // Here the name refers to the eta region that will be displayed/searched using the event plane angle from the opposite region

  TProfile2D *h2_v2ScanTpc = new TProfile2D("h2_v2ScanTpc", "<cos("+ORDER_N_STR+"(#phi^{TPC} - #psi^{EPD}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
					      12, -2, 0, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2ScanEpd = new TProfile2D("h2_v2ScanEpd", "<cos("+ORDER_N_STR+"(#phi^{EPD} - #psi^{TPC}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
					      12, -6.0, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2ScanEpdTpcA = new TProfile2D("h2_v2ScanEpdTpcA", "<cos("+ORDER_N_STR+"(#phi^{EPD} - #psi^{TPC,A}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
						  12, -6.0, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  TProfile2D *h2_v2ScanEpdTpcB = new TProfile2D("h2_v2ScanEpdTpcB", "<cos("+ORDER_N_STR+"(#phi^{EPD} - #psi^{TPC,B}_{"+ORDER_M_STR+"}))>;#eta;Centrality (%)", 
						  12, -6.0, -2.3, CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  h2_v2ScanTpc->SetStats(0);
  h2_v2ScanEpd->SetStats(0);
  h2_v2ScanEpdTpcA->SetStats(0);
  h2_v2ScanEpdTpcB->SetStats(0);

  // The indices here are equivalent to the corresponding centrality ID
  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};

  // The indices here are opposite to the corresponding centrality ID (array is backward)
  //const char *centralityBins[16] = {"0-5", "5-10", "10-15" "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80"};


  TH2D *h2_psiEpdATpcA = new TH2D("h2_psiEpdATpcA", "#psi^{EPD,A} vs #psi^{TPC,A} (Order "+ORDER_M_STR+");#psi^{TPC}_{A};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcA = new TH2D("h2_psiEpdBTpcA", "#psi^{EPD,B} vs #psi^{TPC,A} (Order "+ORDER_M_STR+");#psi^{TPC}_{A};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdCTpcA = new TH2D("h2_psiEpdCTpcA", "#psi^{EPD,C} vs #psi^{TPC,A} (Order "+ORDER_M_STR+");#psi^{TPC}_{A};#psi^{EPD}_{C}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdDTpcA = new TH2D("h2_psiEpdDTpcA", "#psi^{EPD,D} vs #psi^{TPC,A} (Order "+ORDER_M_STR+");#psi^{TPC}_{A};#psi^{EPD}_{D}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);

  TH2D *h2_psiEpdATpcB = new TH2D("h2_psiEpdATpcB", "#psi^{EPD,A} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{EPD}_{A}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdBTpcB = new TH2D("h2_psiEpdBTpcB", "#psi^{EPD,B} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{EPD}_{B}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdCTpcB = new TH2D("h2_psiEpdCTpcB", "#psi^{EPD,C} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{EPD}_{C}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  TH2D *h2_psiEpdDTpcB = new TH2D("h2_psiEpdDTpcB", "#psi^{EPD,D} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{EPD}_{D}", 
				  200, -PSI_BOUNDS, PSI_BOUNDS, 200, -PSI_BOUNDS, PSI_BOUNDS);
  
  TH2D *h2_psiTpcATpcB = new TH2D("h2_psiTpcATpcB", "#psi^{TPC,A} vs #psi^{TPC,B} (Order "+ORDER_M_STR+");#psi^{TPC}_{B};#psi^{TPC}_{A}", 
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

	  h_track_check->Fill(trackSections[1], 1);

	  TVector3 mom_vec  = picoTrack->pMom();
	  Double_t d_dEdx   = picoTrack->dEdx();
	  Double_t d_charge = picoTrack->charge();
	  Double_t d_mom    = picoTrack->pPtot();
	  Double_t d_pT     = picoTrack->pPt();
	  Double_t d_px     = picoTrack->pMom().x();
	  Double_t d_py     = picoTrack->pMom().y();
	  Double_t d_pz     = picoTrack->pMom().z();
	  Double_t d_phi    = mom_vec.Phi();
	  Double_t d_eta    = mom_vec.Eta();
	  Double_t d_TPCnSigmaPion   = picoTrack->nSigmaPion();
	  Double_t d_TPCnSigmaKaon   = picoTrack->nSigmaKaon();
	  Double_t d_TPCnSigmaProton = picoTrack->nSigmaProton();

	  // Get event planes from the TPC here before the TOF cut
	  if (d_charge != 0)
	    {
	      //h_eta_TPC_s->Fill(d_eta - Y_MID);
	      h2_phi_vs_eta_TPC->Fill(d_eta, d_phi);

	      eventInfo.nTracksTpc++;
	      if (d_eta > Y_MID)        // Account for Q vector sign change past mid-rapidity.
		{
		  eventInfo.XnTpc += d_pT * TMath::Cos(ORDER_M * d_phi);
		  eventInfo.YnTpc += d_pT * TMath::Sin(ORDER_M * d_phi);
		}
	      else if (d_eta < Y_MID)
		{
		  eventInfo.XnTpc -= d_pT * TMath::Cos(ORDER_M * d_phi);
		  eventInfo.YnTpc -= d_pT * TMath::Sin(ORDER_M * d_phi);
		}


	      if (d_eta > MIN_TPC_ETA_CUT && d_eta < AGAP_TPC_ETA_CUT)          // TPC A  (sign change happens later)
		{
		  eventInfo.nTracksTpcA++;
		  eventInfo.XnTpcA += d_pT * TMath::Cos(ORDER_M * d_phi);
		  eventInfo.YnTpcA += d_pT * TMath::Sin(ORDER_M * d_phi);
		  eventInfo.phiValuesTpcA.push_back(d_phi);
		  eventInfo.etaValuesTpcA.push_back(d_eta);
		}
	      else if (d_eta > /*GAPB_TPC_ETA_CUT*/-0.5 && d_eta < MAX_TPC_ETA_CUT)     // TPC B
		{
		  eventInfo.nTracksTpcB++;
		  eventInfo.XnTpcB += d_pT * TMath::Cos(ORDER_M * d_phi);
		  eventInfo.YnTpcB += d_pT * TMath::Sin(ORDER_M * d_phi);
		  eventInfo.phiValuesTpcB.push_back(d_phi);
		  eventInfo.etaValuesTpcB.push_back(d_eta);
		}
	    }// End TPC event planes


	  //=========================================================
	  //          TOF Beta Cuts
	  //=========================================================
	  StPicoBTofPidTraits *trait;
	  Double_t d_tofBeta;
	  Double_t d_m2;
	  Bool_t tofTrack = picoTrack->isTofTrack();

	  if (tofTrack)
	    {
	      trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());
	      d_tofBeta = trait->btofBeta();
	      d_m2 = d_mom*d_mom*( (1 / (d_tofBeta*d_tofBeta)) - 1);
	      
	      if (d_tofBeta < 0.01) continue;
	    }
	  //=========================================================
	  //          End TOF Beta Cuts
	  //=========================================================

	  h_track_check->Fill(trackSections[2], 1);


	  // Fill histos and save important event info in the custom struct type
	  if (d_charge != 0)
	    {
	      if (tofTrack)
		{
		  h2_betap->Fill(d_charge * d_mom, 1/d_tofBeta);
		  h2_m2_qp->Fill(d_charge * d_mom, d_mom * d_mom * (1/(d_tofBeta*d_tofBeta) - 1));
		  h2_m2_vs_qpT->Fill(d_charge * d_pT, d_m2);

		  h2_pi_m2_vs_TPC_nsig->Fill(d_TPCnSigmaPion, d_m2);
		  h2_ka_m2_vs_TPC_nsig->Fill(d_TPCnSigmaKaon, d_m2);
		  h2_pr_m2_vs_TPC_nsig->Fill(d_TPCnSigmaProton, d_m2);
		}

	      h2_dEdx_vs_qp->Fill(d_charge * d_mom, d_dEdx);

	      h2_nSig_vs_qp_pi->Fill(d_charge * d_mom, d_TPCnSigmaPion);
	      h2_nSig_vs_qp_ka->Fill(d_charge * d_mom, d_TPCnSigmaKaon);
	      h2_nSig_vs_qp_pr->Fill(d_charge * d_mom, d_TPCnSigmaProton);

	      h2_nSig_vs_qpT_pi->Fill(d_charge * d_pT, d_TPCnSigmaPion);
	      h2_nSig_vs_qpT_ka->Fill(d_charge * d_pT, d_TPCnSigmaKaon);
	      h2_nSig_vs_qpT_pr->Fill(d_charge * d_pT, d_TPCnSigmaProton);

	      //=========================================================
	      //          PID Cuts
	      //=========================================================
	      Bool_t pion = false;
	      Bool_t kaon = false;
	      Bool_t proton = (d_TPCnSigmaProton > -2) && (d_TPCnSigmaProton < 2);

	      if (tofTrack)
		{
		  pion = (d_TPCnSigmaPion > -3) && (d_TPCnSigmaPion < 3) && (d_m2 > -0.1) && (d_m2 < 0.1);
		  kaon = (d_TPCnSigmaKaon > -3) && (d_TPCnSigmaKaon < 3) && (d_m2 > 0.15) && (d_m2 < 0.34);
		}
	    
	      if (!pion && !kaon && !proton) continue;

	      if (pion && proton) { proton = false; }
	      if (kaon && proton) { proton = false; }
	      if (pion && kaon) continue;
	      //=========================================================
	      //          END PID Cuts
	      //=========================================================

	      h_track_check->Fill(trackSections[3], 1);

	      Double_t d_m0_pi = 0.1396;   //Rest masses
	      Double_t d_m0_ka = 0.4937;
	      Double_t d_m0_pr = 0.9383;
	      Double_t d_rapidity;

	      if (pion)
		{
		  if (d_charge > 0) 
		    {
		      d_rapidity = rapidity(d_px, d_py, d_pz, d_m0_pi);

		      h2_pT_vs_yCM_pp->Fill(d_rapidity - Y_MID, d_pT);
			  
		      if (d_rapidity - Y_MID > 0.0 && d_rapidity - Y_MID < 0.545 && d_pT >= 0.18 && d_pT <= 1.6)
			{
			  fillRawSpect(d_px, d_py, d_pz, d_m0_pi, h_pp_dndy, h_pp_dndm, h2_pp_MvsY);
			  h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  h2_y_vs_eta_pp->Fill(d_eta, d_rapidity);
			  h_pp_pT->Fill(d_pT);
			  h_pp_mom->Fill(d_mom);

			  eventInfo.phiValuesPionP.push_back(d_phi);
			}

		      /*
		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pp->Fill(d_eta, d_rapidity); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pp->Fill(d_eta, d_rapidity); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pp->Fill(d_eta, d_rapidity); }
		      */
		    }
		  else if (d_charge < 0) 
		    {
		      d_rapidity = rapidity(d_px, d_py, d_pz, d_m0_pi);

		      h2_pT_vs_yCM_pm->Fill(d_rapidity - Y_MID, d_pT);

		      if (d_rapidity - Y_MID > 0.0 && d_rapidity - Y_MID < 0.545 && d_pT >= 0.18 && d_pT <= 1.6)
			{
			  fillRawSpect(d_px, d_py, d_pz, d_m0_pi, h_pm_dndy, h_pm_dndm, h2_pm_MvsY);
			  h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  h2_y_vs_eta_pm->Fill(d_eta, d_rapidity);
			  h_pm_pT->Fill(d_pT);
			  h_pm_mom->Fill(d_mom);

			  eventInfo.phiValuesPionM.push_back(d_phi);
			}

		      /*
		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pm->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pi)); }
		      */
		    }
		}
	      else if (kaon)
		{
		  if (d_charge > 0) 
		    {
		      d_rapidity = rapidity(d_px, d_py, d_pz, d_m0_ka);

		      h2_pT_vs_yCM_kp->Fill(d_rapidity - Y_MID, d_pT);

		      if (d_rapidity - Y_MID > 0.0 && d_rapidity - Y_MID < 0.545 && d_pT >= 0.4 && d_pT <= 1.6)
			{
			  fillRawSpect(d_px, d_py, d_pz, d_m0_ka, h_kp_dndy, h_kp_dndm, h2_kp_MvsY);
			  h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  h2_y_vs_eta_kp->Fill(d_eta, d_rapidity);
			  h_kp_pT->Fill(d_pT);
			  h_kp_mom->Fill(d_mom);

			  eventInfo.phiValuesKaonP.push_back(d_phi);
			}

		      /*
		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_kp->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      */
		    }
		  else if (d_charge < 0)		 
		    {
		      d_rapidity = rapidity(d_px, d_py, d_pz, d_m0_ka);

		      h2_pT_vs_yCM_km->Fill(d_rapidity - Y_MID, d_pT);

		      if (d_rapidity - Y_MID > 0.0 && d_rapidity - Y_MID < 0.545 && d_pT >= 0.4 && d_pT <= 1.6)
			{
			  fillRawSpect(d_px, d_py, d_pz, d_m0_ka, h_km_dndy, h_km_dndm, h2_km_MvsY);
			  h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  h2_y_vs_eta_km->Fill(d_eta, d_rapidity);
			  h_km_pT->Fill(d_pT);
			  h_km_mom->Fill(d_mom);

			  eventInfo.phiValuesKaonM.push_back(d_phi);
			}

		      /*
		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_km->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_ka)); }
		      */
		    }
		}
	      else if (proton)
		{
		  if (d_charge > 0) 
		    {
		      d_rapidity = rapidity(d_px, d_py, d_pz, d_m0_pr);

		      h2_pT_vs_yCM_pr->Fill(d_rapidity - Y_MID, d_pT);

		      if (d_rapidity - Y_MID > 0.0 && d_rapidity - Y_MID < 0.545 && d_pT >= 0.4 && d_pT <= 2.0)
			{
			  fillRawSpect(d_px, d_py, d_pz, d_m0_pr, h_pr_dndy, h_pr_dndm, h2_pr_MvsY);
			  h2_y_vs_eta->Fill(d_eta, d_rapidity);
			  h2_y_vs_eta_pr->Fill(d_eta, d_rapidity);
			  h_pr_pT->Fill(d_pT);
			  h_pr_mom->Fill(d_mom);

			  eventInfo.phiValuesProton.push_back(d_phi);
			}

		      /*
		      if (d_pT > 0.5 && d_pT < 1.0) { h2_y_vs_eta_pt0p5to1_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr)); }
		      if (d_pT > 1.0 && d_pT < 1.5) { h2_y_vs_eta_pt1to1p5_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr)); }
		      if (d_pT > 1.5 && d_pT < 2.0) { h2_y_vs_eta_pt1p5to2_pr->Fill(d_eta, rapidity(d_px, d_py, d_pz, d_m0_pr)); }
		      */
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
      int tileRow;
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
	  tileRow = epdHit->row();
	  tileEta = tileVector.Eta();
	  tilePhi = tileVector.Phi();
	  //tileWeight = (epdHit->nMIP() > EPD_THRESHOLD) ? ( (epdHit->nMIP() > EPD_MAX_WEIGHT)?EPD_MAX_WEIGHT:epdHit->nMIP() ) : 0;
	  tileWeight = (epdHit->nMIP() > EPD_THRESHOLD) ? 1 : 0;

	  p2_pp_vs_eta->Fill(tileEta, tileSector, tileWeight);
	  h_tileWeights->Fill(tileWeight);
	  h2_phi_vs_eta_EPD->Fill(tileEta, tilePhi);

	  eventInfo.nHitsEpd++;
	  eventInfo.XnEpd += tileWeight * TMath::Cos(ORDER_M * tilePhi);
	  eventInfo.YnEpd += tileWeight * TMath::Sin(ORDER_M * tilePhi);


	  //if (tileRow <= 4)  // Sub A
	  if (tileRow == 2 || tileRow == 3)  // Sub A
	    {
	      eventInfo.nhitsEpdA++;
	      eventInfo.XnEpdA += tileWeight * TMath::Cos(ORDER_M * tilePhi);
	      eventInfo.YnEpdA += tileWeight * TMath::Sin(ORDER_M * tilePhi);
	      eventInfo.phiValuesEpdA.push_back(tilePhi);
	      eventInfo.etaValuesEpdA.push_back(tileEta);
	    }
	  //else if (tileRow >= 5 && tileRow <= 8)  // Sub B
	  else if (tileRow == 6 || tileRow == 7)  // Sub B
	    {
	      eventInfo.nhitsEpdB++;
	      eventInfo.XnEpdB += tileWeight * TMath::Cos(ORDER_M * tilePhi);
	      eventInfo.YnEpdB += tileWeight * TMath::Sin(ORDER_M * tilePhi);
	      eventInfo.phiValuesEpdB.push_back(tilePhi);
	      eventInfo.etaValuesEpdB.push_back(tileEta);
	    }
	  //else if (tileRow >= 9 && tileRow <= 12)  // Sub C
	  else if (tileRow == 10 || tileRow == 11)  // Sub C
	    {
	      eventInfo.nhitsEpdC++;
	      eventInfo.XnEpdC += tileWeight * TMath::Cos(ORDER_M * tilePhi);
	      eventInfo.YnEpdC += tileWeight * TMath::Sin(ORDER_M * tilePhi);
	      eventInfo.phiValuesEpdC.push_back(tilePhi);
	      eventInfo.etaValuesEpdC.push_back(tileEta);
	    }
	  //else if (tileRow >= 13 && tileRow <= 16)  // Sub D
	  else if (tileRow == 14 || tileRow == 15)  // Sub D
	    {
	      eventInfo.nhitsEpdD++;
	      eventInfo.XnEpdD += tileWeight * TMath::Cos(ORDER_M * tilePhi);
	      eventInfo.YnEpdD += tileWeight * TMath::Sin(ORDER_M * tilePhi);
	      eventInfo.phiValuesEpdD.push_back(tilePhi);
	      eventInfo.etaValuesEpdD.push_back(tileEta);
	    }

	}
      delete epdGeom;
      //=========================================================
      //            END EPD STUFF
      //=========================================================


      //if (eventInfo.nTracksTpc  < MIN_TRACKS) continue;
      //if (eventInfo.nTracksTpcA < MIN_TRACKS) continue;
      if (eventInfo.nTracksTpcB < MIN_TRACKS) continue;
      if (eventInfo.nHitsEpd    < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdA   < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdB   < MIN_TRACKS) continue;
      if (eventInfo.nhitsEpdC   < MIN_TRACKS+2) continue;
      if (eventInfo.nhitsEpdD   < MIN_TRACKS+2) continue;
      

      checkZeroQ(eventInfo);
      if (eventInfo.badEvent) continue;


      // RAW SUB-EVENT PLANE ANGLES //
      if (ORDER_M == 1)           // Q vectors must change sign past mid-rapidity; I think this is for 1st order event-planes only. Full TPC already takes this into account.
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
	  eventInfo.XnEpdD *= -1.0;
	  eventInfo.YnEpdD *= -1.0;
	}
      eventInfo.psiTpc  = TMath::ATan2(eventInfo.YnTpc,  eventInfo.XnTpc)  / ORDER_M;
      eventInfo.psiTpcA = TMath::ATan2(eventInfo.YnTpcA, eventInfo.XnTpcA) / ORDER_M;
      eventInfo.psiTpcB = TMath::ATan2(eventInfo.YnTpcB, eventInfo.XnTpcB) / ORDER_M;
      eventInfo.psiEpd  = TMath::ATan2(eventInfo.YnEpd,  eventInfo.XnEpd)  / ORDER_M;
      eventInfo.psiEpdA = TMath::ATan2(eventInfo.YnEpdA, eventInfo.XnEpdA) / ORDER_M;
      eventInfo.psiEpdB = TMath::ATan2(eventInfo.YnEpdB, eventInfo.XnEpdB) / ORDER_M;
      eventInfo.psiEpdC = TMath::ATan2(eventInfo.YnEpdC, eventInfo.XnEpdC) / ORDER_M;
      eventInfo.psiEpdD = TMath::ATan2(eventInfo.YnEpdD, eventInfo.XnEpdD) / ORDER_M;


      // Fill eta/phi distributions here since this is past all possible cuts.
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcA.size(); i++) 
	{ 
	  h_eta_s->Fill(eventInfo.etaValuesTpcA.at(i) - Y_MID); 
	  h_eta_TPC_s->Fill(eventInfo.etaValuesTpcA.at(i) - Y_MID); 
	} 
      for (unsigned int i = 0; i < eventInfo.etaValuesTpcB.size(); i++) 
	{ 
	  h_eta_s->Fill(eventInfo.etaValuesTpcB.at(i) - Y_MID); 
	  h_eta_TPC_s->Fill(eventInfo.etaValuesTpcB.at(i) - Y_MID); 
	} 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdA.size(); i++) { h_eta_s->Fill(eventInfo.etaValuesEpdA.at(i) - Y_MID); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdB.size(); i++) { h_eta_s->Fill(eventInfo.etaValuesEpdB.at(i) - Y_MID); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdC.size(); i++) { h_eta_s->Fill(eventInfo.etaValuesEpdC.at(i) - Y_MID); } 
      for (unsigned int i = 0; i < eventInfo.etaValuesEpdD.size(); i++) { h_eta_s->Fill(eventInfo.etaValuesEpdD.at(i) - Y_MID); } 


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
      h_psiTpcA_RAW->Fill(eventInfo.psiTpcA);
      h_psiTpcB_RAW->Fill(eventInfo.psiTpcB);
      h_psiEpd_RAW->Fill(eventInfo.psiEpd);
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

      h_psiTpc_RC  = new TH1D("h_psiTpc_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", TPC);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiTpcA_RC = new TH1D("h_psiTpcA_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", TPC A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiTpcB_RC = new TH1D("h_psiTpcB_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", TPC B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpd_RC  = new TH1D("h_psiEpd_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdA_RC = new TH1D("h_psiEpdA_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD A);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdB_RC = new TH1D("h_psiEpdB_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD B);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC_RC = new TH1D("h_psiEpdC_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD C);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdD_RC = new TH1D("h_psiEpdD_RC", "Re-centered Event Plane Angles (m = "+ORDER_M_STR+", EPD D);#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

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
	  v_events.at(i).psiTpc  = TMath::ATan2(v_events.at(i).YnTpc,  v_events.at(i).XnTpc)  / ORDER_M; 
	  v_events.at(i).psiTpcA = TMath::ATan2(v_events.at(i).YnTpcA, v_events.at(i).XnTpcA) / ORDER_M; 
	  v_events.at(i).psiTpcB = TMath::ATan2(v_events.at(i).YnTpcB, v_events.at(i).XnTpcB) / ORDER_M; 
	  v_events.at(i).psiEpd  = TMath::ATan2(v_events.at(i).YnEpd,  v_events.at(i).XnEpd)  / ORDER_M;
	  v_events.at(i).psiEpdA = TMath::ATan2(v_events.at(i).YnEpdA, v_events.at(i).XnEpdA) / ORDER_M; 
	  v_events.at(i).psiEpdB = TMath::ATan2(v_events.at(i).YnEpdB, v_events.at(i).XnEpdB) / ORDER_M; 
	  v_events.at(i).psiEpdC = TMath::ATan2(v_events.at(i).YnEpdC, v_events.at(i).XnEpdC) / ORDER_M; 
	  v_events.at(i).psiEpdD = TMath::ATan2(v_events.at(i).YnEpdD, v_events.at(i).XnEpdD) / ORDER_M; 

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_M);
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_M);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_M);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_M);
	  v_events.at(i).psiEpdA = angleShift(v_events.at(i).psiEpdA, ORDER_M);
	  v_events.at(i).psiEpdB = angleShift(v_events.at(i).psiEpdB, ORDER_M);
	  v_events.at(i).psiEpdC = angleShift(v_events.at(i).psiEpdC, ORDER_M);
	  v_events.at(i).psiEpdD = angleShift(v_events.at(i).psiEpdD, ORDER_M);

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
	      p_sinAvgsTpc->Fill(j,  TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiTpc));
	      p_cosAvgsTpc->Fill(j,  TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiTpc));
	      p_sinAvgsTpcA->Fill(j, TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiTpcA));
	      p_cosAvgsTpcA->Fill(j, TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiTpcA));
	      p_sinAvgsTpcB->Fill(j, TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiTpcB));
	      p_cosAvgsTpcB->Fill(j, TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiTpcB));
	      p_sinAvgsEpd->Fill(j,  TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpd));
	      p_cosAvgsEpd->Fill(j,  TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpd));
	      p_sinAvgsEpdA->Fill(j, TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdA));
	      p_cosAvgsEpdA->Fill(j, TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdA));
	      p_sinAvgsEpdB->Fill(j, TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdB));
	      p_cosAvgsEpdB->Fill(j, TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdB));
	      p_sinAvgsEpdC->Fill(j, TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdC));
	      p_cosAvgsEpdC->Fill(j, TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdC));
	      p_sinAvgsEpdD->Fill(j, TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdD));
	      p_cosAvgsEpdD->Fill(j, TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdD));
	    }

	}// End loop over v_events
      std::cout << "Bad Events after re-centering: " << badEvents << std::endl;
    }
  //=========================================================
  //          End Re-centering
  //=========================================================



  if (RUN_ITERATION == 2)
    {
      //=========================================================
      //          Event Plane Angle Shifting
      //=========================================================
      std::cout << "Performing event plane angle shifting..." << std::endl;

      Int_t numOfEvents = v_events.size();

      h_psiTpc_FLAT  = new TH1D("h_psiTpc_FLAT", "Flattened Event Plane Angle (TPC, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiTpcA_FLAT = new TH1D("h_psiTpcA_FLAT", "Flattened Event Plane Angle (TPC A, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiTpcB_FLAT = new TH1D("h_psiTpcB_FLAT", "Flattened Event Plane Angle (TPC B, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);      
      h_psiEpd_FLAT  = new TH1D("h_psiEpd_FLAT", "Flattened Event Plane Angle (EPD, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdA_FLAT = new TH1D("h_psiEpdA_FLAT", "Flattened Event Plane Angle (EPD A, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdB_FLAT = new TH1D("h_psiEpdB_FLAT", "Flattened Event Plane Angle (EPD B, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdC_FLAT = new TH1D("h_psiEpdC_FLAT", "Flattened Event Plane Angle (EPD C, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);
      h_psiEpdD_FLAT = new TH1D("h_psiEpdD_FLAT", "Flattened Event Plane Angle (EPD D, order "+ORDER_M_STR+");#psi_{"+ORDER_M_STR+"};Events", 400, -PSI_BOUNDS, PSI_BOUNDS);

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
	      jthSinAvg_Tpc  = p_sinAvgsTpc_INPUT->GetBinContent(j);
	      jthCosAvg_Tpc  = p_cosAvgsTpc_INPUT->GetBinContent(j);
	      jthSinAvg_TpcA = p_sinAvgsTpcA_INPUT->GetBinContent(j);
	      jthCosAvg_TpcA = p_cosAvgsTpcA_INPUT->GetBinContent(j);
	      jthSinAvg_TpcB = p_sinAvgsTpcB_INPUT->GetBinContent(j);
	      jthCosAvg_TpcB = p_cosAvgsTpcB_INPUT->GetBinContent(j);
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

	      psiTpc_delta  += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_Tpc * TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiTpc) 
							      +jthCosAvg_Tpc * TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiTpc));
	      psiTpcA_delta += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_TpcA * TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiTpcA) 
							      +jthCosAvg_TpcA * TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiTpcA));
	      psiTpcB_delta += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_TpcB * TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiTpcB) 
							      +jthCosAvg_TpcB * TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiTpcB));
	      psiEpd_delta  += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_Epd * TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpd)
							      +jthCosAvg_Epd * TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpd));
	      psiEpdA_delta += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_EpdA*TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdA)
							      +jthCosAvg_EpdA*TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdA));
	      psiEpdB_delta += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_EpdB*TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdB)
							      +jthCosAvg_EpdB*TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdB));
	      psiEpdC_delta += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_EpdC*TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdC)
							      +jthCosAvg_EpdC*TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdC));
	      psiEpdD_delta += (2.0/((Double_t)j*ORDER_M)) * (-jthSinAvg_EpdD*TMath::Cos((Double_t)j * ORDER_M * v_events.at(i).psiEpdD)
							      +jthCosAvg_EpdD*TMath::Sin((Double_t)j * ORDER_M * v_events.at(i).psiEpdD));
	    }


	  v_events.at(i).psiTpc  += psiTpc_delta;
	  v_events.at(i).psiTpcA += psiTpcA_delta;
	  v_events.at(i).psiTpcB += psiTpcB_delta;
	  v_events.at(i).psiEpd  += psiEpd_delta;
	  v_events.at(i).psiEpdA += psiEpdA_delta;
	  v_events.at(i).psiEpdB += psiEpdB_delta;
	  v_events.at(i).psiEpdC += psiEpdC_delta;
	  v_events.at(i).psiEpdD += psiEpdD_delta;

	  v_events.at(i).psiTpc  = angleShift(v_events.at(i).psiTpc,  ORDER_M);
	  v_events.at(i).psiTpcA = angleShift(v_events.at(i).psiTpcA, ORDER_M);
	  v_events.at(i).psiTpcB = angleShift(v_events.at(i).psiTpcB, ORDER_M);
	  v_events.at(i).psiEpd  = angleShift(v_events.at(i).psiEpd,  ORDER_M);
	  v_events.at(i).psiEpdA = angleShift(v_events.at(i).psiEpdA, ORDER_M);
	  v_events.at(i).psiEpdB = angleShift(v_events.at(i).psiEpdB, ORDER_M);
	  v_events.at(i).psiEpdC = angleShift(v_events.at(i).psiEpdC, ORDER_M);
	  v_events.at(i).psiEpdD = angleShift(v_events.at(i).psiEpdD, ORDER_M);

	  h_psiTpc_FLAT->Fill(v_events.at(i).psiTpc);
	  h_psiTpcA_FLAT->Fill(v_events.at(i).psiTpcA);
	  h_psiTpcB_FLAT->Fill(v_events.at(i).psiTpcB);
	  h_psiEpd_FLAT->Fill(v_events.at(i).psiEpd);
	  h_psiEpdA_FLAT->Fill(v_events.at(i).psiEpdA);
	  h_psiEpdB_FLAT->Fill(v_events.at(i).psiEpdB);
	  h_psiEpdC_FLAT->Fill(v_events.at(i).psiEpdC);
	  h_psiEpdD_FLAT->Fill(v_events.at(i).psiEpdD);
	  //=========================================================
	  //          End Event Plane Angle Shifting
	  //=========================================================


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
	  p_TpcAB->Fill(v_events.at(i).centID,    TMath::Cos(ORDER_N * (v_events.at(i).psiTpcA - v_events.at(i).psiTpcB)));

	  p_TpcAEpdA->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdA)));
	  p_TpcAEpdB->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdB)));
	  p_TpcAEpdC->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdC)));
	  p_TpcAEpdD->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcA - v_events.at(i).psiEpdD)));
	  p_TpcBEpdA->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdA)));
	  p_TpcBEpdB->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdB)));
	  p_TpcBEpdC->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdC)));
	  p_TpcBEpdD->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiTpcB - v_events.at(i).psiEpdD)));

	  p_EpdAEpdB->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiEpdA - v_events.at(i).psiEpdB)));
	  p_EpdAEpdC->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiEpdA - v_events.at(i).psiEpdC)));
	  p_EpdAEpdD->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiEpdA - v_events.at(i).psiEpdD)));
	  p_EpdBEpdC->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiEpdB - v_events.at(i).psiEpdC)));
	  p_EpdBEpdD->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiEpdB - v_events.at(i).psiEpdD)));
	  p_EpdCEpdD->Fill(v_events.at(i).centID, TMath::Cos(ORDER_N * (v_events.at(i).psiEpdC - v_events.at(i).psiEpdD)));
	  //


	  //=========================================================
	  //          v_n Scan Plots
	  //=========================================================

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

	      h2_v2ScanEpd->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpc)));
	      h2_v2ScanEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcA)));
	      h2_v2ScanEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsB; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdB.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdB.at(k);

	      h2_v2ScanEpd->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpc)));
	      h2_v2ScanEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcA)));
	      h2_v2ScanEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsC; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdC.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdC.at(k);

	      h2_v2ScanEpd->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpc)));
	      h2_v2ScanEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcA)));
	      h2_v2ScanEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }
	  for (int k = 0; k < epdHitsD; k++)
	    {
	      phiEpd = v_events.at(i).phiValuesEpdD.at(k);
	      etaEpd = v_events.at(i).etaValuesEpdD.at(k);

	      h2_v2ScanEpd->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpc)));
	      h2_v2ScanEpdTpcA->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcA)));
	      h2_v2ScanEpdTpcB->Fill(etaEpd, centralityID, TMath::Cos(ORDER_M * (phiEpd - psiTpcB)));
	      //h2_phiSearchEpd->Fill(phiEpd, centralityID);
	    }

	  for (int k = 0; k < tpcTracksA; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcA.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcA.at(k);

	      h2_v2ScanTpc->Fill(etaTpc, centralityID, TMath::Cos(ORDER_M * (phiTpc - psiEpd)));
	      //h2_phiSearchTpc->Fill(phiTpc, centralityID);
	    }
	  for (int k = 0; k < tpcTracksB; k++)
	    {
	      phiTpc = v_events.at(i).phiValuesTpcB.at(k);
	      etaTpc = v_events.at(i).etaValuesTpcB.at(k);

	      h2_v2ScanTpc->Fill(etaTpc, centralityID, TMath::Cos(ORDER_M * (phiTpc - psiEpd)));
	      //h2_phiSearchTpc->Fill(phiTpc, centralityID);
	    }

	  //=========================================================
	  //          End v_n Scan Plots
	  //=========================================================

	  /*
	  //=========================================================
	  //        RAW Flow Calculations
	  //=========================================================
	  Double_t cosTerm = 0;
	  Double_t jthWeight = 0;
	  Double_t jthPhi = 0;
	  Int_t centID = v_events.at(i).centID;

	  //if (centID < 4) continue;  // ONLY LOOKING AT CENTRALITY 60% AND LOWER

	  // v2 from EPD E
	  for (Int_t j = 0; j < v_events.at(i).nHitsEpdE; j++)  // Loop through the j number of EPD E hits
	    {
	      jthWeight = v_events.at(i).tileWeightsEpdE.at(j);
	      jthPhi = v_events.at(i).phiValuesEpdE.at(j);

	      Double_t newXn = v_events.at(i).XnEpdE - jthWeight * TMath::Cos(ORDER_M * jthPhi);   // For event i, remove the jth particle from event plane
	      Double_t newYn = v_events.at(i).YnEpdE - jthWeight * TMath::Sin(ORDER_M * jthPhi);
	      Double_t newPsi = TMath::ATan2(newYn, newXn) / ORDER_M;
	      //v_events.at(i).eventPlanesEpdE.push_back(newPsi);
	      h_psiEpdE_NoAuto->Fill(newPsi);

	      // Add contribution to v_n from the jth particle using the event plane that omits the jth particle:
	      p_vn_EpdE->Fill(centID, TMath::Cos(ORDER_N * (jthPhi - newPsi)));
	    }



	  Double_t phi = 0;
	  Double_t psi = v_events.at(i).psiEpdE;

	  // v2 from EPD F
	  for (UInt_t j = 0; j < v_events.at(i).phiValuesEpdF.size(); j++)
	    {
	      phi = v_events.at(i).phiValuesEpdF.at(j);

	      p_vn_EpdF->Fill(centID, TMath::Cos(ORDER_N * (phi - psi)));
	    }

	  // v2 from TPC B
	  for (UInt_t j = 0; j < v_events.at(i).phiValuesTpcB.size(); j++)
	    {
	      phi = v_events.at(i).phiValuesTpcB.at(j);

	      p_vn_TpcB->Fill(centID, TMath::Cos(ORDER_N * (phi - psi)));
	    }



	  // Pi+
	  for (UInt_t j = 0; j < v_events.at(i).phiValuesPionP.size(); j++)
	    {
	      phi = v_events.at(i).phiValuesPionP.at(j);

	      p_vn_pp->Fill(centID, TMath::Cos(ORDER_N * (phi - psi)));
	    }

	  // Pi-
	  for (UInt_t j = 0; j < v_events.at(i).phiValuesPionM.size(); j++)
	    {
	      phi = v_events.at(i).phiValuesPionM.at(j);

	      p_vn_pm->Fill(centID, TMath::Cos(ORDER_N * (phi - psi)));
	    }

	  // K+
	  for (UInt_t j = 0; j < v_events.at(i).phiValuesKaonP.size(); j++)
	    {
	      phi = v_events.at(i).phiValuesKaonP.at(j);

	      p_vn_kp->Fill(centID, TMath::Cos(ORDER_N * (phi - psi)));
	    }

	  // K-
	  for (UInt_t j = 0; j < v_events.at(i).phiValuesKaonM.size(); j++)
	    {
	      phi = v_events.at(i).phiValuesKaonM.at(j);

	      p_vn_km->Fill(centID, TMath::Cos(ORDER_N * (phi - psi)));
	    }

	  // Proton
	  for (UInt_t j = 0; j < v_events.at(i).phiValuesProton.size(); j++)
	    {
	      phi = v_events.at(i).phiValuesProton.at(j);

	      p_vn_pr->Fill(centID, TMath::Cos(ORDER_N * (phi - psi)));
	    }
	  //=========================================================
	  //            End RAW Flow Calculations
	  //=========================================================
	  */

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
      h2_v2ScanTpc->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]); 
      h2_v2ScanEpd->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_v2ScanEpdTpcA->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
      h2_v2ScanEpdTpcB->GetYaxis()->SetBinLabel(i, centralityBins[labelIndex]);
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
}//End FlowAnalyzer()





      /* REMOVING AUTO-CORRELATIONS */
      /*
	for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	{
	Double_t newXn = v_events.at(i).Xn - v_events.at(i).pTValues.at(j) * TMath::Cos(ORDER_M * v_events.at(i).phiValues.at(j));
	Double_t newYn = v_events.at(i).Yn - v_events.at(i).pTValues.at(j) * TMath::Sin(ORDER_M * v_events.at(i).phiValues.at(j));
	Double_t newPsi = TMath::ATan2(newYn, newXn) / ORDER_M;
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



	      /*
	      Bool_t pion   = (d_pT >= 1.0) ? 
		((d_m2 > -0.1) && (d_m2 < 0.1) && (d_TPCnSigmaPion > -2) && (d_TPCnSigmaPion < 2)) :
		((d_TPCnSigmaPion > -2) && (d_TPCnSigmaPion < 2));
	      Bool_t kaon   = (d_pT >= 1.0) ? 
 		((d_m2 > 0.15) && (d_m2 < 0.34) && (d_TPCnSigmaKaon > -2) && (d_TPCnSigmaKaon < 2)) :
		((d_TPCnSigmaKaon > -2) && (d_TPCnSigmaKaon < 2));
	      Bool_t proton = (d_pT >= 1.0) ? 
		((d_m2 > 0.65) && (d_m2 < 1.11) && (d_TPCnSigmaProton > -2) && (d_TPCnSigmaProton < 2)) :
		((d_TPCnSigmaProton > -2) && (d_TPCnSigmaProton < 2));
	      */
