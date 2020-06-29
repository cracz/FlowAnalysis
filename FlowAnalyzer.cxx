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
#include "TString.h"
#include "TSystem.h"

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


// Custom type to hold important info for every good event
struct Event
{
  bool badEvent;      // Flag for marking events to ignore
  Int_t nTracks;      // Number of GOOD tracks in the event
  Int_t primTracks;   // Number of primary tracks before track cuts (used for centrality)
  Int_t centID;
  Double_t Xn;
  Double_t Yn;
  Double_t psi;       // Overall EP angle without removing autocorrelations
  Double_t psi_delta;   // Correction factor
  std::vector<Double_t> phiValues;   // Azimuthal angles for all particles in the event
  std::vector<Double_t> etaValues;  // Eta values for all particles in the event
  std::vector<Double_t> pTValues;   // pT values for all particles in the event
  std::vector<Double_t> psiValues;  // EP angles after removing autocorrelations

  Int_t nTracksA;
  Double_t XnA;
  Double_t YnA;
  Double_t psiA;       // Subevent A EP angle without removing autocorrelations
  Double_t psiA_delta;
  std::vector<Double_t> phiValuesA;  // Subevent A azimuthal angles
  std::vector<Double_t> pTValuesA;   // Subevent A pT values
  std::vector<Double_t> psiValuesA; // Subevent A EP angles after removing autocorrelations

  Int_t nTracksB;
  Double_t XnB;
  Double_t YnB;
  Double_t psiB;       // Subevent B EP angle without removing autocorrelations  
  Double_t psiB_delta;
  std::vector<Double_t> phiValuesB;  // Subevent B azimuthal angles
  std::vector<Double_t> pTValuesB;   // Subevent B pT values
  std::vector<Double_t> psiValuesB; // Subevent B EP angles after removing autocorrelations

  Int_t nTracksC;
  Double_t XnC;
  Double_t YnC;
  Double_t psiC;       // Subevent C EP angle without removing autocorrelations  
  Double_t psiC_delta;
  std::vector<Double_t> phiValuesC;  // Subevent C azimuthal angles
  std::vector<Double_t> pTValuesC;   // Subevent C pT values
  std::vector<Double_t> psiValuesC; // Subevent C EP angles after removing autocorrelations

  void reset()
  {
    badEvent  = false;  //Reset all values in the struct to reuse
    nTracks = 0; 
    primTracks = 0;
    centID = -5;
    Xn      = 0;
    Yn      = 0;
    psi     = 0;       
    psi_delta = 0; 
    phiValues.clear();
    etaValues.clear();
    pTValues.clear(); 
    psiValues.clear();

    nTracksA = 0;
    XnA      = 0;
    YnA      = 0;
    psiA     = 0;       
    psiA_delta = 0; 
    phiValuesA.clear();
    pTValuesA.clear(); 
    psiValuesA.clear();

    nTracksB = 0;
    XnB      = 0;
    YnB      = 0;
    psiB     = 0;       
    psiB_delta = 0; 
    phiValuesB.clear();
    pTValuesB.clear(); 
    psiValuesB.clear();

    nTracksC = 0;
    XnC      = 0;
    YnC      = 0;
    psiC     = 0;       
    psiC_delta = 0; 
    phiValuesC.clear();
    pTValuesC.clear(); 
    psiValuesC.clear();
  }
};


TH1D* recenter(TH1D* rawDist, std::vector<Event> &events, TString component);
TH1D* recenterSub(TH1D* rawDist, std::vector<Event> &events, TString component, TString subEvent);
TH1D* flatten(std::vector<Event> &events, const Int_t terms, const Int_t order_n);
TH1D* flattenSub(std::vector<Event> &events, const Int_t terms, const Int_t order_n, TString subEvent);


void FlowAnalyzer(TString inFile, TString jobID)
{
  std::cout << "Initializing..." << std::endl;

  if (gSystem->AccessPathName(inFile)) { std::cout << "Error reading input file!" << std::endl; return;}

  //=========================================================
  //          Some Controls
  //=========================================================
  const Double_t order_n = 2.0;
  TString order_n_str; order_n_str.Form("%d", order_n);

  //const Double_t etaCut = -1.01;    // 2 subevent plane method

  const Double_t etaCut1 = -0.5;     // 3 subevent plane method
  const Double_t etaCut2 = -1.0;
  //=========================================================
  //          
  //=========================================================



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
  cout << "EpdHit Branch returned found = " << found << endl;
  picoChain->SetBranchAddress("EpdHit",&epdHits);



  StEpdEpFinder *epdEpFinder = new StEpdEpFinder(16, "StEpdEpFinderCorrectionHistograms_OUTPUT_"+jobID+".root", "StEpdEpFinderCorrectionHistograms_INPUT.root");
  epdEpFinder->SetEpdHitFormat(2);       // format=0/1/2 for StEpdHit/StMuEpdHit/StPicoEpdHit
  epdEpFinder->SetnMipThreshold(0.3);    // recommended, but variable
  epdEpFinder->SetMaxTileWeight(1);      // recommended, but variable


  
  // OUTPUT FILE
  TString outFile = jobID.Append(".picoDst.result.root");
  TFile *outputFile = new TFile(outFile,"recreate");
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

  TH1D *h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", 100, 190, 210);

  TH1D *h_pT    = new TH1D("h_pT", "Particle p_{T};p_{T} (GeV);Particles", 100, 0, 5);
  TH1D *h_eta   = new TH1D("h_eta", "Particle #eta;#eta;Particles", 600, -6, 1);
  TH1D *h_eta_A = new TH1D("h_eta_A", "Particle #eta (Group A);#eta;Particles", 600, -2, 1);
  TH1D *h_eta_B = new TH1D("h_eta_B", "Particle #eta (Group B);#eta;Particles", 600, -2, 1);
  TH1D *h_eta_C = new TH1D("h_eta_C", "Particle #eta (Group C);#eta;Particles", 600, -2, 1);

  TH1D *h_Xn   = new TH1D("h_Xn", "X_n Distribution;X_n;Events", 100, -25, 25);
  TH1D *h_Yn   = new TH1D("h_Yn", "Y_n Distribution;Y_n;Events", 100, -25, 25);
  TH1D *h_Xn_A = new TH1D("h_Xn_A", "X_n Distribution (sub A);X_n;Events", 100, -25, 25);
  TH1D *h_Yn_A = new TH1D("h_Yn_A", "Y_n Distribution (sub A);Y_n;Events", 100, -25, 25);
  TH1D *h_Xn_B = new TH1D("h_Xn_B", "X_n Distribution (sub B);X_n;Events", 100, -25, 25);
  TH1D *h_Yn_B = new TH1D("h_Yn_B", "Y_n Distribution (sub B);Y_n;Events", 100, -25, 25);
  TH1D *h_Xn_C = new TH1D("h_Xn_C", "X_n Distribution (sub C);X_n;Events", 100, -25, 25);
  TH1D *h_Yn_C = new TH1D("h_Yn_C", "Y_n Distribution (sub C);Y_n;Events", 100, -25, 25);
  TH1D *h_psi  = new TH1D("h_psi", "Event Plane Angles (order "+order_n_str+");#psi_{"+order_n_str+"};Events", 400, -4, 4);
  TH1D *h_psi_A= new TH1D("h_psi_A", "Event Plane Angles (n = "+order_n_str+", sub A);#psi_{"+order_n_str+"};Events", 400, -4, 4);
  TH1D *h_psi_B= new TH1D("h_psi_B", "Event Plane Angles (n = "+order_n_str+", sub B);#psi_{"+order_n_str+"};Events", 400, -4, 4);
  TH1D *h_psi_C= new TH1D("h_psi_C", "Event Plane Angles (n = "+order_n_str+", sub C);#psi_{"+order_n_str+"};Events", 400, -4, 4);
  TH1D *h_v2Plot = new TH1D("h_v2Plot", "Plot to Retrieve v_{2};cos(2(#phi - #psi_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);

  TH1D *h_abCorr_n2 = new TH1D("h_abCorr_n2", "2nd Order Subevent a-b Correlations;cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_acCorr_n2 = new TH1D("h_acCorr_n2", "2nd Order Subevent a-c Correlations;cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{c}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_bcCorr_n2 = new TH1D("h_bcCorr_n2", "2nd Order Subevent b-c Correlations;cos(2(#psi^{b}_{"+order_n_str+"} - #psi^{c}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);

  /*
  TH1D *h_resol2 = new TH1D("h_resol2", "Plot to Retrieve 2nd Order Event Plane Resolution;cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent00 = new TH1D("h_resol2_cent00", "Plot to Retrieve 2nd Order Event Plane Resolution (75%-80% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent01 = new TH1D("h_resol2_cent01", "Plot to Retrieve 2nd Order Event Plane Resolution (70%-75% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent02 = new TH1D("h_resol2_cent02", "Plot to Retrieve 2nd Order Event Plane Resolution (65%-70% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent03 = new TH1D("h_resol2_cent03", "Plot to Retrieve 2nd Order Event Plane Resolution (60%-65% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent04 = new TH1D("h_resol2_cent04", "Plot to Retrieve 2nd Order Event Plane Resolution (55%-60% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent05 = new TH1D("h_resol2_cent05", "Plot to Retrieve 2nd Order Event Plane Resolution (50%-55% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent06 = new TH1D("h_resol2_cent06", "Plot to Retrieve 2nd Order Event Plane Resolution (45%-50% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent07 = new TH1D("h_resol2_cent07", "Plot to Retrieve 2nd Order Event Plane Resolution (40%-45% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent08 = new TH1D("h_resol2_cent08", "Plot to Retrieve 2nd Order Event Plane Resolution (35%-40% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent09 = new TH1D("h_resol2_cent09", "Plot to Retrieve 2nd Order Event Plane Resolution (30%-35% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent10 = new TH1D("h_resol2_cent10", "Plot to Retrieve 2nd Order Event Plane Resolution (25%-30% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent11 = new TH1D("h_resol2_cent11", "Plot to Retrieve 2nd Order Event Plane Resolution (20%-25% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent12 = new TH1D("h_resol2_cent12", "Plot to Retrieve 2nd Order Event Plane Resolution (15%-20% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent13 = new TH1D("h_resol2_cent13", "Plot to Retrieve 2nd Order Event Plane Resolution (10%-15% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent14 = new TH1D("h_resol2_cent14", "Plot to Retrieve 2nd Order Event Plane Resolution (5%-10% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  TH1D *h_resol2_cent15 = new TH1D("h_resol2_cent15", "Plot to Retrieve 2nd Order Event Plane Resolution (0%-5% cent.);cos(2(#psi^{a}_{"+order_n_str+"} - #psi^{b}_{"+order_n_str+"}));Particles", 150, -1.5, 1.5);
  */

  TH2D *h2_beta_p = new TH2D("h2_beta_p","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_m2_p   = new TH2D("h2_m2_p", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 500, -3, 3, 500, -0.1, 15);
  TH2D *h2_trans_vtx = new TH2D("h2_trans_vtx","Transverse Position of Primary Vertex;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);


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

      h_zvtx->Fill(d_zvtx);

      Bool_t b_bad_xvtx = ( (d_xvtx < -1.0) || (d_xvtx > 1.0)  );
      Bool_t b_bad_yvtx = ( (d_yvtx < -2.6) || (d_yvtx > -1.4) );
      Bool_t b_bad_zvtx = ( (d_zvtx < 200)  || (d_zvtx > 201.5));

      if (b_bad_zvtx) continue;

      h2_trans_vtx->Fill(d_xvtx, d_yvtx);

      if (b_bad_xvtx || b_bad_yvtx) continue;
      //=========================================================
      //      END Z-VTX Selection
      //=========================================================

      h_event_check->Fill(eventSections[2], 1);


      Int_t nTracks = dst->numberOfTracks();
      if (nTracks <= 5) continue;                // Require > 5 tracks in each event (This is repeated below after the cuts)

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

	  if (b_bad_hits || b_bad_dEdx || b_bad_tracking) continue;

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
	      eventInfo.nTracks++;
	      eventInfo.phiValues.push_back(d_phi);
	      eventInfo.pTValues.push_back(d_pT);
	      eventInfo.etaValues.push_back(d_eta);

	      h_pT->Fill(d_pT);
	      h_eta->Fill(d_eta);
	      //h2_beta_p->Fill(d_charge*d_mom, 1/d_tofBeta);
	      //h2_m2_p->Fill(d_charge*d_mom, d_mom*d_mom*(1/(d_tofBeta*d_tofBeta) - 1));

	      eventInfo.Xn += d_pT * TMath::Cos(order_n * d_phi);
	      eventInfo.Yn += d_pT * TMath::Sin(order_n * d_phi);
	    }
	}//End track loop


      for (UInt_t i = 0; i < eventInfo.etaValues.size(); i++)
	{
	  if (eventInfo.etaValues.at(i) < 0 && eventInfo.etaValues.at(i) >= etaCut1)
	    {
	      eventInfo.phiValuesA.push_back(eventInfo.phiValues.at(i));
	      eventInfo.pTValuesA.push_back(eventInfo.pTValues.at(i));

	      eventInfo.XnA += eventInfo.pTValues.at(i) * TMath::Cos(order_n * eventInfo.phiValues.at(i));
	      eventInfo.YnA += eventInfo.pTValues.at(i) * TMath::Sin(order_n * eventInfo.phiValues.at(i));

	      h_eta_A->Fill(eventInfo.etaValues.at(i));
	      eventInfo.nTracksA++;
	    }
	  else if (eventInfo.etaValues.at(i) < etaCut1 && eventInfo.etaValues.at(i) >= etaCut2)
	    {
	      eventInfo.phiValuesB.push_back(eventInfo.phiValues.at(i));
	      eventInfo.pTValuesB.push_back(eventInfo.pTValues.at(i));

	      eventInfo.XnB += eventInfo.pTValues.at(i) * TMath::Cos(order_n * eventInfo.phiValues.at(i));
	      eventInfo.YnB += eventInfo.pTValues.at(i) * TMath::Sin(order_n * eventInfo.phiValues.at(i));

	      h_eta_B->Fill(eventInfo.etaValues.at(i));
	      eventInfo.nTracksB++;
	    }
	  else if (eventInfo.etaValues.at(i) < etaCut2 && eventInfo.etaValues.at(i) >= -1.5)
	    {
	      eventInfo.phiValuesC.push_back(eventInfo.phiValues.at(i));
	      eventInfo.pTValuesC.push_back(eventInfo.pTValues.at(i));

	      eventInfo.XnC += eventInfo.pTValues.at(i) * TMath::Cos(order_n * eventInfo.phiValues.at(i));
	      eventInfo.YnC += eventInfo.pTValues.at(i) * TMath::Sin(order_n * eventInfo.phiValues.at(i));

	      h_eta_C->Fill(eventInfo.etaValues.at(i));
	      eventInfo.nTracksC++;
	    }
	}

      if (eventInfo.nTracks < 18) continue;                  // Make sure there are at least 6 GOOD tracks in each subevent
      if (eventInfo.nTracksA < 6) continue;
      if (eventInfo.nTracksB < 6) continue;
      if (eventInfo.nTracksC < 6) continue;

      eventInfo.psi  = TMath::ATan2(eventInfo.Yn,  eventInfo.Xn)  / order_n;
      eventInfo.psiA = TMath::ATan2(eventInfo.YnA, eventInfo.XnA) / order_n;
      eventInfo.psiB = TMath::ATan2(eventInfo.YnB, eventInfo.XnB) / order_n;
      eventInfo.psiC = TMath::ATan2(eventInfo.YnC, eventInfo.XnC) / order_n;

      h_Xn->Fill(eventInfo.Xn);
      h_Yn->Fill(eventInfo.Yn);
      h_psi->Fill(eventInfo.psi);

      h_Xn_A->Fill(eventInfo.XnA);
      h_Yn_A->Fill(eventInfo.YnA);
      h_psi_A->Fill(eventInfo.psiA);

      h_Xn_B->Fill(eventInfo.XnB);
      h_Yn_B->Fill(eventInfo.YnB);
      h_psi_B->Fill(eventInfo.psiB);

      h_Xn_C->Fill(eventInfo.XnC);
      h_Yn_C->Fill(eventInfo.YnC);
      h_psi_C->Fill(eventInfo.psiC);


      // ASSIGN CENTRALITY ID
      if     ( eventInfo.primTracks >=   3 && eventInfo.primTracks <=   4 ) eventInfo.centID =  0;
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
      else if( eventInfo.primTracks >= 141 && eventInfo.primTracks <= 195 ) eventInfo.centID = 15;

      if (eventInfo.centID == -5) continue;


      //=========================================================
      //                EPD STUFF
      //=========================================================
      StEpdEpInfo result = epdEpFinder->Results(epdHits,pVtx,eventInfo.centID);
      StEpdGeom *epdGeom = new StEpdGeom();
      StPicoEpdHit *epdHit;
      TVector3 tileVector;     // Vector from vertex to center of tile that was hit

      for (int iEpdHit = 0; iEpdHit < epdHits->GetEntries(); iEpdHit++)
	{
	  epdHit = (StPicoEpdHit*)(epdHits->At(iEpdHit));
	  int tileID = epdHit->id();
	  tileVector = epdGeom->TileCenter(tileID) - pVtx;

	  double tileEta = tileVector.Eta();

	  h_eta->Fill(tileEta);
	}

      delete epdGeom;
      //=========================================================
      //            END EPD STUFF
      //=========================================================


      v_events.push_back(eventInfo);   // Store this event with all of its attributes
    }//End event loop

  TH1D *h_xvtx = h2_trans_vtx->ProjectionX();
  TH1D *h_yvtx = h2_trans_vtx->ProjectionY();

  Int_t numOfEvents = v_events.size();   // Number of events that passed the cuts

  //=========================================================
  //          Re-centering (Xn, Yn) Distributions
  //=========================================================

  TH1D* h_Xn_s = recenter(h_Xn, v_events, "X");   // Shifted distributions. Vectors are now updated with the recentered values
  TH1D* h_Yn_s = recenter(h_Yn, v_events, "Y");
  TH1D* h_Xn_A_s = recenterSub(h_Xn_A, v_events, "X", "A");
  TH1D* h_Yn_A_s = recenterSub(h_Yn_A, v_events, "Y", "A");
  TH1D* h_Xn_B_s = recenterSub(h_Xn_B, v_events, "X", "B");
  TH1D* h_Yn_B_s = recenterSub(h_Yn_B, v_events, "Y", "B");
  TH1D* h_Xn_C_s = recenterSub(h_Xn_C, v_events, "X", "C");
  TH1D* h_Yn_C_s = recenterSub(h_Yn_C, v_events, "Y", "C");

  Int_t badEvents = 0;

  for (Int_t i = 0; i < numOfEvents; i++)         // Recalculating the event plane angles
    {
      if (v_events.at(i).Xn == 0 && v_events.at(i).Yn == 0) 
	{
	  v_events.at(i).badEvent = true;         // Mark bad events to ignore later
	  badEvents++;
	  continue;
	}
      else if (v_events.at(i).XnA == 0 && v_events.at(i).YnA == 0) 
	{
	  v_events.at(i).badEvent = true;
	  badEvents++;
	  continue;
	}
      else if (v_events.at(i).XnB == 0 && v_events.at(i).YnB == 0) 
	{
	  v_events.at(i).badEvent = true;
	  badEvents++;
	  continue;
	}
      else if (v_events.at(i).XnC == 0 && v_events.at(i).YnC == 0) 
	{
	  v_events.at(i).badEvent = true;
	  badEvents++;
	  continue;
	}
      else
	{ 
	  v_events.at(i).psi  = TMath::ATan2(v_events.at(i).Yn,  v_events.at(i).Xn)  / order_n; 
	  v_events.at(i).psiA = TMath::ATan2(v_events.at(i).YnA, v_events.at(i).XnA) / order_n; 
	  v_events.at(i).psiB = TMath::ATan2(v_events.at(i).YnB, v_events.at(i).XnB) / order_n; 
	  v_events.at(i).psiC = TMath::ATan2(v_events.at(i).YnC, v_events.at(i).XnC) / order_n; 

	  /* REMOVING AUTO-CORRELATIONS */

	  for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	    {
	      Double_t newXn = v_events.at(i).Xn - v_events.at(i).pTValues.at(j) * TMath::Cos(order_n * v_events.at(i).phiValues.at(j));
	      Double_t newYn = v_events.at(i).Yn - v_events.at(i).pTValues.at(j) * TMath::Sin(order_n * v_events.at(i).phiValues.at(j));
	      Double_t newPsi = TMath::ATan2(newYn, newXn) / order_n;
	      v_events.at(i).psiValues.push_back(newPsi);
	    }	  
	}
    }

  std::cout << "Bad Events after recentering: " << badEvents << std::endl;

  //=========================================================
  //          Flattening Event Plane Angle
  //=========================================================

  TH1D *h_psi_flat  = flatten(v_events, 20, order_n);
  TH1D *h_psiA_flat = flattenSub(v_events, 20, order_n, "A");
  TH1D *h_psiB_flat = flattenSub(v_events, 20, order_n, "B");
  TH1D *h_psiC_flat = flattenSub(v_events, 20, order_n, "C");

  
  //=========================================================
  //          Use A Histogram To Get v_2
  //=========================================================

  Double_t cosTerm = 0;
  Double_t phi = 0;
  Double_t psi = 0;
  Double_t abCorrTerm_n2 = 0;
  Double_t acCorrTerm_n2 = 0;
  Double_t bcCorrTerm_n2 = 0;
  Double_t psiA = 0;
  Double_t psiB = 0;
  Double_t psiC = 0;

  // CORRELATIONS AND FLOW COEFFICIENTS
  for (Int_t i = 0; i < numOfEvents; i++)
    {
      // SUBEVENT CORRELATIONS
      psiA = v_events.at(i).psiA;
      psiB = v_events.at(i).psiB;
      psiC = v_events.at(i).psiC;

      abCorrTerm_n2 = TMath::Cos(2 * (psiA - psiB));
      acCorrTerm_n2 = TMath::Cos(2 * (psiA - psiC));
      bcCorrTerm_n2 = TMath::Cos(2 * (psiB - psiC));
      h_abCorr_n2->Fill(abCorrTerm_n2);
      h_acCorr_n2->Fill(acCorrTerm_n2);
      h_bcCorr_n2->Fill(bcCorrTerm_n2);

      /*
      switch (v_events.at(i).centID)
	{
	case 0: h_resol2_cent00->Fill(corrTerm); break;
	case 1: h_resol2_cent01->Fill(corrTerm); break;
	case 2: h_resol2_cent02->Fill(corrTerm); break;
	case 3: h_resol2_cent03->Fill(corrTerm); break;
	case 4: h_resol2_cent04->Fill(corrTerm); break;
	case 5: h_resol2_cent05->Fill(corrTerm); break;
	case 6: h_resol2_cent06->Fill(corrTerm); break;
	case 7: h_resol2_cent07->Fill(corrTerm); break;
	case 8: h_resol2_cent08->Fill(corrTerm); break;
	case 9: h_resol2_cent09->Fill(corrTerm); break;
	case 10: h_resol2_cent10->Fill(corrTerm); break;
	case 11: h_resol2_cent11->Fill(corrTerm); break;
	case 12: h_resol2_cent12->Fill(corrTerm); break;
	case 13: h_resol2_cent13->Fill(corrTerm); break;
	case 14: h_resol2_cent14->Fill(corrTerm); break;
	case 15: h_resol2_cent15->Fill(corrTerm); break;
	}
      */

      // GET V2 WITHOUT AUTOCORRELATIONS HERE
      for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	{
	  phi = v_events.at(i).phiValues.at(j);
	  psi = v_events.at(i).psiValues.at(j);  // This psi was calculated without particle 'j'
	  
	  cosTerm = TMath::Cos(2 * (phi - psi));
	  h_v2Plot->Fill(cosTerm);
	}
    }

  outputFile->cd();

  h_event_check ->Write();
  h_track_check ->Write();
  h_nhits       ->Write();
  h_nhits_fit   ->Write();
  h_nhits_dEdx  ->Write();
  h_xvtx        ->Write();
  h_yvtx        ->Write();
  h_zvtx        ->Write();
  h_pT          ->Write();
  h_eta         ->Write();
  h_eta_A       ->Write();
  h_eta_B       ->Write();
  h_eta_C       ->Write();
  h_Xn          ->Write();
  h_Yn          ->Write();
  h_Xn_A        ->Write();
  h_Yn_A        ->Write();
  h_Xn_B        ->Write();
  h_Yn_B        ->Write();
  h_Xn_C        ->Write();
  h_Yn_C        ->Write();
  h_Xn_s        ->Write();
  h_Yn_s        ->Write();  
  h_Xn_A_s      ->Write();
  h_Yn_A_s      ->Write();
  h_Xn_B_s      ->Write();
  h_Yn_B_s      ->Write();
  h_Xn_C_s      ->Write();
  h_Yn_C_s      ->Write();
  h_psi         ->Write();
  h_psi_A       ->Write();
  h_psi_B       ->Write();
  h_psi_C       ->Write();
  h_psi_flat    ->Write();
  h_psiA_flat   ->Write();
  h_psiB_flat   ->Write();
  h_psiC_flat   ->Write();
  h_v2Plot      ->Write();
  h_abCorr_n2   ->Write();
  h_acCorr_n2   ->Write();
  h_bcCorr_n2   ->Write();

  /*
  h_resol2      ->Write();
  h_resol2_cent00      ->Write();
  h_resol2_cent01      ->Write();
  h_resol2_cent02      ->Write();
  h_resol2_cent03      ->Write();
  h_resol2_cent04      ->Write();
  h_resol2_cent05      ->Write();
  h_resol2_cent06      ->Write();
  h_resol2_cent07      ->Write();
  h_resol2_cent08      ->Write();
  h_resol2_cent09      ->Write();
  h_resol2_cent10      ->Write();
  h_resol2_cent11      ->Write();
  h_resol2_cent12      ->Write();
  h_resol2_cent13      ->Write();
  h_resol2_cent14      ->Write();
  h_resol2_cent15      ->Write();
  */

  h2_beta_p     ->Write();
  h2_m2_p       ->Write();
  h2_trans_vtx  ->Write();

  gROOT->GetListOfFiles()->Remove(outputFile);
  outputFile->Close();

  picoReader->Finish();

  std::cout << "Done!" << std::endl;
}//End EventPlane



////
// This function goes through the (Xn,Yn) values in the events vector 
//and subtracts the average value found from the raw distributions. 
//The (Xn,Yn) values are also updated within the events vector.
////
TH1D* recenter(TH1D* rawDist, std::vector<Event> &events, TString component)
{
  // Copy histogram parameters from the raw distribution
  TString name   = rawDist->GetName(); name = name + "_s";
  TString title  = rawDist->GetTitle(); title = title + " (recentered)";
  TString xLabel = rawDist->GetXaxis()->GetTitle();
  TString yLabel = rawDist->GetYaxis()->GetTitle();
  TString fullTitle = title + ";" + xLabel + ";" + yLabel;

  Int_t bins    = rawDist->GetNbinsX();
  Double_t min  = rawDist->GetXaxis()->GetXmin();
  Double_t max  = rawDist->GetXaxis()->GetXmax();
  //Double_t avg  = rawDist->GetMean();

  TH1D* h_shifted = new TH1D(name, fullTitle, bins, min, max);

  Int_t numOfEvents = events.size();

  // Loop to shift (Xn,Yn) distributions by average values and fill shifted histogram
  if (component.EqualTo("x") || component.EqualTo("X")) 
    {
      for (Int_t i = 0; i < numOfEvents; i++)
	{
	  //events.at(i).Xn -= -0.17149495475033577;             // 2 subevents
	  //events.at(i).Xn -= -0.2906029561873486;              // 3 subevents
	  events.at(i).Xn -= -0.28449056558920105;               // 3 subevents, pT cut
	  //events.at(i).Xn -=                // 3 subevents, pT cut, no TOF beta cut
	  h_shifted->Fill(events.at(i).Xn);
	}
    }
  else if (component.EqualTo("y") || component.EqualTo("Y")) 
    {
      for (Int_t i = 0; i < numOfEvents; i++)
	{
	  //events.at(i).Yn -= -0.12203934435287563;              // 2 subevents
	  //events.at(i).Yn -= -0.20888164534551301;              // 3 subevents
	  events.at(i).Yn -= -0.20357347770740941;                // 3 subevents, pT cut
	  //events.at(i).Yn -=                 // 3 subevents, pT cut, no TOF beta cut
	  h_shifted->Fill(events.at(i).Yn);
	}
    }

  return h_shifted;
}//End recenter()



TH1D* recenterSub(TH1D* rawDist, std::vector<Event> &events, TString component, TString subEvent)
{
  // Copy histogram parameters from the raw distribution
  TString name   = rawDist->GetName(); name = name + "_s";
  TString title  = rawDist->GetTitle(); title = title + " (recentered)";
  TString xLabel = rawDist->GetXaxis()->GetTitle();
  TString yLabel = rawDist->GetYaxis()->GetTitle();
  TString fullTitle = title + ";" + xLabel + ";" + yLabel;

  Int_t bins    = rawDist->GetNbinsX();
  Double_t min  = rawDist->GetXaxis()->GetXmin();
  Double_t max  = rawDist->GetXaxis()->GetXmax();

  Int_t numOfEvents = events.size();

  TH1D* h_shifted = new TH1D(name, fullTitle, bins, min, max);

  // Loop to shift (Xn,Yn) distributions by average values and fill shifted histogram
  if (subEvent.EqualTo("a") || subEvent.EqualTo("A"))
    {
      if (component.EqualTo("x") || component.EqualTo("X")) 
	{
	  for (Int_t i = 0; i < numOfEvents; i++)
	    {
	      //events.at(i).XnA -= -0.07639361131478721;    // Auto eta cuts
	      //events.at(i).XnA -= -0.082509858448988699;    // Straight eta cut, 2 subevents
	      //events.at(i).XnA -= -0.048029071218188758;    // 3 subevents
	      events.at(i).XnA -= -0.047839796373211604;     // 3 subevents, pT cut
	      //events.at(i).XnA -=      // 3 subevents, pT cut, no TOF beta cut
	      h_shifted->Fill(events.at(i).XnA);
	    }
	}      
      else if (component.EqualTo("y") || component.EqualTo("Y")) 
	{
	  for (Int_t i = 0; i < numOfEvents; i++)
	    {
	      //events.at(i).YnA -= 0.052665569322193691;    // Auto eta cuts
	      //events.at(i).YnA -= 0.05781093198280992;    // Straight eta cut, 2 subevents
	      //events.at(i).YnA -= -0.15330719872348675;   // 3 subevents
	      events.at(i).YnA -= -0.15320538376331011;     // 3 subevents, pT cut
	      //events.at(i).YnA -=      // 3 subevents, pT cut, no TOF beta cut
	      h_shifted->Fill(events.at(i).YnA);
	    }
	}
    }
  else if (subEvent.EqualTo("b") || subEvent.EqualTo("B"))
    {
      if (component.EqualTo("x") || component.EqualTo("X")) 
	{
	  for (Int_t i = 0; i < numOfEvents; i++)
	    {
	      //events.at(i).XnB -= -0.091769164703734307;    // Auto eta cuts
	      //events.at(i).XnB -= -0.08898482140291171;   // Straight eta cut, 2 subevents
	      //events.at(i).XnB -= -0.094491708306180172;  // 3 subevents
	      events.at(i).XnB -= -0.092434588910897383;    // 3 subevents, pT cut
	      //events.at(i).XnB -=     // 3 subevents, pT cut, no TOF beta cut
	      h_shifted->Fill(events.at(i).XnB);
	    }
	}      
      else if (component.EqualTo("y") || component.EqualTo("Y")) 
	{
	  for (Int_t i = 0; i < numOfEvents; i++)
	    {
	      //events.at(i).YnB -= -0.17590862668749332;    // Auto eta cuts
	      //events.at(i).YnB -= -0.1798503754306236;    // Straight eta cut, 2 subevents
	      //events.at(i).YnB -= -0.14915417370216158;   // 3 subevents
	      events.at(i).YnB -= -0.14804355679219941;     // 3 subevents, pT cut
	      //events.at(i).YnB -=      // 3 subevents, pT cut, no TOF beta cut
	      h_shifted->Fill(events.at(i).YnB);
	    }
	}
    }
  else if (subEvent.EqualTo("c") || subEvent.EqualTo("C"))
    {
      if (component.EqualTo("x") || component.EqualTo("X")) 
	{
	  for (Int_t i = 0; i < numOfEvents; i++)
	    {
	      //events.at(i).XnC -= -0.14808192683948013;  // 3 subevents
	      events.at(i).XnC -= -0.14421595436210466;     // 3 subevents, pT cut
	      //events.at(i).XnC -=      // 3 subevents, pT cut, no TOF beta cut
	      h_shifted->Fill(events.at(i).XnC);
	    }
	}      
      else if (component.EqualTo("y") || component.EqualTo("Y")) 
	{
	  for (Int_t i = 0; i < numOfEvents; i++)
	    {
	      //events.at(i).YnC -= 0.093576529699855707;  // 3 subevents
	      events.at(i).YnC -= 0.097671928736901739;     // 3 subevents, pT cut
	      //events.at(i).YnC -=      // 3 subevents, pT cut, no TOF beta cut
	      h_shifted->Fill(events.at(i).YnC);
	    }
	}
    }

  return h_shifted;
}//End recenterSub()



////
// This function applies flattening corrections to event plane angle values.
//The number of correction terms is determined by the "terms" parameter, and 
//the "order_n" is just for the axis label on the histogram. The correction
//term is saved in each event, and the event plane angle is updated
//in each event.
////
TH1D* flatten(std::vector<Event> &events, const Int_t terms, const Int_t order_n)
{
  TString order_n_str; order_n_str.Form("%d", order_n);

  TH1D *h_psi_flat = new TH1D("h_psi_flat", "Flattened Event Plane Angle (order "+order_n_str+");#psi_{"+order_n_str+"};Events", 400, -4, 4);

  // Get actual number of good events //
  Int_t numOfEvents = events.size();
  Int_t numOfGoodEvents = 0;
  Int_t numOfBadEvents  = 0;

  for (Int_t i = 0; i < numOfEvents; i++)
    { if (events.at(i).badEvent == true) { numOfBadEvents++; } }

  numOfGoodEvents = numOfEvents - numOfBadEvents;
  ////


  // Get averages //
  Double_t sin_psi_avgs[terms];      // Arrays that contain each average parameter for the correction terms
  Double_t cos_psi_avgs[terms];

  for (Int_t i = 0; i < terms; i++)  // Loop through each term in the correction to build averages
    {
      sin_psi_avgs[i] = 0;  // elements of variable length arrays must be initialized manually
      cos_psi_avgs[i] = 0;

      Int_t j = i + 1;      // index for the sum in deltaPsi 

      for (Int_t k = 0; k < numOfEvents; k++)
	{
	  if (events.at(k).badEvent == true) { continue; }

	  sin_psi_avgs[i] += TMath::Sin(j * order_n * events.at(k).psi);     // First getting the sums for the averages
	  cos_psi_avgs[i] += TMath::Cos(j * order_n * events.at(k).psi);
	}

      sin_psi_avgs[i] = sin_psi_avgs[i]/numOfGoodEvents;
      cos_psi_avgs[i] = cos_psi_avgs[i]/numOfGoodEvents;
    }
  ////


  // Get corrected event plane angles //
  for (Int_t k = 0; k < numOfEvents; k++)  // Loop over psi values to correct them
    {
      if (events.at(k).badEvent == true) { continue; }

      for (Int_t i = 0; i < terms; i++)    // Build the correction term
	{
	  Int_t j = i + 1;

	  events.at(k).psi_delta += ((Double_t)2/(j*order_n))*(-sin_psi_avgs[i]*TMath::Cos(j * order_n * events.at(k).psi) 
								 + cos_psi_avgs[i]*TMath::Sin(j * order_n * events.at(k).psi));
	}

      events.at(k).psi += events.at(k).psi_delta;

      if (events.at(k).psi > TMath::Pi()/order_n) { events.at(k).psi -= TMath::TwoPi()/order_n; }        // Must maintain the angles' periodicity of 2pi/n
      else if (events.at(k).psi < -TMath::Pi()/order_n) { events.at(k).psi += TMath::TwoPi()/order_n; }

      h_psi_flat->Fill(events.at(k).psi);
    }
  ////

  return h_psi_flat;
}// End flatten()



TH1D* flattenSub(std::vector<Event> &events, const Int_t terms, const Int_t order_n, TString subEvent)
{
  TString order_n_str; order_n_str.Form("%d", order_n);

  TH1D *h_psi_flat = new TH1D("h_psi"+subEvent+"_flat", "Flattened Sub-Event Plane "+subEvent+" Angle (order "+order_n_str+");#psi_{"+order_n_str+"};Events", 400, -4, 4);

  // Get actual number of good events //
  Int_t numOfEvents = events.size();
  Int_t numOfGoodEvents = 0;
  Int_t numOfBadEvents  = 0;

  for (Int_t i = 0; i < numOfEvents; i++)
    { if (events.at(i).badEvent == true) { numOfBadEvents++; } }

  numOfGoodEvents = numOfEvents - numOfBadEvents;
  ////

  if (subEvent.EqualTo("a") || subEvent.EqualTo("A"))
    {
      // Get averages //
      Double_t sin_psiA_avgs[terms];      // Arrays that contain each average parameter for the correction terms
      Double_t cos_psiA_avgs[terms];

      for (Int_t i = 0; i < terms; i++)  // Loop through each term in the correction to build averages
	{
	  sin_psiA_avgs[i] = 0;  // elements of variable length arrays must be initialized manually
	  cos_psiA_avgs[i] = 0;

	  Int_t j = i + 1;      // index for the sum in deltaPsi 

	  for (Int_t k = 0; k < numOfEvents; k++)
	    {
	      if (events.at(k).badEvent == true) { continue; }

	      sin_psiA_avgs[i] += TMath::Sin(j * order_n * events.at(k).psiA);     // First getting the sums for the averages
	      cos_psiA_avgs[i] += TMath::Cos(j * order_n * events.at(k).psiA);
	    }

	  sin_psiA_avgs[i] = sin_psiA_avgs[i]/numOfGoodEvents;
	  cos_psiA_avgs[i] = cos_psiA_avgs[i]/numOfGoodEvents;
	}
      ////


      // Get corrected event plane angles //
      for (Int_t k = 0; k < numOfEvents; k++)  // Loop over psi values to correct them
	{
	  if (events.at(k).badEvent == true) { continue; }

	  for (Int_t i = 0; i < terms; i++)    // Build the correction term
	    {
	      Int_t j = i + 1;

	      events.at(k).psiA_delta += ((Double_t)2/(j*order_n))*(-sin_psiA_avgs[i]*TMath::Cos(j * order_n * events.at(k).psiA) 
								   + cos_psiA_avgs[i]*TMath::Sin(j * order_n * events.at(k).psiA));
	    }

	  events.at(k).psiA += events.at(k).psiA_delta;

	  if (events.at(k).psiA > TMath::Pi()/order_n) { events.at(k).psiA -= TMath::TwoPi()/order_n; }        // Must maintain the angles' periodicity of 2pi/n
	  else if (events.at(k).psiA < -TMath::Pi()/order_n) { events.at(k).psiA += TMath::TwoPi()/order_n; }

	  h_psi_flat->Fill(events.at(k).psiA);
	}
      ////
    }  
  else if (subEvent.EqualTo("b") || subEvent.EqualTo("B"))  
    {
      // Get averages //
      Double_t sin_psiB_avgs[terms];      // Arrays that contain each average parameter for the correction terms
      Double_t cos_psiB_avgs[terms];

      for (Int_t i = 0; i < terms; i++)  // Loop through each term in the correction to build averages
	{
	  sin_psiB_avgs[i] = 0;  // elements of variable length arrays must be initialized manually
	  cos_psiB_avgs[i] = 0;

	  Int_t j = i + 1;      // index for the sum in deltaPsi 

	  for (Int_t k = 0; k < numOfEvents; k++)
	    {
	      if (events.at(k).badEvent == true) { continue; }

	      sin_psiB_avgs[i] += TMath::Sin(j * order_n * events.at(k).psiB);     // First getting the sums for the averages
	      cos_psiB_avgs[i] += TMath::Cos(j * order_n * events.at(k).psiB);
	    }

	  sin_psiB_avgs[i] = sin_psiB_avgs[i]/numOfGoodEvents;
	  cos_psiB_avgs[i] = cos_psiB_avgs[i]/numOfGoodEvents;
	}
      ////


      // Get corrected event plane angles //
      for (Int_t k = 0; k < numOfEvents; k++)  // Loop over psi values to correct them
	{
	  if (events.at(k).badEvent == true) { continue; }

	  for (Int_t i = 0; i < terms; i++)    // Build the correction term
	    {
	      Int_t j = i + 1;

	      events.at(k).psiB_delta += ((Double_t)2/(j*order_n))*(-sin_psiB_avgs[i]*TMath::Cos(j * order_n * events.at(k).psiB) 
								   + cos_psiB_avgs[i]*TMath::Sin(j * order_n * events.at(k).psiB));
	    }

	  events.at(k).psiB += events.at(k).psiB_delta;

	  if (events.at(k).psiB > TMath::Pi()/order_n) { events.at(k).psiB -= TMath::TwoPi()/order_n; }        // Must maintain the angles' periodicity of 2pi/n
	  else if (events.at(k).psiB < -TMath::Pi()/order_n) { events.at(k).psiB += TMath::TwoPi()/order_n; }

	  h_psi_flat->Fill(events.at(k).psiB);
	}
      ////
    }
  else if (subEvent.EqualTo("c") || subEvent.EqualTo("C"))  
    {
      // Get averages //
      Double_t sin_psiC_avgs[terms];      // Arrays that contain each average parameter for the correction terms
      Double_t cos_psiC_avgs[terms];

      for (Int_t i = 0; i < terms; i++)  // Loop through each term in the correction to build averages
	{
	  sin_psiC_avgs[i] = 0;  // elements of variable length arrays must be initialized manually
	  cos_psiC_avgs[i] = 0;

	  Int_t j = i + 1;      // index for the sum in deltaPsi 

	  for (Int_t k = 0; k < numOfEvents; k++)
	    {
	      if (events.at(k).badEvent == true) { continue; }

	      sin_psiC_avgs[i] += TMath::Sin(j * order_n * events.at(k).psiC);     // First getting the sums for the averages
	      cos_psiC_avgs[i] += TMath::Cos(j * order_n * events.at(k).psiC);
	    }

	  sin_psiC_avgs[i] = sin_psiC_avgs[i]/numOfGoodEvents;
	  cos_psiC_avgs[i] = cos_psiC_avgs[i]/numOfGoodEvents;
	}
      ////


      // Get corrected event plane angles //
      for (Int_t k = 0; k < numOfEvents; k++)  // Loop over psi values to correct them
	{
	  if (events.at(k).badEvent == true) { continue; }

	  for (Int_t i = 0; i < terms; i++)    // Build the correction term
	    {
	      Int_t j = i + 1;

	      events.at(k).psiC_delta += ((Double_t)2/(j*order_n))*(-sin_psiC_avgs[i]*TMath::Cos(j * order_n * events.at(k).psiC) 
								   + cos_psiC_avgs[i]*TMath::Sin(j * order_n * events.at(k).psiC));
	    }

	  events.at(k).psiC += events.at(k).psiC_delta;

	  if (events.at(k).psiC > TMath::Pi()/order_n) { events.at(k).psiC -= TMath::TwoPi()/order_n; }        // Must maintain the angles' periodicity of 2pi/n
	  else if (events.at(k).psiC < -TMath::Pi()/order_n) { events.at(k).psiC += TMath::TwoPi()/order_n; }

	  h_psi_flat->Fill(events.at(k).psiC);
	}
      ////
    }

  return h_psi_flat;
}// End flattenSub()
