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
  bool badEvent;  // Flag for marking events to ignore
  Int_t nTracks;  // Number of particles in the event
  Double_t Xn;
  Double_t Yn;
  Double_t psi;
  Double_t psi_delta;   // Correction factor
  std::vector<Double_t> phiValues;  // Azimuthal angles for all particles in the event
};


TH1D* recenter(TH1D* rawDist, std::vector<Event> &events, TString component);
TH1D* flatten(std::vector<Event> &events, const Int_t terms, const Int_t order_n);


void FlowAnalyzer(TString inFile, TString outFile)
{
  std::cout << "Initializing..." << std::endl;
  
  if (gSystem->AccessPathName(inFile)) { std::cout << "Error reading input file!" << std::endl; return;}

  //=========================================================
  //          Event Plane Order To Analyze
  //=========================================================
  Int_t order_n = 2;
  TString order_n_str; order_n_str.Form("%d", order_n);
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
  
  // OUTPUT FILE
  outFile.Append(".picoDst.result.root");
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

  TH1D *h_pT   = new TH1D("h_pT", "Particle p_{T};p_{T} (GeV);Particles", 100, 0, 5);
  TH1D *h_eta  = new TH1D("h_eta", "Particle #eta;#eta;Particles", 150, -2, 1);

  TH1D *h_Xn   = new TH1D("h_Xn", "X_n Distribution;X_n;Events", 100, -25, 25);
  TH1D *h_Yn   = new TH1D("h_Yn", "Y_n Distribution;Y_n;Events", 100, -25, 25);
  TH1D *h_psi  = new TH1D("h_psi", "Event Plane Angles (order "+order_n_str+");#psi_{"+order_n_str+"};Events", 400, -4, 4);
  TH1D *h_v2Plot = new TH1D("h_v2Plot", "Plot to retrieve v_{2};cos(2(#phi - #psi_{n}));Particles", 150, -1.5, 1.5);

  TH2D *h2_beta_p = new TH2D("h2_beta_p","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  TH2D *h2_m2_p   = new TH2D("h2_m2_p", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 500, -3, 3, 500, -0.1, 15);
  TH2D *h2_trans_vtx = new TH2D("h2_trans_vtx","Transverse Position of Primary Vertex;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);


  Event eventInfo;
  std::vector<Event> v_events;  // Vector of all events and their info using the custom struct

  std::vector<UInt_t> triggerIDs;

  // EVENT LOOP
  for (Long64_t ievent = 0; ievent < events2read; ievent++)
    {
      eventInfo.badEvent = false;  //Reset all values in the eventInfo struct to reuse
      eventInfo.nTracks = 0;
      eventInfo.Xn = 0;
      eventInfo.Yn = 0;
      eventInfo.psi = 0;
      eventInfo.psi_delta = 0;
      eventInfo.phiValues.clear();


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

      for (UInt_t i = 0; i < triggerIDs.size(); i++) { if (triggerIDs[i] == 620052) {b_bad_trig = false;} }

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
      if (nTracks <= 5) continue;                // Require > 5 tracks in each event

      // TRACK LOOP OVER PRIMARY TRACKS
      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
	{            
	  StPicoTrack *picoTrack = dst->track(iTrk);            
	  if(picoTrack == NULL) continue;
	  if(!picoTrack->isPrimary()) continue;  // Require primary tracks

	  h_track_check->Fill(trackSections[0], 1);

	  //=========================================================
	  //          TOF Beta Cuts
	  //=========================================================

	  if (!picoTrack->isTofTrack()) continue;  // Only TOF tracks

	  StPicoBTofPidTraits *trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());

	  Double_t d_tofBeta = trait->btofBeta();

	  if (d_tofBeta < 0.01 /*|| d_tofBeta >= 1*/) continue;
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

	      h_pT->Fill(d_pT);
	      h_eta->Fill(d_eta);
	      h2_beta_p->Fill(d_charge*d_mom, 1/d_tofBeta);
	      h2_m2_p->Fill(d_charge*d_mom, d_mom*d_mom*(1/(d_tofBeta*d_tofBeta) - 1));

	      eventInfo.Xn += d_pT * TMath::Cos(order_n * d_phi);
	      eventInfo.Yn += d_pT * TMath::Sin(order_n * d_phi);
	    }
	}//End track loop
      
      if (eventInfo.Xn == 0 && eventInfo.Yn == 0) continue;    // Cut out events with no flow vector

      eventInfo.psi = TMath::ATan2(eventInfo.Yn, eventInfo.Xn) / order_n;
      
      h_psi->Fill(eventInfo.psi);
      h_Xn->Fill(eventInfo.Xn);
      h_Yn->Fill(eventInfo.Yn);

      v_events.push_back(eventInfo);
    }//End event loop

  TH1D *h_xvtx = h2_trans_vtx->ProjectionX();
  TH1D *h_yvtx = h2_trans_vtx->ProjectionY();

  Int_t numOfEvents = v_events.size();   // Number of events that passed the cuts

  //=========================================================
  //          Re-centering (Xn, Yn) Distributions
  //=========================================================

  TH1D* h_Xn_s = recenter(h_Xn, v_events, "x");   // Shifted distributions. Both vectors are now updated with the recentered values
  TH1D* h_Yn_s = recenter(h_Yn, v_events, "y");

  for (Int_t i = 0; i < numOfEvents; i++)         // Recalculating the event plane angles
    {
      if (v_events.at(i).Xn == 0 && v_events.at(i).Yn == 0) 
	{
	  v_events.at(i).badEvent = true;         // Mark bad events to ignore later
	  continue;
	}
      else
	{ v_events.at(i).psi = TMath::ATan2(v_events.at(i).Yn, v_events.at(i).Xn) / order_n; }
    }


  //=========================================================
  //          Flattening Event Plane Angle
  //=========================================================

  TH1D *h_psi_flat = flatten(v_events, 20, order_n);

  
  //=========================================================
  //          Use A Histogram To Get v_2
  //=========================================================

  Double_t cosTerm = 0;
  Double_t phi = 0;
  Double_t psi = 0;
  Int_t totalTracks = 0;

  // Get histogram of cos(2(phi - psi)); the average value is (average?) v_2
  for (Int_t i = 0; i < numOfEvents; i++)
    {
      totalTracks += v_events.at(i).nTracks;
      psi = v_events.at(i).psi;

      for (Int_t j = 0; j < v_events.at(i).nTracks; j++)
	{
	  phi = v_events.at(i).phiValues.at(j);
	  
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
  h_Xn          ->Write();
  h_Yn          ->Write();
  h_Xn_s        ->Write();
  h_Yn_s        ->Write();  
  h_psi         ->Write();
  h_psi_flat    ->Write();
  h_v2Plot      ->Write();

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
	  events.at(i).Xn -= -0.26214811775775954;                     // 2000 picoDsts: -0.26214811775775954
	  h_shifted->Fill(events.at(i).Xn);
	}
    }
  else if (component.EqualTo("y") || component.EqualTo("Y")) 
    {
      for (Int_t i = 0; i < numOfEvents; i++)
	{
	  events.at(i).Yn -= -0.39964040264945927;                     // 2000 picoDsts: -0.39964040264945927
	  h_shifted->Fill(events.at(i).Yn);
	}
    }

  return h_shifted;
}//End recenter()



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
