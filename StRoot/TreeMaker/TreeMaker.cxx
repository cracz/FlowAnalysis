// This is an version of Yang Wu's (UCR) tree maker adapted for my own (Cameron Racz) analysis.


// Load header files
#include "TreeMaker.h"
#include "TF1.h"

// Bichsel header
#include "StRoot/StBichsel/Bichsel.h"

// Config header
#include "../ConfigReader/ConfigReader.h"

// Bichsel Function Functions
Double_t bichselZ(Double_t *x,Double_t *par) 
{
  Double_t pove   = x[0];
  Double_t poverm = pove/par[0];
  return TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(poverm),par[1]));
}

Double_t bichsel70(Double_t *x,Double_t *par) 
{
  Double_t pove   = x[0];
  Double_t poverm = pove/par[0];
  return TMath::Exp(Bichsel::Instance()->GetI70(TMath::Log10(poverm),par[1]));
}

// Class implementation in CINT
ClassImp(TreeMaker)

// Constructor
TreeMaker::TreeMaker(StPicoDstMaker* Maker, std::string configFileName, TString JobId, Int_t EventsNumber, Double_t inputParameter) : StMaker()
{
  // Initialize parameters
  Y_MID   = -1.135;
  mPicoDstMaker = Maker;
  configs.read(configFileName);
  JobIdName = JobId;
  cutTest = inputParameter;
  JobIdName.Append(".root"); // Name output file by assigned Job ID
}

// Destructor
TreeMaker::~TreeMaker() {}

Int_t TreeMaker::Init()
{
  outputFile = new TFile(JobIdName,"recreate");
  
  // New Tree
  FxtTree = new TTree("Autree","TTree to hold FXT events and tracks");
  FxtTree->Branch("runId",&tree_runId,"runId/I");
  FxtTree->Branch("eventId",&tree_eventId,"eventId/I");
  FxtTree->Branch("bField",&tree_bField,"bField/F");
  FxtTree->Branch("Vx",&tree_Vx,"Vx/F");
  FxtTree->Branch("Vy",&tree_Vy,"Vy/F");
  FxtTree->Branch("Vz",&tree_Vz,"Vz/F");
  FxtTree->Branch("centrality",&tree_centrality,"centrality/s");
  FxtTree->Branch("tracknumber",&tree_tracknumber,"tracknumber/s");
  FxtTree->Branch("PID",tree_PID,"PID[tracknumber]/S");
  FxtTree->Branch("Charge",tree_Charge,"Charge[tracknumber]/S");
  FxtTree->Branch("Px",tree_Px,"Px[tracknumber]/F");
  FxtTree->Branch("Py",tree_Py,"Py[tracknumber]/F");
  FxtTree->Branch("Pz",tree_Pz,"Pz[tracknumber]/F");
  FxtTree->Branch("nEPDhits",&tree_nEPDhits,"nEPDhits/s");
  FxtTree->Branch("EPDid",tree_EPDid,"EPDid[nEPDhits]/S");
  FxtTree->Branch("EPDnMip",tree_EPDnMip,"EPDnMip[nEPDhits]/F");

  // QA plots
  h_eventCheck = new TH1D("h_eventCheck","Event number after each cut;;Events", 6, 0, 6);
  h_eventCheck->GetXaxis()->SetBinLabel(1,"No Cuts");
  h_eventCheck->GetXaxis()->SetBinLabel(2,"Bad Run Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(3,"Minbias Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(4,"V_{Z} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(5,"V_{r} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(6,"Centrality Cut");
    
  h_trackCheck = new TH1D("h_trackCheck","Track number after each cut;;Tracks", 3, 0, 3);
  h_trackCheck->GetXaxis()->SetBinLabel(1,"Event Cuts Only");
  h_trackCheck->GetXaxis()->SetBinLabel(2,"QA Cuts");
  h_trackCheck->GetXaxis()->SetBinLabel(3,"PID Cuts");
    
  h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", CENT_BINS, FIRST_CENT, FIRST_CENT + CENT_BINS);
  h_centralities->GetXaxis()->SetTitle("Centrality bin");
  h_centralities->GetYaxis()->SetTitle("# of events");
    
  h_zvtx =new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", 100, 190, 210);

  h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
    
  h_refmult = new TH1D("h_refmult","Reference multiplicity",1001,-0.5,1000.5);
  h_refmult->GetXaxis()->SetTitle("RefMult");
  h_refmult->GetXaxis()->SetTitle("# of events");
    
  h_tofmult = new TH1D("h_tofmult","TOF multiplicity",1001,-0.5,1000.5);
  h_tofmult->GetXaxis()->SetTitle("TofMult");
  h_tofmult->GetXaxis()->SetTitle("# of events");
    
  h_trackmult = new TH1D("h_trackmult","Actual track multiplicity",1501,-0.5,1500.5);
  h_trackmult->GetXaxis()->SetTitle("TrackMult");
  h_trackmult->GetXaxis()->SetTitle("# of events");
    
  h2_refmult_vs_trackmult = new TH2D("h2_refmult_vs_trackmult","RefMult vs. Actual track multiplicity;TrackMult;RefMult",1501,-0.5,1500.5,1001,-0.5,1000.5);
    
  h2_tofmult_vs_trackmult = new TH2D("h2_tofmult_vs_trackmult","TofMult vs. Actual track multiplicity;TrackMult;TofMult",1501,-0.5,1500.5,1001,-0.5,1000.5);
    
  h2_tofmult_vs_refmult = new TH2D("h2_tofmult_vs_refmult","TofMult vs. RefMult;RefMult;TofMult",1001,-0.5,1000.5,1001,-0.5,1000.5);

  h_nhits       = new TH1D("h_nhits", "nHits;Number of hits;Tracks", 50, 0, 50);
  h_nhits_dEdx  = new TH1D("h_nhits_dEdx","nHitsdEdx;Number of hits;Tracks", 50, 0, 50);
  h_nhits_ratio = new TH1D("h_nhits_ratio","nhitsFit/nhitsPoss;nhitsFit/nhitsPoss;Tracks",200,0.0,2.0);    
  h_DCA         = new TH1D("h_DCA","DCA (cm);DCA (cm);Tracks",100,0.0,10.0);
    
  h_pT = new TH1D("h_pT","p_{T};p_{T};Tracks",1000,0.0,10.0);
  h_eta = new TH1D("h_eta","#eta;#eta;Tracks",500,-6.0,4.0);
  h_phi = new TH1D("h_phi","#phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h2_dEdx_vs_qp = new TH2D("h2_dEdx_vs_qp", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_half = new TH2D("h2_dEdx_vs_qp_half", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, 0, 4, 500, 0, 12);
  h2_beta_vs_qp = new TH2D("h2_beta_vs_qp","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_m2_vs_qp = new TH2D("h2_m2_vs_qp", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);

 
  h_pT_pp = new TH1D("h_pT_pp","#pi^{+} p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_pm = new TH1D("h_pT_pm","#pi^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_kp = new TH1D("h_pT_kp","K^{+} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_km = new TH1D("h_pT_km","K^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_pr = new TH1D("h_pT_pr","Proton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_de = new TH1D("h_pT_de","Deuteron p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_tr = new TH1D("h_pT_tr","Triton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  
  h_eta_pp = new TH1D("h_eta_pp","#pi^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_pm = new TH1D("h_eta_pm","#pi^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_kp = new TH1D("h_eta_kp","K^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_km = new TH1D("h_eta_km","K^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_pr = new TH1D("h_eta_pr","Proton #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_de = new TH1D("h_eta_de","Deuteron #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_tr = new TH1D("h_eta_tr","Triton #eta;#eta;Tracks",500,-5.0,5.0);

  h_dndy_pp = new TH1D("h_dndy_pp", "#pi^{+} Raw Rapidity Spectrum;y;dN/dy", 40, -2, 0);
  h_dndy_pm = new TH1D("h_dndy_pm", "#pi^{-} Raw Rapidity Spectrum;y;dN/dy", 40, -2, 0);
  h_dndy_kp = new TH1D("h_dndy_kp", "K^{+} Raw Rapidity Spectrum;y;dN/dy",   40, -2, 0);
  h_dndy_km = new TH1D("h_dndy_km", "K^{-} Raw Rapidity Spectrum;y;dN/dy",   40, -2, 0);
  h_dndy_pr = new TH1D("h_dndy_pr", "Proton Raw Rapidity Spectrum;y;dN/dy",  40, -2, 0);
  h_dndy_de = new TH1D("h_dndy_de", "Deuteron Raw Rapidity Spectrum;y;dN/dy",  40, -2, 0);
  h_dndy_tr = new TH1D("h_dndy_tr", "Triton Raw Rapidity Spectrum;y;dN/dy",  40, -2, 0);
  
  h_phi_pp = new TH1D("h_phi_pp","#pi^{+} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_pm = new TH1D("h_phi_pm","#pi^{-} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_kp = new TH1D("h_phi_kp","K^{+} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_km = new TH1D("h_phi_km","K^{-} #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_pr = new TH1D("h_phi_pr","Proton #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_de = new TH1D("h_phi_de","Deuteron #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  h_phi_tr = new TH1D("h_phi_tr","Triton #phi (Radian);#phi;Tracks",1000,-1.5*PI,1.5*PI);
  
  h_mult_pp = new TH1D("h_mult_pp","#pi^{+} track multiplicity;#pi^{+} Mult;Events",1001,-0.5,1000.5);
  h_mult_pm = new TH1D("h_mult_pm","#pi^{-} track multiplicity;#pi^{-} Mult;Events",1001,-0.5,1000.5);
  h_mult_kp = new TH1D("h_mult_kp","K^{#plus} track multiplicity;K^{+} Mult;Events",1001,-0.5,1000.5);
  h_mult_km = new TH1D("h_mult_km","K^{-} track multiplicity;K^{-} Mult;Events",1001,-0.5,1000.5);
  h_mult_pr = new TH1D("h_mult_pr","Proton track multiplicity;Proton Mult;Events",1001,-0.5,1000.5);
  h_mult_de = new TH1D("h_mult_de","Deuteron track multiplicity;Deuteron Mult;Events",1001,-0.5,1000.5);
  h_mult_tr = new TH1D("h_mult_tr","Triton track multiplicity;Triton Mult;Events",1001,-0.5,1000.5);
  
  h2_dEdx_vs_qp_pp = new TH2D("h2_dEdx_vs_qp_pp", "#pi^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_pm = new TH2D("h2_dEdx_vs_qp_pm", "#pi^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_kp = new TH2D("h2_dEdx_vs_qp_kp", "K^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_km = new TH2D("h2_dEdx_vs_qp_km", "K^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_pr = new TH2D("h2_dEdx_vs_qp_pr", "Proton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_de = new TH2D("h2_dEdx_vs_qp_de", "Deuteron dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_tr = new TH2D("h2_dEdx_vs_qp_tr", "Triton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  
  h2_beta_vs_qp_pp = new TH2D("h2_beta_vs_qp_pp","#pi^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_pm = new TH2D("h2_beta_vs_qp_pm","#pi^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_kp = new TH2D("h2_beta_vs_qp_kp","K^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_km = new TH2D("h2_beta_vs_qp_km","K^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_pr = new TH2D("h2_beta_vs_qp_pr","Proton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_de = new TH2D("h2_beta_vs_qp_de","Deuteron 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_tr = new TH2D("h2_beta_vs_qp_tr","Triton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
 
  h2_m2_vs_qp_pp = new TH2D("h2_m2_vs_qp_pp", "#pi^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_pm = new TH2D("h2_m2_vs_qp_pm", "#pi^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_kp = new TH2D("h2_m2_vs_qp_kp", "K^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_km = new TH2D("h2_m2_vs_qp_km", "K^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_pr = new TH2D("h2_m2_vs_qp_pr", "Proton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_de = new TH2D("h2_m2_vs_qp_de", "Deuteron m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_tr = new TH2D("h2_m2_vs_qp_tr", "Triton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  
  h2_pT_vs_yCM_pp = new TH2D("h2_pT_vs_yCM_pp", "#pi^{+};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_pm = new TH2D("h2_pT_vs_yCM_pm", "#pi^{-};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_kp = new TH2D("h2_pT_vs_yCM_kp", "K^{+};y-y_{mid};p_{T} (GeV/c)",   300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_km = new TH2D("h2_pT_vs_yCM_km", "K^{-};y-y_{mid};p_{T} (GeV/c)",   300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_pr = new TH2D("h2_pT_vs_yCM_pr", "Proton;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_de = new TH2D("h2_pT_vs_yCM_de", "Deuteron;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_tr = new TH2D("h2_pT_vs_yCM_tr", "Triton;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);

  return kStOK;
}

Int_t TreeMaker::Finish()
{
  outputFile->Write();
  return kStOK;
}

Int_t TreeMaker::Make()
{
  //=========================================================
  //          Bichsel Function Setup
  //=========================================================
  Double_t log2dx = 1.0;
  Double_t xStart = 0.01;
  Double_t xStop  = 3.0;
  Int_t npx = 10000;
  //                      Mass  log2(dx)
  Double_t params[2] = {  1.0,   log2dx  };

  params[0] = D_M0_DE;
  TF1 *bichselZ_de = new TF1(Form("BichselZ_de_log2dx_%i",(int)log2dx),bichselZ,xStart,xStop,2);
  if (!bichselZ_de) { std::cout << "De function error" << std::endl; return 1; }
  bichselZ_de->SetParameters(params); 
  bichselZ_de->SetNpx(npx);

  params[0] = D_M0_TR;
  TF1 *bichselZ_tr = new TF1(Form("BichselZ_tr_log2dx_%i",(int)log2dx),bichselZ,xStart,xStop,2);
  if (!bichselZ_tr) { std::cout << "Tr function error" << std::endl; return 1; }
  bichselZ_tr->SetParameters(params); 
  bichselZ_tr->SetNpx(npx);
  //=========================================================
  //          END Bichsel Function Setup
  //=========================================================


  h_eventCheck->Fill(0); // Count # of events before any cuts

  StPicoEvent* event = mPicoDstMaker->picoDst()->event(); // Get Event pointer
    
  if(event) // Ensure event pointer is not NULL
    { 
      if( IsGoodRun(event->runId()) )  // Select only good runs
	{ 
	  h_eventCheck->Fill(1); // Count # of good run events

	  /*
	  Bool_t b_good_trig = false;
	  std::vector<UInt_t> triggerIDs = event->triggerIds();
	  for(UInt_t i = 0; i < triggerIDs.size(); i++)
	    { if(triggerIDs.at(i) == 620052) b_good_trig = true; }
	  */
	  Bool_t b_good_trig = true;


	  if(b_good_trig) // Minbias cut
	    { 
	      h_eventCheck->Fill(2); // Count # of events after minbias trigger id cut
	      
	      TVector3 pVtx = event->primaryVertex();
	      Double_t d_xvtx = pVtx.x();
	      Double_t d_yvtx = pVtx.y();
	      Double_t d_zvtx = pVtx.z();
	      Double_t d_rvtx = TMath::Sqrt(d_xvtx * d_xvtx + (d_yvtx + 2) * (d_yvtx + 2));
	      
	      h_zvtx->Fill(pVtx.Z());
	      h2_trans_vtx->Fill(pVtx.X(),pVtx.Y());
	      
	      if(d_zvtx > configs.z_vtx_low && d_zvtx < configs.z_vtx_high) // Vz cut
		{ 
		  h_eventCheck->Fill(3);    // Count # of events after Vz cut
		  
		  if(d_rvtx < configs.r_vtx) // Vr cut
		    { 
		      h_eventCheck->Fill(4); // Count # of events after Vr cut
		      
		      Int_t nTracks = mPicoDstMaker->picoDst()->numberOfTracks(); // Get the number of Primary Tracks
		      Int_t primTracks = 0;
		      Int_t N_pp = 0;
		      Int_t N_pm = 0;
		      Int_t N_kp = 0;
		      Int_t N_km = 0;
		      Int_t N_pr = 0;
		      Int_t N_de = 0;
		      Int_t N_tr = 0;

		      // 1st Primary Tracks loop to get centrality starts
		      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
			{
			  StPicoTrack* track = mPicoDstMaker->picoDst()->track(iTrk); // Get Track pointer
			  if(!track->isPrimary()) continue; // Only Primary Tracks
			  primTracks++;
			  h_nhits->Fill(track->nHits());
			  h_nhits_dEdx->Fill(track->nHitsDedx());
			  h_nhits_ratio->Fill((Double_t)track->nHitsFit()/track->nHitsPoss());
			  h_DCA->Fill(track->gDCA(pVtx).Mag());
			} // 1st Primary tracks loop ends

		      h_refmult->Fill(event->refMultHalfEast());
		      h_trackmult->Fill(primTracks);
		      h_tofmult->Fill(event->nBTOFMatch());
		      h2_refmult_vs_trackmult->Fill(primTracks,event->refMultHalfEast());
		      h2_tofmult_vs_trackmult->Fill(primTracks,event->nBTOFMatch());
		      h2_tofmult_vs_refmult->Fill(event->refMultHalfEast(),event->nBTOFMatch());


		      // GET CENTRALITY
		      Int_t centrality = -99;

		      // 3.0 GeV FXT  --  From Zachary Sweger Nov 11, 2020
		      if (configs.sqrt_s_NN == 3.0)
			{
			  if     ( primTracks >=   5 && primTracks <=   6 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   7 && primTracks <=   8 ) centrality =  1;
			  else if( primTracks >=   9 && primTracks <=  11 ) centrality =  2;
			  else if( primTracks >=  12 && primTracks <=  15 ) centrality =  3;
			  else if( primTracks >=  16 && primTracks <=  20 ) centrality =  4;
			  else if( primTracks >=  21 && primTracks <=  25 ) centrality =  5;
			  else if( primTracks >=  26 && primTracks <=  32 ) centrality =  6;
			  else if( primTracks >=  33 && primTracks <=  40 ) centrality =  7;
			  else if( primTracks >=  41 && primTracks <=  49 ) centrality =  8;
			  else if( primTracks >=  50 && primTracks <=  59 ) centrality =  9;
			  else if( primTracks >=  60 && primTracks <=  71 ) centrality = 10;
			  else if( primTracks >=  72 && primTracks <=  85 ) centrality = 11;
			  else if( primTracks >=  86 && primTracks <= 100 ) centrality = 12;
			  else if( primTracks >= 101 && primTracks <= 118 ) centrality = 13;
			  else if( primTracks >= 119 && primTracks <= 141 ) centrality = 14;
			  else if( primTracks >= 142 && primTracks <= 195 ) centrality = 15;  // 0% - 5% (Central)
			}

		      // 7.2 GeV FXT
		      else if (configs.sqrt_s_NN == 7.2)
			{
			  if     ( primTracks >=   2 && primTracks <=   3 ) centrality =  0;  // 75% - 80% (Peripheral)
			  else if( primTracks >=   4 && primTracks <=   5 ) centrality =  1;
			  else if( primTracks >=   6 && primTracks <=   8 ) centrality =  2;
			  else if( primTracks >=   9 && primTracks <=  11 ) centrality =  3;
			  else if( primTracks >=  12 && primTracks <=  15 ) centrality =  4;
			  else if( primTracks >=  16 && primTracks <=  21 ) centrality =  5;
			  else if( primTracks >=  22 && primTracks <=  29 ) centrality =  6;
			  else if( primTracks >=  30 && primTracks <=  38 ) centrality =  7;
			  else if( primTracks >=  39 && primTracks <=  49 ) centrality =  8;
			  else if( primTracks >=  50 && primTracks <=  63 ) centrality =  9;
			  else if( primTracks >=  64 && primTracks <=  79 ) centrality = 10;
			  else if( primTracks >=  80 && primTracks <=  99 ) centrality = 11;
			  else if( primTracks >= 100 && primTracks <= 123 ) centrality = 12;
			  else if( primTracks >= 124 && primTracks <= 153 ) centrality = 13;
			  else if( primTracks >= 154 && primTracks <= 190 ) centrality = 14;
			  else if( primTracks >= 191 && primTracks <= 240 ) centrality = 15;  // 0% - 5% (Central)
			}


		      // Select only good centrality valued events
		      if(centrality >= FIRST_CENT)
			{
			  h_eventCheck->Fill(5); // Count # of min-bias events after all cuts
			  h_centralities->Fill(centrality);

			  // Prepare New Tree parameters
			  tree_Vx = d_xvtx;
			  tree_Vy = d_yvtx;
			  tree_Vz = d_zvtx;
			  tree_runId   = event->runId();
			  tree_eventId = event->eventId();
			  tree_bField  = event->bField();
			  tree_centrality = centrality;
			  /*tree_tracknumber = primTracks;*/
			  Int_t realTrackIndex = 0;


			  tree_nEPDhits = mPicoDstMaker->picoDst()->numberOfEpdHits();
			  if(tree_nEPDhits)
			    {
			      for(Int_t hit=0; hit<tree_nEPDhits; hit++)
				{
				  StPicoEpdHit *pEpdHit = mPicoDstMaker->picoDst()->epdHit(hit);
				  if(!pEpdHit) continue;
				  if(!pEpdHit->isGood()) continue;
				  tree_EPDid[hit]   = pEpdHit->id();
				  tree_EPDnMip[hit] = pEpdHit->nMIP();
				}
			    }
			    
			  // 2nd Loop through tracks to fill histograms
			  for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
			    {
			      StPicoTrack* track = mPicoDstMaker->picoDst()->track(iTrk); // Get Track pointer
			      if(!track->isPrimary()) continue;

			      h_trackCheck->Fill(0);  // Only event cuts
			      
			      //=========================================================
			      //          Track QA Cuts
			      //=========================================================
			      bool b_bad_hits  = ( track->nHits() < configs.nHits );
			      bool b_bad_dEdx  = ( track->nHitsDedx() <= configs.nHits_dEdx );
			      bool b_bad_ratio = ( ((double)track->nHitsFit() / (double)track->nHitsPoss()) <= configs.nHits_ratio );
			      bool b_bad_DCA   = ( track->gDCA(d_xvtx,d_yvtx,d_zvtx) >= configs.dca );

			      if (b_bad_hits || b_bad_dEdx || b_bad_ratio || b_bad_DCA) continue;
			      //=========================================================
			      //          End Track QA Cuts
			      //=========================================================
			      
			      h_trackCheck->Fill(1); // After QA Cuts
			      
			      Double_t d_pT = track->pPt();
			      Double_t d_px = track->pMom().X();
			      Double_t d_py = track->pMom().Y();
			      Double_t d_pz = track->pMom().Z();
			      Double_t d_eta = track->pMom().Eta();
			      Double_t d_phi = track->pMom().Phi();
			      Double_t d_mom = track->pMom().Mag();
			      Double_t d_dEdx = track->dEdx();
			      Short_t  s_charge = track->charge();
			      Double_t d_TPCnSigmaPion   = track->nSigmaPion();
			      Double_t d_TPCnSigmaKaon   = track->nSigmaKaon();
			      Double_t d_TPCnSigmaProton = track->nSigmaProton();
			      Double_t d_zDeuteron = (s_charge > 0) ? TMath::Log(d_dEdx / bichselZ_de->Eval(d_mom)) : -999.0;
			      Double_t d_zTriton   = (s_charge > 0) ? TMath::Log(d_dEdx / bichselZ_tr->Eval(d_mom)) : -999.0;

			      h_phi->Fill(d_phi);
			      h_eta->Fill(d_eta);
			      h_pT->Fill(d_pT);
			      h2_dEdx_vs_qp->Fill(s_charge * d_mom, track->dEdx());
			      
			      Bool_t tofTrack = false;
			      Double_t d_tofBeta = -999.0;
			      Double_t d_m2 = 0.0;
			      Int_t trackTofIndex = track->bTofPidTraitsIndex();
			      if(trackTofIndex >= 0)
				d_tofBeta = mPicoDstMaker->picoDst()->btofPidTraits(trackTofIndex)->btofBeta();
			      
			      if(d_tofBeta != -999.0)
				{
				  tofTrack = true;
				  d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
				  h2_beta_vs_qp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp->Fill(s_charge*d_mom, d_m2);
				}

			      //=========================================================
			      //          PID Cuts
			      //=========================================================
			      Bool_t pion   = false;
			      Bool_t kaon   = false;
			      Bool_t proton = (d_TPCnSigmaProton > configs.nSig_pr_low) && (d_TPCnSigmaProton < configs.nSig_pr_high) && (s_charge > 0);
			      Bool_t deuteron = false;
			      Bool_t triton   = false;
			      //Bool_t deuteron = (d_zDeuteron > configs.z_de_low) && (d_zDeuteron < configs.z_de_high);
			      //Bool_t triton   = (d_zTriton > configs.z_tr_low) && (d_zTriton < configs.z_tr_high);
			      // d_zDeuteron and d_zDeuteron already ensure that s_charge > 0

			      if (tofTrack)
				{
				  pion = (d_TPCnSigmaPion > configs.nSig_pi_low) &&
				         (d_TPCnSigmaPion < configs.nSig_pi_high) &&
				         (d_m2 > configs.m2_pi_low) &&
				         (d_m2 < configs.m2_pi_high);

				  kaon = (d_TPCnSigmaKaon > configs.nSig_ka_low) &&
				         (d_TPCnSigmaKaon < configs.nSig_ka_high) &&
				         (d_m2 > configs.m2_ka_low) &&
				         (d_m2 < configs.m2_ka_high);

				  deuteron = (d_zDeuteron > configs.z_de_low) &&
				             (d_zDeuteron < configs.z_de_high) &&
				             (d_m2 > configs.m2_de_low) &&
				             (d_m2 < configs.m2_de_high);

				  triton   = (d_zTriton > configs.z_tr_low) &&
				             (d_zTriton < configs.z_tr_high) &&
				             (d_m2 > configs.m2_tr_low) &&
				             (d_m2 < configs.m2_tr_high);
				}

			      //if (deuteron) h2_dEdx_vs_qp_half_postZdCut->Fill(s_charge * d_mom, d_dEdx);
			      //if (triton)   h2_dEdx_vs_qp_half_postZtCut->Fill(s_charge * d_mom, d_dEdx);
	    
			      //if (!pion && !kaon && !proton) continue;

			      if (pion     && proton) { proton = false; }
			      if (kaon     && proton) { proton = false; }
			      if (deuteron && proton) { proton = false; }
			      if (triton   && proton) { proton = false; }

			      //if (pion && deuteron) { deuteron = false; }
			      //if (pion && triton)   { triton   = false; }
			      //if (kaon && deuteron) { deuteron = false; }
			      //if (kaon && triton)   { triton   = false; }

			      /*
				if (deuteron && proton) 
				{ 
				if (TMath::Abs(d_zDeuteron) < TMath::Abs(d_TPCnSigmaProton)) { proton = false; }
				else if (TMath::Abs(d_zDeuteron) == TMath::Abs(d_TPCnSigmaProton)) { proton = false; deuteron = false; }
				else { deuteron = false; }
				}
				if (triton && proton)
				{
				if (TMath::Abs(d_zTriton) < TMath::Abs(d_TPCnSigmaProton)) { proton = false; }
				else if (TMath::Abs(d_zTriton) == TMath::Abs(d_TPCnSigmaProton)) { proton = false; triton = false; }
				else { triton = false; }
				}
				if (deuteron && triton)
				{
				if (TMath::Abs(d_zDeuteron) < TMath::Abs(d_zTriton)) { triton = false; }
				else if (TMath::Abs(d_zDeuteron) == TMath::Abs(d_zTriton)) { triton = false; deuteron = false; }
				else { deuteron = false; }
				}
			      */
			      //=========================================================
			      //          END PID Cuts
			      //=========================================================

			      if (pion || kaon || proton || deuteron || triton) h_trackCheck->Fill(2);

			      Double_t mRapidity = 999.0;

			      if(pion) // PID Pions
				{ 
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_PI);
				  
				  if(s_charge > 0)
				    {
				      N_pp++;
				      h_eta_pp->Fill(d_eta);
				      h_phi_pp->Fill(d_phi);
				      h_pT_pp->Fill(d_pT);
				      h_dndy_pp->Fill(mRapidity);
				      h2_pT_vs_yCM_pp->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_pp->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_pp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_pp->Fill(s_charge*d_mom, d_m2);
				    }
				  else if(s_charge < 0)
				    {
				      N_pm++;
				      h_eta_pm->Fill(d_eta);
				      h_phi_pm->Fill(d_phi);
				      h_pT_pm->Fill(d_pT);
				      h_dndy_pm->Fill(mRapidity);
				      h2_pT_vs_yCM_pm->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_pm->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_pm->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_pm->Fill(s_charge*d_mom, d_m2);
				    }
				}
			      else if(kaon) // PID Kaons
				{
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_KA);
				  
				  if(s_charge > 0)
				    {
				      N_kp++;
				      h_eta_kp->Fill(d_eta);
				      h_phi_kp->Fill(d_phi);
				      h_pT_kp->Fill(d_pT);
				      h_dndy_kp->Fill(mRapidity);
				      h2_pT_vs_yCM_kp->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_kp->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_kp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_kp->Fill(s_charge*d_mom, d_m2);
				    }
				  if(s_charge < 0)
				    {
				      N_km++;
				      h_eta_km->Fill(d_eta);
				      h_phi_km->Fill(d_phi);
				      h_pT_km->Fill(d_pT);
				      h_dndy_km->Fill(mRapidity);
				      h2_pT_vs_yCM_km->Fill(mRapidity - Y_MID, d_pT);
				      h2_dEdx_vs_qp_km->Fill(s_charge*d_mom, track->dEdx());
				      h2_beta_vs_qp_km->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				      h2_m2_vs_qp_km->Fill(s_charge*d_mom, d_m2);
				    }
				}
			      else if(proton) // PID Proton
				{
				  N_pr++;
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_PR);
				  
				  h_eta_pr->Fill(d_eta);
				  h_phi_pr->Fill(d_phi);
				  h_pT_pr->Fill(d_pT);
				  h_dndy_pr->Fill(mRapidity);
				  h2_pT_vs_yCM_pr->Fill(mRapidity - Y_MID, d_pT);
				  h2_dEdx_vs_qp_pr->Fill(s_charge*d_mom, track->dEdx());
				  h2_beta_vs_qp_pr->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp_pr->Fill(s_charge*d_mom, d_m2);
				}		    
			      else if(deuteron) // PID Deuteron
				{
				  N_de++;
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_DE);
				  
				  h_eta_de->Fill(d_eta);
				  h_phi_de->Fill(d_phi);
				  h_pT_de->Fill(d_pT);
				  h_dndy_de->Fill(mRapidity);
				  h2_pT_vs_yCM_de->Fill(mRapidity - Y_MID, d_pT);
				  h2_dEdx_vs_qp_de->Fill(s_charge*d_mom, track->dEdx());
				  h2_beta_vs_qp_de->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp_de->Fill(s_charge*d_mom, d_m2);
				}		    
			      else if(triton) // PID Triton
				{
				  N_tr++;
				  mRapidity = getRapidity(d_px, d_py, d_pz, D_M0_TR);
				  
				  h_eta_tr->Fill(d_eta);
				  h_phi_tr->Fill(d_phi);
				  h_pT_tr->Fill(d_pT);
				  h_dndy_tr->Fill(mRapidity);
				  h2_pT_vs_yCM_tr->Fill(mRapidity - Y_MID, d_pT);
				  h2_dEdx_vs_qp_tr->Fill(s_charge*d_mom, track->dEdx());
				  h2_beta_vs_qp_tr->Fill(s_charge*d_mom, 1.0/d_tofBeta);
				  h2_m2_vs_qp_tr->Fill(s_charge*d_mom, d_m2);
				}		    

			      tree_Charge[realTrackIndex] = s_charge;
			      if(pion)          tree_PID[realTrackIndex] = 0;
			      else if(kaon)     tree_PID[realTrackIndex] = 1;
			      else if(proton)   tree_PID[realTrackIndex] = 2;
			      else if(deuteron) tree_PID[realTrackIndex] = 3;
			      else if(triton)   tree_PID[realTrackIndex] = 4;
			      else              tree_PID[realTrackIndex] = -1;

			      tree_Px[realTrackIndex] = d_px;
			      tree_Py[realTrackIndex] = d_py;
			      tree_Pz[realTrackIndex] = d_pz;
			      realTrackIndex += 1;
			    } // 2nd track loop ends
			    
			  h_mult_pp->Fill(N_pp);
			  h_mult_pm->Fill(N_pm);
			  h_mult_kp->Fill(N_kp);
			  h_mult_km->Fill(N_km);
			  h_mult_pr->Fill(N_pr);
			  h_mult_de->Fill(N_de);
			  h_mult_tr->Fill(N_tr);
			  //tree_tracknumber = N_pr + N_pp + N_pm + N_kp + N_km + N_de + N_tr;
			  tree_tracknumber = realTrackIndex + 1;
			  FxtTree->Fill();			    
			} // Good centrality event selection ends
                    } // Vr cut ends
                } // Vz cut ends
            } // Minbias cut ends
        } // Good run cut ends
    } // Event pointer not NULL cut ends
  return kStOK;
}

Bool_t TreeMaker::IsGoodRun(Int_t runnumber)
{
  // From Ben Kimelman Nov 6, 2020
  Int_t badRunList_3p0GeV[24] = {19151029, 19151045, 19152001, 19152078, 19153023, 19153032, 19153065, 19154012, 19154013, 19154014, 19154015, 19154016, 
    19154017, 19154018, 19154019, 19154020, 19154021, 19154022, 19154023, 19154024, 19154026, 19154046, 19154051, 19154056};

  Bool_t b_good_run = true;
  /*
    if (configs.sqrt_s_NN == 3.0)
    { for (Int_t i = 0; i < 24; i++) { if (event->runId() == badRunList_3p0GeV[i]) {b_bad_run = true; break;} } }
  */
  for (Int_t i = 0; i < 24; i++) { if (runnumber == badRunList_3p0GeV[i]) {b_good_run = false; break;} } 
  
  return b_good_run;
}


Double_t TreeMaker::getRapidity(Double_t px, Double_t py, Double_t pz, Double_t mass)
{
  Double_t rapidity, energy, momentum = 0;
  momentum = TMath::Sqrt(px*px + py*py + pz*pz);
  energy   = TMath::Sqrt(momentum*momentum + mass*mass);
  rapidity = TMath::ATanH(pz/energy);
  return rapidity;
}
