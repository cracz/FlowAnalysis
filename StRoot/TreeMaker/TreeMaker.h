// Preprocessor to avoid loading header file multiple times
#ifndef TreeMaker_def
#define TreeMaker_def

// Load C/C++ header files
#include <stdio.h>
#include <iostream>
#include <vector>

// Load ROOT header files
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

// Load STARLibrary header files
#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"

// Configuration file reader
#include "../ConfigReader/ConfigReader.h"

class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

const Double_t PI = TMath::Pi();
const Int_t Ncentralities =  16;
const Int_t N_track       = 195;
const Int_t CENT_BINS  = 16;
const Int_t FIRST_CENT = 16 - CENT_BINS;
const Double_t D_M0_PI = 0.139571;   //Rest masses
const Double_t D_M0_KA = 0.493677;
const Double_t D_M0_PR = 0.938272;
const Double_t D_M0_DE = 1.875613;   // Deuteron
const Double_t D_M0_TR = 2.808921;   // Triton


class TreeMaker : public StMaker
{
private:
  StPicoDstMaker* mPicoDstMaker;
  ConfigReader    configs;
  TString         JobIdName;
  Double_t        cutTest;
  TFile*          outputFile;
  Double_t        Y_MID;
    
  TTree*              FxtTree;
  Int_t               tree_runId;
  Int_t               tree_eventId;
  Float_t             tree_bField;
  Float_t             tree_Vx;
  Float_t             tree_Vy;
  Float_t             tree_Vz;
  UShort_t            tree_centrality;
  UShort_t            tree_tracknumber;
  Short_t             tree_PID[N_track];
  Short_t             tree_Charge[N_track];
  Float_t             tree_Px[N_track];
  Float_t             tree_Py[N_track];
  Float_t             tree_Pz[N_track];
  UShort_t            tree_nEPDhits;
  Short_t             tree_EPDid[744];
  Float_t             tree_EPDnMip[744];

    
  TH1D*                h_eventCheck;
  TH1D*                h_trackCheck;
  TH1D*                h_centralities;
  TH1D*                h_zvtx;
  TH1D*                h_trackmult;
  TH1D*                h_refmult;
  TH1D*                h_tofmult;
  TH2D*                h2_trans_vtx;
  TH2D*                h2_trans_vtx_cut;
  TH2D*                h2_refmult_vs_trackmult;
  TH2D*                h2_tofmult_vs_trackmult;
  TH2D*                h2_tofmult_vs_refmult;
    
  TH1D*                h_nhits_dEdx;
  TH1D*                h_nhits;
  TH1D*                h_nhits_ratio;
  TH1D*                h_DCA;
  TH1D*                h_pT;
  TH1D*                h_eta;
  TH1D*                h_phi;
  TH2D*                h2_dEdx_vs_qp;
  TH2D*                h2_dEdx_vs_qp_half;
  TH2D*                h2_beta_vs_qp;
  TH2D*                h2_m2_vs_qp;

  TH1D*                h_mult_pp;
  TH1D*                h_mult_pm;
  TH1D*                h_mult_kp;
  TH1D*                h_mult_km;  
  TH1D*                h_mult_pr;
  TH1D*                h_mult_de;
  TH1D*                h_mult_tr;

  TH1D*                h_pT_pp;
  TH1D*                h_pT_pm;
  TH1D*                h_pT_kp;
  TH1D*                h_pT_km;
  TH1D*                h_pT_pr;
  TH1D*                h_pT_de;
  TH1D*                h_pT_tr;

  TH1D*                h_dndy_pp;
  TH1D*                h_dndy_pm;
  TH1D*                h_dndy_kp;
  TH1D*                h_dndy_km;
  TH1D*                h_dndy_pr;
  TH1D*                h_dndy_de;
  TH1D*                h_dndy_tr;

  TH1D*                h_eta_pp;
  TH1D*                h_eta_pm;
  TH1D*                h_eta_kp;
  TH1D*                h_eta_km;
  TH1D*                h_eta_pr;
  TH1D*                h_eta_de;
  TH1D*                h_eta_tr;

  TH1D*                h_phi_pp;
  TH1D*                h_phi_pm;
  TH1D*                h_phi_kp;
  TH1D*                h_phi_km;
  TH1D*                h_phi_pr;
  TH1D*                h_phi_de;
  TH1D*                h_phi_tr;
  
  TH2D*                h2_pT_vs_yCM_pp;
  TH2D*                h2_pT_vs_yCM_pm;
  TH2D*                h2_pT_vs_yCM_kp;
  TH2D*                h2_pT_vs_yCM_km;
  TH2D*                h2_pT_vs_yCM_pr;
  TH2D*                h2_pT_vs_yCM_de;
  TH2D*                h2_pT_vs_yCM_tr;
  
  TH2D*                h2_dEdx_vs_qp_pp;
  TH2D*                h2_dEdx_vs_qp_pm;
  TH2D*                h2_dEdx_vs_qp_kp;
  TH2D*                h2_dEdx_vs_qp_km;
  TH2D*                h2_dEdx_vs_qp_pr;
  TH2D*                h2_dEdx_vs_qp_de;
  TH2D*                h2_dEdx_vs_qp_tr;
  
  TH2D*                h2_beta_vs_qp_pm;
  TH2D*                h2_beta_vs_qp_kp;
  TH2D*                h2_beta_vs_qp_km;
  TH2D*                h2_beta_vs_qp_pr;
  TH2D*                h2_beta_vs_qp_pp;
  TH2D*                h2_beta_vs_qp_de;
  TH2D*                h2_beta_vs_qp_tr;

  TH2D*                h2_m2_vs_qp_pp;
  TH2D*                h2_m2_vs_qp_pm;
  TH2D*                h2_m2_vs_qp_kp;
  TH2D*                h2_m2_vs_qp_km;
  TH2D*                h2_m2_vs_qp_pr;
  TH2D*                h2_m2_vs_qp_de;
  TH2D*                h2_m2_vs_qp_tr;
  
public:
  TreeMaker(StPicoDstMaker* Maker, std::string configFileName, TString JobId, Int_t EventsNumber, Double_t inputParameter);
  virtual ~TreeMaker();
  Int_t    Init();
  Int_t    Make();
  Int_t    Finish();
  Bool_t   IsGoodRun(Int_t runNumber, Double_t sqrt_s_NN);
  Double_t getRapidity(Double_t px, Double_t py, Double_t pz, Double_t mass);
  ClassDef(TreeMaker,1) // Class title
};

// End of preprocessor
#endif
