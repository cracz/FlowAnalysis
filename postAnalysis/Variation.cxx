#include <iostream>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "Variation.h"
#include "PlotUtils.h"

using namespace PlotUtils;

// Constructor for the "normal" results.
Variation::Variation(TString prefix, TString order_n_str)
{
  ID = prefix;
  fileName = prefix + ".picoDst.result.combined.root";

  file = TFile::Open(fileName);

  initialize(order_n_str);
  fixAttributes(order_n_str);
}

// Constructor for the variations to be combined with the normal results.
Variation::Variation(TString prefix, TString order_n_str, Variation *normalData)
{
  ID = prefix;
  fileName = prefix + ".picoDst.result.combined.root";

  file = TFile::Open(fileName);

  initialize(order_n_str);
  combine(normalData);
  stdDevs();
}

// Destructor
Variation::~Variation()
{
  delete h_vn_pp;
  delete h_vn_pm;
  delete h_vn_kp;
  delete h_vn_km;
  delete h_vn_pr;
  delete h_vn_pp_ext;
  delete h_vn_pm_ext;
  delete h_vn_kp_ext;
  delete h_vn_km_ext;
  delete h_vn_pr_ext;
  delete h_vn_pr_for;
  delete h_vn_yCM_00to10_pp;
  delete h_vn_yCM_10to40_pp;
  delete h_vn_yCM_40to60_pp;
  delete h_vn_yCM_00to10_pm;
  delete h_vn_yCM_10to40_pm;
  delete h_vn_yCM_40to60_pm;
  delete h_vn_yCM_00to10_kp;
  delete h_vn_yCM_10to40_kp;
  delete h_vn_yCM_40to60_kp;
  delete h_vn_yCM_00to10_km;
  delete h_vn_yCM_10to40_km;
  delete h_vn_yCM_40to60_km;
  delete h_vn_yCM_00to10_pr;
  delete h_vn_yCM_10to40_pr;
  delete h_vn_yCM_40to60_pr;
  delete h_vn_yCM_00to10_pr_symm;
  delete h_vn_yCM_10to40_pr_symm;
  delete h_vn_yCM_40to60_pr_symm;
  delete h_vn_pT_00to10_pp;
  delete h_vn_pT_10to40_pp;
  delete h_vn_pT_40to60_pp;
  delete h_vn_pT_00to10_pm;
  delete h_vn_pT_10to40_pm;
  delete h_vn_pT_40to60_pm;
  delete h_vn_pT_00to10_kp;
  delete h_vn_pT_10to40_kp;
  delete h_vn_pT_40to60_kp;
  delete h_vn_pT_00to10_km;
  delete h_vn_pT_10to40_km;
  delete h_vn_pT_40to60_km;
  delete h_vn_pT_00to10_pr;
  delete h_vn_pT_10to40_pr;
  delete h_vn_pT_40to60_pr;
  file->Close();
}

// Initialize all histograms
void Variation::initialize(TString order_n_str)
{
  TProfile *p_vn_pp = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr = (TProfile*)file->Get("p_vn_pr");
  p_vn_kp->Rebin();
  p_vn_km->Rebin();
  
  h_vn_pp = p_vn_pp->ProjectionX((TString)p_vn_pp->GetName() +"_"+ ID);
  h_vn_pm = p_vn_pm->ProjectionX((TString)p_vn_pm->GetName() +"_"+ ID);
  h_vn_kp = p_vn_kp->ProjectionX((TString)p_vn_kp->GetName() +"_"+ ID);
  h_vn_km = p_vn_km->ProjectionX((TString)p_vn_km->GetName() +"_"+ ID);
  h_vn_pr = p_vn_pr->ProjectionX((TString)p_vn_pr->GetName() +"_"+ ID); 
  h_vn_pp = flipHisto(h_vn_pp);
  h_vn_pm = flipHisto(h_vn_pm);
  h_vn_kp = flipHisto(h_vn_kp);
  h_vn_km = flipHisto(h_vn_km);
  h_vn_pr = flipHisto(h_vn_pr);
  h_vn_pp = trimCentralityPlot(h_vn_pp);
  h_vn_pm = trimCentralityPlot(h_vn_pm);
  h_vn_kp = trimCentralityPlot(h_vn_kp);
  h_vn_km = trimCentralityPlot(h_vn_km);
  h_vn_pr = trimCentralityPlot(h_vn_pr);
  
  TProfile *p_vn_pp_ext = (TProfile*)file->Get("p_vn_pp_ext");
  TProfile *p_vn_pm_ext = (TProfile*)file->Get("p_vn_pm_ext");
  TProfile *p_vn_kp_ext = (TProfile*)file->Get("p_vn_kp_ext");
  TProfile *p_vn_km_ext = (TProfile*)file->Get("p_vn_km_ext");
  TProfile *p_vn_pr_ext = (TProfile*)file->Get("p_vn_pr_ext");
  p_vn_kp_ext->Rebin();
  p_vn_km_ext->Rebin();
  
  h_vn_pp_ext = p_vn_pp_ext->ProjectionX((TString)p_vn_pp_ext->GetName() +"_"+ ID);
  h_vn_pm_ext = p_vn_pm_ext->ProjectionX((TString)p_vn_pm_ext->GetName() +"_"+ ID);
  h_vn_kp_ext = p_vn_kp_ext->ProjectionX((TString)p_vn_kp_ext->GetName() +"_"+ ID);
  h_vn_km_ext = p_vn_km_ext->ProjectionX((TString)p_vn_km_ext->GetName() +"_"+ ID);
  h_vn_pr_ext = p_vn_pr_ext->ProjectionX((TString)p_vn_pr_ext->GetName() +"_"+ ID);
  h_vn_pp_ext = flipHisto(h_vn_pp_ext);
  h_vn_pm_ext = flipHisto(h_vn_pm_ext);
  h_vn_kp_ext = flipHisto(h_vn_kp_ext);
  h_vn_km_ext = flipHisto(h_vn_km_ext);
  h_vn_pr_ext = flipHisto(h_vn_pr_ext);
  h_vn_pp_ext = trimCentralityPlot(h_vn_pp_ext);
  h_vn_pm_ext = trimCentralityPlot(h_vn_pm_ext);
  h_vn_kp_ext = trimCentralityPlot(h_vn_kp_ext);
  h_vn_km_ext = trimCentralityPlot(h_vn_km_ext);
  h_vn_pr_ext = trimCentralityPlot(h_vn_pr_ext);

  TProfile *p_vn_pr_for = (TProfile*)file->Get("p_vn_pr_for");
  h_vn_pr_for = p_vn_pr_for->ProjectionX((TString)p_vn_pr_for->GetName() +"_"+ ID);
  h_vn_pr_for = flipHisto(h_vn_pr_for);
  h_vn_pr_for = trimCentralityPlot(h_vn_pr_for);

  
  //=== vn vs rapidity stuff
  TProfile2D *p2_vn_yCM_cent_pp = (TProfile2D*)file->Get("p2_vn_yCM_cent_pp");
  TProfile2D *p2_vn_yCM_cent_pm = (TProfile2D*)file->Get("p2_vn_yCM_cent_pm");
  TProfile2D *p2_vn_yCM_cent_kp = (TProfile2D*)file->Get("p2_vn_yCM_cent_kp");
  TProfile2D *p2_vn_yCM_cent_km = (TProfile2D*)file->Get("p2_vn_yCM_cent_km");
  TProfile2D *p2_vn_yCM_cent_pr = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr");
  TProfile2D *p2_vn_yCM_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr_symmetry");

  p2_vn_yCM_cent_kp->RebinY();
  p2_vn_yCM_cent_km->RebinY();


  TProfile *p_vn_yCM_00to10_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_00to10_pp", 15, 16);
  TProfile *p_vn_yCM_10to40_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_10to40_pp", 9, 14);
  TProfile *p_vn_yCM_40to60_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_40to60_pp", 5, 8);
  
  TProfile *p_vn_yCM_00to10_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_00to10_pm", 15, 16);
  TProfile *p_vn_yCM_10to40_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_10to40_pm", 9, 14);
  TProfile *p_vn_yCM_40to60_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_40to60_pm", 5, 8);

  TProfile *p_vn_yCM_00to10_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_00to10_kp", 15, 16);
  TProfile *p_vn_yCM_10to40_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_10to40_kp", 9, 14);
  TProfile *p_vn_yCM_40to60_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_40to60_kp", 5, 8);

  TProfile *p_vn_yCM_00to10_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_00to10_km", 15, 16);
  TProfile *p_vn_yCM_10to40_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_10to40_km", 9, 14);
  TProfile *p_vn_yCM_40to60_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_40to60_km", 5, 8);

  TProfile *p_vn_yCM_00to10_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_00to10_pr", 15, 16);
  TProfile *p_vn_yCM_10to40_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_10to40_pr", 9, 14);
  TProfile *p_vn_yCM_40to60_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_40to60_pr", 5, 8);

  TProfile *p_vn_yCM_00to10_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_00to10_pr_symm", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_10to40_pr_symm", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_40to60_pr_symm", 5, 8);

  
  h_vn_yCM_00to10_pp = new TH1D("h_vn_yCM_00to10_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_10to40_pp = new TH1D("h_vn_yCM_10to40_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_40to60_pp = new TH1D("h_vn_yCM_40to60_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  h_vn_yCM_00to10_pm = new TH1D("h_vn_yCM_00to10_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_10to40_pm = new TH1D("h_vn_yCM_10to40_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_40to60_pm = new TH1D("h_vn_yCM_40to60_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  h_vn_yCM_00to10_kp = new TH1D("h_vn_yCM_00to10_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  h_vn_yCM_10to40_kp = new TH1D("h_vn_yCM_10to40_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  h_vn_yCM_40to60_kp = new TH1D("h_vn_yCM_40to60_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  h_vn_yCM_00to10_km = new TH1D("h_vn_yCM_00to10_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  h_vn_yCM_10to40_km = new TH1D("h_vn_yCM_10to40_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  h_vn_yCM_40to60_km = new TH1D("h_vn_yCM_40to60_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  h_vn_yCM_00to10_pr = new TH1D("h_vn_yCM_00to10_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_10to40_pr = new TH1D("h_vn_yCM_10to40_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_40to60_pr = new TH1D("h_vn_yCM_40to60_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  h_vn_yCM_00to10_pr_symm = new TH1D("h_vn_yCM_00to10_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_10to40_pr_symm = new TH1D("h_vn_yCM_10to40_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_40to60_pr_symm = new TH1D("h_vn_yCM_40to60_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  
  // Convert profiles to histograms
  h_vn_yCM_00to10_pp = p_vn_yCM_00to10_pp->ProjectionX();
  h_vn_yCM_10to40_pp = p_vn_yCM_10to40_pp->ProjectionX();
  h_vn_yCM_40to60_pp = p_vn_yCM_40to60_pp->ProjectionX();
  h_vn_yCM_00to10_pp = trimRapidityPlot(h_vn_yCM_00to10_pp);
  h_vn_yCM_10to40_pp = trimRapidityPlot(h_vn_yCM_10to40_pp);
  h_vn_yCM_40to60_pp = trimRapidityPlot(h_vn_yCM_40to60_pp);
  
  h_vn_yCM_00to10_pm = p_vn_yCM_00to10_pm->ProjectionX();
  h_vn_yCM_10to40_pm = p_vn_yCM_10to40_pm->ProjectionX();
  h_vn_yCM_40to60_pm = p_vn_yCM_40to60_pm->ProjectionX();
  h_vn_yCM_00to10_pm = trimRapidityPlot(h_vn_yCM_00to10_pm);
  h_vn_yCM_10to40_pm = trimRapidityPlot(h_vn_yCM_10to40_pm);
  h_vn_yCM_40to60_pm = trimRapidityPlot(h_vn_yCM_40to60_pm);
  
  h_vn_yCM_00to10_kp = p_vn_yCM_00to10_kp->ProjectionX();
  h_vn_yCM_10to40_kp = p_vn_yCM_10to40_kp->ProjectionX();
  h_vn_yCM_40to60_kp = p_vn_yCM_40to60_kp->ProjectionX();
  h_vn_yCM_00to10_kp = trimRapidityPlot(h_vn_yCM_00to10_kp);
  h_vn_yCM_10to40_kp = trimRapidityPlot(h_vn_yCM_10to40_kp);
  h_vn_yCM_40to60_kp = trimRapidityPlot(h_vn_yCM_40to60_kp);
  
  h_vn_yCM_00to10_km = p_vn_yCM_00to10_km->ProjectionX();
  h_vn_yCM_10to40_km = p_vn_yCM_10to40_km->ProjectionX();
  h_vn_yCM_40to60_km = p_vn_yCM_40to60_km->ProjectionX();
  h_vn_yCM_00to10_km = trimRapidityPlot(h_vn_yCM_00to10_km);
  h_vn_yCM_10to40_km = trimRapidityPlot(h_vn_yCM_10to40_km);
  h_vn_yCM_40to60_km = trimRapidityPlot(h_vn_yCM_40to60_km);
  
  h_vn_yCM_00to10_pr = p_vn_yCM_00to10_pr->ProjectionX();
  h_vn_yCM_10to40_pr = p_vn_yCM_10to40_pr->ProjectionX();
  h_vn_yCM_40to60_pr = p_vn_yCM_40to60_pr->ProjectionX();

  h_vn_yCM_00to10_pr_symm = p_vn_yCM_00to10_pr_symm->ProjectionX();
  h_vn_yCM_10to40_pr_symm = p_vn_yCM_10to40_pr_symm->ProjectionX();
  h_vn_yCM_40to60_pr_symm = p_vn_yCM_40to60_pr_symm->ProjectionX();



  //=== vn vs pT stuff
  TProfile2D *p2_vn_pT_cent_pp = (TProfile2D*)file->Get("p2_vn_pT_cent_pp");
  TProfile2D *p2_vn_pT_cent_pm = (TProfile2D*)file->Get("p2_vn_pT_cent_pm");
  TProfile2D *p2_vn_pT_cent_kp = (TProfile2D*)file->Get("p2_vn_pT_cent_kp");
  TProfile2D *p2_vn_pT_cent_km = (TProfile2D*)file->Get("p2_vn_pT_cent_km");
  TProfile2D *p2_vn_pT_cent_pr = (TProfile2D*)file->Get("p2_vn_pT_cent_pr");
  //TProfile2D *p2_vn_pT_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_pT_cent_pr_symmetry");

  p2_vn_pT_cent_kp->RebinY();
  p2_vn_pT_cent_km->RebinY();


  TProfile *p_vn_pT_00to10_pp = p2_vn_pT_cent_pp->ProfileY("p_vn_pT_00to10_pp", 15, 16);
  TProfile *p_vn_pT_10to40_pp = p2_vn_pT_cent_pp->ProfileY("p_vn_pT_10to40_pp", 9, 14);
  TProfile *p_vn_pT_40to60_pp = p2_vn_pT_cent_pp->ProfileY("p_vn_pT_40to60_pp", 5, 8);
  
  TProfile *p_vn_pT_00to10_pm = p2_vn_pT_cent_pm->ProfileY("p_vn_pT_00to10_pm", 15, 16);
  TProfile *p_vn_pT_10to40_pm = p2_vn_pT_cent_pm->ProfileY("p_vn_pT_10to40_pm", 9, 14);
  TProfile *p_vn_pT_40to60_pm = p2_vn_pT_cent_pm->ProfileY("p_vn_pT_40to60_pm", 5, 8);

  TProfile *p_vn_pT_00to10_kp = p2_vn_pT_cent_kp->ProfileY("p_vn_pT_00to10_kp", 15, 16);
  TProfile *p_vn_pT_10to40_kp = p2_vn_pT_cent_kp->ProfileY("p_vn_pT_10to40_kp", 9, 14);
  TProfile *p_vn_pT_40to60_kp = p2_vn_pT_cent_kp->ProfileY("p_vn_pT_40to60_kp", 5, 8);

  TProfile *p_vn_pT_00to10_km = p2_vn_pT_cent_km->ProfileY("p_vn_pT_00to10_km", 15, 16);
  TProfile *p_vn_pT_10to40_km = p2_vn_pT_cent_km->ProfileY("p_vn_pT_10to40_km", 9, 14);
  TProfile *p_vn_pT_40to60_km = p2_vn_pT_cent_km->ProfileY("p_vn_pT_40to60_km", 5, 8);

  TProfile *p_vn_pT_00to10_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_00to10_pr", 15, 16);
  TProfile *p_vn_pT_10to40_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_10to40_pr", 9, 14);
  TProfile *p_vn_pT_40to60_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_40to60_pr", 5, 8);


  h_vn_pT_00to10_pp = new TH1D("h_vn_pT_00to10_pp", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_10to40_pp = new TH1D("h_vn_pT_10to40_pp", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_40to60_pp = new TH1D("h_vn_pT_40to60_pp", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  h_vn_pT_00to10_pm = new TH1D("h_vn_pT_00to10_pm", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_10to40_pm = new TH1D("h_vn_pT_10to40_pm", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_40to60_pm = new TH1D("h_vn_pT_40to60_pm", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  h_vn_pT_00to10_kp = new TH1D("h_vn_pT_00to10_kp", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  h_vn_pT_10to40_kp = new TH1D("h_vn_pT_10to40_kp", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  h_vn_pT_40to60_kp = new TH1D("h_vn_pT_40to60_kp", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);

  h_vn_pT_00to10_km = new TH1D("h_vn_pT_00to10_km", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  h_vn_pT_10to40_km = new TH1D("h_vn_pT_10to40_km", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  h_vn_pT_40to60_km = new TH1D("h_vn_pT_40to60_km", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);

  h_vn_pT_00to10_pr = new TH1D("h_vn_pT_00to10_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_10to40_pr = new TH1D("h_vn_pT_10to40_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_40to60_pr = new TH1D("h_vn_pT_40to60_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  
  h_vn_pT_00to10_pp = p_vn_pT_00to10_pp->ProjectionX();
  h_vn_pT_10to40_pp = p_vn_pT_10to40_pp->ProjectionX();
  h_vn_pT_40to60_pp = p_vn_pT_40to60_pp->ProjectionX();

  h_vn_pT_00to10_pm = p_vn_pT_00to10_pm->ProjectionX();
  h_vn_pT_10to40_pm = p_vn_pT_10to40_pm->ProjectionX();
  h_vn_pT_40to60_pm = p_vn_pT_40to60_pm->ProjectionX();

  h_vn_pT_00to10_kp = p_vn_pT_00to10_kp->ProjectionX();
  h_vn_pT_10to40_kp = p_vn_pT_10to40_kp->ProjectionX();
  h_vn_pT_40to60_kp = p_vn_pT_40to60_kp->ProjectionX();

  h_vn_pT_00to10_km = p_vn_pT_00to10_km->ProjectionX();
  h_vn_pT_10to40_km = p_vn_pT_10to40_km->ProjectionX();
  h_vn_pT_40to60_km = p_vn_pT_40to60_km->ProjectionX();

  h_vn_pT_00to10_pr = p_vn_pT_00to10_pr->ProjectionX();
  h_vn_pT_10to40_pr = p_vn_pT_10to40_pr->ProjectionX();
  h_vn_pT_40to60_pr = p_vn_pT_40to60_pr->ProjectionX();
}// End initialize()

// Combine this variation with the normal data to get the differences, nSigma_delta, etc.
void Variation::combine(Variation *normalData)
{
  DataPoint point;

  for (int i = 1; i <= h_vn_pp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pp->GetBinContent(i);
      point.variedError = h_vn_pp->GetBinError(i);
      point.normalValue = normalData->h_vn_pp->GetBinContent(i);
      point.normalError = normalData->h_vn_pp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pp.push_back(point);
    }

  for (int i = 1; i <= h_vn_pm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pm->GetBinContent(i);
      point.variedError = h_vn_pm->GetBinError(i);
      point.normalValue = normalData->h_vn_pm->GetBinContent(i);
      point.normalError = normalData->h_vn_pm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pm.push_back(point);
    }

  for (int i = 1; i <= h_vn_kp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_kp->GetBinContent(i);
      point.variedError = h_vn_kp->GetBinError(i);
      point.normalValue = normalData->h_vn_kp->GetBinContent(i);
      point.normalError = normalData->h_vn_kp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_kp.push_back(point);
    }

  for (int i = 1; i <= h_vn_km->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_km->GetBinContent(i);
      point.variedError = h_vn_km->GetBinError(i);
      point.normalValue = normalData->h_vn_km->GetBinContent(i);
      point.normalError = normalData->h_vn_km->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_km.push_back(point);
    }

  for (int i = 1; i <= h_vn_pr->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pr->GetBinContent(i);
      point.variedError = h_vn_pr->GetBinError(i);
      point.normalValue = normalData->h_vn_pr->GetBinContent(i);
      point.normalError = normalData->h_vn_pr->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pr.push_back(point);
    }


  for (int i = 1; i <= h_vn_pp_ext->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pp_ext->GetBinContent(i);
      point.variedError = h_vn_pp_ext->GetBinError(i);
      point.normalValue = normalData->h_vn_pp_ext->GetBinContent(i);
      point.normalError = normalData->h_vn_pp_ext->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pp_ext.push_back(point);
    }

  for (int i = 1; i <= h_vn_pm_ext->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pm_ext->GetBinContent(i);
      point.variedError = h_vn_pm_ext->GetBinError(i);
      point.normalValue = normalData->h_vn_pm_ext->GetBinContent(i);
      point.normalError = normalData->h_vn_pm_ext->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pm_ext.push_back(point);
    }

  for (int i = 1; i <= h_vn_kp_ext->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_kp_ext->GetBinContent(i);
      point.variedError = h_vn_kp_ext->GetBinError(i);
      point.normalValue = normalData->h_vn_kp_ext->GetBinContent(i);
      point.normalError = normalData->h_vn_kp_ext->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_kp_ext.push_back(point);
    }

  for (int i = 1; i <= h_vn_km_ext->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_km_ext->GetBinContent(i);
      point.variedError = h_vn_km_ext->GetBinError(i);
      point.normalValue = normalData->h_vn_km_ext->GetBinContent(i);
      point.normalError = normalData->h_vn_km_ext->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_km_ext.push_back(point);
    }

  for (int i = 1; i <= h_vn_pr_ext->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pr_ext->GetBinContent(i);
      point.variedError = h_vn_pr_ext->GetBinError(i);
      point.normalValue = normalData->h_vn_pr_ext->GetBinContent(i);
      point.normalError = normalData->h_vn_pr_ext->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pr_ext.push_back(point);
    }

  for (int i = 1; i <= h_vn_pr_for->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pr_for->GetBinContent(i);
      point.variedError = h_vn_pr_for->GetBinError(i);
      point.normalValue = normalData->h_vn_pr_for->GetBinContent(i);
      point.normalError = normalData->h_vn_pr_for->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pr_for.push_back(point);
    }

  //== vn vs rapidity - pi+
  for (int i = 1; i <= h_vn_yCM_00to10_pp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_00to10_pp->GetBinContent(i);
      point.variedError = h_vn_yCM_00to10_pp->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_00to10_pp->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_00to10_pp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_00to10_pp.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_10to40_pp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_10to40_pp->GetBinContent(i);
      point.variedError = h_vn_yCM_10to40_pp->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_10to40_pp->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_10to40_pp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_10to40_pp.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_40to60_pp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_40to60_pp->GetBinContent(i);
      point.variedError = h_vn_yCM_40to60_pp->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_40to60_pp->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_40to60_pp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_40to60_pp.push_back(point);
    }
  //
  
  //== vn vs rapidity - pi-
  for (int i = 1; i <= h_vn_yCM_00to10_pm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_00to10_pm->GetBinContent(i);
      point.variedError = h_vn_yCM_00to10_pm->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_00to10_pm->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_00to10_pm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_00to10_pm.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_10to40_pm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_10to40_pm->GetBinContent(i);
      point.variedError = h_vn_yCM_10to40_pm->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_10to40_pm->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_10to40_pm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_10to40_pm.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_40to60_pm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_40to60_pm->GetBinContent(i);
      point.variedError = h_vn_yCM_40to60_pm->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_40to60_pm->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_40to60_pm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_40to60_pm.push_back(point);
    }
  //
  
  //== vn vs rapidity - K+
  for (int i = 1; i <= h_vn_yCM_00to10_kp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_00to10_kp->GetBinContent(i);
      point.variedError = h_vn_yCM_00to10_kp->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_00to10_kp->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_00to10_kp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_00to10_kp.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_10to40_kp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_10to40_kp->GetBinContent(i);
      point.variedError = h_vn_yCM_10to40_kp->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_10to40_kp->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_10to40_kp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_10to40_kp.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_40to60_kp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_40to60_kp->GetBinContent(i);
      point.variedError = h_vn_yCM_40to60_kp->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_40to60_kp->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_40to60_kp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_40to60_kp.push_back(point);
    }
  //

  //== vn vs rapidity - K-
  for (int i = 1; i <= h_vn_yCM_00to10_km->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_00to10_km->GetBinContent(i);
      point.variedError = h_vn_yCM_00to10_km->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_00to10_km->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_00to10_km->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_00to10_km.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_10to40_km->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_10to40_km->GetBinContent(i);
      point.variedError = h_vn_yCM_10to40_km->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_10to40_km->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_10to40_km->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_10to40_km.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_40to60_km->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_40to60_km->GetBinContent(i);
      point.variedError = h_vn_yCM_40to60_km->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_40to60_km->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_40to60_km->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_40to60_km.push_back(point);
    }
  //
  
  //== vn vs rapidity - Proton
  for (int i = 1; i <= h_vn_yCM_00to10_pr->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_00to10_pr->GetBinContent(i);
      point.variedError = h_vn_yCM_00to10_pr->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_00to10_pr->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_00to10_pr->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_00to10_pr.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_10to40_pr->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_10to40_pr->GetBinContent(i);
      point.variedError = h_vn_yCM_10to40_pr->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_10to40_pr->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_10to40_pr->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_10to40_pr.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_40to60_pr->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_40to60_pr->GetBinContent(i);
      point.variedError = h_vn_yCM_40to60_pr->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_40to60_pr->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_40to60_pr->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_40to60_pr.push_back(point);
    }
  //

  //== vn vs rapidity - Proton symmetric across midrapidity
  for (int i = 1; i <= h_vn_yCM_00to10_pr_symm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_00to10_pr_symm->GetBinContent(i);
      point.variedError = h_vn_yCM_00to10_pr_symm->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_00to10_pr_symm->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_00to10_pr_symm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_00to10_pr_symm.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_10to40_pr_symm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_10to40_pr_symm->GetBinContent(i);
      point.variedError = h_vn_yCM_10to40_pr_symm->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_10to40_pr_symm->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_10to40_pr_symm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_10to40_pr_symm.push_back(point);
    }

  for (int i = 1; i <= h_vn_yCM_40to60_pr_symm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_yCM_40to60_pr_symm->GetBinContent(i);
      point.variedError = h_vn_yCM_40to60_pr_symm->GetBinError(i);
      point.normalValue = normalData->h_vn_yCM_40to60_pr_symm->GetBinContent(i);
      point.normalError = normalData->h_vn_yCM_40to60_pr_symm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_yCM_40to60_pr_symm.push_back(point);
    }
  //



  //== vn vs pT - pi+
  for (int i = 1; i <= h_vn_pT_00to10_pp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_00to10_pp->GetBinContent(i);
      point.variedError = h_vn_pT_00to10_pp->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_00to10_pp->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_00to10_pp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_00to10_pp.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_10to40_pp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_10to40_pp->GetBinContent(i);
      point.variedError = h_vn_pT_10to40_pp->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_10to40_pp->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_10to40_pp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_10to40_pp.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_40to60_pp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_40to60_pp->GetBinContent(i);
      point.variedError = h_vn_pT_40to60_pp->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_40to60_pp->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_40to60_pp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_40to60_pp.push_back(point);
    }
  //
  
  //== vn vs pT - pi-
  for (int i = 1; i <= h_vn_pT_00to10_pm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_00to10_pm->GetBinContent(i);
      point.variedError = h_vn_pT_00to10_pm->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_00to10_pm->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_00to10_pm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_00to10_pm.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_10to40_pm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_10to40_pm->GetBinContent(i);
      point.variedError = h_vn_pT_10to40_pm->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_10to40_pm->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_10to40_pm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_10to40_pm.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_40to60_pm->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_40to60_pm->GetBinContent(i);
      point.variedError = h_vn_pT_40to60_pm->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_40to60_pm->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_40to60_pm->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_40to60_pm.push_back(point);
    }
  //
  
  //== vn vs pT - K+
  for (int i = 1; i <= h_vn_pT_00to10_kp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_00to10_kp->GetBinContent(i);
      point.variedError = h_vn_pT_00to10_kp->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_00to10_kp->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_00to10_kp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_00to10_kp.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_10to40_kp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_10to40_kp->GetBinContent(i);
      point.variedError = h_vn_pT_10to40_kp->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_10to40_kp->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_10to40_kp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_10to40_kp.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_40to60_kp->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_40to60_kp->GetBinContent(i);
      point.variedError = h_vn_pT_40to60_kp->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_40to60_kp->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_40to60_kp->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_40to60_kp.push_back(point);
    }
  //

  //== vn vs pT - K-
  for (int i = 1; i <= h_vn_pT_00to10_km->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_00to10_km->GetBinContent(i);
      point.variedError = h_vn_pT_00to10_km->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_00to10_km->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_00to10_km->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_00to10_km.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_10to40_km->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_10to40_km->GetBinContent(i);
      point.variedError = h_vn_pT_10to40_km->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_10to40_km->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_10to40_km->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_10to40_km.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_40to60_km->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_40to60_km->GetBinContent(i);
      point.variedError = h_vn_pT_40to60_km->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_40to60_km->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_40to60_km->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_40to60_km.push_back(point);
    }
  //
  
  //== vn vs pT - Proton
  for (int i = 1; i <= h_vn_pT_00to10_pr->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_00to10_pr->GetBinContent(i);
      point.variedError = h_vn_pT_00to10_pr->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_00to10_pr->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_00to10_pr->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_00to10_pr.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_10to40_pr->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_10to40_pr->GetBinContent(i);
      point.variedError = h_vn_pT_10to40_pr->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_10to40_pr->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_10to40_pr->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_10to40_pr.push_back(point);
    }

  for (int i = 1; i <= h_vn_pT_40to60_pr->GetNbinsX(); i++)
    {
      point.variedValue = h_vn_pT_40to60_pr->GetBinContent(i);
      point.variedError = h_vn_pT_40to60_pr->GetBinError(i);
      point.normalValue = normalData->h_vn_pT_40to60_pr->GetBinContent(i);
      point.normalError = normalData->h_vn_pT_40to60_pr->GetBinError(i);
      point.delta = point.variedValue - point.normalValue;
      point.deltaError = TMath::Sqrt(TMath::Power(point.variedError, 2) - TMath::Power(point.normalError, 2));
      point.nSigma = (point.deltaError == 0.0)?-1:TMath::Abs(point.delta/point.deltaError);
      v_vn_pT_40to60_pr.push_back(point);
    }
  //
}// End combine()

// Get all standard deviations and variances needed for the systematic uncertainties.
void Variation::stdDevs()
{
  TH1D *temp;
  TString histName;
  TH1D *sigmas = new TH1D("sigmas", ";#Delta/#sigma_{#Delta};Bins", 500, 0, 15);
  TH1D *deltaErrors = new TH1D("deltaErrors", ";#sigma_{#Delta};Bins", 500, -0.5, 0.5);
  TFile *detailsFile = new TFile("Details_"+ID+".root", "RECREATE");
  detailsFile->cd();
  
  for (int i = 0; i < v_vn_pp.size(); i++)
    {
      if (v_vn_pp.at(i).nSigma != -1){
	  sigmas->Fill(v_vn_pp.at(i).nSigma);
	  deltaErrors->Fill(v_vn_pp.at(i).deltaError);
	}
      
      histName.Form((TString)h_vn_pp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pp.at(i).variedValue);
      temp->Fill(v_vn_pp.at(i).normalValue);
      temp->Write();
      
      v_vn_pp.at(i).stdDev   = temp->GetStdDev();
      v_vn_pp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pm.size(); i++)
    {
      if (v_vn_pm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pm.at(i).nSigma);
	deltaErrors->Fill(v_vn_pm.at(i).deltaError);
      }

      histName.Form((TString)h_vn_pm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pm.at(i).variedValue);
      temp->Fill(v_vn_pm.at(i).normalValue);
      temp->Write();

      v_vn_pm.at(i).stdDev   = temp->GetStdDev();
      v_vn_pm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_kp.size(); i++)
    {
      if (v_vn_kp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_kp.at(i).nSigma);
	deltaErrors->Fill(v_vn_kp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_kp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_kp.at(i).variedValue);
      temp->Fill(v_vn_kp.at(i).normalValue);
      temp->Write();

      v_vn_kp.at(i).stdDev   = temp->GetStdDev();
      v_vn_kp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_km.size(); i++)
    {
      if (v_vn_km.at(i).nSigma != -1){
	sigmas->Fill(v_vn_km.at(i).nSigma);
	deltaErrors->Fill(v_vn_km.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_km->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_km.at(i).variedValue);
      temp->Fill(v_vn_km.at(i).normalValue);
      temp->Write();

      v_vn_km.at(i).stdDev   = temp->GetStdDev();
      v_vn_km.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pr.size(); i++)
    {
      if (v_vn_pr.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pr.at(i).nSigma);
	deltaErrors->Fill(v_vn_pr.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pr->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pr.at(i).variedValue);
      temp->Fill(v_vn_pr.at(i).normalValue);
      temp->Write();

      v_vn_pr.at(i).stdDev   = temp->GetStdDev();
      v_vn_pr.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }


  for (int i = 0; i < v_vn_pp_ext.size(); i++)
    {
      if (v_vn_pp_ext.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pp_ext.at(i).nSigma);
	deltaErrors->Fill(v_vn_pp_ext.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pp_ext->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pp_ext.at(i).variedValue);
      temp->Fill(v_vn_pp_ext.at(i).normalValue);
      temp->Write();

      v_vn_pp_ext.at(i).stdDev   = temp->GetStdDev();
      v_vn_pp_ext.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pm_ext.size(); i++)
    {
      if (v_vn_pm_ext.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pm_ext.at(i).nSigma);
	deltaErrors->Fill(v_vn_pm_ext.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pm_ext->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pm_ext.at(i).variedValue);
      temp->Fill(v_vn_pm_ext.at(i).normalValue);
      temp->Write();

      v_vn_pm_ext.at(i).stdDev   = temp->GetStdDev();
      v_vn_pm_ext.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_kp_ext.size(); i++)
    {
      if (v_vn_kp_ext.at(i).nSigma != -1){
	sigmas->Fill(v_vn_kp_ext.at(i).nSigma);
	deltaErrors->Fill(v_vn_kp_ext.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_kp_ext->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_kp_ext.at(i).variedValue);
      temp->Fill(v_vn_kp_ext.at(i).normalValue);
      temp->Write();

      v_vn_kp_ext.at(i).stdDev   = temp->GetStdDev();
      v_vn_kp_ext.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_km_ext.size(); i++)
    {
      if (v_vn_km_ext.at(i).nSigma != -1){
	sigmas->Fill(v_vn_km_ext.at(i).nSigma);
	deltaErrors->Fill(v_vn_km_ext.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_km_ext->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_km_ext.at(i).variedValue);
      temp->Fill(v_vn_km_ext.at(i).normalValue);
      temp->Write();

      v_vn_km_ext.at(i).stdDev   = temp->GetStdDev();
      v_vn_km_ext.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pr_ext.size(); i++)
    {
      if (v_vn_pr_ext.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pr_ext.at(i).nSigma);
	deltaErrors->Fill(v_vn_pr_ext.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pr_ext->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pr_ext.at(i).variedValue);
      temp->Fill(v_vn_pr_ext.at(i).normalValue);
      temp->Write();

      v_vn_pr_ext.at(i).stdDev   = temp->GetStdDev();
      v_vn_pr_ext.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pr_for.size(); i++)
    {
      if (v_vn_pr_for.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pr_for.at(i).nSigma);
	deltaErrors->Fill(v_vn_pr_for.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pr_for->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pr_for.at(i).variedValue);
      temp->Fill(v_vn_pr_for.at(i).normalValue);
      temp->Write();

      v_vn_pr_for.at(i).stdDev   = temp->GetStdDev();
      v_vn_pr_for.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  //=== vn vs rapidity - pi+
  for (int i = 0; i < v_vn_yCM_00to10_pp.size(); i++)
    {
      if (v_vn_yCM_00to10_pp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_00to10_pp.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_00to10_pp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_00to10_pp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_00to10_pp.at(i).variedValue);
      temp->Fill(v_vn_yCM_00to10_pp.at(i).normalValue);
      temp->Write();

      v_vn_yCM_00to10_pp.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_00to10_pp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_yCM_10to40_pp.size(); i++)
    {
      if (v_vn_yCM_10to40_pp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_10to40_pp.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_10to40_pp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_10to40_pp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_10to40_pp.at(i).variedValue);
      temp->Fill(v_vn_yCM_10to40_pp.at(i).normalValue);
      temp->Write();

      v_vn_yCM_10to40_pp.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_10to40_pp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_yCM_40to60_pp.size(); i++)
    {
      if (v_vn_yCM_40to60_pp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_40to60_pp.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_40to60_pp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_40to60_pp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_40to60_pp.at(i).variedValue);
      temp->Fill(v_vn_yCM_40to60_pp.at(i).normalValue);
      temp->Write();

      v_vn_yCM_40to60_pp.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_40to60_pp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs rapidity - pi-
  for (int i = 0; i < v_vn_yCM_00to10_pm.size(); i++)
    {
      if (v_vn_yCM_00to10_pm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_00to10_pm.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_00to10_pm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_00to10_pm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_00to10_pm.at(i).variedValue);
      temp->Fill(v_vn_yCM_00to10_pm.at(i).normalValue);
      temp->Write();

      v_vn_yCM_00to10_pm.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_00to10_pm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_yCM_10to40_pm.size(); i++)
    {
      if (v_vn_yCM_10to40_pm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_10to40_pm.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_10to40_pm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_10to40_pm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_10to40_pm.at(i).variedValue);
      temp->Fill(v_vn_yCM_10to40_pm.at(i).normalValue);
      temp->Write();

      v_vn_yCM_10to40_pm.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_10to40_pm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_yCM_40to60_pm.size(); i++)
    {
      if (v_vn_yCM_40to60_pm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_40to60_pm.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_40to60_pm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_40to60_pm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_40to60_pm.at(i).variedValue);
      temp->Fill(v_vn_yCM_40to60_pm.at(i).normalValue);
      temp->Write();

      v_vn_yCM_40to60_pm.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_40to60_pm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs rapidity - K+
  for (int i = 0; i < v_vn_yCM_00to10_kp.size(); i++)
    {
      if (v_vn_yCM_00to10_kp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_00to10_kp.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_00to10_kp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_00to10_kp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_00to10_kp.at(i).variedValue);
      temp->Fill(v_vn_yCM_00to10_kp.at(i).normalValue);
      temp->Write();

      v_vn_yCM_00to10_kp.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_00to10_kp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_yCM_10to40_kp.size(); i++)
    {
      if (v_vn_yCM_10to40_kp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_10to40_kp.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_10to40_kp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_10to40_kp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_10to40_kp.at(i).variedValue);
      temp->Fill(v_vn_yCM_10to40_kp.at(i).normalValue);
      temp->Write();

      v_vn_yCM_10to40_kp.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_10to40_kp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_yCM_40to60_kp.size(); i++)
    {
      if (v_vn_yCM_40to60_kp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_40to60_kp.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_40to60_kp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_40to60_kp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_40to60_kp.at(i).variedValue);
      temp->Fill(v_vn_yCM_40to60_kp.at(i).normalValue);
      temp->Write();

      v_vn_yCM_40to60_kp.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_40to60_kp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs rapidity - K-
  for (int i = 0; i < v_vn_yCM_00to10_km.size(); i++)
    {
      /*
      if (v_vn_yCM_00to10_km.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_00to10_km.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_00to10_km.at(i).deltaError);
	}*/
      
      histName.Form((TString)h_vn_yCM_00to10_km->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_00to10_km.at(i).variedValue);
      temp->Fill(v_vn_yCM_00to10_km.at(i).normalValue);
      temp->Write();

      v_vn_yCM_00to10_km.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_00to10_km.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_yCM_10to40_km.size(); i++)
    {
      if (v_vn_yCM_10to40_km.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_10to40_km.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_10to40_km.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_10to40_km->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_10to40_km.at(i).variedValue);
      temp->Fill(v_vn_yCM_10to40_km.at(i).normalValue);
      temp->Write();

      v_vn_yCM_10to40_km.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_10to40_km.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_yCM_40to60_km.size(); i++)
    {
      /*
      if (v_vn_yCM_40to60_km.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_40to60_km.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_40to60_km.at(i).deltaError);
	}*/
      
      histName.Form((TString)h_vn_yCM_40to60_km->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_40to60_km.at(i).variedValue);
      temp->Fill(v_vn_yCM_40to60_km.at(i).normalValue);
      temp->Write();

      v_vn_yCM_40to60_km.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_40to60_km.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs rapidity - Proton
  for (int i = 0; i < v_vn_yCM_00to10_pr.size(); i++)
    {
      if (v_vn_yCM_00to10_pr.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_00to10_pr.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_00to10_pr.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_00to10_pr->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_00to10_pr.at(i).variedValue);
      temp->Fill(v_vn_yCM_00to10_pr.at(i).normalValue);
      temp->Write();

      v_vn_yCM_00to10_pr.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_00to10_pr.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_yCM_10to40_pr.size(); i++)
    {
      if (v_vn_yCM_10to40_pr.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_10to40_pr.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_10to40_pr.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_10to40_pr->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_10to40_pr.at(i).variedValue);
      temp->Fill(v_vn_yCM_10to40_pr.at(i).normalValue);
      temp->Write();

      v_vn_yCM_10to40_pr.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_10to40_pr.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_yCM_40to60_pr.size(); i++)
    {
      if (v_vn_yCM_40to60_pr.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_40to60_pr.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_40to60_pr.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_40to60_pr->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_40to60_pr.at(i).variedValue);
      temp->Fill(v_vn_yCM_40to60_pr.at(i).normalValue);
      temp->Write();

      v_vn_yCM_40to60_pr.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_40to60_pr.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs rapidity - Proton symmetric across midrapidity
  for (int i = 0; i < v_vn_yCM_00to10_pr_symm.size(); i++)
    {
      if (v_vn_yCM_00to10_pr_symm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_00to10_pr_symm.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_00to10_pr_symm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_00to10_pr_symm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_00to10_pr_symm.at(i).variedValue);
      temp->Fill(v_vn_yCM_00to10_pr_symm.at(i).normalValue);
      temp->Write();

      v_vn_yCM_00to10_pr_symm.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_00to10_pr_symm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_yCM_10to40_pr_symm.size(); i++)
    {
      if (v_vn_yCM_10to40_pr_symm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_10to40_pr_symm.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_10to40_pr_symm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_10to40_pr_symm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_10to40_pr_symm.at(i).variedValue);
      temp->Fill(v_vn_yCM_10to40_pr_symm.at(i).normalValue);
      temp->Write();

      v_vn_yCM_10to40_pr_symm.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_10to40_pr_symm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_yCM_40to60_pr_symm.size(); i++)
    {
      if (v_vn_yCM_40to60_pr_symm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_yCM_40to60_pr_symm.at(i).nSigma);
	deltaErrors->Fill(v_vn_yCM_40to60_pr_symm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_yCM_40to60_pr_symm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_yCM_40to60_pr_symm.at(i).variedValue);
      temp->Fill(v_vn_yCM_40to60_pr_symm.at(i).normalValue);
      temp->Write();

      v_vn_yCM_40to60_pr_symm.at(i).stdDev   = temp->GetStdDev();
      v_vn_yCM_40to60_pr_symm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //


  //=== vn vs pT - pi+
  for (int i = 0; i < v_vn_pT_00to10_pp.size(); i++)
    {
      if (v_vn_pT_00to10_pp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_00to10_pp.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_00to10_pp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_00to10_pp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_00to10_pp.at(i).variedValue);
      temp->Fill(v_vn_pT_00to10_pp.at(i).normalValue);
      temp->Write();

      v_vn_pT_00to10_pp.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_00to10_pp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_pT_10to40_pp.size(); i++)
    {
      if (v_vn_pT_10to40_pp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_10to40_pp.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_10to40_pp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_10to40_pp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_10to40_pp.at(i).variedValue);
      temp->Fill(v_vn_pT_10to40_pp.at(i).normalValue);
      temp->Write();

      v_vn_pT_10to40_pp.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_10to40_pp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pT_40to60_pp.size(); i++)
    {
      if (v_vn_pT_40to60_pp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_40to60_pp.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_40to60_pp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_40to60_pp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_40to60_pp.at(i).variedValue);
      temp->Fill(v_vn_pT_40to60_pp.at(i).normalValue);
      temp->Write();

      v_vn_pT_40to60_pp.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_40to60_pp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs pT - pi-
  for (int i = 0; i < v_vn_pT_00to10_pm.size(); i++)
    {
      if (v_vn_pT_00to10_pm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_00to10_pm.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_00to10_pm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_00to10_pm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_00to10_pm.at(i).variedValue);
      temp->Fill(v_vn_pT_00to10_pm.at(i).normalValue);
      temp->Write();

      v_vn_pT_00to10_pm.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_00to10_pm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_pT_10to40_pm.size(); i++)
    {
      if (v_vn_pT_10to40_pm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_10to40_pm.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_10to40_pm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_10to40_pm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_10to40_pm.at(i).variedValue);
      temp->Fill(v_vn_pT_10to40_pm.at(i).normalValue);
      temp->Write();

      v_vn_pT_10to40_pm.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_10to40_pm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pT_40to60_pm.size(); i++)
    {
      if (v_vn_pT_40to60_pm.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_40to60_pm.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_40to60_pm.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_40to60_pm->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_40to60_pm.at(i).variedValue);
      temp->Fill(v_vn_pT_40to60_pm.at(i).normalValue);
      temp->Write();

      v_vn_pT_40to60_pm.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_40to60_pm.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs pT - K+
  for (int i = 0; i < v_vn_pT_00to10_kp.size(); i++)
    {
      if (v_vn_pT_00to10_kp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_00to10_kp.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_00to10_kp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_00to10_kp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_00to10_kp.at(i).variedValue);
      temp->Fill(v_vn_pT_00to10_kp.at(i).normalValue);
      temp->Write();

      v_vn_pT_00to10_kp.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_00to10_kp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_pT_10to40_kp.size(); i++)
    {
      if (v_vn_pT_10to40_kp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_10to40_kp.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_10to40_kp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_10to40_kp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_10to40_kp.at(i).variedValue);
      temp->Fill(v_vn_pT_10to40_kp.at(i).normalValue);
      temp->Write();

      v_vn_pT_10to40_kp.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_10to40_kp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pT_40to60_kp.size(); i++)
    {
      if (v_vn_pT_40to60_kp.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_40to60_kp.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_40to60_kp.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_40to60_kp->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_40to60_kp.at(i).variedValue);
      temp->Fill(v_vn_pT_40to60_kp.at(i).normalValue);
      temp->Write();

      v_vn_pT_40to60_kp.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_40to60_kp.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs pT - K-
  for (int i = 0; i < v_vn_pT_00to10_km.size(); i++)
    {
      if (v_vn_pT_00to10_km.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_00to10_km.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_00to10_km.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_00to10_km->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_00to10_km.at(i).variedValue);
      temp->Fill(v_vn_pT_00to10_km.at(i).normalValue);
      temp->Write();

      v_vn_pT_00to10_km.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_00to10_km.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_pT_10to40_km.size(); i++)
    {
      if (v_vn_pT_10to40_km.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_10to40_km.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_10to40_km.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_10to40_km->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_10to40_km.at(i).variedValue);
      temp->Fill(v_vn_pT_10to40_km.at(i).normalValue);
      temp->Write();

      v_vn_pT_10to40_km.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_10to40_km.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pT_40to60_km.size(); i++)
    {
      if (v_vn_pT_40to60_km.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_40to60_km.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_40to60_km.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_40to60_km->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_40to60_km.at(i).variedValue);
      temp->Fill(v_vn_pT_40to60_km.at(i).normalValue);
      temp->Write();

      v_vn_pT_40to60_km.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_40to60_km.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  //=== vn vs pT - Proton
  for (int i = 0; i < v_vn_pT_00to10_pr.size(); i++)
    {
      if (v_vn_pT_00to10_pr.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_00to10_pr.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_00to10_pr.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_00to10_pr->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_00to10_pr.at(i).variedValue);
      temp->Fill(v_vn_pT_00to10_pr.at(i).normalValue);
      temp->Write();

      v_vn_pT_00to10_pr.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_00to10_pr.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  
  for (int i = 0; i < v_vn_pT_10to40_pr.size(); i++)
    {
      if (v_vn_pT_10to40_pr.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_10to40_pr.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_10to40_pr.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_10to40_pr->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_10to40_pr.at(i).variedValue);
      temp->Fill(v_vn_pT_10to40_pr.at(i).normalValue);
      temp->Write();

      v_vn_pT_10to40_pr.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_10to40_pr.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }

  for (int i = 0; i < v_vn_pT_40to60_pr.size(); i++)
    {
      if (v_vn_pT_40to60_pr.at(i).nSigma != -1){
	sigmas->Fill(v_vn_pT_40to60_pr.at(i).nSigma);
	deltaErrors->Fill(v_vn_pT_40to60_pr.at(i).deltaError);
      }
      
      histName.Form((TString)h_vn_pT_40to60_pr->GetName()+"_%d", i+1);
      temp = new TH1D(histName, histName, 500, 1, 1);
      temp->Fill(v_vn_pT_40to60_pr.at(i).variedValue);
      temp->Fill(v_vn_pT_40to60_pr.at(i).normalValue);
      temp->Write();

      v_vn_pT_40to60_pr.at(i).stdDev   = temp->GetStdDev();
      v_vn_pT_40to60_pr.at(i).variance = temp->GetStdDev() * temp->GetStdDev();
      delete temp;
    }
  //

  sigmas->Write();
  deltaErrors->Write();
  
  detailsFile->Close();
  delete detailsFile;
}// End stdDevs()

// Use this with the "Normal" Variation to clean up the plots
void Variation::fixAttributes(TString order_n_str)
{
  if (order_n_str == "3")
    {
      Double_t centralityUpperBounds = 0.15;
      Double_t centralityLowerBounds = -0.15;
      
      h_vn_pp->SetMarkerStyle(20);
      h_vn_pp->SetMarkerSize(2.5);
      h_vn_pp->SetMarkerColor(2);
      h_vn_pp->SetLineColor(2);
      h_vn_pp->SetLineWidth(3);
      h_vn_pp->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pp->GetXaxis()->SetNdivisions(210);
      //h_vn_pp->SetMaximum(centralityUpperBounds);
      //h_vn_pp->SetMinimum(centralityLowerBounds);

      h_vn_pm->SetMarkerStyle(20);
      h_vn_pm->SetMarkerSize(2.5);
      h_vn_pm->SetMarkerColor(4);
      h_vn_pm->SetLineColor(4);
      h_vn_pm->SetLineWidth(3);
      h_vn_pm->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pm->GetXaxis()->SetNdivisions(210);
      //h_vn_pm->SetMaximum(centralityUpperBounds);
      //h_vn_pm->SetMinimum(centralityLowerBounds);

      h_vn_kp->SetMarkerStyle(20);
      h_vn_kp->SetMarkerSize(2.5);
      h_vn_kp->SetMarkerColor(2);
      h_vn_kp->SetLineColor(2);
      h_vn_kp->SetLineWidth(3);
      h_vn_kp->GetYaxis()->SetTitleOffset(1.7);
      h_vn_kp->GetXaxis()->SetNdivisions(210);
      //h_vn_kp->SetMaximum(centralityUpperBounds);
      //h_vn_kp->SetMinimum(centralityLowerBounds);

      h_vn_km->SetMarkerStyle(20);
      h_vn_km->SetMarkerSize(2.5);
      h_vn_km->SetMarkerColor(4);
      h_vn_km->SetLineColor(4);
      h_vn_km->SetLineWidth(3);
      h_vn_km->GetYaxis()->SetTitleOffset(1.7);
      h_vn_km->GetXaxis()->SetNdivisions(210);
      //h_vn_km->SetMaximum(centralityUpperBounds);
      //h_vn_km->SetMinimum(centralityLowerBounds);

  
      h_vn_pr->SetTitle(";Centrality (%);v_{"+order_n_str+"}");
      h_vn_pr->SetMarkerStyle(20);
      h_vn_pr->SetMarkerSize(2.5);
      h_vn_pr->SetMarkerColor(2);
      h_vn_pr->SetLineColor(2);
      h_vn_pr->SetLineWidth(3);
      h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pr->GetXaxis()->SetNdivisions(210);
      //h_vn_pr->SetMaximum(centralityUpperBounds);
      //h_vn_pr->SetMinimum(centralityLowerBounds);



      //=== vn vs rapidity
      h_vn_yCM_00to10_pp->SetMarkerStyle(20);
      h_vn_yCM_10to40_pp->SetMarkerStyle(20);
      h_vn_yCM_40to60_pp->SetMarkerStyle(20);
      h_vn_yCM_00to10_pp->SetMarkerColor(2);
      h_vn_yCM_10to40_pp->SetMarkerColor(4);
      h_vn_yCM_40to60_pp->SetMarkerColor(8);
      h_vn_yCM_00to10_pp->SetMarkerSize(2);
      h_vn_yCM_10to40_pp->SetMarkerSize(2);
      h_vn_yCM_40to60_pp->SetMarkerSize(2);
      h_vn_yCM_00to10_pp->SetLineColor(2);
      h_vn_yCM_10to40_pp->SetLineColor(4);
      h_vn_yCM_40to60_pp->SetLineColor(8);
      h_vn_yCM_00to10_pp->SetLineWidth(3);
      h_vn_yCM_10to40_pp->SetLineWidth(3);
      h_vn_yCM_40to60_pp->SetLineWidth(3);

      h_vn_yCM_00to10_pm->SetMarkerStyle(20);
      h_vn_yCM_10to40_pm->SetMarkerStyle(20);
      h_vn_yCM_40to60_pm->SetMarkerStyle(20);
      h_vn_yCM_00to10_pm->SetMarkerColor(2);
      h_vn_yCM_10to40_pm->SetMarkerColor(4);
      h_vn_yCM_40to60_pm->SetMarkerColor(8);
      h_vn_yCM_00to10_pm->SetMarkerSize(2);
      h_vn_yCM_10to40_pm->SetMarkerSize(2);
      h_vn_yCM_40to60_pm->SetMarkerSize(2);
      h_vn_yCM_00to10_pm->SetLineColor(2);
      h_vn_yCM_10to40_pm->SetLineColor(4);
      h_vn_yCM_40to60_pm->SetLineColor(8);
      h_vn_yCM_00to10_pm->SetLineWidth(3);
      h_vn_yCM_10to40_pm->SetLineWidth(3);
      h_vn_yCM_40to60_pm->SetLineWidth(3);

      h_vn_yCM_00to10_kp->SetMarkerStyle(20);
      h_vn_yCM_10to40_kp->SetMarkerStyle(20);
      h_vn_yCM_40to60_kp->SetMarkerStyle(20);
      h_vn_yCM_00to10_kp->SetMarkerColor(2);
      h_vn_yCM_10to40_kp->SetMarkerColor(4);
      h_vn_yCM_40to60_kp->SetMarkerColor(8);
      h_vn_yCM_00to10_kp->SetMarkerSize(2);
      h_vn_yCM_10to40_kp->SetMarkerSize(2);
      h_vn_yCM_40to60_kp->SetMarkerSize(2);
      h_vn_yCM_00to10_kp->SetLineColor(2);
      h_vn_yCM_10to40_kp->SetLineColor(4);
      h_vn_yCM_40to60_kp->SetLineColor(8);
      h_vn_yCM_00to10_kp->SetLineWidth(3);
      h_vn_yCM_10to40_kp->SetLineWidth(3);
      h_vn_yCM_40to60_kp->SetLineWidth(3);

      h_vn_yCM_00to10_km->SetMarkerStyle(20);
      h_vn_yCM_10to40_km->SetMarkerStyle(20);
      h_vn_yCM_40to60_km->SetMarkerStyle(20);
      h_vn_yCM_00to10_km->SetMarkerColor(2);
      h_vn_yCM_10to40_km->SetMarkerColor(4);
      h_vn_yCM_40to60_km->SetMarkerColor(8);
      h_vn_yCM_00to10_km->SetMarkerSize(2);
      h_vn_yCM_10to40_km->SetMarkerSize(2);
      h_vn_yCM_40to60_km->SetMarkerSize(2);
      h_vn_yCM_00to10_km->SetLineColor(2);
      h_vn_yCM_10to40_km->SetLineColor(4);
      h_vn_yCM_40to60_km->SetLineColor(8);
      h_vn_yCM_00to10_km->SetLineWidth(3);
      h_vn_yCM_10to40_km->SetLineWidth(3);
      h_vn_yCM_40to60_km->SetLineWidth(3);

      h_vn_yCM_00to10_pr->SetMarkerStyle(20);
      h_vn_yCM_10to40_pr->SetMarkerStyle(20);
      h_vn_yCM_40to60_pr->SetMarkerStyle(20);
      h_vn_yCM_00to10_pr->SetMarkerColor(2);
      h_vn_yCM_10to40_pr->SetMarkerColor(4);
      h_vn_yCM_40to60_pr->SetMarkerColor(8);
      h_vn_yCM_00to10_pr->SetMarkerSize(2);
      h_vn_yCM_10to40_pr->SetMarkerSize(2);
      h_vn_yCM_40to60_pr->SetMarkerSize(2);
      h_vn_yCM_00to10_pr->SetLineColor(2);
      h_vn_yCM_10to40_pr->SetLineColor(4);
      h_vn_yCM_40to60_pr->SetLineColor(8);
      h_vn_yCM_00to10_pr->SetLineWidth(3);
      h_vn_yCM_10to40_pr->SetLineWidth(3);
      h_vn_yCM_40to60_pr->SetLineWidth(3);

      h_vn_yCM_00to10_pr_symm->SetMarkerStyle(20);
      h_vn_yCM_10to40_pr_symm->SetMarkerStyle(20);
      h_vn_yCM_40to60_pr_symm->SetMarkerStyle(20);
      h_vn_yCM_00to10_pr_symm->SetMarkerColor(2);
      h_vn_yCM_10to40_pr_symm->SetMarkerColor(4);
      h_vn_yCM_40to60_pr_symm->SetMarkerColor(8);
      h_vn_yCM_00to10_pr_symm->SetMarkerSize(2);
      h_vn_yCM_10to40_pr_symm->SetMarkerSize(2);
      h_vn_yCM_40to60_pr_symm->SetMarkerSize(2);
      h_vn_yCM_00to10_pr_symm->SetLineColor(2);
      h_vn_yCM_10to40_pr_symm->SetLineColor(4);
      h_vn_yCM_40to60_pr_symm->SetLineColor(8);
      h_vn_yCM_00to10_pr_symm->SetLineWidth(3);
      h_vn_yCM_10to40_pr_symm->SetLineWidth(3);
      h_vn_yCM_40to60_pr_symm->SetLineWidth(3);


      //=== vn vs pT
      h_vn_pT_00to10_pp->SetMarkerStyle(20);
      h_vn_pT_10to40_pp->SetMarkerStyle(20);
      h_vn_pT_40to60_pp->SetMarkerStyle(20);
      h_vn_pT_00to10_pp->SetMarkerColor(2);
      h_vn_pT_10to40_pp->SetMarkerColor(4);
      h_vn_pT_40to60_pp->SetMarkerColor(8);
      h_vn_pT_00to10_pp->SetMarkerSize(2);
      h_vn_pT_10to40_pp->SetMarkerSize(2);
      h_vn_pT_40to60_pp->SetMarkerSize(2);
      h_vn_pT_00to10_pp->SetLineColor(2);
      h_vn_pT_10to40_pp->SetLineColor(4);
      h_vn_pT_40to60_pp->SetLineColor(8);
      h_vn_pT_00to10_pp->SetLineWidth(3);
      h_vn_pT_10to40_pp->SetLineWidth(3);
      h_vn_pT_40to60_pp->SetLineWidth(3);
  
      h_vn_pT_00to10_pm->SetMarkerStyle(20);
      h_vn_pT_10to40_pm->SetMarkerStyle(20);
      h_vn_pT_40to60_pm->SetMarkerStyle(20);
      h_vn_pT_00to10_pm->SetMarkerColor(2);
      h_vn_pT_10to40_pm->SetMarkerColor(4);
      h_vn_pT_40to60_pm->SetMarkerColor(8);
      h_vn_pT_00to10_pm->SetMarkerSize(2);
      h_vn_pT_10to40_pm->SetMarkerSize(2);
      h_vn_pT_40to60_pm->SetMarkerSize(2);
      h_vn_pT_00to10_pm->SetLineColor(2);
      h_vn_pT_10to40_pm->SetLineColor(4);
      h_vn_pT_40to60_pm->SetLineColor(8);
      h_vn_pT_00to10_pm->SetLineWidth(3);
      h_vn_pT_10to40_pm->SetLineWidth(3);
      h_vn_pT_40to60_pm->SetLineWidth(3);

      h_vn_pT_00to10_kp->SetMarkerStyle(20);
      h_vn_pT_10to40_kp->SetMarkerStyle(20);
      h_vn_pT_40to60_kp->SetMarkerStyle(20);
      h_vn_pT_00to10_kp->SetMarkerColor(2);
      h_vn_pT_10to40_kp->SetMarkerColor(4);
      h_vn_pT_40to60_kp->SetMarkerColor(8);
      h_vn_pT_00to10_kp->SetMarkerSize(2);
      h_vn_pT_10to40_kp->SetMarkerSize(2);
      h_vn_pT_40to60_kp->SetMarkerSize(2);
      h_vn_pT_00to10_kp->SetLineColor(2);
      h_vn_pT_10to40_kp->SetLineColor(4);
      h_vn_pT_40to60_kp->SetLineColor(8);
      h_vn_pT_00to10_kp->SetLineWidth(3);
      h_vn_pT_10to40_kp->SetLineWidth(3);
      h_vn_pT_40to60_kp->SetLineWidth(3);
  
      h_vn_pT_00to10_km->SetMarkerStyle(20);
      h_vn_pT_10to40_km->SetMarkerStyle(20);
      h_vn_pT_40to60_km->SetMarkerStyle(20);
      h_vn_pT_00to10_km->SetMarkerColor(2);
      h_vn_pT_10to40_km->SetMarkerColor(4);
      h_vn_pT_40to60_km->SetMarkerColor(8);
      h_vn_pT_00to10_km->SetMarkerSize(2);
      h_vn_pT_10to40_km->SetMarkerSize(2);
      h_vn_pT_40to60_km->SetMarkerSize(2);
      h_vn_pT_00to10_km->SetLineColor(2);
      h_vn_pT_10to40_km->SetLineColor(4);
      h_vn_pT_40to60_km->SetLineColor(8);
      h_vn_pT_00to10_km->SetLineWidth(3);
      h_vn_pT_10to40_km->SetLineWidth(3);
      h_vn_pT_40to60_km->SetLineWidth(3);
  
      h_vn_pT_00to10_pr->SetMarkerStyle(20);
      h_vn_pT_10to40_pr->SetMarkerStyle(20);
      h_vn_pT_40to60_pr->SetMarkerStyle(20);
      h_vn_pT_00to10_pr->SetMarkerColor(2);
      h_vn_pT_10to40_pr->SetMarkerColor(4);
      h_vn_pT_40to60_pr->SetMarkerColor(8);
      h_vn_pT_00to10_pr->SetMarkerSize(2);
      h_vn_pT_10to40_pr->SetMarkerSize(2);
      h_vn_pT_40to60_pr->SetMarkerSize(2);
      h_vn_pT_00to10_pr->SetLineColor(2);
      h_vn_pT_10to40_pr->SetLineColor(4);
      h_vn_pT_40to60_pr->SetLineColor(8);
      h_vn_pT_00to10_pr->SetLineWidth(3);
      h_vn_pT_10to40_pr->SetLineWidth(3);
      h_vn_pT_40to60_pr->SetLineWidth(3);
    }
}// End fixAttributes()
