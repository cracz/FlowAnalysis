#include "Variation.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "THStack.h"
#include <vector>
#include <iostream>

struct Variances
{
  std::vector<Double_t> v_vn_pp;
  std::vector<Double_t> v_vn_pm;
  std::vector<Double_t> v_vn_kp;
  std::vector<Double_t> v_vn_km;
  std::vector<Double_t> v_vn_pr;
  std::vector<Double_t> v_vn_pp_ext;
  std::vector<Double_t> v_vn_pm_ext;
  std::vector<Double_t> v_vn_kp_ext;
  std::vector<Double_t> v_vn_km_ext;
  std::vector<Double_t> v_vn_pr_ext;
  
  std::vector<Double_t> v_vn_pr_for;

  std::vector<Double_t> v_vn_yCM_00to10_pp;
  std::vector<Double_t> v_vn_yCM_10to40_pp;
  std::vector<Double_t> v_vn_yCM_40to60_pp;
  std::vector<Double_t> v_vn_yCM_00to10_pm;
  std::vector<Double_t> v_vn_yCM_10to40_pm;
  std::vector<Double_t> v_vn_yCM_40to60_pm;
  std::vector<Double_t> v_vn_yCM_00to10_kp;
  std::vector<Double_t> v_vn_yCM_10to40_kp;
  std::vector<Double_t> v_vn_yCM_40to60_kp;
  std::vector<Double_t> v_vn_yCM_00to10_km;
  std::vector<Double_t> v_vn_yCM_10to40_km;
  std::vector<Double_t> v_vn_yCM_40to60_km;
  std::vector<Double_t> v_vn_yCM_00to10_pr;
  std::vector<Double_t> v_vn_yCM_10to40_pr;
  std::vector<Double_t> v_vn_yCM_40to60_pr;
  std::vector<Double_t> v_vn_yCM_00to10_pr_symm;
  std::vector<Double_t> v_vn_yCM_10to40_pr_symm;
  std::vector<Double_t> v_vn_yCM_40to60_pr_symm;

  std::vector<Double_t> v_vn_pT_00to10_pp;
  std::vector<Double_t> v_vn_pT_10to40_pp;
  std::vector<Double_t> v_vn_pT_40to60_pp;
  std::vector<Double_t> v_vn_pT_00to10_pm;
  std::vector<Double_t> v_vn_pT_10to40_pm;
  std::vector<Double_t> v_vn_pT_40to60_pm;
  std::vector<Double_t> v_vn_pT_00to10_kp;
  std::vector<Double_t> v_vn_pT_10to40_kp;
  std::vector<Double_t> v_vn_pT_40to60_kp;
  std::vector<Double_t> v_vn_pT_00to10_km;
  std::vector<Double_t> v_vn_pT_10to40_km;
  std::vector<Double_t> v_vn_pT_40to60_km;
  std::vector<Double_t> v_vn_pT_00to10_pr;
  std::vector<Double_t> v_vn_pT_10to40_pr;
  std::vector<Double_t> v_vn_pT_40to60_pr;


  std::vector<TString> v_vn_pp_st;
  std::vector<TString> v_vn_pm_st;
  std::vector<TString> v_vn_kp_st;
  std::vector<TString> v_vn_km_st;
  std::vector<TString> v_vn_pr_st;
  std::vector<TString> v_vn_pp_ext_st;
  std::vector<TString> v_vn_pm_ext_st;
  std::vector<TString> v_vn_kp_ext_st;
  std::vector<TString> v_vn_km_ext_st;
  std::vector<TString> v_vn_pr_ext_st;
  
  std::vector<TString> v_vn_pr_for_st;

  std::vector<TString> v_vn_yCM_00to10_pp_st;
  std::vector<TString> v_vn_yCM_10to40_pp_st;
  std::vector<TString> v_vn_yCM_40to60_pp_st;
  std::vector<TString> v_vn_yCM_00to10_pm_st;
  std::vector<TString> v_vn_yCM_10to40_pm_st;
  std::vector<TString> v_vn_yCM_40to60_pm_st;
  std::vector<TString> v_vn_yCM_00to10_kp_st;
  std::vector<TString> v_vn_yCM_10to40_kp_st;
  std::vector<TString> v_vn_yCM_40to60_kp_st;
  std::vector<TString> v_vn_yCM_00to10_km_st;
  std::vector<TString> v_vn_yCM_10to40_km_st;
  std::vector<TString> v_vn_yCM_40to60_km_st;
  std::vector<TString> v_vn_yCM_00to10_pr_st;
  std::vector<TString> v_vn_yCM_10to40_pr_st;
  std::vector<TString> v_vn_yCM_40to60_pr_st;
  std::vector<TString> v_vn_yCM_00to10_pr_symm_st;
  std::vector<TString> v_vn_yCM_10to40_pr_symm_st;
  std::vector<TString> v_vn_yCM_40to60_pr_symm_st;

  std::vector<TString> v_vn_pT_00to10_pp_st;
  std::vector<TString> v_vn_pT_10to40_pp_st;
  std::vector<TString> v_vn_pT_40to60_pp_st;
  std::vector<TString> v_vn_pT_00to10_pm_st;
  std::vector<TString> v_vn_pT_10to40_pm_st;
  std::vector<TString> v_vn_pT_40to60_pm_st;
  std::vector<TString> v_vn_pT_00to10_kp_st;
  std::vector<TString> v_vn_pT_10to40_kp_st;
  std::vector<TString> v_vn_pT_40to60_kp_st;
  std::vector<TString> v_vn_pT_00to10_km_st;
  std::vector<TString> v_vn_pT_10to40_km_st;
  std::vector<TString> v_vn_pT_40to60_km_st;
  std::vector<TString> v_vn_pT_00to10_pr_st;
  std::vector<TString> v_vn_pT_10to40_pr_st;
  std::vector<TString> v_vn_pT_40to60_pr_st;
};

void printPointErrors(TGraphErrors *graph);
void outputMaxVariance(std::vector<Variation*> variations);




void finalWithSystematics(TString order_n_str)
{
  //TFile *newFile = new TFile("resultsOfSystematics.root", "RECREATE");
  
  Variation *Normal = new Variation("Normal", order_n_str);

  // Variations are not placed into one array of Variations because some plots may need to
  //omit some of the Variations that don't affect them.
  Variation *tpcEff = new Variation("tpcEff", order_n_str, Normal);
  Variation *Rvtx_low = new Variation("Rvtx_low", order_n_str, Normal);
  Variation *Rvtx_high = new Variation("Rvtx_high", order_n_str, Normal);
  Variation *Zvtx_low = new Variation("Zvtx_low", order_n_str, Normal);
  Variation *Zvtx_high = new Variation("Zvtx_high", order_n_str, Normal);
  Variation *dca_low = new Variation("dca_low", order_n_str, Normal);
  Variation *dca_high = new Variation("dca_high", order_n_str, Normal);
  //Variation *nHits_low = new Variation("nHits_low", order_n_str, Normal);
  Variation *nHits_high = new Variation("nHits_high", order_n_str, Normal);
  Variation *nHitsdEdx_low = new Variation("nHitsdEdx_low", order_n_str, Normal);
  Variation *nHitsdEdx_high = new Variation("nHitsdEdx_high", order_n_str, Normal);
  Variation *tracking_low = new Variation("tracking_low", order_n_str, Normal);
  Variation *tracking_high = new Variation("tracking_high", order_n_str, Normal);
  Variation *nSigPi_low = new Variation("nSigPi_low", order_n_str, Normal);
  Variation *nSigPi_high = new Variation("nSigPi_high", order_n_str, Normal);
  Variation *nSigKa_low = new Variation("nSigKa_low", order_n_str, Normal);
  Variation *nSigKa_high = new Variation("nSigKa_high", order_n_str, Normal);
  Variation *nSigPr_low = new Variation("nSigPr_low", order_n_str, Normal);
  Variation *nSigPr_high = new Variation("nSigPr_high", order_n_str, Normal);
  Variation *m2Pi_low = new Variation("m2Pi_low", order_n_str, Normal);
  Variation *m2Pi_high = new Variation("m2Pi_high", order_n_str, Normal);
  Variation *m2Ka_low = new Variation("m2Ka_low", order_n_str, Normal);
  Variation *m2Ka_high = new Variation("m2Ka_high", order_n_str, Normal);
  Variation *epdC = new Variation("epdC", order_n_str, Normal);
  Variation *yPidPi_low = new Variation("yPidPi_low", order_n_str, Normal);
  Variation *yPidPi_high = new Variation("yPidPi_high", order_n_str, Normal);
  Variation *yPidKa_low = new Variation("yPidKa_low", order_n_str, Normal);
  Variation *yPidKa_high = new Variation("yPidKa_high", order_n_str, Normal);
  Variation *yPidPr_low = new Variation("yPidPr_low", order_n_str, Normal);
  Variation *yPidPr_high = new Variation("yPidPr_high", order_n_str, Normal);
  Variation *ptPidPi_low = new Variation("ptPidPi_low", order_n_str, Normal);
  Variation *ptPidPi_high = new Variation("ptPidPi_high", order_n_str, Normal);
  Variation *ptPidKa_low = new Variation("ptPidKa_low", order_n_str, Normal);
  Variation *ptPidKa_high = new Variation("ptPidKa_high", order_n_str, Normal);
  Variation *ptPidPr_low = new Variation("ptPidPr_low", order_n_str, Normal);
  Variation *ptPidPr_high = new Variation("ptPidPr_high", order_n_str, Normal);
  Variation *yFlow_low = new Variation("yFlow_low", order_n_str, Normal);
  Variation *yFlow_high = new Variation("yFlow_high", order_n_str, Normal);
  Variation *ptFlow_low = new Variation("ptFlow_low", order_n_str, Normal);
  Variation *ptFlow_high = new Variation("ptFlow_high", order_n_str, Normal);

  // Print out which variation had the highest variance for each bin of each histogram
  std::vector<Variation*> variations = {tpcEff,epdC};

  //outputMaxVariance(variations);
  

  //======= CALCULATION OF SYSTEMATIC ERRORS
  Double_t totalVariance = 0;
  Double_t totalTpcEff = 0;
  Double_t totalEpdC = 0;
  Double_t ithBinSysErr;
  Int_t bins;

  //=== pi+ vs centrality
  std::vector<Double_t> v_sys_pp;
  bins = Normal->h_vn_pp->GetNbinsX();
  for (int i = 0; i < bins; i++)      // loop over bins of h_vn_pp, same as elements of v_vn_pp but Normal doesn't have v_vn_pp
    {
      ithBinSysErr = TMath::Sqrt(epdC->v_vn_pp.at(i).variance +
				 nHits_high->v_vn_pp.at(i).variance + nHitsdEdx_low->v_vn_pp.at(i).variance +
				 nHitsdEdx_high->v_vn_pp.at(i).variance + tracking_low->v_vn_pp.at(i).variance +
				 nSigPi_low->v_vn_pp.at(i).variance + nSigPi_high->v_vn_pp.at(i).variance +
				 m2Pi_low->v_vn_pp.at(i).variance + m2Pi_high->v_vn_pp.at(i).variance +
				 ptPidPr_low->v_vn_pp.at(i).variance + ptPidPr_high->v_vn_pp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pp.at(i).variance;
      totalEpdC   += epdC->v_vn_pp.at(i).variance;
      
      v_sys_pp.push_back(ithBinSysErr);
    }// End h_vn_pp loop
  ithBinSysErr = 0;
  //===

  //=== pi- vs centrality
  std::vector<Double_t> v_sys_pm;
  bins = Normal->h_vn_pm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pm.at(i).variance + epdC->v_vn_pm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pm.at(i).variance;
      totalEpdC   += epdC->v_vn_pm.at(i).variance;
      
      v_sys_pm.push_back(ithBinSysErr);
    }// End h_vn_pm loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs centrality
  std::vector<Double_t> v_sys_kp;
  bins = Normal->h_vn_kp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_kp.at(i).variance + epdC->v_vn_kp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_kp.at(i).variance;
      totalEpdC   += epdC->v_vn_kp.at(i).variance;
      
      v_sys_kp.push_back(ithBinSysErr);
    }// End h_vn_kp loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs centrality
  std::vector<Double_t> v_sys_km;
  bins = Normal->h_vn_km->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_km.at(i).variance + epdC->v_vn_km.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_km.at(i).variance;
      totalEpdC   += epdC->v_vn_km.at(i).variance;
      
      v_sys_km.push_back(ithBinSysErr);
    }// End h_vn_km loop
  ithBinSysErr = 0;
  //===

  //=== proton vs centrality
  std::vector<Double_t> v_sys_pr;
  bins = Normal->h_vn_pr->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pr.at(i).variance + epdC->v_vn_pr.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pr.at(i).variance;
      totalEpdC   += epdC->v_vn_pr.at(i).variance;
      
      v_sys_pr.push_back(ithBinSysErr);
    }// End h_vn_pr loop
  ithBinSysErr = 0;
  //===


  //=== pi+ vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_pp;
  bins = Normal->h_vn_yCM_00to10_pp->GetNbinsX();
  for (int i = 0; i < bins; i++)      // loop over bins of h_vn_yCM_00to10_pp, same as elements of v_vn_yCM_00to10_pp but Normal doesn't have v_vn_yCM_00to10_pp
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_00to10_pp.at(i).variance + epdC->v_vn_yCM_00to10_pp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_00to10_pp.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_00to10_pp.at(i).variance;
      
      v_sys_yCM_00to10_pp.push_back(ithBinSysErr);
    }// End h_vn_pp loop
  ithBinSysErr = 0;
  //===

  //=== pi+ vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_pp;
  bins = Normal->h_vn_yCM_10to40_pp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_10to40_pp.at(i).variance + epdC->v_vn_yCM_10to40_pp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_10to40_pp.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_10to40_pp.at(i).variance;
      
      v_sys_yCM_10to40_pp.push_back(ithBinSysErr);
    }// End h_vn_pp loop
  ithBinSysErr = 0;
  //===

  //=== pi+ vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_pp;
  bins = Normal->h_vn_yCM_40to60_pp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_40to60_pp.at(i).variance + epdC->v_vn_yCM_40to60_pp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_40to60_pp.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_40to60_pp.at(i).variance;
      
      v_sys_yCM_40to60_pp.push_back(ithBinSysErr);
    }// End h_vn_pp loop
  ithBinSysErr = 0;
  //===


  //=== pi- vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_pm;
  bins = Normal->h_vn_yCM_00to10_pm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_00to10_pm.at(i).variance + epdC->v_vn_yCM_00to10_pm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_00to10_pm.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_00to10_pm.at(i).variance;
      
      v_sys_yCM_00to10_pm.push_back(ithBinSysErr);
    }// End h_vn_pm loop
  ithBinSysErr = 0;
  //===

  //=== pi- vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_pm;
  bins = Normal->h_vn_yCM_10to40_pm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_10to40_pm.at(i).variance + epdC->v_vn_yCM_10to40_pm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_10to40_pm.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_10to40_pm.at(i).variance;
      
      v_sys_yCM_10to40_pm.push_back(ithBinSysErr);
    }// End h_vn_pm loop
  ithBinSysErr = 0;
  //===

  //=== pi- vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_pm;
  bins = Normal->h_vn_yCM_40to60_pm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_40to60_pm.at(i).variance + epdC->v_vn_yCM_40to60_pm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_40to60_pm.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_40to60_pm.at(i).variance;
      
      v_sys_yCM_40to60_pm.push_back(ithBinSysErr);
    }// End h_vn_pm loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_kp;
  bins = Normal->h_vn_yCM_00to10_kp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_00to10_kp.at(i).variance + epdC->v_vn_yCM_00to10_kp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_00to10_kp.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_00to10_kp.at(i).variance;
      
      v_sys_yCM_00to10_kp.push_back(ithBinSysErr);
    }// End h_vn_kp loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_kp;
  bins = Normal->h_vn_yCM_10to40_kp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_10to40_kp.at(i).variance + epdC->v_vn_yCM_10to40_kp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_10to40_kp.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_10to40_kp.at(i).variance;
      
      v_sys_yCM_10to40_kp.push_back(ithBinSysErr);
    }// End h_vn_kp loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_kp;
  bins = Normal->h_vn_yCM_40to60_kp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_40to60_kp.at(i).variance + epdC->v_vn_yCM_40to60_kp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_40to60_kp.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_40to60_kp.at(i).variance;
      
      v_sys_yCM_40to60_kp.push_back(ithBinSysErr);
    }// End h_vn_kp loop
  ithBinSysErr = 0;
  //===

  
  //=== K- vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_km;
  bins = Normal->h_vn_yCM_00to10_km->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_00to10_km.at(i).variance + epdC->v_vn_yCM_00to10_km.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_00to10_km.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_00to10_km.at(i).variance;
      
      v_sys_yCM_00to10_km.push_back(ithBinSysErr);
    }// End h_vn_km loop
  ithBinSysErr = 0;
  //===

  //=== K- vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_km;
  bins = Normal->h_vn_yCM_10to40_km->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_10to40_km.at(i).variance + epdC->v_vn_yCM_10to40_km.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_10to40_km.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_10to40_km.at(i).variance;
      
      v_sys_yCM_10to40_km.push_back(ithBinSysErr);
    }// End h_vn_km loop
  ithBinSysErr = 0;
  //===

  //=== K- vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_km;
  bins = Normal->h_vn_yCM_40to60_km->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_40to60_km.at(i).variance + epdC->v_vn_yCM_40to60_km.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_40to60_km.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_40to60_km.at(i).variance;
      
      v_sys_yCM_40to60_km.push_back(ithBinSysErr);
    }// End h_vn_km loop
  ithBinSysErr = 0;
  //===


  //=== Proton vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_pr;
  bins = Normal->h_vn_yCM_00to10_pr->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_00to10_pr.at(i).variance + epdC->v_vn_yCM_00to10_pr.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_00to10_pr.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_00to10_pr.at(i).variance;
      
      v_sys_yCM_00to10_pr.push_back(ithBinSysErr);
    }// End h_vn_pr loop
  ithBinSysErr = 0;
  //===

  //=== Proton vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_pr;
  bins = Normal->h_vn_yCM_10to40_pr->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_10to40_pr.at(i).variance + epdC->v_vn_yCM_10to40_pr.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_10to40_pr.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_10to40_pr.at(i).variance;
      
      v_sys_yCM_10to40_pr.push_back(ithBinSysErr);
    }// End h_vn_pr loop
  ithBinSysErr = 0;
  //===

  //=== Proton vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_pr;
  bins = Normal->h_vn_yCM_40to60_pr->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_40to60_pr.at(i).variance + epdC->v_vn_yCM_40to60_pr.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_40to60_pr.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_40to60_pr.at(i).variance;
      
      v_sys_yCM_40to60_pr.push_back(ithBinSysErr);
    }// End h_vn_pr loop
  ithBinSysErr = 0;
  //===


  //=== Proton vs rapidity 0 - 10% symmetric
  std::vector<Double_t> v_sys_yCM_00to10_pr_symm;
  bins = Normal->h_vn_yCM_00to10_pr_symm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_00to10_pr_symm.at(i).variance + epdC->v_vn_yCM_00to10_pr_symm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_00to10_pr_symm.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_00to10_pr_symm.at(i).variance;
      
      v_sys_yCM_00to10_pr_symm.push_back(ithBinSysErr);
    }// End h_vn_pr_symm loop
  ithBinSysErr = 0;
  //===

  //=== Proton vs rapidity 10 - 40% symmetric
  std::vector<Double_t> v_sys_yCM_10to40_pr_symm;
  bins = Normal->h_vn_yCM_10to40_pr_symm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_10to40_pr_symm.at(i).variance + epdC->v_vn_yCM_10to40_pr_symm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_10to40_pr_symm.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_10to40_pr_symm.at(i).variance;
      
      v_sys_yCM_10to40_pr_symm.push_back(ithBinSysErr);
    }// End h_vn_pr_symm loop
  ithBinSysErr = 0;
  //===

  //=== Proton vs rapidity 40 - 60% symmetric
  std::vector<Double_t> v_sys_yCM_40to60_pr_symm;
  bins = Normal->h_vn_yCM_40to60_pr_symm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_yCM_40to60_pr_symm.at(i).variance + epdC->v_vn_yCM_40to60_pr_symm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_yCM_40to60_pr_symm.at(i).variance;
      totalEpdC   += epdC->v_vn_yCM_40to60_pr_symm.at(i).variance;
      
      v_sys_yCM_40to60_pr_symm.push_back(ithBinSysErr);
    }// End h_vn_pr_symm loop
  ithBinSysErr = 0;
  //===




  //=== pi+ vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_pp;
  bins = Normal->h_vn_pT_00to10_pp->GetNbinsX();
  for (int i = 0; i < bins; i++)      // loop over bins of h_vn_pT_00to10_pp, same as elements of v_vn_pT_00to10_pp but Normal doesn't have v_vn_pT_00to10_pp
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_00to10_pp.at(i).variance + epdC->v_vn_pT_00to10_pp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_00to10_pp.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_00to10_pp.at(i).variance;
      
      v_sys_pT_00to10_pp.push_back(ithBinSysErr);
    }// End h_vn_pp loop
  ithBinSysErr = 0;
  //===

  //=== pi+ vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_pp;
  bins = Normal->h_vn_pT_10to40_pp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_10to40_pp.at(i).variance + epdC->v_vn_pT_10to40_pp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_10to40_pp.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_10to40_pp.at(i).variance;
      
      v_sys_pT_10to40_pp.push_back(ithBinSysErr);
    }// End h_vn_pp loop
  ithBinSysErr = 0;
  //===

  //=== pi+ vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_pp;
  bins = Normal->h_vn_pT_40to60_pp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_40to60_pp.at(i).variance + epdC->v_vn_pT_40to60_pp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_40to60_pp.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_40to60_pp.at(i).variance;
      
      v_sys_pT_40to60_pp.push_back(ithBinSysErr);
    }// End h_vn_pp loop
  ithBinSysErr = 0;
  //===


  //=== pi- vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_pm;
  bins = Normal->h_vn_pT_00to10_pm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_00to10_pm.at(i).variance + epdC->v_vn_pT_00to10_pm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_00to10_pm.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_00to10_pm.at(i).variance;
      
      v_sys_pT_00to10_pm.push_back(ithBinSysErr);
    }// End h_vn_pm loop
  ithBinSysErr = 0;
  //===

  //=== pi- vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_pm;
  bins = Normal->h_vn_pT_10to40_pm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_10to40_pm.at(i).variance + epdC->v_vn_pT_10to40_pm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_10to40_pm.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_10to40_pm.at(i).variance;
      
      v_sys_pT_10to40_pm.push_back(ithBinSysErr);
    }// End h_vn_pm loop
  ithBinSysErr = 0;
  //===

  //=== pi- vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_pm;
  bins = Normal->h_vn_pT_40to60_pm->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_40to60_pm.at(i).variance + epdC->v_vn_pT_40to60_pm.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_40to60_pm.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_40to60_pm.at(i).variance;
      
      v_sys_pT_40to60_pm.push_back(ithBinSysErr);
    }// End h_vn_pm loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_kp;
  bins = Normal->h_vn_pT_00to10_kp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_00to10_kp.at(i).variance + epdC->v_vn_pT_00to10_kp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_00to10_kp.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_00to10_kp.at(i).variance;
      
      v_sys_pT_00to10_kp.push_back(ithBinSysErr);
    }// End h_vn_kp loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_kp;
  bins = Normal->h_vn_pT_10to40_kp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_10to40_kp.at(i).variance + epdC->v_vn_pT_10to40_kp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_10to40_kp.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_10to40_kp.at(i).variance;
      
      v_sys_pT_10to40_kp.push_back(ithBinSysErr);
    }// End h_vn_kp loop
  ithBinSysErr = 0;
  //===

  //=== K+ vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_kp;
  bins = Normal->h_vn_pT_40to60_kp->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_40to60_kp.at(i).variance + epdC->v_vn_pT_40to60_kp.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_40to60_kp.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_40to60_kp.at(i).variance;
      
      v_sys_pT_40to60_kp.push_back(ithBinSysErr);
    }// End h_vn_kp loop
  ithBinSysErr = 0;
  //===

  
  //=== K- vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_km;
  bins = Normal->h_vn_pT_00to10_km->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_00to10_km.at(i).variance + epdC->v_vn_pT_00to10_km.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_00to10_km.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_00to10_km.at(i).variance;
      
      v_sys_pT_00to10_km.push_back(ithBinSysErr);
    }// End h_vn_km loop
  ithBinSysErr = 0;
  //===

  //=== K- vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_km;
  bins = Normal->h_vn_pT_10to40_km->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_10to40_km.at(i).variance + epdC->v_vn_pT_10to40_km.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_10to40_km.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_10to40_km.at(i).variance;
      
      v_sys_pT_10to40_km.push_back(ithBinSysErr);
    }// End h_vn_km loop
  ithBinSysErr = 0;
  //===

  //=== K- vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_km;
  bins = Normal->h_vn_pT_40to60_km->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_40to60_km.at(i).variance + epdC->v_vn_pT_40to60_km.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_40to60_km.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_40to60_km.at(i).variance;
      
      v_sys_pT_40to60_km.push_back(ithBinSysErr);
    }// End h_vn_km loop
  ithBinSysErr = 0;
  //===


  //=== Proton vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_pr;
  bins = Normal->h_vn_pT_00to10_pr->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_00to10_pr.at(i).variance + epdC->v_vn_pT_00to10_pr.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_00to10_pr.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_00to10_pr.at(i).variance;
      
      v_sys_pT_00to10_pr.push_back(ithBinSysErr);
    }// End h_vn_pr loop
  ithBinSysErr = 0;
  //===

  //=== Proton vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_pr;
  bins = Normal->h_vn_pT_10to40_pr->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_10to40_pr.at(i).variance + epdC->v_vn_pT_10to40_pr.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_10to40_pr.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_10to40_pr.at(i).variance;
      
      v_sys_pT_10to40_pr.push_back(ithBinSysErr);
    }// End h_vn_pr loop
  ithBinSysErr = 0;
  //===

  //=== Proton vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_pr;
  bins = Normal->h_vn_pT_40to60_pr->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(tpcEff->v_vn_pT_40to60_pr.at(i).variance + epdC->v_vn_pT_40to60_pr.at(i).variance);

      totalVariance  += TMath::Power(ithBinSysErr,2);
      totalTpcEff += tpcEff->v_vn_pT_40to60_pr.at(i).variance;
      totalEpdC   += epdC->v_vn_pT_40to60_pr.at(i).variance;
      
      v_sys_pT_40to60_pr.push_back(ithBinSysErr);
    }// End h_vn_pr loop
  ithBinSysErr = 0;
  //===


  //std::cout << "TPC Efficiency: " << (totalTpcEff/totalVariance) * 100 << "%" << std::endl;
  //std::cout << "EPD C: " << (totalEpdC/totalVariance) * 100 << "%" << std::endl;
  

  // PLOTTING
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
  //canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLeftMargin(0.15);
  canvas->cd();

  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(6);

  TLine *zeroLine = new TLine(0, 0, 60, 0);
  zeroLine->SetLineStyle(9);

  TLine *zeroLine_y = new TLine(0, 0, 1, 0);
  zeroLine_y->SetLineStyle(9);

  TLine *zeroLine_y_pr = new TLine(-1, 0, 1, 0);
  zeroLine_y_pr->SetLineStyle(9);

  TLine *zeroLine_pt = new TLine(0, 0, 2, 0);
  zeroLine_pt->SetLineStyle(9);


  if (order_n_str == "3")
    {      
      TLegend *piLegend = new TLegend(0.67, 0.71, 0.805, 0.87);
      piLegend->AddEntry(Normal->h_vn_pp,"#pi^{+}");
      piLegend->AddEntry(Normal->h_vn_pm,"#pi^{-}");
      piLegend->SetFillColorAlpha(0,0);
      piLegend->SetLineColorAlpha(0,0);

      TLegend *kaLegend = new TLegend(0.7, 0.72, 0.8, 0.87);
      kaLegend->AddEntry(Normal->h_vn_kp,"K^{+}");
      kaLegend->AddEntry(Normal->h_vn_km,"K^{-}");
      kaLegend->SetFillColorAlpha(0,0);
      kaLegend->SetLineColorAlpha(0,0);

      /*
      TLegend *ppExtLegend = new TLegend(0.4, 0.62, 0.7, 0.82);
      ppExtLegend->AddEntry(vn_pp,"#pi^{+}, 0 < y_{CM} < 0.5");
      ppExtLegend->AddEntry(vn_pp_ext,"#pi^{+}, 0.5 < y_{CM} < 1.0");
      ppExtLegend->SetFillColorAlpha(0,0);
      ppExtLegend->SetLineColorAlpha(0,0);

      TLegend *pmExtLegend = new TLegend(0.15, 0.67, 0.45, 0.9);
      pmExtLegend->AddEntry(vn_pm,"#pi^{-}, 0 < y_{CM} < 0.5");
      pmExtLegend->AddEntry(vn_pm_ext,"#pi^{-}, 0.5 < y_{CM} < 1.0");
      pmExtLegend->SetFillColorAlpha(0,0);
      pmExtLegend->SetLineColorAlpha(0,0);

      TLegend *kpExtLegend = new TLegend(0.55, 0.7, 0.85, 0.9);
      kpExtLegend->AddEntry(vn_kp,"K^{+}, 0 < y_{CM} < 0.5");
      kpExtLegend->AddEntry(vn_kp_ext,"K^{+}, 0.5 < y_{CM} < 1.0");
      kpExtLegend->SetFillColorAlpha(0,0);
      kpExtLegend->SetLineColorAlpha(0,0);

      TLegend *kmExtLegend = new TLegend(0.28, 0.68, 0.55, 0.85);
      kmExtLegend->AddEntry(vn_km,"K^{-}, 0 < y_{CM} < 0.5");
      kmExtLegend->AddEntry(vn_km_ext,"K^{-}, 0.5 < y_{CM} < 1.0");
      kmExtLegend->SetFillColorAlpha(0,0);
      kmExtLegend->SetLineColorAlpha(0,0);

      TLegend *prExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      prExtLegend->AddEntry(vn_pr_for,"Proton, -0.5 < y_{CM} < 0");
      prExtLegend->AddEntry(vn_pr,"Proton, 0 < y_{CM} < 0.5");
      prExtLegend->AddEntry(vn_pr_ext,"Proton, 0.5 < y_{CM} < 1.0");
      prExtLegend->SetFillColorAlpha(0,0);
      prExtLegend->SetLineColorAlpha(0,0);

      TLegend *etaLegend = new TLegend(0.65, 0.25, 0.9, 0.45);
      etaLegend->AddEntry(vn_EpdE, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(vn_EpdF, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(vn_TpcB, "TPC -1.0 < #eta < 0");
      */
      
      TLegend *ppLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      ppLegend->AddEntry(Normal->h_vn_yCM_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(Normal->h_vn_yCM_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(Normal->h_vn_yCM_40to60_pp, "40 - 60%");
      ppLegend->SetBorderSize(0);
      ppLegend->SetFillColorAlpha(0,0);

      TLegend *pmLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      pmLegend->AddEntry(Normal->h_vn_yCM_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(Normal->h_vn_yCM_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(Normal->h_vn_yCM_40to60_pm, "40 - 60%");
      pmLegend->SetBorderSize(0);
      pmLegend->SetFillColorAlpha(0,0);

      TLegend *kpLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      kpLegend->AddEntry(Normal->h_vn_yCM_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(Normal->h_vn_yCM_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(Normal->h_vn_yCM_40to60_kp, "40 - 60%");
      kpLegend->SetBorderSize(0);
      kpLegend->SetFillColorAlpha(0,0);

      TLegend *kmLegend = new TLegend(0.18, 0.77, 0.38, 0.87);
      //kmLegend->AddEntry(Normal->h_vn_yCM_00to10_km, "0 - 10%");
      kmLegend->AddEntry(Normal->h_vn_yCM_10to40_km, "10 - 40%");
      //kmLegend->AddEntry(Normal->h_vn_yCM_40to60_km, "40 - 60%");
      kmLegend->SetBorderSize(0);
      kmLegend->SetFillColorAlpha(0,0);

      TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      prLegend->AddEntry(Normal->h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(Normal->h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(Normal->h_vn_yCM_40to60_pr, "40 - 60%");
      prLegend->SetBorderSize(0);
      prLegend->SetFillColorAlpha(0,0);

      

      TLegend *ppPtLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      ppPtLegend->AddEntry(Normal->h_vn_pT_00to10_pp, "0 - 10%");
      ppPtLegend->AddEntry(Normal->h_vn_pT_10to40_pp, "10 - 40%");
      ppPtLegend->AddEntry(Normal->h_vn_pT_40to60_pp, "40 - 60%");
      ppPtLegend->SetBorderSize(0);
      ppPtLegend->SetFillColorAlpha(0,0);

      TLegend *pmPtLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      pmPtLegend->AddEntry(Normal->h_vn_pT_00to10_pm, "0 - 10%");
      pmPtLegend->AddEntry(Normal->h_vn_pT_10to40_pm, "10 - 40%");
      pmPtLegend->AddEntry(Normal->h_vn_pT_40to60_pm, "40 - 60%");
      pmPtLegend->SetBorderSize(0);
      pmPtLegend->SetFillColorAlpha(0,0);

      TLegend *kpPtLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      kpPtLegend->AddEntry(Normal->h_vn_pT_00to10_kp, "0 - 10%");
      kpPtLegend->AddEntry(Normal->h_vn_pT_10to40_kp, "10 - 40%");
      kpPtLegend->AddEntry(Normal->h_vn_pT_40to60_kp, "40 - 60%");
      kpPtLegend->SetBorderSize(0);
      kpPtLegend->SetFillColorAlpha(0,0);

      TLegend *kmPtLegend = new TLegend(0.19, 0.12, 0.39, 0.27);
      //kmPtLegend->AddEntry(Normal->h_vn_pT_00to10_km, "0 - 10%");
      kmPtLegend->AddEntry(Normal->h_vn_pT_10to40_km, "10 - 40%");
      //kmPtLegend->AddEntry(Normal->h_vn_pT_40to60_km, "40 - 60%");
      kmPtLegend->SetBorderSize(0);
      kmPtLegend->SetFillColorAlpha(0,0);

      TLegend *prPtLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      prPtLegend->AddEntry(Normal->h_vn_pT_00to10_pr, "0 - 10%");
      prPtLegend->AddEntry(Normal->h_vn_pT_10to40_pr, "10 - 40%");
      prPtLegend->AddEntry(Normal->h_vn_pT_40to60_pr, "40 - 60%");
      prPtLegend->SetBorderSize(0);
      prPtLegend->SetFillColorAlpha(0,0);


      
      TPaveText *piText = new TPaveText(5, 0.025, 38, 0.07, "NB");
      piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      piText->AddText("0 < y_{CM} < 0.5 GeV");
      piText->AddText("0.18 < p_{T} < 1.6 GeV");
      piText->SetFillColorAlpha(0,0);
      piText->SetLineColorAlpha(0,0);

      TPaveText *kaText = new TPaveText(5, 0.025, 38, 0.07, "NB");
      kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kaText->AddText("0 < y_{CM} < 0.5 GeV");
      kaText->AddText("0.4 < p_{T} < 1.6 GeV");
      kaText->SetFillColorAlpha(0,0);
      kaText->SetLineColorAlpha(0,0);
      
      TPaveText *prText = new TPaveText(5, 0.025, 38, 0.07, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 < p_{T} < 2.0 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);


      TPaveText *ppText = new TPaveText(0.5, 0.07, 0.8, 0.14, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0.18 < p_{T} < 1.6 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);
      ppText->SetTextSize(.04);

      TPaveText *pmText = new TPaveText(0.5, 0.07, 0.8, 0.14, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0.18 < p_{T} < 1.6 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);
      pmText->SetTextSize(.04);
      
      TPaveText *kpText = new TPaveText(0.5, 0.07, 0.8, 0.14, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0.4 < p_{T} < 1.6 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);
      kpText->SetTextSize(.04);
      
      TPaveText *kmText = new TPaveText(0.3, 0.05, 0.7, 0.12, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0.4 < p_{T} < 1.6 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);
      kmText->SetTextSize(.04);
      
      TPaveText *prText_y = new TPaveText(-0.2, 0.02, 0.9, 0.05, "NB");
      prText_y->AddText("Proton");
      prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_y->AddText("0.4 < p_{T} < 2.0 GeV");
      prText_y->SetFillColorAlpha(0,0);
      prText_y->SetLineColorAlpha(0,0);
      prText_y->SetTextSize(.035);
      
      TPaveText *prText_y_symm = new TPaveText(-0.2, 0.02, 0.9, 0.05, "NB");
      prText_y_symm->AddText("Proton");
      prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_y_symm->AddText("1.0 < p_{T} < 2.5 GeV");
      prText_y_symm->SetFillColorAlpha(0,0);
      prText_y_symm->SetLineColorAlpha(0,0);
      prText_y_symm->SetTextSize(.035);
      


      TPaveText *ppPtText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      ppPtText->AddText("#pi^{+}");
      ppPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppPtText->AddText("0 < y_{CM} < 0.5");
      ppPtText->SetFillColorAlpha(0,0);
      ppPtText->SetLineColorAlpha(0,0);
      ppPtText->SetTextSize(.04);

      TPaveText *pmPtText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      pmPtText->AddText("#pi^{-}");
      pmPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmPtText->AddText("0 < y_{CM} < 0.5");
      pmPtText->SetFillColorAlpha(0,0);
      pmPtText->SetLineColorAlpha(0,0);
      pmPtText->SetTextSize(.04);

      TPaveText *kpPtText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      kpPtText->AddText("K^{+}");
      kpPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpPtText->AddText("0 < y_{CM} < 0.5");
      kpPtText->SetFillColorAlpha(0,0);
      kpPtText->SetLineColorAlpha(0,0);
      kpPtText->SetTextSize(.04);

      TPaveText *kmPtText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      kmPtText->AddText("K^{-}");
      kmPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmPtText->AddText("0 < y_{CM} < 0.5");
      kmPtText->SetFillColorAlpha(0,0);
      kmPtText->SetLineColorAlpha(0,0);
      kmPtText->SetTextSize(.04);

      TPaveText *prPtText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      prPtText->AddText("Proton");
      prPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prPtText->AddText("0 < y_{CM} < 0.5");
      prPtText->SetFillColorAlpha(0,0);
      prPtText->SetLineColorAlpha(0,0);
      prPtText->SetTextSize(.04);

      Double_t centralityUpperBound = 0.08;
      Double_t centralityLowerBound = -0.08;
      Double_t rapidityUpperBound = 0.15;
      Double_t rapidityLowerBound = -0.1;
      Double_t rapidityUpperBound_pr = 0.06;
      Double_t rapidityLowerBound_pr = -0.1;
      Double_t ptUpperBound = 0.25;
      Double_t ptLowerBound = -0.25;
      
  
      TGraphErrors *copyWithNewErrors1;
      TGraphErrors *copyWithNewErrors2;
      TGraphErrors *copyWithNewErrors3;
  
      //=== pi+- vs centrality
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pp->Clone());
      for (int i = 0; i < v_sys_pp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pm->Clone());
      for (int i = 0; i < v_sys_pm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pm.at(i)); }

      THStack *piCentralityStack = new THStack("piCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
      piCentralityStack->Add(Normal->h_vn_pp);
      piCentralityStack->Add(Normal->h_vn_pm);

      piCentralityStack->Draw();
      piCentralityStack->SetMinimum(centralityLowerBound);
      piCentralityStack->SetMaximum(centralityUpperBound);
      piCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      piLegend->Draw();
      piText->Draw();
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      canvas->SaveAs("sys_h_vn_pi.png");
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      canvas->Clear();
      //===

      //=== K+- vs centrality
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_kp->Clone());
      for (int i = 0; i < v_sys_kp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_kp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_km->Clone());
      for (int i = 0; i < v_sys_km.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_km.at(i)); }

      THStack *kaCentralityStack = new THStack("kaCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
      kaCentralityStack->Add(Normal->h_vn_kp);
      kaCentralityStack->Add(Normal->h_vn_km);

      kaCentralityStack->Draw();
      kaCentralityStack->SetMinimum(centralityLowerBound);
      kaCentralityStack->SetMaximum(centralityUpperBound);
      kaCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kaLegend->Draw();
      kaText->Draw();
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      canvas->SaveAs("sys_h_vn_ka.png");
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      canvas->Clear();
      //===

      //=== Proton vs centrality
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pr->Clone());
      for (int i = 0; i < v_sys_pr.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pr.at(i)); }

      Normal->h_vn_pr->Draw();
      Normal->h_vn_pr->SetMinimum(centralityLowerBound);
      Normal->h_vn_pr->SetMaximum(centralityUpperBound);
      Normal->h_vn_pr->Draw("E1");
      copyWithNewErrors1->Draw("[]");
      zeroLine->Draw("SAME");
      prText->Draw();
      canvas->SaveAs("sys_h_vn_pr.png");
      delete copyWithNewErrors1;
      canvas->Clear();
      //===


      //=== Pi+ vs rapidity
      THStack *ppRapidityStack = new THStack("ppRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      ppRapidityStack->Add(Normal->h_vn_yCM_00to10_pp);
      ppRapidityStack->Add(Normal->h_vn_yCM_10to40_pp);
      ppRapidityStack->Add(Normal->h_vn_yCM_40to60_pp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pp->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_pp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pp->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pp->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pp.at(i)); }
      
      ppRapidityStack->Draw();
      ppRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      ppRapidityStack->GetXaxis()->SetNdivisions(210);
      ppRapidityStack->SetMaximum(rapidityUpperBound);
      ppRapidityStack->SetMinimum(rapidityLowerBound);
      ppRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      ppLegend->Draw();
      ppText->Draw();
      canvas->SaveAs("sys_ppRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      //=== Pi- vs rapidity
      THStack *pmRapidityStack = new THStack("pmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      pmRapidityStack->Add(Normal->h_vn_yCM_00to10_pm);
      pmRapidityStack->Add(Normal->h_vn_yCM_10to40_pm);
      pmRapidityStack->Add(Normal->h_vn_yCM_40to60_pm);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pm->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_pm.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pm.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pm->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pm.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pm->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pm.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pm.at(i)); }
      
      pmRapidityStack->Draw();
      pmRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      pmRapidityStack->GetXaxis()->SetNdivisions(210);
      pmRapidityStack->SetMaximum(rapidityUpperBound);
      pmRapidityStack->SetMinimum(rapidityLowerBound);
      pmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      pmLegend->Draw();
      pmText->Draw();
      canvas->SaveAs("sys_pmRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===

      //=== K+ vs rapidity
      THStack *kpRapidityStack = new THStack("kpRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      kpRapidityStack->Add(Normal->h_vn_yCM_00to10_kp);
      kpRapidityStack->Add(Normal->h_vn_yCM_10to40_kp);
      kpRapidityStack->Add(Normal->h_vn_yCM_40to60_kp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_kp->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_kp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_kp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_kp->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_kp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_kp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_kp->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_kp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_kp.at(i)); }
      
      kpRapidityStack->Draw();
      kpRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      kpRapidityStack->GetXaxis()->SetNdivisions(210);
      kpRapidityStack->SetMaximum(rapidityUpperBound);
      kpRapidityStack->SetMinimum(rapidityLowerBound);
      kpRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      kpLegend->Draw();
      kpText->Draw();
      canvas->SaveAs("sys_kpRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===

      //=== K- vs rapidity
      THStack *kmRapidityStack = new THStack("kmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      //kmRapidityStack->Add(Normal->h_vn_yCM_00to10_km);
      kmRapidityStack->Add(Normal->h_vn_yCM_10to40_km);
      //kmRapidityStack->Add(Normal->h_vn_yCM_40to60_km);
      /*
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_km->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_km.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_km.at(i)); }
      */
      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_km->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_km.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_km.at(i)); }
      /*
      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_km->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_km.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_km.at(i)); }
      */
      kmRapidityStack->Draw();
      kmRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      kmRapidityStack->GetXaxis()->SetNdivisions(210);
      kmRapidityStack->SetMaximum(rapidityUpperBound);
      kmRapidityStack->SetMinimum(rapidityLowerBound);
      kmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      //copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      //copyWithNewErrors3->Draw("[]");
      kmLegend->Draw();
      kmText->Draw();
      canvas->SaveAs("sys_kmRapidityStack.png");
      canvas->Clear();
      //delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      //delete copyWithNewErrors3;
      //===

      //=== Proton vs rapidity
      THStack *prRapidityStack = new THStack("prRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      prRapidityStack->Add(Normal->h_vn_yCM_00to10_pr);
      prRapidityStack->Add(Normal->h_vn_yCM_10to40_pr);
      prRapidityStack->Add(Normal->h_vn_yCM_40to60_pr);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pr->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_pr.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pr.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pr->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pr.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pr.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pr->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pr.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pr.at(i)); }
      
      prRapidityStack->Draw();
      prRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      prRapidityStack->GetXaxis()->SetNdivisions(210);
      prRapidityStack->SetMaximum(rapidityUpperBound_pr);
      prRapidityStack->SetMinimum(rapidityLowerBound_pr);
      prRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y_pr->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      prLegend->Draw();
      prText_y->Draw();
      canvas->SaveAs("sys_prRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      //=== Proton vs rapidity symmetric across midrapidity
      THStack *prRapidityStack_symm = new THStack("prRapidityStack_symm", ";y-y_{mid};v_{"+order_n_str+"}");
      prRapidityStack_symm->Add(Normal->h_vn_yCM_00to10_pr_symm);
      prRapidityStack_symm->Add(Normal->h_vn_yCM_10to40_pr_symm);
      prRapidityStack_symm->Add(Normal->h_vn_yCM_40to60_pr_symm);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pr_symm->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_pr_symm.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pr_symm.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pr_symm->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pr_symm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pr_symm.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pr_symm->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pr_symm.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pr_symm.at(i)); }
      
      prRapidityStack_symm->Draw();
      prRapidityStack_symm->GetYaxis()->SetTitleOffset(1.9);
      prRapidityStack_symm->GetXaxis()->SetNdivisions(210);
      prRapidityStack_symm->SetMaximum(rapidityUpperBound_pr);
      prRapidityStack_symm->SetMinimum(rapidityLowerBound_pr);
      prRapidityStack_symm->Draw("NOSTACK E1P");
      zeroLine_y_pr->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      prLegend->Draw();
      prText_y_symm->Draw();
      canvas->SaveAs("sys_pr_symmRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===




      //=== Pi+ vs pT
      THStack *ppPtStack = new THStack("ppPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      ppPtStack->Add(Normal->h_vn_pT_00to10_pp);
      ppPtStack->Add(Normal->h_vn_pT_40to60_pp);
      ppPtStack->Add(Normal->h_vn_pT_10to40_pp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_pp->Clone());
      for (int i = 0; i < v_sys_pT_00to10_pp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_pp.at(i));}

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_pp->Clone());
      for (int i = 0; i < v_sys_pT_10to40_pp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_pp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_pp->Clone());
      for (int i = 0; i < v_sys_pT_40to60_pp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_pp.at(i)); }

      printPointErrors(copyWithNewErrors1);
      printPointErrors(copyWithNewErrors2);
      printPointErrors(copyWithNewErrors3);
      
      ppPtStack->Draw();
      ppPtStack->GetYaxis()->SetTitleOffset(1.9);
      ppPtStack->GetXaxis()->SetNdivisions(210);
      ppPtStack->SetMaximum(ptUpperBound);
      ppPtStack->SetMinimum(ptLowerBound);
      ppPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      ppPtLegend->Draw();
      ppPtText->Draw();
      canvas->SaveAs("sys_ppPtStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===



      //=== Pi- vs pT
      THStack *pmPtStack = new THStack("pmPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      pmPtStack->Add(Normal->h_vn_pT_00to10_pm);
      pmPtStack->Add(Normal->h_vn_pT_40to60_pm);
      pmPtStack->Add(Normal->h_vn_pT_10to40_pm);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_pm->Clone());
      for (int i = 0; i < v_sys_pT_00to10_pm.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_pm.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_pm->Clone());
      for (int i = 0; i < v_sys_pT_10to40_pm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_pm.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_pm->Clone());
      for (int i = 0; i < v_sys_pT_40to60_pm.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_pm.at(i)); }


      printPointErrors(copyWithNewErrors1);
      printPointErrors(copyWithNewErrors2);
      printPointErrors(copyWithNewErrors3);

      pmPtStack->Draw();
      pmPtStack->GetYaxis()->SetTitleOffset(1.9);
      pmPtStack->GetXaxis()->SetNdivisions(210);
      pmPtStack->SetMaximum(ptUpperBound);
      pmPtStack->SetMinimum(ptLowerBound);
      pmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      pmPtLegend->Draw();
      pmPtText->Draw();
      canvas->SaveAs("sys_pmPtStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      //=== K+ vs pT
      THStack *kpPtStack = new THStack("kpPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      kpPtStack->Add(Normal->h_vn_pT_00to10_kp);
      kpPtStack->Add(Normal->h_vn_pT_40to60_kp);
      kpPtStack->Add(Normal->h_vn_pT_10to40_kp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_kp->Clone());
      for (int i = 0; i < v_sys_pT_00to10_kp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_kp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_kp->Clone());
      for (int i = 0; i < v_sys_pT_10to40_kp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_kp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_kp->Clone());
      for (int i = 0; i < v_sys_pT_40to60_kp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_kp.at(i)); }

      printPointErrors(copyWithNewErrors1);
      printPointErrors(copyWithNewErrors2);
      printPointErrors(copyWithNewErrors3);

      kpPtStack->Draw();
      kpPtStack->GetYaxis()->SetTitleOffset(1.9);
      kpPtStack->GetXaxis()->SetNdivisions(210);
      kpPtStack->SetMaximum(ptUpperBound);
      kpPtStack->SetMinimum(ptLowerBound);
      kpPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      kpPtLegend->Draw();
      kpPtText->Draw();
      canvas->SaveAs("sys_kpPtStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      //=== K- vs pT
      THStack *kmPtStack = new THStack("kmPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      //kmPtStack->Add(Normal->h_vn_pT_00to10_km);
      //kmPtStack->Add(Normal->h_vn_pT_40to60_km);
      kmPtStack->Add(Normal->h_vn_pT_10to40_km);
      /*
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_km->Clone());
      for (int i = 0; i < v_sys_pT_00to10_km.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_km.at(i)); }
      */
      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_km->Clone());
      for (int i = 0; i < v_sys_pT_10to40_km.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_km.at(i)); }
      /*
      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_km->Clone());
      for (int i = 0; i < v_sys_pT_40to60_km.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_km.at(i)); }
      */

      printPointErrors(copyWithNewErrors2);

      kmPtStack->Draw();
      kmPtStack->GetYaxis()->SetTitleOffset(1.9);
      kmPtStack->GetXaxis()->SetNdivisions(210);
      kmPtStack->SetMaximum(ptUpperBound);
      kmPtStack->SetMinimum(ptLowerBound);
      kmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      //copyWithNewErrors1->Draw("[]");
      //copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      kmPtLegend->Draw();
      kmPtText->Draw();
      canvas->SaveAs("sys_kmPtStack.png");
      canvas->Clear();
      //delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      //delete copyWithNewErrors3;
      //===


      //=== Proton vs pT
      THStack *prPtStack = new THStack("prPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      prPtStack->Add(Normal->h_vn_pT_00to10_pr);
      prPtStack->Add(Normal->h_vn_pT_40to60_pr);
      prPtStack->Add(Normal->h_vn_pT_10to40_pr);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_pr->Clone());
      for (int i = 0; i < v_sys_pT_00to10_pr.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_pr.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_pr->Clone());
      for (int i = 0; i < v_sys_pT_10to40_pr.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_pr.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_pr->Clone());
      for (int i = 0; i < v_sys_pT_40to60_pr.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_pr.at(i)); }

      printPointErrors(copyWithNewErrors1);
      printPointErrors(copyWithNewErrors2);
      printPointErrors(copyWithNewErrors3);
      
      prPtStack->Draw();
      prPtStack->GetYaxis()->SetTitleOffset(1.9);
      prPtStack->GetXaxis()->SetNdivisions(210);
      prPtStack->SetMaximum(ptUpperBound);
      prPtStack->SetMinimum(ptLowerBound);
      prPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      prPtLegend->Draw();
      prPtText->Draw();
      canvas->SaveAs("sys_prPtStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      delete canvas;
    }// End if order_n_str == 3



  //newFile->cd(); // It is important to call this right before writing.


  

  delete Normal;
  delete tpcEff;
  delete Rvtx_low;
  delete Rvtx_high;
  delete Zvtx_low;
  delete Zvtx_high;
  delete dca_low;
  delete dca_high;
  delete nHits_high;
  delete nHitsdEdx_low;
  delete nHitsdEdx_high;
  delete tracking_low;
  delete tracking_high;
  delete nSigPi_low;
  delete nSigPi_high;
  delete nSigKa_low;
  delete nSigKa_high;
  delete nSigPr_low;
  delete nSigPr_high;
  delete m2Pi_low;
  delete m2Pi_high;
  delete m2Ka_low;
  delete m2Ka_high;
  delete epdC;
  delete yPidPi_low;
  delete yPidPi_high;
  delete yPidKa_low;
  delete yPidKa_high;
  delete yPidPr_low;
  delete yPidPr_high;
  delete ptPidPi_low;
  delete ptPidPi_high;
  delete ptPidKa_low;
  delete ptPidKa_high;
  delete ptPidPr_low;
  delete ptPidPr_high;
  delete yFlow_low;
  delete yFlow_high;
  delete ptFlow_low;
  delete ptFlow_high;
  
  //newFile->Close();
  //delete newFile;
}



void printPointErrors(TGraphErrors *graph)
{
  Int_t nPoints = graph->GetN();

  std::cout << graph->GetName() << ":" << std::endl;
  
  for (int i = 0; i < nPoints; i++)
    {
      std::cout << "Bin " << i+1 << ": " << graph->GetErrorY(i) << std::endl;
    }
}


// Maybe unnecessary
void outputMaxVariance(std::vector<Variation*> variations)
{
  Variances maxVariances;
  
  // Output the max delta for each histograms' bins
  for (int i = 0; i < variations.size(); i++)
    {
      if (i == 0)
	{
	  for (int j = 0; j < variations.at(i)->h_vn_pp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pp.push_back(variations.at(i)->v_vn_pp.at(j).variance);
	      maxVariances.v_vn_pp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pm.push_back(variations.at(i)->v_vn_pm.at(j).variance);
	      maxVariances.v_vn_pm_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_kp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_kp.push_back(variations.at(i)->v_vn_kp.at(j).variance);
	      maxVariances.v_vn_kp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_km->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_km.push_back(variations.at(i)->v_vn_km.at(j).variance);
	      maxVariances.v_vn_km_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pr->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pr.push_back(variations.at(i)->v_vn_pr.at(j).variance);
	      maxVariances.v_vn_pr_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pp_ext->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pp_ext.push_back(variations.at(i)->v_vn_pp_ext.at(j).variance);
	      maxVariances.v_vn_pp_ext_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pm_ext->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pm_ext.push_back(variations.at(i)->v_vn_pm_ext.at(j).variance);
	      maxVariances.v_vn_pm_ext_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_kp_ext->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_kp_ext.push_back(variations.at(i)->v_vn_kp_ext.at(j).variance);
	      maxVariances.v_vn_kp_ext_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_km_ext->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_km_ext.push_back(variations.at(i)->v_vn_km_ext.at(j).variance);
	      maxVariances.v_vn_km_ext_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pr_ext->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pr_ext.push_back(variations.at(i)->v_vn_pr_ext.at(j).variance);
	      maxVariances.v_vn_pr_ext_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pr_for->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pr_for.push_back(variations.at(i)->v_vn_pr_for.at(j).variance);
	      maxVariances.v_vn_pr_for_st.push_back(variations.at(i)->ID);
	    }

	  // pi+ vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_00to10_pp.push_back(variations.at(i)->v_vn_yCM_00to10_pp.at(j).variance);
	      maxVariances.v_vn_yCM_00to10_pp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_10to40_pp.push_back(variations.at(i)->v_vn_yCM_10to40_pp.at(j).variance);
	      maxVariances.v_vn_yCM_10to40_pp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_40to60_pp.push_back(variations.at(i)->v_vn_yCM_40to60_pp.at(j).variance);
	      maxVariances.v_vn_yCM_40to60_pp_st.push_back(variations.at(i)->ID);
	    }

	  // pi- vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_00to10_pm.push_back(variations.at(i)->v_vn_yCM_00to10_pm.at(j).variance);
	      maxVariances.v_vn_yCM_00to10_pm_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_10to40_pm.push_back(variations.at(i)->v_vn_yCM_10to40_pm.at(j).variance);
	      maxVariances.v_vn_yCM_10to40_pm_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_40to60_pm.push_back(variations.at(i)->v_vn_yCM_40to60_pm.at(j).variance);
	      maxVariances.v_vn_yCM_40to60_pm_st.push_back(variations.at(i)->ID);
	    }

	  // K+ vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_kp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_00to10_kp.push_back(variations.at(i)->v_vn_yCM_00to10_kp.at(j).variance);
	      maxVariances.v_vn_yCM_00to10_kp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_kp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_10to40_kp.push_back(variations.at(i)->v_vn_yCM_10to40_kp.at(j).variance);
	      maxVariances.v_vn_yCM_10to40_kp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_kp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_40to60_kp.push_back(variations.at(i)->v_vn_yCM_40to60_kp.at(j).variance);
	      maxVariances.v_vn_yCM_40to60_kp_st.push_back(variations.at(i)->ID);
	    }

	  // K- vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_km->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_00to10_km.push_back(variations.at(i)->v_vn_yCM_00to10_km.at(j).variance);
	      maxVariances.v_vn_yCM_00to10_km_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_km->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_10to40_km.push_back(variations.at(i)->v_vn_yCM_10to40_km.at(j).variance);
	      maxVariances.v_vn_yCM_10to40_km_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_km->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_40to60_km.push_back(variations.at(i)->v_vn_yCM_40to60_km.at(j).variance);
	      maxVariances.v_vn_yCM_40to60_km_st.push_back(variations.at(i)->ID);
	    }

	  // pr vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pr->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_00to10_pr.push_back(variations.at(i)->v_vn_yCM_00to10_pr.at(j).variance);
	      maxVariances.v_vn_yCM_00to10_pr_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pr->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_10to40_pr.push_back(variations.at(i)->v_vn_yCM_10to40_pr.at(j).variance);
	      maxVariances.v_vn_yCM_10to40_pr_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pr->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_40to60_pr.push_back(variations.at(i)->v_vn_yCM_40to60_pr.at(j).variance);
	      maxVariances.v_vn_yCM_40to60_pr_st.push_back(variations.at(i)->ID);
	    }

	  // pr vs y symmetric
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pr_symm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_00to10_pr_symm.push_back(variations.at(i)->v_vn_yCM_00to10_pr_symm.at(j).variance);
	      maxVariances.v_vn_yCM_00to10_pr_symm_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pr_symm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_10to40_pr_symm.push_back(variations.at(i)->v_vn_yCM_10to40_pr_symm.at(j).variance);
	      maxVariances.v_vn_yCM_10to40_pr_symm_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pr_symm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_yCM_40to60_pr_symm.push_back(variations.at(i)->v_vn_yCM_40to60_pr_symm.at(j).variance);
	      maxVariances.v_vn_yCM_40to60_pr_symm_st.push_back(variations.at(i)->ID);
	    }

	  // pi+ vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_pp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_00to10_pp.push_back(variations.at(i)->v_vn_pT_00to10_pp.at(j).variance);
	      maxVariances.v_vn_pT_00to10_pp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_pp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_10to40_pp.push_back(variations.at(i)->v_vn_pT_10to40_pp.at(j).variance);
	      maxVariances.v_vn_pT_10to40_pp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_pp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_40to60_pp.push_back(variations.at(i)->v_vn_pT_40to60_pp.at(j).variance);
	      maxVariances.v_vn_pT_40to60_pp_st.push_back(variations.at(i)->ID);
	    }

	  // pi- vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_pm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_00to10_pm.push_back(variations.at(i)->v_vn_pT_00to10_pm.at(j).variance);
	      maxVariances.v_vn_pT_00to10_pm_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_pm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_10to40_pm.push_back(variations.at(i)->v_vn_pT_10to40_pm.at(j).variance);
	      maxVariances.v_vn_pT_10to40_pm_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_pm->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_40to60_pm.push_back(variations.at(i)->v_vn_pT_40to60_pm.at(j).variance);
	      maxVariances.v_vn_pT_40to60_pm_st.push_back(variations.at(i)->ID);
	    }

	  // K+ vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_kp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_00to10_kp.push_back(variations.at(i)->v_vn_pT_00to10_kp.at(j).variance);
	      maxVariances.v_vn_pT_00to10_kp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_kp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_10to40_kp.push_back(variations.at(i)->v_vn_pT_10to40_kp.at(j).variance);
	      maxVariances.v_vn_pT_10to40_kp_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_kp->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_40to60_kp.push_back(variations.at(i)->v_vn_pT_40to60_kp.at(j).variance);
	      maxVariances.v_vn_pT_40to60_kp_st.push_back(variations.at(i)->ID);
	    }

	  // K- vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_km->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_00to10_km.push_back(variations.at(i)->v_vn_pT_00to10_km.at(j).variance);
	      maxVariances.v_vn_pT_00to10_km_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_km->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_10to40_km.push_back(variations.at(i)->v_vn_pT_10to40_km.at(j).variance);
	      maxVariances.v_vn_pT_10to40_km_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_km->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_40to60_km.push_back(variations.at(i)->v_vn_pT_40to60_km.at(j).variance);
	      maxVariances.v_vn_pT_40to60_km_st.push_back(variations.at(i)->ID);
	    }

	  // proton vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_pr->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_00to10_pr.push_back(variations.at(i)->v_vn_pT_00to10_pr.at(j).variance);
	      maxVariances.v_vn_pT_00to10_pr_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_pr->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_10to40_pr.push_back(variations.at(i)->v_vn_pT_10to40_pr.at(j).variance);
	      maxVariances.v_vn_pT_10to40_pr_st.push_back(variations.at(i)->ID);
	    }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_pr->GetNbinsX(); j++)
	    {
	      maxVariances.v_vn_pT_40to60_pr.push_back(variations.at(i)->v_vn_pT_40to60_pr.at(j).variance);
	      maxVariances.v_vn_pT_40to60_pr_st.push_back(variations.at(i)->ID);
	    }
	}// End if (i == 0)
      else
	{
	  for (int j = 0; j < variations.at(i)->h_vn_pp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pp.at(j).variance > maxVariances.v_vn_pp.at(j))
	      {
		maxVariances.v_vn_pp.at(j)    = variations.at(i)->v_vn_pp.at(j).variance;
		maxVariances.v_vn_pp_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_pm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pm.at(j).variance > maxVariances.v_vn_pm.at(j)) 
	      {
		maxVariances.v_vn_pm.at(j)    = variations.at(i)->v_vn_pm.at(j).variance;
		maxVariances.v_vn_pm_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_kp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_kp.at(j).variance > maxVariances.v_vn_kp.at(j)) 
	      {
		maxVariances.v_vn_kp.at(j)    = variations.at(i)->v_vn_kp.at(j).variance;
		maxVariances.v_vn_kp_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_km->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_km.at(j).variance > maxVariances.v_vn_km.at(j)) 
	      {
		maxVariances.v_vn_km.at(j)    = variations.at(i)->v_vn_km.at(j).variance;
		maxVariances.v_vn_km_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_pr->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pr.at(j).variance > maxVariances.v_vn_pr.at(j)) 
	      {
		maxVariances.v_vn_pr.at(j)    = variations.at(i)->v_vn_pr.at(j).variance;
		maxVariances.v_vn_pr_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_pp_ext->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pp_ext.at(j).variance > maxVariances.v_vn_pp_ext.at(j)) 
	      {
		maxVariances.v_vn_pp_ext.at(j)    = variations.at(i)->v_vn_pp_ext.at(j).variance;
		maxVariances.v_vn_pp_ext_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_pm_ext->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pm_ext.at(j).variance > maxVariances.v_vn_pm_ext.at(j)) 
	      {
		maxVariances.v_vn_pm_ext.at(j)    = variations.at(i)->v_vn_pm_ext.at(j).variance;
		maxVariances.v_vn_pm_ext_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_kp_ext->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_kp_ext.at(j).variance > maxVariances.v_vn_kp_ext.at(j)) 
	      {
		maxVariances.v_vn_kp_ext.at(j)    = variations.at(i)->v_vn_kp_ext.at(j).variance;
		maxVariances.v_vn_kp_ext_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_km_ext->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_km_ext.at(j).variance > maxVariances.v_vn_km_ext.at(j)) 
	      {
		maxVariances.v_vn_km_ext.at(j)    = variations.at(i)->v_vn_km_ext.at(j).variance;
		maxVariances.v_vn_km_ext_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_pr_ext->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pr_ext.at(j).variance > maxVariances.v_vn_pr_ext.at(j)) 
	      {
		maxVariances.v_vn_pr_ext.at(j)    = variations.at(i)->v_vn_pr_ext.at(j).variance;
		maxVariances.v_vn_pr_ext_st.at(j) = variations.at(i)->ID;
	      }

	  for (int j = 0; j < variations.at(i)->h_vn_pr_for->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pr_for.at(j).variance > maxVariances.v_vn_pr_for.at(j))
	      {
		maxVariances.v_vn_pr_for.at(j)    = variations.at(i)->v_vn_pr_for.at(j).variance;
		maxVariances.v_vn_pr_for_st.at(j) = variations.at(i)->ID;
	      }

 	  // pi+ vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_00to10_pp.at(j).variance > maxVariances.v_vn_yCM_00to10_pp.at(j))
	      {
		maxVariances.v_vn_yCM_00to10_pp.at(j)    = variations.at(i)->v_vn_yCM_00to10_pp.at(j).variance;
		maxVariances.v_vn_yCM_00to10_pp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_10to40_pp.at(j).variance > maxVariances.v_vn_yCM_10to40_pp.at(j))
	      {
		maxVariances.v_vn_yCM_10to40_pp.at(j)    = variations.at(i)->v_vn_yCM_10to40_pp.at(j).variance;
		maxVariances.v_vn_yCM_10to40_pp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_40to60_pp.at(j).variance > maxVariances.v_vn_yCM_40to60_pp.at(j))
	      {
		maxVariances.v_vn_yCM_40to60_pp.at(j)    = variations.at(i)->v_vn_yCM_40to60_pp.at(j).variance;
		maxVariances.v_vn_yCM_40to60_pp_st.at(j) = variations.at(i)->ID;
	      }

	  // pi- vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_00to10_pm.at(j).variance > maxVariances.v_vn_yCM_00to10_pm.at(j))
	      {
		maxVariances.v_vn_yCM_00to10_pm.at(j)    = variations.at(i)->v_vn_yCM_00to10_pm.at(j).variance;
		maxVariances.v_vn_yCM_00to10_pm_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_10to40_pm.at(j).variance > maxVariances.v_vn_yCM_10to40_pm.at(j))
	      {
		maxVariances.v_vn_yCM_10to40_pm.at(j)    = variations.at(i)->v_vn_yCM_10to40_pm.at(j).variance;
		maxVariances.v_vn_yCM_10to40_pm_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_40to60_pm.at(j).variance > maxVariances.v_vn_yCM_40to60_pm.at(j))
	      {
		maxVariances.v_vn_yCM_40to60_pm.at(j)    = variations.at(i)->v_vn_yCM_40to60_pm.at(j).variance;
		maxVariances.v_vn_yCM_40to60_pm_st.at(j) = variations.at(i)->ID;
	      }

	  // K+ vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_kp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_00to10_kp.at(j).variance > maxVariances.v_vn_yCM_00to10_kp.at(j))
	      {
		maxVariances.v_vn_yCM_00to10_kp.at(j)    = variations.at(i)->v_vn_yCM_00to10_kp.at(j).variance;
		maxVariances.v_vn_yCM_00to10_kp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_kp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_10to40_kp.at(j).variance > maxVariances.v_vn_yCM_10to40_kp.at(j))
	      {
		maxVariances.v_vn_yCM_10to40_kp.at(j)    = variations.at(i)->v_vn_yCM_10to40_kp.at(j).variance;
		maxVariances.v_vn_yCM_10to40_kp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_kp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_40to60_kp.at(j).variance > maxVariances.v_vn_yCM_40to60_kp.at(j))
	      {
		maxVariances.v_vn_yCM_40to60_kp.at(j)    = variations.at(i)->v_vn_yCM_40to60_kp.at(j).variance;
		maxVariances.v_vn_yCM_40to60_kp_st.at(j) = variations.at(i)->ID;
	      }

	  // K- vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_km->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_00to10_km.at(j).variance > maxVariances.v_vn_yCM_00to10_km.at(j))
	      {
		maxVariances.v_vn_yCM_00to10_km.at(j)    = variations.at(i)->v_vn_yCM_00to10_km.at(j).variance;
		maxVariances.v_vn_yCM_00to10_km_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_km->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_10to40_km.at(j).variance > maxVariances.v_vn_yCM_10to40_km.at(j))
	      {
		maxVariances.v_vn_yCM_10to40_km.at(j)    = variations.at(i)->v_vn_yCM_10to40_km.at(j).variance;
		maxVariances.v_vn_yCM_10to40_km_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_km->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_40to60_km.at(j).variance > maxVariances.v_vn_yCM_40to60_km.at(j))
	      {
		maxVariances.v_vn_yCM_40to60_km.at(j)    = variations.at(i)->v_vn_yCM_40to60_km.at(j).variance;
		maxVariances.v_vn_yCM_40to60_km_st.at(j) = variations.at(i)->ID;
	      }

	  // proton vs y
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pr->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_00to10_pr.at(j).variance > maxVariances.v_vn_yCM_00to10_pr.at(j))
	      {
		maxVariances.v_vn_yCM_00to10_pr.at(j)    = variations.at(i)->v_vn_yCM_00to10_pr.at(j).variance;
		maxVariances.v_vn_yCM_00to10_pr_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pr->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_10to40_pr.at(j).variance > maxVariances.v_vn_yCM_10to40_pr.at(j))
	      {
		maxVariances.v_vn_yCM_10to40_pr.at(j)    = variations.at(i)->v_vn_yCM_10to40_pr.at(j).variance;
		maxVariances.v_vn_yCM_10to40_pr_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pr->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_40to60_pr.at(j).variance > maxVariances.v_vn_yCM_40to60_pr.at(j))
	      {
		maxVariances.v_vn_yCM_40to60_pr.at(j)    = variations.at(i)->v_vn_yCM_40to60_pr.at(j).variance;
		maxVariances.v_vn_yCM_40to60_pr_st.at(j) = variations.at(i)->ID;
	      }

	  // proton vs y symmmetric
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_00to10_pr_symm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_00to10_pr_symm.at(j).variance > maxVariances.v_vn_yCM_00to10_pr_symm.at(j))
	      {
		maxVariances.v_vn_yCM_00to10_pr_symm.at(j)    = variations.at(i)->v_vn_yCM_00to10_pr_symm.at(j).variance;
		maxVariances.v_vn_yCM_00to10_pr_symm_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_10to40_pr_symm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_10to40_pr_symm.at(j).variance > maxVariances.v_vn_yCM_10to40_pr_symm.at(j))
	      {
		maxVariances.v_vn_yCM_10to40_pr_symm.at(j)    = variations.at(i)->v_vn_yCM_10to40_pr_symm.at(j).variance;
		maxVariances.v_vn_yCM_10to40_pr_symm_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_yCM_40to60_pr_symm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_yCM_40to60_pr_symm.at(j).variance > maxVariances.v_vn_yCM_40to60_pr_symm.at(j))
	      {
		maxVariances.v_vn_yCM_40to60_pr_symm.at(j)    = variations.at(i)->v_vn_yCM_40to60_pr_symm.at(j).variance;
		maxVariances.v_vn_yCM_40to60_pr_symm_st.at(j) = variations.at(i)->ID;
	      }


	  // pi+ vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_pp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_00to10_pp.at(j).variance > maxVariances.v_vn_pT_00to10_pp.at(j))
	      {
		maxVariances.v_vn_pT_00to10_pp.at(j)    = variations.at(i)->v_vn_pT_00to10_pp.at(j).variance;
		maxVariances.v_vn_pT_00to10_pp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_pp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_10to40_pp.at(j).variance > maxVariances.v_vn_pT_10to40_pp.at(j))
	      {
		maxVariances.v_vn_pT_10to40_pp.at(j)    = variations.at(i)->v_vn_pT_10to40_pp.at(j).variance;
		maxVariances.v_vn_pT_10to40_pp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_pp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_40to60_pp.at(j).variance > maxVariances.v_vn_pT_40to60_pp.at(j))
	      {
		maxVariances.v_vn_pT_40to60_pp.at(j)    = variations.at(i)->v_vn_pT_40to60_pp.at(j).variance;
		maxVariances.v_vn_pT_40to60_pp_st.at(j) = variations.at(i)->ID;
	      }

	  // pi- vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_pm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_00to10_pm.at(j).variance > maxVariances.v_vn_pT_00to10_pm.at(j))
	      {
		maxVariances.v_vn_pT_00to10_pm.at(j)    = variations.at(i)->v_vn_pT_00to10_pm.at(j).variance;
		maxVariances.v_vn_pT_00to10_pm_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_pm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_10to40_pm.at(j).variance > maxVariances.v_vn_pT_10to40_pm.at(j))
	      {
		maxVariances.v_vn_pT_10to40_pm.at(j)    = variations.at(i)->v_vn_pT_10to40_pm.at(j).variance;
		maxVariances.v_vn_pT_10to40_pm_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_pm->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_40to60_pm.at(j).variance > maxVariances.v_vn_pT_40to60_pm.at(j))
	      {
		maxVariances.v_vn_pT_40to60_pm.at(j)    = variations.at(i)->v_vn_pT_40to60_pm.at(j).variance;
		maxVariances.v_vn_pT_40to60_pm_st.at(j) = variations.at(i)->ID;
	      }

	  // K+ vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_kp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_00to10_kp.at(j).variance > maxVariances.v_vn_pT_00to10_kp.at(j))
	      {
		maxVariances.v_vn_pT_00to10_kp.at(j)    = variations.at(i)->v_vn_pT_00to10_kp.at(j).variance;
		maxVariances.v_vn_pT_00to10_kp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_kp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_10to40_kp.at(j).variance > maxVariances.v_vn_pT_10to40_kp.at(j))
	      {
		maxVariances.v_vn_pT_10to40_kp.at(j)    = variations.at(i)->v_vn_pT_10to40_kp.at(j).variance;
		maxVariances.v_vn_pT_10to40_kp_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_kp->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_40to60_kp.at(j).variance > maxVariances.v_vn_pT_40to60_kp.at(j))
	      {
		maxVariances.v_vn_pT_40to60_kp.at(j)    = variations.at(i)->v_vn_pT_40to60_kp.at(j).variance;
		maxVariances.v_vn_pT_40to60_kp_st.at(j) = variations.at(i)->ID;
	      }

	  // K- vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_km->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_00to10_km.at(j).variance > maxVariances.v_vn_pT_00to10_km.at(j))
	      {
		maxVariances.v_vn_pT_00to10_km.at(j)    = variations.at(i)->v_vn_pT_00to10_km.at(j).variance;
		maxVariances.v_vn_pT_00to10_km_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_km->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_10to40_km.at(j).variance > maxVariances.v_vn_pT_10to40_km.at(j))
	      {
		maxVariances.v_vn_pT_10to40_km.at(j)    = variations.at(i)->v_vn_pT_10to40_km.at(j).variance;
		maxVariances.v_vn_pT_10to40_km_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_km->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_40to60_km.at(j).variance > maxVariances.v_vn_pT_40to60_km.at(j))
	      {
		maxVariances.v_vn_pT_40to60_km.at(j)    = variations.at(i)->v_vn_pT_40to60_km.at(j).variance;
		maxVariances.v_vn_pT_40to60_km_st.at(j) = variations.at(i)->ID;
	      }

	  // proton vs pT
	  for (int j = 0; j < variations.at(i)->h_vn_pT_00to10_pr->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_00to10_pr.at(j).variance > maxVariances.v_vn_pT_00to10_pr.at(j))
	      {
		maxVariances.v_vn_pT_00to10_pr.at(j)    = variations.at(i)->v_vn_pT_00to10_pr.at(j).variance;
		maxVariances.v_vn_pT_00to10_pr_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_10to40_pr->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_10to40_pr.at(j).variance > maxVariances.v_vn_pT_10to40_pr.at(j))
	      {
		maxVariances.v_vn_pT_10to40_pr.at(j)    = variations.at(i)->v_vn_pT_10to40_pr.at(j).variance;
		maxVariances.v_vn_pT_10to40_pr_st.at(j) = variations.at(i)->ID;
	      }
	  for (int j = 0; j < variations.at(i)->h_vn_pT_40to60_pr->GetNbinsX(); j++)
	    if (variations.at(i)->v_vn_pT_40to60_pr.at(j).variance > maxVariances.v_vn_pT_40to60_pr.at(j))
	      {
		maxVariances.v_vn_pT_40to60_pr.at(j)    = variations.at(i)->v_vn_pT_40to60_pr.at(j).variance;
		maxVariances.v_vn_pT_40to60_pr_st.at(j) = variations.at(i)->ID;
	      }

	}// End else
    } // End for(variations)

  std::cout << "h_vn_pp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pp_st.at(i) << std::endl; }
  std::cout << std::endl;
  
  std::cout << "h_vn_pm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pm_st.at(i) << std::endl; }
  std::cout << std::endl;
  
  std::cout << "h_vn_kp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_kp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_kp_st.at(i) << std::endl; }
  std::cout << std::endl;
  
  std::cout << "h_vn_km:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_km_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_km_st.at(i) << std::endl; }
  std::cout << std::endl;
 
  std::cout << "h_vn_pr:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pr_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pr_st.at(i) << std::endl; }
  std::cout << std::endl;
  

  std::cout << "h_vn_yCM_00to10_pp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_00to10_pp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_00to10_pp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_10to40_pp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_10to40_pp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_10to40_pp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_40to60_pp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_40to60_pp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_40to60_pp_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_yCM_00to10_pm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_00to10_pm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_00to10_pm_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_10to40_pm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_10to40_pm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_10to40_pm_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_40to60_pm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_40to60_pm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_40to60_pm_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_yCM_00to10_kp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_00to10_kp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_00to10_kp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_10to40_kp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_10to40_kp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_10to40_kp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_40to60_kp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_40to60_kp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_40to60_kp_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_yCM_00to10_km:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_00to10_km_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_00to10_km_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_10to40_km:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_10to40_km_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_10to40_km_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_40to60_km:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_40to60_km_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_40to60_km_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_yCM_00to10_pr:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_00to10_pr_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_00to10_pr_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_10to40_pr:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_10to40_pr_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_10to40_pr_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_40to60_pr:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_40to60_pr_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_40to60_pr_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_yCM_00to10_pr_symm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_00to10_pr_symm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_00to10_pr_symm_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_10to40_pr_symm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_10to40_pr_symm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_10to40_pr_symm_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_yCM_40to60_pr_symm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_yCM_40to60_pr_symm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_yCM_40to60_pr_symm_st.at(i) << std::endl; }
  std::cout << std::endl;



    std::cout << "h_vn_pT_00to10_pp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_00to10_pp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_00to10_pp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_10to40_pp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_10to40_pp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_10to40_pp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_40to60_pp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_40to60_pp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_40to60_pp_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_pT_00to10_pm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_00to10_pm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_00to10_pm_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_10to40_pm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_10to40_pm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_10to40_pm_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_40to60_pm:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_40to60_pm_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_40to60_pm_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_pT_00to10_kp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_00to10_kp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_00to10_kp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_10to40_kp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_10to40_kp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_10to40_kp_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_40to60_kp:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_40to60_kp_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_40to60_kp_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_pT_00to10_km:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_00to10_km_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_00to10_km_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_10to40_km:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_10to40_km_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_10to40_km_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_40to60_km:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_40to60_km_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_40to60_km_st.at(i) << std::endl; }
  std::cout << std::endl;


  std::cout << "h_vn_pT_00to10_pr:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_00to10_pr_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_00to10_pr_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_10to40_pr:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_10to40_pr_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_10to40_pr_st.at(i) << std::endl; }
  std::cout << std::endl;

  std::cout << "h_vn_pT_40to60_pr:" << std::endl;
  for (int i = 0; i < maxVariances.v_vn_pT_40to60_pr_st.size(); i++) { std::cout << "Bin " << i+1 << ": " << maxVariances.v_vn_pT_40to60_pr_st.at(i) << std::endl; }
  std::cout << std::endl;

}//End outputMaxVariances()
