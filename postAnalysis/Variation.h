#ifndef VARIATION_H
#define VARIATION_H

#include <vector>
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"

class Variation
{
 public:
  Variation(TString prefix, TString order_n_str);
  Variation(TString prefix, TString order_n_str, Variation *normalData);
  ~Variation();
  TString ID;

  struct DataPoint
  {
    Double_t variedValue;
    Double_t variedError;
    Double_t normalValue;
    Double_t normalError;
    Double_t delta;
    Double_t deltaError;
    Double_t nSigma;
    Double_t stdDev;
    Double_t variance;
  };


  TH1D *h_vn_pp;
  TH1D *h_vn_pm;
  TH1D *h_vn_kp;
  TH1D *h_vn_km;
  TH1D *h_vn_pr;
  
  TH1D *h_vn_pp_ext;
  TH1D *h_vn_pm_ext;
  TH1D *h_vn_kp_ext;
  TH1D *h_vn_km_ext;
  TH1D *h_vn_pr_ext;
  
  TH1D *h_vn_pr_for;
  
  TH1D *h_vn_yCM_00to10_pp;
  TH1D *h_vn_yCM_10to40_pp;
  TH1D *h_vn_yCM_40to60_pp;
  TH1D *h_vn_yCM_00to10_pm;
  TH1D *h_vn_yCM_10to40_pm;
  TH1D *h_vn_yCM_40to60_pm;
  TH1D *h_vn_yCM_00to10_kp;
  TH1D *h_vn_yCM_10to40_kp;
  TH1D *h_vn_yCM_40to60_kp;
  TH1D *h_vn_yCM_00to10_km;
  TH1D *h_vn_yCM_10to40_km;
  TH1D *h_vn_yCM_40to60_km;
  TH1D *h_vn_yCM_00to10_pr;
  TH1D *h_vn_yCM_10to40_pr;
  TH1D *h_vn_yCM_40to60_pr;
  TH1D *h_vn_yCM_00to10_pr_symm;
  TH1D *h_vn_yCM_10to40_pr_symm;
  TH1D *h_vn_yCM_40to60_pr_symm;

  TH1D *h_vn_pT_00to10_pp;
  TH1D *h_vn_pT_10to40_pp;
  TH1D *h_vn_pT_40to60_pp;
  TH1D *h_vn_pT_00to10_pm;
  TH1D *h_vn_pT_10to40_pm;
  TH1D *h_vn_pT_40to60_pm;
  TH1D *h_vn_pT_00to10_kp;
  TH1D *h_vn_pT_10to40_kp;
  TH1D *h_vn_pT_40to60_kp;
  TH1D *h_vn_pT_00to10_km;
  TH1D *h_vn_pT_10to40_km;
  TH1D *h_vn_pT_40to60_km;
  TH1D *h_vn_pT_00to10_pr;
  TH1D *h_vn_pT_10to40_pr;
  TH1D *h_vn_pT_40to60_pr;


  std::vector<DataPoint> v_vn_pp;
  std::vector<DataPoint> v_vn_pm;
  std::vector<DataPoint> v_vn_kp;
  std::vector<DataPoint> v_vn_km;
  std::vector<DataPoint> v_vn_pr;
  
  std::vector<DataPoint> v_vn_pp_ext;
  std::vector<DataPoint> v_vn_pm_ext;
  std::vector<DataPoint> v_vn_kp_ext;
  std::vector<DataPoint> v_vn_km_ext;
  std::vector<DataPoint> v_vn_pr_ext;
  
  std::vector<DataPoint> v_vn_pr_for;
  
  std::vector<DataPoint> v_vn_yCM_00to10_pp;
  std::vector<DataPoint> v_vn_yCM_10to40_pp;
  std::vector<DataPoint> v_vn_yCM_40to60_pp;
  std::vector<DataPoint> v_vn_yCM_00to10_pm;
  std::vector<DataPoint> v_vn_yCM_10to40_pm;
  std::vector<DataPoint> v_vn_yCM_40to60_pm;
  std::vector<DataPoint> v_vn_yCM_00to10_kp;
  std::vector<DataPoint> v_vn_yCM_10to40_kp;
  std::vector<DataPoint> v_vn_yCM_40to60_kp;
  std::vector<DataPoint> v_vn_yCM_00to10_km;
  std::vector<DataPoint> v_vn_yCM_10to40_km;
  std::vector<DataPoint> v_vn_yCM_40to60_km;
  std::vector<DataPoint> v_vn_yCM_00to10_pr;
  std::vector<DataPoint> v_vn_yCM_10to40_pr;
  std::vector<DataPoint> v_vn_yCM_40to60_pr;
  std::vector<DataPoint> v_vn_yCM_00to10_pr_symm;
  std::vector<DataPoint> v_vn_yCM_10to40_pr_symm;
  std::vector<DataPoint> v_vn_yCM_40to60_pr_symm;

  std::vector<DataPoint> v_vn_pT_00to10_pp;
  std::vector<DataPoint> v_vn_pT_10to40_pp;
  std::vector<DataPoint> v_vn_pT_40to60_pp;
  std::vector<DataPoint> v_vn_pT_00to10_pm;
  std::vector<DataPoint> v_vn_pT_10to40_pm;
  std::vector<DataPoint> v_vn_pT_40to60_pm;
  std::vector<DataPoint> v_vn_pT_00to10_kp;
  std::vector<DataPoint> v_vn_pT_10to40_kp;
  std::vector<DataPoint> v_vn_pT_40to60_kp;
  std::vector<DataPoint> v_vn_pT_00to10_km;
  std::vector<DataPoint> v_vn_pT_10to40_km;
  std::vector<DataPoint> v_vn_pT_40to60_km;
  std::vector<DataPoint> v_vn_pT_00to10_pr;
  std::vector<DataPoint> v_vn_pT_10to40_pr;
  std::vector<DataPoint> v_vn_pT_40to60_pr;

  
 private:
  TString fileName;
  TFile *file;
  void initialize(TString order_n_str);
  void combine(Variation *normalData);
  void stdDevs();
  void fixAttributes(TString order_n_str);
};

#endif
