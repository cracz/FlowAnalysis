#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <string>
#include "TString.h"

class ConfigReader
{
 public:
  ConfigReader();
  ~ConfigReader();

  Int_t epd_max_weight;
  Int_t nHits;
  Int_t dEdx;
  Int_t min_tracks;
  Int_t shift_terms;
  Int_t epdA_inner_row;
  Int_t epdA_outer_row;
  Int_t epdB_inner_row;
  Int_t epdB_outer_row;
  Double_t sqrt_s_NN;
  Double_t order_n;
  Double_t order_m;
  Double_t epd_threshold;
  Double_t tracking;
  Double_t dca;
  Double_t min_abs_tpc_eta;
  Double_t near_abs_tpc_eta;
  Double_t far_abs_tpc_eta;
  Double_t max_abs_tpc_eta;
  Double_t r_vtx;
  Double_t z_vtx_low;
  Double_t z_vtx_high;
  Double_t y_mid;
  Double_t nSig_pi_low;
  Double_t nSig_pi_high;
  Double_t nSig_ka_low;
  Double_t nSig_ka_high;
  Double_t nSig_pr_low;
  Double_t nSig_pr_high;
  Double_t m2_pi_low;
  Double_t m2_pi_high;
  Double_t m2_ka_low;
  Double_t m2_ka_high;
  Double_t y_mid_pi_low_wide;
  Double_t y_mid_pi_low;
  Double_t y_mid_pi_high;
  Double_t y_mid_pi_high_wide;
  Double_t y_mid_ka_low_wide;
  Double_t y_mid_ka_low;
  Double_t y_mid_ka_high;
  Double_t y_mid_ka_high_wide;
  Double_t y_mid_pr_low_wide;
  Double_t y_mid_pr_low;
  Double_t y_mid_pr_high;
  Double_t y_mid_pr_high_wide;
  Double_t pt_pi_low;
  Double_t pt_pi_high;
  Double_t pt_ka_low;
  Double_t pt_ka_high;
  Double_t pt_pr_low_wide;
  Double_t pt_pr_low;
  Double_t pt_pr_high;
  Double_t pt_pr_high_wide;
  TString order_n_str;
  TString order_m_str;

  bool errorFound();
  void notifyError();
  void read(std::string fileName);

 private:
  bool errorFlag;
  TString lastGoodKey;
  TString lastGoodValue;
};

#endif
