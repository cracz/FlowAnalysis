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

  Double_t yCM_pid_pi_low;
  Double_t yCM_pid_pi_high;
  Double_t yCM_flow_pi_low;
  Double_t yCM_flow_pi_high;
  Double_t yCM_ext_flow_pi_low;
  Double_t yCM_ext_flow_pi_high;

  Double_t yCM_pid_ka_low;
  Double_t yCM_pid_ka_high;
  Double_t yCM_flow_ka_low;
  Double_t yCM_flow_ka_high;
  Double_t yCM_ext_flow_ka_low;
  Double_t yCM_ext_flow_ka_high;

  Double_t yCM_pid_pr_low;
  Double_t yCM_pid_pr_high;
  Double_t yCM_flow_pr_low;
  Double_t yCM_flow_pr_high;
  Double_t yCM_dep_flow_pr_low;
  Double_t yCM_dep_flow_pr_high;
  Double_t yCM_ext_flow_pr_low;
  Double_t yCM_ext_flow_pr_high;
  Double_t yCM_sym_flow_pr_low;
  Double_t yCM_sym_flow_pr_high;
  Double_t yCM_for_flow_pr_low;
  Double_t yCM_for_flow_pr_high;

  Double_t pt_pid_pi_low;
  Double_t pt_pid_pi_high;
  Double_t pt_pid_ka_low;
  Double_t pt_pid_ka_high;
  Double_t pt_pid_pr_low;
  Double_t pt_pid_pr_high;
  Double_t pt_flow_pr_low;
  Double_t pt_flow_pr_high;
  Double_t pt_ydep_flow_pr_low;
  Double_t pt_ydep_flow_pr_high;
  Double_t pt_ext_flow_pr_low;
  Double_t pt_ext_flow_pr_high;
  Double_t pt_sym_flow_pr_low;
  Double_t pt_sym_flow_pr_high;
  Double_t pt_for_flow_pr_low;
  Double_t pt_for_flow_pr_high;

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
