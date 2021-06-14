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
