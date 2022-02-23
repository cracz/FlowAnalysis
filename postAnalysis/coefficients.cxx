#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TStyle.h"
#include "PlotUtils.h"

void coefficients(TString jobID, TString order_n_str)
{
  //TH1::SetDefaultSumw2();
  
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {std::cout << "Wrong file!" << std::endl; return;}
  /*
  TFile *resolutionInfo_INPUT = TFile::Open("resolutionInfo_INPUT.root", "READ");
  if(!resolutionInfo_INPUT) { cout << "No resolution file found!" << endl; return; }

  TH1D *h_resolutions = (TH1D*)resolutionInfo_INPUT->Get("h_resolutions");
  TH2D *h2_resolutions = (TH2D*)resolutionInfo_INPUT->Get("h2_resolutions");
  */
  /*
  Double_t resolutionIDs[16] = {};
  Double_t resolutionIDsError[16] = {};
  
  for (int i = 1; i < h_resolutions->GetNbinsX(); i++)
    {
      resolutionIDs[i-1] = h_resolutions->GetBinContent(i);
      resolutionIDsError[i-1] = h_resolutions->GetBinError(i);
    }
  */
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
  //canvas->SetGridx();
  //canvas->SetGridy();
  canvas->SetLogy(0);
  canvas->SetTicks();
  //canvas->SetLeftMargin(0.15);
  canvas->cd();

  // Shaowei's data:
  double cent_label[7]={2.5, 7.5, 15, 25, 35, 45, 55};
  double cent_bin_lows[8]={0, 5, 10, 20, 30, 40, 50, 60};

  double v2val_cent_prp_data[7]={-0.00662406, -0.00906423, -0.0129282, -0.0202661, -0.0315561, -0.0467831, -0.0654673};
  double v2err_cent_prp_data[7]={0.000562772, 0.000228575, 9.699e-05, 8.17342e-05, 9.80836e-05, 0.000138984, 0.000257859};

  double v2val_cent_pip_data[7]={-0.0121035, -0.0200037, -0.0244652, -0.0279546, -0.0312375, -0.035842, -0.0413696};
  double v2err_cent_pip_data[7]={0.00136922, 0.000563976, 0.000240995, 0.000202438, 0.000239846, 0.000333549, 0.000602461};

  double v2val_cent_pim_data[7]={-0.00579955, -0.0125058, -0.0149095, -0.0170165, -0.0194754, -0.0216241, -0.0261989};
  double v2err_cent_pim_data[7]={0.00122163, 0.000505636, 0.000215966, 0.000181535, 0.000215232, 0.000298714, 0.000536848};

  double v2val_cent_kap_data[7]={-0.00341182, -0.0202501, -0.0132731, -0.0183739, -0.0281251, -0.0409518, -0.0722278};
  double v2err_cent_kap_data[7]={0.00604229, 0.00254523, 0.00113305, 0.00101225, 0.0012811, 0.00187359, 0.00345125};

  double v2val_cent_kam_data[7]={-0.0593462, -0.0051962, -0.0205428, -0.0224589, -0.0185285, -0.00424369, -0.0782483};
  double v2err_cent_kam_data[7]={0.0239337, 0.0106676, 0.00477639, 0.00426836, 0.00530658, 0.00769645, 0.0140865};
  ////////

  
  TH1D *sh_cent_pp = new TH1D("sh_cent_pp", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_pp->FillN(7, cent_label, v2val_cent_pip_data);
  for (int i = 1; i <= sh_cent_pp->GetNbinsX(); i++) { sh_cent_pp->SetBinError(i, v2err_cent_pip_data[i-1]); }
  
  TH1D *sh_cent_pm = new TH1D("sh_cent_pm", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_pm->FillN(7, cent_label, v2val_cent_pim_data);
  for (int i = 1; i <= sh_cent_pm->GetNbinsX(); i++) { sh_cent_pm->SetBinError(i, v2err_cent_pim_data[i-1]); }
  
  TH1D *sh_cent_kp = new TH1D("sh_cent_kp", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_kp->FillN(7, cent_label, v2val_cent_kap_data);
  for (int i = 1; i <= sh_cent_kp->GetNbinsX(); i++) { sh_cent_kp->SetBinError(i, v2err_cent_kap_data[i-1]); }
  
  TH1D *sh_cent_km = new TH1D("sh_cent_km", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_km->FillN(7, cent_label, v2val_cent_kam_data);
  for (int i = 1; i <= sh_cent_km->GetNbinsX(); i++) { sh_cent_km->SetBinError(i, v2err_cent_kam_data[i-1]); }
  
  TH1D *sh_cent_pr = new TH1D("sh_cent_pr", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_pr->FillN(7, cent_label, v2val_cent_prp_data);
  for (int i = 1; i <= sh_cent_pr->GetNbinsX(); i++) { sh_cent_pr->SetBinError(i, v2err_cent_prp_data[i-1]); }
  ////

  TProfile *p_vn_EpdA = (TProfile*)file->Get("p_vn_EpdA");
  TProfile *p_vn_EpdB = (TProfile*)file->Get("p_vn_EpdB");
  TProfile *p_vn_TpcB = (TProfile*)file->Get("p_vn_TpcB");
  TProfile *p_vn_Tpc  = (TProfile*)file->Get("p_vn_Tpc_pT_0p2to2");
  TProfile *p_vn_pp   = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm   = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp   = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km   = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr   = (TProfile*)file->Get("p_vn_pr");
  TProfile *p_vn_pr_alt   = (TProfile*)file->Get("p_vn_pr_alt");
  TProfile *p_vn_de   = (TProfile*)file->Get("p_vn_de");
  TProfile *p_vn_tr   = (TProfile*)file->Get("p_vn_tr");
  TProfile *p_vn_de_overA = (TProfile*)file->Get("p_vn_de_overA");
  TProfile *p_vn_tr_overA = (TProfile*)file->Get("p_vn_tr_overA");

  p_vn_kp->Rebin();
  p_vn_km->Rebin();

  TProfile *p_vn_pp_ext = (TProfile*)file->Get("p_vn_pp_ext");
  TProfile *p_vn_pm_ext = (TProfile*)file->Get("p_vn_pm_ext");
  TProfile *p_vn_kp_ext = (TProfile*)file->Get("p_vn_kp_ext");
  TProfile *p_vn_km_ext = (TProfile*)file->Get("p_vn_km_ext");
  TProfile *p_vn_pr_ext = (TProfile*)file->Get("p_vn_pr_ext");
  TProfile *p_vn_de_ext = (TProfile*)file->Get("p_vn_de_ext");
  TProfile *p_vn_tr_ext = (TProfile*)file->Get("p_vn_tr_ext");
  TProfile *p_vn_de_ext_overA = (TProfile*)file->Get("p_vn_de_ext_overA");
  TProfile *p_vn_tr_ext_overA = (TProfile*)file->Get("p_vn_tr_ext_overA");

  p_vn_kp_ext->Rebin();
  p_vn_km_ext->Rebin();

  TProfile *p_vn_pr_for = (TProfile*)file->Get("p_vn_pr_for");


  // Convert profiles to histograms
  TH1D *h_vn_EpdA = p_vn_EpdA->ProjectionX();
  TH1D *h_vn_EpdB = p_vn_EpdB->ProjectionX();
  TH1D *h_vn_TpcB = p_vn_TpcB->ProjectionX();
  TH1D *h_vn_Tpc = p_vn_Tpc->ProjectionX();

  TH1D *h_vn_pp = p_vn_pp->ProjectionX();
  TH1D *h_vn_pm = p_vn_pm->ProjectionX();
  TH1D *h_vn_kp = p_vn_kp->ProjectionX();
  TH1D *h_vn_km = p_vn_km->ProjectionX();
  TH1D *h_vn_pr = p_vn_pr->ProjectionX();
  TH1D *h_vn_pr_alt = p_vn_pr_alt->ProjectionX();
  TH1D *h_vn_de = p_vn_de->ProjectionX();
  TH1D *h_vn_tr = p_vn_tr->ProjectionX();
  TH1D *h_vn_de_overA = p_vn_de_overA->ProjectionX();
  TH1D *h_vn_tr_overA = p_vn_tr_overA->ProjectionX();

  TH1D *h_vn_pp_ext = p_vn_pp_ext->ProjectionX();
  TH1D *h_vn_pm_ext = p_vn_pm_ext->ProjectionX();
  TH1D *h_vn_kp_ext = p_vn_kp_ext->ProjectionX();
  TH1D *h_vn_km_ext = p_vn_km_ext->ProjectionX();
  TH1D *h_vn_pr_ext = p_vn_pr_ext->ProjectionX();
  TH1D *h_vn_de_ext = p_vn_de_ext->ProjectionX();
  TH1D *h_vn_tr_ext = p_vn_tr_ext->ProjectionX();
  TH1D *h_vn_de_ext_overA = p_vn_de_ext_overA->ProjectionX();
  TH1D *h_vn_tr_ext_overA = p_vn_tr_ext_overA->ProjectionX();
  
  TH1D *h_vn_pr_for = p_vn_pr_for->ProjectionX();

  // Flip centrality plots
  h_vn_EpdA = PlotUtils::flipHisto(h_vn_EpdA);
  h_vn_EpdB = PlotUtils::flipHisto(h_vn_EpdB);
  h_vn_TpcB = PlotUtils::flipHisto(h_vn_TpcB);
  h_vn_Tpc  = PlotUtils::flipHisto(h_vn_Tpc);

  h_vn_pp = PlotUtils::flipHisto(h_vn_pp);
  h_vn_pm = PlotUtils::flipHisto(h_vn_pm);
  h_vn_kp = PlotUtils::flipHisto(h_vn_kp);
  h_vn_km = PlotUtils::flipHisto(h_vn_km);
  h_vn_pr = PlotUtils::flipHisto(h_vn_pr);
  h_vn_pr_alt = PlotUtils::flipHisto(h_vn_pr_alt);
  h_vn_de = PlotUtils::flipHisto(h_vn_de);
  h_vn_tr = PlotUtils::flipHisto(h_vn_tr);
  h_vn_de_overA = PlotUtils::flipHisto(h_vn_de_overA);
  h_vn_tr_overA = PlotUtils::flipHisto(h_vn_tr_overA);

  h_vn_pp_ext = PlotUtils::flipHisto(h_vn_pp_ext);
  h_vn_pm_ext = PlotUtils::flipHisto(h_vn_pm_ext);
  h_vn_kp_ext = PlotUtils::flipHisto(h_vn_kp_ext);
  h_vn_km_ext = PlotUtils::flipHisto(h_vn_km_ext);
  h_vn_pr_ext = PlotUtils::flipHisto(h_vn_pr_ext);
  h_vn_de_ext = PlotUtils::flipHisto(h_vn_de_ext);
  h_vn_tr_ext = PlotUtils::flipHisto(h_vn_tr_ext);
  h_vn_de_ext_overA = PlotUtils::flipHisto(h_vn_de_ext_overA);
  h_vn_tr_ext_overA = PlotUtils::flipHisto(h_vn_tr_ext_overA);
  
  h_vn_pr_for = PlotUtils::flipHisto(h_vn_pr_for);

  // Trim and clean up x-axis
  h_vn_EpdA = PlotUtils::trimCentralityPlot(h_vn_EpdA);
  h_vn_EpdB = PlotUtils::trimCentralityPlot(h_vn_EpdB);
  h_vn_TpcB = PlotUtils::trimCentralityPlot(h_vn_TpcB);
  h_vn_Tpc  = PlotUtils::trimCentralityPlot(h_vn_Tpc);

  h_vn_pp = PlotUtils::trimCentralityPlot(h_vn_pp);
  h_vn_pm = PlotUtils::trimCentralityPlot(h_vn_pm);
  h_vn_kp = PlotUtils::trimCentralityPlot(h_vn_kp);
  h_vn_km = PlotUtils::trimCentralityPlot(h_vn_km);
  h_vn_pr = PlotUtils::trimCentralityPlot(h_vn_pr);
  h_vn_pr_alt = PlotUtils::trimCentralityPlot(h_vn_pr_alt);
  h_vn_de = PlotUtils::trimCentralityPlot(h_vn_de);
  h_vn_tr = PlotUtils::trimCentralityPlot(h_vn_tr);
  h_vn_de_overA = PlotUtils::trimCentralityPlot(h_vn_de_overA);
  h_vn_tr_overA = PlotUtils::trimCentralityPlot(h_vn_tr_overA);

  h_vn_pp_ext = PlotUtils::trimCentralityPlot(h_vn_pp_ext);
  h_vn_pm_ext = PlotUtils::trimCentralityPlot(h_vn_pm_ext);
  h_vn_kp_ext = PlotUtils::trimCentralityPlot(h_vn_kp_ext);
  h_vn_km_ext = PlotUtils::trimCentralityPlot(h_vn_km_ext);
  h_vn_pr_ext = PlotUtils::trimCentralityPlot(h_vn_pr_ext);
  h_vn_de_ext = PlotUtils::trimCentralityPlot(h_vn_de_ext);
  h_vn_tr_ext = PlotUtils::trimCentralityPlot(h_vn_tr_ext);
  h_vn_de_ext_overA = PlotUtils::trimCentralityPlot(h_vn_de_ext_overA);
  h_vn_tr_ext_overA = PlotUtils::trimCentralityPlot(h_vn_tr_ext_overA);
  
  h_vn_pr_for = PlotUtils::trimCentralityPlot(h_vn_pr_for);


  TH1D *h_vn_de_scaled = (TH1D*)h_vn_de->Clone();
  h_vn_de_scaled->Scale(1.0/2.0);
  
  TH1D *h_vn_de_ext_scaled = (TH1D*)h_vn_de_ext->Clone();
  h_vn_de_ext_scaled->Scale(1.0/2.0);

  TH1D *h_vn_tr_scaled = (TH1D*)h_vn_tr->Clone();
  h_vn_tr_scaled->Scale(1.0/3.0);
  
  TH1D *h_vn_tr_ext_scaled = (TH1D*)h_vn_tr_ext->Clone();
  h_vn_tr_ext_scaled->Scale(1.0/3.0);

  /*
  TH1D *h_vn_de_scaled = (TH1D*)h_vn_de_overA->Clone();
  h_vn_de_scaled->Scale(1.0/2.0);
  
  TH1D *h_vn_de_ext_scaled = (TH1D*)h_vn_de_ext_overA->Clone();
  h_vn_de_ext_scaled->Scale(1.0/2.0);

  TH1D *h_vn_tr_scaled = (TH1D*)h_vn_tr_overA->Clone();
  h_vn_tr_scaled->Scale(1.0/3.0);
  
  TH1D *h_vn_tr_ext_scaled = (TH1D*)h_vn_tr_ext_overA->Clone();
  h_vn_tr_ext_scaled->Scale(1.0/3.0);
  */
  
  THStack *piCentralityStack = new THStack("piCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kaCentralityStack = new THStack("kaCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");

  THStack *ppExtCentralityStack = new THStack("ppExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *pmExtCentralityStack = new THStack("pmExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kpExtCentralityStack = new THStack("kpExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kmExtCentralityStack = new THStack("kmExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *prExtCentralityStack = new THStack("prExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *deExtCentralityStack = new THStack("deExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *trExtCentralityStack = new THStack("trExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");

  THStack *pdtCentralityStack = new THStack("pdtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *pdtExtCentralityStack = new THStack("pdtExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *pdtScaledCentralityStack = new THStack("pdtScaledCentralityStack", ";Centrality (%);v_{"+order_n_str+"}/A");
  THStack *pdtExtScaledCentralityStack = new THStack("pdtExtScaledCentralityStack", ";Centrality (%);v_{"+order_n_str+"}/A");
  
  THStack *etaRegionStack = new THStack("etaRegionStack", ";Centrality (%);v_{"+order_n_str+"}");


  sh_cent_pp->SetMarkerStyle(25);
  sh_cent_pp->SetMarkerSize(2);
  sh_cent_pp->SetMarkerColor(2);
  sh_cent_pp->SetLineColor(2);

  sh_cent_pm->SetMarkerStyle(25);
  sh_cent_pm->SetMarkerSize(2);
  sh_cent_pm->SetMarkerColor(4);
  sh_cent_pm->SetLineColor(4);

  sh_cent_kp->SetMarkerStyle(25);
  sh_cent_kp->SetMarkerSize(2);
  sh_cent_kp->SetMarkerColor(2);
  sh_cent_kp->SetLineColor(2);

  sh_cent_km->SetMarkerStyle(25);
  sh_cent_km->SetMarkerSize(2);
  sh_cent_km->SetMarkerColor(4);
  sh_cent_km->SetLineColor(4);
  
  sh_cent_pr->SetMarkerStyle(25);
  sh_cent_pr->SetMarkerSize(2);
  sh_cent_pr->SetMarkerColor(2);
  sh_cent_pr->SetLineColor(2);


  
  h_vn_pp->SetMarkerStyle(20);
  h_vn_pp->SetMarkerSize(2.5);
  h_vn_pp->SetMarkerColor(2);
  h_vn_pp->SetLineColor(2);
  h_vn_pp->SetLineWidth(3);
  h_vn_pp->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pm->SetMarkerStyle(20);
  h_vn_pm->SetMarkerSize(2.5);
  h_vn_pm->SetMarkerColor(4);
  h_vn_pm->SetLineColor(4);
  h_vn_pm->SetLineWidth(3);
  h_vn_pm->GetYaxis()->SetTitleOffset(1.7);

  h_vn_kp->SetMarkerStyle(20);
  h_vn_kp->SetMarkerSize(2.5);
  h_vn_kp->SetMarkerColor(2);
  h_vn_kp->SetLineColor(2);
  h_vn_kp->SetLineWidth(3);
  h_vn_kp->GetYaxis()->SetTitleOffset(1.7);

  h_vn_km->SetMarkerStyle(20);
  h_vn_km->SetMarkerSize(2.5);
  h_vn_km->SetMarkerColor(4);
  h_vn_km->SetLineColor(4);
  h_vn_km->SetLineWidth(3);
  h_vn_km->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pr->SetMarkerStyle(20);
  h_vn_pr->SetMarkerSize(2.5);
  h_vn_pr->SetMarkerColor(2);
  h_vn_pr->SetLineColor(2);
  h_vn_pr->SetLineWidth(3);
  h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
  h_vn_pr->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_pr_alt->SetMarkerStyle(20);
  h_vn_pr_alt->SetMarkerSize(2.5);
  h_vn_pr_alt->SetMarkerColor(2);
  h_vn_pr_alt->SetLineColor(2);
  h_vn_pr_alt->SetLineWidth(3);
  h_vn_pr_alt->GetYaxis()->SetTitleOffset(1.7);
  h_vn_pr_alt->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_de->SetMarkerStyle(20);
  h_vn_de->SetMarkerSize(2.5);
  //h_vn_de->SetMarkerColor(2);
  //h_vn_de->SetLineColor(2);
  h_vn_de->SetLineWidth(3);
  h_vn_de->GetYaxis()->SetTitleOffset(1.7);
  h_vn_de->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_tr->SetMarkerStyle(20);
  h_vn_tr->SetMarkerSize(2.5);
  h_vn_tr->SetMarkerColor(44);
  h_vn_tr->SetLineColor(44);
  h_vn_tr->SetLineWidth(3);
  h_vn_tr->GetYaxis()->SetTitleOffset(1.7);
  h_vn_tr->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_de_overA->SetMarkerStyle(20);
  h_vn_de_overA->SetMarkerSize(2.5);
  //h_vn_de_overA->SetMarkerColor(2);
  //h_vn_de_overA->SetLineColor(2);
  h_vn_de_overA->SetLineWidth(3);
  h_vn_de_overA->GetYaxis()->SetTitleOffset(1.7);
  h_vn_de_overA->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_tr_overA->SetMarkerStyle(20);
  h_vn_tr_overA->SetMarkerSize(2.5);
  h_vn_tr_overA->SetMarkerColor(44);
  h_vn_tr_overA->SetLineColor(44);
  h_vn_tr_overA->SetLineWidth(3);
  h_vn_tr_overA->GetYaxis()->SetTitleOffset(1.7);
  h_vn_tr_overA->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_de_scaled->SetMarkerStyle(20);
  h_vn_de_scaled->SetMarkerSize(2.5);
  //h_vn_de_scaled->SetMarkerColor(2);
  //h_vn_de_scaled->SetLineColor(2);
  h_vn_de_scaled->SetLineWidth(3);
  h_vn_de_scaled->GetYaxis()->SetTitleOffset(1.7);
  h_vn_de_scaled->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_tr_scaled->SetMarkerStyle(20);
  h_vn_tr_scaled->SetMarkerSize(2.5);
  h_vn_tr_scaled->SetMarkerColor(44);
  h_vn_tr_scaled->SetLineColor(44);
  h_vn_tr_scaled->SetLineWidth(3);
  h_vn_tr_scaled->GetYaxis()->SetTitleOffset(1.7);
  h_vn_tr_scaled->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_de_ext_scaled->SetMarkerStyle(20);
  h_vn_de_ext_scaled->SetMarkerSize(2.5);
  //h_vn_de_ext_scaled->SetMarkerColor(2);
  //h_vn_de_ext_scaled->SetLineColor(2);
  h_vn_de_ext_scaled->SetLineWidth(3);
  h_vn_de_ext_scaled->GetYaxis()->SetTitleOffset(1.7);
  h_vn_de_ext_scaled->GetYaxis()->SetTitle("v_{"+order_n_str+"}");

  h_vn_tr_ext_scaled->SetMarkerStyle(20);
  h_vn_tr_ext_scaled->SetMarkerSize(2.5);
  h_vn_tr_ext_scaled->SetMarkerColor(44);
  h_vn_tr_ext_scaled->SetLineColor(44);
  h_vn_tr_ext_scaled->SetLineWidth(3);
  h_vn_tr_ext_scaled->GetYaxis()->SetTitleOffset(1.7);
  h_vn_tr_ext_scaled->GetYaxis()->SetTitle("v_{"+order_n_str+"}");


  h_vn_pp_ext->SetMarkerStyle(20);
  h_vn_pp_ext->SetMarkerSize(2.5);
  h_vn_pp_ext->SetLineWidth(3);
  h_vn_pp->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pm_ext->SetMarkerStyle(20);
  h_vn_pm_ext->SetMarkerSize(2.5);
  h_vn_pm_ext->SetLineWidth(3);
  h_vn_pm_ext->GetYaxis()->SetTitleOffset(1.7);

  h_vn_kp_ext->SetMarkerStyle(20);
  h_vn_kp_ext->SetMarkerSize(2.5);
  h_vn_kp_ext->SetLineWidth(3);
  h_vn_kp_ext->GetYaxis()->SetTitleOffset(1.7);

  h_vn_km_ext->SetMarkerStyle(20);
  h_vn_km_ext->SetMarkerSize(2.5);
  h_vn_km_ext->SetLineWidth(3);
  h_vn_km_ext->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pr_ext->SetMarkerStyle(20);
  h_vn_pr_ext->SetMarkerSize(2.5);
  h_vn_pr_ext->SetLineWidth(3);
  h_vn_pr_ext->SetMarkerColor(2);
  h_vn_pr_ext->SetLineColor(2);
  h_vn_pr_ext->GetYaxis()->SetTitleOffset(1.7);

  h_vn_de_ext->SetMarkerStyle(20);
  h_vn_de_ext->SetMarkerSize(2.5);
  h_vn_de_ext->SetLineWidth(3);
  h_vn_de_ext->GetYaxis()->SetTitleOffset(1.7);

  h_vn_tr_ext->SetMarkerStyle(20);
  h_vn_tr_ext->SetMarkerSize(2.5);
  h_vn_tr_ext->SetLineWidth(3);
  h_vn_tr_ext->SetMarkerColor(44);
  h_vn_tr_ext->SetLineColor(44);
  h_vn_tr_ext->GetYaxis()->SetTitleOffset(1.7);

  h_vn_de_ext_overA->SetMarkerStyle(20);
  h_vn_de_ext_overA->SetMarkerSize(2.5);
  h_vn_de_ext_overA->SetLineWidth(3);
  h_vn_de_ext_overA->GetYaxis()->SetTitleOffset(1.7);

  h_vn_tr_ext_overA->SetMarkerStyle(20);
  h_vn_tr_ext_overA->SetMarkerSize(2.5);
  h_vn_tr_ext_overA->SetLineWidth(3);
  h_vn_tr_ext_overA->SetMarkerColor(44);
  h_vn_tr_ext_overA->SetLineColor(44);
  h_vn_tr_ext_overA->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pr_for->SetMarkerStyle(20);
  h_vn_pr_for->SetMarkerSize(2.5);
  h_vn_pr_for->SetMarkerColor(4);
  h_vn_pr_for->SetLineWidth(3);
  h_vn_pr_for->GetYaxis()->SetTitleOffset(1.7);
  
  h_vn_EpdA->SetMarkerStyle(20);
  h_vn_EpdA->SetMarkerSize(2.5);
  h_vn_EpdA->SetMarkerColor(46);
  h_vn_EpdA->SetLineColor(46);
  h_vn_EpdA->SetLineWidth(3);
  h_vn_EpdA->GetYaxis()->SetTitleOffset(1.7);

  h_vn_EpdB->SetMarkerStyle(20);
  h_vn_EpdB->SetMarkerSize(2.5);
  h_vn_EpdB->SetMarkerColor(38);
  h_vn_EpdB->SetLineColor(38);
  h_vn_EpdB->SetLineWidth(3);
  h_vn_EpdB->GetYaxis()->SetTitleOffset(1.7);

  h_vn_TpcB->SetMarkerStyle(20);
  h_vn_TpcB->SetMarkerSize(2.5);
  h_vn_TpcB->SetMarkerColor(8);
  h_vn_TpcB->SetLineColor(8);
  h_vn_TpcB->SetLineWidth(3);
  h_vn_TpcB->GetYaxis()->SetTitleOffset(1.7);

  h_vn_Tpc->SetMarkerStyle(20);
  h_vn_Tpc->SetMarkerSize(2.5);
  h_vn_Tpc->SetLineWidth(3);
  h_vn_Tpc->GetYaxis()->SetTitleOffset(1.7);
  //vn_Tpc->SetMarkerColor(8);
  //vn_Tpc->SetLineColor(8);

  
  piCentralityStack->Add(h_vn_pp);
  piCentralityStack->Add(h_vn_pm);

  kaCentralityStack->Add(h_vn_kp);
  kaCentralityStack->Add(h_vn_km);


  ppExtCentralityStack->Add(h_vn_pp);
  ppExtCentralityStack->Add(h_vn_pp_ext);

  pmExtCentralityStack->Add(h_vn_pm);
  pmExtCentralityStack->Add(h_vn_pm_ext);

  kpExtCentralityStack->Add(h_vn_kp);
  kpExtCentralityStack->Add(h_vn_kp_ext);

  kmExtCentralityStack->Add(h_vn_km);
  kmExtCentralityStack->Add(h_vn_km_ext);

  deExtCentralityStack->Add(h_vn_de);
  deExtCentralityStack->Add(h_vn_de_ext);

  trExtCentralityStack->Add(h_vn_tr);
  trExtCentralityStack->Add(h_vn_tr_ext);

  pdtCentralityStack->Add(h_vn_pr_alt);
  pdtCentralityStack->Add(h_vn_de);
  pdtCentralityStack->Add(h_vn_tr);

  pdtExtCentralityStack->Add(h_vn_pr_ext);
  pdtExtCentralityStack->Add(h_vn_de_ext);
  pdtExtCentralityStack->Add(h_vn_tr_ext);

  pdtScaledCentralityStack->Add(h_vn_pr_alt);
  pdtScaledCentralityStack->Add(h_vn_de_scaled);
  pdtScaledCentralityStack->Add(h_vn_tr_scaled);

  pdtExtScaledCentralityStack->Add(h_vn_pr_ext);
  pdtExtScaledCentralityStack->Add(h_vn_de_ext_scaled);
  pdtExtScaledCentralityStack->Add(h_vn_tr_ext_scaled);

  /*
  pdtScaledCentralityStack->Add(h_vn_pr);
  pdtScaledCentralityStack->Add(h_vn_de_overA);
  pdtScaledCentralityStack->Add(h_vn_tr_overA);

  pdtExtScaledCentralityStack->Add(h_vn_pr_ext);
  pdtExtScaledCentralityStack->Add(h_vn_de_ext_overA);
  pdtExtScaledCentralityStack->Add(h_vn_tr_ext_overA);
  */
  etaRegionStack->Add(h_vn_EpdA);
  etaRegionStack->Add(h_vn_EpdB);
  etaRegionStack->Add(h_vn_TpcB);



  if (order_n_str == "2")
    {
      TLegend *piLegend = new TLegend(0.775, 0.7, 0.9, 0.85);
      piLegend->AddEntry(h_vn_pp,"#pi^{+}");
      piLegend->AddEntry(h_vn_pm,"#pi^{-}");
      piLegend->SetFillColorAlpha(0,0);
      piLegend->SetLineColorAlpha(0,0);

      TLegend *kaLegend = new TLegend(0.275, 0.26, 0.425, 0.39);
      kaLegend->AddEntry(h_vn_kp,"K^{+}");
      kaLegend->AddEntry(h_vn_km,"K^{-}");
      kaLegend->SetFillColorAlpha(0,0);
      kaLegend->SetLineColorAlpha(0,0);


      TLegend *ppExtLegend = new TLegend(0.55, 0.65, 0.85, 0.88);
      ppExtLegend->AddEntry(h_vn_pp,"#pi^{+}, 0 < y_{CM} < 0.5");
      ppExtLegend->AddEntry(h_vn_pp_ext,"#pi^{+}, 0.5 < y_{CM} < 1.0");
      ppExtLegend->SetFillColorAlpha(0,0);
      ppExtLegend->SetLineColorAlpha(0,0);

      TLegend *pmExtLegend = new TLegend(0.55, 0.65, 0.85, 0.88);
      pmExtLegend->AddEntry(h_vn_pm,"#pi^{-}, 0 < y_{CM} < 0.5");
      pmExtLegend->AddEntry(h_vn_pm_ext,"#pi^{-}, 0.5 < y_{CM} < 1.0");
      pmExtLegend->SetFillColorAlpha(0,0);
      pmExtLegend->SetLineColorAlpha(0,0);

      TLegend *kpExtLegend = new TLegend(0.55, 0.68, 0.85, 0.9);
      kpExtLegend->AddEntry(h_vn_kp,"K^{+}, 0 < y_{CM} < 0.5");
      kpExtLegend->AddEntry(h_vn_kp_ext,"K^{+}, 0.5 < y_{CM} < 1.0");
      kpExtLegend->SetFillColorAlpha(0,0);
      kpExtLegend->SetLineColorAlpha(0,0);

      TLegend *kmExtLegend = new TLegend(0.25, 0.13, 0.55, 0.28);
      kmExtLegend->AddEntry(h_vn_km,"K^{-}, 0 < y_{CM} < 0.5");
      kmExtLegend->AddEntry(h_vn_km_ext,"K^{-}, 0.5 < y_{CM} < 1.0");
      kmExtLegend->SetFillColorAlpha(0,0);
      kmExtLegend->SetLineColorAlpha(0,0);

      TLegend *prExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      prExtLegend->AddEntry(h_vn_pr_for,"Proton, -0.5 < y_{CM} < 0");
      prExtLegend->AddEntry(h_vn_pr,"Proton, 0 < y_{CM} < 0.5");
      prExtLegend->AddEntry(h_vn_pr_ext,"Proton, 0.5 < y_{CM} < 1.0");
      prExtLegend->SetFillColorAlpha(0,0);
      prExtLegend->SetLineColorAlpha(0,0);

      TLegend *deExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      deExtLegend->AddEntry(h_vn_de,"Deuteron, 0 < y_{CM} < 0.5");
      deExtLegend->AddEntry(h_vn_de_ext,"Deuteron, 0.5 < y_{CM} < 1.0");
      deExtLegend->SetFillColorAlpha(0,0);
      deExtLegend->SetLineColorAlpha(0,0);

      TLegend *trExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      trExtLegend->AddEntry(h_vn_tr,"Triton, 0 < y_{CM} < 0.5");
      trExtLegend->AddEntry(h_vn_tr_ext,"Triton, 0.5 < y_{CM} < 1.0");
      trExtLegend->SetFillColorAlpha(0,0);
      trExtLegend->SetLineColorAlpha(0,0);

      
      TLegend *etaLegend = new TLegend(0.65, 0.7, 0.9, 0.9);
      etaLegend->AddEntry(h_vn_EpdA, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(h_vn_EpdB, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(h_vn_TpcB, "TPC -1.0 < #eta < 0");

      
      TPaveText *piText = new TPaveText(15, -0.004, 45, 0.008, "NB");
      piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      piText->AddText("0 < y_{CM} < 0.5 GeV");
      piText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      piText->SetFillColorAlpha(0,0);
      piText->SetLineColorAlpha(0,0);

      TPaveText *kaText = new TPaveText(20, -0.1, 50, -0.07, "NB");
      kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kaText->AddText("0 < y_{CM} < 0.5 GeV");
      kaText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      kaText->SetFillColorAlpha(0,0);
      kaText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(5, -0.07, 35, -0.05, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);

      
      canvas->SetLeftMargin(0.13);
      gStyle->SetErrorX(0);
      gStyle->SetOptStat(0);

      TLine *zeroLine = new TLine(0, 0, 60, 0);
      zeroLine->SetLineStyle(9);

  
      piCentralityStack->Draw();
      piCentralityStack->GetXaxis()->SetNdivisions(210);
      piCentralityStack->SetMaximum(0.01);
      piCentralityStack->SetMinimum(-0.05);
      sh_cent_pp->SetMaximum(0.01);
      sh_cent_pp->SetMinimum(-0.05);
      sh_cent_pm->SetMaximum(0.01);
      sh_cent_pm->SetMinimum(-0.05);
      sh_cent_pp->Draw("E1P");
      sh_cent_pm->Draw("E1P SAME");
      piCentralityStack->Draw("NOSTACK E1P SAME");
      zeroLine->Draw("SAME");
      piLegend->Draw();
      piText->Draw();
      canvas->SaveAs(jobID + "_piCentralityStack.png");
      canvas->Clear();

      kaCentralityStack->Draw();
      kaCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kaCentralityStack->GetXaxis()->SetNdivisions(210);
      kaCentralityStack->SetMaximum(0.02);
      kaCentralityStack->SetMinimum(-0.12);
      sh_cent_kp->SetMaximum(0.02);
      sh_cent_kp->SetMinimum(-0.12);
      sh_cent_km->SetMaximum(0.02);
      sh_cent_km->SetMinimum(-0.12);
      sh_cent_kp->Draw("E1P");
      sh_cent_km->Draw("E1P SAME");
      kaCentralityStack->Draw("NOSTACK E1P SAME");
      zeroLine->Draw("SAME");
      kaLegend->Draw();
      kaText->Draw();
      canvas->SaveAs(jobID + "_kaCentralityStack.png");
      canvas->Clear();

      h_vn_pr->SetTitle("");
      h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pr->GetXaxis()->SetNdivisions(210);
      h_vn_pr->Draw("E1P");
      h_vn_pr->SetMaximum(0.01);
      h_vn_pr->SetMinimum(-0.09);
      sh_cent_pr->SetMaximum(0.01);
      sh_cent_pr->SetMinimum(-0.09);
      sh_cent_pr->SetMaximum(0.01);
      sh_cent_pr->SetMinimum(-0.09);
      sh_cent_pr->Draw("E1P SAME");
      zeroLine->Draw("SAME");
      prText->Draw();
      canvas->SaveAs(jobID + "_vn_pr.png");
      canvas->Clear();
     

      h_vn_EpdA->SetMarkerStyle(20);
      h_vn_EpdA->SetMarkerSize(2);
      //h_vn_EpdA->SetMarkerColor(2);
      //h_vn_EpdA->SetLineColor(2);
      h_vn_EpdA->GetXaxis()->SetNdivisions(210);
      h_vn_EpdA->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdA.png");
      canvas->Clear();

      h_vn_EpdB->SetMarkerStyle(20);
      h_vn_EpdB->SetMarkerSize(2);
      //h_vn_EpdB->SetMarkerColor(2);
      //h_vn_EpdB->SetLineColor(2);
      h_vn_EpdB->GetXaxis()->SetNdivisions(210);
      h_vn_EpdB->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdB.png");
      canvas->Clear();

      h_vn_TpcB->SetMarkerStyle(20);
      h_vn_TpcB->SetMarkerSize(2);
      //h_vn_TpcB->SetMarkerColor(2);
      //h_vn_TpcB->SetLineColor(2);
      h_vn_TpcB->GetXaxis()->SetNdivisions(210);
      h_vn_TpcB->GetYaxis()->SetTitleOffset(1.7);
      h_vn_TpcB->Draw("E1P");
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_TpcB.png");
      canvas->Clear();

      h_vn_Tpc->SetMarkerStyle(20);
      h_vn_Tpc->SetMarkerSize(2);
      h_vn_Tpc->GetXaxis()->SetNdivisions(210);
      h_vn_Tpc->GetYaxis()->SetTitleOffset(1.7);
      //h_vn_Tpc->SetMaximum(0.002);
      h_vn_Tpc->Draw("E1P");
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_Tpc.png");
      canvas->Clear();


      ppExtCentralityStack->Draw();
      ppExtCentralityStack->GetXaxis()->SetNdivisions(210);
      ppExtCentralityStack->SetMaximum(0.01);
      ppExtCentralityStack->SetMinimum(-0.05);
      ppExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      ppExtLegend->Draw();
      //ppExtText->Draw();
      canvas->SaveAs(jobID + "_ppExtCentralityStack.png");
      canvas->Clear();

      pmExtCentralityStack->Draw();
      pmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      pmExtCentralityStack->SetMaximum(0.01);
      pmExtCentralityStack->SetMinimum(-0.05);
      pmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pmExtLegend->Draw();
      //pmExtText->Draw();
      canvas->SaveAs(jobID + "_pmExtCentralityStack.png");
      canvas->Clear();

      kpExtCentralityStack->Draw();
      kpExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kpExtCentralityStack->SetMaximum(0.02);
      kpExtCentralityStack->SetMinimum(-0.12);
      kpExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kpExtLegend->Draw();
      //kpExtText->Draw();
      canvas->SaveAs(jobID + "_kpExtCentralityStack.png");
      canvas->Clear();

      kmExtCentralityStack->Draw();
      kmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kmExtCentralityStack->SetMaximum(0.02);
      kmExtCentralityStack->SetMinimum(-0.12);
      kmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kmExtLegend->Draw();
      //kmExtText->Draw();
      canvas->SaveAs(jobID + "_kmExtCentralityStack.png");
      canvas->Clear();


      h_vn_pr->SetMarkerColor(1);
      h_vn_pr->SetLineColor(1);

      prExtCentralityStack->Add(h_vn_pr_for);
      prExtCentralityStack->Add(h_vn_pr);
      prExtCentralityStack->Add(h_vn_pr_ext);
      
      prExtCentralityStack->Draw();
      prExtCentralityStack->GetXaxis()->SetNdivisions(210);
      prExtCentralityStack->SetMaximum(0.03);
      prExtCentralityStack->SetMinimum(-0.09);
      prExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      prExtLegend->Draw();
      //prExtText->Draw();
      canvas->SaveAs(jobID + "_prExtCentralityStack.png");
      canvas->Clear();

      h_vn_pr->SetMarkerColor(2);
      h_vn_pr->SetLineColor(2);


      deExtCentralityStack->Draw();
      deExtCentralityStack->GetXaxis()->SetNdivisions(210);
      deExtCentralityStack->SetMaximum(0.03);
      deExtCentralityStack->SetMinimum(-0.09);
      deExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      deExtLegend->Draw();
      //deExtText->Draw();
      canvas->SaveAs(jobID + "_deExtCentralityStack.png");
      canvas->Clear();

      trExtCentralityStack->Draw();
      trExtCentralityStack->GetXaxis()->SetNdivisions(210);
      trExtCentralityStack->SetMaximum(0.03);
      trExtCentralityStack->SetMinimum(-0.09);
      trExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      trExtLegend->Draw();
      //trExtText->Draw();
      canvas->SaveAs(jobID + "_trExtCentralityStack.png");
      canvas->Clear();

      
      etaRegionStack->Draw();
      etaRegionStack->GetYaxis()->SetTitleOffset(1.7);
      etaRegionStack->GetXaxis()->SetNdivisions(210);
      etaRegionStack->Draw();
      //etaRegionStack->SetMaximum(0.3);
      etaRegionStack->SetMinimum(-0.1);
      etaRegionStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      etaLegend->Draw();
      canvas->SaveAs(jobID + "_etaRegionStack.png");
      canvas->Clear();

      // Zoom in on last one
      etaRegionStack->Draw();
      etaRegionStack->GetYaxis()->SetTitleOffset(1.7);
      etaRegionStack->GetXaxis()->SetNdivisions(210);
      etaRegionStack->Draw();
      etaRegionStack->SetMaximum(0.1);
      etaRegionStack->SetMinimum(-0.1);
      etaRegionStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      etaLegend->Draw();
      canvas->SaveAs(jobID + "_etaRegionStack_zoom.png");
      canvas->Clear();
    }
  else if (order_n_str == "3")
    {
      TLegend *piLegend = new TLegend(0.67, 0.71, 0.805, 0.87);
      piLegend->AddEntry(h_vn_pp,"#pi^{+}");
      piLegend->AddEntry(h_vn_pm,"#pi^{-}");
      piLegend->SetFillColorAlpha(0,0);
      piLegend->SetLineColorAlpha(0,0);

      TLegend *kaLegend = new TLegend(0.7, 0.72, 0.8, 0.87);
      kaLegend->AddEntry(h_vn_kp,"K^{+}");
      kaLegend->AddEntry(h_vn_km,"K^{-}");
      kaLegend->SetFillColorAlpha(0,0);
      kaLegend->SetLineColorAlpha(0,0);

      TLegend *pdtLegend = new TLegend(0.67, 0.71, 0.805, 0.87);
      pdtLegend->AddEntry(h_vn_pr,"p");
      pdtLegend->AddEntry(h_vn_de,"d");
      pdtLegend->AddEntry(h_vn_tr,"t");
      pdtLegend->SetFillColorAlpha(0,0);
      pdtLegend->SetLineColorAlpha(0,0);

      TLegend *prLegend = new TLegend(0.7, 0.7, 0.83, 0.78);
      prLegend->AddEntry(h_vn_pr,"p");
      prLegend->SetFillColorAlpha(0,0);
      prLegend->SetLineColorAlpha(0,0);

      TLegend *deLegend = new TLegend(0.7, 0.7, 0.83, 0.78);
      deLegend->AddEntry(h_vn_de,"d");
      deLegend->SetFillColorAlpha(0,0);
      deLegend->SetLineColorAlpha(0,0);

      TLegend *trLegend = new TLegend(0.7, 0.7, 0.83, 0.78);
      trLegend->AddEntry(h_vn_tr,"t");
      trLegend->SetFillColorAlpha(0,0);
      trLegend->SetLineColorAlpha(0,0);


      TLegend *ppExtLegend = new TLegend(0.4, 0.62, 0.7, 0.82);
      ppExtLegend->AddEntry(h_vn_pp,"#pi^{+}, 0 < y_{CM} < 0.5");
      ppExtLegend->AddEntry(h_vn_pp_ext,"#pi^{+}, 0.5 < y_{CM} < 1.0");
      ppExtLegend->SetFillColorAlpha(0,0);
      ppExtLegend->SetLineColorAlpha(0,0);

      TLegend *pmExtLegend = new TLegend(0.15, 0.67, 0.45, 0.9);
      pmExtLegend->AddEntry(h_vn_pm,"#pi^{-}, 0 < y_{CM} < 0.5");
      pmExtLegend->AddEntry(h_vn_pm_ext,"#pi^{-}, 0.5 < y_{CM} < 1.0");
      pmExtLegend->SetFillColorAlpha(0,0);
      pmExtLegend->SetLineColorAlpha(0,0);

      TLegend *kpExtLegend = new TLegend(0.55, 0.7, 0.85, 0.9);
      kpExtLegend->AddEntry(h_vn_kp,"K^{+}, 0 < y_{CM} < 0.5");
      kpExtLegend->AddEntry(h_vn_kp_ext,"K^{+}, 0.5 < y_{CM} < 1.0");
      kpExtLegend->SetFillColorAlpha(0,0);
      kpExtLegend->SetLineColorAlpha(0,0);

      TLegend *kmExtLegend = new TLegend(0.28, 0.68, 0.55, 0.85);
      kmExtLegend->AddEntry(h_vn_km,"K^{-}, 0 < y_{CM} < 0.5");
      kmExtLegend->AddEntry(h_vn_km_ext,"K^{-}, 0.5 < y_{CM} < 1.0");
      kmExtLegend->SetFillColorAlpha(0,0);
      kmExtLegend->SetLineColorAlpha(0,0);

      TLegend *prExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      prExtLegend->AddEntry(h_vn_pr_for,"Proton, -0.5 < y_{CM} < 0");
      prExtLegend->AddEntry(h_vn_pr,"Proton, 0 < y_{CM} < 0.5");
      prExtLegend->AddEntry(h_vn_pr_ext,"Proton, 0.5 < y_{CM} < 1.0");
      prExtLegend->SetFillColorAlpha(0,0);
      prExtLegend->SetLineColorAlpha(0,0);

      TLegend *deExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      deExtLegend->AddEntry(h_vn_de,"Deuteron, 0 < y_{CM} < 0.5");
      deExtLegend->AddEntry(h_vn_de_ext,"Deuteron, 0.5 < y_{CM} < 1.0");
      deExtLegend->SetFillColorAlpha(0,0);
      deExtLegend->SetLineColorAlpha(0,0);

      TLegend *trExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      trExtLegend->AddEntry(h_vn_tr,"Triton, 0 < y_{CM} < 0.5");
      trExtLegend->AddEntry(h_vn_tr_ext,"Triton, 0.5 < y_{CM} < 1.0");
      trExtLegend->SetFillColorAlpha(0,0);
      trExtLegend->SetLineColorAlpha(0,0);

      TLegend *pdtExtLegend = new TLegend(0.67, 0.71, 0.805, 0.87);
      pdtExtLegend->AddEntry(h_vn_pr_ext,"p");
      pdtExtLegend->AddEntry(h_vn_de_ext,"d");
      pdtExtLegend->AddEntry(h_vn_tr_ext,"t");
      pdtExtLegend->SetFillColorAlpha(0,0);
      pdtExtLegend->SetLineColorAlpha(0,0);

      TLegend *etaLegend = new TLegend(0.65, 0.25, 0.9, 0.45);
      etaLegend->AddEntry(h_vn_EpdA, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(h_vn_EpdB, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(h_vn_TpcB, "TPC -1.0 < #eta < 0");

      

      TPaveText *piText = new TPaveText(5, 0.055, 38, 0.1, "NB");
      piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      piText->AddText("0 < y_{CM} < 0.5 GeV");
      piText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      piText->SetFillColorAlpha(0,0);
      piText->SetLineColorAlpha(0,0);

      TPaveText *kaText = new TPaveText(5, 0.055, 38, 0.1, "NB");
      kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kaText->AddText("0 < y_{CM} < 0.5 GeV");
      kaText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      kaText->SetFillColorAlpha(0,0);
      kaText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(5, 0.055, 38, 0.1, "NB");
      //prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);

      TPaveText *deText = new TPaveText(5, 0.055, 38, 0.1, "NB");
      //deText->AddText("Deuteron");
      deText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      deText->AddText("0.5 < y_{CM} < 1.0 GeV");
      deText->AddText("0.3 #leq p_{T} #leq 1.0 GeV");
      deText->SetFillColorAlpha(0,0);
      deText->SetLineColorAlpha(0,0);

      TPaveText *trText = new TPaveText(5, 0.055, 38, 0.1, "NB");
      //trText->AddText("Triton");
      trText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      trText->AddText("0.5 < y_{CM} < 1.0 GeV");
      trText->AddText("0.3 #leq p_{T} #leq 1.0 GeV");
      trText->SetFillColorAlpha(0,0);
      trText->SetLineColorAlpha(0,0);

      TPaveText *pdtText = new TPaveText(5, -0.12, 38, -0.08, "NB");
      pdtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pdtText->AddText("0.5 < y_{CM} < 1.0 GeV");
      pdtText->AddText("0.3 #leq p_{T}/A #leq 1.0 GeV");
      pdtText->SetFillColorAlpha(0,0);
      pdtText->SetLineColorAlpha(0,0);

      TPaveText *pdtExtText = new TPaveText(5, -0.12, 38, -0.08, "NB");
      pdtExtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pdtExtText->AddText("0.5 < y_{CM} < 1.0 GeV");
      pdtExtText->AddText("0.4 #leq p_{T}/A #leq 2.0 GeV");
      pdtExtText->SetFillColorAlpha(0,0);
      pdtExtText->SetLineColorAlpha(0,0);

      TPaveText *pdtScaledText = new TPaveText(5, -0.08, 38, -0.04, "NB");
      pdtScaledText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pdtScaledText->AddText("0.5 < y_{CM} < 1.0 GeV");
      pdtScaledText->AddText("0.3 #leq p_{T}/A #leq 1.0 GeV");
      pdtScaledText->SetFillColorAlpha(0,0);
      pdtScaledText->SetLineColorAlpha(0,0);

      TPaveText *pdtExtScaledText = new TPaveText(5, -0.12, 38, -0.08, "NB");
      pdtExtScaledText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pdtExtScaledText->AddText("0.5 < y_{CM} < 1.0 GeV");
      pdtExtScaledText->AddText("0.4 #leq p_{T}/A #leq 2.0 GeV");
      pdtExtScaledText->SetFillColorAlpha(0,0);
      pdtExtScaledText->SetLineColorAlpha(0,0);


      canvas->SetLeftMargin(0.13);
      gStyle->SetErrorX(0);
      gStyle->SetOptStat(0);
      gStyle->SetEndErrorSize(6);


      TLine *zeroLine = new TLine(0, 0, 60, 0);
      zeroLine->SetLineStyle(9);

      Double_t centralityUpperBounds = 0.12;
      Double_t centralityLowerBounds = -0.07;
      
      piCentralityStack->Draw();
      piCentralityStack->GetYaxis()->SetTitleOffset(1.8);
      piCentralityStack->GetXaxis()->SetNdivisions(210);
      piCentralityStack->SetMaximum(centralityUpperBounds);
      piCentralityStack->SetMinimum(centralityLowerBounds);
      piCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      piCentralityStack->Draw("NOSTACK E1P SAME");
      piLegend->Draw();
      piText->Draw();
      canvas->SaveAs(jobID + "_piCentralityStack.png");
      canvas->Clear();

      // Zoomed in pions
      piCentralityStack->SetMaximum(0.03);
      piCentralityStack->SetMinimum(-0.02);
      piCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      piCentralityStack->Draw("NOSTACK E1P SAME");
      canvas->SaveAs(jobID + "_piCentralityStack_ZOOM.png");
      canvas->Clear();
      //


      kaCentralityStack->Draw();
      kaCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kaCentralityStack->GetXaxis()->SetNdivisions(210);
      kaCentralityStack->SetMaximum(centralityUpperBounds);
      kaCentralityStack->SetMinimum(centralityLowerBounds);
      kaCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kaCentralityStack->Draw("NOSTACK E1P SAME");
      kaLegend->Draw();
      kaText->Draw();
      canvas->SaveAs(jobID + "_kaCentralityStack.png");
      canvas->Clear();

      
      pdtCentralityStack->Draw();
      pdtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      pdtCentralityStack->GetXaxis()->SetNdivisions(210);
      pdtCentralityStack->SetMaximum(0.07);
      pdtCentralityStack->SetMinimum(-0.15);
      pdtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pdtCentralityStack->Draw("NOSTACK E1P SAME");
      pdtLegend->Draw();
      pdtText->Draw();
      canvas->SaveAs(jobID + "_pdtCentralityStack.png");
      canvas->Clear();
      /*
      pdtExtCentralityStack->Draw();
      pdtExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      pdtExtCentralityStack->GetXaxis()->SetNdivisions(210);
      pdtExtCentralityStack->SetMaximum(0.07);
      pdtExtCentralityStack->SetMinimum(-0.14);
      pdtExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pdtExtCentralityStack->Draw("NOSTACK E1P SAME");
      pdtExtLegend->Draw();
      pdtExtText->Draw();
      canvas->SaveAs(jobID + "_pdtExtCentralityStack.png");
      canvas->Clear();

      pdtExtScaledCentralityStack->Draw();
      pdtExtScaledCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      pdtExtScaledCentralityStack->GetXaxis()->SetNdivisions(210);
      pdtExtScaledCentralityStack->SetMaximum(0.07);
      pdtExtScaledCentralityStack->SetMinimum(-0.14);
      pdtExtScaledCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pdtExtScaledCentralityStack->Draw("NOSTACK E1P SAME");
      pdtExtLegend->Draw();
      pdtExtScaledText->Draw();
      canvas->SaveAs(jobID + "_pdtExtScaledCentralityStack.png");
      canvas->Clear();
      */
      pdtScaledCentralityStack->Draw();
      pdtScaledCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      pdtScaledCentralityStack->GetXaxis()->SetNdivisions(210);
      pdtScaledCentralityStack->SetMaximum(0.07);
      pdtScaledCentralityStack->SetMinimum(-0.1);
      pdtScaledCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pdtScaledCentralityStack->Draw("NOSTACK E1P SAME");
      pdtLegend->Draw();
      pdtScaledText->Draw();
      canvas->SaveAs(jobID + "_pdtScaledCentralityStack.png");
      canvas->Clear();

      h_vn_pr->SetTitle("");
      h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pr->GetXaxis()->SetNdivisions(210);
      h_vn_pr->Draw("E1P");
      h_vn_pr->SetMaximum(centralityUpperBounds);
      h_vn_pr->SetMinimum(centralityLowerBounds);
      zeroLine->Draw("SAME");
      h_vn_pr->Draw("E1P SAME");
      prLegend->Draw();
      prText->Draw();
      canvas->SaveAs(jobID + "_vn_pr.png");
      canvas->Clear();

      //Zoomed in proton
      h_vn_pr->SetMaximum(0.01);
      h_vn_pr->SetMinimum(-0.03);
      h_vn_pr->Draw("E1P");
      zeroLine->Draw("SAME");
      h_vn_pr->Draw("E1P SAME");
      canvas->SaveAs(jobID + "_vn_pr_ZOOM.png");
      canvas->Clear();
      //////


      h_vn_de->SetTitle("");
      h_vn_de->GetYaxis()->SetTitleOffset(1.7);
      h_vn_de->GetXaxis()->SetNdivisions(210);
      h_vn_de->Draw("E1P");
      h_vn_de->SetMaximum(centralityUpperBounds);
      h_vn_de->SetMinimum(centralityLowerBounds);
      zeroLine->Draw("SAME");
      h_vn_de->Draw("E1P SAME");
      deLegend->Draw();
      deText->Draw();
      canvas->SaveAs(jobID + "_vn_de.png");
      canvas->Clear();

      //Zoomed in Deuteron
      h_vn_de->SetMaximum(0.03);
      h_vn_de->SetMinimum(-0.04);
      h_vn_de->Draw("E1P");
      zeroLine->Draw("SAME");
      h_vn_de->Draw("E1P SAME");
      canvas->SaveAs(jobID + "_vn_de_ZOOM.png");
      canvas->Clear();
      //////


      h_vn_tr->SetTitle("");
      h_vn_tr->GetYaxis()->SetTitleOffset(1.7);
      h_vn_tr->GetXaxis()->SetNdivisions(210);
      h_vn_tr->Draw("E1P");
      h_vn_tr->SetMaximum(centralityUpperBounds);
      h_vn_tr->SetMinimum(-0.15);
      zeroLine->Draw("SAME");
      h_vn_tr->Draw("E1P SAME");
      trLegend->Draw();
      trText->Draw();
      canvas->SaveAs(jobID + "_vn_tr.png");
      canvas->Clear();

      //Zoomed in Triton
      h_vn_tr->SetMaximum(0.05);
      h_vn_tr->SetMinimum(-0.15);
      h_vn_tr->Draw("E1P");
      zeroLine->Draw("SAME");
      h_vn_tr->Draw("E1P SAME");
      canvas->SaveAs(jobID + "_vn_tr_ZOOM.png");
      canvas->Clear();
      //////

      h_vn_EpdA->SetMarkerStyle(20);
      h_vn_EpdA->SetMarkerSize(2);
      //h_vn_EpdA->SetMarkerColor(2);
      //h_vn_EpdA->SetLineColor(2);
      h_vn_EpdA->GetYaxis()->SetTitleOffset(1.7);
      h_vn_EpdA->GetXaxis()->SetNdivisions(210);
      h_vn_EpdA->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdA.png");
      canvas->Clear();

      h_vn_EpdB->SetMarkerStyle(20);
      h_vn_EpdB->SetMarkerSize(2);
      //h_vn_EpdB->SetMarkerColor(2);
      //h_vn_EpdB->SetLineColor(2);
      h_vn_EpdB->GetYaxis()->SetTitleOffset(1.8);
      h_vn_EpdB->GetXaxis()->SetNdivisions(210);
      h_vn_EpdB->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdB.png");
      canvas->Clear();

      h_vn_TpcB->SetMarkerStyle(20);
      h_vn_TpcB->SetMarkerSize(2);
      //h_vn_TpcB->SetMarkerColor(2);
      //h_vn_TpcB->SetLineColor(2);
      h_vn_TpcB->GetXaxis()->SetNdivisions(210);
      h_vn_TpcB->GetYaxis()->SetTitleOffset(1.7);
      h_vn_TpcB->Draw("E1P");
      h_vn_TpcB->SetMaximum(0.01);
      h_vn_TpcB->SetMinimum(-0.03);
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_TpcB.png");
      canvas->Clear();

      h_vn_Tpc->SetMarkerStyle(20);
      h_vn_Tpc->SetMarkerSize(2);
      h_vn_Tpc->GetXaxis()->SetNdivisions(210);
      h_vn_Tpc->GetYaxis()->SetTitleOffset(1.7);
      h_vn_Tpc->SetMaximum(0.002);
      h_vn_Tpc->Draw("E1P");
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_Tpc.png");
      canvas->Clear();


      ppExtCentralityStack->Draw();
      ppExtCentralityStack->GetYaxis()->SetTitleOffset(1.8);
      ppExtCentralityStack->GetXaxis()->SetNdivisions(210);
      ppExtCentralityStack->SetMaximum(0.03);
      ppExtCentralityStack->SetMinimum(-0.015);
      ppExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      ppExtCentralityStack->Draw("NOSTACK E1P SAME");
      ppExtLegend->Draw();
      //ppExtText->Draw();
      canvas->SaveAs(jobID + "_ppExtCentralityStack.png");
      canvas->Clear();

      pmExtCentralityStack->Draw();
      pmExtCentralityStack->GetYaxis()->SetTitleOffset(1.8);
      pmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      pmExtCentralityStack->SetMaximum(0.02);
      pmExtCentralityStack->SetMinimum(-0.015);
      pmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pmExtCentralityStack->Draw("NOSTACK E1P SAME");
      pmExtLegend->Draw();
      //pmExtText->Draw();
      canvas->SaveAs(jobID + "_pmExtCentralityStack.png");
      canvas->Clear();

      kpExtCentralityStack->Draw();
      kpExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kpExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kpExtCentralityStack->SetMaximum(0.2);
      kpExtCentralityStack->SetMinimum(-0.2);
      kpExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kpExtCentralityStack->Draw("NOSTACK E1P SAME");
      kpExtLegend->Draw();
      //kpExtText->Draw();
      canvas->SaveAs(jobID + "_kpExtCentralityStack.png");
      canvas->Clear();

      kmExtCentralityStack->Draw();
      kmExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kmExtCentralityStack->SetMaximum(0.4);
      kmExtCentralityStack->SetMinimum(-0.2);
      kmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kmExtCentralityStack->Draw("NOSTACK E1P SAME");
      kmExtLegend->Draw();
      //kmExtText->Draw();
      canvas->SaveAs(jobID + "_kmExtCentralityStack.png");
      canvas->Clear();

      h_vn_pr->SetMarkerColor(1);
      h_vn_pr->SetLineColor(1);

      prExtCentralityStack->Add(h_vn_pr_for);
      prExtCentralityStack->Add(h_vn_pr);
      prExtCentralityStack->Add(h_vn_pr_ext);

      prExtCentralityStack->Draw();
      prExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      prExtCentralityStack->GetXaxis()->SetNdivisions(210);
      prExtCentralityStack->SetMaximum(0.04);
      prExtCentralityStack->SetMinimum(-0.08);
      prExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      prExtCentralityStack->Draw("NOSTACK E1P SAME");
      prExtLegend->Draw();
      //prExtText->Draw();
      canvas->SaveAs(jobID + "_prExtCentralityStack.png");
      canvas->Clear();

      h_vn_pr->SetMarkerColor(2);
      h_vn_pr->SetLineColor(2);


      deExtCentralityStack->Draw();
      deExtCentralityStack->GetXaxis()->SetNdivisions(210);
      deExtCentralityStack->SetMaximum(0.03);
      deExtCentralityStack->SetMinimum(-0.09);
      deExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      deExtLegend->Draw();
      //deExtText->Draw();
      canvas->SaveAs(jobID + "_deExtCentralityStack.png");
      canvas->Clear();

      trExtCentralityStack->Draw();
      trExtCentralityStack->GetXaxis()->SetNdivisions(210);
      trExtCentralityStack->SetMaximum(0.03);
      trExtCentralityStack->SetMinimum(-0.09);
      trExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      trExtLegend->Draw();
      //trExtText->Draw();
      canvas->SaveAs(jobID + "_trExtCentralityStack.png");
      canvas->Clear();

      etaRegionStack->Draw();
      etaRegionStack->GetYaxis()->SetTitleOffset(1.7);
      etaRegionStack->GetXaxis()->SetNdivisions(210);
      etaRegionStack->Draw();
      etaRegionStack->SetMaximum(0.05);
      //etaRegionStack->SetMinimum(-0.1);
      etaRegionStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      etaRegionStack->Draw("NOSTACK E1P SAME");
      etaLegend->Draw();
      canvas->SaveAs(jobID + "_etaRegionStack.png");
      canvas->Clear();
    }

  //resolutionInfo_INPUT->Close();
  file->Close();
}





/*
void applyResolution(TH1D *histogram, Double_t resolution, Double_t resolutionError)
{
  for (int i = 1; i < histogram->GetNbinsX(); i++)
    {
      Double_t rawBinContent = histogram->GetBinContent(i);
      if (rawBinContent == 0.0) continue;
      
      Double_t rawBinError = histogram->GetBinError(i);

      Double_t newBinContent = rawBinContent/resolution;
      Double_t newBinError = newBinContent * TMath::Sqrt( pow(rawBinError/rawBinContent, 2) + pow(resolutionError/resolution, 2) );

      histogram->SetBinContent(i, newBinContent);
      histogram->SetBinError(i, newBinError);
    }
}
*/
