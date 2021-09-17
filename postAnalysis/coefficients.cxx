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
  canvas->SetGridx();
  canvas->SetGridy();
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

  TProfile *p_vn_EpdE = (TProfile*)file->Get("p_vn_EpdE");
  TProfile *p_vn_EpdF = (TProfile*)file->Get("p_vn_EpdF");
  TProfile *p_vn_TpcB = (TProfile*)file->Get("p_vn_TpcB");
  TProfile *p_vn_Tpc  = (TProfile*)file->Get("p_vn_Tpc_pT_0p2to2");
  TProfile *p_vn_pp   = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm   = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp   = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km   = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr   = (TProfile*)file->Get("p_vn_pr");

  p_vn_kp->Rebin();
  p_vn_km->Rebin();

  TProfile *p_vn_pp_ext = (TProfile*)file->Get("p_vn_pp_ext");
  TProfile *p_vn_pm_ext = (TProfile*)file->Get("p_vn_pm_ext");
  TProfile *p_vn_kp_ext = (TProfile*)file->Get("p_vn_kp_ext");
  TProfile *p_vn_km_ext = (TProfile*)file->Get("p_vn_km_ext");
  TProfile *p_vn_pr_ext = (TProfile*)file->Get("p_vn_pr_ext");

  p_vn_kp_ext->Rebin();
  p_vn_km_ext->Rebin();

  TProfile *p_vn_pr_for = (TProfile*)file->Get("p_vn_pr_for");


  // Convert profiles to histograms
  TH1D *h_vn_EpdE = p_vn_EpdE->ProjectionX();
  TH1D *h_vn_EpdF = p_vn_EpdF->ProjectionX();
  TH1D *h_vn_TpcB = p_vn_TpcB->ProjectionX();
  TH1D *h_vn_Tpc = p_vn_Tpc->ProjectionX();

  TH1D *h_vn_pp = p_vn_pp->ProjectionX();
  TH1D *h_vn_pm = p_vn_pm->ProjectionX();
  TH1D *h_vn_kp = p_vn_kp->ProjectionX();
  TH1D *h_vn_km = p_vn_km->ProjectionX();
  TH1D *h_vn_pr = p_vn_pr->ProjectionX();

  TH1D *h_vn_pp_ext = p_vn_pp_ext->ProjectionX();
  TH1D *h_vn_pm_ext = p_vn_pm_ext->ProjectionX();
  TH1D *h_vn_kp_ext = p_vn_kp_ext->ProjectionX();
  TH1D *h_vn_km_ext = p_vn_km_ext->ProjectionX();
  TH1D *h_vn_pr_ext = p_vn_pr_ext->ProjectionX();
  
  TH1D *h_vn_pr_for = p_vn_pr_for->ProjectionX();

  // Flip centrality plots
  h_vn_EpdE = PlotUtils::flipHisto(h_vn_EpdE);
  h_vn_EpdF = PlotUtils::flipHisto(h_vn_EpdF);
  h_vn_TpcB = PlotUtils::flipHisto(h_vn_TpcB);
  h_vn_Tpc  = PlotUtils::flipHisto(h_vn_Tpc);

  h_vn_pp = PlotUtils::flipHisto(h_vn_pp);
  h_vn_pm = PlotUtils::flipHisto(h_vn_pm);
  h_vn_kp = PlotUtils::flipHisto(h_vn_kp);
  h_vn_km = PlotUtils::flipHisto(h_vn_km);
  h_vn_pr = PlotUtils::flipHisto(h_vn_pr);

  h_vn_pp_ext = PlotUtils::flipHisto(h_vn_pp_ext);
  h_vn_pm_ext = PlotUtils::flipHisto(h_vn_pm_ext);
  h_vn_kp_ext = PlotUtils::flipHisto(h_vn_kp_ext);
  h_vn_km_ext = PlotUtils::flipHisto(h_vn_km_ext);
  h_vn_pr_ext = PlotUtils::flipHisto(h_vn_pr_ext);
  
  h_vn_pr_for = PlotUtils::flipHisto(h_vn_pr_for);

  // Trim and clean up x-axis
  h_vn_EpdE = PlotUtils::trimCentralityPlot(h_vn_EpdE);
  h_vn_EpdF = PlotUtils::trimCentralityPlot(h_vn_EpdF);
  h_vn_TpcB = PlotUtils::trimCentralityPlot(h_vn_TpcB);
  h_vn_Tpc  = PlotUtils::trimCentralityPlot(h_vn_Tpc);

  h_vn_pp = PlotUtils::trimCentralityPlot(h_vn_pp);
  h_vn_pm = PlotUtils::trimCentralityPlot(h_vn_pm);
  h_vn_kp = PlotUtils::trimCentralityPlot(h_vn_kp);
  h_vn_km = PlotUtils::trimCentralityPlot(h_vn_km);
  h_vn_pr = PlotUtils::trimCentralityPlot(h_vn_pr);

  h_vn_pp_ext = PlotUtils::trimCentralityPlot(h_vn_pp_ext);
  h_vn_pm_ext = PlotUtils::trimCentralityPlot(h_vn_pm_ext);
  h_vn_kp_ext = PlotUtils::trimCentralityPlot(h_vn_kp_ext);
  h_vn_km_ext = PlotUtils::trimCentralityPlot(h_vn_km_ext);
  h_vn_pr_ext = PlotUtils::trimCentralityPlot(h_vn_pr_ext);
  
  h_vn_pr_for = PlotUtils::trimCentralityPlot(h_vn_pr_for);

  
  THStack *piCentralityStack = new THStack("piCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kaCentralityStack = new THStack("kaCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");

  THStack *ppExtCentralityStack = new THStack("ppExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *pmExtCentralityStack = new THStack("pmExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kpExtCentralityStack = new THStack("kpExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kmExtCentralityStack = new THStack("kmExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *prExtCentralityStack = new THStack("prExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");

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
  h_vn_pr_ext->GetYaxis()->SetTitleOffset(1.7);

  h_vn_pr_for->SetMarkerStyle(20);
  h_vn_pr_for->SetMarkerSize(2.5);
  h_vn_pr_for->SetMarkerColor(4);
  h_vn_pr_for->SetLineWidth(3);
  h_vn_pr_for->GetYaxis()->SetTitleOffset(1.7);
  
  h_vn_EpdE->SetMarkerStyle(20);
  h_vn_EpdE->SetMarkerSize(2.5);
  h_vn_EpdE->SetMarkerColor(46);
  h_vn_EpdE->SetLineColor(46);
  h_vn_EpdE->SetLineWidth(3);
  h_vn_EpdE->GetYaxis()->SetTitleOffset(1.7);

  h_vn_EpdF->SetMarkerStyle(20);
  h_vn_EpdF->SetMarkerSize(2.5);
  h_vn_EpdF->SetMarkerColor(38);
  h_vn_EpdF->SetLineColor(38);
  h_vn_EpdF->SetLineWidth(3);
  h_vn_EpdF->GetYaxis()->SetTitleOffset(1.7);

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

  prExtCentralityStack->Add(h_vn_pr_for);
  prExtCentralityStack->Add(h_vn_pr);
  prExtCentralityStack->Add(h_vn_pr_ext);

  etaRegionStack->Add(h_vn_EpdE);
  etaRegionStack->Add(h_vn_EpdF);
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

      
      TLegend *etaLegend = new TLegend(0.65, 0.7, 0.9, 0.9);
      etaLegend->AddEntry(h_vn_EpdE, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(h_vn_EpdF, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(h_vn_TpcB, "TPC -1.0 < #eta < 0");

      
      TPaveText *piText = new TPaveText(15, -0.004, 45, 0.008, "NB");
      piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      piText->AddText("0 < y_{CM} < 0.5 GeV");
      piText->AddText("0.18 < p_{T} < 1.6 GeV");
      piText->SetFillColorAlpha(0,0);
      piText->SetLineColorAlpha(0,0);

      TPaveText *kaText = new TPaveText(20, -0.1, 50, -0.07, "NB");
      kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kaText->AddText("0 < y_{CM} < 0.5 GeV");
      kaText->AddText("0.18 < p_{T} < 1.6 GeV");
      kaText->SetFillColorAlpha(0,0);
      kaText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(5, -0.07, 35, -0.05, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 < p_{T} < 2.0 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);

      
      canvas->SetLeftMargin(0.13);
      canvas->SetGridx();
      canvas->SetGridy();
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
     

      h_vn_EpdE->SetMarkerStyle(20);
      h_vn_EpdE->SetMarkerSize(2);
      //h_vn_EpdE->SetMarkerColor(2);
      //h_vn_EpdE->SetLineColor(2);
      h_vn_EpdE->GetXaxis()->SetNdivisions(210);
      h_vn_EpdE->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdE.png");
      canvas->Clear();

      h_vn_EpdF->SetMarkerStyle(20);
      h_vn_EpdF->SetMarkerSize(2);
      //h_vn_EpdF->SetMarkerColor(2);
      //h_vn_EpdF->SetLineColor(2);
      h_vn_EpdF->GetXaxis()->SetNdivisions(210);
      h_vn_EpdF->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdF.png");
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

      TLegend *etaLegend = new TLegend(0.65, 0.25, 0.9, 0.45);
      etaLegend->AddEntry(h_vn_EpdE, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(h_vn_EpdF, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(h_vn_TpcB, "TPC -1.0 < #eta < 0");

      

      TPaveText *piText = new TPaveText(5, 0.025, 38, 0.07, "NB");
      piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      piText->AddText("0 < y_{CM} < 0.5 GeV");
      piText->AddText("0.18 < p_{T} < 1.6 GeV");
      piText->SetFillColorAlpha(0,0);
      piText->SetLineColorAlpha(0,0);

      TPaveText *kaText = new TPaveText(5, 0.025, 38, 0.07, "NB");
      kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kaText->AddText("0 < y_{CM} < 0.5 GeV");
      kaText->AddText("0.18 < p_{T} < 1.6 GeV");
      kaText->SetFillColorAlpha(0,0);
      kaText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(5, 0.025, 38, 0.07, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 < p_{T} < 2.0 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);


      canvas->SetLeftMargin(0.13);
      canvas->SetGridx();
      canvas->SetGridy();
      gStyle->SetErrorX(0);
      gStyle->SetOptStat(0);
      gStyle->SetEndErrorSize(6);


      TLine *zeroLine = new TLine(0, 0, 60, 0);
      zeroLine->SetLineStyle(9);

      Double_t centralityUpperBounds = 0.08;
      Double_t centralityLowerBounds = -0.08;
      
      piCentralityStack->Draw();
      piCentralityStack->GetYaxis()->SetTitleOffset(1.8);
      piCentralityStack->GetXaxis()->SetNdivisions(210);
      piCentralityStack->SetMaximum(centralityUpperBounds);
      piCentralityStack->SetMinimum(centralityLowerBounds);
      piCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      piLegend->Draw();
      piText->Draw();
      canvas->SaveAs(jobID + "_piCentralityStack.png");
      canvas->Clear();

      kaCentralityStack->Draw();
      kaCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kaCentralityStack->GetXaxis()->SetNdivisions(210);
      kaCentralityStack->SetMaximum(centralityUpperBounds);
      kaCentralityStack->SetMinimum(centralityLowerBounds);
      kaCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kaLegend->Draw();
      kaText->Draw();
      canvas->SaveAs(jobID + "_kaCentralityStack.png");
      canvas->Clear();

      h_vn_pr->SetTitle("");
      h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pr->GetXaxis()->SetNdivisions(210);
      h_vn_pr->Draw("E1P");
      h_vn_pr->SetMaximum(centralityUpperBounds);
      h_vn_pr->SetMinimum(centralityLowerBounds);
      zeroLine->Draw("SAME");
      prText->Draw();
      canvas->SaveAs(jobID + "_vn_pr.png");
      canvas->Clear();

      h_vn_EpdE->SetMarkerStyle(20);
      h_vn_EpdE->SetMarkerSize(2);
      //h_vn_EpdE->SetMarkerColor(2);
      //h_vn_EpdE->SetLineColor(2);
      h_vn_EpdE->GetYaxis()->SetTitleOffset(1.7);
      h_vn_EpdE->GetXaxis()->SetNdivisions(210);
      h_vn_EpdE->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdE.png");
      canvas->Clear();

      h_vn_EpdF->SetMarkerStyle(20);
      h_vn_EpdF->SetMarkerSize(2);
      //h_vn_EpdF->SetMarkerColor(2);
      //h_vn_EpdF->SetLineColor(2);
      h_vn_EpdF->GetYaxis()->SetTitleOffset(1.8);
      h_vn_EpdF->GetXaxis()->SetNdivisions(210);
      h_vn_EpdF->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdF.png");
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
      kmExtLegend->Draw();
      //kmExtText->Draw();
      canvas->SaveAs(jobID + "_kmExtCentralityStack.png");
      canvas->Clear();

      prExtCentralityStack->Draw();
      prExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      prExtCentralityStack->GetXaxis()->SetNdivisions(210);
      prExtCentralityStack->SetMaximum(0.04);
      prExtCentralityStack->SetMinimum(-0.08);
      prExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      prExtLegend->Draw();
      //prExtText->Draw();
      canvas->SaveAs(jobID + "_prExtCentralityStack.png");
      canvas->Clear();
      
      etaRegionStack->Draw();
      etaRegionStack->GetYaxis()->SetTitleOffset(1.7);
      etaRegionStack->GetXaxis()->SetNdivisions(210);
      etaRegionStack->Draw();
      etaRegionStack->SetMaximum(0.05);
      //etaRegionStack->SetMinimum(-0.1);
      etaRegionStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
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
