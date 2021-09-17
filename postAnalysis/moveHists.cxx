#include "PlotUtils.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"

void moveHists(TString fileName = "Normal.picoDst.result.combined.root")
{
  TFile *file = TFile::Open(fileName);

  TProfile *p_vn_Tpc_pT_0p2to2 = (TProfile*)file->Get("p_vn_Tpc_pT_0p2to2");
    
  TProfile *p_vn_pp = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr = (TProfile*)file->Get("p_vn_pr");

  TProfile2D *p2_vn_yCM_cent_pp = (TProfile2D*)file->Get("p2_vn_yCM_cent_pp");
  TProfile2D *p2_vn_yCM_cent_pm = (TProfile2D*)file->Get("p2_vn_yCM_cent_pm");
  TProfile2D *p2_vn_yCM_cent_kp = (TProfile2D*)file->Get("p2_vn_yCM_cent_kp");
  TProfile2D *p2_vn_yCM_cent_km = (TProfile2D*)file->Get("p2_vn_yCM_cent_km");
  TProfile2D *p2_vn_yCM_cent_pr = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr");
  TProfile2D *p2_vn_yCM_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr_symmetry");

  p2_vn_yCM_cent_kp->RebinY();
  p2_vn_yCM_cent_km->RebinY();

  TProfile2D *p2_vn_pT_cent_pp = (TProfile2D*)file->Get("p2_vn_pT_cent_pp");
  TProfile2D *p2_vn_pT_cent_pm = (TProfile2D*)file->Get("p2_vn_pT_cent_pm");
  TProfile2D *p2_vn_pT_cent_kp = (TProfile2D*)file->Get("p2_vn_pT_cent_kp");
  TProfile2D *p2_vn_pT_cent_km = (TProfile2D*)file->Get("p2_vn_pT_cent_km");
  TProfile2D *p2_vn_pT_cent_pr = (TProfile2D*)file->Get("p2_vn_pT_cent_pr");

  p2_vn_pT_cent_kp->RebinY();
  p2_vn_pT_cent_km->RebinY();

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



  TH1D *h_vn_Tpc_pT_0p2to2 = p_vn_Tpc_pT_0p2to2->ProjectionX();
  TH1D *h_vn_pp = p_vn_pp->ProjectionX();
  TH1D *h_vn_pm = p_vn_pm->ProjectionX();
  TH1D *h_vn_kp = p_vn_kp->ProjectionX();
  TH1D *h_vn_km = p_vn_km->ProjectionX();
  TH1D *h_vn_pr = p_vn_pr->ProjectionX();

  TH1D *h_vn_yCM_00to10_pp = p_vn_yCM_00to10_pp->ProjectionX();
  TH1D *h_vn_yCM_10to40_pp = p_vn_yCM_10to40_pp->ProjectionX();
  TH1D *h_vn_yCM_40to60_pp = p_vn_yCM_40to60_pp->ProjectionX();
  
  TH1D *h_vn_yCM_00to10_pm = p_vn_yCM_00to10_pm->ProjectionX();
  TH1D *h_vn_yCM_10to40_pm = p_vn_yCM_10to40_pm->ProjectionX();
  TH1D *h_vn_yCM_40to60_pm = p_vn_yCM_40to60_pm->ProjectionX();

  TH1D *h_vn_yCM_00to10_kp = p_vn_yCM_00to10_kp->ProjectionX();
  TH1D *h_vn_yCM_10to40_kp = p_vn_yCM_10to40_kp->ProjectionX();
  TH1D *h_vn_yCM_40to60_kp = p_vn_yCM_40to60_kp->ProjectionX();

  TH1D *h_vn_yCM_00to10_km = p_vn_yCM_00to10_km->ProjectionX();
  TH1D *h_vn_yCM_10to40_km = p_vn_yCM_10to40_km->ProjectionX();
  TH1D *h_vn_yCM_40to60_km = p_vn_yCM_40to60_km->ProjectionX();

  TH1D *h_vn_yCM_00to10_pr = p_vn_yCM_00to10_pr->ProjectionX();
  TH1D *h_vn_yCM_10to40_pr = p_vn_yCM_10to40_pr->ProjectionX();
  TH1D *h_vn_yCM_40to60_pr = p_vn_yCM_40to60_pr->ProjectionX();

  TH1D *h_vn_yCM_00to10_pr_symm = p_vn_yCM_00to10_pr_symm->ProjectionX();
  TH1D *h_vn_yCM_10to40_pr_symm = p_vn_yCM_10to40_pr_symm->ProjectionX();
  TH1D *h_vn_yCM_40to60_pr_symm = p_vn_yCM_40to60_pr_symm->ProjectionX();


  TH1D *h_vn_pT_00to10_pp = p_vn_pT_00to10_pp->ProjectionX();
  TH1D *h_vn_pT_10to40_pp = p_vn_pT_10to40_pp->ProjectionX();
  TH1D *h_vn_pT_40to60_pp = p_vn_pT_40to60_pp->ProjectionX();
  
  TH1D *h_vn_pT_00to10_pm = p_vn_pT_00to10_pm->ProjectionX();
  TH1D *h_vn_pT_10to40_pm = p_vn_pT_10to40_pm->ProjectionX();
  TH1D *h_vn_pT_40to60_pm = p_vn_pT_40to60_pm->ProjectionX();

  TH1D *h_vn_pT_00to10_kp = p_vn_pT_00to10_kp->ProjectionX();
  TH1D *h_vn_pT_10to40_kp = p_vn_pT_10to40_kp->ProjectionX();
  TH1D *h_vn_pT_40to60_kp = p_vn_pT_40to60_kp->ProjectionX();

  TH1D *h_vn_pT_00to10_km = p_vn_pT_00to10_km->ProjectionX();
  TH1D *h_vn_pT_10to40_km = p_vn_pT_10to40_km->ProjectionX();
  TH1D *h_vn_pT_40to60_km = p_vn_pT_40to60_km->ProjectionX();

  TH1D *h_vn_pT_00to10_pr = p_vn_pT_00to10_pr->ProjectionX();
  TH1D *h_vn_pT_10to40_pr = p_vn_pT_10to40_pr->ProjectionX();
  TH1D *h_vn_pT_40to60_pr = p_vn_pT_40to60_pr->ProjectionX();


  h_vn_Tpc_pT_0p2to2 = PlotUtils::flipHisto(h_vn_Tpc_pT_0p2to2);
  
  h_vn_pp = PlotUtils::flipHisto(h_vn_pp);
  h_vn_pm = PlotUtils::flipHisto(h_vn_pm);
  h_vn_kp = PlotUtils::flipHisto(h_vn_kp);
  h_vn_km = PlotUtils::flipHisto(h_vn_km);
  h_vn_pr = PlotUtils::flipHisto(h_vn_pr);

  h_vn_Tpc_pT_0p2to2 = PlotUtils::trimCentralityPlot(h_vn_Tpc_pT_0p2to2);

  h_vn_pp = PlotUtils::trimCentralityPlot(h_vn_pp);
  h_vn_pm = PlotUtils::trimCentralityPlot(h_vn_pm);
  h_vn_kp = PlotUtils::trimCentralityPlot(h_vn_kp);
  h_vn_km = PlotUtils::trimCentralityPlot(h_vn_km);
  h_vn_pr = PlotUtils::trimCentralityPlot(h_vn_pr);
  
  h_vn_yCM_00to10_pp = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_pp);
  h_vn_yCM_10to40_pp = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_pp);
  h_vn_yCM_40to60_pp = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_pp);
  
  h_vn_yCM_00to10_pm = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_pm);
  h_vn_yCM_10to40_pm = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_pm);
  h_vn_yCM_40to60_pm = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_pm);

  h_vn_yCM_00to10_kp = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_kp);
  h_vn_yCM_10to40_kp = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_kp);
  h_vn_yCM_40to60_kp = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_kp);

  h_vn_yCM_00to10_km = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_km);
  h_vn_yCM_10to40_km = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_km);
  h_vn_yCM_40to60_km = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_km);


  TFile *newFile = new TFile("v3_results.root", "RECREATE");
  newFile->cd();

  h_vn_Tpc_pT_0p2to2->Write();

  h_vn_pp->Write();
  h_vn_pm->Write();
  h_vn_kp->Write();
  h_vn_km->Write();
  h_vn_pr->Write();

  h_vn_yCM_00to10_pp->Write();
  h_vn_yCM_10to40_pp->Write();
  h_vn_yCM_40to60_pp->Write();

  h_vn_yCM_00to10_pm->Write();
  h_vn_yCM_10to40_pm->Write();
  h_vn_yCM_40to60_pm->Write();
  
  h_vn_yCM_00to10_kp->Write();
  h_vn_yCM_10to40_kp->Write();
  h_vn_yCM_40to60_kp->Write();
  
  h_vn_yCM_00to10_km->Write();
  h_vn_yCM_10to40_km->Write();
  h_vn_yCM_40to60_km->Write();
  
  h_vn_yCM_00to10_pr->Write();
  h_vn_yCM_10to40_pr->Write();
  h_vn_yCM_40to60_pr->Write();

  h_vn_yCM_00to10_pr_symm->Write();
  h_vn_yCM_10to40_pr_symm->Write();
  h_vn_yCM_40to60_pr_symm->Write();

  h_vn_pT_00to10_pp->Write();
  h_vn_pT_10to40_pp->Write();
  h_vn_pT_40to60_pp->Write();

  h_vn_pT_00to10_pm->Write();
  h_vn_pT_10to40_pm->Write();
  h_vn_pT_40to60_pm->Write();
  
  h_vn_pT_00to10_kp->Write();
  h_vn_pT_10to40_kp->Write();
  h_vn_pT_40to60_kp->Write();
  
  h_vn_pT_00to10_km->Write();
  h_vn_pT_10to40_km->Write();
  h_vn_pT_40to60_km->Write();
  
  h_vn_pT_00to10_pr->Write();
  h_vn_pT_10to40_pr->Write();
  h_vn_pT_40to60_pr->Write();


  file->Close();
  newFile->Close();
}
