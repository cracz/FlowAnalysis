void acceptanceCuts(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TProfile2D *p2_pp_vs_eta = (TProfile2D*)file->Get("p2_pp_vs_eta");
  TH2D *h2_phi_vs_eta_EPD = (TH2D*)file->Get("h2_phi_vs_eta_EPD");
  
  TH2D *h2_pT_vs_yCM_pp = (TH2D*)file->Get("h2_pT_vs_yCM_pp");
  TH2D *h2_pT_vs_yCM_pm = (TH2D*)file->Get("h2_pT_vs_yCM_pm");
  TH2D *h2_pT_vs_yCM_kp = (TH2D*)file->Get("h2_pT_vs_yCM_kp");
  TH2D *h2_pT_vs_yCM_km = (TH2D*)file->Get("h2_pT_vs_yCM_km");
  TH2D *h2_pT_vs_yCM_pr = (TH2D*)file->Get("h2_pT_vs_yCM_pr");
  h2_pT_vs_yCM_pp->SetTitle("");
  h2_pT_vs_yCM_pm->SetTitle("");
  h2_pT_vs_yCM_kp->SetTitle("");
  h2_pT_vs_yCM_km->SetTitle("");
  h2_pT_vs_yCM_pr->SetTitle("");

  Double_t yCM_low_pp  = 0.0;
  Double_t yCM_high_pp = 0.545;
  Double_t pT_low_pp   = 0.18;
  Double_t pT_high_pp  = 1.6;

  Double_t yCM_low_pm  = 0.0;
  Double_t yCM_high_pm = 0.545;
  Double_t pT_low_pm   = 0.18;
  Double_t pT_high_pm  = 1.6;

  Double_t yCM_low_kp  = 0.0;
  Double_t yCM_high_kp = 0.545;
  Double_t pT_low_kp   = 0.4;
  Double_t pT_high_kp  = 1.6;

  Double_t yCM_low_km  = 0.0;
  Double_t yCM_high_km = 0.545;
  Double_t pT_low_km   = 0.4;
  Double_t pT_high_km  = 1.6;

  Double_t yCM_low_pr  = 0.0;
  Double_t yCM_high_pr = 0.545;
  Double_t pT_low_pr   = 0.4;
  Double_t pT_high_pr  = 2.0;

  TLine *y_mid = new TLine(0, 0, 0, 2.5);
  y_mid->SetLineColor(kRed);
  y_mid->SetLineWidth(4);

  TLine *y_mid_pr = new TLine(0, 0, 0, 3);
  y_mid_pr->SetLineColor(kRed);
  y_mid_pr->SetLineWidth(4);

  TLine *left_pp = new TLine(yCM_low_pp, pT_low_pp, yCM_low_pp, pT_high_pp);
  TLine *right_pp = new TLine(yCM_high_pp, pT_low_pp, yCM_high_pp, pT_high_pp);
  TLine *top_pp = new TLine(yCM_low_pp, pT_high_pp, yCM_high_pp, pT_high_pp);
  TLine *bottom_pp = new TLine(yCM_low_pp, pT_low_pp, yCM_high_pp, pT_low_pp);
  left_pp->SetLineWidth(4);
  right_pp->SetLineWidth(4);
  top_pp->SetLineWidth(4);
  bottom_pp->SetLineWidth(4);
  
  TLine *left_pm = new TLine(yCM_low_pm, pT_low_pm, yCM_low_pm, pT_high_pm);
  TLine *right_pm = new TLine(yCM_high_pm, pT_low_pm, yCM_high_pm, pT_high_pm);
  TLine *top_pm = new TLine(yCM_low_pm, pT_high_pm, yCM_high_pm, pT_high_pm);
  TLine *bottom_pm = new TLine(yCM_low_pm, pT_low_pm, yCM_high_pm, pT_low_pm);
  left_pm->SetLineWidth(4);
  right_pm->SetLineWidth(4);
  top_pm->SetLineWidth(4);
  bottom_pm->SetLineWidth(4);

  TLine *left_kp = new TLine(yCM_low_kp, pT_low_kp, yCM_low_kp, pT_high_kp);
  TLine *right_kp = new TLine(yCM_high_kp, pT_low_kp, yCM_high_kp, pT_high_kp);
  TLine *top_kp = new TLine(yCM_low_kp, pT_high_kp, yCM_high_kp, pT_high_kp);
  TLine *bottom_kp = new TLine(yCM_low_kp, pT_low_kp, yCM_high_kp, pT_low_kp);
  left_kp->SetLineWidth(4);
  right_kp->SetLineWidth(4);
  top_kp->SetLineWidth(4);
  bottom_kp->SetLineWidth(4);

  TLine *left_km = new TLine(yCM_low_km, pT_low_km, yCM_low_km, pT_high_km);
  TLine *right_km = new TLine(yCM_high_km, pT_low_km, yCM_high_km, pT_high_km);
  TLine *top_km = new TLine(yCM_low_km, pT_high_km, yCM_high_km, pT_high_km);
  TLine *bottom_km = new TLine(yCM_low_km, pT_low_km, yCM_high_km, pT_low_km);
  left_km->SetLineWidth(4);
  right_km->SetLineWidth(4);
  top_km->SetLineWidth(4);
  bottom_km->SetLineWidth(4);

  TLine *left_pr = new TLine(yCM_low_pr, pT_low_pr, yCM_low_pr, pT_high_pr);
  TLine *right_pr = new TLine(yCM_high_pr, pT_low_pr, yCM_high_pr, pT_high_pr);
  TLine *top_pr = new TLine(yCM_low_pr, pT_high_pr, yCM_high_pr, pT_high_pr);
  TLine *bottom_pr = new TLine(yCM_low_pr, pT_low_pr, yCM_high_pr, pT_low_pr);
  left_pr->SetLineWidth(4);
  right_pr->SetLineWidth(4);
  top_pr->SetLineWidth(4);
  bottom_pr->SetLineWidth(4);

  TLine *pp_vs_eta_cutoff = new TLine(-5.1, 0.5, -5.1, 12.5);
  pp_vs_eta_cutoff->SetLineWidth(3);
  pp_vs_eta_cutoff->SetLineColor(kRed);

  TLine *phi_vs_eta_cutoff = new TLine(-5.1, -4.0, -5.1, 4.0);
  phi_vs_eta_cutoff->SetLineWidth(3);
  phi_vs_eta_cutoff->SetLineColor(kRed);

  TCanvas *canvas = new TCanvas("canvas", "canvas", 875, 675);
  canvas->SetRightMargin(0.12);
  canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLogz();
  canvas->cd();

  gStyle->SetOptStat(0);

  TPaveText *text_pp = new TPaveText(-0.9, 2.1, -0.6, 2.4);
  text_pp->SetFillColorAlpha(0, 0);
  text_pp->AddText("#pi^{+}");

  TPaveText *text_pm = new TPaveText(-0.9, 2.1, -0.6, 2.4);
  text_pm->SetFillColorAlpha(0, 0);
  text_pm->AddText("#pi^{-}");

  TPaveText *text_kp = new TPaveText(-0.9, 2.1, -0.6, 2.4);
  text_kp->SetFillColorAlpha(0, 0);
  text_kp->AddText("K^{+}");

  TPaveText *text_km = new TPaveText(-0.9, 2.1, -0.6, 2.4);
  text_km->SetFillColorAlpha(0, 0);
  text_km->AddText("K^{-}");

  TPaveText *text_pr = new TPaveText(-0.9, 2.1, -0.6, 2.4);
  text_pr->SetFillColorAlpha(0, 0);
  text_pr->AddText("p");
  
  h2_pT_vs_yCM_pp->Draw("colz");
  y_mid->Draw("SAME");
  left_pp->Draw("SAME");
  right_pp->Draw("SAME");
  top_pp->Draw("SAME");
  bottom_pp->Draw("SAME");
  text_pp->Draw("SAME");
  canvas->SaveAs("Acceptance_pp.png");
  canvas->Clear();

  h2_pT_vs_yCM_pm->Draw("colz");
  y_mid->Draw("SAME");
  left_pm->Draw("SAME");
  right_pm->Draw("SAME");
  top_pm->Draw("SAME");
  bottom_pm->Draw("SAME");
  text_pm->Draw("SAME");
  canvas->SaveAs("Acceptance_pm.png");
  canvas->Clear();

  h2_pT_vs_yCM_kp->Draw("colz");
  y_mid->Draw("SAME");
  left_kp->Draw("SAME");
  right_kp->Draw("SAME");
  top_kp->Draw("SAME");
  bottom_kp->Draw("SAME");
  text_kp->Draw("SAME");
  canvas->SaveAs("Acceptance_kp.png");
  canvas->Clear();

  h2_pT_vs_yCM_km->Draw("colz");
  y_mid->Draw("SAME");
  left_km->Draw("SAME");
  right_km->Draw("SAME");
  top_km->Draw("SAME");
  bottom_km->Draw("SAME");
  text_km->Draw("SAME");
  canvas->SaveAs("Acceptance_km.png");
  canvas->Clear();

  h2_pT_vs_yCM_pr->Draw("colz");
  y_mid_pr->Draw("SAME");
  left_pr->Draw("SAME");
  right_pr->Draw("SAME");
  top_pr->Draw("SAME");
  bottom_pr->Draw("SAME");
  text_pr->Draw("SAME");
  canvas->SaveAs("Acceptance_pr.png");
  canvas->Clear();

  canvas->SetLogz(0);
  p2_pp_vs_eta->Draw("colz");
  pp_vs_eta_cutoff->Draw("SAME");
  canvas->SaveAs("epdAcceptance_TnMIP.png");
  canvas->Clear();

  canvas->SetLogz();
  h2_phi_vs_eta_EPD->Draw("colz");
  phi_vs_eta_cutoff->Draw("SAME");
  canvas->SaveAs("epdAcceptance_phi.png");
  canvas->Clear();

  delete canvas;
  file->Close();
}
