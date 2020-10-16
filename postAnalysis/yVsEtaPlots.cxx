void yVsEtaPlots(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 875, 675);

  TH1D *h_y_eta_asymm_pp = new TH1D("h_y_eta_asymm_pp", "#pi^{+} y-#eta Asymmetry;p_{T} (GeV);<y> - <#eta>", 5, 0, 2.5);//#frac{<y> - <#eta>}{<y> + <#eta>}", 5, 0, 2.5);
  TH1D *h_y_eta_asymm_pm = new TH1D("h_y_eta_asymm_pm", "#pi^{-} y-#eta Asymmetry;p_{T} (GeV);<y> - <#eta>", 5, 0, 2.5);//#frac{<y> - <#eta>}{<y> + <#eta>}", 5, 0, 2.5);
  TH1D *h_y_eta_asymm_kp = new TH1D("h_y_eta_asymm_kp", "K^{+} y-#eta Asymmetry;p_{T} (GeV);<y> - <#eta>", 5, 0, 2.5);//#frac{<y> - <#eta>}{<y> + <#eta>}", 5, 0, 2.5);
  TH1D *h_y_eta_asymm_km = new TH1D("h_y_eta_asymm_km", "K^{-} y-#eta Asymmetry;p_{T} (GeV);<y> - <#eta>", 5, 0, 2.5);//#frac{<y> - <#eta>}{<y> + <#eta>}", 5, 0, 2.5);
  TH1D *h_y_eta_asymm_pr = new TH1D("h_y_eta_asymm_pr", "Proton y-#eta Asymmetry;p_{T} (GeV);<y> - <#eta>", 5, 0, 2.5);//#frac{<y> - <#eta>}{<y> + <#eta>}", 5, 0, 2.5);

  TObjArray *array_pp = new TObjArray();
  array_pp->Add((TH1D*)file->Get("h2_y_vs_eta_pt0p5to1_pp"));
  array_pp->Add((TH1D*)file->Get("h2_y_vs_eta_pt1to1p5_pp"));
  array_pp->Add((TH1D*)file->Get("h2_y_vs_eta_pt1p5to2_pp"));
  array_pp->SetOwner(kTRUE);

  TObjArray *array_pm = new TObjArray();
  array_pm->Add((TH1D*)file->Get("h2_y_vs_eta_pt0p5to1_pm"));
  array_pm->Add((TH1D*)file->Get("h2_y_vs_eta_pt1to1p5_pm"));
  array_pm->Add((TH1D*)file->Get("h2_y_vs_eta_pt1p5to2_pm"));
  array_pm->SetOwner(kTRUE);

  TObjArray *array_kp = new TObjArray();
  array_kp->Add((TH1D*)file->Get("h2_y_vs_eta_pt0p5to1_kp"));
  array_kp->Add((TH1D*)file->Get("h2_y_vs_eta_pt1to1p5_kp"));
  array_kp->Add((TH1D*)file->Get("h2_y_vs_eta_pt1p5to2_kp"));
  array_kp->SetOwner(kTRUE);

  TObjArray *array_km = new TObjArray();
  array_km->Add((TH1D*)file->Get("h2_y_vs_eta_pt0p5to1_km"));
  array_km->Add((TH1D*)file->Get("h2_y_vs_eta_pt1to1p5_km"));
  array_km->Add((TH1D*)file->Get("h2_y_vs_eta_pt1p5to2_km"));
  array_km->SetOwner(kTRUE);

  TObjArray *array_pr = new TObjArray();
  array_pr->Add((TH1D*)file->Get("h2_y_vs_eta_pt0p5to1_pr"));
  array_pr->Add((TH1D*)file->Get("h2_y_vs_eta_pt1to1p5_pr"));
  array_pr->Add((TH1D*)file->Get("h2_y_vs_eta_pt1p5to2_pr"));
  array_pr->SetOwner(kTRUE);

  for (int i = 2; i <= 4; i++)
  {
    Double_t y_mean_pp         = ((TH1D*)array_pp->At(i-2))->GetMean(2);
    Double_t y_mean_error_pp   = ((TH1D*)array_pp->At(i-2))->GetMean(12);
    Double_t eta_mean_pp       = ((TH1D*)array_pp->At(i-2))->GetMean(1);
    Double_t eta_mean_error_pp = ((TH1D*)array_pp->At(i-2))->GetMean(11);
    
    Double_t y_mean_pm         = ((TH1D*)array_pm->At(i-2))->GetMean(2);
    Double_t y_mean_error_pm   = ((TH1D*)array_pm->At(i-2))->GetMean(12);
    Double_t eta_mean_pm       = ((TH1D*)array_pm->At(i-2))->GetMean(1);
    Double_t eta_mean_error_pm = ((TH1D*)array_pm->At(i-2))->GetMean(11);

    Double_t y_mean_kp         = ((TH1D*)array_kp->At(i-2))->GetMean(2);
    Double_t y_mean_error_kp   = ((TH1D*)array_kp->At(i-2))->GetMean(12);
    Double_t eta_mean_kp       = ((TH1D*)array_kp->At(i-2))->GetMean(1);
    Double_t eta_mean_error_kp = ((TH1D*)array_kp->At(i-2))->GetMean(11);
    
    Double_t y_mean_km         = ((TH1D*)array_km->At(i-2))->GetMean(2);
    Double_t y_mean_error_km   = ((TH1D*)array_km->At(i-2))->GetMean(12);
    Double_t eta_mean_km       = ((TH1D*)array_km->At(i-2))->GetMean(1);
    Double_t eta_mean_error_km = ((TH1D*)array_km->At(i-2))->GetMean(11);

    Double_t y_mean_pr         = ((TH1D*)array_pr->At(i-2))->GetMean(2);
    Double_t y_mean_error_pr   = ((TH1D*)array_pr->At(i-2))->GetMean(12);
    Double_t eta_mean_pr       = ((TH1D*)array_pr->At(i-2))->GetMean(1);
    Double_t eta_mean_error_pr = ((TH1D*)array_pr->At(i-2))->GetMean(11);
    /*
    // NORMALIZED CALCULATIONS    
    Double_t asymm_pp = (y_mean_pp - eta_mean_pp) / (y_mean_pp + eta_mean_pp);
    Double_t num_pp = y_mean_pp - eta_mean_pp;
    Double_t den_pp = y_mean_pp + eta_mean_pp;
    Double_t num_error_pp = TMath::Sqrt(y_mean_error_pp*y_mean_error_pp + eta_mean_error_pp*eta_mean_error_pp);
    Double_t den_error_pp = TMath::Sqrt(y_mean_error_pp*y_mean_error_pp + eta_mean_error_pp*eta_mean_error_pp);
    Double_t asymm_error_pp = asymm_pp * TMath::Sqrt( (num_error_pp/num_pp)*(num_error_pp/num_pp) + (den_error_pp/den_pp)*(den_error_pp/den_pp) );

    Double_t asymm_pm = (y_mean_pm - eta_mean_pm) / (y_mean_pm + eta_mean_pm);
    Double_t num_pm = y_mean_pm - eta_mean_pm;
    Double_t den_pm = y_mean_pm + eta_mean_pm;
    Double_t num_error_pm = TMath::Sqrt(y_mean_error_pm*y_mean_error_pm + eta_mean_error_pm*eta_mean_error_pm);
    Double_t den_error_pm = TMath::Sqrt(y_mean_error_pm*y_mean_error_pm + eta_mean_error_pm*eta_mean_error_pm);
    Double_t asymm_error_pm = asymm_pm * TMath::Sqrt( (num_error_pm/num_pm)*(num_error_pm/num_pm) + (den_error_pm/den_pm)*(den_error_pm/den_pm) );

    Double_t asymm_kp = (y_mean_kp - eta_mean_kp) / (y_mean_kp + eta_mean_kp);
    Double_t num_kp = y_mean_kp - eta_mean_kp;
    Double_t den_kp = y_mean_kp + eta_mean_kp;
    Double_t num_error_kp = TMath::Sqrt(y_mean_error_kp*y_mean_error_kp + eta_mean_error_kp*eta_mean_error_kp);
    Double_t den_error_kp = TMath::Sqrt(y_mean_error_kp*y_mean_error_kp + eta_mean_error_kp*eta_mean_error_kp);
    Double_t asymm_error_kp = asymm_kp * TMath::Sqrt( (num_error_kp/num_kp)*(num_error_kp/num_kp) + (den_error_kp/den_kp)*(den_error_kp/den_kp) );

    Double_t asymm_km = (y_mean_km - eta_mean_km) / (y_mean_km + eta_mean_km);
    Double_t num_km = y_mean_km - eta_mean_km;
    Double_t den_km = y_mean_km + eta_mean_km;
    Double_t num_error_km = TMath::Sqrt(y_mean_error_km*y_mean_error_km + eta_mean_error_km*eta_mean_error_km);
    Double_t den_error_km = TMath::Sqrt(y_mean_error_km*y_mean_error_km + eta_mean_error_km*eta_mean_error_km);
    Double_t asymm_error_km = asymm_km * TMath::Sqrt( (num_error_km/num_km)*(num_error_km/num_km) + (den_error_km/den_km)*(den_error_km/den_km) );

    Double_t asymm_pr = (y_mean_pr - eta_mean_pr) / (y_mean_pr + eta_mean_pr);
    Double_t num_pr = y_mean_pr - eta_mean_pr;
    Double_t den_pr = y_mean_pr + eta_mean_pr;
    Double_t num_error_pr = TMath::Sqrt(y_mean_error_pr*y_mean_error_pr + eta_mean_error_pr*eta_mean_error_pr);
    Double_t den_error_pr = TMath::Sqrt(y_mean_error_pr*y_mean_error_pr + eta_mean_error_pr*eta_mean_error_pr);
    Double_t asymm_error_pr = asymm_pr * TMath::Sqrt( (num_error_pr/num_pr)*(num_error_pr/num_pr) + (den_error_pr/den_pr)*(den_error_pr/den_pr) );
    */

    // NON-NORMALIZED CALCULATIONS
    Double_t asymm_pp = (y_mean_pp - eta_mean_pp);
    Double_t asymm_error_pp = TMath::Sqrt(y_mean_error_pp*y_mean_error_pp + eta_mean_error_pp*eta_mean_error_pp);

    Double_t asymm_pm = (y_mean_pm - eta_mean_pm);
    Double_t asymm_error_pm = TMath::Sqrt(y_mean_error_pm*y_mean_error_pm + eta_mean_error_pm*eta_mean_error_pm);

    Double_t asymm_kp = (y_mean_kp - eta_mean_kp);
    Double_t asymm_error_kp = TMath::Sqrt(y_mean_error_kp*y_mean_error_kp + eta_mean_error_kp*eta_mean_error_kp);

    Double_t asymm_km = (y_mean_km - eta_mean_km);
    Double_t asymm_error_km = TMath::Sqrt(y_mean_error_km*y_mean_error_km + eta_mean_error_km*eta_mean_error_km);

    Double_t asymm_pr = (y_mean_pr - eta_mean_pr);
    Double_t asymm_error_pr = TMath::Sqrt(y_mean_error_pr*y_mean_error_pr + eta_mean_error_pr*eta_mean_error_pr);

    h_y_eta_asymm_pp->SetBinContent(i, asymm_pp);
    h_y_eta_asymm_pm->SetBinContent(i, asymm_pm);
    h_y_eta_asymm_kp->SetBinContent(i, asymm_kp);
    h_y_eta_asymm_km->SetBinContent(i, asymm_km);
    h_y_eta_asymm_pr->SetBinContent(i, asymm_pr);

    h_y_eta_asymm_pp->SetBinError(i, asymm_error_pp);
    h_y_eta_asymm_pm->SetBinError(i, asymm_error_pm);
    h_y_eta_asymm_kp->SetBinError(i, asymm_error_kp);
    h_y_eta_asymm_km->SetBinError(i, asymm_error_km);
    h_y_eta_asymm_pr->SetBinError(i, asymm_error_pr);
  }

  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLeftMargin(0.15);
  canvas->cd();

  gStyle->SetOptStat(0);

  h_y_eta_asymm_pp->GetYaxis()->SetTitleOffset(1.8);
  h_y_eta_asymm_pp->Draw("E1");
  canvas->SaveAs("h_y_eta_asymm_pp.png");
  canvas->Clear();

  h_y_eta_asymm_pm->GetYaxis()->SetTitleOffset(1.8);
  h_y_eta_asymm_pm->Draw("E1");
  canvas->SaveAs("h_y_eta_asymm_pm.png");
  canvas->Clear();

  h_y_eta_asymm_kp->GetYaxis()->SetTitleOffset(1.8);
  h_y_eta_asymm_kp->Draw("E1");
  canvas->SaveAs("h_y_eta_asymm_kp.png");
  canvas->Clear();

  h_y_eta_asymm_km->GetYaxis()->SetTitleOffset(1.8);
  h_y_eta_asymm_km->Draw("E1");
  canvas->SaveAs("h_y_eta_asymm_km.png");
  canvas->Clear();

  h_y_eta_asymm_pr->GetYaxis()->SetTitleOffset(1.8);
  h_y_eta_asymm_pr->Draw("E1");
  canvas->SaveAs("h_y_eta_asymm_pr.png");
  canvas->Clear();



  THStack *stack = new THStack("stack","TPC y-#eta Asymmetries;p_{T} (GeV);<y> - <#eta>");//#frac{<y> - <#eta>}{<y> + <#eta>}");
  h_y_eta_asymm_pp->SetMarkerStyle(21);
  h_y_eta_asymm_pm->SetMarkerStyle(21);
  h_y_eta_asymm_kp->SetMarkerStyle(21);
  h_y_eta_asymm_km->SetMarkerStyle(21);
  h_y_eta_asymm_pr->SetMarkerStyle(21);

  h_y_eta_asymm_pp->SetMarkerColor(kGreen);
  h_y_eta_asymm_pm->SetMarkerColor(kCyan);
  h_y_eta_asymm_kp->SetMarkerColor(kBlue);
  h_y_eta_asymm_km->SetMarkerColor(kMagenta);
  h_y_eta_asymm_pr->SetMarkerColor(kRed);

  h_y_eta_asymm_pp->SetLineColor(kGreen);
  h_y_eta_asymm_pm->SetLineColor(kCyan);
  h_y_eta_asymm_kp->SetLineColor(kBlue);
  h_y_eta_asymm_km->SetLineColor(kMagenta);
  h_y_eta_asymm_pr->SetLineColor(kRed);

  stack->Add(h_y_eta_asymm_pp);
  stack->Add(h_y_eta_asymm_pm);
  stack->Add(h_y_eta_asymm_kp);
  stack->Add(h_y_eta_asymm_km);
  stack->Add(h_y_eta_asymm_pr);

  TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85);
  legend->SetHeader("Particles", "C");
  legend->AddEntry(h_y_eta_asymm_pp,"#pi^{+}");
  legend->AddEntry(h_y_eta_asymm_pm,"#pi^{-}");
  legend->AddEntry(h_y_eta_asymm_kp,"K^{+}");
  legend->AddEntry(h_y_eta_asymm_km,"K^{-}");
  legend->AddEntry(h_y_eta_asymm_pr,"Protons");
  
  TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 875, 675);
  canvas2->SetGridx();
  canvas2->SetGridy();
  canvas2->cd();
  canvas2->SetLeftMargin(0.13);
  stack->Draw();
  stack->GetYaxis()->SetTitleOffset(1.3);
  canvas2->Clear();
  stack->Draw("NOSTACK E1P");
  legend->Draw();
  canvas2->SaveAs("all_asymmetries.png");
  canvas2->Clear();
  delete canvas2;

  canvas->Clear();
  file->Close();
}
