void eventPlanes(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 2400, 900);
  //canvas->SetGridx();
  //canvas->SetGridy();
  //canvas->SetLeftMargin(0.15);
  //canvas->SetTopMargin(0.15);
  canvas->Divide(3,1);


  TH1D *h_psiEpdE_RAW = (TH1D*)file->Get("h_psiEpdE_RAW");
  TH1D *h_psiEpdF_RAW = (TH1D*)file->Get("h_psiEpdF_RAW");
  TH1D *h_psiTpcB_RAW = (TH1D*)file->Get("h_psiTpcB_RAW");
  h_psiEpdE_RAW->SetLineWidth(3);
  h_psiEpdF_RAW->SetLineWidth(3);
  h_psiTpcB_RAW->SetLineWidth(3);

  TH1D *h_psiEpdE_RC = (TH1D*)file->Get("h_psiEpdE_RC");
  TH1D *h_psiEpdF_RC = (TH1D*)file->Get("h_psiEpdF_RC");
  TH1D *h_psiTpcB_RC = (TH1D*)file->Get("h_psiTpcB_RC");
  h_psiEpdE_RC->SetLineColor(kBlue);
  h_psiEpdF_RC->SetLineColor(kBlue);
  h_psiTpcB_RC->SetLineColor(kBlue);
  h_psiEpdE_RC->SetLineWidth(3);
  h_psiEpdF_RC->SetLineWidth(3);
  h_psiTpcB_RC->SetLineWidth(3);

  
  TH1D *h_psiEpdE_FLAT = (TH1D*)file->Get("h_psiEpdE_FLAT");
  TH1D *h_psiEpdF_FLAT = (TH1D*)file->Get("h_psiEpdF_FLAT");
  TH1D *h_psiTpcB_FLAT = (TH1D*)file->Get("h_psiTpcB_FLAT");
  h_psiEpdE_FLAT->SetLineColor(kRed);
  h_psiEpdF_FLAT->SetLineColor(kRed);
  h_psiTpcB_FLAT->SetLineColor(kRed);
  h_psiEpdE_FLAT->SetLineWidth(3);
  h_psiEpdF_FLAT->SetLineWidth(3);
  h_psiTpcB_FLAT->SetLineWidth(3);


  THStack *stackEpdE = new THStack("stackEpdE", "Epd E;#psi;");
  stackEpdE->Add(h_psiEpdE_RAW);
  stackEpdE->Add(h_psiEpdE_RC);
  stackEpdE->Add(h_psiEpdE_FLAT);
  
  THStack *stackEpdF = new THStack("stackEpdF", "Epd F;#psi;");
  stackEpdF->Add(h_psiEpdF_RAW);
  stackEpdF->Add(h_psiEpdF_RC);
  stackEpdF->Add(h_psiEpdF_FLAT);

  THStack *stackTpcB = new THStack("stackTpcB", "TPC B;#psi;");
  stackTpcB->Add(h_psiTpcB_RAW);
  stackTpcB->Add(h_psiTpcB_RC);
  stackTpcB->Add(h_psiTpcB_FLAT);

  TLegend *legend = new TLegend(0.2, 0.65, 0.9, 0.9);
  legend->AddEntry(h_psiEpdE_RAW, "Raw");
  legend->AddEntry(h_psiEpdE_RC, "Recentered");
  legend->AddEntry(h_psiEpdE_FLAT, "Recentered and shifted");
  legend->SetFillColorAlpha(0,0);
  legend->SetLineColorAlpha(0,0);
  
  canvas->cd(1);
  gPad->SetTicky();
  gPad->SetRightMargin(0);
  stackEpdE->Draw();
  stackEpdE->SetMinimum(5e5);
  stackEpdE->SetMaximum(1e6);
  stackEpdE->Draw("NOSTACK");
  legend->Draw();
  
  canvas->cd(2);
  gPad->SetTicky();
  gPad->SetLeftMargin(0);
  gPad->SetRightMargin(0);
  stackEpdF->Draw();
  stackEpdF->SetMinimum(5e5);
  stackEpdF->SetMaximum(1e6);
  stackEpdF->Draw("NOSTACK");

  canvas->cd(3);  
  gPad->SetTicky();
  gPad->SetLeftMargin(0);
  stackTpcB->Draw();
  stackTpcB->SetMinimum(5e5);
  stackTpcB->SetMaximum(1e6);
  stackTpcB->Draw("NOSTACK");

  canvas->SaveAs("psiCombined.png");
  file->Close();
}
