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


  TH1D *h_psiEpdA_RAW = (TH1D*)file->Get("h_psiEpdA_RAW");
  TH1D *h_psiEpdB_RAW = (TH1D*)file->Get("h_psiEpdB_RAW");
  TH1D *h_psiTpcB_RAW = (TH1D*)file->Get("h_psiTpcB_RAW");
  h_psiEpdA_RAW->SetLineWidth(3);
  h_psiEpdB_RAW->SetLineWidth(3);
  h_psiTpcB_RAW->SetLineWidth(3);

  TH1D *h_psiEpdA_RC = (TH1D*)file->Get("h_psiEpdA_RC");
  TH1D *h_psiEpdB_RC = (TH1D*)file->Get("h_psiEpdB_RC");
  TH1D *h_psiTpcB_RC = (TH1D*)file->Get("h_psiTpcB_RC");
  h_psiEpdA_RC->SetLineColor(kBlue);
  h_psiEpdB_RC->SetLineColor(kBlue);
  h_psiTpcB_RC->SetLineColor(kBlue);
  h_psiEpdA_RC->SetLineWidth(3);
  h_psiEpdB_RC->SetLineWidth(3);
  h_psiTpcB_RC->SetLineWidth(3);

  
  TH1D *h_psiEpdA_FLAT = (TH1D*)file->Get("h_psiEpdA_FLAT");
  TH1D *h_psiEpdB_FLAT = (TH1D*)file->Get("h_psiEpdB_FLAT");
  TH1D *h_psiTpcB_FLAT = (TH1D*)file->Get("h_psiTpcB_FLAT");
  h_psiEpdA_FLAT->SetLineColor(kRed);
  h_psiEpdB_FLAT->SetLineColor(kRed);
  h_psiTpcB_FLAT->SetLineColor(kRed);
  h_psiEpdA_FLAT->SetLineWidth(3);
  h_psiEpdB_FLAT->SetLineWidth(3);
  h_psiTpcB_FLAT->SetLineWidth(3);


  THStack *stackEpdA = new THStack("stackEpdA", "Epd E;#psi;");
  stackEpdA->Add(h_psiEpdA_RAW);
  stackEpdA->Add(h_psiEpdA_RC);
  stackEpdA->Add(h_psiEpdA_FLAT);
  
  THStack *stackEpdB = new THStack("stackEpdB", "Epd F;#psi;");
  stackEpdB->Add(h_psiEpdB_RAW);
  stackEpdB->Add(h_psiEpdB_RC);
  stackEpdB->Add(h_psiEpdB_FLAT);

  THStack *stackTpcB = new THStack("stackTpcB", "TPC B;#psi;");
  stackTpcB->Add(h_psiTpcB_RAW);
  stackTpcB->Add(h_psiTpcB_RC);
  stackTpcB->Add(h_psiTpcB_FLAT);

  TLegend *legend = new TLegend(0.15, 0.7, 0.8, 0.9);
  legend->AddEntry(h_psiEpdA_RAW, "Raw");
  legend->AddEntry(h_psiEpdA_RC, "Recentered");
  legend->AddEntry(h_psiEpdA_FLAT, "Recentered and shifted");
  legend->SetFillColorAlpha(0,0);
  legend->SetLineColorAlpha(0,0);
  
  canvas->cd(1);
  gPad->SetTicky();
  gPad->SetRightMargin(0);
  stackEpdA->Draw();
  stackEpdA->SetMinimum(3e5);
  stackEpdA->SetMaximum(1e6);
  stackEpdA->Draw("NOSTACK");
  legend->Draw();
  
  canvas->cd(2);
  gPad->SetTicky();
  gPad->SetLeftMargin(0);
  gPad->SetRightMargin(0);
  stackEpdB->Draw();
  stackEpdB->SetMinimum(3e5);
  stackEpdB->SetMaximum(1e6);
  stackEpdB->Draw("NOSTACK");

  canvas->cd(3);  
  gPad->SetTicky();
  gPad->SetLeftMargin(0);
  stackTpcB->Draw();
  stackTpcB->SetMinimum(3e5);
  stackTpcB->SetMaximum(1e6);
  stackTpcB->Draw("NOSTACK");

  canvas->SaveAs(jobID+"_psiCombined.png");
  file->Close();
}
