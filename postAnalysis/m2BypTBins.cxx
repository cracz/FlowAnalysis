void m2BypTBins(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 875, 675);
  canvas->SetGrid();
  canvas->SetLogy();
  gStyle->SetOptStat(0);

  Double_t low_pT_values[7]  = {1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2};
  Double_t high_pT_values[7] = {1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

  TH2D *h2_m2_vs_qpT = (TH2D*)file->Get("h2_m2_vs_qpT");

  for(int i = 0; i < 7; i++)
    {
      Int_t low_pT_bin = h2_m2_vs_qpT->GetXaxis()->FindBin(low_pT_values[i]);
      Int_t high_pT_bin = h2_m2_vs_qpT->GetXaxis()->FindBin(high_pT_values[i]);

      Double_t low_pT  = h2_m2_vs_qpT->GetXaxis()->GetBinLowEdge(low_pT_bin);
      Double_t high_pT = h2_m2_vs_qpT->GetXaxis()->GetBinLowEdge(high_pT_bin);

      TString low_pT_str;
      TString high_pT_str;
      low_pT_str.Form("%1.3f", low_pT);
      high_pT_str.Form("%1.3f", high_pT);

      TH1D *h_m2 = h2_m2_vs_qpT->ProjectionY("h_m2", low_pT_bin, high_pT_bin);
      h_m2->GetYaxis()->SetTitle("Tracks");

      TPaveText *text = new TPaveText(0.2, 0.6, 0.4, 0.9, "NDC");
      text->SetFillColorAlpha(0,0);
      text->AddText(low_pT_str + " < p_{T} < " + high_pT_str);
  
      h_m2->Draw();
      text->Draw("SAME");
      canvas->SaveAs("h_m2_for_pT_"+low_pT_str+"_to_"+high_pT_str+".png");
      canvas->Clear();
    }

  file->Close();
}
