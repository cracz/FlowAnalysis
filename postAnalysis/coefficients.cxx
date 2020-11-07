void coefficients(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 800);
  canvas->SetGridx();
  canvas->SetGridy();
  //canvas->SetLeftMargin(0.15);
  canvas->cd();

  TProfile *p_vn_EpdE = (TProfile*)file->Get("p_vn_EpdE");
  TProfile *p_vn_EpdF = (TProfile*)file->Get("p_vn_EpdF");
  TProfile *p_vn_TpcB = (TProfile*)file->Get("p_vn_TpcB");
  TProfile *p_vn_pp = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr = (TProfile*)file->Get("p_vn_pr");

  TH1D *h_vn_EpdE = new TH1D("h_vn_EpdE", ";Centrality;v_{2}", 16, 0, 16);
  TH1D *h_vn_EpdF = new TH1D("h_vn_EpdF", ";Centrality;v_{2}", 16, 0, 16);
  TH1D *h_vn_TpcB = new TH1D("h_vn_TpcB", ";Centrality;v_{2}", 16, 0, 16);
  TH1D *h_vn_pp = new TH1D("h_vn_pp", ";Centrality;v_{2}", 16, 0, 16);
  TH1D *h_vn_pm = new TH1D("h_vn_pm", ";Centrality;v_{2}", 16, 0, 16);
  TH1D *h_vn_kp = new TH1D("h_vn_kp", ";Centrality;v_{2}", 16, 0, 16);
  TH1D *h_vn_km = new TH1D("h_vn_km", ";Centrality;v_{2}", 16, 0, 16);
  TH1D *h_vn_pr = new TH1D("h_vn_pr", ";Centrality;v_{2}", 16, 0, 16);


  for (int i = 1; i <= 16; i++)
    {
      h_vn_EpdE->SetBinContent(i, p_vn_EpdE->GetBinContent(i));
      h_vn_EpdE->SetBinError(i, p_vn_EpdE->GetBinError(i));

      h_vn_EpdF->SetBinContent(i, p_vn_EpdF->GetBinContent(i));
      h_vn_EpdF->SetBinError(i, p_vn_EpdF->GetBinError(i));

      h_vn_TpcB->SetBinContent(i, p_vn_TpcB->GetBinContent(i));
      h_vn_TpcB->SetBinError(i, p_vn_TpcB->GetBinError(i));

      h_vn_pp->SetBinContent(i, p_vn_pp->GetBinContent(i));
      h_vn_pp->SetBinError(i, p_vn_pp->GetBinError(i));

      h_vn_pm->SetBinContent(i, p_vn_pm->GetBinContent(i));
      h_vn_pm->SetBinError(i, p_vn_pm->GetBinError(i));

      h_vn_kp->SetBinContent(i, p_vn_kp->GetBinContent(i));
      h_vn_kp->SetBinError(i, p_vn_kp->GetBinError(i));

      h_vn_km->SetBinContent(i, p_vn_km->GetBinContent(i));
      h_vn_km->SetBinError(i, p_vn_km->GetBinError(i));

      h_vn_pr->SetBinContent(i, p_vn_pr->GetBinContent(i));
      h_vn_pr->SetBinError(i, p_vn_pr->GetBinError(i));
    }


  TFile *resolutionInfo_INPUT = TFile::Open("resolutionInfo_INPUT.root", "READ");
  if(!resolutionInfo_INPUT) { cout << "No resolution file found!" << endl; return; }

  TH1D *h_resolutions = (TH1D*)resolutionInfo_INPUT->Get("h_resolutions");

  h_vn_EpdE->Divide(h_vn_EpdE, h_resolutions);
  h_vn_EpdF->Divide(h_vn_EpdF, h_resolutions);
  h_vn_TpcB->Divide(h_vn_TpcB, h_resolutions);
  h_vn_pp->Divide(h_vn_pp, h_resolutions);
  h_vn_pm->Divide(h_vn_pm, h_resolutions);
  h_vn_kp->Divide(h_vn_kp, h_resolutions);
  h_vn_km->Divide(h_vn_km, h_resolutions);
  h_vn_pr->Divide(h_vn_pr, h_resolutions);
  

  Int_t centBins    = h_vn_TpcB->GetNbinsX();
  Int_t firstCentID = h_vn_TpcB->GetBinLowEdge(1);
  Int_t lastCentID  = h_vn_TpcB->GetBinLowEdge(h_vn_TpcB->GetNbinsX());

  TH1D *h_vn_EpdE_flip = new TH1D("h_vn_EpdE_flip",h_vn_EpdE->GetTitle(),centBins,0,centBins);
  h_vn_EpdE_flip->GetXaxis()->SetTitle((TString)h_vn_EpdE->GetXaxis()->GetTitle()+" (%)");
  h_vn_EpdE_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_EpdE->GetYaxis()->GetTitle());

  TH1D *h_vn_EpdF_flip = new TH1D("h_vn_EpdF_flip",h_vn_EpdF->GetTitle(),centBins,0,centBins);
  h_vn_EpdF_flip->GetXaxis()->SetTitle((TString)h_vn_EpdF->GetXaxis()->GetTitle()+" (%)");
  h_vn_EpdF_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_EpdF->GetYaxis()->GetTitle());

  TH1D *h_vn_TpcB_flip = new TH1D("h_vn_TpcB_flip",h_vn_TpcB->GetTitle(),centBins,0,centBins);
  h_vn_TpcB_flip->GetXaxis()->SetTitle((TString)h_vn_TpcB->GetXaxis()->GetTitle()+" (%)");
  h_vn_TpcB_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_TpcB->GetYaxis()->GetTitle());
  
  TH1D *h_vn_pp_flip = new TH1D("h_vn_pp_flip",h_vn_pp->GetTitle(),centBins,0,centBins);
  h_vn_pp_flip->GetXaxis()->SetTitle((TString)h_vn_pp->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_pp->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_flip = new TH1D("h_vn_pm_flip",h_vn_pm->GetTitle(),centBins,0,centBins);
  h_vn_pm_flip->GetXaxis()->SetTitle((TString)h_vn_pm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_pm->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_flip = new TH1D("h_vn_kp_flip",h_vn_kp->GetTitle(),centBins,0,centBins);
  h_vn_kp_flip->GetXaxis()->SetTitle((TString)h_vn_kp->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_kp->GetYaxis()->GetTitle());

  TH1D *h_vn_km_flip = new TH1D("h_vn_km_flip",h_vn_km->GetTitle(),centBins,0,centBins);
  h_vn_km_flip->GetXaxis()->SetTitle((TString)h_vn_km->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_km->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_flip = new TH1D("h_vn_pr_flip",h_vn_pr->GetTitle(),centBins,0,centBins);
  h_vn_pr_flip->GetXaxis()->SetTitle((TString)h_vn_pr->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_flip->GetYaxis()->SetTitle("v_{2}");//h_vn_pr->GetYaxis()->GetTitle());



  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};


  std::vector<TString> newBinLabels;

  // Get list of bin labels, but flipped
  for (int i = lastCentID; i >= firstCentID; i--) { newBinLabels.push_back(centralityBins[i]); }

  // Flip the bin contents into the new histograms
  int j = 1;
  for (int i = centBins; i >= 1; i--)
    {
      h_vn_EpdE_flip->SetBinContent(j, h_vn_EpdE->GetBinContent(i));
      h_vn_EpdE_flip->SetBinError(j, h_vn_EpdE->GetBinError(i));

      h_vn_EpdF_flip->SetBinContent(j, h_vn_EpdF->GetBinContent(i));
      h_vn_EpdF_flip->SetBinError(j, h_vn_EpdF->GetBinError(i));

      h_vn_TpcB_flip->SetBinContent(j, h_vn_TpcB->GetBinContent(i));
      h_vn_TpcB_flip->SetBinError(j, h_vn_TpcB->GetBinError(i));

      h_vn_pp_flip->SetBinContent(j, h_vn_pp->GetBinContent(i));
      h_vn_pp_flip->SetBinError(j, h_vn_pp->GetBinError(i));

      h_vn_pm_flip->SetBinContent(j, h_vn_pm->GetBinContent(i));
      h_vn_pm_flip->SetBinError(j, h_vn_pm->GetBinError(i));

      h_vn_kp_flip->SetBinContent(j, h_vn_kp->GetBinContent(i));
      h_vn_kp_flip->SetBinError(j, h_vn_kp->GetBinError(i));

      h_vn_km_flip->SetBinContent(j, h_vn_km->GetBinContent(i));
      h_vn_km_flip->SetBinError(j, h_vn_km->GetBinError(i));

      h_vn_pr_flip->SetBinContent(j, h_vn_pr->GetBinContent(i));
      h_vn_pr_flip->SetBinError(j, h_vn_pr->GetBinError(i));

      
      j++;
    }

  // Put the bin labels on the new histograms  MAYBE UNNECESSARY
  /*
  for (int i = 1; i <= centBins; i++)
    {
      h_vn_EpdE_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_vn_EpdF_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_vn_TpcB_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_vn_pp_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_vn_pm_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_vn_kp_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_vn_km_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_vn_pr_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
    }
  */



  canvas->SetLogy(0);
  canvas->SetTicks();

  TH1D *vn_EpdE = new TH1D("vn_EpdE", ";Centrality (%);v_{2}", 16, 0, 80);   // Use these to change the x axis without fighting with TAxis.
  TH1D *vn_EpdF = new TH1D("vn_EpdF", ";Centrality (%);v_{2}", 16, 0, 80);
  TH1D *vn_TpcB = new TH1D("vn_TpcB", ";Centrality (%);v_{2}", 16, 0, 80);
  TH1D *vn_pp = new TH1D("vn_pp", ";Centrality (%);v_{2}", 16, 0, 80);
  TH1D *vn_pm = new TH1D("vn_pm", ";Centrality (%);v_{2}", 16, 0, 80);
  TH1D *vn_kp = new TH1D("vn_kp", ";Centrality (%);v_{2}", 16, 0, 80);
  TH1D *vn_km = new TH1D("vn_km", ";Centrality (%);v_{2}", 16, 0, 80);
  TH1D *vn_pr = new TH1D("vn_pr", ";Centrality (%);v_{2}", 16, 0, 80);

  for (int i = 1; i <= 16; i++)
    {
      vn_EpdE->SetBinContent(i, h_vn_EpdE_flip->GetBinContent(i));
      vn_EpdE->SetBinError(i, h_vn_EpdE_flip->GetBinError(i));

      vn_EpdF->SetBinContent(i, h_vn_EpdF_flip->GetBinContent(i));
      vn_EpdF->SetBinError(i, h_vn_EpdF_flip->GetBinError(i));

      vn_TpcB->SetBinContent(i, h_vn_TpcB_flip->GetBinContent(i));
      vn_TpcB->SetBinError(i, h_vn_TpcB_flip->GetBinError(i));

      vn_pp->SetBinContent(i, h_vn_pp_flip->GetBinContent(i));
      vn_pp->SetBinError(i, h_vn_pp_flip->GetBinError(i));

      vn_pm->SetBinContent(i, h_vn_pm_flip->GetBinContent(i));
      vn_pm->SetBinError(i, h_vn_pm_flip->GetBinError(i));

      vn_kp->SetBinContent(i, h_vn_kp_flip->GetBinContent(i));
      vn_kp->SetBinError(i, h_vn_kp_flip->GetBinError(i));

      vn_km->SetBinContent(i, h_vn_km_flip->GetBinContent(i));
      vn_km->SetBinError(i, h_vn_km_flip->GetBinError(i));

      vn_pr->SetBinContent(i, h_vn_pr_flip->GetBinContent(i));
      vn_pr->SetBinError(i, h_vn_pr_flip->GetBinError(i));
    }
  

  THStack *piStack = new THStack("piStack", ";Centrality (%);v_{2}");
  THStack *kaStack = new THStack("kaStack", ";Centrality (%);v_{2}");

  vn_pp->SetMarkerStyle(20);
  vn_pp->SetMarkerColor(2);
  vn_pp->SetLineColor(2);

  vn_pm->SetMarkerStyle(20);
  vn_pm->SetMarkerColor(4);
  vn_pm->SetLineColor(4);

  vn_kp->SetMarkerStyle(20);
  vn_kp->SetMarkerColor(2);
  vn_kp->SetLineColor(2);

  vn_km->SetMarkerStyle(20);
  vn_km->SetMarkerColor(4);
  vn_km->SetLineColor(4);

  vn_pr->SetMarkerStyle(20);
  vn_pr->SetMarkerColor(2);
  vn_pr->SetLineColor(2);

  piStack->Add(vn_pp);
  piStack->Add(vn_pm);

  kaStack->Add(vn_kp);
  kaStack->Add(vn_km);

  TLegend *piLegend = new TLegend(0.775, 0.75, 0.9, 0.9);
  piLegend->AddEntry(vn_pp,"#pi^{+}");
  piLegend->AddEntry(vn_pm,"#pi^{-}");

  TLegend *kaLegend = new TLegend(0.775, 0.75, 0.9, 0.9);
  kaLegend->AddEntry(vn_kp,"K^{+}");
  kaLegend->AddEntry(vn_km,"K^{-}");

  canvas->SetLeftMargin(0.12);
  canvas->SetGrid(0);
  gStyle->SetErrorX(0);

  TLine *zeroLine = new TLine(0, 0, 80, 0);
  zeroLine->SetLineStyle(9);

  
  piStack->Draw();
  piStack->GetXaxis()->SetNdivisions(210);
  piStack->SetMaximum(0.01);
  piStack->SetMinimum(-0.12);
  piStack->Draw("NOSTACK E1P");
  zeroLine->Draw("SAME");
  piLegend->Draw();
  canvas->SaveAs(jobID + "_piStack.png");
  canvas->Clear();

  kaStack->Draw();
  kaStack->GetXaxis()->SetNdivisions(210);
  kaStack->Draw("NOSTACK E1P");
  zeroLine->Draw("SAME");
  kaLegend->Draw();
  canvas->SaveAs(jobID + "_kaStack.png");
  canvas->Clear();

  vn_pr->SetTitle("");
  vn_pr->GetXaxis()->SetNdivisions(210);
  vn_pr->Draw("E1P");
  zeroLine->Draw("SAME");
  canvas->SaveAs(jobID + "_vn_pr.png");
  canvas->Clear();

  vn_EpdE->SetMarkerStyle(20);
  vn_EpdE->SetMarkerColor(2);
  vn_EpdE->SetLineColor(2);
  vn_EpdE->GetXaxis()->SetNdivisions(210);
  vn_EpdE->Draw("E1P");
  canvas->SaveAs(jobID + "_vn_EpdE.png");
  canvas->Clear();

  vn_EpdF->SetMarkerStyle(20);
  vn_EpdF->SetMarkerColor(2);
  vn_EpdF->SetLineColor(2);
  vn_EpdF->GetXaxis()->SetNdivisions(210);
  vn_EpdF->Draw("E1P");
  canvas->SaveAs(jobID + "_vn_EpdF.png");
  canvas->Clear();

  vn_TpcB->SetMarkerStyle(20);
  vn_TpcB->SetMarkerColor(2);
  vn_TpcB->SetLineColor(2);
  vn_TpcB->GetXaxis()->SetNdivisions(210);
  vn_TpcB->Draw("E1P");
  canvas->SaveAs(jobID + "_vn_TpcB.png");
  canvas->Clear();


  resolutionInfo_INPUT->Close();
  file->Close();
}
