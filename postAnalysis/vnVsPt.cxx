void vnVsPt(TString jobID, TString order_n_str)
{
  //TH1::SetDefaultSumw2();
  
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 800);
  //canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLeftMargin(0.15);
  canvas->cd();
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);



  TProfile2D *p2_vn_pT_cent_pp = (TProfile2D*)file->Get("p2_vn_pT_cent_pp");
  TProfile2D *p2_vn_pT_cent_pm = (TProfile2D*)file->Get("p2_vn_pT_cent_pm");
  TProfile2D *p2_vn_pT_cent_kp = (TProfile2D*)file->Get("p2_vn_pT_cent_kp");
  TProfile2D *p2_vn_pT_cent_km = (TProfile2D*)file->Get("p2_vn_pT_cent_km");
  TProfile2D *p2_vn_pT_cent_pr = (TProfile2D*)file->Get("p2_vn_pT_cent_pr");
  //TProfile2D *p2_vn_pT_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_pT_cent_pr_symmetry");

  p2_vn_pT_cent_kp->RebinY();
  p2_vn_pT_cent_km->RebinY();


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


  TH1D *h_vn_pT_00to10_pp = new TH1D("h_vn_pT_00to10_pp", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pp = new TH1D("h_vn_pT_10to40_pp", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pp = new TH1D("h_vn_pT_40to60_pp", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  TH1D *h_vn_pT_00to10_pm = new TH1D("h_vn_pT_00to10_pm", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pm = new TH1D("h_vn_pT_10to40_pm", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pm = new TH1D("h_vn_pT_40to60_pm", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  TH1D *h_vn_pT_00to10_kp = new TH1D("h_vn_pT_00to10_kp", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_pT_10to40_kp = new TH1D("h_vn_pT_10to40_kp", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_pT_40to60_kp = new TH1D("h_vn_pT_40to60_kp", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);

  TH1D *h_vn_pT_00to10_km = new TH1D("h_vn_pT_00to10_km", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_pT_10to40_km = new TH1D("h_vn_pT_10to40_km", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);
  TH1D *h_vn_pT_40to60_km = new TH1D("h_vn_pT_40to60_km", ";p_{T} (GeV);v_{"+order_n_str+"}", 5, 0, 2);

  TH1D *h_vn_pT_00to10_pr = new TH1D("h_vn_pT_00to10_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pr = new TH1D("h_vn_pT_10to40_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pr = new TH1D("h_vn_pT_40to60_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);


  h_vn_pT_00to10_kp = p_vn_pT_00to10_kp->ProjectionX();
  h_vn_pT_10to40_kp = p_vn_pT_10to40_kp->ProjectionX();
  h_vn_pT_40to60_kp = p_vn_pT_40to60_kp->ProjectionX();
  h_vn_pT_00to10_km = p_vn_pT_00to10_km->ProjectionX();
  h_vn_pT_10to40_km = p_vn_pT_10to40_km->ProjectionX();
  h_vn_pT_40to60_km = p_vn_pT_40to60_km->ProjectionX();
  
  h_vn_pT_00to10_pp = p_vn_pT_00to10_pp->ProjectionX();
  h_vn_pT_10to40_pp = p_vn_pT_10to40_pp->ProjectionX();
  h_vn_pT_40to60_pp = p_vn_pT_40to60_pp->ProjectionX();

  h_vn_pT_00to10_pm = p_vn_pT_00to10_pm->ProjectionX();
  h_vn_pT_10to40_pm = p_vn_pT_10to40_pm->ProjectionX();
  h_vn_pT_40to60_pm = p_vn_pT_40to60_pm->ProjectionX();

  h_vn_pT_00to10_pr = p_vn_pT_00to10_pr->ProjectionX();
  h_vn_pT_10to40_pr = p_vn_pT_10to40_pr->ProjectionX();
  h_vn_pT_40to60_pr = p_vn_pT_40to60_pr->ProjectionX();


  THStack *ppPtStack   = new THStack("ppPtStack", ";p_{T};v_{"+order_n_str+"}");
  THStack *pmPtStack   = new THStack("pmPtStack", ";p_{T};v_{"+order_n_str+"}");
  THStack *kpPtStack   = new THStack("kpPtStack", ";p_{T};v_{"+order_n_str+"}");
  THStack *kmPtStack   = new THStack("kmPtStack", ";p_{T};v_{"+order_n_str+"}");
  THStack *prPtStack   = new THStack("prPtStack", ";p_{T};v_{"+order_n_str+"}");


  h_vn_pT_00to10_pp->SetMarkerStyle(20);
  h_vn_pT_10to40_pp->SetMarkerStyle(20);
  h_vn_pT_40to60_pp->SetMarkerStyle(20);

  h_vn_pT_00to10_pp->SetMarkerColor(2);
  h_vn_pT_10to40_pp->SetMarkerColor(4);
  h_vn_pT_40to60_pp->SetMarkerColor(8);
  h_vn_pT_00to10_pp->SetMarkerSize(2);
  h_vn_pT_10to40_pp->SetMarkerSize(2);
  h_vn_pT_40to60_pp->SetMarkerSize(2);
  h_vn_pT_00to10_pp->SetLineColor(2);
  h_vn_pT_10to40_pp->SetLineColor(4);
  h_vn_pT_40to60_pp->SetLineColor(8);
  
  h_vn_pT_00to10_pm->SetMarkerStyle(20);
  h_vn_pT_10to40_pm->SetMarkerStyle(20);
  h_vn_pT_40to60_pm->SetMarkerStyle(20);
  h_vn_pT_00to10_pm->SetMarkerColor(2);
  h_vn_pT_10to40_pm->SetMarkerColor(4);
  h_vn_pT_40to60_pm->SetMarkerColor(8);
  h_vn_pT_00to10_pm->SetMarkerSize(2);
  h_vn_pT_10to40_pm->SetMarkerSize(2);
  h_vn_pT_40to60_pm->SetMarkerSize(2);
  h_vn_pT_00to10_pm->SetLineColor(2);
  h_vn_pT_10to40_pm->SetLineColor(4);
  h_vn_pT_40to60_pm->SetLineColor(8);
  
  h_vn_pT_00to10_kp->SetMarkerStyle(20);
  h_vn_pT_10to40_kp->SetMarkerStyle(20);
  h_vn_pT_40to60_kp->SetMarkerStyle(20);
  h_vn_pT_00to10_kp->SetMarkerColor(2);
  h_vn_pT_10to40_kp->SetMarkerColor(4);
  h_vn_pT_40to60_kp->SetMarkerColor(8);
  h_vn_pT_00to10_kp->SetMarkerSize(2);
  h_vn_pT_10to40_kp->SetMarkerSize(2);
  h_vn_pT_40to60_kp->SetMarkerSize(2);
  h_vn_pT_00to10_kp->SetLineColor(2);
  h_vn_pT_10to40_kp->SetLineColor(4);
  h_vn_pT_40to60_kp->SetLineColor(8);
  
  h_vn_pT_00to10_km->SetMarkerStyle(20);
  h_vn_pT_10to40_km->SetMarkerStyle(20);
  h_vn_pT_40to60_km->SetMarkerStyle(20);
  h_vn_pT_00to10_km->SetMarkerColor(2);
  h_vn_pT_10to40_km->SetMarkerColor(4);
  h_vn_pT_40to60_km->SetMarkerColor(8);
  h_vn_pT_00to10_km->SetMarkerSize(2);
  h_vn_pT_10to40_km->SetMarkerSize(2);
  h_vn_pT_40to60_km->SetMarkerSize(2);
  h_vn_pT_00to10_km->SetLineColor(2);
  h_vn_pT_10to40_km->SetLineColor(4);
  h_vn_pT_40to60_km->SetLineColor(8);
  
  h_vn_pT_00to10_pr->SetMarkerStyle(20);
  h_vn_pT_10to40_pr->SetMarkerStyle(20);
  h_vn_pT_40to60_pr->SetMarkerStyle(20);
  h_vn_pT_00to10_pr->SetMarkerColor(2);
  h_vn_pT_10to40_pr->SetMarkerColor(4);
  h_vn_pT_40to60_pr->SetMarkerColor(8);
  h_vn_pT_00to10_pr->SetMarkerSize(2);
  h_vn_pT_10to40_pr->SetMarkerSize(2);
  h_vn_pT_40to60_pr->SetMarkerSize(2);
  h_vn_pT_00to10_pr->SetLineColor(2);
  h_vn_pT_10to40_pr->SetLineColor(4);
  h_vn_pT_40to60_pr->SetLineColor(8);
  

  if (order_n_str == "2")
    {
      ppPtStack->Add(h_vn_pT_00to10_pp);
      ppPtStack->Add(h_vn_pT_10to40_pp);
      ppPtStack->Add(h_vn_pT_40to60_pp);

      pmPtStack->Add(h_vn_pT_00to10_pm);
      pmPtStack->Add(h_vn_pT_10to40_pm);
      pmPtStack->Add(h_vn_pT_40to60_pm);

      kpPtStack->Add(h_vn_pT_00to10_kp);
      kpPtStack->Add(h_vn_pT_10to40_kp);
      kpPtStack->Add(h_vn_pT_40to60_kp);

      kmPtStack->Add(h_vn_pT_10to40_km);

      prPtStack->Add(h_vn_pT_00to10_pr);
      prPtStack->Add(h_vn_pT_10to40_pr);
      prPtStack->Add(h_vn_pT_40to60_pr);


      ppPtStack->Draw();
      ppPtStack->GetYaxis()->SetTitleOffset(1.7);
      ppPtStack->GetXaxis()->SetNdivisions(210);
      ppPtStack->Draw();
      //ppPtStack->SetMaximum(0.03);
      //ppPtStack->SetMinimum(-0.05);
      ppPtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //ppLegend->Draw();
      //ppText->Draw();
      canvas->SaveAs(jobID + "_ppPtStack.png");
      canvas->Clear();

      pmPtStack->Draw();
      pmPtStack->GetYaxis()->SetTitleOffset(1.7);
      pmPtStack->GetXaxis()->SetNdivisions(210);
      pmPtStack->Draw();
      //pmPtStack->SetMaximum(0.03);
      //pmPtStack->SetMinimum(-0.04);
      pmPtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //pmLegend->Draw();
      //pmText->Draw();
      canvas->SaveAs(jobID + "_pmPtStack.png");
      canvas->Clear();

      kpPtStack->Draw();
      kpPtStack->GetYaxis()->SetTitleOffset(1.7);
      kpPtStack->GetXaxis()->SetNdivisions(210);
      kpPtStack->Draw();
      //kpPtStack->SetMaximum(0.08);
      //kpPtStack->SetMinimum(-0.13);
      kpPtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //kpLegend->Draw();
      //kpText->Draw();
      canvas->SaveAs(jobID + "_kpPtStack.png");
      canvas->Clear();

      kmPtStack->Draw();
      kmPtStack->GetYaxis()->SetTitleOffset(1.7);
      kmPtStack->GetXaxis()->SetNdivisions(210);
      kmPtStack->Draw();
      //kmPtStack->SetMaximum(0.02);
      //kmPtStack->SetMinimum(-0.06);
      kmPtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //kmLegend->Draw();
      //kmText->Draw();
      canvas->SaveAs(jobID + "_kmPtStack.png");
      canvas->Clear();

      prPtStack->Draw();
      prPtStack->GetYaxis()->SetTitleOffset(1.7);
      prPtStack->GetXaxis()->SetNdivisions(210);
      prPtStack->Draw();
      //prPtStack->SetMaximum(0.1);
      //prPtStack->SetMinimum(-0.08);
      prPtStack->Draw("NOSTACK E1P");
      //zeroLine_y->Draw("SAME");
      //prLegend->Draw();
      //prText_y->Draw();
      canvas->SaveAs(jobID + "_prPtStack.png");
      canvas->Clear();
    }
  else if (order_n_str == "3")
    {
      TLine *zeroLine_pt = new TLine(0, 0, 2, 0);
      zeroLine_pt->SetLineStyle(9);

      
      ppPtStack->Add(h_vn_pT_00to10_pp);
      ppPtStack->Add(h_vn_pT_10to40_pp);
      ppPtStack->Add(h_vn_pT_40to60_pp);

      pmPtStack->Add(h_vn_pT_00to10_pm);
      pmPtStack->Add(h_vn_pT_10to40_pm);
      pmPtStack->Add(h_vn_pT_40to60_pm);

      kpPtStack->Add(h_vn_pT_00to10_kp);
      kpPtStack->Add(h_vn_pT_10to40_kp);
      kpPtStack->Add(h_vn_pT_40to60_kp);

      kmPtStack->Add(h_vn_pT_10to40_km);

      prPtStack->Add(h_vn_pT_00to10_pr);
      prPtStack->Add(h_vn_pT_10to40_pr);
      prPtStack->Add(h_vn_pT_40to60_pr);


      ppPtStack->Draw();
      ppPtStack->GetYaxis()->SetTitleOffset(1.7);
      ppPtStack->GetXaxis()->SetNdivisions(210);
      ppPtStack->Draw();
      //ppPtStack->SetMaximum(0.03);
      ppPtStack->SetMinimum(-0.1);
      ppPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      //ppLegend->Draw();
      //ppText->Draw();
      canvas->SaveAs(jobID + "_ppPtStack.png");
      canvas->Clear();

      pmPtStack->Draw();
      pmPtStack->GetYaxis()->SetTitleOffset(1.7);
      pmPtStack->GetXaxis()->SetNdivisions(210);
      pmPtStack->Draw();
      //pmPtStack->SetMaximum(0.03);
      //pmPtStack->SetMinimum(-0.04);
      pmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      //pmLegend->Draw();
      //pmText->Draw();
      canvas->SaveAs(jobID + "_pmPtStack.png");
      canvas->Clear();

      kpPtStack->Draw();
      kpPtStack->GetYaxis()->SetTitleOffset(1.7);
      kpPtStack->GetXaxis()->SetNdivisions(210);
      kpPtStack->Draw();
      //kpPtStack->SetMaximum(0.08);
      kpPtStack->SetMinimum(-0.1);
      kpPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      //kpLegend->Draw();
      //kpText->Draw();
      canvas->SaveAs(jobID + "_kpPtStack.png");
      canvas->Clear();

      kmPtStack->Draw();
      kmPtStack->GetYaxis()->SetTitleOffset(1.7);
      kmPtStack->GetXaxis()->SetNdivisions(210);
      kmPtStack->Draw();
      //kmPtStack->SetMaximum(0.02);
      //kmPtStack->SetMinimum(-0.06);
      kmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      //kmLegend->Draw();
      //kmText->Draw();
      canvas->SaveAs(jobID + "_kmPtStack.png");
      canvas->Clear();

      prPtStack->Draw();
      prPtStack->GetYaxis()->SetTitleOffset(1.7);
      prPtStack->GetXaxis()->SetNdivisions(210);
      prPtStack->Draw();
      prPtStack->SetMaximum(0.02);
      prPtStack->SetMinimum(-0.11);
      prPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      //prLegend->Draw();
      //prText_y->Draw();
      canvas->SaveAs(jobID + "_prPtStack.png");
      canvas->Clear();
    }

  file->Close();
  
}
