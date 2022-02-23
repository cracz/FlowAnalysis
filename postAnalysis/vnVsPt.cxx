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
  canvas->SetLeftMargin(0.14);
  canvas->cd();
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(6);


  TProfile2D *p2_vn_pT_cent_pp = (TProfile2D*)file->Get("p2_vn_pT_cent_pp");
  TProfile2D *p2_vn_pT_cent_pm = (TProfile2D*)file->Get("p2_vn_pT_cent_pm");
  TProfile2D *p2_vn_pT_cent_kp = (TProfile2D*)file->Get("p2_vn_pT_cent_kp");
  TProfile2D *p2_vn_pT_cent_km = (TProfile2D*)file->Get("p2_vn_pT_cent_km");
  TProfile2D *p2_vn_pT_cent_pr = (TProfile2D*)file->Get("p2_vn_pT_cent_pr");
  TProfile2D *p2_vn_pT_cent_pr_alt = (TProfile2D*)file->Get("p2_vn_pT_cent_pr_alt");
  //TProfile2D *p2_vn_pT_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_pT_cent_pr_symmetry");
  TProfile2D *p2_vn_pT_cent_de = (TProfile2D*)file->Get("p2_vn_pT_cent_de");
  TProfile2D *p2_vn_pT_cent_tr = (TProfile2D*)file->Get("p2_vn_pT_cent_tr");

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

  TProfile *p_vn_pT_00to10_pr_alt = p2_vn_pT_cent_pr_alt->ProfileY("p_vn_pT_00to10_pr_alt", 15, 16);
  TProfile *p_vn_pT_10to40_pr_alt = p2_vn_pT_cent_pr_alt->ProfileY("p_vn_pT_10to40_pr_alt", 9, 14);
  TProfile *p_vn_pT_40to60_pr_alt = p2_vn_pT_cent_pr_alt->ProfileY("p_vn_pT_40to60_pr_alt", 5, 8);

  TProfile *p_vn_pT_00to10_de = p2_vn_pT_cent_de->ProfileY("p_vn_pT_00to10_de", 15, 16);
  TProfile *p_vn_pT_10to40_de = p2_vn_pT_cent_de->ProfileY("p_vn_pT_10to40_de", 9, 14);
  TProfile *p_vn_pT_40to60_de = p2_vn_pT_cent_de->ProfileY("p_vn_pT_40to60_de", 5, 8);

  TProfile *p_vn_pT_00to10_tr = p2_vn_pT_cent_tr->ProfileY("p_vn_pT_00to10_tr", 15, 16);
  TProfile *p_vn_pT_10to40_tr = p2_vn_pT_cent_tr->ProfileY("p_vn_pT_10to40_tr", 9, 14);
  TProfile *p_vn_pT_40to60_tr = p2_vn_pT_cent_tr->ProfileY("p_vn_pT_40to60_tr", 5, 8);


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

  TH1D *h_vn_pT_00to10_pr_alt = new TH1D("h_vn_pT_00to10_pr_alt", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);
  TH1D *h_vn_pT_10to40_pr_alt = new TH1D("h_vn_pT_10to40_pr_alt", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);
  TH1D *h_vn_pT_40to60_pr_alt = new TH1D("h_vn_pT_40to60_pr_alt", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);

  TH1D *h_vn_pT_00to10_pr = new TH1D("h_vn_pT_00to10_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_10to40_pr = new TH1D("h_vn_pT_10to40_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  TH1D *h_vn_pT_40to60_pr = new TH1D("h_vn_pT_40to60_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);

  TH1D *h_vn_pT_00to10_de = new TH1D("h_vn_pT_00to10_de", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);
  TH1D *h_vn_pT_10to40_de = new TH1D("h_vn_pT_10to40_de", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);
  TH1D *h_vn_pT_40to60_de = new TH1D("h_vn_pT_40to60_de", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);

  TH1D *h_vn_pT_00to10_tr = new TH1D("h_vn_pT_00to10_tr", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);
  TH1D *h_vn_pT_10to40_tr = new TH1D("h_vn_pT_10to40_tr", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);
  TH1D *h_vn_pT_40to60_tr = new TH1D("h_vn_pT_40to60_tr", ";p_{T} (GeV);v_{"+order_n_str+"}", 15, 0, 2.5);


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

  h_vn_pT_00to10_pr_alt = p_vn_pT_00to10_pr_alt->ProjectionX();
  h_vn_pT_10to40_pr_alt = p_vn_pT_10to40_pr_alt->ProjectionX();
  h_vn_pT_40to60_pr_alt = p_vn_pT_40to60_pr_alt->ProjectionX();

  h_vn_pT_00to10_pr = p_vn_pT_00to10_pr->ProjectionX();
  h_vn_pT_10to40_pr = p_vn_pT_10to40_pr->ProjectionX();
  h_vn_pT_40to60_pr = p_vn_pT_40to60_pr->ProjectionX();

  h_vn_pT_00to10_de = p_vn_pT_00to10_de->ProjectionX();
  h_vn_pT_10to40_de = p_vn_pT_10to40_de->ProjectionX();
  h_vn_pT_40to60_de = p_vn_pT_40to60_de->ProjectionX();

  h_vn_pT_00to10_tr = p_vn_pT_00to10_tr->ProjectionX();
  h_vn_pT_10to40_tr = p_vn_pT_10to40_tr->ProjectionX();
  h_vn_pT_40to60_tr = p_vn_pT_40to60_tr->ProjectionX();


  THStack *ppPtStack   = new THStack("ppPtStack", ";p_{T} GeV;v_{"+order_n_str+"}");
  THStack *pmPtStack   = new THStack("pmPtStack", ";p_{T} GeV;v_{"+order_n_str+"}");
  THStack *kpPtStack   = new THStack("kpPtStack", ";p_{T} GeV;v_{"+order_n_str+"}");
  THStack *kmPtStack   = new THStack("kmPtStack", ";p_{T} GeV;v_{"+order_n_str+"}");
  THStack *prPtStack   = new THStack("prPtStack", ";p_{T} GeV;v_{"+order_n_str+"}");
  THStack *prPtStack_alt   = new THStack("prPtStack_alt", ";p_{T} GeV;v_{"+order_n_str+"}");
  THStack *dePtStack   = new THStack("dePtStack", ";p_{T} GeV;v_{"+order_n_str+"}");
  THStack *trPtStack   = new THStack("trPtStack", ";p_{T} GeV;v_{"+order_n_str+"}");


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

  h_vn_pT_00to10_pr_alt->SetMarkerStyle(20);
  h_vn_pT_10to40_pr_alt->SetMarkerStyle(20);
  h_vn_pT_40to60_pr_alt->SetMarkerStyle(20);
  h_vn_pT_00to10_pr_alt->SetMarkerColor(2);
  h_vn_pT_10to40_pr_alt->SetMarkerColor(4);
  h_vn_pT_40to60_pr_alt->SetMarkerColor(8);
  h_vn_pT_00to10_pr_alt->SetMarkerSize(2);
  h_vn_pT_10to40_pr_alt->SetMarkerSize(2);
  h_vn_pT_40to60_pr_alt->SetMarkerSize(2);
  h_vn_pT_00to10_pr_alt->SetLineColor(2);
  h_vn_pT_10to40_pr_alt->SetLineColor(4);
  h_vn_pT_40to60_pr_alt->SetLineColor(8);
  
  h_vn_pT_00to10_de->SetMarkerStyle(20);
  h_vn_pT_10to40_de->SetMarkerStyle(20);
  h_vn_pT_40to60_de->SetMarkerStyle(20);
  h_vn_pT_00to10_de->SetMarkerColor(2);
  h_vn_pT_10to40_de->SetMarkerColor(4);
  h_vn_pT_40to60_de->SetMarkerColor(8);
  h_vn_pT_00to10_de->SetMarkerSize(2);
  h_vn_pT_10to40_de->SetMarkerSize(2);
  h_vn_pT_40to60_de->SetMarkerSize(2);
  h_vn_pT_00to10_de->SetLineColor(2);
  h_vn_pT_10to40_de->SetLineColor(4);
  h_vn_pT_40to60_de->SetLineColor(8);

  h_vn_pT_00to10_tr->SetMarkerStyle(20);
  h_vn_pT_10to40_tr->SetMarkerStyle(20);
  h_vn_pT_40to60_tr->SetMarkerStyle(20);
  h_vn_pT_00to10_tr->SetMarkerColor(2);
  h_vn_pT_10to40_tr->SetMarkerColor(4);
  h_vn_pT_40to60_tr->SetMarkerColor(8);
  h_vn_pT_00to10_tr->SetMarkerSize(2);
  h_vn_pT_10to40_tr->SetMarkerSize(2);
  h_vn_pT_40to60_tr->SetMarkerSize(2);
  h_vn_pT_00to10_tr->SetLineColor(2);
  h_vn_pT_10to40_tr->SetLineColor(4);
  h_vn_pT_40to60_tr->SetLineColor(8);


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

      TLine *zeroLine_pdt = new TLine(0, 0, 2.5, 0);
      zeroLine_pdt->SetLineStyle(9);

      
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

      prPtStack_alt->Add(h_vn_pT_00to10_pr_alt);
      prPtStack_alt->Add(h_vn_pT_10to40_pr_alt);
      prPtStack_alt->Add(h_vn_pT_40to60_pr_alt);

      dePtStack->Add(h_vn_pT_00to10_de);
      dePtStack->Add(h_vn_pT_10to40_de);
      dePtStack->Add(h_vn_pT_40to60_de);

      trPtStack->Add(h_vn_pT_00to10_tr);
      trPtStack->Add(h_vn_pT_10to40_tr);
      trPtStack->Add(h_vn_pT_40to60_tr);


      TLegend *ppLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      ppLegend->AddEntry(h_vn_pT_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(h_vn_pT_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(h_vn_pT_40to60_pp, "40 - 60%");
      ppLegend->SetBorderSize(0);
      ppLegend->SetFillColorAlpha(0,0);

      TLegend *pmLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      pmLegend->AddEntry(h_vn_pT_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(h_vn_pT_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(h_vn_pT_40to60_pm, "40 - 60%");
      pmLegend->SetBorderSize(0);
      pmLegend->SetFillColorAlpha(0,0);

      TLegend *kpLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      kpLegend->AddEntry(h_vn_pT_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(h_vn_pT_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(h_vn_pT_40to60_kp, "40 - 60%");
      kpLegend->SetBorderSize(0);
      kpLegend->SetFillColorAlpha(0,0);

      TLegend *kmLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      //kmLegend->AddEntry(h_vn_pT_00to10_km, "0 - 10%");
      kmLegend->AddEntry(h_vn_pT_10to40_km, "10 - 40%");
      //kmLegend->AddEntry(h_vn_pT_40to60_km, "40 - 60%");
      kmLegend->SetBorderSize(0);
      kmLegend->SetFillColorAlpha(0,0);

      TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      prLegend->AddEntry(h_vn_pT_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_pT_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_pT_40to60_pr, "40 - 60%");
      prLegend->SetBorderSize(0);
      prLegend->SetFillColorAlpha(0,0);

      TLegend *deLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      deLegend->AddEntry(h_vn_pT_00to10_de, "0 - 10%");
      deLegend->AddEntry(h_vn_pT_10to40_de, "10 - 40%");
      deLegend->AddEntry(h_vn_pT_40to60_de, "40 - 60%");
      deLegend->SetBorderSize(0);
      deLegend->SetFillColorAlpha(0,0);

      TLegend *trLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      trLegend->AddEntry(h_vn_pT_00to10_tr, "0 - 10%");
      trLegend->AddEntry(h_vn_pT_10to40_tr, "10 - 40%");
      trLegend->AddEntry(h_vn_pT_40to60_tr, "40 - 60%");
      trLegend->SetBorderSize(0);
      trLegend->SetFillColorAlpha(0,0);

      

      TPaveText *ppText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0 < y_{CM} < 0.5 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);

      TPaveText *pmText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0 < y_{CM} < 0.5 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);

      TPaveText *kpText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0 < y_{CM} < 0.5 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);

      TPaveText *kmText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0 < y_{CM} < 0.5 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);

      TPaveText *prText_alt = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      prText_alt->AddText("Proton");
      prText_alt->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_alt->AddText("0.5 < y_{CM} < 1.0 GeV");
      prText_alt->SetFillColorAlpha(0,0);
      prText_alt->SetLineColorAlpha(0,0);

      TPaveText *deText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      deText->AddText("Deuteron");
      deText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      deText->AddText("0.5 < y_{CM} < 1.0 GeV");
      deText->SetFillColorAlpha(0,0);
      deText->SetLineColorAlpha(0,0);

      TPaveText *trText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      trText->AddText("Triton");
      trText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      trText->AddText("0.5 < y_{CM} < 1.0 GeV");
      trText->SetFillColorAlpha(0,0);
      trText->SetLineColorAlpha(0,0);


      Double_t ptUpperBound = 0.25;
      Double_t ptLowerBound = -0.25;


      ppPtStack->Draw();
      ppPtStack->GetYaxis()->SetTitleOffset(1.7);
      ppPtStack->GetXaxis()->SetNdivisions(210);
      ppPtStack->Draw();
      ppPtStack->SetMaximum(ptUpperBound);
      ppPtStack->SetMinimum(ptLowerBound);
      ppPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      ppPtStack->Draw("NOSTACK E1P SAME");
      ppLegend->Draw();
      ppText->Draw();
      canvas->SaveAs(jobID + "_ppPtStack.png");
      canvas->Clear();

      pmPtStack->Draw();
      pmPtStack->GetYaxis()->SetTitleOffset(1.7);
      pmPtStack->GetXaxis()->SetNdivisions(210);
      pmPtStack->Draw();
      pmPtStack->SetMaximum(ptUpperBound);
      pmPtStack->SetMinimum(ptLowerBound);
      pmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      pmPtStack->Draw("NOSTACK E1P SAME");
      pmLegend->Draw();
      pmText->Draw();
      canvas->SaveAs(jobID + "_pmPtStack.png");
      canvas->Clear();

      kpPtStack->Draw();
      kpPtStack->GetYaxis()->SetTitleOffset(1.7);
      kpPtStack->GetXaxis()->SetNdivisions(210);
      kpPtStack->Draw();
      kpPtStack->SetMaximum(ptUpperBound);
      kpPtStack->SetMinimum(ptLowerBound);
      kpPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      kpPtStack->Draw("NOSTACK E1P SAME");
      kpLegend->Draw();
      kpText->Draw();
      canvas->SaveAs(jobID + "_kpPtStack.png");
      canvas->Clear();
      
      kmPtStack->Draw();
      kmPtStack->GetYaxis()->SetTitleOffset(1.7);
      kmPtStack->GetXaxis()->SetNdivisions(210);
      kmPtStack->Draw();
      kmPtStack->SetMaximum(ptUpperBound);
      kmPtStack->SetMinimum(ptLowerBound);
      kmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      kmPtStack->Draw("NOSTACK E1P SAME");
      kmLegend->Draw();
      kmText->Draw();
      canvas->SaveAs(jobID + "_kmPtStack.png");
      canvas->Clear();
      gStyle->SetErrorX(0);
      
      prPtStack->Draw();
      prPtStack->GetYaxis()->SetTitleOffset(1.7);
      prPtStack->GetXaxis()->SetNdivisions(210);
      prPtStack->Draw();
      prPtStack->SetMaximum(ptUpperBound);
      prPtStack->SetMinimum(ptLowerBound);
      prPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      prPtStack->Draw("NOSTACK E1P SAME");
      prLegend->Draw();
      prText->Draw();
      canvas->SaveAs(jobID + "_prPtStack.png");
      canvas->Clear();


      prPtStack_alt->Draw();
      prPtStack_alt->GetYaxis()->SetTitleOffset(1.7);
      prPtStack_alt->GetXaxis()->SetNdivisions(210);
      prPtStack_alt->Draw();
      prPtStack_alt->SetMaximum(ptUpperBound);
      prPtStack_alt->SetMinimum(ptLowerBound);
      prPtStack_alt->Draw("NOSTACK E1P");
      zeroLine_pdt->Draw("SAME");
      prPtStack_alt->Draw("NOSTACK E1P SAME");
      prLegend->Draw();
      prText_alt->Draw();
      canvas->SaveAs(jobID + "_prPtStack_alt.png");
      canvas->Clear();

      dePtStack->Draw();
      dePtStack->GetYaxis()->SetTitleOffset(1.7);
      dePtStack->GetXaxis()->SetNdivisions(210);
      dePtStack->Draw();
      dePtStack->SetMaximum(ptUpperBound);
      dePtStack->SetMinimum(ptLowerBound);
      dePtStack->Draw("NOSTACK E1P");
      zeroLine_pdt->Draw("SAME");
      dePtStack->Draw("NOSTACK E1P SAME");
      deLegend->Draw();
      deText->Draw();
      canvas->SaveAs(jobID + "_dePtStack.png");
      canvas->Clear();

      trPtStack->Draw();
      trPtStack->GetYaxis()->SetTitleOffset(1.7);
      trPtStack->GetXaxis()->SetNdivisions(210);
      trPtStack->Draw();
      trPtStack->SetMaximum(ptUpperBound);
      trPtStack->SetMinimum(ptLowerBound);
      trPtStack->Draw("NOSTACK E1P");
      zeroLine_pdt->Draw("SAME");
      trPtStack->Draw("NOSTACK E1P SAME");
      trLegend->Draw();
      trText->Draw();
      canvas->SaveAs(jobID + "_trPtStack.png");
      canvas->Clear();
    }

  file->Close();
  
}
