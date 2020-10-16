#include "TObjArray.h"

void subsEFpostAnalysis(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  ////////
  // PLOTTING EVERYTHING FROM THE FILE
  ////////

  TH1F  *hist;
  TString name;
  TIter next(file->GetListOfKeys());
  TKey  *key;

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 875, 675);

  /*
  while ((key = (TKey*)next()))
    {
      canvas->SetCanvasSize(875, 675);   // reset canvas size
      
      hist = (TH1F*)key->ReadObj();
      name = hist->GetName();

      if (hist->GetEntries() == 0) { continue; }

      canvas->SetLogy(0);    //Reset y axis to linear
      canvas->SetGrid();     //Draw grid lines
      gStyle->SetOptStat(1); //Reset stat box
      //gStyle->SetPalette();  //Reset color palette

      Bool_t errors = false; //Draw error bars
      
      TClass *cl = gROOT->GetClass(key->GetClassName());

      if (cl->InheritsFrom("TH2"))                       //Keep this if/else order!! A TH2 still inherits from TH1.
	{
	  canvas->SetRightMargin(0.14);
	  canvas->SetLeftMargin(0.12);
	  hist->SetTitleOffset(1.2, "y");
	  gStyle->SetOptStat(11);
	  //gStyle->SetPalette(53); //kDarkBodyRadiator
	  //gStyle->SetPalette(57); //kBird Doesn't work
	  canvas->SetLogz();

	  if (name.Contains("psi") || name.Contains("Search") || name.Contains("MvsY"))
	    {
	      canvas->SetLogz(0); // Remove log z
	    }
	  else if (name.Contains("vtx"))
	    {
	      //hist->GetZaxis()->SetRangeUser(1,10000);
	      canvas->SetCanvasSize(700,700);
	    }
	  else if (name.Contains("y_vs_eta_pt"))
	    {
	      //gStyle->SetOptStat(211);
	      continue;
	    }
	  else if (name.Contains("y_vs_eta"))
	    {
	      gStyle->SetOptStat(211);
	    }

	  hist->Draw("COLZ");
	  canvas->Update();
	}
      else if (cl->InheritsFrom("TH1"))
	{
	  hist->SetTitleOffset(1.2, "y");

	  if (name.Contains("sinAvgs") || name.Contains("cosAvgs") || name.Contains("Xn") || name.Contains("Yn"))
	    {
	      continue;
	    }
	  else if (name.Contains("vtx") || name.Contains("pT") || name.Contains("nhits") || name.Contains("mom")
		   || name.Contains("primTracks") || name.Contains("dndm"))
	    {
	      canvas->SetLogy();
	    }
	  else if (name.Contains("p_") && !name.Contains("pp") && !name.Contains("kp"))
	    {
	      canvas->SetLeftMargin(0.15);
	      hist->SetTitleOffset(1.8, "y");
	      gStyle->SetOptStat(0);
	    }
	  else if (name.Contains("TOF_beta"))
	    {
	      canvas->SetLogy();
	      hist->GetYaxis()->SetRangeUser(0.1,10E+7);
	      gStyle->SetOptStat(11);
	    }
	  else if (name.Contains("check"))
	    {
	      hist->SetFillColorAlpha(4,0.6);
	    }
	  else if (name.Contains("dndm")) 
	    { 
	      hist->SetMinimum(0.1);
	      canvas->SetLogy();
	      errors = true;
	    }


	  if (errors) 
	    hist->Draw("E1");
	  else 
	    hist->Draw();

	  canvas->Update();
	}
      else
	continue;
      
      canvas->SaveAs(jobID+"_"+name+".png");
    }

  canvas->Clear();
*/
  ////////
  // END PLOTTING SECTION
  ////////


  ////////
  // PLOTTING m^2 BY p_T BINS
  ////////
  canvas->SetGrid();
  canvas->SetLogy();
  gStyle->SetOptStat(0);

  Double_t low_pT_values[4]  = {1.3, 1.6, 1.9, 2.2};
  Double_t high_pT_values[4] = {1.6, 1.9, 2.2, 2.5};

  TH2D *h2_m2_vs_qpT = (TH2D*)file->Get("h2_m2_vs_qpT");

  for(int i = 0; i < 4; i++)
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

  gStyle->SetOptStat(1);		      
  ////////
  // END PLOTTING m^2 BY p_T BINS
  ////////


  
  /*
  ////////
  // RAPIDITY AND PSEUDO-RAPIDITY PLOTS
  ////////
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
  /*
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
*/
  /*
  TH1D *h_centralities = (TH1D*)file->Get("h_centralities");
  
  ////////
  // PLOTTING CORRELATIONS AND RESOLUTIONS
  ////////

  TFile *resolutionInfo_INPUT = new TFile("resolutionInfo_INPUT.root", "RECREATE");
  
  TProfile *h_EpdEEpdF = (TProfile*)file->Get("p_EpdEEpdF");

  TProfile *h_TpcBEpdE = (TProfile*)file->Get("p_TpcBEpdE");
  TProfile *h_TpcBEpdF = (TProfile*)file->Get("p_TpcBEpdF");

  
  Int_t centBins    = h_EpdEEpdF->GetNbinsX();
  Int_t firstCentID = h_EpdEEpdF->GetBinLowEdge(1);
  Int_t lastCentID  = h_EpdEEpdF->GetBinLowEdge(h_EpdEEpdF->GetNbinsX());

  
  TH1D *h_EpdEEpdF_flip = new TH1D("h_EpdEEpdF_flip",h_EpdEEpdF->GetTitle(),centBins,1,centBins+1);
  h_EpdEEpdF_flip->GetXaxis()->SetTitle((TString)h_EpdEEpdF->GetXaxis()->GetTitle()+" (%)");
  h_EpdEEpdF_flip->GetYaxis()->SetTitle(h_EpdEEpdF->GetYaxis()->GetTitle());

  TH1D *h_TpcBEpdE_flip = new TH1D("h_TpcBEpdE_flip",h_TpcBEpdE->GetTitle(),centBins,1,centBins+1);
  h_TpcBEpdE_flip->GetXaxis()->SetTitle((TString)h_TpcBEpdE->GetXaxis()->GetTitle()+" (%)");
  h_TpcBEpdE_flip->GetYaxis()->SetTitle(h_TpcBEpdE->GetYaxis()->GetTitle());
  TH1D *h_TpcBEpdF_flip = new TH1D("h_TpcBEpdF_flip",h_TpcBEpdF->GetTitle(),centBins,1,centBins+1);
  h_TpcBEpdF_flip->GetXaxis()->SetTitle((TString)h_TpcBEpdF->GetXaxis()->GetTitle()+" (%)");
  h_TpcBEpdF_flip->GetYaxis()->SetTitle(h_TpcBEpdF->GetYaxis()->GetTitle());

  TH1D *h_centralities_flip = new TH1D("h_centralities_flip",h_centralities->GetTitle(),centBins,1,centBins+1);
  h_centralities_flip->GetXaxis()->SetTitle((TString)h_centralities->GetXaxis()->GetTitle()+" (%)");
  h_centralities_flip->GetYaxis()->SetTitle(h_centralities->GetYaxis()->GetTitle());


  // Make the possible resolution plots
  TH1D *h_resolEvsF = new TH1D("h_resolEvsF","EPD E vs EPD F and TPC B;Centrality (%);R^{E}_{21}",centBins,1,centBins+1);
  TH1D *h_resolFvsE = new TH1D("h_resolFvsE","EPD F vs EPD E and TPC B;Centrality (%);R^{F}_{21}",centBins,1,centBins+1);
  TH1D *h_resolTpcB = new TH1D("h_resolTpcB","TPC B vs EPD E and EPD F;Centrality (%);R^{B}_{21}",centBins,1,centBins+1);

  TH1D *h_resolution = new TH1D("h_resolution","TITLE;Centrality;R_{21}",centBins,1,centBins+1);
  
  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};

  std::vector<TString> newBinLabels;

  // Get list of bin labels, but flipped
  for (int i = lastCentID; i >= firstCentID; i--) { newBinLabels.push_back(centralityBins[i]); }

  // Flip the bin contents into the new histograms
  int j = 1;
  for (int i = centBins; i >= 1; i--)
    {
      h_EpdEEpdF_flip->SetBinContent(j, h_EpdEEpdF->GetBinContent(i));
      h_EpdEEpdF_flip->SetBinError(j, h_EpdEEpdF->GetBinError(i));

      h_TpcBEpdE_flip->SetBinContent(j, h_TpcBEpdE->GetBinContent(i));
      h_TpcBEpdE_flip->SetBinError(j, h_TpcBEpdE->GetBinError(i));
      h_TpcBEpdF_flip->SetBinContent(j, h_TpcBEpdF->GetBinContent(i));
      h_TpcBEpdF_flip->SetBinError(j, h_TpcBEpdF->GetBinError(i));

      h_centralities_flip->SetBinContent(j, h_centralities->GetBinContent(i));

      j++;
    }

  
  Double_t EpdEEpdF;
  Double_t TpcBEpdE;
  Double_t TpcBEpdF;

  Double_t dEpdEEpdF;
  Double_t dTpcBEpdE;
  Double_t dTpcBEpdF;

  Double_t R_EvsF;
  Double_t R_FvsE;
  Double_t R_TpcB;
  Double_t dR_EvsF;
  Double_t dR_FvsE;
  Double_t dR_TpcB;

  Double_t EpdEEpdF_save;
  Double_t TpcBEpdE_save;
  Double_t TpcBEpdF_save;

  Double_t dEpdEEpdF_save;
  Double_t dTpcBEpdE_save;
  Double_t dTpcBEpdF_save;

  Double_t R_EvsF_save;
  Double_t dR_EvsF_save;
      
  // Fill resolution plots
  for (int i = 1; i <= centBins; i++)
    {
      EpdEEpdF_save = h_EpdEEpdF->GetBinContent(i);      
      TpcBEpdE_save = h_TpcBEpdE->GetBinContent(i);
      TpcBEpdF_save = h_TpcBEpdF->GetBinContent(i);

      dEpdEEpdF_save = h_EpdEEpdF->GetBinError(i);
      dTpcBEpdE_save = h_TpcBEpdE->GetBinError(i);
      dTpcBEpdF_save = h_TpcBEpdF->GetBinError(i);

      R_EvsF_save  = TMath::Sqrt( (EpdEEpdF_save * TpcBEpdE_save) / TpcBEpdF_save );
      dR_EvsF_save = R_EvsF_save * TMath::Sqrt((dEpdEEpdF_save/(2*EpdEEpdF_save))*(dEpdEEpdF_save/(2*EpdEEpdF_save)) +
					       (dTpcBEpdE_save/(2*TpcBEpdE_save))*(dTpcBEpdE_save/(2*TpcBEpdE_save)) +
					       (dTpcBEpdF_save/(2*TpcBEpdF_save))*(dTpcBEpdF_save/(2*TpcBEpdF_save)));
	    
      
      EpdEEpdF = h_EpdEEpdF_flip->GetBinContent(i);      
      TpcBEpdE = h_TpcBEpdE_flip->GetBinContent(i);
      TpcBEpdF = h_TpcBEpdF_flip->GetBinContent(i);

      dEpdEEpdF = h_EpdEEpdF_flip->GetBinError(i);
      dTpcBEpdE = h_TpcBEpdE_flip->GetBinError(i);
      dTpcBEpdF = h_TpcBEpdF_flip->GetBinError(i);

      R_EvsF = TMath::Sqrt( (EpdEEpdF * TpcBEpdE) / TpcBEpdF );
      R_FvsE = TMath::Sqrt( (EpdEEpdF * TpcBEpdF) / TpcBEpdE );
      R_TpcB = TMath::Sqrt( (TpcBEpdE * TpcBEpdF) / EpdEEpdF );

      dR_EvsF = R_EvsF * TMath::Sqrt((dEpdEEpdF/(2*EpdEEpdF))*(dEpdEEpdF/(2*EpdEEpdF)) +
				     (dTpcBEpdE/(2*TpcBEpdE))*(dTpcBEpdE/(2*TpcBEpdE)) +
				     (dTpcBEpdF/(2*TpcBEpdF))*(dTpcBEpdF/(2*TpcBEpdF)));

      dR_FvsE = R_FvsE * TMath::Sqrt((dEpdEEpdF/(2*EpdEEpdF))*(dEpdEEpdF/(2*EpdEEpdF)) +
				     (dTpcBEpdE/(2*TpcBEpdE))*(dTpcBEpdE/(2*TpcBEpdE)) +
				     (dTpcBEpdF/(2*TpcBEpdF))*(dTpcBEpdF/(2*TpcBEpdF)));

      dR_TpcB = R_TpcB * TMath::Sqrt((dTpcBEpdE/(2*TpcBEpdE))*(dTpcBEpdE/(2*TpcBEpdE)) +
				     (dTpcBEpdF/(2*TpcBEpdF))*(dTpcBEpdF/(2*TpcBEpdF)) +
				     (dEpdEEpdF/(2*EpdEEpdF))*(dEpdEEpdF/(2*EpdEEpdF)));

      if(TMath::IsNaN(R_EvsF)) { R_EvsF = 0; dR_EvsF = 0; }
      if(TMath::IsNaN(R_FvsE)) { R_FvsE = 0; dR_FvsE = 0; }
      if(TMath::IsNaN(R_TpcB)) { R_TpcB = 0; dR_TpcB = 0; }
      if(TMath::IsNaN(R_EvsF_save)) { R_EvsF_save = 0; dR_EvsF_save = 0.05; }      

      h_resolEvsF->SetBinContent(i, R_EvsF);
      h_resolEvsF->SetBinError(i, dR_EvsF);

      h_resolFvsE->SetBinContent(i, R_FvsE);
      h_resolFvsE->SetBinError(i, dR_FvsE);

      h_resolTpcB->SetBinContent(i, R_TpcB);
      h_resolTpcB->SetBinError(i, dR_TpcB);

      
      h_resolution->SetBinContent(i, R_EvsF_save);
      h_resolution->SetBinError(i, dR_EvsF_save);
    }

  h_resolution->Write();
  
  // Put the bin labels on the new histograms
  for (int i = 1; i <= centBins; i++)
    {
      h_resolEvsF->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_resolFvsE->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_resolTpcB->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));

      h_centralities_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
    }
  // END RESOLUTIONS

  
  //TCanvas *canvas = new TCanvas("canvas","Canvas",900,700);
  canvas->SetGridx();
  canvas->SetGridy();
  //canvas->SetLeftMargin(0.15);
  canvas->cd();

  gStyle->SetOptStat(0);

  THStack *stack = new THStack("stack", ";Centrality (%);R_{21}");

  h_resolEvsF->SetMarkerStyle(20);
  h_resolFvsE->SetMarkerStyle(20);
  h_resolTpcB->SetMarkerStyle(20);

  h_resolEvsF->SetMarkerSize(1.5);
  h_resolFvsE->SetMarkerSize(1.5);
  h_resolTpcB->SetMarkerSize(1.5);

  h_resolEvsF->SetMarkerColor(kBlue-7);
  h_resolFvsE->SetMarkerColor(kRed-7);
  h_resolTpcB->SetMarkerColor(kGreen-3);

  h_resolEvsF->SetLineColor(kBlue-7);
  h_resolFvsE->SetLineColor(kRed-7);
  h_resolTpcB->SetLineColor(kGreen-3);

  stack->Add(h_resolEvsF);
  stack->Add(h_resolFvsE);
  stack->Add(h_resolTpcB);

  TLegend *legend = new TLegend(0.74, 0.84, 0.98, 0.98);
  legend->AddEntry(h_resolEvsF,"EPD E vs EPD F, TPC B");
  legend->AddEntry(h_resolFvsE,"EPD F vs EPD E, TPC B");
  legend->AddEntry(h_resolTpcB,"TPC B vs EPD E, EPD F");

  stack->Draw("NOSTACK E1P");
  legend->Draw();
  canvas->SaveAs("stack.png");
  canvas->Clear();

  canvas->SetLogy();
  h_centralities_flip->Draw();
  canvas->SaveAs("h_centralities_flip.png");
  canvas->Clear();
  */
  
  /*
  h_resolAvsB->Draw();
  canvas->SaveAs("h_resolAvsB.png");
  canvas->Clear();

  h_resolAvsC->Draw();
  canvas->SaveAs("h_resolAvsC.png");
  canvas->Clear();

  h_resolAvsD->Draw();
  canvas->SaveAs("h_resolAvsD.png");
  canvas->Clear();


  h_resolBvsA->Draw();
  canvas->SaveAs("h_resolBvsA.png");
  canvas->Clear();

  h_resolBvsC->Draw();
  canvas->SaveAs("h_resolBvsC.png");
  canvas->Clear();

  h_resolBvsD->Draw();
  canvas->SaveAs("h_resolBvsD.png");
  canvas->Clear();


  h_resolCvsA->Draw();
  canvas->SaveAs("h_resolCvsA.png");
  canvas->Clear();

  h_resolCvsB->Draw();
  canvas->SaveAs("h_resolCvsB.png");
  canvas->Clear();

  h_resolCvsD->Draw();
  canvas->SaveAs("h_resolCvsD.png");
  canvas->Clear();


  h_resolDvsA->Draw();
  canvas->SaveAs("h_resolDvsA.png");
  canvas->Clear();

  h_resolDvsB->Draw();
  canvas->SaveAs("h_resolDvsB.png");
  canvas->Clear();

  h_resolDvsC->Draw();
  canvas->SaveAs("h_resolDvsC.png");
  canvas->Clear();
  */
  ////////
  // END RESOLUTIONS SECTION
  ////////
  
  canvas->Close();

  //resolutionInfo_INPUT->Close();
  file->Close();
  
  delete canvas;
  delete file;

  cout << "Done!" << endl;
}

