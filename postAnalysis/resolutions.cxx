void resolutions(TString jobID, TString order_n_str)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  TFile *resolutionInfo_INPUT = new TFile("resolutionInfo_INPUT.root", "RECREATE");

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 875, 675);
  canvas->SetGridx();
  canvas->SetGridy();
  //canvas->SetLeftMargin(0.15);
  canvas->cd();

  TH1D *h_centralities = (TH1D*)file->Get("h_centralities");
  
  TProfile *h_EpdEEpdF = (TProfile*)file->Get("p_EpdEEpdF");

  TProfile *h_TpcBEpdE = (TProfile*)file->Get("p_TpcBEpdE");
  TProfile *h_TpcBEpdF = (TProfile*)file->Get("p_TpcBEpdF");

  
  Int_t centBins    = h_EpdEEpdF->GetNbinsX();
  Int_t firstCentID = h_EpdEEpdF->GetBinLowEdge(1);
  Int_t lastCentID  = h_EpdEEpdF->GetBinLowEdge(h_EpdEEpdF->GetNbinsX());

  
  TH1D *h_EpdEEpdF_flip = new TH1D("h_EpdEEpdF_flip",h_EpdEEpdF->GetTitle(),centBins,0,centBins);
  h_EpdEEpdF_flip->GetXaxis()->SetTitle((TString)h_EpdEEpdF->GetXaxis()->GetTitle()+" (%)");
  h_EpdEEpdF_flip->GetYaxis()->SetTitle(h_EpdEEpdF->GetYaxis()->GetTitle());

  TH1D *h_TpcBEpdE_flip = new TH1D("h_TpcBEpdE_flip",h_TpcBEpdE->GetTitle(),centBins,0,centBins);
  h_TpcBEpdE_flip->GetXaxis()->SetTitle((TString)h_TpcBEpdE->GetXaxis()->GetTitle()+" (%)");
  h_TpcBEpdE_flip->GetYaxis()->SetTitle(h_TpcBEpdE->GetYaxis()->GetTitle());
  TH1D *h_TpcBEpdF_flip = new TH1D("h_TpcBEpdF_flip",h_TpcBEpdF->GetTitle(),centBins,0,centBins);
  h_TpcBEpdF_flip->GetXaxis()->SetTitle((TString)h_TpcBEpdF->GetXaxis()->GetTitle()+" (%)");
  h_TpcBEpdF_flip->GetYaxis()->SetTitle(h_TpcBEpdF->GetYaxis()->GetTitle());

  TH1D *h_centralities_flip = new TH1D("h_centralities_flip",h_centralities->GetTitle(),centBins,0,centBins);
  h_centralities_flip->GetXaxis()->SetTitle((TString)h_centralities->GetXaxis()->GetTitle()+" (%)");
  h_centralities_flip->GetYaxis()->SetTitle(h_centralities->GetYaxis()->GetTitle());



  // Make the possible resolution plots
  TH1D *h_resolEvsF = new TH1D("h_resolEvsF","EPD E vs EPD F and TPC B;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);
  TH1D *h_resolFvsE = new TH1D("h_resolFvsE","EPD F vs EPD E and TPC B;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);
  TH1D *h_resolTpcB = new TH1D("h_resolTpcB","TPC B vs EPD E and EPD F;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);

  TH1D *h_resolutions = new TH1D("h_resolutions","EPD E Resolutions;Centrality;R_{"+order_n_str+"1}",centBins,0,centBins);
  TH2D *h2_resolutions = new TH2D("h2_resolutions","EPD E Resolutions;Centrality;y-y_{mid}",centBins,0,centBins, 20, -1, 1);
  
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
      EpdEEpdF_save = h_EpdEEpdF->GetBinContent(i);  //Don't use the flipped values here in the saved histogram!      
      TpcBEpdE_save = h_TpcBEpdE->GetBinContent(i);  // We need the centrality ID's in order, not the centrality percentages.
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

      if(TMath::IsNaN(R_EvsF) || dR_EvsF > 0.1) { R_EvsF = 0; dR_EvsF = 0; }
      if(TMath::IsNaN(R_FvsE) || dR_FvsE > 0.1) { R_FvsE = 0; dR_FvsE = 0; }
      if(TMath::IsNaN(R_TpcB) || dR_TpcB > 0.1) { R_TpcB = 0; dR_TpcB = 0; }
      //if(TMath::IsNaN(R_EvsF_save)) { R_EvsF_save = 0; dR_EvsF_save = 0.1; }      

      h_resolEvsF->SetBinContent(i, R_EvsF);
      h_resolEvsF->SetBinError(i, dR_EvsF);

      h_resolFvsE->SetBinContent(i, R_FvsE);
      h_resolFvsE->SetBinError(i, dR_FvsE);

      h_resolTpcB->SetBinContent(i, R_TpcB);
      h_resolTpcB->SetBinError(i, dR_TpcB);

      if(!TMath::IsNaN(R_EvsF_save))
	{
	  h_resolutions->SetBinContent(i, R_EvsF_save);
	  h_resolutions->SetBinError(i, dR_EvsF_save);

	  for (int j = 11; j <= 20; j++)
	    {
	      h2_resolutions->SetBinContent(i, j, R_EvsF_save);
	      h2_resolutions->SetBinError(i, j, dR_EvsF_save);
	    }
	}
    }

  h_resolutions->Write();
  h2_resolutions->Write();
  
  // Put the bin labels on the new histograms
  for (int i = 1; i <= centBins; i++)
    {
      h_resolEvsF->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_resolFvsE->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_resolTpcB->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));

      h_centralities_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
    }
  // END RESOLUTIONS

  
  gStyle->SetOptStat(0);

  THStack *stack = new THStack("stack", ";Centrality (%);R_{"+order_n_str+"1}");

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

  TLegend *legend = new TLegend(0.65, 0.75, 0.9, 0.9);
  legend->AddEntry(h_resolEvsF,"EPD E vs EPD F, TPC B");
  legend->AddEntry(h_resolFvsE,"EPD F vs EPD E, TPC B");
  legend->AddEntry(h_resolTpcB,"TPC B vs EPD E, EPD F");

  canvas->SetTicks();
  stack->Draw("NOSTACK E1P");
  legend->Draw();
  canvas->SaveAs(jobID + "_resolutions.png");
  canvas->Clear();

  TLegend *legend2 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legend2->AddEntry(h_resolEvsF,"EPD E vs EPD F, TPC B");

  h_resolEvsF->SetMaximum(0.7);
  h_resolEvsF->SetTitle("");
  h_resolEvsF->Draw("E1P");
  legend2->Draw();
  canvas->SaveAs(jobID + "_resolutionEonly.png");
  canvas->Clear();

  canvas->SetTicks(0);
  canvas->SetLogy();
  h_centralities_flip->Draw();
  canvas->SaveAs(jobID + "_h_centralities_flip.png");
  canvas->Clear();

  

  
/*  
  canvas->SetLogy(0);
  h_resolEvsF->Draw();
  canvas->SaveAs("h_resolEvsF.png");
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

  resolutionInfo_INPUT->Close();
  file->Close();
}
