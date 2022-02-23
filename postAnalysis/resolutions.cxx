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
  
  TProfile *h_EpdAEpdB = (TProfile*)file->Get("p_EpdAEpdB");

  TProfile *h_TpcBEpdA = (TProfile*)file->Get("p_TpcBEpdA");
  TProfile *h_TpcBEpdB = (TProfile*)file->Get("p_TpcBEpdB");

  
  Int_t centBins    = h_EpdAEpdB->GetNbinsX();
  Int_t firstCentID = h_EpdAEpdB->GetBinLowEdge(1);
  Int_t lastCentID  = h_EpdAEpdB->GetBinLowEdge(h_EpdAEpdB->GetNbinsX());

  
  TH1D *h_EpdAEpdB_flip = new TH1D("h_EpdAEpdB_flip",h_EpdAEpdB->GetTitle(),centBins,0,centBins);
  h_EpdAEpdB_flip->GetXaxis()->SetTitle((TString)h_EpdAEpdB->GetXaxis()->GetTitle()+" (%)");
  h_EpdAEpdB_flip->GetYaxis()->SetTitle(h_EpdAEpdB->GetYaxis()->GetTitle());

  TH1D *h_TpcBEpdA_flip = new TH1D("h_TpcBEpdA_flip",h_TpcBEpdA->GetTitle(),centBins,0,centBins);
  h_TpcBEpdA_flip->GetXaxis()->SetTitle((TString)h_TpcBEpdA->GetXaxis()->GetTitle()+" (%)");
  h_TpcBEpdA_flip->GetYaxis()->SetTitle(h_TpcBEpdA->GetYaxis()->GetTitle());
  TH1D *h_TpcBEpdB_flip = new TH1D("h_TpcBEpdB_flip",h_TpcBEpdB->GetTitle(),centBins,0,centBins);
  h_TpcBEpdB_flip->GetXaxis()->SetTitle((TString)h_TpcBEpdB->GetXaxis()->GetTitle()+" (%)");
  h_TpcBEpdB_flip->GetYaxis()->SetTitle(h_TpcBEpdB->GetYaxis()->GetTitle());

  TH1D *h_centralities_flip = new TH1D("h_centralities_flip",h_centralities->GetTitle(),centBins,0,centBins);
  h_centralities_flip->GetXaxis()->SetTitle((TString)h_centralities->GetXaxis()->GetTitle()+" (%)");
  h_centralities_flip->GetYaxis()->SetTitle(h_centralities->GetYaxis()->GetTitle());



  // Make the possible resolution plots
  TH1D *h_resolAvsB = new TH1D("h_resolAvsB","EPD A vs EPD B and TPC B;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);
  TH1D *h_resolBvsA = new TH1D("h_resolBvsA","EPD B vs EPD A and TPC B;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);
  TH1D *h_resolTpcB = new TH1D("h_resolTpcB","TPC B vs EPD A and EPD B;Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);
  TH1D *h_resolAvsB_2sub = new TH1D("h_resolAvsB_2sub","EPD A vs EPD B (2 sub-event);Centrality (%);R_{"+order_n_str+"1}",centBins,0,centBins);

  TH1D *h_resolutions = new TH1D("h_resolutions","EPD A Resolutions;Centrality;R_{"+order_n_str+"1}",centBins,0,centBins);
  TH2D *h2_resolutions = new TH2D("h2_resolutions","EPD A Resolutions;Centrality;y-y_{mid}",centBins,0,centBins, 20, -1, 1);
  
  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};

  std::vector<TString> newBinLabels;

  // Get list of bin labels, but flipped
  for (int i = lastCentID; i >= firstCentID; i--) { newBinLabels.push_back(centralityBins[i]); }

  // Flip the bin contents into the new histograms
  int j = 1;
  for (int i = centBins; i >= 1; i--)
    {
      h_EpdAEpdB_flip->SetBinContent(j, h_EpdAEpdB->GetBinContent(i));
      h_EpdAEpdB_flip->SetBinError(j, h_EpdAEpdB->GetBinError(i));

      h_TpcBEpdA_flip->SetBinContent(j, h_TpcBEpdA->GetBinContent(i));
      h_TpcBEpdA_flip->SetBinError(j, h_TpcBEpdA->GetBinError(i));
      h_TpcBEpdB_flip->SetBinContent(j, h_TpcBEpdB->GetBinContent(i));
      h_TpcBEpdB_flip->SetBinError(j, h_TpcBEpdB->GetBinError(i));

      h_centralities_flip->SetBinContent(j, h_centralities->GetBinContent(i));
      
      j++;
    }

  
  Double_t EpdAEpdB;
  Double_t TpcBEpdA;
  Double_t TpcBEpdB;

  Double_t dEpdAEpdB;
  Double_t dTpcBEpdA;
  Double_t dTpcBEpdB;

  Double_t R_AvsB;
  Double_t R_BvsA;
  Double_t R_TpcB;
  Double_t dR_AvsB;
  Double_t dR_BvsA;
  Double_t dR_TpcB;

  Double_t R_AvsB_2sub;
  Double_t dR_AvsB_2sub;

  Double_t EpdAEpdB_save;
  Double_t TpcBEpdA_save;
  Double_t TpcBEpdB_save;

  Double_t dEpdAEpdB_save;
  Double_t dTpcBEpdA_save;
  Double_t dTpcBEpdB_save;

  Double_t R_AvsB_save;
  Double_t dR_AvsB_save;
      
  // Fill resolution plots
  for (int i = 1; i <= centBins; i++)
    {
      EpdAEpdB_save = h_EpdAEpdB->GetBinContent(i);  //Don't use the flipped values here in the saved histogram!      
      TpcBEpdA_save = h_TpcBEpdA->GetBinContent(i);  // We need the centrality ID's in order, not the centrality percentages.
      TpcBEpdB_save = h_TpcBEpdB->GetBinContent(i);

      dEpdAEpdB_save = h_EpdAEpdB->GetBinError(i);
      dTpcBEpdA_save = h_TpcBEpdA->GetBinError(i);
      dTpcBEpdB_save = h_TpcBEpdB->GetBinError(i);

      R_AvsB_save  = TMath::Sqrt( (EpdAEpdB_save * TpcBEpdA_save) / TpcBEpdB_save );
      dR_AvsB_save = R_AvsB_save * TMath::Sqrt((dEpdAEpdB_save/(2*EpdAEpdB_save))*(dEpdAEpdB_save/(2*EpdAEpdB_save)) +
					       (dTpcBEpdA_save/(2*TpcBEpdA_save))*(dTpcBEpdA_save/(2*TpcBEpdA_save)) +
					       (dTpcBEpdB_save/(2*TpcBEpdB_save))*(dTpcBEpdB_save/(2*TpcBEpdB_save)));
	    
      
      EpdAEpdB = h_EpdAEpdB_flip->GetBinContent(i);      
      TpcBEpdA = h_TpcBEpdA_flip->GetBinContent(i);
      TpcBEpdB = h_TpcBEpdB_flip->GetBinContent(i);

      dEpdAEpdB = h_EpdAEpdB_flip->GetBinError(i);
      dTpcBEpdA = h_TpcBEpdA_flip->GetBinError(i);
      dTpcBEpdB = h_TpcBEpdB_flip->GetBinError(i);

      R_AvsB = TMath::Sqrt( (EpdAEpdB * TpcBEpdA) / TpcBEpdB );
      R_BvsA = TMath::Sqrt( (EpdAEpdB * TpcBEpdB) / TpcBEpdA );
      R_TpcB = TMath::Sqrt( (TpcBEpdA * TpcBEpdB) / EpdAEpdB );
      if (EpdAEpdB > 0)
	{
	  R_AvsB_2sub = TMath::Sqrt( EpdAEpdB );
	  dR_AvsB_2sub = (R_AvsB_2sub / 2) * ( dEpdAEpdB / EpdAEpdB);
	}

      dR_AvsB = R_AvsB * TMath::Sqrt((dEpdAEpdB/(2*EpdAEpdB))*(dEpdAEpdB/(2*EpdAEpdB)) +
				     (dTpcBEpdA/(2*TpcBEpdA))*(dTpcBEpdA/(2*TpcBEpdA)) +
				     (dTpcBEpdB/(2*TpcBEpdB))*(dTpcBEpdB/(2*TpcBEpdB)));

      dR_BvsA = R_BvsA * TMath::Sqrt((dEpdAEpdB/(2*EpdAEpdB))*(dEpdAEpdB/(2*EpdAEpdB)) +
				     (dTpcBEpdA/(2*TpcBEpdA))*(dTpcBEpdA/(2*TpcBEpdA)) +
				     (dTpcBEpdB/(2*TpcBEpdB))*(dTpcBEpdB/(2*TpcBEpdB)));

      dR_TpcB = R_TpcB * TMath::Sqrt((dTpcBEpdA/(2*TpcBEpdA))*(dTpcBEpdA/(2*TpcBEpdA)) +
				     (dTpcBEpdB/(2*TpcBEpdB))*(dTpcBEpdB/(2*TpcBEpdB)) +
				     (dEpdAEpdB/(2*EpdAEpdB))*(dEpdAEpdB/(2*EpdAEpdB)));


      if(TMath::IsNaN(R_AvsB) || dR_AvsB > 0.1) { R_AvsB = 0; dR_AvsB = 0; }
      if(TMath::IsNaN(R_BvsA) || dR_BvsA > 0.1) { R_BvsA = 0; dR_BvsA = 0; }
      if(TMath::IsNaN(R_TpcB) || dR_TpcB > 0.1) { R_TpcB = 0; dR_TpcB = 0; }
      //if(TMath::IsNaN(R_AvsB_save)) { R_AvsB_save = 0; dR_AvsB_save = 0.1; }      

      h_resolAvsB->SetBinContent(i, R_AvsB);
      h_resolAvsB->SetBinError(i, dR_AvsB);

      h_resolBvsA->SetBinContent(i, R_BvsA);
      h_resolBvsA->SetBinError(i, dR_BvsA);

      h_resolTpcB->SetBinContent(i, R_TpcB);
      h_resolTpcB->SetBinError(i, dR_TpcB);

      if (EpdAEpdB > 0)
	{
	  h_resolAvsB_2sub->SetBinContent(i, R_AvsB_2sub);
	  h_resolAvsB_2sub->SetBinError(i, dR_AvsB_2sub);
	}

      if(!TMath::IsNaN(R_AvsB_save))
	{
	  h_resolutions->SetBinContent(i, R_AvsB_save);
	  h_resolutions->SetBinError(i, dR_AvsB_save);

	  for (int j = 11; j <= 20; j++)
	    {
	      h2_resolutions->SetBinContent(i, j, R_AvsB_save);
	      h2_resolutions->SetBinError(i, j, dR_AvsB_save);
	    }
	}
    }

  h_resolutions->Write();
  h2_resolutions->Write();
  
  // Put the bin labels on the new histograms
  for (int i = 1; i <= centBins; i++)
    {
      h_resolAvsB->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_resolBvsA->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_resolTpcB->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_resolAvsB_2sub->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
      h_centralities_flip->GetXaxis()->SetBinLabel(i, newBinLabels.at(i-1));
    }
  // END RESOLUTIONS

  
  gStyle->SetOptStat(0);

  THStack *stack = new THStack("stack", ";Centrality (%);R_{"+order_n_str+"1}");

  h_resolAvsB->SetMarkerStyle(20);
  h_resolBvsA->SetMarkerStyle(20);
  h_resolTpcB->SetMarkerStyle(20);
  h_resolAvsB_2sub->SetMarkerStyle(20);

  h_resolAvsB->SetMarkerSize(1.5);
  h_resolBvsA->SetMarkerSize(1.5);
  h_resolTpcB->SetMarkerSize(1.5);
  h_resolAvsB_2sub->SetMarkerSize(1.5);

  h_resolAvsB->SetMarkerColor(kBlue-7);
  h_resolBvsA->SetMarkerColor(kRed-7);
  h_resolTpcB->SetMarkerColor(kGreen-3);

  h_resolAvsB->SetLineColor(kBlue-7);
  h_resolBvsA->SetLineColor(kRed-7);
  h_resolTpcB->SetLineColor(kGreen-3);

  stack->Add(h_resolAvsB);
  stack->Add(h_resolBvsA);
  stack->Add(h_resolTpcB);

  TLegend *legend = new TLegend(0.7, 0.75, 0.9, 0.9);
  legend->AddEntry(h_resolAvsB,"EPD A");
  legend->AddEntry(h_resolBvsA,"EPD B");
  legend->AddEntry(h_resolTpcB,"TPC B");

  canvas->SetTicks();
  stack->Draw("NOSTACK E1P");
  legend->Draw();
  canvas->SaveAs(jobID + "_resolutions.png");
  canvas->Clear();

  TLegend *legend2 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legend2->AddEntry(h_resolAvsB,"EPD A vs EPD B, TPC B");

  h_resolAvsB->SetMaximum(0.7);
  h_resolAvsB->SetTitle("");
  h_resolAvsB->Draw("E1P");
  legend2->Draw();
  canvas->SaveAs(jobID + "_resolutionAonly.png");
  canvas->Clear();

  h_resolAvsB_2sub->SetMaximum(0.7);
  h_resolAvsB_2sub->SetTitle("");
  h_resolAvsB_2sub->Draw("E1P");
  //legend2->Draw();
  canvas->SaveAs(jobID + "_resolution2SubOnly.png");
  canvas->Clear();

  canvas->SetTicks(0);
  canvas->SetLogy();
  h_centralities_flip->Draw();
  canvas->SaveAs(jobID + "_h_centralities_flip.png");
  canvas->Clear();

  

  
/*  
  canvas->SetLogy(0);
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

  resolutionInfo_INPUT->Close();
  file->Close();
}
