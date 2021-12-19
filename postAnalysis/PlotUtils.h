#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include <iostream>
#include "TH1D.h"

namespace PlotUtils
{
  TH1D* flipHisto(TH1D* histo)
  {
    TH1D *h_flipped = (TH1D*)histo->Clone((TString)histo->GetName() + "_flip");
    Int_t bins = histo->GetNbinsX();
    
    Int_t j = 1;
    for (int i = bins; i >= 1; i--)
      {
	h_flipped->SetBinContent(j, histo->GetBinContent(i));
	h_flipped->SetBinError(j, histo->GetBinError(i));
	j++;
      }

    return h_flipped;
  };// End flipHisto

    // Cut off everything above 60% centrality
  TH1D* trimCentralityPlot(TH1D *histo)
  {
    TH1D *h_trimmed = (TH1D*)histo->Clone();
    Int_t oldBins = histo->GetNbinsX();
    Int_t newBins = (oldBins == 16) ? 12 : 6; // Usually 12, 6 for rebinned kaons.
    h_trimmed->SetBins(newBins, 0, 60);
    h_trimmed->GetXaxis()->SetTitle("Centrality (%)");

    for (int i = 1; i < newBins; i++)
      {
	h_trimmed->SetBinContent(i, histo->GetBinContent(i));
	h_trimmed->SetBinError(i, histo->GetBinError(i));
      }

    return h_trimmed;
  };// End trimCentralityPlot()

    // Cut off everything below y_cm = 0
  TH1D* trimRapidityPlot(TH1D* histo)
  {
    TH1D *h_trimmed = (TH1D*)histo->Clone();
    Int_t oldBins = histo->GetNbinsX();
    Int_t newBins = (oldBins == 20) ? 10 : 5; //Cut Nbins in half, 5 for rebinned kaons.
    Int_t binOffset = newBins; // old and new are offset by this amount
    h_trimmed->SetBins(newBins, 0, 1);
    h_trimmed->GetXaxis()->SetTitle("y-y_{mid}");

    for (int i = 1; i <= newBins; i++)
      {
	h_trimmed->SetBinContent(i, histo->GetBinContent(i+binOffset));
	h_trimmed->SetBinError(i, histo->GetBinError(i+binOffset));
      }

    return h_trimmed;
  };// End trimRapidityPlot()

  // Apply resolutions to standard v_n vs centrality plot (histo)
  // USE THIS BEFORE FLIPPING OR TRIMMING THE ORIGINAL HISTOGRAM!!!
  TH1D* applyResolutions(TH1D* histo, TH1D* resolutions)
  {
    TH1D *newHisto = (TH1D*)histo->Clone((TString)histo->GetName() + "_resCor");
    int bins = histo->GetNbinsX();

    for (int i = 1; i <= bins; i++)
      {
	if (resolutions->GetBinContent(i) == 0.0) continue;

	Double_t originalValue = histo->GetBinContent(i);
	Double_t originalError = histo->GetBinError(i);
	Double_t currentResolution = resolutions->GetBinContent(i);
	Double_t currentResolutionError = resolutions->GetBinError(i);

	Double_t newValue = originalValue / currentResolution;
	Double_t newError = newValue * TMath::Sqrt( TMath::Power(originalError/originalValue, 2) +
						    TMath::Power(currentResolutionError/currentResolution, 2) );

	newHisto->SetBinContent(i, newValue);
	newHisto->SetBinError(i, newError);
      }

    return newHisto;
  }// End applyResolutions()

  // Apply resolutions to 2D v_n plot with y vs centrality or p_T vs centrality (y-axis vs x-axis)
  // USE THIS BEFORE PROJECTING INTO WIDER CENTRALITY BINS!!!
  TH2D* applyResolutions2D(TH2D* histo, TH1D* resolutions)
  {
    TH2D *newHisto = (TH2D*)histo->Clone((TString)histo->GetName() + "_resCor");
    int nBinsX = histo->GetNbinsX();
    int nBinsY = histo->GetNbinsY();

    for (int xBin = 1; xBin <= nBinsX; xBin++)
      {
	if (resolutions->GetBinContent(xBin) == 0.0) continue;

	Double_t currentResolution = resolutions->GetBinContent(xBin);
	Double_t currentResolutionError = resolutions->GetBinError(xBin);
	
	for (int yBin = 1; yBin <= nBinsY; yBin++)
	  {
	    if (histo->GetBinContent(xBin, yBin) == 0.0) continue;

	    Double_t originalValue = histo->GetBinContent(xBin, yBin);
	    Double_t originalError = histo->GetBinError(xBin, yBin);

	    Double_t newValue = originalValue / currentResolution;
	    Double_t newError = newValue * TMath::Sqrt( TMath::Power(originalError/originalValue, 2) +
							TMath::Power(currentResolutionError/currentResolution, 2) );

	    newHisto->SetBinContent(xBin, yBin, newValue);
	    newHisto->SetBinError(xBin, yBin, newError);
	  }
      }

    return newHisto;
  }// End applyResolutions2D()
}

#endif
