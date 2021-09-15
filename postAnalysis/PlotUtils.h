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
}

#endif
