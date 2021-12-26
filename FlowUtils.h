#ifndef FLOWUTILS_H
#define FLOWUTILS_H

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TMath.h"


namespace FlowUtils
{
  const Int_t I_BAD_VALUE    = -999;
  const Double_t D_BAD_VALUE = -999.0;

  // Custom type to hold info for each good particle
  struct Particle
  {
    Double_t phi;
    Double_t eta;
    Double_t pT;
    Double_t weight;
    Double_t rapidity;

    Bool_t ppTag;
    Bool_t pmTag;
    Bool_t kpTag;
    Bool_t kmTag;
    Bool_t prTag;
    Bool_t deTag;
    Bool_t trTag;

    Bool_t isInTpcA;
    Bool_t isInTpcB;
    Bool_t isInEpdA;
    Bool_t isInEpdB;

    void reset()
    {
      phi = D_BAD_VALUE;
      eta = D_BAD_VALUE;
      weight = 0;
      rapidity = D_BAD_VALUE;

      ppTag = false;
      pmTag = false;
      kpTag = false;
      kmTag = false;
      prTag = false;
      deTag = false;
      trTag = false;

      isInTpcA = false;
      isInTpcB = false;
      isInEpdA = false;
      isInEpdB = false;
    }
  };




  // Custom type to hold important info for every good event
  struct Event
  {
    std::vector<Particle> tpcParticles;
    std::vector<Particle> epdParticles;

    Double_t XnTpc;
    Double_t YnTpc;
    Double_t psiTpc; 
    Double_t XnTpcA;
    Double_t YnTpcA;
    Double_t psiTpcA;       // Overall EP angle without removing autocorrelations
    Double_t XnTpcB;
    Double_t YnTpcB;
    Double_t psiTpcB;
    Double_t XnEpd;
    Double_t YnEpd;
    Double_t psiEpd;
    Double_t XnEpdA;
    Double_t YnEpdA;
    Double_t psiEpdA;
    Double_t XnEpdB;
    Double_t YnEpdB;
    Double_t psiEpdB;

    Int_t nTracksTpc;
    Int_t nTracksTpcA;      // Number of GOOD tracks in the sub-event
    Int_t nTracksTpcB;
    Int_t nHitsEpd;
    Int_t nHitsEpdA;
    Int_t nHitsEpdB;

    bool badEvent;      // Flag for marking events to ignore
    Int_t centID;
    Int_t primTracks;   // Number of primary tracks before track cuts (used for centrality)

    void reset()//Reset all values in the struct to reuse
    {
      std::vector<Particle>().swap(tpcParticles);
      std::vector<Particle>().swap(epdParticles);

      XnTpc = 0;
      YnTpc = 0;
      psiTpc = D_BAD_VALUE;        //Just some number to use that is out of bounds
      XnTpcA = 0;
      YnTpcA = 0;
      psiTpcA = D_BAD_VALUE;
      XnTpcB = 0;
      YnTpcB = 0;
      psiTpcB = D_BAD_VALUE;
      XnEpd = 0;
      YnEpd = 0;
      psiEpd = D_BAD_VALUE;
      XnEpdA = 0;
      YnEpdA = 0;
      psiEpdA = D_BAD_VALUE;
      XnEpdB = 0;
      YnEpdB = 0;
      psiEpdB = D_BAD_VALUE;

      nTracksTpc = 0;
      nTracksTpcA = 0;
      nTracksTpcB = 0;
      nHitsEpd = 0;
      nHitsEpdA = 0;
      nHitsEpdB = 0;

      badEvent  = false;
      primTracks = 0;
      centID = I_BAD_VALUE;
    }
  };



  
  ////////
  //   Calculates the event plane angle in every subevent.
  ////////
  void getAllPsi(Event &eventInfo, Double_t order_m)
  {
    eventInfo.psiTpc  = TMath::ATan2(eventInfo.YnTpc,  eventInfo.XnTpc)  / order_m;
    eventInfo.psiTpcA = TMath::ATan2(eventInfo.YnTpcA, eventInfo.XnTpcA) / order_m;
    eventInfo.psiTpcB = TMath::ATan2(eventInfo.YnTpcB, eventInfo.XnTpcB) / order_m;
    eventInfo.psiEpd  = TMath::ATan2(eventInfo.YnEpd,  eventInfo.XnEpd)  / order_m;
    eventInfo.psiEpdA = TMath::ATan2(eventInfo.YnEpdA, eventInfo.XnEpdA) / order_m;
    eventInfo.psiEpdB = TMath::ATan2(eventInfo.YnEpdB, eventInfo.XnEpdB) / order_m;
  }

  ////////
  //   Moves event plane angles back into the -pi to pi period.
  ////////
  Double_t angleShift(Double_t angle, Int_t order)
  {
    if (angle < -TMath::Pi()/(Double_t)order) { angle += TMath::TwoPi()/(Double_t)order; }
    else if (angle >  TMath::Pi()/(Double_t)order) { angle -= TMath::TwoPi()/(Double_t)order; }
    return angle;
  }

  ////////
  //   Moves all event plane angles into the -pi to pi period.
  ////////
  void setAllPeriods(Event &eventInfo, Double_t order_m)
  {
    eventInfo.psiTpc  = angleShift(eventInfo.psiTpc,  order_m);
    eventInfo.psiTpcA = angleShift(eventInfo.psiTpcA, order_m);
    eventInfo.psiTpcB = angleShift(eventInfo.psiTpcB, order_m);
    eventInfo.psiEpd  = angleShift(eventInfo.psiEpd,  order_m);
    eventInfo.psiEpdA = angleShift(eventInfo.psiEpdA, order_m);
    eventInfo.psiEpdB = angleShift(eventInfo.psiEpdB, order_m);
  }

  ////////
  //   Checks if an event has any flow vectors equal to zero. Updates the event's member variable "badEvent".
  ////////
  void checkZeroQ(Event &event)
  {
    if (event.XnTpc == 0 && event.YnTpc == 0) { event.badEvent = true; }
    //else if (event.XnTpcA == 0 && event.YnTpcA == 0) { event.badEvent = true; }
    if (event.XnTpcB == 0 && event.YnTpcB == 0) { event.badEvent = true; }
    else if (event.XnEpd  == 0 && event.YnEpd  == 0) { event.badEvent = true; }
    else if (event.XnEpdA == 0 && event.YnEpdA == 0) { event.badEvent = true; }
    else if (event.XnEpdB == 0 && event.YnEpdB == 0) { event.badEvent = true; }
  }


  ////////
  //   Using px, py, pz, and rest mass, return rapidity
  ////////
  Double_t rapidity(Double_t px, Double_t py, Double_t pz, Double_t mass)
  {
    Double_t rapidity, energy, momentum = 0;
    momentum = TMath::Sqrt(px*px + py*py + pz*pz);
    energy   = TMath::Sqrt(momentum*momentum + mass*mass);
    rapidity = TMath::ATanH(pz/energy);
    return rapidity;
  }

  Double_t pseudorapidity(Double_t px, Double_t py, Double_t pz)
  {
    Double_t momentum = TMath::Sqrt(px*px + py*py + pz*pz);
    return 0.5 * TMath::Log( (momentum + pz)/(momentum - pz) );
  }

  ////////
  //   Using px, py, pz, and rest mass, return transverse mass
  ////////
  Double_t transMass(Double_t px, Double_t py, Double_t mass) 
  { return TMath::Sqrt(mass*mass + px*px + py*py); }


  Double_t transMomentum(Double_t px, Double_t py)
  { return TMath::Sqrt(px*px + py*py); }

 
  Double_t totalMomentum(Double_t px, Double_t py, Double_t pz)
  {
    return TMath::Sqrt(px*px + py*py + pz*pz);
  }


  Double_t phi(Double_t px, Double_t py)
  { return TMath::ATan2(py, px); }

  Int_t epdSector(Short_t tileID)
  { return TMath::Abs( tileID/100 ); }


  Int_t epdRow(Short_t tileID)
  {
    Short_t absID = TMath::Abs(tileID);
    Short_t tileNum = absID % 100;
    Int_t rowNum = 0;
    
    if (tileNum == 1) rowNum = 1;
    else if (tileNum == 2 || tileNum == 3) rowNum = 2;
    else if (tileNum == 4 || tileNum == 5) rowNum = 3;
    else if (tileNum == 6 || tileNum == 7) rowNum = 4;
    else if (tileNum == 8 || tileNum == 9) rowNum = 5;
    else if (tileNum == 10 || tileNum == 11) rowNum = 6;
    else if (tileNum == 12 || tileNum == 13) rowNum = 7;
    else if (tileNum == 14 || tileNum == 15) rowNum = 8;
    else if (tileNum == 16 || tileNum == 17) rowNum = 9;
    else if (tileNum == 18 || tileNum == 19) rowNum = 10;
    else if (tileNum == 20 || tileNum == 21) rowNum = 11;
    else if (tileNum == 22 || tileNum == 23) rowNum = 12;
    else if (tileNum == 24 || tileNum == 25) rowNum = 13;
    else if (tileNum == 26 || tileNum == 27) rowNum = 14;
    else if (tileNum == 28 || tileNum == 29) rowNum = 15;
    else if (tileNum == 30 || tileNum == 31) rowNum = 16;

    return rowNum;
  }

  ////////
  //   Using px, py, pz, and rest mass, fill histograms raw of dN/dy,
  // and dN/dmT (shifted left by m0), and a 2D histogram of mT-m0 vs y.
  ////////
  void fillRawSpect(Double_t px, Double_t py, Double_t pz, Double_t mass, TH1D *dndy, TH1D *dndm, TH2D *MvsY)
  {
    Double_t y  = rapidity(px, py, pz, mass);
    Double_t mT = transMass(px, py, mass);
    Double_t M  = mT - mass;
    dndy->Fill(y);
    dndm->Fill(M);
    MvsY->Fill(y,M, 1/(TMath::TwoPi() * mT));
  }


  Double_t getTpcEff(Double_t yCM, Double_t pT, TH2D *h2_ratio)
  {
    Int_t xBin = h2_ratio->GetXaxis()->FindBin(yCM);
    Int_t yBin = h2_ratio->GetYaxis()->FindBin(pT);
    Double_t efficiency = h2_ratio->GetBinContent(xBin, yBin);
    return (efficiency == 0) ? -1 : efficiency;
  }


  ////////
  //   Recenters the flow vectors of every subevent region in an event using the averages 
  //  found over all events in the full dataset and then recalculates event plane angles.
  ////////
  void recenterQ(Event &eventInfo, TFile *correctionInputFile, Double_t order_m)
  {
    TH1D *h_XnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_XnTpc");
    TH1D *h_XnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcA");
    TH1D *h_XnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_XnTpcB");
    TH1D *h_XnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_XnEpd");
    TH1D *h_XnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdA");
    TH1D *h_XnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_XnEpdB");

    TH1D *h_YnTpc_INPUT  = (TH1D*)correctionInputFile->Get("h_YnTpc");
    TH1D *h_YnTpcA_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcA");
    TH1D *h_YnTpcB_INPUT = (TH1D*)correctionInputFile->Get("h_YnTpcB");
    TH1D *h_YnEpd_INPUT  = (TH1D*)correctionInputFile->Get("h_YnEpd");
    TH1D *h_YnEpdA_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdA");
    TH1D *h_YnEpdB_INPUT = (TH1D*)correctionInputFile->Get("h_YnEpdB");

    Double_t d_XnTpc_Avg  = h_XnTpc_INPUT->GetMean();
    Double_t d_XnTpcA_Avg = h_XnTpcA_INPUT->GetMean();
    Double_t d_XnTpcB_Avg = h_XnTpcB_INPUT->GetMean();
    Double_t d_XnEpd_Avg  = h_XnEpd_INPUT->GetMean();
    Double_t d_XnEpdA_Avg = h_XnEpdA_INPUT->GetMean();
    Double_t d_XnEpdB_Avg = h_XnEpdB_INPUT->GetMean();

    Double_t d_YnTpc_Avg  = h_YnTpc_INPUT->GetMean();
    Double_t d_YnTpcA_Avg = h_YnTpcA_INPUT->GetMean();
    Double_t d_YnTpcB_Avg = h_YnTpcB_INPUT->GetMean();
    Double_t d_YnEpd_Avg  = h_YnEpd_INPUT->GetMean();
    Double_t d_YnEpdA_Avg = h_YnEpdA_INPUT->GetMean();
    Double_t d_YnEpdB_Avg = h_YnEpdB_INPUT->GetMean();


    eventInfo.XnTpc  -= d_XnTpc_Avg;
    eventInfo.XnTpcA -= d_XnTpcA_Avg;
    eventInfo.XnTpcB -= d_XnTpcB_Avg;
    eventInfo.XnEpd  -= d_XnEpd_Avg;
    eventInfo.XnEpdA -= d_XnEpdA_Avg;
    eventInfo.XnEpdB -= d_XnEpdB_Avg;

    eventInfo.YnTpc  -= d_YnTpc_Avg;
    eventInfo.YnTpcA -= d_YnTpcA_Avg;
    eventInfo.YnTpcB -= d_YnTpcB_Avg;
    eventInfo.YnEpd  -= d_YnEpd_Avg;
    eventInfo.YnEpdA -= d_YnEpdA_Avg;
    eventInfo.YnEpdB -= d_YnEpdB_Avg;

    checkZeroQ(eventInfo);

    getAllPsi(eventInfo, order_m);
    setAllPeriods(eventInfo, order_m);
  }// End recenterQ()


  ////////
  //   Performs the event-by-event shifting described in the Poskanzer paper to flatten the event 
  //  plane angle distributions of each subevent.
  ////////
  void shiftPsi(Event &eventInfo, TFile *correctionInputFile, Double_t order_m, Int_t shiftTerms)
  {
    TProfile *p_sinAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsTpc");
    TProfile *p_cosAvgsTpc_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsTpc");
    TProfile *p_sinAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcA");
    TProfile *p_cosAvgsTpcA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcA");
    TProfile *p_sinAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsTpcB");
    TProfile *p_cosAvgsTpcB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsTpcB");
    TProfile *p_sinAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_sinAvgsEpd");
    TProfile *p_cosAvgsEpd_INPUT  = (TProfile*)correctionInputFile->Get("p_cosAvgsEpd");
    TProfile *p_sinAvgsEpdA_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdA");
    TProfile *p_cosAvgsEpdA_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdA");
    TProfile *p_sinAvgsEpdB_INPUT = (TProfile*)correctionInputFile->Get("p_sinAvgsEpdB");
    TProfile *p_cosAvgsEpdB_INPUT = (TProfile*)correctionInputFile->Get("p_cosAvgsEpdB");


    // Get corrected event plane angles //

    Double_t psiTpc_delta  = 0;
    Double_t psiTpcA_delta = 0;
    Double_t psiTpcB_delta = 0;
    Double_t psiEpd_delta  = 0;
    Double_t psiEpdA_delta = 0;
    Double_t psiEpdB_delta = 0;

    Double_t jthSinAvg_Tpc  = 0;
    Double_t jthCosAvg_Tpc  = 0;
    Double_t jthSinAvg_TpcA = 0;
    Double_t jthCosAvg_TpcA = 0;
    Double_t jthSinAvg_TpcB = 0;
    Double_t jthCosAvg_TpcB = 0;
    Double_t jthSinAvg_Epd  = 0;
    Double_t jthCosAvg_Epd  = 0;
    Double_t jthSinAvg_EpdA = 0;
    Double_t jthCosAvg_EpdA = 0;
    Double_t jthSinAvg_EpdB = 0;
    Double_t jthCosAvg_EpdB = 0;


    for (Int_t j = 1; j <= shiftTerms; j++)    // Build the correction sums
      {
	jthSinAvg_Tpc  = p_sinAvgsTpc_INPUT->GetBinContent(j);
	jthCosAvg_Tpc  = p_cosAvgsTpc_INPUT->GetBinContent(j);
	jthSinAvg_TpcA = p_sinAvgsTpcA_INPUT->GetBinContent(j);
	jthCosAvg_TpcA = p_cosAvgsTpcA_INPUT->GetBinContent(j);
	jthSinAvg_TpcB = p_sinAvgsTpcB_INPUT->GetBinContent(j);
	jthCosAvg_TpcB = p_cosAvgsTpcB_INPUT->GetBinContent(j);
	jthSinAvg_Epd  = p_sinAvgsEpd_INPUT->GetBinContent(j);
	jthCosAvg_Epd  = p_cosAvgsEpd_INPUT->GetBinContent(j);
	jthSinAvg_EpdA = p_sinAvgsEpdA_INPUT->GetBinContent(j);
	jthCosAvg_EpdA = p_cosAvgsEpdA_INPUT->GetBinContent(j);
	jthSinAvg_EpdB = p_sinAvgsEpdB_INPUT->GetBinContent(j);
	jthCosAvg_EpdB = p_cosAvgsEpdB_INPUT->GetBinContent(j);

	psiTpc_delta  += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_Tpc * TMath::Cos((Double_t)j * order_m * eventInfo.psiTpc) 
							+jthCosAvg_Tpc * TMath::Sin((Double_t)j * order_m * eventInfo.psiTpc));
	psiTpcA_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_TpcA * TMath::Cos((Double_t)j * order_m * eventInfo.psiTpcA) 
							+jthCosAvg_TpcA * TMath::Sin((Double_t)j * order_m * eventInfo.psiTpcA));
	psiTpcB_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_TpcB * TMath::Cos((Double_t)j * order_m * eventInfo.psiTpcB) 
							+jthCosAvg_TpcB * TMath::Sin((Double_t)j * order_m * eventInfo.psiTpcB));
	psiEpd_delta  += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_Epd * TMath::Cos((Double_t)j * order_m * eventInfo.psiEpd)
							+jthCosAvg_Epd * TMath::Sin((Double_t)j * order_m * eventInfo.psiEpd));
	psiEpdA_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_EpdA * TMath::Cos((Double_t)j * order_m * eventInfo.psiEpdA)
							+jthCosAvg_EpdA * TMath::Sin((Double_t)j * order_m * eventInfo.psiEpdA));
	psiEpdB_delta += (2.0/((Double_t)j*order_m)) * (-jthSinAvg_EpdB * TMath::Cos((Double_t)j * order_m * eventInfo.psiEpdB)
							+jthCosAvg_EpdB * TMath::Sin((Double_t)j * order_m * eventInfo.psiEpdB));
      }

    // Shift event plane angles
    eventInfo.psiTpc  += psiTpc_delta;
    eventInfo.psiTpcA += psiTpcA_delta;
    eventInfo.psiTpcB += psiTpcB_delta;
    eventInfo.psiEpd  += psiEpd_delta;
    eventInfo.psiEpdA += psiEpdA_delta;
    eventInfo.psiEpdB += psiEpdB_delta;

    // Keep angles in the correct period
    setAllPeriods(eventInfo, order_m);
  }// End shiftPsi()



}

#endif
