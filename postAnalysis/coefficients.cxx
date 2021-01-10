void coefficients(TString jobID, TString order_n_str)
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

  TProfile *p_vn_yCM_00to10_pp = (TProfile*)file->Get("p_vn_yCM_00to10_pp");
  TProfile *p_vn_yCM_10to40_pp = (TProfile*)file->Get("p_vn_yCM_10to40_pp");
  TProfile *p_vn_yCM_40to60_pp = (TProfile*)file->Get("p_vn_yCM_40to60_pp");
  TProfile *p_vn_yCM_00to60_pp = (TProfile*)file->Get("p_vn_yCM_00to60_pp");

  TProfile *p_vn_yCM_00to10_pm = (TProfile*)file->Get("p_vn_yCM_00to10_pm");
  TProfile *p_vn_yCM_10to40_pm = (TProfile*)file->Get("p_vn_yCM_10to40_pm");
  TProfile *p_vn_yCM_40to60_pm = (TProfile*)file->Get("p_vn_yCM_40to60_pm");
  TProfile *p_vn_yCM_00to60_pm = (TProfile*)file->Get("p_vn_yCM_00to60_pm");

  TProfile *p_vn_yCM_00to10_kp = (TProfile*)file->Get("p_vn_yCM_00to10_kp");
  TProfile *p_vn_yCM_10to40_kp = (TProfile*)file->Get("p_vn_yCM_10to40_kp");
  TProfile *p_vn_yCM_40to60_kp = (TProfile*)file->Get("p_vn_yCM_40to60_kp");
  TProfile *p_vn_yCM_00to60_kp = (TProfile*)file->Get("p_vn_yCM_00to60_kp");

  TProfile *p_vn_yCM_00to10_km = (TProfile*)file->Get("p_vn_yCM_00to10_km");
  TProfile *p_vn_yCM_10to40_km = (TProfile*)file->Get("p_vn_yCM_10to40_km");
  TProfile *p_vn_yCM_40to60_km = (TProfile*)file->Get("p_vn_yCM_40to60_km");
  TProfile *p_vn_yCM_00to60_km = (TProfile*)file->Get("p_vn_yCM_00to60_km");

  TProfile *p_vn_yCM_00to10_pr = (TProfile*)file->Get("p_vn_yCM_00to10_pr");
  TProfile *p_vn_yCM_10to40_pr = (TProfile*)file->Get("p_vn_yCM_10to40_pr");
  TProfile *p_vn_yCM_40to60_pr = (TProfile*)file->Get("p_vn_yCM_40to60_pr");
  TProfile *p_vn_yCM_00to60_pr = (TProfile*)file->Get("p_vn_yCM_00to60_pr");




  TH1D *h_vn_EpdE = new TH1D("h_vn_EpdE", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_EpdF = new TH1D("h_vn_EpdF", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_TpcB = new TH1D("h_vn_TpcB", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_pp = new TH1D("h_vn_pp", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_pm = new TH1D("h_vn_pm", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_kp = new TH1D("h_vn_kp", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_km = new TH1D("h_vn_km", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_pr = new TH1D("h_vn_pr", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);


  TH1D *h_vn_yCM_00to10_pp = new TH1D("h_vn_yCM_00to10_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pp = new TH1D("h_vn_yCM_10to40_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pp = new TH1D("h_vn_yCM_40to60_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pp = new TH1D("h_vn_yCM_00to60_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pm = new TH1D("h_vn_yCM_00to10_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pm = new TH1D("h_vn_yCM_10to40_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pm = new TH1D("h_vn_yCM_40to60_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pm = new TH1D("h_vn_yCM_00to60_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_kp = new TH1D("h_vn_yCM_00to10_kp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_kp = new TH1D("h_vn_yCM_10to40_kp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_kp = new TH1D("h_vn_yCM_40to60_kp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_kp = new TH1D("h_vn_yCM_00to60_kp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_km = new TH1D("h_vn_yCM_00to10_km", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_km = new TH1D("h_vn_yCM_10to40_km", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_km = new TH1D("h_vn_yCM_40to60_km", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_km = new TH1D("h_vn_yCM_00to60_km", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pr = new TH1D("h_vn_yCM_00to10_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr = new TH1D("h_vn_yCM_10to40_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr = new TH1D("h_vn_yCM_40to60_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pr = new TH1D("h_vn_yCM_00to60_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  // Convert profiles to histograms
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

  // Convert profiles to histograms
  for (int i = 1; i <= 20; i++)
    {
      h_vn_yCM_00to10_pp->SetBinContent(i, p_vn_yCM_00to10_pp->GetBinContent(i));
      h_vn_yCM_00to10_pp->SetBinError(i, p_vn_yCM_00to10_pp->GetBinError(i));
      h_vn_yCM_10to40_pp->SetBinContent(i, p_vn_yCM_10to40_pp->GetBinContent(i));
      h_vn_yCM_10to40_pp->SetBinError(i, p_vn_yCM_10to40_pp->GetBinError(i));
      h_vn_yCM_40to60_pp->SetBinContent(i, p_vn_yCM_40to60_pp->GetBinContent(i));
      h_vn_yCM_40to60_pp->SetBinError(i, p_vn_yCM_40to60_pp->GetBinError(i));
      h_vn_yCM_00to60_pp->SetBinContent(i, p_vn_yCM_00to60_pp->GetBinContent(i));
      h_vn_yCM_00to60_pp->SetBinError(i, p_vn_yCM_00to60_pp->GetBinError(i));

      h_vn_yCM_00to10_pm->SetBinContent(i, p_vn_yCM_00to10_pm->GetBinContent(i));
      h_vn_yCM_00to10_pm->SetBinError(i, p_vn_yCM_00to10_pm->GetBinError(i));
      h_vn_yCM_10to40_pm->SetBinContent(i, p_vn_yCM_10to40_pm->GetBinContent(i));
      h_vn_yCM_10to40_pm->SetBinError(i, p_vn_yCM_10to40_pm->GetBinError(i));
      h_vn_yCM_40to60_pm->SetBinContent(i, p_vn_yCM_40to60_pm->GetBinContent(i));
      h_vn_yCM_40to60_pm->SetBinError(i, p_vn_yCM_40to60_pm->GetBinError(i));
      h_vn_yCM_00to60_pm->SetBinContent(i, p_vn_yCM_00to60_pm->GetBinContent(i));
      h_vn_yCM_00to60_pm->SetBinError(i, p_vn_yCM_00to60_pm->GetBinError(i));

      h_vn_yCM_00to10_kp->SetBinContent(i, p_vn_yCM_00to10_kp->GetBinContent(i));
      h_vn_yCM_00to10_kp->SetBinError(i, p_vn_yCM_00to10_kp->GetBinError(i));
      h_vn_yCM_10to40_kp->SetBinContent(i, p_vn_yCM_10to40_kp->GetBinContent(i));
      h_vn_yCM_10to40_kp->SetBinError(i, p_vn_yCM_10to40_kp->GetBinError(i));
      h_vn_yCM_40to60_kp->SetBinContent(i, p_vn_yCM_40to60_kp->GetBinContent(i));
      h_vn_yCM_40to60_kp->SetBinError(i, p_vn_yCM_40to60_kp->GetBinError(i));
      h_vn_yCM_00to60_kp->SetBinContent(i, p_vn_yCM_00to60_kp->GetBinContent(i));
      h_vn_yCM_00to60_kp->SetBinError(i, p_vn_yCM_00to60_kp->GetBinError(i));

      h_vn_yCM_00to10_km->SetBinContent(i, p_vn_yCM_00to10_km->GetBinContent(i));
      h_vn_yCM_00to10_km->SetBinError(i, p_vn_yCM_00to10_km->GetBinError(i));
      h_vn_yCM_10to40_km->SetBinContent(i, p_vn_yCM_10to40_km->GetBinContent(i));
      h_vn_yCM_10to40_km->SetBinError(i, p_vn_yCM_10to40_km->GetBinError(i));
      h_vn_yCM_40to60_km->SetBinContent(i, p_vn_yCM_40to60_km->GetBinContent(i));
      h_vn_yCM_40to60_km->SetBinError(i, p_vn_yCM_40to60_km->GetBinError(i));
      h_vn_yCM_00to60_km->SetBinContent(i, p_vn_yCM_00to60_km->GetBinContent(i));
      h_vn_yCM_00to60_km->SetBinError(i, p_vn_yCM_00to60_km->GetBinError(i));

      h_vn_yCM_00to10_pr->SetBinContent(i, p_vn_yCM_00to10_pr->GetBinContent(i));
      h_vn_yCM_00to10_pr->SetBinError(i, p_vn_yCM_00to10_pr->GetBinError(i));
      h_vn_yCM_10to40_pr->SetBinContent(i, p_vn_yCM_10to40_pr->GetBinContent(i));
      h_vn_yCM_10to40_pr->SetBinError(i, p_vn_yCM_10to40_pr->GetBinError(i));
      h_vn_yCM_40to60_pr->SetBinContent(i, p_vn_yCM_40to60_pr->GetBinContent(i));
      h_vn_yCM_40to60_pr->SetBinError(i, p_vn_yCM_40to60_pr->GetBinError(i));
      h_vn_yCM_00to60_pr->SetBinContent(i, p_vn_yCM_00to60_pr->GetBinContent(i));
      h_vn_yCM_00to60_pr->SetBinError(i, p_vn_yCM_00to60_pr->GetBinError(i));
    }

  /*
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
  */

  Int_t centBins    = h_vn_TpcB->GetNbinsX();
  Int_t firstCentID = h_vn_TpcB->GetBinLowEdge(1);
  Int_t lastCentID  = h_vn_TpcB->GetBinLowEdge(h_vn_TpcB->GetNbinsX());

  TH1D *h_vn_EpdE_flip = new TH1D("h_vn_EpdE_flip",h_vn_EpdE->GetTitle(),centBins,0,centBins);
  h_vn_EpdE_flip->GetXaxis()->SetTitle((TString)h_vn_EpdE->GetXaxis()->GetTitle()+" (%)");
  h_vn_EpdE_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_EpdE->GetYaxis()->GetTitle());

  TH1D *h_vn_EpdF_flip = new TH1D("h_vn_EpdF_flip",h_vn_EpdF->GetTitle(),centBins,0,centBins);
  h_vn_EpdF_flip->GetXaxis()->SetTitle((TString)h_vn_EpdF->GetXaxis()->GetTitle()+" (%)");
  h_vn_EpdF_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_EpdF->GetYaxis()->GetTitle());

  TH1D *h_vn_TpcB_flip = new TH1D("h_vn_TpcB_flip",h_vn_TpcB->GetTitle(),centBins,0,centBins);
  h_vn_TpcB_flip->GetXaxis()->SetTitle((TString)h_vn_TpcB->GetXaxis()->GetTitle()+" (%)");
  h_vn_TpcB_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_TpcB->GetYaxis()->GetTitle());
  
  TH1D *h_vn_pp_flip = new TH1D("h_vn_pp_flip",h_vn_pp->GetTitle(),centBins,0,centBins);
  h_vn_pp_flip->GetXaxis()->SetTitle((TString)h_vn_pp->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pp->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_flip = new TH1D("h_vn_pm_flip",h_vn_pm->GetTitle(),centBins,0,centBins);
  h_vn_pm_flip->GetXaxis()->SetTitle((TString)h_vn_pm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pm->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_flip = new TH1D("h_vn_kp_flip",h_vn_kp->GetTitle(),centBins,0,centBins);
  h_vn_kp_flip->GetXaxis()->SetTitle((TString)h_vn_kp->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_kp->GetYaxis()->GetTitle());

  TH1D *h_vn_km_flip = new TH1D("h_vn_km_flip",h_vn_km->GetTitle(),centBins,0,centBins);
  h_vn_km_flip->GetXaxis()->SetTitle((TString)h_vn_km->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_km->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_flip = new TH1D("h_vn_pr_flip",h_vn_pr->GetTitle(),centBins,0,centBins);
  h_vn_pr_flip->GetXaxis()->SetTitle((TString)h_vn_pr->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pr->GetYaxis()->GetTitle());


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

  // Use these to change the x axis of centrality plots without fighting with TAxis.
  TH1D *vn_EpdE = new TH1D("vn_EpdE", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_EpdF = new TH1D("vn_EpdF", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_TpcB = new TH1D("vn_TpcB", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_pp = new TH1D("vn_pp", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_pm = new TH1D("vn_pm", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_kp = new TH1D("vn_kp", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_km = new TH1D("vn_km", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_pr = new TH1D("vn_pr", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);

  for (int i = 1; i <= 12; i++)
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
  

  THStack *piCentralityStack = new THStack("piCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kaCentralityStack = new THStack("kaCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");

  THStack *etaRegionStack    = new THStack("etaRegionStack", ";Centrality (%);v_{"+order_n_str+"}");
    
  THStack *ppRapidityStack   = new THStack("ppRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *pmRapidityStack   = new THStack("pmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *kpRapidityStack   = new THStack("kpRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *kmRapidityStack   = new THStack("kmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *prRapidityStack   = new THStack("prRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");

  vn_pp->SetMarkerStyle(20);
  vn_pp->SetMarkerSize(2);
  vn_pp->SetMarkerColor(2);
  vn_pp->SetLineColor(2);

  vn_pm->SetMarkerStyle(20);
  vn_pm->SetMarkerSize(2);
  vn_pm->SetMarkerColor(4);
  vn_pm->SetLineColor(4);

  vn_kp->SetMarkerStyle(20);
  vn_kp->SetMarkerSize(2);
  vn_kp->SetMarkerColor(2);
  vn_kp->SetLineColor(2);

  vn_km->SetMarkerStyle(20);
  vn_km->SetMarkerSize(2);
  vn_km->SetMarkerColor(4);
  vn_km->SetLineColor(4);

  vn_pr->SetMarkerStyle(20);
  vn_pr->SetMarkerSize(2);
  vn_pr->SetMarkerColor(2);
  vn_pr->SetLineColor(2);

  
  vn_EpdE->SetMarkerStyle(20);
  vn_EpdE->SetMarkerSize(2);
  vn_EpdE->SetMarkerColor(46);
  vn_EpdE->SetLineColor(46);

  vn_EpdF->SetMarkerStyle(20);
  vn_EpdF->SetMarkerSize(2);
  vn_EpdF->SetMarkerColor(38);
  vn_EpdF->SetLineColor(38);

  vn_TpcB->SetMarkerStyle(20);
  vn_TpcB->SetMarkerSize(2);
  vn_TpcB->SetMarkerColor(8);
  vn_TpcB->SetLineColor(8);


  h_vn_yCM_00to10_pp->SetMarkerStyle(20);
  h_vn_yCM_10to40_pp->SetMarkerStyle(20);
  h_vn_yCM_40to60_pp->SetMarkerStyle(20);
  h_vn_yCM_00to60_pp->SetMarkerStyle(20);
  h_vn_yCM_00to10_pp->SetMarkerColor(2);
  h_vn_yCM_10to40_pp->SetMarkerColor(4);
  h_vn_yCM_40to60_pp->SetMarkerColor(3);
  h_vn_yCM_00to60_pp->SetMarkerColor(4);
  h_vn_yCM_00to10_pp->SetMarkerSize(2);
  h_vn_yCM_10to40_pp->SetMarkerSize(2);
  h_vn_yCM_40to60_pp->SetMarkerSize(2);
  h_vn_yCM_00to60_pp->SetMarkerSize(2);
  h_vn_yCM_00to10_pp->SetLineColor(2);
  h_vn_yCM_10to40_pp->SetLineColor(4);
  h_vn_yCM_40to60_pp->SetLineColor(3);
  h_vn_yCM_00to60_pp->SetLineColor(4);

  h_vn_yCM_00to10_pm->SetMarkerStyle(20);
  h_vn_yCM_10to40_pm->SetMarkerStyle(20);
  h_vn_yCM_40to60_pm->SetMarkerStyle(20);
  h_vn_yCM_00to60_pm->SetMarkerStyle(20);
  h_vn_yCM_00to10_pm->SetMarkerColor(2);
  h_vn_yCM_10to40_pm->SetMarkerColor(4);
  h_vn_yCM_40to60_pm->SetMarkerColor(3);
  h_vn_yCM_00to60_pm->SetMarkerColor(4);
  h_vn_yCM_00to10_pm->SetMarkerSize(2);
  h_vn_yCM_10to40_pm->SetMarkerSize(2);
  h_vn_yCM_40to60_pm->SetMarkerSize(2);
  h_vn_yCM_00to60_pm->SetMarkerSize(2);
  h_vn_yCM_00to10_pm->SetLineColor(2);
  h_vn_yCM_10to40_pm->SetLineColor(4);
  h_vn_yCM_40to60_pm->SetLineColor(3);
  h_vn_yCM_00to60_pm->SetLineColor(4);

  h_vn_yCM_00to10_kp->SetMarkerStyle(20);
  h_vn_yCM_10to40_kp->SetMarkerStyle(20);
  h_vn_yCM_40to60_kp->SetMarkerStyle(20);
  h_vn_yCM_00to60_kp->SetMarkerStyle(20);
  h_vn_yCM_00to10_kp->SetMarkerColor(2);
  h_vn_yCM_10to40_kp->SetMarkerColor(4);
  h_vn_yCM_40to60_kp->SetMarkerColor(3);
  h_vn_yCM_00to60_kp->SetMarkerColor(4);
  h_vn_yCM_00to10_kp->SetMarkerSize(2);
  h_vn_yCM_10to40_kp->SetMarkerSize(2);
  h_vn_yCM_40to60_kp->SetMarkerSize(2);
  h_vn_yCM_00to60_kp->SetMarkerSize(2);
  h_vn_yCM_00to10_kp->SetLineColor(2);
  h_vn_yCM_10to40_kp->SetLineColor(4);
  h_vn_yCM_40to60_kp->SetLineColor(3);
  h_vn_yCM_00to60_kp->SetLineColor(4);

  h_vn_yCM_00to10_km->SetMarkerStyle(20);
  h_vn_yCM_10to40_km->SetMarkerStyle(20);
  h_vn_yCM_40to60_km->SetMarkerStyle(20);
  h_vn_yCM_00to60_km->SetMarkerStyle(20);
  h_vn_yCM_00to10_km->SetMarkerColor(2);
  h_vn_yCM_10to40_km->SetMarkerColor(4);
  h_vn_yCM_40to60_km->SetMarkerColor(3);
  h_vn_yCM_00to60_km->SetMarkerColor(4);
  h_vn_yCM_00to10_km->SetMarkerSize(2);
  h_vn_yCM_10to40_km->SetMarkerSize(2);
  h_vn_yCM_40to60_km->SetMarkerSize(2);
  h_vn_yCM_00to60_km->SetMarkerSize(2);
  h_vn_yCM_00to10_km->SetLineColor(2);
  h_vn_yCM_10to40_km->SetLineColor(4);
  h_vn_yCM_40to60_km->SetLineColor(3);
  h_vn_yCM_00to60_km->SetLineColor(4);

  h_vn_yCM_00to10_pr->SetMarkerStyle(20);
  h_vn_yCM_10to40_pr->SetMarkerStyle(20);
  h_vn_yCM_40to60_pr->SetMarkerStyle(20);
  h_vn_yCM_00to60_pr->SetMarkerStyle(20);
  h_vn_yCM_00to10_pr->SetMarkerColor(2);
  h_vn_yCM_10to40_pr->SetMarkerColor(4);
  h_vn_yCM_40to60_pr->SetMarkerColor(3);
  h_vn_yCM_00to60_pr->SetMarkerColor(4);
  h_vn_yCM_00to10_pr->SetMarkerSize(2);
  h_vn_yCM_10to40_pr->SetMarkerSize(2);
  h_vn_yCM_40to60_pr->SetMarkerSize(2);
  h_vn_yCM_00to60_pr->SetMarkerSize(2);
  h_vn_yCM_00to10_pr->SetLineColor(2);
  h_vn_yCM_10to40_pr->SetLineColor(4);
  h_vn_yCM_40to60_pr->SetLineColor(3);
  h_vn_yCM_00to60_pr->SetLineColor(4);

  
  piCentralityStack->Add(vn_pp);
  piCentralityStack->Add(vn_pm);

  kaCentralityStack->Add(vn_kp);
  kaCentralityStack->Add(vn_km);

  etaRegionStack->Add(vn_EpdE);
  etaRegionStack->Add(vn_EpdF);
  etaRegionStack->Add(vn_TpcB);

  ppRapidityStack->Add(h_vn_yCM_00to10_pp);
  ppRapidityStack->Add(h_vn_yCM_10to40_pp);
  ppRapidityStack->Add(h_vn_yCM_40to60_pp);

  pmRapidityStack->Add(h_vn_yCM_00to10_pm);
  pmRapidityStack->Add(h_vn_yCM_10to40_pm);
  pmRapidityStack->Add(h_vn_yCM_40to60_pm);

  kpRapidityStack->Add(h_vn_yCM_00to10_kp);
  kpRapidityStack->Add(h_vn_yCM_10to40_kp);
  kpRapidityStack->Add(h_vn_yCM_40to60_kp);

  kmRapidityStack->Add(h_vn_yCM_00to10_km);
  kmRapidityStack->Add(h_vn_yCM_10to40_km);
  kmRapidityStack->Add(h_vn_yCM_40to60_km);

  prRapidityStack->Add(h_vn_yCM_00to10_pr);
  prRapidityStack->Add(h_vn_yCM_10to40_pr);
  prRapidityStack->Add(h_vn_yCM_40to60_pr);

  if (order_n_str == "2")
    {
      TLegend *piLegend = new TLegend(0.775, 0.625, 0.9, 0.775);
      piLegend->AddEntry(vn_pp,"#pi^{+}");
      piLegend->AddEntry(vn_pm,"#pi^{-}");
      piLegend->SetFillColorAlpha(0,0);
      piLegend->SetLineColorAlpha(0,0);

      TLegend *kaLegend = new TLegend(0.275, 0.32, 0.425, 0.45);
      kaLegend->AddEntry(vn_kp,"K^{+}");
      kaLegend->AddEntry(vn_km,"K^{-}");
      kaLegend->SetFillColorAlpha(0,0);
      kaLegend->SetLineColorAlpha(0,0);

      TLegend *etaLegend = new TLegend(0.65, 0.7, 0.9, 0.9);
      etaLegend->AddEntry(vn_EpdE, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(vn_EpdF, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(vn_TpcB, "TPC -1.0 < #eta < 0");

      TLegend *ppLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      ppLegend->AddEntry(h_vn_yCM_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(h_vn_yCM_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(h_vn_yCM_40to60_pp, "40 - 60%");

      TLegend *pmLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      pmLegend->AddEntry(h_vn_yCM_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(h_vn_yCM_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(h_vn_yCM_40to60_pm, "40 - 60%");

      TLegend *kpLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      kpLegend->AddEntry(h_vn_yCM_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(h_vn_yCM_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(h_vn_yCM_40to60_kp, "40 - 60%");

      TLegend *kmLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      kmLegend->AddEntry(h_vn_yCM_00to10_km, "0 - 10%");
      kmLegend->AddEntry(h_vn_yCM_10to40_km, "10 - 40%");
      kmLegend->AddEntry(h_vn_yCM_40to60_km, "40 - 60%");

      TLegend *prLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      prLegend->AddEntry(h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_yCM_40to60_pr, "40 - 60%");

      canvas->SetLeftMargin(0.13);
      canvas->SetGridx();
      canvas->SetGridy();
      gStyle->SetErrorX(0);
      gStyle->SetOptStat(0);

      TLine *zeroLine = new TLine(0, 0, 60, 0);
      zeroLine->SetLineStyle(9);

      TLine *zeroLine_y = new TLine(-1, 0, 1, 0);
      zeroLine_y->SetLineStyle(9);

  
      piCentralityStack->Draw();
      piCentralityStack->GetXaxis()->SetNdivisions(210);
      piCentralityStack->SetMaximum(0.01);
      piCentralityStack->SetMinimum(-0.05);
      piCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      piLegend->Draw();
      canvas->SaveAs(jobID + "_piCentralityStack.png");
      canvas->Clear();

      kaCentralityStack->Draw();
      kaCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kaCentralityStack->GetXaxis()->SetNdivisions(210);
      kaCentralityStack->SetMaximum(0.02);
      kaCentralityStack->SetMinimum(-0.12);
      kaCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kaLegend->Draw();
      canvas->SaveAs(jobID + "_kaCentralityStack.png");
      canvas->Clear();

      vn_pr->SetTitle("");
      vn_pr->GetYaxis()->SetTitleOffset(1.7);
      vn_pr->GetXaxis()->SetNdivisions(210);
      vn_pr->Draw("E1P");
      vn_pr->SetMaximum(0.01);
      vn_pr->SetMinimum(-0.09);
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_pr.png");
      canvas->Clear();

      vn_EpdE->SetMarkerStyle(20);
      vn_EpdE->SetMarkerSize(2);
      //vn_EpdE->SetMarkerColor(2);
      //vn_EpdE->SetLineColor(2);
      vn_EpdE->GetXaxis()->SetNdivisions(210);
      vn_EpdE->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdE.png");
      canvas->Clear();

      vn_EpdF->SetMarkerStyle(20);
      vn_EpdF->SetMarkerSize(2);
      //vn_EpdF->SetMarkerColor(2);
      //vn_EpdF->SetLineColor(2);
      vn_EpdF->GetXaxis()->SetNdivisions(210);
      vn_EpdF->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdF.png");
      canvas->Clear();

      vn_TpcB->SetMarkerStyle(20);
      vn_TpcB->SetMarkerSize(2);
      //vn_TpcB->SetMarkerColor(2);
      //vn_TpcB->SetLineColor(2);
      vn_TpcB->GetXaxis()->SetNdivisions(210);
      vn_TpcB->GetYaxis()->SetTitleOffset(1.7);
      vn_TpcB->Draw("E1P");
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_TpcB.png");
      canvas->Clear();

      
      etaRegionStack->Draw();
      etaRegionStack->GetYaxis()->SetTitleOffset(1.7);
      etaRegionStack->GetXaxis()->SetNdivisions(210);
      etaRegionStack->Draw();
      //etaRegionStack->SetMaximum(0.3);
      etaRegionStack->SetMinimum(-0.1);
      etaRegionStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      etaLegend->Draw();
      canvas->SaveAs(jobID + "_etaRegionStack.png");
      canvas->Clear();


      ppRapidityStack->Draw();
      ppRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      ppRapidityStack->GetXaxis()->SetNdivisions(210);
      ppRapidityStack->Draw();
      ppRapidityStack->SetMaximum(0.03);
      ppRapidityStack->SetMinimum(-0.05);
      ppRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      ppLegend->Draw();
      canvas->SaveAs(jobID + "_ppRapidityStack.png");
      canvas->Clear();

      pmRapidityStack->Draw();
      pmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      pmRapidityStack->GetXaxis()->SetNdivisions(210);
      pmRapidityStack->Draw();
      pmRapidityStack->SetMaximum(0.03);
      pmRapidityStack->SetMinimum(-0.05);
      pmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      pmLegend->Draw();
      canvas->SaveAs(jobID + "_pmRapidityStack.png");
      canvas->Clear();

      kpRapidityStack->Draw();
      kpRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kpRapidityStack->GetXaxis()->SetNdivisions(210);
      kpRapidityStack->Draw();
      kpRapidityStack->SetMaximum(0.03);
      kpRapidityStack->SetMinimum(-0.05);
      kpRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kpLegend->Draw();
      canvas->SaveAs(jobID + "_kpRapidityStack.png");
      canvas->Clear();

      kmRapidityStack->Draw();
      kmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kmRapidityStack->GetXaxis()->SetNdivisions(210);
      kmRapidityStack->Draw();
      kmRapidityStack->SetMaximum(0.03);
      kmRapidityStack->SetMinimum(-0.05);
      kmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kmLegend->Draw();
      canvas->SaveAs(jobID + "_kmRapidityStack.png");
      canvas->Clear();

      prRapidityStack->Draw();
      prRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      prRapidityStack->GetXaxis()->SetNdivisions(210);
      prRapidityStack->Draw();
      prRapidityStack->SetMaximum(0.1);
      prRapidityStack->SetMinimum(-0.1);
      prRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      prLegend->Draw();
      canvas->SaveAs(jobID + "_prRapidityStack.png");
      canvas->Clear();
    }
  else if (order_n_str == "3")
    {
      TLegend *piLegend = new TLegend(0.675, 0.625, 0.8, 0.775);
      piLegend->AddEntry(vn_pp,"#pi^{+}");
      piLegend->AddEntry(vn_pm,"#pi^{-}");
      piLegend->SetFillColorAlpha(0,0);
      piLegend->SetLineColorAlpha(0,0);

      TLegend *kaLegend = new TLegend(0.275, 0.15, 0.425, 0.27);
      kaLegend->AddEntry(vn_kp,"K^{+}");
      kaLegend->AddEntry(vn_km,"K^{-}");
      kaLegend->SetFillColorAlpha(0,0);
      kaLegend->SetLineColorAlpha(0,0);


      TLegend *etaLegend = new TLegend(0.65, 0.25, 0.9, 0.45);
      etaLegend->AddEntry(vn_EpdE, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(vn_EpdF, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(vn_TpcB, "TPC -1.0 < #eta < 0");


      TLegend *ppLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      ppLegend->AddEntry(h_vn_yCM_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(h_vn_yCM_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(h_vn_yCM_40to60_pp, "40 - 60%");

      TLegend *pmLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      pmLegend->AddEntry(h_vn_yCM_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(h_vn_yCM_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(h_vn_yCM_40to60_pm, "40 - 60%");

      TLegend *kpLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      kpLegend->AddEntry(h_vn_yCM_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(h_vn_yCM_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(h_vn_yCM_40to60_kp, "40 - 60%");

      TLegend *kmLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      kmLegend->AddEntry(h_vn_yCM_00to10_km, "0 - 10%");
      kmLegend->AddEntry(h_vn_yCM_10to40_km, "10 - 40%");
      kmLegend->AddEntry(h_vn_yCM_40to60_km, "40 - 60%");

      TLegend *prLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      prLegend->AddEntry(h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_yCM_40to60_pr, "40 - 60%");

      canvas->SetLeftMargin(0.13);
      canvas->SetGridx();
      canvas->SetGridy();
      gStyle->SetErrorX(0);
      gStyle->SetOptStat(0);

      TLine *zeroLine = new TLine(0, 0, 60, 0);
      zeroLine->SetLineStyle(9);

      TLine *zeroLine_y = new TLine(-1, 0, 1, 0);
      zeroLine_y->SetLineStyle(9);

  
      piCentralityStack->Draw();
      piCentralityStack->GetYaxis()->SetTitleOffset(1.8);
      piCentralityStack->GetXaxis()->SetNdivisions(210);
      piCentralityStack->SetMaximum(0.015);
      piCentralityStack->SetMinimum(-0.015);
      piCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      piLegend->Draw();
      canvas->SaveAs(jobID + "_piCentralityStack.png");
      canvas->Clear();

      kaCentralityStack->Draw();
      kaCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kaCentralityStack->GetXaxis()->SetNdivisions(210);
      kaCentralityStack->SetMaximum(0.1);
      kaCentralityStack->SetMinimum(-0.16);
      kaCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kaLegend->Draw();
      canvas->SaveAs(jobID + "_kaCentralityStack.png");
      canvas->Clear();

      vn_pr->SetTitle("");
      vn_pr->GetYaxis()->SetTitleOffset(1.7);
      vn_pr->GetXaxis()->SetNdivisions(210);
      vn_pr->Draw("E1P");
      vn_pr->SetMaximum(0.01);
      vn_pr->SetMinimum(-0.03);
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_pr.png");
      canvas->Clear();

      vn_EpdE->SetMarkerStyle(20);
      vn_EpdE->SetMarkerSize(2);
      //vn_EpdE->SetMarkerColor(2);
      //vn_EpdE->SetLineColor(2);
      vn_EpdE->GetYaxis()->SetTitleOffset(1.7);
      vn_EpdE->GetXaxis()->SetNdivisions(210);
      vn_EpdE->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdE.png");
      canvas->Clear();

      vn_EpdF->SetMarkerStyle(20);
      vn_EpdF->SetMarkerSize(2);
      //vn_EpdF->SetMarkerColor(2);
      //vn_EpdF->SetLineColor(2);
      vn_EpdF->GetYaxis()->SetTitleOffset(1.8);
      vn_EpdF->GetXaxis()->SetNdivisions(210);
      vn_EpdF->Draw("E1P");
      canvas->SaveAs(jobID + "_vn_EpdF.png");
      canvas->Clear();

      vn_TpcB->SetMarkerStyle(20);
      vn_TpcB->SetMarkerSize(2);
      //vn_TpcB->SetMarkerColor(2);
      //vn_TpcB->SetLineColor(2);
      vn_TpcB->GetXaxis()->SetNdivisions(210);
      vn_TpcB->GetYaxis()->SetTitleOffset(1.7);
      vn_TpcB->Draw("E1P");
      vn_TpcB->SetMaximum(0.01);
      vn_TpcB->SetMinimum(-0.03);
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_TpcB.png");
      canvas->Clear();


      etaRegionStack->Draw();
      etaRegionStack->GetYaxis()->SetTitleOffset(1.7);
      etaRegionStack->GetXaxis()->SetNdivisions(210);
      etaRegionStack->Draw();
      etaRegionStack->SetMaximum(0.05);
      //etaRegionStack->SetMinimum(-0.1);
      etaRegionStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      etaLegend->Draw();
      canvas->SaveAs(jobID + "_etaRegionStack.png");
      canvas->Clear();


      ppRapidityStack->Draw();
      ppRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      ppRapidityStack->GetXaxis()->SetNdivisions(210);
      ppRapidityStack->Draw();
      ppRapidityStack->SetMaximum(0.02);
      ppRapidityStack->SetMinimum(-0.01);
      ppRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      ppLegend->Draw();
      canvas->SaveAs(jobID + "_ppRapidityStack.png");
      canvas->Clear();

      pmRapidityStack->Draw();
      pmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      pmRapidityStack->GetXaxis()->SetNdivisions(210);
      pmRapidityStack->Draw();
      pmRapidityStack->SetMaximum(0.01);
      pmRapidityStack->SetMinimum(-0.02);
      pmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      pmLegend->Draw();
      canvas->SaveAs(jobID + "_pmRapidityStack.png");
      canvas->Clear();

      kpRapidityStack->Draw();
      kpRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kpRapidityStack->GetXaxis()->SetNdivisions(210);
      kpRapidityStack->Draw();
      kpRapidityStack->SetMaximum(0.03);
      kpRapidityStack->SetMinimum(-0.05);
      kpRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kpLegend->Draw();
      canvas->SaveAs(jobID + "_kpRapidityStack.png");
      canvas->Clear();

      kmRapidityStack->Draw();
      kmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kmRapidityStack->GetXaxis()->SetNdivisions(210);
      kmRapidityStack->Draw();
      kmRapidityStack->SetMaximum(0.03);
      kmRapidityStack->SetMinimum(-0.05);
      kmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kmLegend->Draw();
      canvas->SaveAs(jobID + "_kmRapidityStack.png");
      canvas->Clear();

      prRapidityStack->Draw();
      prRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      prRapidityStack->GetXaxis()->SetNdivisions(210);
      prRapidityStack->Draw();
      prRapidityStack->SetMaximum(0.02);
      prRapidityStack->SetMinimum(-0.04);
      prRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      prLegend->Draw();
      canvas->SaveAs(jobID + "_prRapidityStack.png");
      canvas->Clear();
    }

  //resolutionInfo_INPUT->Close();
  file->Close();
}
