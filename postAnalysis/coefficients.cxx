void applyResolution(TH1D *histogram, Double_t resolution, Double_t resolutionError)
{
  for (int i = 1; i < histogram->GetNbinsX(); i++)
    {
      Double_t rawBinContent = histogram->GetBinContent(i);
      if (rawBinContent == 0.0) continue;
      
      Double_t rawBinError = histogram->GetBinError(i);

      Double_t newBinContent = rawBinContent/resolution;
      Double_t newBinError = newBinContent * TMath::Sqrt( pow(rawBinError/rawBinContent, 2) + pow(resolutionError/resolution, 2) );

      histogram->SetBinContent(i, newBinContent);
      histogram->SetBinError(i, newBinError);
    }
}



void coefficients(TString jobID, TString order_n_str)
{
  //TH1::SetDefaultSumw2();
  
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}
  /*
  TFile *resolutionInfo_INPUT = TFile::Open("resolutionInfo_INPUT.root", "READ");
  if(!resolutionInfo_INPUT) { cout << "No resolution file found!" << endl; return; }

  TH1D *h_resolutions = (TH1D*)resolutionInfo_INPUT->Get("h_resolutions");
  TH2D *h2_resolutions = (TH2D*)resolutionInfo_INPUT->Get("h2_resolutions");
  */
  /*
  Double_t resolutionIDs[16] = {};
  Double_t resolutionIDsError[16] = {};
  
  for (int i = 1; i < h_resolutions->GetNbinsX(); i++)
    {
      resolutionIDs[i-1] = h_resolutions->GetBinContent(i);
      resolutionIDsError[i-1] = h_resolutions->GetBinError(i);
    }
  */
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 800);
  canvas->SetGridx();
  canvas->SetGridy();
  //canvas->SetLeftMargin(0.15);
  canvas->cd();

  // Shaowei's data:
  double cent_label[7]={2.5, 7.5, 15, 25, 35, 45, 55};
  double cent_bin_lows[8]={0, 5, 10, 20, 30, 40, 50, 60};

  double v2val_cent_prp_data[7]={-0.00662406, -0.00906423, -0.0129282, -0.0202661, -0.0315561, -0.0467831, -0.0654673};
  double v2err_cent_prp_data[7]={0.000562772, 0.000228575, 9.699e-05, 8.17342e-05, 9.80836e-05, 0.000138984, 0.000257859};

  double v2val_cent_pip_data[7]={-0.0121035, -0.0200037, -0.0244652, -0.0279546, -0.0312375, -0.035842, -0.0413696};
  double v2err_cent_pip_data[7]={0.00136922, 0.000563976, 0.000240995, 0.000202438, 0.000239846, 0.000333549, 0.000602461};

  double v2val_cent_pim_data[7]={-0.00579955, -0.0125058, -0.0149095, -0.0170165, -0.0194754, -0.0216241, -0.0261989};
  double v2err_cent_pim_data[7]={0.00122163, 0.000505636, 0.000215966, 0.000181535, 0.000215232, 0.000298714, 0.000536848};

  double v2val_cent_kap_data[7]={-0.00341182, -0.0202501, -0.0132731, -0.0183739, -0.0281251, -0.0409518, -0.0722278};
  double v2err_cent_kap_data[7]={0.00604229, 0.00254523, 0.00113305, 0.00101225, 0.0012811, 0.00187359, 0.00345125};

  double v2val_cent_kam_data[7]={-0.0593462, -0.0051962, -0.0205428, -0.0224589, -0.0185285, -0.00424369, -0.0782483};
  double v2err_cent_kam_data[7]={0.0239337, 0.0106676, 0.00477639, 0.00426836, 0.00530658, 0.00769645, 0.0140865};

  /*
  //{-0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05};
  double y_prp_data[10]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
  double v2val_y_prp_data[10]={0.0374998, 0.0196273, 0.00722679, -0.00228046, -0.00941107, -0.0142688, -0.01755, -0.0192977, -0.0208411, -0.0221447};
  double v2err_y_prp_data[10]={0.00018097, 0.000154358, 0.000138407, 0.000133027, 0.000130472, 0.000130127, 0.00012812, 0.000127518, 0.000125894, 0.000124516};

  double y_pip_data[10]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
  double v2val_y_pip_data[10]={-0.0311928, -0.031139, -0.0300568, -0.028663, -0.0281598, -0.0276927, -0.0272488, -0.0267038, -0.0262154, -0.0270011};
  double v2err_y_pip_data[10]={0.000474489, 0.000431578, 0.000398728, 0.000371556, 0.000350646, 0.000330814, 0.000319569, 0.000310681, 0.000307709, 0.000304301};

  double y_pim_data[10]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
  double v2val_y_pim_data[10]={-0.0204974, -0.0213052, -0.0186106, -0.0185628, -0.017647, -0.0168865, -0.0170284, -0.0160957, -0.0160634, -0.0165181};
  double v2err_y_pim_data[10]={0.000436394, 0.000392559, 0.000360703, 0.000334713, 0.000314866, 0.000296397, 0.000286299, 0.000278582, 0.000276067, 0.000273199};

  //{-0.9, -0.7, -0.5, -0.3, -0.1};
  double y_kap_data[5]={0.15, 0.35, 0.55, 0.75, 0.95};
  double v2val_y_kap_data[5]={-0.0364406, -0.0285906, -0.0251823, -0.0176064, -0.0169012};
  double v2err_y_kap_data[5]={0.00296857, 0.00201864, 0.00182397, 0.00162463, 0.00150312};

  double y_kam_data[5]={0.15, 0.35, 0.55, 0.75, 0.95};
  double v2val_y_kam_data[5]={-0.0250802, -0.0234748, -0.0324117, -0.0128608, -0.0212328};
  double v2err_y_kam_data[5]={0.0114492, 0.00788562, 0.00580941, 0.00480216, 0.00424148};
  */

  //double y_prp_data[10]={-0.05, -0.15, -0.25, -0.35, -0.45, -0.55, -0.65, -0.75, -0.85, -0.95};
  double y_prp_data[10]={0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
  double v2val_y_prp_data[10]={-0.0221447, -0.0208411, -0.0192977, -0.01755, -0.0142688, -0.00941107, -0.00228046, 0.00722679, 0.0196273, 0.0374998};
  double v2err_y_prp_data[10]={0.000124516, 0.000125894, 0.000127518, 0.00012812, 0.000130127, 0.000130472, 0.000133027, 0.000138407, 0.000154358, 0.00018097};

  double y_pip_data[10]={0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
  double v2val_y_pip_data[10]={ -0.0270011, -0.0262154, -0.0267038, -0.0272488, -0.0276927, -0.0281598, -0.028663, -0.0300568, -0.031139, -0.0311928};
  double v2err_y_pip_data[10]={0.000304301, 0.000307709, 0.000310681, 0.000319569, 0.000330814, 0.000350646, 0.000371556, 0.000398728, 0.000431578, 0.000474489};

  double y_pim_data[10]={0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
  double v2val_y_pim_data[10]={-0.0165181, -0.0160634, -0.0160957, -0.0170284, -0.0168865, -0.017647, -0.0185628, -0.0186106, -0.0213052, -0.0204974,};
  double v2err_y_pim_data[10]={ 0.000273199, 0.000276067, 0.000278582, 0.000286299, 0.000296397, 0.000314866, 0.000334713, 0.000360703, 0.000392559, 0.000436394};

  //double y_kap_data[5]={0.15, 0.35, 0.55, 0.75, 0.95};
  double y_kap_data[5]={0.1,0.3, 0.5, 0.7, 0.9};
  double v2val_y_kap_data[5]={-0.0169012, -0.0176064, -0.0251823, -0.0285906, -0.0364406};
  double v2err_y_kap_data[5]={0.00150312, 0.00162463, 0.00182397, 0.00201864, 0.00296857};

  //double y_kam_data[5]={0.15, 0.35, 0.55, 0.75, 0.95};
  double y_kam_data[5]={0.1,0.3, 0.5, 0.7, 0.9};
  double v2val_y_kam_data[5]={-0.0212328, -0.0128608, -0.0324117, -0.0234748, -0.0250802};
  double v2err_y_kam_data[5]={0.00424148, 0.00480216, 0.00580941, 0.00788562, 0.0114492};
  ////////

  
  TH1D *sh_cent_pp = new TH1D("sh_cent_pp", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_pp->FillN(7, cent_label, v2val_cent_pip_data);
  for (int i = 1; i <= sh_cent_pp->GetNbinsX(); i++) { sh_cent_pp->SetBinError(i, v2err_cent_pip_data[i-1]); }
  
  TH1D *sh_cent_pm = new TH1D("sh_cent_pm", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_pm->FillN(7, cent_label, v2val_cent_pim_data);
  for (int i = 1; i <= sh_cent_pm->GetNbinsX(); i++) { sh_cent_pm->SetBinError(i, v2err_cent_pim_data[i-1]); }
  
  TH1D *sh_cent_kp = new TH1D("sh_cent_kp", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_kp->FillN(7, cent_label, v2val_cent_kap_data);
  for (int i = 1; i <= sh_cent_kp->GetNbinsX(); i++) { sh_cent_kp->SetBinError(i, v2err_cent_kap_data[i-1]); }
  
  TH1D *sh_cent_km = new TH1D("sh_cent_km", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_km->FillN(7, cent_label, v2val_cent_kam_data);
  for (int i = 1; i <= sh_cent_km->GetNbinsX(); i++) { sh_cent_km->SetBinError(i, v2err_cent_kam_data[i-1]); }
  
  TH1D *sh_cent_pr = new TH1D("sh_cent_pr", ";Centrality (%); v_{"+order_n_str+"}", 7, cent_bin_lows);
  sh_cent_pr->FillN(7, cent_label, v2val_cent_prp_data);
  for (int i = 1; i <= sh_cent_pr->GetNbinsX(); i++) { sh_cent_pr->SetBinError(i, v2err_cent_prp_data[i-1]); }


  TH1D *sh_y_pp = new TH1D("sh_y_pp", ";y-y_{mid};v_{2}", 20, -1, 1);
  sh_y_pp->FillN(10, y_pip_data, v2val_y_pip_data);
  for (int i = 0; i < 10; i++) { sh_y_pp->SetBinError(i+11, v2err_y_pip_data[i]); }
  
  TH1D *sh_y_pm = new TH1D("sh_y_pm", ";y-y_{mid};v_{2}", 20, -1, 1);
  sh_y_pm->FillN(10, y_pim_data, v2val_y_pim_data);
  for (int i = 0; i < 10; i++) { sh_y_pm->SetBinError(i+11, v2err_y_pim_data[i]); }

  TH1D *sh_y_pr = new TH1D("sh_y_pr", ";y-y_{mid};v_{2}", 20, -1, 1);
  sh_y_pr->FillN(10, y_prp_data, v2val_y_prp_data);
  for (int i = 0; i < 10; i++) { sh_y_pr->SetBinError(i+11, v2err_y_prp_data[i]); }
  
  TH1D *sh_y_kp = new TH1D("sh_y_kp", ";y-y_{mid};v_{2}", 10, -1, 1);
  sh_y_kp->FillN(5, y_kap_data, v2val_y_kap_data);
  for (int i = 0; i < 5; i++) { sh_y_kp->SetBinError(i+6, v2err_y_kap_data[i]); }
  
  TH1D *sh_y_km = new TH1D("sh_y_km", ";y-y_{mid};v_{2}", 10, -1, 1);
  sh_y_km->FillN(5, y_kam_data, v2val_y_kam_data);
  for (int i = 0; i < 5; i++) { sh_y_km->SetBinError(i+6, v2err_y_kam_data[i]); }
  ////

  TProfile *p_vn_EpdE = (TProfile*)file->Get("p_vn_EpdE");
  TProfile *p_vn_EpdF = (TProfile*)file->Get("p_vn_EpdF");
  TProfile *p_vn_TpcB = (TProfile*)file->Get("p_vn_TpcB");
  TProfile *p_vn_Tpc = (TProfile*)file->Get("p_vn_Tpc_pT_0p2to2");
  TProfile *p_vn_pp = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr = (TProfile*)file->Get("p_vn_pr");

  p_vn_kp->Rebin();
  p_vn_km->Rebin();

  TProfile *p_vn_pp_ext = (TProfile*)file->Get("p_vn_pp_ext");
  TProfile *p_vn_pm_ext = (TProfile*)file->Get("p_vn_pm_ext");
  TProfile *p_vn_kp_ext = (TProfile*)file->Get("p_vn_kp_ext");
  TProfile *p_vn_km_ext = (TProfile*)file->Get("p_vn_km_ext");
  TProfile *p_vn_pr_ext = (TProfile*)file->Get("p_vn_pr_ext");

  p_vn_kp_ext->Rebin();
  p_vn_km_ext->Rebin();

  TProfile2D *p2_vn_yCM_cent_pp = (TProfile2D*)file->Get("p2_vn_yCM_cent_pp");
  TProfile2D *p2_vn_yCM_cent_pm = (TProfile2D*)file->Get("p2_vn_yCM_cent_pm");
  TProfile2D *p2_vn_yCM_cent_kp = (TProfile2D*)file->Get("p2_vn_yCM_cent_kp");
  TProfile2D *p2_vn_yCM_cent_km = (TProfile2D*)file->Get("p2_vn_yCM_cent_km");
  TProfile2D *p2_vn_yCM_cent_pr = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr");

  p2_vn_yCM_cent_kp->RebinY();
  p2_vn_yCM_cent_km->RebinY();
  
  /*
  // RESOLUTION APPLICATION
  TH1D *h_vn_EpdE = p_vn_EpdE->ProjectionX();
  TH1D *h_vn_EpdF = p_vn_EpdF->ProjectionX();
  TH1D *h_vn_TpcB = p_vn_TpcB->ProjectionX();
  TH1D *h_vn_pp = p_vn_pp->ProjectionX();
  TH1D *h_vn_pm = p_vn_pm->ProjectionX();
  TH1D *h_vn_kp = p_vn_kp->ProjectionX();
  TH1D *h_vn_km = p_vn_km->ProjectionX();
  TH1D *h_vn_pr = p_vn_pr->ProjectionX();

  TH1D *h_vn_pp_ext = p_vn_pp_ext->ProjectionX();
  TH1D *h_vn_pm_ext = p_vn_pm_ext->ProjectionX();
  TH1D *h_vn_kp_ext = p_vn_kp_ext->ProjectionX();
  TH1D *h_vn_km_ext = p_vn_km_ext->ProjectionX();
  TH1D *h_vn_pr_ext = p_vn_pr_ext->ProjectionX();

  TH2D *h2_vn_yCM_cent_pp = p2_vn_yCM_cent_pp->ProjectionXY();
  TH2D *h2_vn_yCM_cent_pm = p2_vn_yCM_cent_pm->ProjectionXY();
  TH2D *h2_vn_yCM_cent_kp = p2_vn_yCM_cent_kp->ProjectionXY();
  TH2D *h2_vn_yCM_cent_km = p2_vn_yCM_cent_km->ProjectionXY();
  TH2D *h2_vn_yCM_cent_pr = p2_vn_yCM_cent_pr->ProjectionXY();

  h_vn_EpdE->Divide(h_resolutions);
  h_vn_EpdF->Divide(h_resolutions);
  h_vn_TpcB->Divide(h_resolutions);
  h_vn_pp->Divide(h_resolutions);
  h_vn_pm->Divide(h_resolutions);
  h_vn_kp->Divide(h_resolutions);
  h_vn_km->Divide(h_resolutions);
  h_vn_pr->Divide(h_resolutions);

  h_vn_pp_ext->Divide(h_resolutions);
  h_vn_pm_ext->Divide(h_resolutions);
  h_vn_kp_ext->Divide(h_resolutions);
  h_vn_km_ext->Divide(h_resolutions);
  h_vn_pr_ext->Divide(h_resolutions);
  
  h2_vn_yCM_cent_pp->Divide(h2_resolutions);
  h2_vn_yCM_cent_pm->Divide(h2_resolutions);
  h2_vn_yCM_cent_kp->Divide(h2_resolutions);
  h2_vn_yCM_cent_km->Divide(h2_resolutions);
  h2_vn_yCM_cent_pr->Divide(h2_resolutions);
  */
  ////

  
  /*
  TH1D *pp_cent15 = p2_vn_yCM_cent_pp->ProfileX("pp_cent15", 16, 16)->ProjectionX();
  TH1D *pp_cent14 = p2_vn_yCM_cent_pp->ProfileX("pp_cent14", 15, 15)->ProjectionX();
  TH1D *pp_cent13 = p2_vn_yCM_cent_pp->ProfileX("pp_cent13", 14, 14)->ProjectionX();
  TH1D *pp_cent12 = p2_vn_yCM_cent_pp->ProfileX("pp_cent12", 13, 13)->ProjectionX();
  TH1D *pp_cent11 = p2_vn_yCM_cent_pp->ProfileX("pp_cent11", 12, 12)->ProjectionX();
  TH1D *pp_cent10 = p2_vn_yCM_cent_pp->ProfileX("pp_cent10", 11, 11)->ProjectionX();
  TH1D *pp_cent9 = p2_vn_yCM_cent_pp->ProfileX("pp_cent9", 10, 10)->ProjectionX();
  TH1D *pp_cent8 = p2_vn_yCM_cent_pp->ProfileX("pp_cent8", 9, 9)->ProjectionX();
  TH1D *pp_cent7 = p2_vn_yCM_cent_pp->ProfileX("pp_cent7", 8, 8)->ProjectionX();
  TH1D *pp_cent6 = p2_vn_yCM_cent_pp->ProfileX("pp_cent6", 7, 7)->ProjectionX();
  TH1D *pp_cent5 = p2_vn_yCM_cent_pp->ProfileX("pp_cent5", 6, 6)->ProjectionX();
  TH1D *pp_cent4 = p2_vn_yCM_cent_pp->ProfileX("pp_cent3", 5, 5)->ProjectionX();

  std::cout << "pp_cent10 bin 18 content: " << pp_cent10->GetBinContent(18) << std::endl;

  applyResolution(pp_cent15, resolutionIDs[15], resolutionIDsError[15]);
  applyResolution(pp_cent14, resolutionIDs[14], resolutionIDsError[14]);
  applyResolution(pp_cent13, resolutionIDs[13], resolutionIDsError[13]);
  applyResolution(pp_cent12, resolutionIDs[12], resolutionIDsError[12]);
  applyResolution(pp_cent11, resolutionIDs[11], resolutionIDsError[11]);
  applyResolution(pp_cent10, resolutionIDs[10], resolutionIDsError[10]);
  applyResolution(pp_cent9, resolutionIDs[9], resolutionIDsError[9]);
  applyResolution(pp_cent8, resolutionIDs[8], resolutionIDsError[8]);
  applyResolution(pp_cent7, resolutionIDs[7], resolutionIDsError[7]);
  applyResolution(pp_cent6, resolutionIDs[6], resolutionIDsError[6]);
  applyResolution(pp_cent5, resolutionIDs[5], resolutionIDsError[5]);
  applyResolution(pp_cent4, resolutionIDs[4], resolutionIDsError[4]);
  
  std::cout << "pp_cent10 bin 18 content: " << pp_cent10->GetBinContent(18) << std::endl;
  */

  /*
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
  */
  TProfile *p_vn_yCM_00to10_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_00to10_pp", 15, 16);
  TProfile *p_vn_yCM_10to40_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_10to40_pp", 9, 14);
  TProfile *p_vn_yCM_40to60_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_40to60_pp", 5, 8);
  //TProfile *p_vn_yCM_00to60_pp = (TProfile*)file->Get("p_vn_yCM_00to60_pp");
  
  TProfile *p_vn_yCM_00to10_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_00to10_pm", 15, 16);
  TProfile *p_vn_yCM_10to40_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_10to40_pm", 9, 14);
  TProfile *p_vn_yCM_40to60_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_40to60_pm", 5, 8);
  //TProfile *p_vn_yCM_00to60_pm = (TProfile*)file->Get("p_vn_yCM_00to60_pm");

  TProfile *p_vn_yCM_00to10_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_00to10_kp", 15, 16);
  TProfile *p_vn_yCM_10to40_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_10to40_kp", 9, 14);
  TProfile *p_vn_yCM_40to60_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_40to60_kp", 5, 8);
  //TProfile *p_vn_yCM_00to60_kp = (TProfile*)file->Get("p_vn_yCM_00to60_kp");

  TProfile *p_vn_yCM_00to10_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_00to10_km", 15, 16);
  TProfile *p_vn_yCM_10to40_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_10to40_km", 9, 14);
  TProfile *p_vn_yCM_40to60_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_40to60_km", 5, 8);
  //TProfile *p_vn_yCM_00to60_km = (TProfile*)file->Get("p_vn_yCM_00to60_km");

  TProfile *p_vn_yCM_00to10_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_00to10_pr", 15, 16);
  TProfile *p_vn_yCM_10to40_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_10to40_pr", 9, 14);
  TProfile *p_vn_yCM_40to60_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_40to60_pr", 5, 8);
  //TProfile *p_vn_yCM_00to60_pr = (TProfile*)file->Get("p_vn_yCM_00to60_pr");


  TH1D *h_vn_yCM_00to10_pp = new TH1D("h_vn_yCM_00to10_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pp = new TH1D("h_vn_yCM_10to40_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pp = new TH1D("h_vn_yCM_40to60_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pp = new TH1D("h_vn_yCM_00to60_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pm = new TH1D("h_vn_yCM_00to10_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pm = new TH1D("h_vn_yCM_10to40_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pm = new TH1D("h_vn_yCM_40to60_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pm = new TH1D("h_vn_yCM_00to60_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_kp = new TH1D("h_vn_yCM_00to10_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_kp = new TH1D("h_vn_yCM_10to40_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_kp = new TH1D("h_vn_yCM_40to60_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_00to60_kp = new TH1D("h_vn_yCM_00to60_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_km = new TH1D("h_vn_yCM_00to10_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_km = new TH1D("h_vn_yCM_10to40_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_km = new TH1D("h_vn_yCM_40to60_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_00to60_km = new TH1D("h_vn_yCM_00to60_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_pr = new TH1D("h_vn_yCM_00to10_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr = new TH1D("h_vn_yCM_10to40_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr = new TH1D("h_vn_yCM_40to60_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pr = new TH1D("h_vn_yCM_00to60_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  /*
  for (int i = 11; i <= 20; i++)
    {
      // 0 to 10%
      h_vn_yCM_00to10_pp->SetBinContent(i, (h2_vn_yCM_cent_pp->GetBinContent(15, i) + h2_vn_yCM_cent_pp->GetBinContent(16, i)) /2);
      h_vn_yCM_00to10_pp->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pp->GetBinError(15, i) / 2, 2) + pow(h2_vn_yCM_cent_pp->GetBinError(16, i) / 2, 2) ));

      h_vn_yCM_00to10_pm->SetBinContent(i, (h2_vn_yCM_cent_pm->GetBinContent(15, i) + h2_vn_yCM_cent_pm->GetBinContent(16, i)) /2);
      h_vn_yCM_00to10_pm->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pm->GetBinError(15, i) / 2, 2) + pow(h2_vn_yCM_cent_pm->GetBinError(16, i) / 2, 2) ));

      h_vn_yCM_00to10_kp->SetBinContent(i, (h2_vn_yCM_cent_kp->GetBinContent(15, i) + h2_vn_yCM_cent_kp->GetBinContent(16, i)) /2);
      h_vn_yCM_00to10_kp->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_kp->GetBinError(15, i) / 2, 2) + pow(h2_vn_yCM_cent_kp->GetBinError(16, i) / 2, 2) ));

      h_vn_yCM_00to10_km->SetBinContent(i, (h2_vn_yCM_cent_km->GetBinContent(15, i) + h2_vn_yCM_cent_km->GetBinContent(16, i)) /2);
      h_vn_yCM_00to10_km->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_km->GetBinError(15, i) / 2, 2) + pow(h2_vn_yCM_cent_km->GetBinError(16, i) / 2, 2) ));

      h_vn_yCM_00to10_pr->SetBinContent(i, (h2_vn_yCM_cent_pr->GetBinContent(15, i) + h2_vn_yCM_cent_pr->GetBinContent(16, i)) /2);
      h_vn_yCM_00to10_pr->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pr->GetBinError(15, i) / 2, 2) + pow(h2_vn_yCM_cent_pr->GetBinError(16, i) / 2, 2) ));

      // 10 to 40%
      h_vn_yCM_10to40_pp->SetBinContent(i, (h2_vn_yCM_cent_pp->GetBinContent(9, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(10, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(11, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(12, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(13, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(14, i)) /6);
      h_vn_yCM_10to40_pp->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pp->GetBinError(9, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(10, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(11, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(12, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(13, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(14, i) / 6, 2) ));

      h_vn_yCM_10to40_pm->SetBinContent(i, (h2_vn_yCM_cent_pm->GetBinContent(9, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(10, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(11, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(12, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(13, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(14, i)) /6);
      h_vn_yCM_10to40_pm->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pm->GetBinError(9, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(10, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(11, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(12, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(13, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(14, i) / 6, 2) ));

      h_vn_yCM_10to40_kp->SetBinContent(i, (h2_vn_yCM_cent_kp->GetBinContent(9, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(10, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(11, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(12, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(13, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(14, i)) /6);
      h_vn_yCM_10to40_kp->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_kp->GetBinError(9, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(10, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(11, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(12, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(13, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(14, i) / 6, 2) ));

      h_vn_yCM_10to40_km->SetBinContent(i, (h2_vn_yCM_cent_km->GetBinContent(9, i) +
					    h2_vn_yCM_cent_km->GetBinContent(10, i) +
					    h2_vn_yCM_cent_km->GetBinContent(11, i) +
					    h2_vn_yCM_cent_km->GetBinContent(12, i) +
					    h2_vn_yCM_cent_km->GetBinContent(13, i) +
					    h2_vn_yCM_cent_km->GetBinContent(14, i)) /6);
      h_vn_yCM_10to40_km->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_km->GetBinError(9, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(10, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(11, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(12, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(13, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(14, i) / 6, 2) ));

      h_vn_yCM_10to40_pr->SetBinContent(i, (h2_vn_yCM_cent_pr->GetBinContent(9, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(10, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(11, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(12, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(13, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(14, i)) /6);
      h_vn_yCM_10to40_pr->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pr->GetBinError(9, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(10, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(11, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(12, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(13, i) / 6, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(14, i) / 6, 2) ));

      

      // 40 to 60%
      h_vn_yCM_40to60_pp->SetBinContent(i, (h2_vn_yCM_cent_pp->GetBinContent(4, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(5, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(6, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(7, i) +
					    h2_vn_yCM_cent_pp->GetBinContent(8, i)) /5);
      h_vn_yCM_40to60_pp->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pp->GetBinError(4, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(5, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(6, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(7, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pp->GetBinError(8, i) / 5, 2) ));

      h_vn_yCM_40to60_pm->SetBinContent(i, (h2_vn_yCM_cent_pm->GetBinContent(4, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(5, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(6, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(7, i) +
					    h2_vn_yCM_cent_pm->GetBinContent(8, i)) /5);
      h_vn_yCM_40to60_pm->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pm->GetBinError(4, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(5, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(6, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(7, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pm->GetBinError(8, i) / 5, 2) ));

      h_vn_yCM_40to60_kp->SetBinContent(i, (h2_vn_yCM_cent_kp->GetBinContent(4, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(5, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(6, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(7, i) +
					    h2_vn_yCM_cent_kp->GetBinContent(8, i)) /5);
      h_vn_yCM_40to60_kp->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_kp->GetBinError(4, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(5, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(6, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(7, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_kp->GetBinError(8, i) / 5, 2) ));

      h_vn_yCM_40to60_km->SetBinContent(i, (h2_vn_yCM_cent_km->GetBinContent(4, i) +
					    h2_vn_yCM_cent_km->GetBinContent(5, i) +
					    h2_vn_yCM_cent_km->GetBinContent(6, i) +
					    h2_vn_yCM_cent_km->GetBinContent(7, i) +
					    h2_vn_yCM_cent_km->GetBinContent(8, i)) /5);
      h_vn_yCM_40to60_km->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_km->GetBinError(4, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(5, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(6, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(7, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_km->GetBinError(8, i) / 5, 2) ));

      h_vn_yCM_40to60_pr->SetBinContent(i, (h2_vn_yCM_cent_pr->GetBinContent(4, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(5, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(6, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(7, i) +
					    h2_vn_yCM_cent_pr->GetBinContent(8, i)) /5);
      h_vn_yCM_40to60_pr->SetBinError(i, TMath::Sqrt( pow(h2_vn_yCM_cent_pr->GetBinError(4, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(5, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(6, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(7, i) / 5, 2) +
						      pow(h2_vn_yCM_cent_pr->GetBinError(8, i) / 5, 2) ));

    }
  */
  /*
  h_vn_yCM_00to10_pp->Draw();
  h_vn_yCM_00to10_pp->SetMinimum(-0.05);
  h_vn_yCM_00to10_pp->SetMaximum(0.03);
  h_vn_yCM_00to10_pp->Draw("E1P");
  return;
  */


  TH1D *h_vn_EpdE = new TH1D("h_vn_EpdE", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_EpdF = new TH1D("h_vn_EpdF", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_TpcB = new TH1D("h_vn_TpcB", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_Tpc = new TH1D("h_vn_Tpc", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_pp = new TH1D("h_vn_pp", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_pm = new TH1D("h_vn_pm", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_kp = new TH1D("h_vn_kp", ";Centrality;v_{"+order_n_str+"}", 8, 0, 16);
  TH1D *h_vn_km = new TH1D("h_vn_km", ";Centrality;v_{"+order_n_str+"}", 8, 0, 16);
  TH1D *h_vn_pr = new TH1D("h_vn_pr", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);

  TH1D *h_vn_EpdE_ext = new TH1D("h_vn_EpdE_ext", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_EpdF_ext = new TH1D("h_vn_EpdF_ext", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_TpcB_ext = new TH1D("h_vn_TpcB_ext", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_pp_ext = new TH1D("h_vn_pp_ext", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_pm_ext = new TH1D("h_vn_pm_ext", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);
  TH1D *h_vn_kp_ext = new TH1D("h_vn_kp_ext", ";Centrality;v_{"+order_n_str+"}", 8, 0, 16);
  TH1D *h_vn_km_ext = new TH1D("h_vn_km_ext", ";Centrality;v_{"+order_n_str+"}", 8, 0, 16);
  TH1D *h_vn_pr_ext = new TH1D("h_vn_pr_ext", ";Centrality;v_{"+order_n_str+"}", 16, 0, 16);



  /*
  h_vn_yCM_00to10_pp->Add(pp_cent15, pp_cent14);
  
  h_vn_yCM_10to40_pp->Add(pp_cent13, pp_cent12);
  h_vn_yCM_10to40_pp->Add(h_vn_yCM_10to40_pp, pp_cent11);
  h_vn_yCM_10to40_pp->Add(h_vn_yCM_10to40_pp, pp_cent10);
  h_vn_yCM_10to40_pp->Add(h_vn_yCM_10to40_pp, pp_cent9);
  h_vn_yCM_10to40_pp->Add(h_vn_yCM_10to40_pp, pp_cent8);

  
  
  h_vn_yCM_40to60_pp->Add(pp_cent7, pp_cent6);
  h_vn_yCM_40to60_pp->Add(h_vn_yCM_40to60_pp, pp_cent5);
  h_vn_yCM_40to60_pp->Add(h_vn_yCM_40to60_pp, pp_cent4);
  */

  //mirrored plots
  TH1D *h_vn_yCM_00to10_pp_mirror = new TH1D("h_vn_yCM_00to10_pp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pp_mirror = new TH1D("h_vn_yCM_10to40_pp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pp_mirror = new TH1D("h_vn_yCM_40to60_pp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pp_mirror = new TH1D("h_vn_yCM_00to60_pp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pm_mirror = new TH1D("h_vn_yCM_00to10_pm_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pm_mirror = new TH1D("h_vn_yCM_10to40_pm_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pm_mirror = new TH1D("h_vn_yCM_40to60_pm_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pm_mirror = new TH1D("h_vn_yCM_00to60_pm_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_kp_mirror = new TH1D("h_vn_yCM_00to10_kp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_kp_mirror = new TH1D("h_vn_yCM_10to40_kp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_kp_mirror = new TH1D("h_vn_yCM_40to60_kp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_00to60_kp_mirror = new TH1D("h_vn_yCM_00to60_kp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_km_mirror = new TH1D("h_vn_yCM_00to10_km_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_km_mirror = new TH1D("h_vn_yCM_10to40_km_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_km_mirror = new TH1D("h_vn_yCM_40to60_km_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_00to60_km_mirror = new TH1D("h_vn_yCM_00to60_km_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_pr_mirror = new TH1D("h_vn_yCM_00to10_pr_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr_mirror = new TH1D("h_vn_yCM_10to40_pr_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr_mirror = new TH1D("h_vn_yCM_40to60_pr_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_00to60_pr_mirror = new TH1D("h_vn_yCM_00to60_pr_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  // Convert profiles to histograms

  for (int i = 1; i <= 8; i++)     // KAONS
    {
      h_vn_kp->SetBinContent(i, p_vn_kp->GetBinContent(i));
      h_vn_kp->SetBinError(i, p_vn_kp->GetBinError(i));

      h_vn_km->SetBinContent(i, p_vn_km->GetBinContent(i));
      h_vn_km->SetBinError(i, p_vn_km->GetBinError(i));

      h_vn_kp_ext->SetBinContent(i, p_vn_kp_ext->GetBinContent(i));
      h_vn_kp_ext->SetBinError(i, p_vn_kp_ext->GetBinError(i));

      h_vn_km_ext->SetBinContent(i, p_vn_km_ext->GetBinContent(i));
      h_vn_km_ext->SetBinError(i, p_vn_km_ext->GetBinError(i));
    }
  
  for (int i = 1; i <= 16; i++)
    {
      h_vn_EpdE->SetBinContent(i, p_vn_EpdE->GetBinContent(i));
      h_vn_EpdE->SetBinError(i, p_vn_EpdE->GetBinError(i));

      h_vn_EpdF->SetBinContent(i, p_vn_EpdF->GetBinContent(i));
      h_vn_EpdF->SetBinError(i, p_vn_EpdF->GetBinError(i));

      h_vn_TpcB->SetBinContent(i, p_vn_TpcB->GetBinContent(i));
      h_vn_TpcB->SetBinError(i, p_vn_TpcB->GetBinError(i));

      h_vn_Tpc->SetBinContent(i, p_vn_Tpc->GetBinContent(i));
      h_vn_Tpc->SetBinError(i, p_vn_Tpc->GetBinError(i));

      h_vn_pp->SetBinContent(i, p_vn_pp->GetBinContent(i));
      h_vn_pp->SetBinError(i, p_vn_pp->GetBinError(i));

      h_vn_pm->SetBinContent(i, p_vn_pm->GetBinContent(i));
      h_vn_pm->SetBinError(i, p_vn_pm->GetBinError(i));

      
      h_vn_pr->SetBinContent(i, p_vn_pr->GetBinContent(i));
      h_vn_pr->SetBinError(i, p_vn_pr->GetBinError(i));

      // Extended rapidity plots
      h_vn_pp_ext->SetBinContent(i, p_vn_pp_ext->GetBinContent(i));
      h_vn_pp_ext->SetBinError(i, p_vn_pp_ext->GetBinError(i));

      h_vn_pm_ext->SetBinContent(i, p_vn_pm_ext->GetBinContent(i));
      h_vn_pm_ext->SetBinError(i, p_vn_pm_ext->GetBinError(i));


      h_vn_pr_ext->SetBinContent(i, p_vn_pr_ext->GetBinContent(i));
      h_vn_pr_ext->SetBinError(i, p_vn_pr_ext->GetBinError(i));
    }

  // Convert profiles to histograms and make mirrored plots

  for (int i = 1; i <= 10; i++)     // KAONS
    {
      int j = 0;  // mirrored bin
      
      switch(i)
	{
	case 6: j = 5; break;
	case 7: j = 4; break;
	case 8: j = 3; break;
	case 9: j = 2; break;
	case 10: j = 1; break;
	}

      h_vn_yCM_00to10_kp->SetBinContent(i, p_vn_yCM_00to10_kp->GetBinContent(i));
      h_vn_yCM_00to10_kp->SetBinError(i, p_vn_yCM_00to10_kp->GetBinError(i));
      h_vn_yCM_10to40_kp->SetBinContent(i, p_vn_yCM_10to40_kp->GetBinContent(i));
      h_vn_yCM_10to40_kp->SetBinError(i, p_vn_yCM_10to40_kp->GetBinError(i));
      h_vn_yCM_40to60_kp->SetBinContent(i, p_vn_yCM_40to60_kp->GetBinContent(i));
      h_vn_yCM_40to60_kp->SetBinError(i, p_vn_yCM_40to60_kp->GetBinError(i));
      //h_vn_yCM_00to60_kp->SetBinContent(i, p_vn_yCM_00to60_kp->GetBinContent(i));
      //h_vn_yCM_00to60_kp->SetBinError(i, p_vn_yCM_00to60_kp->GetBinError(i));

      h_vn_yCM_00to10_km->SetBinContent(i, p_vn_yCM_00to10_km->GetBinContent(i));
      h_vn_yCM_00to10_km->SetBinError(i, p_vn_yCM_00to10_km->GetBinError(i));
      h_vn_yCM_10to40_km->SetBinContent(i, p_vn_yCM_10to40_km->GetBinContent(i));
      h_vn_yCM_10to40_km->SetBinError(i, p_vn_yCM_10to40_km->GetBinError(i));
      h_vn_yCM_40to60_km->SetBinContent(i, p_vn_yCM_40to60_km->GetBinContent(i));
      h_vn_yCM_40to60_km->SetBinError(i, p_vn_yCM_40to60_km->GetBinError(i));
      //h_vn_yCM_00to60_km->SetBinContent(i, p_vn_yCM_00to60_km->GetBinContent(i));
      //h_vn_yCM_00to60_km->SetBinError(i, p_vn_yCM_00to60_km->GetBinError(i));

      if (j!=0)
	{
	  h_vn_yCM_00to10_kp_mirror->SetBinContent(j, h_vn_yCM_00to10_kp->GetBinContent(i));
	  h_vn_yCM_00to10_kp_mirror->SetBinError(j, h_vn_yCM_00to10_kp->GetBinError(i));
	  h_vn_yCM_10to40_kp_mirror->SetBinContent(j, h_vn_yCM_10to40_kp->GetBinContent(i));
	  h_vn_yCM_10to40_kp_mirror->SetBinError(j, h_vn_yCM_10to40_kp->GetBinError(i));
	  h_vn_yCM_40to60_kp_mirror->SetBinContent(j, h_vn_yCM_40to60_kp->GetBinContent(i));
	  h_vn_yCM_40to60_kp_mirror->SetBinError(j, h_vn_yCM_40to60_kp->GetBinError(i));
	  //h_vn_yCM_00to60_kp_mirror->SetBinContent(j, h_vn_yCM_00to60_kp->GetBinContent(i));
	  //h_vn_yCM_00to60_kp_mirror->SetBinError(j, h_vn_yCM_00to60_kp->GetBinError(i));

	  h_vn_yCM_00to10_km_mirror->SetBinContent(j, h_vn_yCM_00to10_km->GetBinContent(i));
	  h_vn_yCM_00to10_km_mirror->SetBinError(j, h_vn_yCM_00to10_km->GetBinError(i));
	  h_vn_yCM_10to40_km_mirror->SetBinContent(j, h_vn_yCM_10to40_km->GetBinContent(i));
	  h_vn_yCM_10to40_km_mirror->SetBinError(j, h_vn_yCM_10to40_km->GetBinError(i));
	  h_vn_yCM_40to60_km_mirror->SetBinContent(j, h_vn_yCM_40to60_km->GetBinContent(i));
	  h_vn_yCM_40to60_km_mirror->SetBinError(j, h_vn_yCM_40to60_km->GetBinError(i));
	  //h_vn_yCM_00to60_km_mirror->SetBinContent(j, h_vn_yCM_00to60_km->GetBinContent(i));
	  //h_vn_yCM_00to60_km_mirror->SetBinError(j, h_vn_yCM_00to60_km->GetBinError(i));
	}
    }
  
  for (int i = 1; i <= 20; i++)
    {
      int j = 0;  // mirrored bin
      
      switch(i)
	{
	case 11: j = 10; break;
	case 12: j = 9; break;
	case 13: j = 8; break;
	case 14: j = 7; break;
	case 15: j = 6; break;
	case 16: j = 5; break;
	case 17: j = 4; break;
	case 18: j = 3; break;
	case 19: j = 2; break;
	case 20: j = 1; break;
	}

      h_vn_yCM_00to10_pp->SetBinContent(i, p_vn_yCM_00to10_pp->GetBinContent(i));
      h_vn_yCM_00to10_pp->SetBinError(i, p_vn_yCM_00to10_pp->GetBinError(i));
      h_vn_yCM_10to40_pp->SetBinContent(i, p_vn_yCM_10to40_pp->GetBinContent(i));
      h_vn_yCM_10to40_pp->SetBinError(i, p_vn_yCM_10to40_pp->GetBinError(i));
      h_vn_yCM_40to60_pp->SetBinContent(i, p_vn_yCM_40to60_pp->GetBinContent(i));
      h_vn_yCM_40to60_pp->SetBinError(i, p_vn_yCM_40to60_pp->GetBinError(i));
      //h_vn_yCM_00to60_pp->SetBinContent(i, p_vn_yCM_00to60_pp->GetBinContent(i));
      //h_vn_yCM_00to60_pp->SetBinError(i, p_vn_yCM_00to60_pp->GetBinError(i));

      h_vn_yCM_00to10_pm->SetBinContent(i, p_vn_yCM_00to10_pm->GetBinContent(i));
      h_vn_yCM_00to10_pm->SetBinError(i, p_vn_yCM_00to10_pm->GetBinError(i));
      h_vn_yCM_10to40_pm->SetBinContent(i, p_vn_yCM_10to40_pm->GetBinContent(i));
      h_vn_yCM_10to40_pm->SetBinError(i, p_vn_yCM_10to40_pm->GetBinError(i));
      h_vn_yCM_40to60_pm->SetBinContent(i, p_vn_yCM_40to60_pm->GetBinContent(i));
      h_vn_yCM_40to60_pm->SetBinError(i, p_vn_yCM_40to60_pm->GetBinError(i));
      //h_vn_yCM_00to60_pm->SetBinContent(i, p_vn_yCM_00to60_pm->GetBinContent(i));
      //h_vn_yCM_00to60_pm->SetBinError(i, p_vn_yCM_00to60_pm->GetBinError(i));


      h_vn_yCM_00to10_pr->SetBinContent(i, p_vn_yCM_00to10_pr->GetBinContent(i));
      h_vn_yCM_00to10_pr->SetBinError(i, p_vn_yCM_00to10_pr->GetBinError(i));
      h_vn_yCM_10to40_pr->SetBinContent(i, p_vn_yCM_10to40_pr->GetBinContent(i));
      h_vn_yCM_10to40_pr->SetBinError(i, p_vn_yCM_10to40_pr->GetBinError(i));
      h_vn_yCM_40to60_pr->SetBinContent(i, p_vn_yCM_40to60_pr->GetBinContent(i));
      h_vn_yCM_40to60_pr->SetBinError(i, p_vn_yCM_40to60_pr->GetBinError(i));
      //h_vn_yCM_00to60_pr->SetBinContent(i, p_vn_yCM_00to60_pr->GetBinContent(i));
      //h_vn_yCM_00to60_pr->SetBinError(i, p_vn_yCM_00to60_pr->GetBinError(i));


      if(j != 0)
	{
	  h_vn_yCM_00to10_pp_mirror->SetBinContent(j, h_vn_yCM_00to10_pp->GetBinContent(i));
	  h_vn_yCM_00to10_pp_mirror->SetBinError(j, h_vn_yCM_00to10_pp->GetBinError(i));
	  h_vn_yCM_10to40_pp_mirror->SetBinContent(j, h_vn_yCM_10to40_pp->GetBinContent(i));
	  h_vn_yCM_10to40_pp_mirror->SetBinError(j, h_vn_yCM_10to40_pp->GetBinError(i));
	  h_vn_yCM_40to60_pp_mirror->SetBinContent(j, h_vn_yCM_40to60_pp->GetBinContent(i));
	  h_vn_yCM_40to60_pp_mirror->SetBinError(j, h_vn_yCM_40to60_pp->GetBinError(i));
	  //h_vn_yCM_00to60_pp_mirror->SetBinContent(j, h_vn_yCM_00to60_pp->GetBinContent(i));
	  //h_vn_yCM_00to60_pp_mirror->SetBinError(j, h_vn_yCM_00to60_pp->GetBinError(i));

	  h_vn_yCM_00to10_pm_mirror->SetBinContent(j, h_vn_yCM_00to10_pm->GetBinContent(i));
	  h_vn_yCM_00to10_pm_mirror->SetBinError(j, h_vn_yCM_00to10_pm->GetBinError(i));
	  h_vn_yCM_10to40_pm_mirror->SetBinContent(j, h_vn_yCM_10to40_pm->GetBinContent(i));
	  h_vn_yCM_10to40_pm_mirror->SetBinError(j, h_vn_yCM_10to40_pm->GetBinError(i));
	  h_vn_yCM_40to60_pm_mirror->SetBinContent(j, h_vn_yCM_40to60_pm->GetBinContent(i));
	  h_vn_yCM_40to60_pm_mirror->SetBinError(j, h_vn_yCM_40to60_pm->GetBinError(i));
	  //h_vn_yCM_00to60_pm_mirror->SetBinContent(j, h_vn_yCM_00to60_pm->GetBinContent(i));
	  //h_vn_yCM_00to60_pm_mirror->SetBinError(j, h_vn_yCM_00to60_pm->GetBinError(i));


	  h_vn_yCM_00to10_pr_mirror->SetBinContent(j, h_vn_yCM_00to10_pr->GetBinContent(i));
	  h_vn_yCM_00to10_pr_mirror->SetBinError(j, h_vn_yCM_00to10_pr->GetBinError(i));
	  h_vn_yCM_10to40_pr_mirror->SetBinContent(j, h_vn_yCM_10to40_pr->GetBinContent(i));
	  h_vn_yCM_10to40_pr_mirror->SetBinError(j, h_vn_yCM_10to40_pr->GetBinError(i));
	  h_vn_yCM_40to60_pr_mirror->SetBinContent(j, h_vn_yCM_40to60_pr->GetBinContent(i));
	  h_vn_yCM_40to60_pr_mirror->SetBinError(j, h_vn_yCM_40to60_pr->GetBinError(i));
	  //h_vn_yCM_00to60_pr_mirror->SetBinContent(j, h_vn_yCM_00to60_pr->GetBinContent(i));
	  //h_vn_yCM_00to60_pr_mirror->SetBinError(j, h_vn_yCM_00to60_pr->GetBinError(i));
	}
    }


  /*
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

  TH1D *h_vn_Tpc_flip = new TH1D("h_vn_Tpc_flip",h_vn_Tpc->GetTitle(),centBins,0,centBins);
  h_vn_Tpc_flip->GetXaxis()->SetTitle((TString)h_vn_Tpc->GetXaxis()->GetTitle()+" (%)");
  h_vn_Tpc_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_Tpc->GetYaxis()->GetTitle());

  TH1D *h_vn_pp_flip = new TH1D("h_vn_pp_flip",h_vn_pp->GetTitle(),centBins,0,centBins);
  h_vn_pp_flip->GetXaxis()->SetTitle((TString)h_vn_pp->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pp->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_flip = new TH1D("h_vn_pm_flip",h_vn_pm->GetTitle(),centBins,0,centBins);
  h_vn_pm_flip->GetXaxis()->SetTitle((TString)h_vn_pm->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pm->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_flip = new TH1D("h_vn_kp_flip",h_vn_kp->GetTitle(),centBins / 2,0,centBins);
  h_vn_kp_flip->GetXaxis()->SetTitle((TString)h_vn_kp->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_kp->GetYaxis()->GetTitle());

  TH1D *h_vn_km_flip = new TH1D("h_vn_km_flip",h_vn_km->GetTitle(),centBins / 2,0,centBins);
  h_vn_km_flip->GetXaxis()->SetTitle((TString)h_vn_km->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_km->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_flip = new TH1D("h_vn_pr_flip",h_vn_pr->GetTitle(),centBins,0,centBins);
  h_vn_pr_flip->GetXaxis()->SetTitle((TString)h_vn_pr->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pr->GetYaxis()->GetTitle());

  // Extended rapidity plots
  TH1D *h_vn_pp_ext_flip = new TH1D("h_vn_pp_ext_flip",h_vn_pp_ext->GetTitle(),centBins,0,centBins);
  h_vn_pp_ext_flip->GetXaxis()->SetTitle((TString)h_vn_pp_ext->GetXaxis()->GetTitle()+" (%)");
  h_vn_pp_ext_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pp_ext->GetYaxis()->GetTitle());

  TH1D *h_vn_pm_ext_flip = new TH1D("h_vn_pm_ext_flip",h_vn_pm_ext->GetTitle(),centBins,0,centBins);
  h_vn_pm_ext_flip->GetXaxis()->SetTitle((TString)h_vn_pm_ext->GetXaxis()->GetTitle()+" (%)");
  h_vn_pm_ext_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pm_ext->GetYaxis()->GetTitle());

  TH1D *h_vn_kp_ext_flip = new TH1D("h_vn_kp_ext_flip",h_vn_kp_ext->GetTitle(),centBins / 2,0,centBins);
  h_vn_kp_ext_flip->GetXaxis()->SetTitle((TString)h_vn_kp_ext->GetXaxis()->GetTitle()+" (%)");
  h_vn_kp_ext_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_kp_ext->GetYaxis()->GetTitle());

  TH1D *h_vn_km_ext_flip = new TH1D("h_vn_km_ext_flip",h_vn_km_ext->GetTitle(),centBins / 2,0,centBins);
  h_vn_km_ext_flip->GetXaxis()->SetTitle((TString)h_vn_km_ext->GetXaxis()->GetTitle()+" (%)");
  h_vn_km_ext_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_km_ext->GetYaxis()->GetTitle());

  TH1D *h_vn_pr_ext_flip = new TH1D("h_vn_pr_ext_flip",h_vn_pr_ext->GetTitle(),centBins,0,centBins);
  h_vn_pr_ext_flip->GetXaxis()->SetTitle((TString)h_vn_pr_ext->GetXaxis()->GetTitle()+" (%)");
  h_vn_pr_ext_flip->GetYaxis()->SetTitle("v_{"+order_n_str+"}");//h_vn_pr_ext->GetYaxis()->GetTitle());


  const char *centralityBins[16] = {"75-80", "70-75", "65-70", "60-65", "55-60", "50-55", "45-50", "40-45", "35-40", "30-35", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5"};


  std::vector<TString> newBinLabels;

  // Get list of bin labels, but flipped
  for (int i = lastCentID; i >= firstCentID; i--) { newBinLabels.push_back(centralityBins[i]); }

  // Flip the bin contents into the new histograms
  int j = 1;
  for (int i = centBins / 2; i >= 1; i--)     // KAONS
    {
      h_vn_kp_flip->SetBinContent(j, h_vn_kp->GetBinContent(i));
      h_vn_kp_flip->SetBinError(j, h_vn_kp->GetBinError(i));

      h_vn_km_flip->SetBinContent(j, h_vn_km->GetBinContent(i));
      h_vn_km_flip->SetBinError(j, h_vn_km->GetBinError(i));

      h_vn_kp_ext_flip->SetBinContent(j, h_vn_kp_ext->GetBinContent(i));
      h_vn_kp_ext_flip->SetBinError(j, h_vn_kp_ext->GetBinError(i));

      h_vn_km_ext_flip->SetBinContent(j, h_vn_km_ext->GetBinContent(i));
      h_vn_km_ext_flip->SetBinError(j, h_vn_km_ext->GetBinError(i));

      j++;
    }

  j = 1;
  for (int i = centBins; i >= 1; i--)
    {
      h_vn_EpdE_flip->SetBinContent(j, h_vn_EpdE->GetBinContent(i));
      h_vn_EpdE_flip->SetBinError(j, h_vn_EpdE->GetBinError(i));

      h_vn_EpdF_flip->SetBinContent(j, h_vn_EpdF->GetBinContent(i));
      h_vn_EpdF_flip->SetBinError(j, h_vn_EpdF->GetBinError(i));

      h_vn_TpcB_flip->SetBinContent(j, h_vn_TpcB->GetBinContent(i));
      h_vn_TpcB_flip->SetBinError(j, h_vn_TpcB->GetBinError(i));

      h_vn_Tpc_flip->SetBinContent(j, h_vn_Tpc->GetBinContent(i));
      h_vn_Tpc_flip->SetBinError(j, h_vn_Tpc->GetBinError(i));

      h_vn_pp_flip->SetBinContent(j, h_vn_pp->GetBinContent(i));
      h_vn_pp_flip->SetBinError(j, h_vn_pp->GetBinError(i));

      h_vn_pm_flip->SetBinContent(j, h_vn_pm->GetBinContent(i));
      h_vn_pm_flip->SetBinError(j, h_vn_pm->GetBinError(i));


      h_vn_pr_flip->SetBinContent(j, h_vn_pr->GetBinContent(i));
      h_vn_pr_flip->SetBinError(j, h_vn_pr->GetBinError(i));


      h_vn_pp_ext_flip->SetBinContent(j, h_vn_pp_ext->GetBinContent(i));
      h_vn_pp_ext_flip->SetBinError(j, h_vn_pp_ext->GetBinError(i));

      h_vn_pm_ext_flip->SetBinContent(j, h_vn_pm_ext->GetBinContent(i));
      h_vn_pm_ext_flip->SetBinError(j, h_vn_pm_ext->GetBinError(i));


      h_vn_pr_ext_flip->SetBinContent(j, h_vn_pr_ext->GetBinContent(i));
      h_vn_pr_ext_flip->SetBinError(j, h_vn_pr_ext->GetBinError(i));

      
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
  TH1D *vn_Tpc = new TH1D("vn_Tpc", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_pp = new TH1D("vn_pp", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_pm = new TH1D("vn_pm", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_kp = new TH1D("vn_kp", ";Centrality (%);v_{"+order_n_str+"}", 6, 0, 60);
  TH1D *vn_km = new TH1D("vn_km", ";Centrality (%);v_{"+order_n_str+"}", 6, 0, 60);
  TH1D *vn_pr = new TH1D("vn_pr", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_pp_ext = new TH1D("vn_pp_ext", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_pm_ext = new TH1D("vn_pm_ext", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);
  TH1D *vn_kp_ext = new TH1D("vn_kp_ext", ";Centrality (%);v_{"+order_n_str+"}", 6, 0, 60);
  TH1D *vn_km_ext = new TH1D("vn_km_ext", ";Centrality (%);v_{"+order_n_str+"}", 6, 0, 60);
  TH1D *vn_pr_ext = new TH1D("vn_pr_ext", ";Centrality (%);v_{"+order_n_str+"}", 12, 0, 60);

  for (int i = 1; i <= 6; i++)     // KAONS
    {
      vn_kp->SetBinContent(i, h_vn_kp_flip->GetBinContent(i));
      vn_kp->SetBinError(i, h_vn_kp_flip->GetBinError(i));

      vn_km->SetBinContent(i, h_vn_km_flip->GetBinContent(i));
      vn_km->SetBinError(i, h_vn_km_flip->GetBinError(i));

      vn_kp_ext->SetBinContent(i, h_vn_kp_ext_flip->GetBinContent(i));
      vn_kp_ext->SetBinError(i, h_vn_kp_ext_flip->GetBinError(i));

      vn_km_ext->SetBinContent(i, h_vn_km_ext_flip->GetBinContent(i));
      vn_km_ext->SetBinError(i, h_vn_km_ext_flip->GetBinError(i));
    }

  for (int i = 1; i <= 12; i++)
    {
      vn_EpdE->SetBinContent(i, h_vn_EpdE_flip->GetBinContent(i));
      vn_EpdE->SetBinError(i, h_vn_EpdE_flip->GetBinError(i));

      vn_EpdF->SetBinContent(i, h_vn_EpdF_flip->GetBinContent(i));
      vn_EpdF->SetBinError(i, h_vn_EpdF_flip->GetBinError(i));

      vn_TpcB->SetBinContent(i, h_vn_TpcB_flip->GetBinContent(i));
      vn_TpcB->SetBinError(i, h_vn_TpcB_flip->GetBinError(i));

      vn_Tpc->SetBinContent(i, h_vn_Tpc_flip->GetBinContent(i));
      vn_Tpc->SetBinError(i, h_vn_Tpc_flip->GetBinError(i));

      vn_pp->SetBinContent(i, h_vn_pp_flip->GetBinContent(i));
      vn_pp->SetBinError(i, h_vn_pp_flip->GetBinError(i));

      vn_pm->SetBinContent(i, h_vn_pm_flip->GetBinContent(i));
      vn_pm->SetBinError(i, h_vn_pm_flip->GetBinError(i));

      

      vn_pr->SetBinContent(i, h_vn_pr_flip->GetBinContent(i));
      vn_pr->SetBinError(i, h_vn_pr_flip->GetBinError(i));

      vn_pp_ext->SetBinContent(i, h_vn_pp_ext_flip->GetBinContent(i));
      vn_pp_ext->SetBinError(i, h_vn_pp_ext_flip->GetBinError(i));

      vn_pm_ext->SetBinContent(i, h_vn_pm_ext_flip->GetBinContent(i));
      vn_pm_ext->SetBinError(i, h_vn_pm_ext_flip->GetBinError(i));

      

      vn_pr_ext->SetBinContent(i, h_vn_pr_ext_flip->GetBinContent(i));
      vn_pr_ext->SetBinError(i, h_vn_pr_ext_flip->GetBinError(i));
    }
  

  
  THStack *piCentralityStack = new THStack("piCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kaCentralityStack = new THStack("kaCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");

  THStack *ppExtCentralityStack = new THStack("ppExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *pmExtCentralityStack = new THStack("pmExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kpExtCentralityStack = new THStack("kpExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *kmExtCentralityStack = new THStack("kmExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  THStack *prExtCentralityStack = new THStack("prExtCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");

  THStack *etaRegionStack    = new THStack("etaRegionStack", ";Centrality (%);v_{"+order_n_str+"}");
    
  THStack *ppRapidityStack   = new THStack("ppRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *pmRapidityStack   = new THStack("pmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *kpRapidityStack   = new THStack("kpRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *kmRapidityStack   = new THStack("kmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *prRapidityStack   = new THStack("prRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");


  sh_cent_pp->SetMarkerStyle(25);
  sh_cent_pp->SetMarkerSize(2);
  sh_cent_pp->SetMarkerColor(2);
  sh_cent_pp->SetLineColor(2);

  sh_cent_pm->SetMarkerStyle(25);
  sh_cent_pm->SetMarkerSize(2);
  sh_cent_pm->SetMarkerColor(4);
  sh_cent_pm->SetLineColor(4);

  sh_cent_kp->SetMarkerStyle(25);
  sh_cent_kp->SetMarkerSize(2);
  sh_cent_kp->SetMarkerColor(2);
  sh_cent_kp->SetLineColor(2);

  sh_cent_km->SetMarkerStyle(25);
  sh_cent_km->SetMarkerSize(2);
  sh_cent_km->SetMarkerColor(4);
  sh_cent_km->SetLineColor(4);
  
  sh_cent_pr->SetMarkerStyle(25);
  sh_cent_pr->SetMarkerSize(2);
  sh_cent_pr->SetMarkerColor(2);
  sh_cent_pr->SetLineColor(2);

  sh_y_pp->SetMarkerStyle(25);
  sh_y_pp->SetMarkerSize(2);
  sh_y_pp->SetMarkerColor(4);
  sh_y_pp->SetLineColor(4);

  sh_y_pm->SetMarkerStyle(25);
  sh_y_pm->SetMarkerSize(2);
  sh_y_pm->SetMarkerColor(4);
  sh_y_pm->SetLineColor(4);

  sh_y_pr->SetMarkerStyle(25);
  sh_y_pr->SetMarkerSize(2);
  sh_y_pr->SetMarkerColor(4);
  sh_y_pr->SetLineColor(4);
      
  sh_y_kp->SetMarkerStyle(25);
  sh_y_kp->SetMarkerSize(2);
  sh_y_kp->SetMarkerColor(4);
  sh_y_kp->SetLineColor(4);

  sh_y_km->SetMarkerStyle(25);
  sh_y_km->SetMarkerSize(2);
  sh_y_km->SetMarkerColor(4);
  sh_y_km->SetLineColor(4);

  
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


  vn_pp_ext->SetMarkerStyle(20);
  vn_pp_ext->SetMarkerSize(2);

  vn_pm_ext->SetMarkerStyle(20);
  vn_pm_ext->SetMarkerSize(2);

  vn_kp_ext->SetMarkerStyle(20);
  vn_kp_ext->SetMarkerSize(2);

  vn_km_ext->SetMarkerStyle(20);
  vn_km_ext->SetMarkerSize(2);

  vn_pr_ext->SetMarkerStyle(20);
  vn_pr_ext->SetMarkerSize(2);

  
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

  vn_Tpc->SetMarkerStyle(20);
  vn_Tpc->SetMarkerSize(2);
  //vn_Tpc->SetMarkerColor(8);
  //vn_Tpc->SetLineColor(8);

  h_vn_yCM_00to10_pp->SetMarkerStyle(20);
  h_vn_yCM_10to40_pp->SetMarkerStyle(20);
  h_vn_yCM_40to60_pp->SetMarkerStyle(20);
  h_vn_yCM_00to60_pp->SetMarkerStyle(20);
  h_vn_yCM_00to10_pp->SetMarkerColor(2);
  h_vn_yCM_10to40_pp->SetMarkerColor(4);
  h_vn_yCM_40to60_pp->SetMarkerColor(8);
  h_vn_yCM_00to60_pp->SetMarkerColor(4);
  h_vn_yCM_00to10_pp->SetMarkerSize(2);
  h_vn_yCM_10to40_pp->SetMarkerSize(2);
  h_vn_yCM_40to60_pp->SetMarkerSize(2);
  h_vn_yCM_00to60_pp->SetMarkerSize(2);
  h_vn_yCM_00to10_pp->SetLineColor(2);
  h_vn_yCM_10to40_pp->SetLineColor(4);
  h_vn_yCM_40to60_pp->SetLineColor(8);
  h_vn_yCM_00to60_pp->SetLineColor(4);

  h_vn_yCM_00to10_pm->SetMarkerStyle(20);
  h_vn_yCM_10to40_pm->SetMarkerStyle(20);
  h_vn_yCM_40to60_pm->SetMarkerStyle(20);
  h_vn_yCM_00to60_pm->SetMarkerStyle(20);
  h_vn_yCM_00to10_pm->SetMarkerColor(2);
  h_vn_yCM_10to40_pm->SetMarkerColor(4);
  h_vn_yCM_40to60_pm->SetMarkerColor(8);
  h_vn_yCM_00to60_pm->SetMarkerColor(4);
  h_vn_yCM_00to10_pm->SetMarkerSize(2);
  h_vn_yCM_10to40_pm->SetMarkerSize(2);
  h_vn_yCM_40to60_pm->SetMarkerSize(2);
  h_vn_yCM_00to60_pm->SetMarkerSize(2);
  h_vn_yCM_00to10_pm->SetLineColor(2);
  h_vn_yCM_10to40_pm->SetLineColor(4);
  h_vn_yCM_40to60_pm->SetLineColor(8);
  h_vn_yCM_00to60_pm->SetLineColor(4);

  h_vn_yCM_00to10_kp->SetMarkerStyle(20);
  h_vn_yCM_10to40_kp->SetMarkerStyle(20);
  h_vn_yCM_40to60_kp->SetMarkerStyle(20);
  h_vn_yCM_00to60_kp->SetMarkerStyle(20);
  h_vn_yCM_00to10_kp->SetMarkerColor(2);
  h_vn_yCM_10to40_kp->SetMarkerColor(4);
  h_vn_yCM_40to60_kp->SetMarkerColor(8);
  h_vn_yCM_00to60_kp->SetMarkerColor(4);
  h_vn_yCM_00to10_kp->SetMarkerSize(2);
  h_vn_yCM_10to40_kp->SetMarkerSize(2);
  h_vn_yCM_40to60_kp->SetMarkerSize(2);
  h_vn_yCM_00to60_kp->SetMarkerSize(2);
  h_vn_yCM_00to10_kp->SetLineColor(2);
  h_vn_yCM_10to40_kp->SetLineColor(4);
  h_vn_yCM_40to60_kp->SetLineColor(8);
  h_vn_yCM_00to60_kp->SetLineColor(4);

  h_vn_yCM_00to10_km->SetMarkerStyle(20);
  h_vn_yCM_10to40_km->SetMarkerStyle(20);
  h_vn_yCM_40to60_km->SetMarkerStyle(20);
  h_vn_yCM_00to60_km->SetMarkerStyle(20);
  h_vn_yCM_00to10_km->SetMarkerColor(2);
  h_vn_yCM_10to40_km->SetMarkerColor(4);
  h_vn_yCM_40to60_km->SetMarkerColor(8);
  h_vn_yCM_00to60_km->SetMarkerColor(4);
  h_vn_yCM_00to10_km->SetMarkerSize(2);
  h_vn_yCM_10to40_km->SetMarkerSize(2);
  h_vn_yCM_40to60_km->SetMarkerSize(2);
  h_vn_yCM_00to60_km->SetMarkerSize(2);
  h_vn_yCM_00to10_km->SetLineColor(2);
  h_vn_yCM_10to40_km->SetLineColor(4);
  h_vn_yCM_40to60_km->SetLineColor(8);
  h_vn_yCM_00to60_km->SetLineColor(4);

  h_vn_yCM_00to10_pr->SetMarkerStyle(20);
  h_vn_yCM_10to40_pr->SetMarkerStyle(20);
  h_vn_yCM_40to60_pr->SetMarkerStyle(20);
  h_vn_yCM_00to60_pr->SetMarkerStyle(20);
  h_vn_yCM_00to10_pr->SetMarkerColor(2);
  h_vn_yCM_10to40_pr->SetMarkerColor(4);
  h_vn_yCM_40to60_pr->SetMarkerColor(8);
  h_vn_yCM_00to60_pr->SetMarkerColor(4);
  h_vn_yCM_00to10_pr->SetMarkerSize(2);
  h_vn_yCM_10to40_pr->SetMarkerSize(2);
  h_vn_yCM_40to60_pr->SetMarkerSize(2);
  h_vn_yCM_00to60_pr->SetMarkerSize(2);
  h_vn_yCM_00to10_pr->SetLineColor(2);
  h_vn_yCM_10to40_pr->SetLineColor(4);
  h_vn_yCM_40to60_pr->SetLineColor(8);
  h_vn_yCM_00to60_pr->SetLineColor(4);

  //mirrored plots
  h_vn_yCM_00to10_pp_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_pp_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_pp_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to60_pp_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_pp_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_pp_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_pp_mirror->SetMarkerColor(8);
  h_vn_yCM_00to60_pp_mirror->SetMarkerColor(4);
  h_vn_yCM_00to10_pp_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_pp_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_pp_mirror->SetMarkerSize(2);
  h_vn_yCM_00to60_pp_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_pp_mirror->SetLineColor(2);
  h_vn_yCM_10to40_pp_mirror->SetLineColor(4);
  h_vn_yCM_40to60_pp_mirror->SetLineColor(8);
  h_vn_yCM_00to60_pp_mirror->SetLineColor(4);

  h_vn_yCM_00to10_pm_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_pm_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_pm_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to60_pm_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_pm_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_pm_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_pm_mirror->SetMarkerColor(8);
  h_vn_yCM_00to60_pm_mirror->SetMarkerColor(4);
  h_vn_yCM_00to10_pm_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_pm_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_pm_mirror->SetMarkerSize(2);
  h_vn_yCM_00to60_pm_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_pm_mirror->SetLineColor(2);
  h_vn_yCM_10to40_pm_mirror->SetLineColor(4);
  h_vn_yCM_40to60_pm_mirror->SetLineColor(8);
  h_vn_yCM_00to60_pm_mirror->SetLineColor(4);

  h_vn_yCM_00to10_kp_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_kp_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_kp_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to60_kp_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_kp_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_kp_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_kp_mirror->SetMarkerColor(8);
  h_vn_yCM_00to60_kp_mirror->SetMarkerColor(4);
  h_vn_yCM_00to10_kp_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_kp_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_kp_mirror->SetMarkerSize(2);
  h_vn_yCM_00to60_kp_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_kp_mirror->SetLineColor(2);
  h_vn_yCM_10to40_kp_mirror->SetLineColor(4);
  h_vn_yCM_40to60_kp_mirror->SetLineColor(8);
  h_vn_yCM_00to60_kp_mirror->SetLineColor(4);

  h_vn_yCM_00to10_km_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_km_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_km_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to60_km_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_km_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_km_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_km_mirror->SetMarkerColor(8);
  h_vn_yCM_00to60_km_mirror->SetMarkerColor(4);
  h_vn_yCM_00to10_km_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_km_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_km_mirror->SetMarkerSize(2);
  h_vn_yCM_00to60_km_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_km_mirror->SetLineColor(2);
  h_vn_yCM_10to40_km_mirror->SetLineColor(4);
  h_vn_yCM_40to60_km_mirror->SetLineColor(8);
  h_vn_yCM_00to60_km_mirror->SetLineColor(4);

  h_vn_yCM_00to10_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to60_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_pr_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_pr_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_pr_mirror->SetMarkerColor(8);
  h_vn_yCM_00to60_pr_mirror->SetMarkerColor(4);
  h_vn_yCM_00to10_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_00to60_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_pr_mirror->SetLineColor(2);
  h_vn_yCM_10to40_pr_mirror->SetLineColor(4);
  h_vn_yCM_40to60_pr_mirror->SetLineColor(8);
  h_vn_yCM_00to60_pr_mirror->SetLineColor(4);

  
  piCentralityStack->Add(vn_pp);
  piCentralityStack->Add(vn_pm);

  kaCentralityStack->Add(vn_kp);
  kaCentralityStack->Add(vn_km);


  ppExtCentralityStack->Add(vn_pp);
  ppExtCentralityStack->Add(vn_pp_ext);

  pmExtCentralityStack->Add(vn_pm);
  pmExtCentralityStack->Add(vn_pm_ext);

  kpExtCentralityStack->Add(vn_kp);
  kpExtCentralityStack->Add(vn_kp_ext);

  kmExtCentralityStack->Add(vn_km);
  kmExtCentralityStack->Add(vn_km_ext);

  prExtCentralityStack->Add(vn_pr);
  prExtCentralityStack->Add(vn_pr_ext);


  etaRegionStack->Add(vn_EpdE);
  etaRegionStack->Add(vn_EpdF);
  etaRegionStack->Add(vn_TpcB);



  if (order_n_str == "2")
    {
      ppRapidityStack->Add(h_vn_yCM_00to10_pp);
      ppRapidityStack->Add(h_vn_yCM_10to40_pp);
      ppRapidityStack->Add(h_vn_yCM_40to60_pp);

      pmRapidityStack->Add(h_vn_yCM_00to10_pm);
      pmRapidityStack->Add(h_vn_yCM_10to40_pm);
      pmRapidityStack->Add(h_vn_yCM_40to60_pm);

      kpRapidityStack->Add(h_vn_yCM_00to10_kp);
      kpRapidityStack->Add(h_vn_yCM_10to40_kp);
      kpRapidityStack->Add(h_vn_yCM_40to60_kp);

      //kmRapidityStack->Add(h_vn_yCM_00to10_km);
      kmRapidityStack->Add(h_vn_yCM_10to40_km);
      //kmRapidityStack->Add(h_vn_yCM_40to60_km);

      prRapidityStack->Add(h_vn_yCM_00to10_pr);
      prRapidityStack->Add(h_vn_yCM_10to40_pr);
      prRapidityStack->Add(h_vn_yCM_40to60_pr);


      TLegend *piLegend = new TLegend(0.775, 0.7, 0.9, 0.85);
      piLegend->AddEntry(vn_pp,"#pi^{+}");
      piLegend->AddEntry(vn_pm,"#pi^{-}");
      piLegend->SetFillColorAlpha(0,0);
      piLegend->SetLineColorAlpha(0,0);

      TLegend *kaLegend = new TLegend(0.275, 0.26, 0.425, 0.39);
      kaLegend->AddEntry(vn_kp,"K^{+}");
      kaLegend->AddEntry(vn_km,"K^{-}");
      kaLegend->SetFillColorAlpha(0,0);
      kaLegend->SetLineColorAlpha(0,0);


      TLegend *ppExtLegend = new TLegend(0.55, 0.65, 0.85, 0.88);
      ppExtLegend->AddEntry(vn_pp,"#pi^{+}, 0 < y_{CM} < 0.5");
      ppExtLegend->AddEntry(vn_pp_ext,"#pi^{+}, 0.5 < y_{CM} < 1.0");
      ppExtLegend->SetFillColorAlpha(0,0);
      ppExtLegend->SetLineColorAlpha(0,0);

      TLegend *pmExtLegend = new TLegend(0.55, 0.65, 0.85, 0.88);
      pmExtLegend->AddEntry(vn_pm,"#pi^{-}, 0 < y_{CM} < 0.5");
      pmExtLegend->AddEntry(vn_pm_ext,"#pi^{-}, 0.5 < y_{CM} < 1.0");
      pmExtLegend->SetFillColorAlpha(0,0);
      pmExtLegend->SetLineColorAlpha(0,0);

      TLegend *kpExtLegend = new TLegend(0.55, 0.68, 0.85, 0.9);
      kpExtLegend->AddEntry(vn_kp,"K^{+}, 0 < y_{CM} < 0.5");
      kpExtLegend->AddEntry(vn_kp_ext,"K^{+}, 0.5 < y_{CM} < 1.0");
      kpExtLegend->SetFillColorAlpha(0,0);
      kpExtLegend->SetLineColorAlpha(0,0);

      TLegend *kmExtLegend = new TLegend(0.25, 0.13, 0.55, 0.28);
      kmExtLegend->AddEntry(vn_km,"K^{-}, 0 < y_{CM} < 0.5");
      kmExtLegend->AddEntry(vn_km_ext,"K^{-}, 0.5 < y_{CM} < 1.0");
      kmExtLegend->SetFillColorAlpha(0,0);
      kmExtLegend->SetLineColorAlpha(0,0);

      TLegend *prExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      prExtLegend->AddEntry(vn_pr,"Proton, 0 < y_{CM} < 0.5");
      prExtLegend->AddEntry(vn_pr_ext,"Proton, 0.5 < y_{CM} < 1.0");
      prExtLegend->SetFillColorAlpha(0,0);
      prExtLegend->SetLineColorAlpha(0,0);

      
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

      //TLegend *kpLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      TLegend *kpLegend = new TLegend(0.13, 0.75, 0.35, 0.9);
      kpLegend->AddEntry(h_vn_yCM_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(h_vn_yCM_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(h_vn_yCM_40to60_kp, "40 - 60%");

      TLegend *kmLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      //kmLegend->AddEntry(h_vn_yCM_00to10_km, "0 - 10%");
      kmLegend->AddEntry(h_vn_yCM_10to40_km, "10 - 40%");
      //kmLegend->AddEntry(h_vn_yCM_40to60_km, "40 - 60%");

      TLegend *prLegend = new TLegend(0.7, 0.75, 0.9, 0.9);
      prLegend->AddEntry(h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_yCM_40to60_pr, "40 - 60%");


      
      TPaveText *piText = new TPaveText(15, -0.004, 45, 0.008, "NB");
      piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      piText->AddText("0 < y_{CM} < 0.5 GeV");
      piText->AddText("0.18 < p_{T} < 1.633 GeV");
      piText->SetFillColorAlpha(0,0);
      piText->SetLineColorAlpha(0,0);

      TPaveText *kaText = new TPaveText(20, -0.1, 50, -0.07, "NB");
      kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kaText->AddText("0 < y_{CM} < 0.5 GeV");
      kaText->AddText("0.18 < p_{T} < 1.62 GeV");
      kaText->SetFillColorAlpha(0,0);
      kaText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(5, -0.07, 35, -0.05, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 < p_{T} < 2.0314 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);

      
      TPaveText *ppText = new TPaveText(-0.77, 0.008, 0.2, 0.028, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0.18 < p_{T} < 1.633 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);

      TPaveText *pmText = new TPaveText(-0.77, 0.008, 0.2, 0.028, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0.18 < p_{T} < 1.633 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);

      TPaveText *kpText = new TPaveText(-0.4, 0.03, 0.6, 0.075, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0.18 < p_{T} < 1.62 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);

      TPaveText *kmText = new TPaveText(-0.75, 0.0, 0.15, 0.018, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0.18 < p_{T} < 1.62 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);

      TPaveText *prText_y = new TPaveText(-0.75, 0.038, 0.15, 0.083, "NB");
      prText_y->AddText("Proton");
      prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_y->AddText("0.4 < p_{T} < 2.0314 GeV");
      prText_y->SetFillColorAlpha(0,0);
      prText_y->SetLineColorAlpha(0,0);

      
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
      sh_cent_pp->SetMaximum(0.01);
      sh_cent_pp->SetMinimum(-0.05);
      sh_cent_pm->SetMaximum(0.01);
      sh_cent_pm->SetMinimum(-0.05);
      sh_cent_pp->Draw("E1P");
      sh_cent_pm->Draw("E1P SAME");
      piCentralityStack->Draw("NOSTACK E1P SAME");
      zeroLine->Draw("SAME");
      piLegend->Draw();
      piText->Draw();
      canvas->SaveAs(jobID + "_piCentralityStack.png");
      canvas->Clear();

      kaCentralityStack->Draw();
      kaCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kaCentralityStack->GetXaxis()->SetNdivisions(210);
      kaCentralityStack->SetMaximum(0.02);
      kaCentralityStack->SetMinimum(-0.12);
      sh_cent_kp->SetMaximum(0.02);
      sh_cent_kp->SetMinimum(-0.12);
      sh_cent_km->SetMaximum(0.02);
      sh_cent_km->SetMinimum(-0.12);
      sh_cent_kp->Draw("E1P");
      sh_cent_km->Draw("E1P SAME");
      kaCentralityStack->Draw("NOSTACK E1P SAME");
      zeroLine->Draw("SAME");
      kaLegend->Draw();
      kaText->Draw();
      canvas->SaveAs(jobID + "_kaCentralityStack.png");
      canvas->Clear();

      vn_pr->SetTitle("");
      vn_pr->GetYaxis()->SetTitleOffset(1.7);
      vn_pr->GetXaxis()->SetNdivisions(210);
      vn_pr->Draw("E1P");
      vn_pr->SetMaximum(0.01);
      vn_pr->SetMinimum(-0.09);
      sh_cent_pr->SetMaximum(0.01);
      sh_cent_pr->SetMinimum(-0.09);
      sh_cent_pr->SetMaximum(0.01);
      sh_cent_pr->SetMinimum(-0.09);
      sh_cent_pr->Draw("E1P SAME");
      zeroLine->Draw("SAME");
      prText->Draw();
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

      vn_Tpc->SetMarkerStyle(20);
      vn_Tpc->SetMarkerSize(2);
      vn_Tpc->GetXaxis()->SetNdivisions(210);
      vn_Tpc->GetYaxis()->SetTitleOffset(1.7);
      //vn_Tpc->SetMaximum(0.002);
      vn_Tpc->Draw("E1P");
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_Tpc.png");
      canvas->Clear();


      ppExtCentralityStack->Draw();
      ppExtCentralityStack->GetXaxis()->SetNdivisions(210);
      ppExtCentralityStack->SetMaximum(0.01);
      ppExtCentralityStack->SetMinimum(-0.05);
      ppExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      ppExtLegend->Draw();
      //ppExtText->Draw();
      canvas->SaveAs(jobID + "_ppExtCentralityStack.png");
      canvas->Clear();

      pmExtCentralityStack->Draw();
      pmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      pmExtCentralityStack->SetMaximum(0.01);
      pmExtCentralityStack->SetMinimum(-0.05);
      pmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pmExtLegend->Draw();
      //pmExtText->Draw();
      canvas->SaveAs(jobID + "_pmExtCentralityStack.png");
      canvas->Clear();

      kpExtCentralityStack->Draw();
      kpExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kpExtCentralityStack->SetMaximum(0.02);
      kpExtCentralityStack->SetMinimum(-0.12);
      kpExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kpExtLegend->Draw();
      //kpExtText->Draw();
      canvas->SaveAs(jobID + "_kpExtCentralityStack.png");
      canvas->Clear();

      kmExtCentralityStack->Draw();
      kmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kmExtCentralityStack->SetMaximum(0.02);
      kmExtCentralityStack->SetMinimum(-0.12);
      kmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kmExtLegend->Draw();
      //kmExtText->Draw();
      canvas->SaveAs(jobID + "_kmExtCentralityStack.png");
      canvas->Clear();

      prExtCentralityStack->Draw();
      prExtCentralityStack->GetXaxis()->SetNdivisions(210);
      prExtCentralityStack->SetMaximum(0.03);
      prExtCentralityStack->SetMinimum(-0.09);
      prExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      prExtLegend->Draw();
      //prExtText->Draw();
      canvas->SaveAs(jobID + "_prExtCentralityStack.png");
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

      // Zoom in on last one
      etaRegionStack->Draw();
      etaRegionStack->GetYaxis()->SetTitleOffset(1.7);
      etaRegionStack->GetXaxis()->SetNdivisions(210);
      etaRegionStack->Draw();
      etaRegionStack->SetMaximum(0.1);
      etaRegionStack->SetMinimum(-0.1);
      etaRegionStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      etaLegend->Draw();
      canvas->SaveAs(jobID + "_etaRegionStack_zoom.png");
      canvas->Clear();


      ppRapidityStack->Draw();
      ppRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      ppRapidityStack->GetXaxis()->SetNdivisions(210);
      ppRapidityStack->Draw();
      ppRapidityStack->SetMaximum(0.03);
      ppRapidityStack->SetMinimum(-0.05);
      ppRapidityStack->Draw("NOSTACK E1P");
      sh_y_pp->Draw("E1P SAME");
      zeroLine_y->Draw("SAME");
      ppLegend->Draw();
      ppText->Draw();
      canvas->SaveAs(jobID + "_ppRapidityStack.png");
      canvas->Clear();

      pmRapidityStack->Draw();
      pmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      pmRapidityStack->GetXaxis()->SetNdivisions(210);
      pmRapidityStack->Draw();
      pmRapidityStack->SetMaximum(0.03);
      pmRapidityStack->SetMinimum(-0.04);
      pmRapidityStack->Draw("NOSTACK E1P");
      sh_y_pm->Draw("E1P SAME");
      zeroLine_y->Draw("SAME");
      pmLegend->Draw();
      pmText->Draw();
      canvas->SaveAs(jobID + "_pmRapidityStack.png");
      canvas->Clear();

      kpRapidityStack->Draw();
      kpRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kpRapidityStack->GetXaxis()->SetNdivisions(210);
      kpRapidityStack->Draw();
      kpRapidityStack->SetMaximum(0.08);
      kpRapidityStack->SetMinimum(-0.13);
      kpRapidityStack->Draw("NOSTACK E1P");
      sh_y_kp->Draw("E1P SAME");
      zeroLine_y->Draw("SAME");
      kpLegend->Draw();
      kpText->Draw();
      canvas->SaveAs(jobID + "_kpRapidityStack.png");
      canvas->Clear();

      kmRapidityStack->Draw();
      kmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kmRapidityStack->GetXaxis()->SetNdivisions(210);
      kmRapidityStack->Draw();
      kmRapidityStack->SetMaximum(0.02);
      kmRapidityStack->SetMinimum(-0.06);
      kmRapidityStack->Draw("NOSTACK E1P");
      sh_y_km->Draw("E1P SAME");
      zeroLine_y->Draw("SAME");
      kmLegend->Draw();
      kmText->Draw();
      canvas->SaveAs(jobID + "_kmRapidityStack.png");
      canvas->Clear();

      prRapidityStack->Draw();
      prRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      prRapidityStack->GetXaxis()->SetNdivisions(210);
      prRapidityStack->Draw();
      prRapidityStack->SetMaximum(0.1);
      prRapidityStack->SetMinimum(-0.08);
      prRapidityStack->Draw("NOSTACK E1P");
      sh_y_pr->Draw("E1P SAME");
      zeroLine_y->Draw("SAME");
      prLegend->Draw();
      prText_y->Draw();
      canvas->SaveAs(jobID + "_prRapidityStack.png");
      canvas->Clear();      

    }
  else if (order_n_str == "3")
    {
      ppRapidityStack->Add(h_vn_yCM_00to10_pp);
      //ppRapidityStack->Add(h_vn_yCM_00to10_pp_mirror);
      ppRapidityStack->Add(h_vn_yCM_10to40_pp);
      //ppRapidityStack->Add(h_vn_yCM_10to40_pp_mirror);
      ppRapidityStack->Add(h_vn_yCM_40to60_pp);
      //ppRapidityStack->Add(h_vn_yCM_40to60_pp_mirror);

      pmRapidityStack->Add(h_vn_yCM_00to10_pm);
      //pmRapidityStack->Add(h_vn_yCM_00to10_pm_mirror);
      pmRapidityStack->Add(h_vn_yCM_10to40_pm);
      //pmRapidityStack->Add(h_vn_yCM_10to40_pm_mirror);
      pmRapidityStack->Add(h_vn_yCM_40to60_pm);
      //pmRapidityStack->Add(h_vn_yCM_40to60_pm_mirror);

      kpRapidityStack->Add(h_vn_yCM_00to10_kp);
      //kpRapidityStack->Add(h_vn_yCM_00to10_kp_mirror);
      kpRapidityStack->Add(h_vn_yCM_10to40_kp);
      //kpRapidityStack->Add(h_vn_yCM_10to40_kp_mirror);
      kpRapidityStack->Add(h_vn_yCM_40to60_kp);
      //kpRapidityStack->Add(h_vn_yCM_40to60_kp_mirror);

      //kmRapidityStack->Add(h_vn_yCM_00to10_km);
      kmRapidityStack->Add(h_vn_yCM_10to40_km);
      //kmRapidityStack->Add(h_vn_yCM_10to40_km_mirror);
      //kmRapidityStack->Add(h_vn_yCM_40to60_km);

      prRapidityStack->Add(h_vn_yCM_00to10_pr);
      //prRapidityStack->Add(h_vn_yCM_00to10_pr_mirror);
      prRapidityStack->Add(h_vn_yCM_10to40_pr);
      //prRapidityStack->Add(h_vn_yCM_10to40_pr_mirror);
      prRapidityStack->Add(h_vn_yCM_40to60_pr);
      //prRapidityStack->Add(h_vn_yCM_40to60_pr_mirror);

      TLegend *piLegend = new TLegend(0.67, 0.71, 0.805, 0.87);
      piLegend->AddEntry(vn_pp,"#pi^{+}");
      piLegend->AddEntry(vn_pm,"#pi^{-}");
      piLegend->SetFillColorAlpha(0,0);
      piLegend->SetLineColorAlpha(0,0);

      TLegend *kaLegend = new TLegend(0.55, 0.72, 0.65, 0.87);
      kaLegend->AddEntry(vn_kp,"K^{+}");
      kaLegend->AddEntry(vn_km,"K^{-}");
      kaLegend->SetFillColorAlpha(0,0);
      kaLegend->SetLineColorAlpha(0,0);


      TLegend *ppExtLegend = new TLegend(0.4, 0.62, 0.7, 0.82);
      ppExtLegend->AddEntry(vn_pp,"#pi^{+}, 0 < y_{CM} < 0.5");
      ppExtLegend->AddEntry(vn_pp_ext,"#pi^{+}, 0.5 < y_{CM} < 1.0");
      ppExtLegend->SetFillColorAlpha(0,0);
      ppExtLegend->SetLineColorAlpha(0,0);

      TLegend *pmExtLegend = new TLegend(0.15, 0.67, 0.45, 0.9);
      pmExtLegend->AddEntry(vn_pm,"#pi^{-}, 0 < y_{CM} < 0.5");
      pmExtLegend->AddEntry(vn_pm_ext,"#pi^{-}, 0.5 < y_{CM} < 1.0");
      pmExtLegend->SetFillColorAlpha(0,0);
      pmExtLegend->SetLineColorAlpha(0,0);

      TLegend *kpExtLegend = new TLegend(0.55, 0.7, 0.85, 0.9);
      kpExtLegend->AddEntry(vn_kp,"K^{+}, 0 < y_{CM} < 0.5");
      kpExtLegend->AddEntry(vn_kp_ext,"K^{+}, 0.5 < y_{CM} < 1.0");
      kpExtLegend->SetFillColorAlpha(0,0);
      kpExtLegend->SetLineColorAlpha(0,0);

      TLegend *kmExtLegend = new TLegend(0.28, 0.68, 0.55, 0.85);
      kmExtLegend->AddEntry(vn_km,"K^{-}, 0 < y_{CM} < 0.5");
      kmExtLegend->AddEntry(vn_km_ext,"K^{-}, 0.5 < y_{CM} < 1.0");
      kmExtLegend->SetFillColorAlpha(0,0);
      kmExtLegend->SetLineColorAlpha(0,0);

      TLegend *prExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      prExtLegend->AddEntry(vn_pr,"Proton, 0 < y_{CM} < 0.5");
      prExtLegend->AddEntry(vn_pr_ext,"Proton, 0.5 < y_{CM} < 1.0");
      prExtLegend->SetFillColorAlpha(0,0);
      prExtLegend->SetLineColorAlpha(0,0);

      

      TLegend *etaLegend = new TLegend(0.65, 0.25, 0.9, 0.45);
      etaLegend->AddEntry(vn_EpdE, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(vn_EpdF, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(vn_TpcB, "TPC -1.0 < #eta < 0");


      TLegend *ppLegend = new TLegend(0.16, 0.75, 0.36, 0.9);
      ppLegend->AddEntry(h_vn_yCM_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(h_vn_yCM_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(h_vn_yCM_40to60_pp, "40 - 60%");

      TLegend *pmLegend = new TLegend(0.16, 0.75, 0.36, 0.9);
      pmLegend->AddEntry(h_vn_yCM_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(h_vn_yCM_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(h_vn_yCM_40to60_pm, "40 - 60%");

      TLegend *kpLegend = new TLegend(0.16, 0.75, 0.36, 0.9);
      kpLegend->AddEntry(h_vn_yCM_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(h_vn_yCM_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(h_vn_yCM_40to60_kp, "40 - 60%");

      TLegend *kmLegend = new TLegend(0.16, 0.8, 0.36, 0.9);
      //kmLegend->AddEntry(h_vn_yCM_00to10_km, "0 - 10%");
      kmLegend->AddEntry(h_vn_yCM_10to40_km, "10 - 40%");
      //kmLegend->AddEntry(h_vn_yCM_40to60_km, "40 - 60%");

      TLegend *prLegend = new TLegend(0.16, 0.75, 0.36, 0.9);
      prLegend->AddEntry(h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_yCM_40to60_pr, "40 - 60%");


      TPaveText *piText = new TPaveText(12, 0.012, 38, 0.0185, "NB");
      piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      piText->AddText("0 < y_{CM} < 0.5 GeV");
      piText->AddText("0.18 < p_{T} < 1.633 GeV");
      piText->SetFillColorAlpha(0,0);
      piText->SetLineColorAlpha(0,0);

      TPaveText *kaText = new TPaveText(1, 0.1, 32, 0.185, "NB");
      kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kaText->AddText("0 < y_{CM} < 0.5 GeV");
      kaText->AddText("0.18 < p_{T} < 1.62 GeV");
      kaText->SetFillColorAlpha(0,0);
      kaText->SetLineColorAlpha(0,0);

      TPaveText *prText = new TPaveText(10, -0.015, 40, -0.025, "NB");
      prText->AddText("Proton");
      prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText->AddText("0 < y_{CM} < 0.5 GeV");
      prText->AddText("0.4 < p_{T} < 2.0314 GeV");
      prText->SetFillColorAlpha(0,0);
      prText->SetLineColorAlpha(0,0);

      TPaveText *ppText = new TPaveText(-0.4, 0.03, 0.6, 0.04, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0.18 < p_{T} < 1.633 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);

      TPaveText *pmText = new TPaveText(-0.4, 0.022, 0.6, 0.035, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0.18 < p_{T} < 1.633 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);

      TPaveText *kpText = new TPaveText(-0.4, 0.13, 0.6, 0.19, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0.18 < p_{T} < 1.62 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);

      TPaveText *kmText = new TPaveText(-0.4, 0.12, 0.6, 0.18, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0.18 < p_{T} < 1.62 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);

      TPaveText *prText_y = new TPaveText(-0.4, 0.005, 0.6, 0.02, "NB");
      prText_y->AddText("Proton");
      prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_y->AddText("0.4 < p_{T} < 2.0314 GeV");
      prText_y->SetFillColorAlpha(0,0);
      prText_y->SetLineColorAlpha(0,0);

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
      piCentralityStack->SetMaximum(0.02);
      piCentralityStack->SetMinimum(-0.015);
      piCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      piLegend->Draw();
      piText->Draw();
      canvas->SaveAs(jobID + "_piCentralityStack.png");
      canvas->Clear();

      kaCentralityStack->Draw();
      kaCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kaCentralityStack->GetXaxis()->SetNdivisions(210);
      kaCentralityStack->SetMaximum(0.2);
      kaCentralityStack->SetMinimum(-0.15);
      kaCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kaLegend->Draw();
      kaText->Draw();
      canvas->SaveAs(jobID + "_kaCentralityStack.png");
      canvas->Clear();

      vn_pr->SetTitle("");
      vn_pr->GetYaxis()->SetTitleOffset(1.7);
      vn_pr->GetXaxis()->SetNdivisions(210);
      vn_pr->Draw("E1P");
      vn_pr->SetMaximum(0.01);
      vn_pr->SetMinimum(-0.03);
      zeroLine->Draw("SAME");
      prText->Draw();
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

      vn_Tpc->SetMarkerStyle(20);
      vn_Tpc->SetMarkerSize(2);
      vn_Tpc->GetXaxis()->SetNdivisions(210);
      vn_Tpc->GetYaxis()->SetTitleOffset(1.7);
      vn_Tpc->SetMaximum(0.002);
      vn_Tpc->Draw("E1P");
      zeroLine->Draw("SAME");
      canvas->SaveAs(jobID + "_vn_Tpc.png");
      canvas->Clear();


      ppExtCentralityStack->Draw();
      ppExtCentralityStack->GetYaxis()->SetTitleOffset(1.8);
      ppExtCentralityStack->GetXaxis()->SetNdivisions(210);
      ppExtCentralityStack->SetMaximum(0.03);
      ppExtCentralityStack->SetMinimum(-0.015);
      ppExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      ppExtLegend->Draw();
      //ppExtText->Draw();
      canvas->SaveAs(jobID + "_ppExtCentralityStack.png");
      canvas->Clear();

      pmExtCentralityStack->Draw();
      pmExtCentralityStack->GetYaxis()->SetTitleOffset(1.8);
      pmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      pmExtCentralityStack->SetMaximum(0.02);
      pmExtCentralityStack->SetMinimum(-0.015);
      pmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      pmExtLegend->Draw();
      //pmExtText->Draw();
      canvas->SaveAs(jobID + "_pmExtCentralityStack.png");
      canvas->Clear();

      kpExtCentralityStack->Draw();
      kpExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kpExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kpExtCentralityStack->SetMaximum(0.2);
      kpExtCentralityStack->SetMinimum(-0.2);
      kpExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kpExtLegend->Draw();
      //kpExtText->Draw();
      canvas->SaveAs(jobID + "_kpExtCentralityStack.png");
      canvas->Clear();

      kmExtCentralityStack->Draw();
      kmExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      kmExtCentralityStack->GetXaxis()->SetNdivisions(210);
      kmExtCentralityStack->SetMaximum(0.4);
      kmExtCentralityStack->SetMinimum(-0.2);
      kmExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      kmExtLegend->Draw();
      //kmExtText->Draw();
      canvas->SaveAs(jobID + "_kmExtCentralityStack.png");
      canvas->Clear();

      prExtCentralityStack->Draw();
      prExtCentralityStack->GetYaxis()->SetTitleOffset(1.7);
      prExtCentralityStack->GetXaxis()->SetNdivisions(210);
      prExtCentralityStack->SetMaximum(0.01);
      prExtCentralityStack->SetMinimum(-0.08);
      prExtCentralityStack->Draw("NOSTACK E1P");
      zeroLine->Draw("SAME");
      prExtLegend->Draw();
      //prExtText->Draw();
      canvas->SaveAs(jobID + "_prExtCentralityStack.png");
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
      ppRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      ppRapidityStack->GetXaxis()->SetNdivisions(210);
      ppRapidityStack->Draw();
      ppRapidityStack->SetMaximum(0.045);
      ppRapidityStack->SetMinimum(-0.01);
      ppRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      ppLegend->Draw();
      ppText->Draw();
      canvas->SaveAs(jobID + "_ppRapidityStack.png");
      canvas->Clear();

      pmRapidityStack->Draw();
      pmRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      pmRapidityStack->GetXaxis()->SetNdivisions(210);
      pmRapidityStack->Draw();
      pmRapidityStack->SetMaximum(0.04);
      pmRapidityStack->SetMinimum(-0.02);
      pmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      pmLegend->Draw();
      pmText->Draw();
      canvas->SaveAs(jobID + "_pmRapidityStack.png");
      canvas->Clear();

      kpRapidityStack->Draw();
      kpRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kpRapidityStack->GetXaxis()->SetNdivisions(210);
      kpRapidityStack->Draw();
      kpRapidityStack->SetMaximum(0.2);
      kpRapidityStack->SetMinimum(-0.12);
      kpRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kpLegend->Draw();
      kpText->Draw();
      canvas->SaveAs(jobID + "_kpRapidityStack.png");
      canvas->Clear();

      kmRapidityStack->Draw();
      kmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kmRapidityStack->GetXaxis()->SetNdivisions(210);
      kmRapidityStack->Draw();
      kmRapidityStack->SetMaximum(0.2);
      kmRapidityStack->SetMinimum(-0.1);
      kmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kmLegend->Draw();
      kmText->Draw();
      canvas->SaveAs(jobID + "_kmRapidityStack.png");
      canvas->Clear();

      prRapidityStack->Draw();
      prRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      prRapidityStack->GetXaxis()->SetNdivisions(210);
      prRapidityStack->Draw();
      prRapidityStack->SetMaximum(0.025);
      prRapidityStack->SetMinimum(-0.055);
      prRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      prLegend->Draw();
      prText_y->Draw();
      canvas->SaveAs(jobID + "_prRapidityStack.png");
      canvas->Clear();
    }

  //resolutionInfo_INPUT->Close();
  file->Close();
}
