#include "PlotUtils.h"

void vnVsY(TString jobID, TString order_n_str)
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
  gStyle->SetEndErrorSize(6);

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


  TProfile2D *p2_vn_yCM_cent_pp = (TProfile2D*)file->Get("p2_vn_yCM_cent_pp");
  TProfile2D *p2_vn_yCM_cent_pm = (TProfile2D*)file->Get("p2_vn_yCM_cent_pm");
  TProfile2D *p2_vn_yCM_cent_kp = (TProfile2D*)file->Get("p2_vn_yCM_cent_kp");
  TProfile2D *p2_vn_yCM_cent_km = (TProfile2D*)file->Get("p2_vn_yCM_cent_km");
  TProfile2D *p2_vn_yCM_cent_pr = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr");
  TProfile2D *p2_vn_yCM_cent_pr_alt = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr_alt");
  TProfile2D *p2_vn_yCM_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr_symmetry");
  TProfile2D *p2_vn_yCM_cent_de = (TProfile2D*)file->Get("p2_vn_yCM_cent_de");
  TProfile2D *p2_vn_yCM_cent_tr = (TProfile2D*)file->Get("p2_vn_yCM_cent_tr");

  p2_vn_yCM_cent_kp->RebinY();
  p2_vn_yCM_cent_km->RebinY();


  TProfile *p_vn_yCM_00to10_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_00to10_pp", 15, 16);
  TProfile *p_vn_yCM_10to40_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_10to40_pp", 9, 14);
  TProfile *p_vn_yCM_40to60_pp = p2_vn_yCM_cent_pp->ProfileY("p_vn_yCM_40to60_pp", 5, 8);
  
  TProfile *p_vn_yCM_00to10_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_00to10_pm", 15, 16);
  TProfile *p_vn_yCM_10to40_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_10to40_pm", 9, 14);
  TProfile *p_vn_yCM_40to60_pm = p2_vn_yCM_cent_pm->ProfileY("p_vn_yCM_40to60_pm", 5, 8);

  TProfile *p_vn_yCM_00to10_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_00to10_kp", 15, 16);
  TProfile *p_vn_yCM_10to40_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_10to40_kp", 9, 14);
  TProfile *p_vn_yCM_40to60_kp = p2_vn_yCM_cent_kp->ProfileY("p_vn_yCM_40to60_kp", 5, 8);

  TProfile *p_vn_yCM_00to10_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_00to10_km", 15, 16);
  TProfile *p_vn_yCM_10to40_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_10to40_km", 9, 14);
  TProfile *p_vn_yCM_40to60_km = p2_vn_yCM_cent_km->ProfileY("p_vn_yCM_40to60_km", 5, 8);

  TProfile *p_vn_yCM_00to10_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_00to10_pr", 15, 16);
  TProfile *p_vn_yCM_10to40_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_10to40_pr", 9, 14);
  TProfile *p_vn_yCM_40to60_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_40to60_pr", 5, 8);

  TProfile *p_vn_yCM_00to10_pr_alt = p2_vn_yCM_cent_pr_alt->ProfileY("p_vn_yCM_00to10_pr_alt", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_alt = p2_vn_yCM_cent_pr_alt->ProfileY("p_vn_yCM_10to40_pr_alt", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_alt = p2_vn_yCM_cent_pr_alt->ProfileY("p_vn_yCM_40to60_pr_alt", 5, 8);

  TProfile *p_vn_yCM_00to10_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_00to10_pr_symm", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_10to40_pr_symm", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_40to60_pr_symm", 5, 8);

  TProfile *p_vn_yCM_00to10_de = p2_vn_yCM_cent_de->ProfileY("p_vn_yCM_00to10_de", 15, 16);
  TProfile *p_vn_yCM_10to40_de = p2_vn_yCM_cent_de->ProfileY("p_vn_yCM_10to40_de", 9, 14);
  TProfile *p_vn_yCM_40to60_de = p2_vn_yCM_cent_de->ProfileY("p_vn_yCM_40to60_de", 5, 8);

  TProfile *p_vn_yCM_00to10_tr = p2_vn_yCM_cent_tr->ProfileY("p_vn_yCM_00to10_tr", 15, 16);
  TProfile *p_vn_yCM_10to40_tr = p2_vn_yCM_cent_tr->ProfileY("p_vn_yCM_10to40_tr", 9, 14);
  TProfile *p_vn_yCM_40to60_tr = p2_vn_yCM_cent_tr->ProfileY("p_vn_yCM_40to60_tr", 5, 8);

  
  TH1D *h_vn_yCM_00to10_pp = new TH1D("h_vn_yCM_00to10_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pp = new TH1D("h_vn_yCM_10to40_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pp = new TH1D("h_vn_yCM_40to60_pp", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pm = new TH1D("h_vn_yCM_00to10_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pm = new TH1D("h_vn_yCM_10to40_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pm = new TH1D("h_vn_yCM_40to60_pm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_kp = new TH1D("h_vn_yCM_00to10_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_kp = new TH1D("h_vn_yCM_10to40_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_kp = new TH1D("h_vn_yCM_40to60_kp", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_km = new TH1D("h_vn_yCM_00to10_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_km = new TH1D("h_vn_yCM_10to40_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_km = new TH1D("h_vn_yCM_40to60_km", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_pr = new TH1D("h_vn_yCM_00to10_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr = new TH1D("h_vn_yCM_10to40_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr = new TH1D("h_vn_yCM_40to60_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pr_alt = new TH1D("h_vn_yCM_00to10_pr_alt", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr_alt = new TH1D("h_vn_yCM_10to40_pr_alt", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr_alt = new TH1D("h_vn_yCM_40to60_pr_alt", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pr_symm = new TH1D("h_vn_yCM_00to10_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr_symm = new TH1D("h_vn_yCM_10to40_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr_symm = new TH1D("h_vn_yCM_40to60_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_de = new TH1D("h_vn_yCM_00to10_de", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_de = new TH1D("h_vn_yCM_10to40_de", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_de = new TH1D("h_vn_yCM_40to60_de", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_tr = new TH1D("h_vn_yCM_00to10_tr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_tr = new TH1D("h_vn_yCM_10to40_tr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_tr = new TH1D("h_vn_yCM_40to60_tr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);


    //mirrored plots
  TH1D *h_vn_yCM_00to10_pp_mirror = new TH1D("h_vn_yCM_00to10_pp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pp_mirror = new TH1D("h_vn_yCM_10to40_pp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pp_mirror = new TH1D("h_vn_yCM_40to60_pp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_pm_mirror = new TH1D("h_vn_yCM_00to10_pm_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pm_mirror = new TH1D("h_vn_yCM_10to40_pm_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pm_mirror = new TH1D("h_vn_yCM_40to60_pm_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  TH1D *h_vn_yCM_00to10_kp_mirror = new TH1D("h_vn_yCM_00to10_kp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_kp_mirror = new TH1D("h_vn_yCM_10to40_kp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_kp_mirror = new TH1D("h_vn_yCM_40to60_kp_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_km_mirror = new TH1D("h_vn_yCM_00to10_km_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_10to40_km_mirror = new TH1D("h_vn_yCM_10to40_km_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);
  TH1D *h_vn_yCM_40to60_km_mirror = new TH1D("h_vn_yCM_40to60_km_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 10, -1, 1);

  TH1D *h_vn_yCM_00to10_pr_mirror = new TH1D("h_vn_yCM_00to10_pr_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_10to40_pr_mirror = new TH1D("h_vn_yCM_10to40_pr_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  TH1D *h_vn_yCM_40to60_pr_mirror = new TH1D("h_vn_yCM_40to60_pr_mirror", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);



  // Convert profiles to histograms
  h_vn_yCM_00to10_kp = p_vn_yCM_00to10_kp->ProjectionX();
  h_vn_yCM_10to40_kp = p_vn_yCM_10to40_kp->ProjectionX();
  h_vn_yCM_40to60_kp = p_vn_yCM_40to60_kp->ProjectionX();
  h_vn_yCM_00to10_km = p_vn_yCM_00to10_km->ProjectionX();
  h_vn_yCM_10to40_km = p_vn_yCM_10to40_km->ProjectionX();
  h_vn_yCM_40to60_km = p_vn_yCM_40to60_km->ProjectionX();
  
  h_vn_yCM_00to10_pp = p_vn_yCM_00to10_pp->ProjectionX();
  h_vn_yCM_10to40_pp = p_vn_yCM_10to40_pp->ProjectionX();
  h_vn_yCM_40to60_pp = p_vn_yCM_40to60_pp->ProjectionX();

  h_vn_yCM_00to10_pm = p_vn_yCM_00to10_pm->ProjectionX();
  h_vn_yCM_10to40_pm = p_vn_yCM_10to40_pm->ProjectionX();
  h_vn_yCM_40to60_pm = p_vn_yCM_40to60_pm->ProjectionX();

  h_vn_yCM_00to10_pr = p_vn_yCM_00to10_pr->ProjectionX();
  h_vn_yCM_10to40_pr = p_vn_yCM_10to40_pr->ProjectionX();
  h_vn_yCM_40to60_pr = p_vn_yCM_40to60_pr->ProjectionX();

  h_vn_yCM_00to10_pr_alt = p_vn_yCM_00to10_pr_alt->ProjectionX();
  h_vn_yCM_10to40_pr_alt = p_vn_yCM_10to40_pr_alt->ProjectionX();
  h_vn_yCM_40to60_pr_alt = p_vn_yCM_40to60_pr_alt->ProjectionX();

  h_vn_yCM_00to10_pr_symm = p_vn_yCM_00to10_pr_symm->ProjectionX();
  h_vn_yCM_10to40_pr_symm = p_vn_yCM_10to40_pr_symm->ProjectionX();
  h_vn_yCM_40to60_pr_symm = p_vn_yCM_40to60_pr_symm->ProjectionX();

  h_vn_yCM_00to10_de = p_vn_yCM_00to10_de->ProjectionX();
  h_vn_yCM_10to40_de = p_vn_yCM_10to40_de->ProjectionX();
  h_vn_yCM_40to60_de = p_vn_yCM_40to60_de->ProjectionX();

  h_vn_yCM_00to10_tr = p_vn_yCM_00to10_tr->ProjectionX();
  h_vn_yCM_10to40_tr = p_vn_yCM_10to40_tr->ProjectionX();
  h_vn_yCM_40to60_tr = p_vn_yCM_40to60_tr->ProjectionX();


  // Make mirrored plots
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

      if (j!=0)
	{
	  h_vn_yCM_00to10_kp_mirror->SetBinContent(j, h_vn_yCM_00to10_kp->GetBinContent(i));
	  h_vn_yCM_00to10_kp_mirror->SetBinError(j, h_vn_yCM_00to10_kp->GetBinError(i));
	  h_vn_yCM_10to40_kp_mirror->SetBinContent(j, h_vn_yCM_10to40_kp->GetBinContent(i));
	  h_vn_yCM_10to40_kp_mirror->SetBinError(j, h_vn_yCM_10to40_kp->GetBinError(i));
	  h_vn_yCM_40to60_kp_mirror->SetBinContent(j, h_vn_yCM_40to60_kp->GetBinContent(i));
	  h_vn_yCM_40to60_kp_mirror->SetBinError(j, h_vn_yCM_40to60_kp->GetBinError(i));

	  h_vn_yCM_00to10_km_mirror->SetBinContent(j, h_vn_yCM_00to10_km->GetBinContent(i));
	  h_vn_yCM_00to10_km_mirror->SetBinError(j, h_vn_yCM_00to10_km->GetBinError(i));
	  h_vn_yCM_10to40_km_mirror->SetBinContent(j, h_vn_yCM_10to40_km->GetBinContent(i));
	  h_vn_yCM_10to40_km_mirror->SetBinError(j, h_vn_yCM_10to40_km->GetBinError(i));
	  h_vn_yCM_40to60_km_mirror->SetBinContent(j, h_vn_yCM_40to60_km->GetBinContent(i));
	  h_vn_yCM_40to60_km_mirror->SetBinError(j, h_vn_yCM_40to60_km->GetBinError(i));
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
      
      if(j != 0)
	{
	  h_vn_yCM_00to10_pp_mirror->SetBinContent(j, h_vn_yCM_00to10_pp->GetBinContent(i));
	  h_vn_yCM_00to10_pp_mirror->SetBinError(j, h_vn_yCM_00to10_pp->GetBinError(i));
	  h_vn_yCM_10to40_pp_mirror->SetBinContent(j, h_vn_yCM_10to40_pp->GetBinContent(i));
	  h_vn_yCM_10to40_pp_mirror->SetBinError(j, h_vn_yCM_10to40_pp->GetBinError(i));
	  h_vn_yCM_40to60_pp_mirror->SetBinContent(j, h_vn_yCM_40to60_pp->GetBinContent(i));
	  h_vn_yCM_40to60_pp_mirror->SetBinError(j, h_vn_yCM_40to60_pp->GetBinError(i));

	  h_vn_yCM_00to10_pm_mirror->SetBinContent(j, h_vn_yCM_00to10_pm->GetBinContent(i));
	  h_vn_yCM_00to10_pm_mirror->SetBinError(j, h_vn_yCM_00to10_pm->GetBinError(i));
	  h_vn_yCM_10to40_pm_mirror->SetBinContent(j, h_vn_yCM_10to40_pm->GetBinContent(i));
	  h_vn_yCM_10to40_pm_mirror->SetBinError(j, h_vn_yCM_10to40_pm->GetBinError(i));
	  h_vn_yCM_40to60_pm_mirror->SetBinContent(j, h_vn_yCM_40to60_pm->GetBinContent(i));
	  h_vn_yCM_40to60_pm_mirror->SetBinError(j, h_vn_yCM_40to60_pm->GetBinError(i));

	  h_vn_yCM_00to10_pr_mirror->SetBinContent(j, h_vn_yCM_00to10_pr->GetBinContent(i));
	  h_vn_yCM_00to10_pr_mirror->SetBinError(j, h_vn_yCM_00to10_pr->GetBinError(i));
	  h_vn_yCM_10to40_pr_mirror->SetBinContent(j, h_vn_yCM_10to40_pr->GetBinContent(i));
	  h_vn_yCM_10to40_pr_mirror->SetBinError(j, h_vn_yCM_10to40_pr->GetBinError(i));
	  h_vn_yCM_40to60_pr_mirror->SetBinContent(j, h_vn_yCM_40to60_pr->GetBinContent(i));
	  h_vn_yCM_40to60_pr_mirror->SetBinError(j, h_vn_yCM_40to60_pr->GetBinError(i));
	}
    }


  THStack *ppRapidityStack   = new THStack("ppRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *pmRapidityStack   = new THStack("pmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *kpRapidityStack   = new THStack("kpRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *kmRapidityStack   = new THStack("kmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *prRapidityStack   = new THStack("prRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *prRapidityStack_alt = new THStack("prRapidityStack_alt", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *prRapidityStack_symm = new THStack("prRapidityStack_symm", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *deRapidityStack   = new THStack("deRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  THStack *trRapidityStack   = new THStack("trRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");


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


  h_vn_yCM_00to10_pp->SetMarkerStyle(20);
  h_vn_yCM_10to40_pp->SetMarkerStyle(20);
  h_vn_yCM_40to60_pp->SetMarkerStyle(20);
  h_vn_yCM_00to10_pp->SetMarkerColor(2);
  h_vn_yCM_10to40_pp->SetMarkerColor(4);
  h_vn_yCM_40to60_pp->SetMarkerColor(8);
  h_vn_yCM_00to10_pp->SetMarkerSize(2);
  h_vn_yCM_10to40_pp->SetMarkerSize(2);
  h_vn_yCM_40to60_pp->SetMarkerSize(2);
  h_vn_yCM_00to10_pp->SetLineColor(2);
  h_vn_yCM_10to40_pp->SetLineColor(4);
  h_vn_yCM_40to60_pp->SetLineColor(8);

  h_vn_yCM_00to10_pm->SetMarkerStyle(20);
  h_vn_yCM_10to40_pm->SetMarkerStyle(20);
  h_vn_yCM_40to60_pm->SetMarkerStyle(20);
  h_vn_yCM_00to10_pm->SetMarkerColor(2);
  h_vn_yCM_10to40_pm->SetMarkerColor(4);
  h_vn_yCM_40to60_pm->SetMarkerColor(8);
  h_vn_yCM_00to10_pm->SetMarkerSize(2);
  h_vn_yCM_10to40_pm->SetMarkerSize(2);
  h_vn_yCM_40to60_pm->SetMarkerSize(2);
  h_vn_yCM_00to10_pm->SetLineColor(2);
  h_vn_yCM_10to40_pm->SetLineColor(4);
  h_vn_yCM_40to60_pm->SetLineColor(8);

  h_vn_yCM_00to10_kp->SetMarkerStyle(20);
  h_vn_yCM_10to40_kp->SetMarkerStyle(20);
  h_vn_yCM_40to60_kp->SetMarkerStyle(20);
  h_vn_yCM_00to10_kp->SetMarkerColor(2);
  h_vn_yCM_10to40_kp->SetMarkerColor(4);
  h_vn_yCM_40to60_kp->SetMarkerColor(8);
  h_vn_yCM_00to10_kp->SetMarkerSize(2);
  h_vn_yCM_10to40_kp->SetMarkerSize(2);
  h_vn_yCM_40to60_kp->SetMarkerSize(2);
  h_vn_yCM_00to10_kp->SetLineColor(2);
  h_vn_yCM_10to40_kp->SetLineColor(4);
  h_vn_yCM_40to60_kp->SetLineColor(8);

  h_vn_yCM_00to10_km->SetMarkerStyle(20);
  h_vn_yCM_10to40_km->SetMarkerStyle(20);
  h_vn_yCM_40to60_km->SetMarkerStyle(20);
  h_vn_yCM_00to10_km->SetMarkerColor(2);
  h_vn_yCM_10to40_km->SetMarkerColor(4);
  h_vn_yCM_40to60_km->SetMarkerColor(8);
  h_vn_yCM_00to10_km->SetMarkerSize(2);
  h_vn_yCM_10to40_km->SetMarkerSize(2);
  h_vn_yCM_40to60_km->SetMarkerSize(2);
  h_vn_yCM_00to10_km->SetLineColor(2);
  h_vn_yCM_10to40_km->SetLineColor(4);
  h_vn_yCM_40to60_km->SetLineColor(8);

  h_vn_yCM_00to10_pr->SetMarkerStyle(20);
  h_vn_yCM_10to40_pr->SetMarkerStyle(20);
  h_vn_yCM_40to60_pr->SetMarkerStyle(20);
  h_vn_yCM_00to10_pr->SetMarkerColor(2);
  h_vn_yCM_10to40_pr->SetMarkerColor(4);
  h_vn_yCM_40to60_pr->SetMarkerColor(8);
  h_vn_yCM_00to10_pr->SetMarkerSize(2);
  h_vn_yCM_10to40_pr->SetMarkerSize(2);
  h_vn_yCM_40to60_pr->SetMarkerSize(2);
  h_vn_yCM_00to10_pr->SetLineColor(2);
  h_vn_yCM_10to40_pr->SetLineColor(4);
  h_vn_yCM_40to60_pr->SetLineColor(8);

  h_vn_yCM_00to10_pr_alt->SetMarkerStyle(20);
  h_vn_yCM_10to40_pr_alt->SetMarkerStyle(20);
  h_vn_yCM_40to60_pr_alt->SetMarkerStyle(20);
  h_vn_yCM_00to10_pr_alt->SetMarkerColor(2);
  h_vn_yCM_10to40_pr_alt->SetMarkerColor(4);
  h_vn_yCM_40to60_pr_alt->SetMarkerColor(8);
  h_vn_yCM_00to10_pr_alt->SetMarkerSize(2);
  h_vn_yCM_10to40_pr_alt->SetMarkerSize(2);
  h_vn_yCM_40to60_pr_alt->SetMarkerSize(2);
  h_vn_yCM_00to10_pr_alt->SetLineColor(2);
  h_vn_yCM_10to40_pr_alt->SetLineColor(4);
  h_vn_yCM_40to60_pr_alt->SetLineColor(8);

  h_vn_yCM_00to10_pr_symm->SetMarkerStyle(20);
  h_vn_yCM_10to40_pr_symm->SetMarkerStyle(20);
  h_vn_yCM_40to60_pr_symm->SetMarkerStyle(20);
  h_vn_yCM_00to10_pr_symm->SetMarkerColor(2);
  h_vn_yCM_10to40_pr_symm->SetMarkerColor(4);
  h_vn_yCM_40to60_pr_symm->SetMarkerColor(8);
  h_vn_yCM_00to10_pr_symm->SetMarkerSize(2);
  h_vn_yCM_10to40_pr_symm->SetMarkerSize(2);
  h_vn_yCM_40to60_pr_symm->SetMarkerSize(2);
  h_vn_yCM_00to10_pr_symm->SetLineColor(2);
  h_vn_yCM_10to40_pr_symm->SetLineColor(4);
  h_vn_yCM_40to60_pr_symm->SetLineColor(8);

  h_vn_yCM_00to10_de->SetMarkerStyle(20);
  h_vn_yCM_10to40_de->SetMarkerStyle(20);
  h_vn_yCM_40to60_de->SetMarkerStyle(20);
  h_vn_yCM_00to10_de->SetMarkerColor(2);
  h_vn_yCM_10to40_de->SetMarkerColor(4);
  h_vn_yCM_40to60_de->SetMarkerColor(8);
  h_vn_yCM_00to10_de->SetMarkerSize(2);
  h_vn_yCM_10to40_de->SetMarkerSize(2);
  h_vn_yCM_40to60_de->SetMarkerSize(2);
  h_vn_yCM_00to10_de->SetLineColor(2);
  h_vn_yCM_10to40_de->SetLineColor(4);
  h_vn_yCM_40to60_de->SetLineColor(8);

  h_vn_yCM_00to10_tr->SetMarkerStyle(20);
  h_vn_yCM_10to40_tr->SetMarkerStyle(20);
  h_vn_yCM_40to60_tr->SetMarkerStyle(20);
  h_vn_yCM_00to10_tr->SetMarkerColor(2);
  h_vn_yCM_10to40_tr->SetMarkerColor(4);
  h_vn_yCM_40to60_tr->SetMarkerColor(8);
  h_vn_yCM_00to10_tr->SetMarkerSize(2);
  h_vn_yCM_10to40_tr->SetMarkerSize(2);
  h_vn_yCM_40to60_tr->SetMarkerSize(2);
  h_vn_yCM_00to10_tr->SetLineColor(2);
  h_vn_yCM_10to40_tr->SetLineColor(4);
  h_vn_yCM_40to60_tr->SetLineColor(8);



  //mirrored plots
  h_vn_yCM_00to10_pp_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_pp_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_pp_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_pp_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_pp_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_pp_mirror->SetMarkerColor(8);
  h_vn_yCM_00to10_pp_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_pp_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_pp_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_pp_mirror->SetLineColor(2);
  h_vn_yCM_10to40_pp_mirror->SetLineColor(4);
  h_vn_yCM_40to60_pp_mirror->SetLineColor(8);

  h_vn_yCM_00to10_pm_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_pm_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_pm_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_pm_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_pm_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_pm_mirror->SetMarkerColor(8);
  h_vn_yCM_00to10_pm_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_pm_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_pm_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_pm_mirror->SetLineColor(2);
  h_vn_yCM_10to40_pm_mirror->SetLineColor(4);
  h_vn_yCM_40to60_pm_mirror->SetLineColor(8);

  h_vn_yCM_00to10_kp_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_kp_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_kp_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_kp_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_kp_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_kp_mirror->SetMarkerColor(8);
  h_vn_yCM_00to10_kp_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_kp_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_kp_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_kp_mirror->SetLineColor(2);
  h_vn_yCM_10to40_kp_mirror->SetLineColor(4);
  h_vn_yCM_40to60_kp_mirror->SetLineColor(8);

  h_vn_yCM_00to10_km_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_km_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_km_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_km_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_km_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_km_mirror->SetMarkerColor(8);
  h_vn_yCM_00to10_km_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_km_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_km_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_km_mirror->SetLineColor(2);
  h_vn_yCM_10to40_km_mirror->SetLineColor(4);
  h_vn_yCM_40to60_km_mirror->SetLineColor(8);

  h_vn_yCM_00to10_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_10to40_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_40to60_pr_mirror->SetMarkerStyle(24);
  h_vn_yCM_00to10_pr_mirror->SetMarkerColor(2);
  h_vn_yCM_10to40_pr_mirror->SetMarkerColor(4);
  h_vn_yCM_40to60_pr_mirror->SetMarkerColor(8);
  h_vn_yCM_00to10_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_10to40_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_40to60_pr_mirror->SetMarkerSize(2);
  h_vn_yCM_00to10_pr_mirror->SetLineColor(2);
  h_vn_yCM_10to40_pr_mirror->SetLineColor(4);
  h_vn_yCM_40to60_pr_mirror->SetLineColor(8);



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


      TPaveText *ppText = new TPaveText(-0.77, 0.008, 0.2, 0.028, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);

      TPaveText *pmText = new TPaveText(-0.77, 0.008, 0.2, 0.028, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);

      TPaveText *kpText = new TPaveText(-0.4, 0.03, 0.6, 0.075, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);

      TPaveText *kmText = new TPaveText(-0.75, 0.0, 0.15, 0.018, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);

      TPaveText *prText_y = new TPaveText(-0.75, 0.038, 0.15, 0.083, "NB");
      prText_y->AddText("Proton");
      prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_y->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
      prText_y->SetFillColorAlpha(0,0);
      prText_y->SetLineColorAlpha(0,0);


      TLine *zeroLine_y = new TLine(-1, 0, 1, 0);
      zeroLine_y->SetLineStyle(9);

      
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
      // Trim and clean up x-axis
      h_vn_yCM_00to10_kp = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_kp);
      h_vn_yCM_10to40_kp = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_kp);
      h_vn_yCM_40to60_kp = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_kp);
      h_vn_yCM_00to10_km = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_km);
      h_vn_yCM_10to40_km = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_km);
      h_vn_yCM_40to60_km = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_km);
  
      h_vn_yCM_00to10_pp = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_pp);
      h_vn_yCM_10to40_pp = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_pp);
      h_vn_yCM_40to60_pp = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_pp);

      h_vn_yCM_00to10_pm = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_pm);
      h_vn_yCM_10to40_pm = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_pm);
      h_vn_yCM_40to60_pm = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_pm);

      h_vn_yCM_00to10_pr_alt = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_pr_alt);
      h_vn_yCM_10to40_pr_alt = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_pr_alt);
      h_vn_yCM_40to60_pr_alt = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_pr_alt);

      h_vn_yCM_00to10_de = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_de);
      h_vn_yCM_10to40_de = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_de);
      h_vn_yCM_40to60_de = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_de);

      h_vn_yCM_00to10_tr = PlotUtils::trimRapidityPlot(h_vn_yCM_00to10_tr);
      h_vn_yCM_10to40_tr = PlotUtils::trimRapidityPlot(h_vn_yCM_10to40_tr);
      h_vn_yCM_40to60_tr = PlotUtils::trimRapidityPlot(h_vn_yCM_40to60_tr);

      
      ppRapidityStack->Add(h_vn_yCM_00to10_pp);
      ppRapidityStack->Add(h_vn_yCM_10to40_pp);
      ppRapidityStack->Add(h_vn_yCM_40to60_pp);

      pmRapidityStack->Add(h_vn_yCM_00to10_pm);
      pmRapidityStack->Add(h_vn_yCM_10to40_pm);
      pmRapidityStack->Add(h_vn_yCM_40to60_pm);

      kpRapidityStack->Add(h_vn_yCM_00to10_kp);
      kpRapidityStack->Add(h_vn_yCM_10to40_kp);
      kpRapidityStack->Add(h_vn_yCM_40to60_kp);

      kmRapidityStack->Add(h_vn_yCM_10to40_km);

      prRapidityStack->Add(h_vn_yCM_00to10_pr);
      prRapidityStack->Add(h_vn_yCM_10to40_pr);
      prRapidityStack->Add(h_vn_yCM_40to60_pr);

      prRapidityStack_alt->Add(h_vn_yCM_00to10_pr_alt);
      prRapidityStack_alt->Add(h_vn_yCM_10to40_pr_alt);
      prRapidityStack_alt->Add(h_vn_yCM_40to60_pr_alt);

      prRapidityStack_symm->Add(h_vn_yCM_00to10_pr_symm);
      prRapidityStack_symm->Add(h_vn_yCM_10to40_pr_symm);
      prRapidityStack_symm->Add(h_vn_yCM_40to60_pr_symm);

      deRapidityStack->Add(h_vn_yCM_00to10_de);
      deRapidityStack->Add(h_vn_yCM_10to40_de);
      deRapidityStack->Add(h_vn_yCM_40to60_de);

      trRapidityStack->Add(h_vn_yCM_00to10_tr);
      trRapidityStack->Add(h_vn_yCM_10to40_tr);
      trRapidityStack->Add(h_vn_yCM_40to60_tr);

      /*
      TFile *newFile = new TFile("v3_vs_yCM.root", "RECREATE");
      newFile->cd();

      h_vn_yCM_00to10_pr->Write();
      h_vn_yCM_10to40_pr->Write();
      h_vn_yCM_40to60_pr->Write();
      h_vn_yCM_00to10_pr_symm->Write();
      h_vn_yCM_10to40_pr_symm->Write();
      h_vn_yCM_40to60_pr_symm->Write();

      newFile->Close();
      */

      TLegend *ppLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      ppLegend->AddEntry(h_vn_yCM_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(h_vn_yCM_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(h_vn_yCM_40to60_pp, "40 - 60%");
      ppLegend->SetBorderSize(0);
      ppLegend->SetFillColorAlpha(0,0);

      TLegend *pmLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      pmLegend->AddEntry(h_vn_yCM_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(h_vn_yCM_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(h_vn_yCM_40to60_pm, "40 - 60%");
      pmLegend->SetBorderSize(0);
      pmLegend->SetFillColorAlpha(0,0);

      TLegend *kpLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      kpLegend->AddEntry(h_vn_yCM_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(h_vn_yCM_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(h_vn_yCM_40to60_kp, "40 - 60%");
      kpLegend->SetBorderSize(0);
      kpLegend->SetFillColorAlpha(0,0);

      TLegend *kmLegend = new TLegend(0.18, 0.77, 0.38, 0.87);
      //kmLegend->AddEntry(h_vn_yCM_00to10_km, "0 - 10%");
      kmLegend->AddEntry(h_vn_yCM_10to40_km, "10 - 40%");
      //kmLegend->AddEntry(h_vn_yCM_40to60_km, "40 - 60%");
      kmLegend->SetBorderSize(0);
      kmLegend->SetFillColorAlpha(0,0);
      /*
      TLegend *prLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      prLegend->AddEntry(h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_yCM_40to60_pr, "40 - 60%");
      prLegend->SetBorderSize(0);
      prLegend->SetFillColorAlpha(0,0);
      */
      TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      prLegend->AddEntry(h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(h_vn_yCM_40to60_pr, "40 - 60%");
      prLegend->SetBorderSize(0);
      prLegend->SetFillColorAlpha(0,0);

      TLegend *deLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      deLegend->AddEntry(h_vn_yCM_00to10_de, "0 - 10%");
      deLegend->AddEntry(h_vn_yCM_10to40_de, "10 - 40%");
      deLegend->AddEntry(h_vn_yCM_40to60_de, "40 - 60%");
      deLegend->SetBorderSize(0);
      deLegend->SetFillColorAlpha(0,0);

      TLegend *trLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      trLegend->AddEntry(h_vn_yCM_00to10_tr, "0 - 10%");
      trLegend->AddEntry(h_vn_yCM_10to40_tr, "10 - 40%");
      trLegend->AddEntry(h_vn_yCM_40to60_tr, "40 - 60%");
      trLegend->SetBorderSize(0);
      trLegend->SetFillColorAlpha(0,0);

      
      
      TPaveText *ppText = new TPaveText(0.5, 0.09, 0.8, 0.16, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);
      ppText->SetTextSize(.04);

      TPaveText *pmText = new TPaveText(0.5, 0.09, 0.8, 0.16, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);
      pmText->SetTextSize(.04);
 
      TPaveText *kpText = new TPaveText(0.5, 0.1, 0.8, 0.17, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);
      kpText->SetTextSize(.04);
 
      TPaveText *kmText = new TPaveText(0.3, 0.05, 0.7, 0.12, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0.18 #leq p_{T} #leq 1.6 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);
      kmText->SetTextSize(.04);
      
      TPaveText *prText_y = new TPaveText(-0.2, 0.02, 0.9, 0.05, "NB");
      prText_y->AddText("Proton");
      prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_y->AddText("0.4 #leq p_{T} #leq 2.0 GeV");
      prText_y->SetFillColorAlpha(0,0);
      prText_y->SetLineColorAlpha(0,0);
      prText_y->SetTextSize(.035);
 
      TPaveText *prText_y_symm = new TPaveText(-0.2, 0.02, 0.9, 0.05, "NB");
      prText_y_symm->AddText("Proton");
      prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_y_symm->AddText("1.0 #leq p_{T} #leq 2.5 GeV");
      prText_y_symm->SetFillColorAlpha(0,0);
      prText_y_symm->SetLineColorAlpha(0,0);
      prText_y_symm->SetTextSize(.035);

      TPaveText *prText_alt = new TPaveText(0.2, 0.02, 0.4, 0.1, "NB");
      prText_alt->AddText("p");
      prText_alt->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prText_alt->AddText("1.0 #leq p_{T} #leq 2.5 GeV");
      prText_alt->SetFillColorAlpha(0,0);
      prText_alt->SetLineColorAlpha(0,0);
      prText_alt->SetTextSize(.035);

      TPaveText *deText = new TPaveText(0.2, -0.17, 0.4, -0.1, "NB");
      deText->AddText("d");
      deText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      deText->AddText("1.0 #leq p_{T} #leq 2.5 GeV");
      deText->SetFillColorAlpha(0,0);
      deText->SetLineColorAlpha(0,0);
      deText->SetTextSize(.035);

      TPaveText *trText = new TPaveText(0.2, -0.17, 0.4, -0.1, "NB");
      trText->AddText("t");
      trText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      trText->AddText("1.0 #leq p_{T} #leq 2.5 GeV");
      trText->SetFillColorAlpha(0,0);
      trText->SetLineColorAlpha(0,0);
      trText->SetTextSize(.035);

      TLine *zeroLine_y = new TLine(0, 0, 1, 0);
      zeroLine_y->SetLineStyle(9);

      TLine *zeroLine_y_pr = new TLine(-1, 0, 1, 0);
      zeroLine_y_pr->SetLineStyle(9);

      Double_t rapidityUpperBound = 0.18;
      Double_t rapidityLowerBound = -0.15;
      Double_t rapidityUpperBound_pr = 0.06;
      Double_t rapidityLowerBound_pr = -0.1;
      Double_t rapidityUpperBound_pdt = 0.15;
      Double_t rapidityLowerBound_pdt = -0.2;
      
      ppRapidityStack->Draw();
      ppRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      ppRapidityStack->GetXaxis()->SetNdivisions(210);
      ppRapidityStack->Draw();
      ppRapidityStack->SetMaximum(rapidityUpperBound);
      ppRapidityStack->SetMinimum(rapidityLowerBound);
      ppRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      ppRapidityStack->Draw("NOSTACK E1P SAME");
      ppLegend->Draw();
      ppText->Draw();
      canvas->SaveAs(jobID + "_ppRapidityStack.png");
      canvas->Clear();

      pmRapidityStack->Draw();
      pmRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      pmRapidityStack->GetXaxis()->SetNdivisions(210);
      pmRapidityStack->Draw();
      pmRapidityStack->SetMaximum(rapidityUpperBound);
      pmRapidityStack->SetMinimum(rapidityLowerBound);
      pmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      pmRapidityStack->Draw("NOSTACK E1P SAME");
      pmLegend->Draw();
      pmText->Draw();
      canvas->SaveAs(jobID + "_pmRapidityStack.png");
      canvas->Clear();

      kpRapidityStack->Draw();
      kpRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kpRapidityStack->GetXaxis()->SetNdivisions(210);
      kpRapidityStack->Draw();
      kpRapidityStack->SetMaximum(rapidityUpperBound);
      kpRapidityStack->SetMinimum(rapidityLowerBound);
      kpRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kpRapidityStack->Draw("NOSTACK E1P SAME");
      kpLegend->Draw();
      kpText->Draw();
      canvas->SaveAs(jobID + "_kpRapidityStack.png");
      canvas->Clear();

      kmRapidityStack->Draw();
      kmRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      kmRapidityStack->GetXaxis()->SetNdivisions(210);
      kmRapidityStack->Draw();
      kmRapidityStack->SetMaximum(rapidityUpperBound);
      kmRapidityStack->SetMinimum(rapidityLowerBound);
      kmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      kmRapidityStack->Draw("NOSTACK E1P SAME");
      kmLegend->Draw();
      kmText->Draw();
      canvas->SaveAs(jobID + "_kmRapidityStack.png");
      canvas->Clear();

      prRapidityStack->Draw();
      prRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      prRapidityStack->GetXaxis()->SetNdivisions(210);
      prRapidityStack->Draw();
      prRapidityStack->SetMaximum(rapidityUpperBound_pr);
      prRapidityStack->SetMinimum(rapidityLowerBound_pr);
      prRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y_pr->Draw("SAME");
      prRapidityStack->Draw("NOSTACK E1P SAME");
      prLegend->Draw();
      prText_y->Draw();
      canvas->SaveAs(jobID + "_prRapidityStack.png");
      canvas->Clear();


      prRapidityStack_alt->Draw();
      prRapidityStack_alt->GetYaxis()->SetTitleOffset(1.7);
      prRapidityStack_alt->GetXaxis()->SetNdivisions(210);
      prRapidityStack_alt->Draw();
      prRapidityStack_alt->SetMaximum(rapidityUpperBound_pdt);
      prRapidityStack_alt->SetMinimum(rapidityLowerBound_pdt);
      prRapidityStack_alt->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      prRapidityStack_alt->Draw("NOSTACK E1P SAME");
      prLegend->Draw();
      prText_alt->Draw();
      canvas->SaveAs(jobID + "_prRapidityStack_alt.png");
      canvas->Clear();


      prRapidityStack_symm->Draw();
      prRapidityStack_symm->GetYaxis()->SetTitleOffset(1.7);
      prRapidityStack_symm->GetXaxis()->SetNdivisions(210);
      prRapidityStack_symm->Draw();
      prRapidityStack_symm->SetMaximum(rapidityUpperBound_pr);
      prRapidityStack_symm->SetMinimum(rapidityLowerBound_pr);
      prRapidityStack_symm->Draw("NOSTACK E1P");
      zeroLine_y_pr->Draw("SAME");
      prRapidityStack_symm->Draw("NOSTACK E1P SAME");
      prLegend->Draw();
      prText_y_symm->Draw();
      canvas->SaveAs(jobID + "_prRapidityStack_symm.png");
      canvas->Clear();

      deRapidityStack->Draw();
      deRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      deRapidityStack->GetXaxis()->SetNdivisions(210);
      deRapidityStack->Draw();
      deRapidityStack->SetMaximum(rapidityUpperBound_pdt);
      deRapidityStack->SetMinimum(rapidityLowerBound_pdt);
      deRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      deRapidityStack->Draw("NOSTACK E1P SAME");
      deLegend->Draw();
      deText->Draw();
      canvas->SaveAs(jobID + "_deRapidityStack.png");
      canvas->Clear();

      trRapidityStack->Draw();
      trRapidityStack->GetYaxis()->SetTitleOffset(1.7);
      trRapidityStack->GetXaxis()->SetNdivisions(210);
      trRapidityStack->Draw();
      trRapidityStack->SetMaximum(rapidityUpperBound_pdt);
      trRapidityStack->SetMinimum(rapidityLowerBound_pdt);
      trRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      trRapidityStack->Draw("NOSTACK E1P SAME");
      trLegend->Draw();
      trText->Draw();
      canvas->SaveAs(jobID + "_trRapidityStack.png");
      canvas->Clear();
    }

  file->Close();
}
