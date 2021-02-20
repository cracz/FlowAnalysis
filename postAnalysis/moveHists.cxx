void moveHists(TString fileName = "86E594D916F5AA37F98DE84D418CB64A.picoDst.result.combined.root")
{
  TFile *file = TFile::Open(fileName);

  TProfile *p_vn_yCM_00to10_pp = (TProfile*)file->Get("p_vn_yCM_00to10_pp");
  TProfile *p_vn_yCM_10to40_pp = (TProfile*)file->Get("p_vn_yCM_10to40_pp");
  TProfile *p_vn_yCM_40to60_pp = (TProfile*)file->Get("p_vn_yCM_40to60_pp");

  TProfile *p_vn_yCM_00to10_pm = (TProfile*)file->Get("p_vn_yCM_00to10_pm");
  TProfile *p_vn_yCM_10to40_pm = (TProfile*)file->Get("p_vn_yCM_10to40_pm");
  TProfile *p_vn_yCM_40to60_pm = (TProfile*)file->Get("p_vn_yCM_40to60_pm");
  
  TProfile *p_vn_yCM_00to10_kp = (TProfile*)file->Get("p_vn_yCM_00to10_kp");
  TProfile *p_vn_yCM_10to40_kp = (TProfile*)file->Get("p_vn_yCM_10to40_kp");
  TProfile *p_vn_yCM_40to60_kp = (TProfile*)file->Get("p_vn_yCM_40to60_kp");
  
  TProfile *p_vn_yCM_00to10_km = (TProfile*)file->Get("p_vn_yCM_00to10_km");
  TProfile *p_vn_yCM_10to40_km = (TProfile*)file->Get("p_vn_yCM_10to40_km");
  TProfile *p_vn_yCM_40to60_km = (TProfile*)file->Get("p_vn_yCM_40to60_km");
  
  TProfile *p_vn_yCM_00to10_pr = (TProfile*)file->Get("p_vn_yCM_00to10_pr");
  TProfile *p_vn_yCM_10to40_pr = (TProfile*)file->Get("p_vn_yCM_10to40_pr");
  TProfile *p_vn_yCM_40to60_pr = (TProfile*)file->Get("p_vn_yCM_40to60_pr");

  TFile *newFile = new TFile("v2_vs_yCM.root", "RECREATE");
  newFile->cd();

  p_vn_yCM_00to10_pp->Write();
  p_vn_yCM_10to40_pp->Write();
  p_vn_yCM_40to60_pp->Write();

  p_vn_yCM_00to10_pm->Write();
  p_vn_yCM_10to40_pm->Write();
  p_vn_yCM_40to60_pm->Write();
  
  p_vn_yCM_00to10_kp->Write();
  p_vn_yCM_10to40_kp->Write();
  p_vn_yCM_40to60_kp->Write();
  
  p_vn_yCM_00to10_km->Write();
  p_vn_yCM_10to40_km->Write();
  p_vn_yCM_40to60_km->Write();
  
  p_vn_yCM_00to10_pr->Write();
  p_vn_yCM_10to40_pr->Write();
  p_vn_yCM_40to60_pr->Write();


  file->Close();
  newFile->Close();
}
