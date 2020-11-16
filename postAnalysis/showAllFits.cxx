void usage();

void showAllFits(TString jobID, TString fitType)
{
  if (jobID.EqualTo("")) { usage(); return; }
  if (!fitType.EqualTo("E") && !fitType.EqualTo("L") && !fitType.EqualTo("B")) { usage(); return; }
  
  TString intFileName = jobID + ".picoDst.result."+fitType+"yields.root";
  TFile *intFile = TFile::Open(intFileName, "read");
  if (!intFile) { std::cout << "File with integrated spectra not found!" << std::endl; return; }

  TH1D *pp_int = (TH1D*)intFile->Get("pp_yield");
  TH1D *pm_int = (TH1D*)intFile->Get("pm_yield");
  TH1D *kp_int = (TH1D*)intFile->Get("kp_yield");
  TH1D *km_int = (TH1D*)intFile->Get("km_yield");
  TH1D *pr_int = (TH1D*)intFile->Get("pr_yield");

  TObjArray *intArray = new TObjArray();
  intArray->Add(pp_int);
  intArray->Add(pm_int);
  intArray->Add(kp_int);
  intArray->Add(km_int);
  intArray->Add(pr_int);
  intArray->Compress();


  TCanvas *canvas = new TCanvas("canvas", "Title", 875, 675);
  canvas->SetGrid();
  canvas->SetLogy();
  canvas->SetTicks();

  gStyle->SetOptFit(1);

  for (Int_t j = 0; j < intArray->GetEntriesFast(); j++)  //Loop over particle types
    {
      TString prefix = ((TH1D*)intArray->At(j))->GetName();

      for (Int_t k = 1; k <= ((TH1D*)intArray->At(j))->GetNbinsX(); k++)  //Loop over bin numbers
	{
	  TString binStr;
	  binStr.Form("%d", k);

	  TString histName = prefix + "_bin" + binStr;
	  TString funcName = prefix + "_func" + binStr;

	  TH1D *hist = (TH1D*)intFile->Get(histName);
	  if (!hist) continue;

	  TF1 *func = (TF1*)intFile->Get(funcName);
	  func->SetLineColor(2);

	  hist->Draw("E1");
	  func->Draw("SAME");

	  canvas->Update();
	  canvas->SaveAs(histName+".png");
	  canvas->Clear();
	}
    }

  canvas->Close();
  intFile->Close();
  delete intArray;
  delete intFile;
}//End showAllFits()


void usage()
{
  std::cout << "Usage: root -l -b -q 'showAllFits.cxx(\"[job ID]\",\"[fit type]\")'" << std::endl
	    << "Fit types: E = exponential, L = Levy, B = Blast-wave" << std::endl;
}

