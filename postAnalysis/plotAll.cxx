void plotAll(TString jobID)
{
  if (!jobID) { std::cout << "Supply a job ID!" << std::endl; return; }
  TString fileName = jobID + ".picoDst.result.combined.root";

  TFile *file = TFile::Open(fileName);
  if(!file) {cout << "Wrong file!" << endl; return;}

  ////////
  // PLOTTING EVERYTHING FROM THE FILE
  ////////

  TH1F  *hist;
  TString name;
  TIter next(file->GetListOfKeys());
  TKey  *key;

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 875, 675);


  while ((key = (TKey*)next()))
    {
      canvas->SetCanvasSize(875, 675);   // reset canvas size
      canvas->SetBottomMargin(0.1);
      
      hist = (TH1F*)key->ReadObj();
      name = hist->GetName();

      if (hist->GetEntries() == 0) { continue; }

      canvas->SetLogy(0);    //Reset y axis to linear
      canvas->SetGrid();     //Draw grid lines
      gStyle->SetOptStat(1); //Reset stat box
      //gStyle->SetPalette();  //Reset color palette

      Bool_t errors = false; //Draw error bars
      
      TClass *cl = gROOT->GetClass(key->GetClassName());

      if (cl->InheritsFrom("TH2"))                       //Keep this if/else order!! A TH2 still inherits from TH1.
	{
	  canvas->SetRightMargin(0.14);
	  canvas->SetLeftMargin(0.12);
	  hist->SetTitleOffset(1.2, "y");
	  gStyle->SetOptStat(11);
	  //gStyle->SetPalette(53); //kDarkBodyRadiator
	  //gStyle->SetPalette(57); //kBird Doesn't work
	  canvas->SetLogz();

	  if (name.Contains("psi") || name.Contains("Scan") || name.Contains("MvsY") ||
	      name.Contains("pp_vs_eta") || name.Contains("p2"))
	    {
	      canvas->SetLogz(0); // Remove log z
	    }
	  else if (name.Contains("vtx"))
	    {
	      //hist->GetZaxis()->SetRangeUser(1,10000);
	      canvas->SetCanvasSize(700,700);
	    }
	  else if (name.Contains("y_vs_eta_pt"))
	    {
	      //gStyle->SetOptStat(211);
	      continue;
	    }
	  else if (name.Contains("y_vs_eta"))
	    {
	      gStyle->SetOptStat(211);
	    }

	  hist->Draw("COLZ");
	  canvas->Update();
	}
      else if (cl->InheritsFrom("TH1"))
	{
	  hist->SetTitleOffset(1.2, "y");

	  if (name.Contains("sinAvgs") || name.Contains("cosAvgs") || name.Contains("Xn") || name.Contains("Yn"))
	    {
	      continue;
	    }
	  else if (name.Contains("vtx") || name.Contains("pT") || name.Contains("nhits") || name.Contains("mom")
		   || name.Contains("primTracks") || name.Contains("dndm") || name.Contains("tofBeta"))
	    {
	      canvas->SetLogy();
	    }
	  else if (name.Contains("p_") && !name.Contains("pp") && !name.Contains("kp"))
	    {
	      canvas->SetLeftMargin(0.15);
	      hist->SetTitleOffset(1.8, "y");
	      gStyle->SetOptStat(0);
	    }
	  else if (name.Contains("TOF_beta"))
	    {
	      canvas->SetLogy();
	      hist->GetYaxis()->SetRangeUser(0.1,10E+7);
	      gStyle->SetOptStat(11);
	    }
	  else if (name.Contains("Check"))
	    {
	      hist->SetFillColorAlpha(4,0.6);
	      hist->GetXaxis()->SetLabelSize(0.07);
	      if (name.Contains("trackCheck")) canvas->SetBottomMargin(0.15);
	    }
	  else if (name.Contains("dndm")) 
	    { 
	      hist->SetMinimum(0.1);
	      canvas->SetLogy();
	      errors = true;
	    }

	  
	  if (errors) 
	    hist->Draw("E1");
	  else 
	    hist->Draw("HIST");

	  canvas->Update();
	}
      else
	continue;
      
      canvas->SaveAs(jobID+"_"+name+".png");
    }

  canvas->Clear();
  file->Close();
}
