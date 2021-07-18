#include <iostream>
#include <fstream>
#include <sstream>
#include "ConfigReader.h"


ConfigReader::ConfigReader() 
{ errorFlag = false; }

ConfigReader::~ConfigReader()
{}

bool ConfigReader::errorFound()
{ return errorFlag; }

void ConfigReader::notifyError()
{
  std::cout << std::endl;
  std::cout << "There was an error in reading the config file." << std::endl;

  if (lastGoodKey != "")
    { std::cout << "The last key and value successflly read were: " << lastGoodKey << ", " << lastGoodValue << std::endl; }
  else
    { std::cout << "There were no keys or values read successfully." << std::endl; }
}

void ConfigReader::read(std::string fileName)
{
  const TString intValKeys[9] = {"epd_max_weight", "nHits", "dEdx", "min_tracks", 
				 "shift_terms", "epdA_inner_row", "epdA_outer_row", 
				 "epdB_inner_row", "epdB_outer_row"};
  const TString doubleValKeys[44] = {"sqrt_s_NN", "order_n", "order_m", "epd_threshold", 
				     "tracking", "dca", "min_abs_tpc_eta", 
				     "near_abs_tpc_eta", "far_abs_tpc_eta", 
				     "max_abs_tpc_eta", "r_vtx", 
				     "z_vtx_low", "z_vtx_high", "y_mid", "nSig_pi_low", 
				     "nSig_pi_high", "nSig_ka_low", "nSig_ka_high", 
				     "nSig_pr_low", "nSig_pr_high", "m2_pi_low", 
				     "m2_pi_high", "m2_ka_low", "m2_ka_high", 
				     "y_mid_pi_low_wide", "y_mid_pi_low", 
				     "y_mid_pi_high", "y_mid_pi_high_wide", 
				     "y_mid_ka_low_wide", "y_mid_ka_low", 
				     "y_mid_ka_high", "y_mid_ka_high_wide", 
				     "y_mid_pr_low_wide", "y_mid_pr_low", 
				     "y_mid_pr_high", "y_mid_pr_high_wide", 
				     "pt_pi_low", "pt_pi_high", 
				     "pt_ka_low", "pt_ka_high", 
				     "pt_pr_low_wide", "pt_pr_low", 
				     "pt_pr_high","pt_pr_high_wide"};

  int intValKeysLength = sizeof(intValKeys)/sizeof(intValKeys[0]);
  int doubleValKeysLength = sizeof(doubleValKeys)/sizeof(doubleValKeys[0]);


  std::ifstream inputStream(fileName.c_str());

  std::string line;
  std::getline(inputStream, line);  // Get the first line

  // Loop over lines of input in current file
  while (inputStream.good())
    {
      // Skip the text lines and the empty lines
      if (line[0] == '#' || line.empty())
	{
	  std::getline(inputStream, line);
	  continue;
	}

      // Split string by delimeter '='
      std::string delimeter = "=";
      size_t pos = line.find(delimeter);

      if (pos == std::string::npos) 
	{ 
	  std::cout << "Missing \'=\' in a line of the config file." << std::endl; 
	  errorFlag = true;
	  break;
	}

      TString key   = line.substr(0, pos);
      TString value = line.substr(pos+1, std::string::npos);
      
      lastGoodKey   = key;
      lastGoodValue = value;


      Int_t value_int;
      Double_t value_double;
      bool keyFound = false;

      // Check if the key that was just read is a known cut parameter and what type the value is.
      for (int i = 0; i < intValKeysLength; i++)
	{
	  if (key == intValKeys[i])
	    {
	      value_int = value.Atoi(); 
	      keyFound = true;
	      break;
	    }
	}

      if (keyFound) // Known parameter and it's an integer value
	{
	  if (key == "epd_max_weight") { epd_max_weight = value_int; }
	  else if (key == "nHits") { nHits = value_int; }
	  else if (key == "dEdx") { dEdx = value_int; }
	  else if (key == "min_tracks") { min_tracks = value_int; }
	  else if (key == "shift_terms") { shift_terms = value_int; }
	  else if (key == "epdA_inner_row") { epdA_inner_row = value_int; }
	  else if (key == "epdA_outer_row") { epdA_outer_row = value_int; }
	  else if (key == "epdB_inner_row") { epdB_inner_row = value_int; }
	  else if (key == "epdB_outer_row") { epdB_outer_row = value_int; }
	}
      else
	{
	  for (int i = 0; i < doubleValKeysLength; i++)
	    {
	      if (key == doubleValKeys[i])
		{
		  value_double = value.Atof(); 
		  keyFound = true;
		  break;
		}
	    }
	}

      if (keyFound) // Known parameter and it's a double value
	{
	  if (key == "sqrt_s_NN") { sqrt_s_NN = value_double; }
	  else if (key == "order_n") { order_n = value_double; order_n_str.Form("%d", (Int_t)value_double); }
	  else if (key == "order_m") { order_m = value_double; order_m_str.Form("%d", (Int_t)value_double); }
	  else if (key == "epd_threshold") { epd_threshold = value_double; }
	  else if (key == "tracking") { tracking = value_double; }
	  else if (key == "dca") { dca = value_double; }
	  else if (key == "min_abs_tpc_eta") { min_abs_tpc_eta = value_double; }
	  else if (key == "near_abs_tpc_eta") { near_abs_tpc_eta = value_double; }
	  else if (key == "far_abs_tpc_eta") { far_abs_tpc_eta = value_double; }
	  else if (key == "max_abs_tpc_eta") { max_abs_tpc_eta = value_double; }
	  else if (key == "r_vtx") { r_vtx = value_double; }
	  else if (key == "z_vtx_low") { z_vtx_low = value_double; }
	  else if (key == "z_vtx_high") { z_vtx_high = value_double; }
	  else if (key == "y_mid") { y_mid = value_double; }
	  else if (key == "nSig_pi_low") { nSig_pi_low = value_double; }
	  else if (key == "nSig_pi_high") { nSig_pi_high = value_double; }
	  else if (key == "nSig_ka_low") { nSig_ka_low = value_double; }
	  else if (key == "nSig_ka_high") { nSig_ka_high = value_double; }
	  else if (key == "nSig_pr_low") { nSig_pr_low = value_double; }
	  else if (key == "nSig_pr_high") { nSig_pr_high = value_double; }
	  else if (key == "m2_pi_low") { m2_pi_low = value_double; }
	  else if (key == "m2_pi_high") { m2_pi_high = value_double; }
	  else if (key == "m2_ka_low") { m2_ka_low = value_double; }
	  else if (key == "m2_ka_high") { m2_ka_high = value_double; }
	  else if (key == "y_mid_pi_low_wide") { y_mid_pi_low_wide = value_double; }
	  else if (key == "y_mid_pi_low") { y_mid_pi_low = value_double; }
	  else if (key == "y_mid_pi_high") { y_mid_pi_high = value_double; }
	  else if (key == "y_mid_pi_high_wide") { y_mid_pi_high_wide = value_double; }
	  else if (key == "y_mid_ka_low_wide") { y_mid_ka_low_wide = value_double; }
	  else if (key == "y_mid_ka_low") { y_mid_ka_low = value_double; }
	  else if (key == "y_mid_ka_high") { y_mid_ka_high = value_double; }
	  else if (key == "y_mid_ka_high_wide") { y_mid_ka_high_wide = value_double; }
	  else if (key == "y_mid_pr_low_wide") { y_mid_pr_low_wide = value_double; }
	  else if (key == "y_mid_pr_low") { y_mid_pr_low = value_double; }
	  else if (key == "y_mid_pr_high") { y_mid_pr_high = value_double; }
	  else if (key == "y_mid_pr_high_wide") { y_mid_pr_high_wide = value_double; }
	  else if (key == "pt_pi_low") { pt_pi_low = value_double; }
	  else if (key == "pt_pi_high") { pt_pi_high = value_double; }
	  else if (key == "pt_ka_low") { pt_ka_low = value_double; }
	  else if (key == "pt_ka_high") { pt_ka_high = value_double; }
	  else if (key == "pt_pr_low_wide") { pt_pr_low_wide = value_double; }
	  else if (key == "pt_pr_low") { pt_pr_low = value_double; }
	  else if (key == "pt_pr_high") { pt_pr_high = value_double; }
	  else if (key == "pt_pr_high_wide") { pt_pr_high_wide = value_double; }
	}
      else
	{ 
	  std::cout << "Unknown config key: " << key << std::endl;
	  errorFlag = true;
	  break;
	}

      std::getline(inputStream, line); // Get the next line
    }// End while(inputStream.good())

  if (errorFlag) { notifyError(); }
}// End function read()
