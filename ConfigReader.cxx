#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "ConfigReader.h"

// I have to initialize the maps like this because c++98 sucks (or maybe just the CINT does)
void ConfigReader::initialize()
{
  intValCuts["minbias"] = -999;
  intValCuts["epd_max_weight"] = -999;
  intValCuts["nHits"] = -999;
  intValCuts["dEdx"] = -999;
  intValCuts["min_tracks"] = -999;
  intValCuts["shift_terms"] = -999;
  intValCuts["epdA_inner_row"] = -999;
  intValCuts["epdA_outer_row"] = -999;
  intValCuts["epdB_inner_row"] = -999;
  intValCuts["epdB_outer_row"] = -999;

  dblValCuts["sqrt_s_NN"] = -999.0;
  dblValCuts["order_n"] = -999.0; 
  dblValCuts["order_m"] = -999.0; 
  dblValCuts["epd_threshold"] = -999.0; 
  dblValCuts["tracking"] = -999.0; 
  dblValCuts["dca"] = -999.0; 
  dblValCuts["min_abs_tpc_eta"] = -999.0; 
  dblValCuts["near_abs_tpc_eta"] = -999.0; 
  dblValCuts["far_abs_tpc_eta"] = -999.0; 
  dblValCuts["max_abs_tpc_eta"] = -999.0; 
  dblValCuts["r_vtx"] = -999.0; 
  dblValCuts["z_vtx_low"] = -999.0; 
  dblValCuts["z_vtx_high"] = -999.0; 
  dblValCuts["y_mid"] = -999.0; 
  dblValCuts["nSig_pi_low"] = -999.0; 
  dblValCuts["nSig_pi_high"] = -999.0; 
  dblValCuts["nSig_ka_low"] = -999.0; 
  dblValCuts["nSig_ka_high"] = -999.0; 
  dblValCuts["nSig_pr_low"] = -999.0; 
  dblValCuts["nSig_pr_high"] = -999.0; 
  dblValCuts["m2_pi_low"] = -999.0; 
  dblValCuts["m2_pi_high"] = -999.0; 
  dblValCuts["m2_ka_low"] = -999.0; 
  dblValCuts["m2_ka_high"] = -999.0; 
  dblValCuts["yCM_pid_pi_low"] = -999.0;
  dblValCuts["yCM_pid_pi_high"] = -999.0;
  dblValCuts["yCM_flow_pi_low"] = -999.0;
  dblValCuts["yCM_flow_pi_high"] = -999.0;
  dblValCuts["yCM_ext_flow_pi_low"] = -999.0;
  dblValCuts["yCM_ext_flow_pi_high"] = -999.0;
  dblValCuts["yCM_pid_ka_low"] = -999.0;
  dblValCuts["yCM_pid_ka_high"] = -999.0;
  dblValCuts["yCM_flow_ka_low"] = -999.0;
  dblValCuts["yCM_flow_ka_high"] = -999.0;
  dblValCuts["yCM_ext_flow_ka_low"] = -999.0;
  dblValCuts["yCM_ext_flow_ka_high"] = -999.0;
  dblValCuts["yCM_pid_pr_low"] = -999.0;
  dblValCuts["yCM_pid_pr_high"] = -999.0;
  dblValCuts["yCM_flow_pr_low"] = -999.0;
  dblValCuts["yCM_flow_pr_high"] = -999.0;
  dblValCuts["yCM_dep_flow_pr_low"] = -999.0;
  dblValCuts["yCM_dep_flow_pr_high"] = -999.0;
  dblValCuts["yCM_ext_flow_pr_low"] = -999.0;
  dblValCuts["yCM_ext_flow_pr_high"] = -999.0;
  dblValCuts["yCM_sym_flow_pr_low"] = -999.0;
  dblValCuts["yCM_sym_flow_pr_high"] = -999.0;
  dblValCuts["yCM_for_flow_pr_low"] = -999.0;
  dblValCuts["yCM_for_flow_pr_high"] = -999.0;
  dblValCuts["pt_pid_pi_low"] = -999.0; 
  dblValCuts["pt_pid_pi_high"] = -999.0; 
  dblValCuts["pt_pid_ka_low"] = -999.0; 
  dblValCuts["pt_pid_ka_high"] = -999.0; 
  dblValCuts["pt_pid_pr_low"] = -999.0; 
  dblValCuts["pt_pid_pr_high"] = -999.0; 
  dblValCuts["pt_flow_pr_low"] = -999.0; 
  dblValCuts["pt_flow_pr_high"] = -999.0; 
  dblValCuts["pt_ydep_flow_pr_low"] = -999.0;
  dblValCuts["pt_ydep_flow_pr_high"] = -999.0;
  dblValCuts["pt_ext_flow_pr_low"] = -999.0;
  dblValCuts["pt_ext_flow_pr_high"] = -999.0;
  dblValCuts["pt_sym_flow_pr_low"] = -999.0;
  dblValCuts["pt_sym_flow_pr_high"] = -999.0;
  dblValCuts["pt_for_flow_pr_low"] = -999.0;
  dblValCuts["pt_for_flow_pr_high"] = -999.0;
}

void ConfigReader::setAllCuts()
{
  minbias = intValCuts["minbias"];
  epd_max_weight = intValCuts["epd_max_weight"];
  nHits = intValCuts["nHits"];
  dEdx = intValCuts["dEdx"];
  min_tracks = intValCuts["min_tracks"];
  shift_terms = intValCuts["shift_terms"];
  epdA_inner_row = intValCuts["epdA_inner_row"];
  epdA_outer_row = intValCuts["epdA_outer_row"];
  epdB_inner_row = intValCuts["epdB_inner_row"];
  epdB_outer_row = intValCuts["epdB_outer_row"];

  sqrt_s_NN = dblValCuts["sqrt_s_NN"];
  order_n = dblValCuts["order_n"]; 
  order_m = dblValCuts["order_m"]; 
  epd_threshold = dblValCuts["epd_threshold"]; 
  tracking = dblValCuts["tracking"]; 
  dca = dblValCuts["dca"]; 
  min_abs_tpc_eta = dblValCuts["min_abs_tpc_eta"]; 
  near_abs_tpc_eta = dblValCuts["near_abs_tpc_eta"]; 
  far_abs_tpc_eta = dblValCuts["far_abs_tpc_eta"]; 
  max_abs_tpc_eta = dblValCuts["max_abs_tpc_eta"]; 
  r_vtx = dblValCuts["r_vtx"]; 
  z_vtx_low = dblValCuts["z_vtx_low"]; 
  z_vtx_high = dblValCuts["z_vtx_high"]; 
  y_mid = dblValCuts["y_mid"]; 
  nSig_pi_low = dblValCuts["nSig_pi_low"]; 
  nSig_pi_high = dblValCuts["nSig_pi_high"]; 
  nSig_ka_low = dblValCuts["nSig_ka_low"]; 
  nSig_ka_high = dblValCuts["nSig_ka_high"]; 
  nSig_pr_low = dblValCuts["nSig_pr_low"]; 
  nSig_pr_high = dblValCuts["nSig_pr_high"]; 
  m2_pi_low = dblValCuts["m2_pi_low"]; 
  m2_pi_high = dblValCuts["m2_pi_high"]; 
  m2_ka_low = dblValCuts["m2_ka_low"]; 
  m2_ka_high = dblValCuts["m2_ka_high"]; 
  yCM_pid_pi_low = dblValCuts["yCM_pid_pi_low"];
  yCM_pid_pi_high = dblValCuts["yCM_pid_pi_high"];
  yCM_flow_pi_low = dblValCuts["yCM_flow_pi_low"];
  yCM_flow_pi_high = dblValCuts["yCM_flow_pi_high"];
  yCM_ext_flow_pi_low = dblValCuts["yCM_ext_flow_pi_low"];
  yCM_ext_flow_pi_high = dblValCuts["yCM_ext_flow_pi_high"];
  yCM_pid_ka_low = dblValCuts["yCM_pid_ka_low"];
  yCM_pid_ka_high = dblValCuts["yCM_pid_ka_high"];
  yCM_flow_ka_low = dblValCuts["yCM_flow_ka_low"];
  yCM_flow_ka_high = dblValCuts["yCM_flow_ka_high"];
  yCM_ext_flow_ka_low = dblValCuts["yCM_ext_flow_ka_low"];
  yCM_ext_flow_ka_high = dblValCuts["yCM_ext_flow_ka_high"];
  yCM_pid_pr_low = dblValCuts["yCM_pid_pr_low"];
  yCM_pid_pr_high = dblValCuts["yCM_pid_pr_high"];
  yCM_flow_pr_low = dblValCuts["yCM_flow_pr_low"];
  yCM_flow_pr_high = dblValCuts["yCM_flow_pr_high"];
  yCM_dep_flow_pr_low = dblValCuts["yCM_dep_flow_pr_low"];
  yCM_dep_flow_pr_high = dblValCuts["yCM_dep_flow_pr_high"];
  yCM_ext_flow_pr_low = dblValCuts["yCM_ext_flow_pr_low"];
  yCM_ext_flow_pr_high = dblValCuts["yCM_ext_flow_pr_high"];
  yCM_sym_flow_pr_low = dblValCuts["yCM_sym_flow_pr_low"];
  yCM_sym_flow_pr_high = dblValCuts["yCM_sym_flow_pr_high"];
  yCM_for_flow_pr_low = dblValCuts["yCM_for_flow_pr_low"];
  yCM_for_flow_pr_high = dblValCuts["yCM_for_flow_pr_high"];
  pt_pid_pi_low = dblValCuts["pt_pid_pi_low"]; 
  pt_pid_pi_high = dblValCuts["pt_pid_pi_high"]; 
  pt_pid_ka_low = dblValCuts["pt_pid_ka_low"]; 
  pt_pid_ka_high = dblValCuts["pt_pid_ka_high"]; 
  pt_pid_pr_low = dblValCuts["pt_pid_pr_low"]; 
  pt_pid_pr_high = dblValCuts["pt_pid_pr_high"]; 
  pt_flow_pr_low = dblValCuts["pt_flow_pr_low"]; 
  pt_flow_pr_high = dblValCuts["pt_flow_pr_high"]; 
  pt_ydep_flow_pr_low = dblValCuts["pt_ydep_flow_pr_low"];
  pt_ydep_flow_pr_high = dblValCuts["pt_ydep_flow_pr_high"];
  pt_ext_flow_pr_low = dblValCuts["pt_ext_flow_pr_low"];
  pt_ext_flow_pr_high = dblValCuts["pt_ext_flow_pr_high"];
  pt_sym_flow_pr_low = dblValCuts["pt_sym_flow_pr_low"];
  pt_sym_flow_pr_high = dblValCuts["pt_sym_flow_pr_high"];
  pt_for_flow_pr_low = dblValCuts["pt_for_flow_pr_low"];
  pt_for_flow_pr_high = dblValCuts["pt_for_flow_pr_high"];
}

ConfigReader::ConfigReader() 
{ 
  initialize();
  errorFlag = false; 
}

ConfigReader::~ConfigReader()
{}

bool ConfigReader::errorFound()
{ return errorFlag; }

void ConfigReader::notifyError()
{
  std::cout << std::endl;
  std::cout << "There was an error in reading the config file." << std::endl;

  if (lastKey != "")
    { std::cout << "The last key and value read were: " << lastKey << ", " << lastValue << std::endl; }
  else
    { std::cout << "There were no keys or values read successfully." << std::endl; }
}

void ConfigReader::read(std::string fileName)
{
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

      std::string key   = line.substr(0, pos);
      std::string value = line.substr(pos+1, std::string::npos);
      
      lastKey   = key;
      lastValue = value;

      try
	{ intValCuts.at(key) = std::atoi(value.c_str()); }
      catch (...)//(const std::out_of_range& oorInt) 
	{
	  try
	    { dblValCuts.at(key) = std::atof(value.c_str()); }
	  catch (...)//(const std::out_of_range& oorDbl) 
	    {
	      std::cout << "Unknown config key: " << key << std::endl;
	      errorFlag = true;
	      break;	      
	    }
	}

      std::getline(inputStream, line); // Get the next line
    }// End while(inputStream.good())

  setAllCuts();

  if (errorFlag) { notifyError(); }
}// End function read()
