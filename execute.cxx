#include <iostream>
#include <string>
#include "TROOT.h"
#include "TString.h"

void execute(TString fileList, TString jobID, TString configFile)
{
  gROOT->ProcessLine(".L StRoot/StPicoEvent/libStPicoDst.so");
  gROOT->ProcessLine(".L StRoot/StEpdUtil/libStEpdUtil.so");
  gROOT->ProcessLine(".L ConfigReader_cxx.so");
  //gROOT->ProcessLine(".x FlowAnalyzer.cxx+\(\"$FILELIST\",\"$JOBID\",\"config_3p0GeV.txt\"\)");

  TString command = ".x FlowAnalyzer.cxx+\(\""+fileList+"\",\""+jobID+"\",\""+configFile+"\"\)";

  gROOT->ProcessLine(command);
}
