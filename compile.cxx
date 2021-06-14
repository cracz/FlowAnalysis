/**
 * \brief Helper script to compile the FlowAnalyzer.cxx macro
 *
 * This macro is meant to be used with 3.85 GeV production data, 
 * so before running this compile script, first run the commmand
 *     starver SL19b
 *
 **/

void compile()
{
  gROOT->ProcessLine(".L StRoot/StPicoEvent/libStPicoDst.so");
  gROOT->ProcessLine(".L StRoot/StEpdUtil/libStEpdUtil.so");
  gROOT->ProcessLine(".L ConfigReader_cxx.so");
  gROOT->ProcessLine(".L FlowAnalyzer.cxx++");
}
