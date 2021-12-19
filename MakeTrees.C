#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;

Int_t MakeTrees(Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir, TString JobIdName, std::string configFileName, Double_t testParameter)
{
  // Load libraries
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
    
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StEpdUtil");
  gSystem->Load("StPileupUtil");
  gSystem->Load("TreeMaker.so");
  gSystem->Load("ConfigReader.so");

  std::cout << "config file name: " << configFileName << std::endl;
  
  // List of member links in the chain
  StChain* chain = new StChain;
  StPicoDstMaker* picoMaker = new StPicoDstMaker(2,InputFileList,"picoDst");
  TreeMaker* doAnalysis = new TreeMaker(picoMaker,configFileName,JobIdName,nEvents,testParameter);
  
  if (nEvents == 0)  nEvents = 100000000 ;       // Take all events in nFiles if nEvents = 0
    
  // Loop over the links in the chain
  if(chain->Init() == kStErr)
    {
      cout<<"chain->Init();"<<endl;
      return;
    }
    
  int total = picoMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents > total) nEvents = total;
    
  for (Int_t i = 0; i < nEvents; i++)
    {

      if(i % 10000 == 0) cout << "Working on eventNumber " << i << endl;
      
      chain->Clear();
      int iret = chain->Make(i);      
      if (iret) { cout << "Bad return code!" << iret << endl; break;}

      total++;
    }
    
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!" << endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
    
  // Cleanup
  delete doAnalysis;
  delete picoMaker;
  delete chain;
  return 0;
}
