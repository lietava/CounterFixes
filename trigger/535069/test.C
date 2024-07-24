#include <iostream>
#include <TFile.h>
#include <TH1.h>

void test()
{
  TFile file("AnalysisResults_fullrun.root");
  if(!file.IsOpen()) {
    std::cout << "File AnalysisResults.root can not be opned" << std::endl;
    return;
  }
  // Extract the histograms
  TH1* hist1 = dynamic_cast<TH1*>(file.Get("central-event-filter-task/scalers/mFiltered;1")); // Replace with the correct path
  if(hist1 == nullptr) {
	  std::cout << " Can not find labels" << std::endl;
    return;
  }
  TH1* mFilterCounters = nullptr;  
  mFilterCounters = hist1;
  int n = mFilterCounters->GetNbinsX();
  for(int i = 0; i < n; i++) {
    cout << hist1->GetBinContent(i) << std::endl;
  }
  return;
}
