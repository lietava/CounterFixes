void fillChainFromAO2D(TChain &chain, const TString &fileName)
{
  TFile file(fileName);
  for (auto key : *file.GetListOfKeys())
  {
    TString keyName = key->GetName();
    if (keyName.Contains("DF_"))
    {
      chain.Add((fileName + "/" + keyName + "/" + chain.GetName()).Data());
    }
  }
}
void analMultFromTree(TString dataTRG = "dataTRG/513390/AO2D.root")
{
  //ROOT::EnableImplicitMT();
  //
  ROOT::RDF::RResultPtr<TH1D> dataCentHist;
  ROOT::RDF::RResultPtr<TH1D> dataFT0MHist;

  TChain dataChain("O2npcasctablent");
  fillChainFromAO2D(dataChain, dataTRG);
  ROOT::RDataFrame dataDF(dataChain);
  auto dataFilteredDF = dataDF.Filter("fToiMask & 0x1 && fNoSameBunchPileup");
  dataCentHist =  dataFilteredDF.Histo1D({"Cent","Cent",4000,0.,10.},"fCentFT0M");
  dataFT0MHist =  dataFilteredDF.Histo1D({"FT0M","FT0M",7000,0.,7000.},"fMultFT0M");
  TFile output("outputMultTree.root", "recreate");
  dataCentHist->Write();
  dataFT0MHist->Write();
  //
  std::map<int, TH1*> centPerRunHist;
  int lastRun = 0;
  dataFilteredDF.Foreach([&lastRun, &centPerRunHist] (int runNumber, float cent, float mult) {
    //std::cout << runNumber << " " << cent << " " << mult << std::endl;
    if(lastRun != runNumber) {
      if(centPerRunHist.count(runNumber) == 0) {
        std::string histName = "cent" + to_string(runNumber);
        std::cout << "Creating:" << histName << std::endl;
        //std::cout << runNumber << ",";
        centPerRunHist[runNumber] = new TH1F(histName.c_str(), histName.c_str(), 4000, 0., 10.);
        lastRun = runNumber;
      }
    }
    //std::cout << "Filling:" << lastRun << std::endl;
    centPerRunHist[lastRun]->Fill(cent);

  },{ "fRunNumber", "fCentFT0M", "fMultFT0M"});

  for(auto const& h: centPerRunHist) {
    h.second->Write();
  }
}
