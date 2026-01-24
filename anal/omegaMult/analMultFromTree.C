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
  //auto dataFilteredDF = dataDF.Filter("fToiMask & 0x1");
  TFile output("outputMultTree.root", "recreate");
  // dataCentHist =  dataFilteredDF.Histo1D({"Cent","Cent",4000,0.,10.},"fCentFT0M");
  // dataFT0MHist =  dataFilteredDF.Histo1D({"FT0M","FT0M",7000,0.,7000.},"fMultFT0M");
  // dataCentHist->Write();
  // dataFT0MHist->Write();
  auto h2CentVsFT0M = dataFilteredDF.Histo2D(
    { "hCentVsFT0M", "Centrality vs FT0M multiplicity",
      4000, 0., 10.,     // X bins for Centrality (Cent)
      7000, 0., 7000. }, // Y bins for FT0M multiplicity
    "fCentFT0M", "fMultFT0M"
  );
  auto h2CentVsNTrk = dataFilteredDF.Histo2D(
    { "hCentVsNTrk", "Centrality vs FT0M multiplicity",
      4000, 0., 10.,     // X bins for Centrality (Cent)
      100, 0., 100. }, // Y bins for FT0M multiplicity
    "fCentFT0M", "fMultNTracksGlobal"
  );
  auto pCentVsMeanMult = h2CentVsFT0M->ProfileX("pCentVsMeanMult");
  auto pCentVsNTrk = h2CentVsNTrk->ProfileX("pCentVsNTrk");
  h2CentVsFT0M->Write();
  pCentVsMeanMult->Write();
  h2CentVsNTrk->Write();
  pCentVsNTrk->Write();
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
