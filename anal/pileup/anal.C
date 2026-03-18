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
void anal(TString dataFilename = "AO2D.root")  // hyperloop
//void anal(TString dataFilename = "../omegaProdDerived/AO2D.root")
{
  ROOT::EnableImplicitMT();
  TChain dataChain(Form("O2nppileuptable"));
  fillChainFromAO2D(dataChain, dataFilename);
  ROOT::RDataFrame dataDF0(dataChain);
  auto c1 = dataDF0.Filter("fNumContrib").Count();
  std::cout << *c1 << std::endl;
  //
  //th1d->DrawClone();
  TFile output("output.root", "recreate");
  auto h2 = dataDF0.Histo2D({"h2", "Run vs NumContrib",
                    10000, 550000, 560000,
                    200, 0, 200},
                   "fRunNumber",
                   "fNumContrib");
  TH1D* hx = h2->ProjectionX("hx");
  h2->Write();
  hx->Write();
  auto dataDF = dataDF0.Filter("fMultFT0M > 000");
  //auto bcVec = dataDF.Filter("fRunNumber == 551843").Take<ULong64_t>("fGlobalBC");
  //
  auto bcVec = dataDF.Take<ULong64_t>("fGlobalBC");
  auto tracksVec = dataDF.Take<uint16_t>("fNumContrib");
  auto ft0mVec   = dataDF.Take<float>("fMultFT0M");
  std::unordered_map<ULong64_t, std::vector<size_t>> bcMap;
  //
  std::map<ULong64_t, int> bcCounts;
  int i = 0;
  for (auto bc : *bcVec) {
    bcMap[bc].push_back(i);
    bcCounts[bc]++;
    i++;
  }
  std::cout << bcVec->size() << ":" << bcCounts.size() << std::endl;
  TH1I* hNcollPerBC = new TH1I("hNcollPerBC", "Collisions per BC;N collisions in BC;Number of BCs", 10, 0.5, 10.5);
  for (auto const& [bc, n] : bcCounts) {
    hNcollPerBC->Fill(n);
  }
  TF1* fPois = new TF1("fPois", "[0]*TMath::Poisson(x,[1])", 0, 10);
  fPois->SetParameters(hNcollPerBC->Integral(), hNcollPerBC->GetMean());
  hNcollPerBC->Fit(fPois, "R");
  hNcollPerBC->Write();
  //
  auto hconts = std::make_unique<TH2D>(
    "NumContrib_sameBC2",
    "Contrib1 vs Contrib2 for collisions with exactly two same GlobalBC",
    100, 0, 100,
    100, 0, 100);
  auto hft0vstrk = std::make_unique<TH2D>(
    "hFT0MvsTracks_nopileup",
    "FT0M vs tracks for NC with one collision",
    1000, 0, 10000,
    100, 0, 100);
  for (const auto& [bc, idxs] : bcMap) {
    if (idxs.size() == 2) {
      hconts->Fill((*tracksVec)[idxs[0]], (*tracksVec)[idxs[1]]);
    }
    if (idxs.size() == 1) {
      hft0vstrk->Fill((*ft0mVec)[idxs[0]], (*tracksVec)[idxs[0]]);
    }
  }
  hft0vstrk->Write();
  hconts->Write();
  auto ht0vstrkP = hft0vstrk->ProfileX();
  ht0vstrkP->Write();
  //
 
}
