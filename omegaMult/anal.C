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
void anal(TString dataFilename = "dataMB/AO2D.root", TString normalisationFilename = "dataMB/AnalysisResults.root")
{
    ROOT::EnableImplicitMT();
    std::string names[2]{"", "nt"};
    int i = 1;
    TFile output("output.root", "recreate");
    auto dir = output.mkdir(Form("Omega%s", names[i].c_str()));
    dir->cd();
    TChain dataChain(Form("O2npcasctable%s", names[i].c_str()));
    fillChainFromAO2D(dataChain, dataFilename);
    ROOT::RDataFrame dataDF(dataChain);
    auto c1 = dataDF.Filter("fSel8").Count();
    auto c2 = dataDF.Filter("fMultFT0C").Count();
    auto c3 = dataDF.Filter("fMultFT0A").Count();
    auto c4 = dataDF.Filter("fProtonPt").Count();
    //
    //auto h2 = dataDF.Histo1D("fMultFT0C");
    auto h2 = dataDF.Filter("fMultFT0C >= 0").Histo1D("fMultFT0C");
    auto h3 = dataDF.Filter("fMultFT0A >= 0").Histo1D("fMultFT0A");

    //auto th1d = dataDF.Fill<float>(TH1D("th1d", "th1d", 64, 0, 128), {"fMultFT0A"});

    h3->DrawClone();
    //th1d->DrawClone();
    std::cout << *c1 << " " << *c2 << " " << *c3 << " " << *c4 <<std::endl;
}
