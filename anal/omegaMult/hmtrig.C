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
void hmtrig()
{
    //TString path = "../omegaProdDerived/";
    TString path = "dataTRG/559609/";
    //TString normalisationFilename = "dataTRG/557047/AnalysisResults.root";
    TString normalisationFilename = path + "AnalysisResults.root";
    TFile normalisationFile(normalisationFilename.Data());
    ZorroSummary *zorroSummary = (ZorroSummary *)normalisationFile.Get("non-prompt-cascade-task/ZorroSummary");
    std::vector<std::string> tois;
    if(zorroSummary == nullptr) {
        std::cout << "Can not get zorro" << std::endl;
        exit(1);
    } else {
        std::string text = zorroSummary->getTOInames();
        tois = o2::utils::Str::tokenize(text, ',');
        std::cout << "TOI:" <<  text << ":" << tois.size() << std::endl;
    }
    //
    TString fileName = path + "AO2D.root";
    TChain dataChain("O2npcasctablent");
    fillChainFromAO2D(dataChain, fileName);
    ROOT::RDataFrame dataDF(dataChain);
    for(int i = 0; i < tois.size(); i++) {
        uint32_t mask = 1 << i;
        auto nsel = dataDF.Filter(Form("fToiMask & %u", mask)).Count();
        std::cout << tois[i] << " " << *nsel << std::endl;
    }
    auto h1 = dataDF.Filter("fToiMask & 0x40")
           .Histo1D({"hMult", "Mult distribution", 10000, 0, 10000}, "fMultFT0M");
    auto h2 = dataDF.Filter("fToiMask & 0x40")
           .Histo1D({"hTracks", "Tracks distribution", 100, 0, 100}, "fMultNTracksGlobal");
    TFile outfile("outtest.root","recreate");
    outfile.cd();
    h1->Write();
    h2->Write();
}
